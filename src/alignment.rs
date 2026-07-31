use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};

use noodles::core::Region;
use noodles::sam::alignment::{Record, RecordBuf};
use noodles::sam::header::record::value::map::read_group::tag;
use noodles::{bam, bgzf, cram, fasta, sam};
use rsomics_bamio::raw::{self, RawRecord, RawRecordEncoder};
use rsomics_bamio::{IndexedAlignmentReader, open_indexed_alignment};

use crate::{CallError, Result, SampleMap, SampleMapBuilder, SampleSelection};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignmentInput {
    source_id: u32,
    path: PathBuf,
    sample_name: Box<str>,
    ignore_read_groups: bool,
}

impl AlignmentInput {
    pub fn new(source_id: u32, path: impl Into<PathBuf>, sample_name: impl Into<Box<str>>) -> Self {
        Self {
            source_id,
            path: path.into(),
            sample_name: sample_name.into(),
            ignore_read_groups: false,
        }
    }

    pub fn ignore_read_groups(mut self, ignore: bool) -> Self {
        self.ignore_read_groups = ignore;
        self
    }

    pub fn source_id(&self) -> u32 {
        self.source_id
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceSequence {
    name: Box<[u8]>,
    length: u64,
}

impl ReferenceSequence {
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    pub fn length(&self) -> u64 {
        self.length
    }
}

pub struct AlignmentSet {
    readers: Vec<SourceReader>,
    pending: BinaryHeap<PendingRecord>,
    references: Box<[ReferenceSequence]>,
    samples: SampleMap,
}

impl AlignmentSet {
    pub fn open(
        inputs: impl IntoIterator<Item = AlignmentInput>,
        reference: Option<&Path>,
        selection: SampleSelection,
    ) -> Result<Self> {
        let repository = reference.map(reference_repository).transpose()?;
        let mut readers = Vec::new();
        let mut metadata = AlignmentMetadataBuilder::new(selection);

        for input in inputs {
            let reader = SourceReader::open(input, repository.clone())?;
            metadata.add_source(
                reader.source_id,
                &reader.path,
                &reader.header,
                reader.sample_name.clone(),
                reader.ignore_read_groups,
            )?;
            readers.push(reader);
        }

        if readers.is_empty() {
            return Err(CallError::MissingAlignmentInputs);
        }
        let (references, samples) = metadata.finish()?;

        let mut set = Self {
            readers,
            pending: BinaryHeap::new(),
            references,
            samples,
        };
        for source_index in 0..set.readers.len() {
            set.fill_source(source_index)?;
        }
        Ok(set)
    }

    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.references
    }

    pub fn samples(&self) -> &SampleMap {
        &self.samples
    }

    pub fn next_record(&mut self) -> Result<Option<(u32, RawRecord)>> {
        let Some(pending) = self.pending.pop() else {
            return Ok(None);
        };
        let source_id = self.readers[pending.source_index].source_id;
        self.fill_source(pending.source_index)?;
        Ok(Some((source_id, pending.record)))
    }

    fn fill_source(&mut self, source_index: usize) -> Result<()> {
        if let Some(record) = self.readers[source_index].next_record()? {
            self.pending.push(PendingRecord::new(source_index, record));
        }
        Ok(())
    }
}

pub(crate) struct IndexedAlignmentSet {
    readers: Vec<IndexedSourceReader>,
    references: Box<[ReferenceSequence]>,
    samples: SampleMap,
}

impl IndexedAlignmentSet {
    pub(crate) fn open(
        inputs: impl IntoIterator<Item = AlignmentInput>,
        reference: Option<&Path>,
        selection: SampleSelection,
    ) -> Result<Self> {
        let mut readers = Vec::new();
        let mut metadata = AlignmentMetadataBuilder::new(selection);

        for input in inputs {
            let mut reader = open_indexed_alignment(&input.path, reference)
                .map_err(|error| input_error(&input.path, error))?;
            let header = reader
                .read_header()
                .map_err(|error| input_error(&input.path, error))?;
            metadata.add_source(
                input.source_id,
                &input.path,
                &header,
                input.sample_name.clone(),
                input.ignore_read_groups,
            )?;
            readers.push(IndexedSourceReader {
                source_id: input.source_id,
                path: input.path,
                header,
                reader,
            });
        }

        if readers.is_empty() {
            return Err(CallError::MissingAlignmentInputs);
        }
        let (references, samples) = metadata.finish()?;
        Ok(Self {
            readers,
            references,
            samples,
        })
    }

    pub(crate) fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.references
    }

    pub(crate) fn samples(&self) -> &SampleMap {
        &self.samples
    }

    pub(crate) fn visit_region(
        &mut self,
        region: &Region,
        mut visit: impl FnMut(&SampleMap, u32, RawRecord) -> Result<()>,
    ) -> Result<()> {
        let samples = &self.samples;
        let mut readers = self
            .readers
            .iter_mut()
            .map(|source| RegionSourceReader::new(source, region))
            .collect::<Result<Vec<_>>>()?;
        let mut pending = BinaryHeap::new();
        for source_index in 0..readers.len() {
            fill_region_source(&mut readers, &mut pending, source_index)?;
        }
        while let Some(record) = pending.pop() {
            let source_id = readers[record.source_index].source_id;
            fill_region_source(&mut readers, &mut pending, record.source_index)?;
            visit(samples, source_id, record.record)?;
        }
        Ok(())
    }
}

struct IndexedSourceReader {
    source_id: u32,
    path: PathBuf,
    header: sam::Header,
    reader: IndexedAlignmentReader,
}

struct RegionSourceReader<'a> {
    source_id: u32,
    path: &'a Path,
    header: &'a sam::Header,
    records: Box<dyn Iterator<Item = io::Result<Box<dyn Record>>> + 'a>,
    encoder: RawRecordEncoder,
}

impl<'a> RegionSourceReader<'a> {
    fn new(source: &'a mut IndexedSourceReader, region: &Region) -> Result<Self> {
        let records = source
            .reader
            .query(&source.header, region)
            .map_err(|error| input_error(&source.path, error))?;
        Ok(Self {
            source_id: source.source_id,
            path: &source.path,
            header: &source.header,
            records: Box::new(records),
            encoder: RawRecordEncoder::default(),
        })
    }

    fn next_record(&mut self) -> Result<Option<RawRecord>> {
        let Some(record) = self.records.next() else {
            return Ok(None);
        };
        let record = record.map_err(|error| input_error(self.path, error))?;
        self.encoder
            .encode(self.header, record.as_ref())
            .map(Some)
            .map_err(|error| input_error(self.path, error))
    }
}

fn fill_region_source(
    readers: &mut [RegionSourceReader<'_>],
    pending: &mut BinaryHeap<PendingRecord>,
    source_index: usize,
) -> Result<()> {
    if let Some(record) = readers[source_index].next_record()? {
        pending.push(PendingRecord::new(source_index, record));
    }
    Ok(())
}

struct AlignmentMetadataBuilder {
    references: Option<Vec<ReferenceSequence>>,
    samples: SampleMapBuilder,
}

impl AlignmentMetadataBuilder {
    fn new(selection: SampleSelection) -> Self {
        Self {
            references: None,
            samples: SampleMapBuilder::new(selection),
        }
    }

    fn add_source(
        &mut self,
        source_id: u32,
        path: &Path,
        header: &sam::Header,
        sample_name: Box<str>,
        ignore_read_groups: bool,
    ) -> Result<()> {
        let current_references = reference_sequences(header);
        if let Some(expected) = &self.references {
            if expected != &current_references {
                return Err(CallError::ReferenceDictionaryMismatch(
                    path.display().to_string(),
                ));
            }
        } else {
            self.references = Some(current_references);
        }

        let read_groups = header
            .read_groups()
            .iter()
            .filter_map(|(id, group)| {
                group
                    .other_fields()
                    .get(&tag::SAMPLE)
                    .map(|sample| (id, sample))
            })
            .map(|(id, sample)| {
                let sample = std::str::from_utf8(sample.as_ref())
                    .map_err(|error| input_error(path, error))?;
                Ok((id.to_vec(), sample.to_owned()))
            })
            .collect::<Result<Vec<_>>>()?;
        self.samples
            .add_source(source_id, sample_name, ignore_read_groups, read_groups)
    }

    fn finish(self) -> Result<(Box<[ReferenceSequence]>, SampleMap)> {
        Ok((
            self.references.unwrap().into_boxed_slice(),
            self.samples.finish()?,
        ))
    }
}

struct SourceReader {
    source_id: u32,
    path: PathBuf,
    sample_name: Box<str>,
    ignore_read_groups: bool,
    header: sam::Header,
    inner: SourceInner,
}

impl SourceReader {
    fn open(input: AlignmentInput, repository: Option<fasta::Repository>) -> Result<Self> {
        let (format, compression) = detect_source(&input.path)?;
        let file = File::open(&input.path).map_err(|error| input_error(&input.path, error))?;
        let mut inner = match (format, compression) {
            (Format::Sam, Compression::None) => SourceInner::Sam {
                reader: sam::io::Reader::new(BufReader::new(file)),
                record: RecordBuf::default(),
                encoder: RawRecordEncoder::default(),
            },
            (Format::Sam, Compression::Bgzf) => SourceInner::SamGz {
                reader: sam::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file))),
                record: RecordBuf::default(),
                encoder: RawRecordEncoder::default(),
            },
            (Format::Bam, Compression::None) => {
                SourceInner::BamRaw(bam::io::Reader::from(BufReader::new(file)))
            }
            (Format::Bam, Compression::Bgzf) => {
                SourceInner::Bam(bam::io::Reader::new(BufReader::new(file)))
            }
            (Format::Cram, Compression::None) => {
                let repository = repository.unwrap_or_default();
                let reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repository.clone())
                    .build_from_reader(BufReader::new(file));
                SourceInner::Cram {
                    reader,
                    repository,
                    container: cram::io::reader::Container::default(),
                    records: VecDeque::new(),
                    encoder: RawRecordEncoder::default(),
                }
            }
            (Format::Cram, Compression::Bgzf) => unreachable!(),
        };
        let header = inner
            .read_header()
            .map_err(|error| input_error(&input.path, error))?;

        Ok(Self {
            source_id: input.source_id,
            path: input.path,
            sample_name: input.sample_name,
            ignore_read_groups: input.ignore_read_groups,
            header,
            inner,
        })
    }

    fn next_record(&mut self) -> Result<Option<RawRecord>> {
        self.inner
            .next_record(&self.header)
            .map_err(|error| input_error(&self.path, error))
    }
}

enum SourceInner {
    Sam {
        reader: sam::io::Reader<BufReader<File>>,
        record: RecordBuf,
        encoder: RawRecordEncoder,
    },
    SamGz {
        reader: sam::io::Reader<bgzf::io::Reader<BufReader<File>>>,
        record: RecordBuf,
        encoder: RawRecordEncoder,
    },
    Bam(bam::io::Reader<bgzf::io::Reader<BufReader<File>>>),
    BamRaw(bam::io::Reader<BufReader<File>>),
    Cram {
        reader: cram::io::Reader<BufReader<File>>,
        repository: fasta::Repository,
        container: cram::io::reader::Container,
        records: VecDeque<RecordBuf>,
        encoder: RawRecordEncoder,
    },
}

impl SourceInner {
    fn read_header(&mut self) -> std::io::Result<sam::Header> {
        match self {
            Self::Sam { reader, .. } => reader.read_header(),
            Self::SamGz { reader, .. } => reader.read_header(),
            Self::Bam(reader) => reader.read_header(),
            Self::BamRaw(reader) => reader.read_header(),
            Self::Cram { reader, .. } => reader.read_header(),
        }
    }

    fn next_record(&mut self, header: &sam::Header) -> std::io::Result<Option<RawRecord>> {
        match self {
            Self::Sam {
                reader,
                record,
                encoder,
            } => {
                if reader.read_record_buf(header, record)? == 0 {
                    Ok(None)
                } else {
                    encode_record(encoder, header, record).map(Some)
                }
            }
            Self::SamGz {
                reader,
                record,
                encoder,
            } => {
                if reader.read_record_buf(header, record)? == 0 {
                    Ok(None)
                } else {
                    encode_record(encoder, header, record).map(Some)
                }
            }
            Self::Bam(reader) => read_bam_record(reader.get_mut()),
            Self::BamRaw(reader) => read_bam_record(reader.get_mut()),
            Self::Cram {
                reader,
                repository,
                container,
                records,
                encoder,
            } => {
                while records.is_empty() {
                    if reader.read_container(container)? == 0 {
                        return Ok(None);
                    }
                    let compression_header = container.compression_header()?;
                    for slice in container.slices() {
                        let slice = slice?;
                        let (core, external) = slice.decode_blocks()?;
                        let decoded = slice.records(
                            repository.clone(),
                            header,
                            &compression_header,
                            &core,
                            &external,
                        )?;
                        for record in decoded {
                            records
                                .push_back(RecordBuf::try_from_alignment_record(header, &record)?);
                        }
                    }
                }
                encode_record(encoder, header, &records.pop_front().unwrap()).map(Some)
            }
        }
    }
}

fn encode_record(
    encoder: &mut RawRecordEncoder,
    header: &sam::Header,
    record: &dyn sam::alignment::Record,
) -> io::Result<RawRecord> {
    encoder
        .encode(header, record)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
}

struct PendingRecord {
    source_index: usize,
    coordinate: (i32, i32),
    record: RawRecord,
}

impl PendingRecord {
    fn new(source_index: usize, record: RawRecord) -> Self {
        let coordinate = if record.flags() & 0x04 != 0 {
            (i32::MAX, i32::MAX)
        } else {
            (record.reference_sequence_id(), record.alignment_start())
        };
        Self {
            source_index,
            coordinate,
            record,
        }
    }
}

impl PartialEq for PendingRecord {
    fn eq(&self, other: &Self) -> bool {
        (self.coordinate, self.source_index) == (other.coordinate, other.source_index)
    }
}

impl Eq for PendingRecord {}

impl PartialOrd for PendingRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PendingRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.coordinate, other.source_index).cmp(&(self.coordinate, self.source_index))
    }
}

#[derive(Clone, Copy)]
enum Format {
    Sam,
    Bam,
    Cram,
}

#[derive(Clone, Copy)]
enum Compression {
    None,
    Bgzf,
}

fn detect_source(path: &Path) -> Result<(Format, Compression)> {
    let mut source = BufReader::new(File::open(path).map_err(|error| input_error(path, error))?);
    let mut magic = [0; 4];
    source
        .read_exact(&mut magic)
        .map_err(|error| input_error(path, error))?;
    if magic == *b"CRAM" {
        return Ok((Format::Cram, Compression::None));
    }
    if magic == *b"BAM\x01" {
        return Ok((Format::Bam, Compression::None));
    }
    if magic[..2] != [0x1f, 0x8b] {
        return Ok((Format::Sam, Compression::None));
    }

    let file = File::open(path).map_err(|error| input_error(path, error))?;
    let mut reader = bgzf::io::Reader::new(file);
    reader
        .read_exact(&mut magic)
        .map_err(|error| input_error(path, error))?;
    Ok((
        if magic == *b"BAM\x01" {
            Format::Bam
        } else {
            Format::Sam
        },
        Compression::Bgzf,
    ))
}

fn reference_repository(path: &Path) -> Result<fasta::Repository> {
    fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .map_err(|error| input_error(path, error))
}

fn reference_sequences(header: &sam::Header) -> Vec<ReferenceSequence> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| ReferenceSequence {
            name: name.to_vec().into_boxed_slice(),
            length: usize::from(reference.length()) as u64,
        })
        .collect()
}

fn read_bam_record(reader: &mut impl Read) -> std::io::Result<Option<RawRecord>> {
    let mut record = RawRecord::default();
    match raw::read_record(reader, &mut record) {
        Ok(0) => Ok(None),
        Ok(_) => Ok(Some(record)),
        Err(error) => Err(std::io::Error::new(std::io::ErrorKind::InvalidData, error)),
    }
}

fn input_error(path: &Path, error: impl std::fmt::Display) -> CallError {
    CallError::AlignmentInput {
        path: path.display().to_string(),
        message: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use noodles::sam::alignment::io::Write as _;
    use tempfile::NamedTempFile;

    use super::*;

    const FIRST: &str = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20\n\
@RG\tID:rg1\tSM:S1\n\
r1\t0\tchr1\t3\t60\t1M\t*\t0\t0\tA\tI\tRG:Z:rg1\n";
    const SECOND: &str = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20\n\
@RG\tID:rg2\tSM:S2\n\
r2\t0\tchr1\t2\t60\t1M\t*\t0\t0\tG\tI\tRG:Z:rg2\n";

    fn sam_file(contents: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        file
    }

    fn bam_file(contents: &str) -> NamedTempFile {
        let mut reader = sam::io::Reader::new(contents.as_bytes());
        let header = reader.read_header().unwrap();
        let records = reader
            .record_bufs(&header)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        let file = NamedTempFile::new().unwrap();
        let mut writer = bam::io::Writer::new(file.reopen().unwrap());
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
        file
    }

    fn raw_bam_file(contents: &str) -> NamedTempFile {
        let mut reader = sam::io::Reader::new(contents.as_bytes());
        let header = reader.read_header().unwrap();
        let records = reader
            .record_bufs(&header)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        let file = NamedTempFile::new().unwrap();
        let mut writer = bam::io::Writer::from(file.reopen().unwrap());
        writer.write_header(&header).unwrap();
        for record in records {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        file
    }

    fn bgzf_sam_file(contents: &str) -> NamedTempFile {
        let file = NamedTempFile::new().unwrap();
        let mut writer = bgzf::io::Writer::new(file.reopen().unwrap());
        writer.write_all(contents.as_bytes()).unwrap();
        writer.try_finish().unwrap();
        file
    }

    #[test]
    fn merges_sources_and_builds_samples() {
        let first = sam_file(FIRST);
        let second = bam_file(SECOND);
        let inputs = [
            AlignmentInput::new(7, first.path(), "first"),
            AlignmentInput::new(9, second.path(), "second"),
        ];
        let mut set = AlignmentSet::open(inputs, None, SampleSelection::default()).unwrap();

        assert_eq!(
            set.reference_sequences(),
            [ReferenceSequence {
                name: Box::from(&b"chr1"[..]),
                length: 20
            }]
        );
        assert_eq!(
            set.samples()
                .samples()
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            ["S1", "S2"]
        );

        let (source, record) = set.next_record().unwrap().unwrap();
        assert_eq!(
            (source, record.name(), record.alignment_start()),
            (9, b"r2".as_slice(), 1)
        );
        assert_eq!(set.samples().sample_index(source, &record), Ok(Some(1)));

        let (source, record) = set.next_record().unwrap().unwrap();
        assert_eq!(
            (source, record.name(), record.alignment_start()),
            (7, b"r1".as_slice(), 2)
        );
        assert_eq!(set.samples().sample_index(source, &record), Ok(Some(0)));
        assert!(set.next_record().unwrap().is_none());
    }

    #[test]
    fn rejects_reference_dictionary_mismatch() {
        let first = sam_file(FIRST);
        let second = sam_file(&SECOND.replace("LN:20", "LN:21"));
        let inputs = [
            AlignmentInput::new(7, first.path(), "first"),
            AlignmentInput::new(9, second.path(), "second"),
        ];
        assert!(matches!(
            AlignmentSet::open(inputs, None, SampleSelection::default()),
            Err(CallError::ReferenceDictionaryMismatch(_))
        ));
    }

    #[test]
    fn reads_raw_bam_and_bgzf_sam() {
        for file in [raw_bam_file(FIRST), bgzf_sam_file(FIRST)] {
            let mut set = AlignmentSet::open(
                [AlignmentInput::new(1, file.path(), "first")],
                None,
                SampleSelection::default(),
            )
            .unwrap();
            let (_, record) = set.next_record().unwrap().unwrap();
            assert_eq!(
                (record.name(), record.alignment_start()),
                (b"r1".as_slice(), 2)
            );
            assert!(set.next_record().unwrap().is_none());
        }
    }

    #[test]
    fn retains_every_record_in_a_cram_container() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let reference = fixtures.join("alignment-reference.fa");
        let mut set = AlignmentSet::open(
            [AlignmentInput::new(
                1,
                fixtures.join("alignment.cram"),
                "first",
            )],
            Some(&reference),
            SampleSelection::default(),
        )
        .unwrap();

        let (_, first) = set.next_record().unwrap().unwrap();
        let (_, second) = set.next_record().unwrap().unwrap();
        assert_eq!(
            (first.name(), first.alignment_start()),
            (b"r1".as_slice(), 1)
        );
        assert_eq!(
            (second.name(), second.alignment_start()),
            (b"r2".as_slice(), 2)
        );
        assert!(set.next_record().unwrap().is_none());
    }
}
