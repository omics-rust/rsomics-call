use std::fs::File;
use std::io::{Read, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

use noodles::{
    core::{Position, Region},
    fasta,
};
use rsomics_pileup::{BaqOptions, PileupEngine, PileupOptions};

use crate::{
    AlignmentInput, AlignmentSet, CallError, CalledSite, CalledVariantWriter, CalledVcfSchema,
    IndelLikelihoodConfig, IndelSiteBuilder, LikelihoodSite, LikelihoodVariantReader, Nucleotide,
    ReferenceSequence, Result, SampleMap, SampleSelection, SnpLikelihoodConfig, SnpSiteBuilder,
    VariantOutputFormat,
};

pub struct SnpLikelihoodRun {
    alignments: AlignmentSet,
    reference: ReferenceCache,
    pileup: PileupEngine,
    sites: SnpSiteBuilder,
    indels: Option<IndelSiteBuilder>,
    baq: Option<BaqRun>,
}

#[derive(Clone, Copy)]
struct BaqRun {
    mode: BaqMode,
    maximum_read_len: usize,
    options: BaqOptions,
}

#[derive(Clone, Copy)]
enum BaqMode {
    Full,
    Partial,
}

impl SnpLikelihoodRun {
    pub fn open(
        inputs: impl IntoIterator<Item = AlignmentInput>,
        reference: impl AsRef<Path>,
        selection: SampleSelection,
        pileup_options: PileupOptions,
        likelihood_config: SnpLikelihoodConfig,
    ) -> Result<Self> {
        let reference = reference.as_ref();
        let alignments = AlignmentSet::open(inputs, Some(reference), selection)?;
        let lengths = alignments
            .reference_sequences()
            .iter()
            .map(ReferenceSequence::length);
        let pileup = PileupEngine::new(lengths, pileup_options);
        let sites = SnpSiteBuilder::new(alignments.samples().samples().len(), likelihood_config)?;
        let reference = ReferenceCache::open(reference, alignments.reference_sequences())?;
        Ok(Self {
            alignments,
            reference,
            pileup,
            sites,
            indels: None,
            baq: None,
        })
    }

    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        self.alignments.reference_sequences()
    }

    pub fn samples(&self) -> &SampleMap {
        self.alignments.samples()
    }

    pub fn with_full_baq(mut self, maximum_read_len: usize, redo: bool) -> Self {
        self.set_baq(BaqMode::Full, maximum_read_len, redo);
        self
    }

    pub fn with_partial_baq(mut self, maximum_read_len: usize, redo: bool) -> Self {
        self.set_baq(BaqMode::Partial, maximum_read_len, redo);
        self
    }

    pub fn with_indels(mut self, config: IndelLikelihoodConfig) -> Result<Self> {
        self.indels = Some(IndelSiteBuilder::new(
            self.alignments.samples().samples().len(),
            config,
        )?);
        Ok(self)
    }

    fn set_baq(&mut self, mode: BaqMode, maximum_read_len: usize, redo: bool) {
        self.baq = Some(BaqRun {
            mode,
            maximum_read_len,
            options: BaqOptions {
                adjust_qualities: true,
                extended: true,
                redo,
            },
        });
    }

    pub fn run(mut self, mut emit: impl FnMut(LikelihoodSite) -> Result<()>) -> Result<()> {
        while let Some((source_id, record)) = self.alignments.next_record()? {
            self.pileup.push_with_source(source_id, record)?;
            drain_sites(
                &mut self.pileup,
                &mut self.reference,
                &mut self.sites,
                &mut self.indels,
                self.alignments.samples(),
                self.baq,
                &mut emit,
            )?;
        }
        self.pileup.finish()?;
        drain_sites(
            &mut self.pileup,
            &mut self.reference,
            &mut self.sites,
            &mut self.indels,
            self.alignments.samples(),
            self.baq,
            &mut emit,
        )
    }

    pub fn run_called(
        self,
        mut call: impl FnMut(&LikelihoodSite) -> Result<CalledSite>,
        mut emit: impl FnMut(CalledSite) -> Result<()>,
    ) -> Result<()> {
        self.run(|site| emit(call(&site)?))
    }
}

pub fn run_likelihood_calls<R, W>(
    mut reader: LikelihoodVariantReader<R>,
    writer: W,
    output_schema: CalledVcfSchema,
    output_format: VariantOutputFormat,
    mut call: impl FnMut(&LikelihoodSite) -> Result<CalledSite>,
) -> Result<W>
where
    R: Read,
    W: Write,
{
    let mut writer = CalledVariantWriter::new(writer, output_schema, output_format)?;
    while let Some(site) = reader.read_site()? {
        let record = reader.record_number();
        let called = call(&site).map_err(|error| CallError::LikelihoodCallRecord {
            record,
            source: Box::new(error),
        })?;
        writer.write_site(&called)?;
    }
    writer.finish()
}

fn drain_sites(
    pileup: &mut PileupEngine,
    reference: &mut ReferenceCache,
    sites: &mut SnpSiteBuilder,
    indels: &mut Option<IndelSiteBuilder>,
    samples: &SampleMap,
    baq: Option<BaqRun>,
    emit: &mut impl FnMut(LikelihoodSite) -> Result<()>,
) -> Result<()> {
    pileup.drain_with(|context| {
        if let Some(baq) = baq {
            let mut fetch_reference = |reference_id, range, buffer: &mut Vec<u8>| {
                buffer.extend_from_slice(reference.sequence(reference_id, range)?);
                Ok::<_, CallError>(())
            };
            match baq.mode {
                BaqMode::Full => context.apply_full_baq(
                    baq.maximum_read_len,
                    baq.options,
                    &mut fetch_reference,
                )?,
                BaqMode::Partial => context.apply_partial_baq(
                    baq.maximum_read_len,
                    baq.options,
                    &mut fetch_reference,
                )?,
            }
        }
        let column = context.column();
        let reference_base = reference.base(column.reference_id(), column.position())?;
        let site = sites.build(&column, reference_base, |source_id, record| {
            samples.sample_index(source_id, record)
        })?;
        emit(site)?;
        if let Some(indels) = indels {
            let reference_length = reference.length(column.reference_id())?;
            let mut fetch_reference = |range, buffer: &mut Vec<u8>| {
                buffer.extend_from_slice(reference.sequence(column.reference_id(), range)?);
                Ok::<_, CallError>(())
            };
            if let Some(site) = indels.build(
                &column,
                reference_length,
                |source_id, record| samples.sample_index(source_id, record),
                &mut fetch_reference,
            )? {
                emit(site)?;
            }
        }
        Ok(())
    })
}

struct ReferenceCache {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    references: Box<[(Box<[u8]>, u64)]>,
    reference_id: Option<usize>,
    sequence_start: usize,
    sequence: Vec<u8>,
}

impl ReferenceCache {
    fn open(path: &Path, references: &[ReferenceSequence]) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| reference_error(path, error))?;
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            references: references
                .iter()
                .map(|reference| {
                    (
                        reference.name().to_vec().into_boxed_slice(),
                        reference.length(),
                    )
                })
                .collect(),
            reference_id: None,
            sequence_start: 0,
            sequence: Vec::new(),
        })
    }

    fn base(&mut self, reference_id: i32, position: i64) -> Result<Nucleotide> {
        let position =
            usize::try_from(position).map_err(|error| reference_error(&self.path, error))?;
        let end = position
            .checked_add(1)
            .ok_or_else(|| reference_error(&self.path, "reference position overflows"))?;
        let base = self.sequence(reference_id, position..end)?[0];
        Ok(match base.to_ascii_uppercase() {
            b'A' => Nucleotide::A,
            b'C' => Nucleotide::C,
            b'G' => Nucleotide::G,
            b'T' => Nucleotide::T,
            _ => Nucleotide::N,
        })
    }

    fn length(&self, reference_id: i32) -> Result<usize> {
        let reference_id =
            usize::try_from(reference_id).map_err(|error| reference_error(&self.path, error))?;
        let length = self
            .references
            .get(reference_id)
            .map(|(_, length)| *length)
            .ok_or_else(|| reference_error(&self.path, "reference ID is absent"))?;
        usize::try_from(length).map_err(|error| reference_error(&self.path, error))
    }

    fn sequence(&mut self, reference_id: i32, range: Range<usize>) -> Result<&[u8]> {
        const CHUNK_SIZE: usize = 1024 * 1024;

        let reference_id =
            usize::try_from(reference_id).map_err(|error| reference_error(&self.path, error))?;
        let (name, length) = self
            .references
            .get(reference_id)
            .map(|(name, length)| (name.to_vec(), *length))
            .ok_or_else(|| reference_error(&self.path, "reference ID is absent"))?;
        let length = usize::try_from(length).map_err(|error| reference_error(&self.path, error))?;
        if range.start >= range.end || range.end > length {
            return Err(reference_error(
                &self.path,
                "requested reference range is invalid",
            ));
        }
        if self.reference_id != Some(reference_id)
            || range.start < self.sequence_start
            || range.end > self.sequence_start + self.sequence.len()
        {
            let start = range.start / CHUNK_SIZE * CHUNK_SIZE;
            let chunk_end = start
                .checked_add(CHUNK_SIZE)
                .map_or(length, |end| end.min(length));
            let end = chunk_end.max(range.end);
            let interval_start = Position::try_from(start + 1)
                .map_err(|error| reference_error(&self.path, error))?;
            let interval_end =
                Position::try_from(end).map_err(|error| reference_error(&self.path, error))?;
            let record = self
                .reader
                .query(&Region::new(name, interval_start..=interval_end))
                .map_err(|error| reference_error(&self.path, error))?;
            self.sequence.clear();
            self.sequence.extend_from_slice(record.sequence().as_ref());
            self.reference_id = Some(reference_id);
            self.sequence_start = start;
        }
        let start = range
            .start
            .checked_sub(self.sequence_start)
            .ok_or_else(|| reference_error(&self.path, "invalid reference cache position"))?;
        let end = range
            .end
            .checked_sub(self.sequence_start)
            .ok_or_else(|| reference_error(&self.path, "invalid reference cache position"))?;
        self.sequence
            .get(start..end)
            .ok_or_else(|| reference_error(&self.path, "reference range is outside the cache"))
    }
}

fn reference_error(path: &Path, error: impl std::fmt::Display) -> CallError {
    CallError::ReferenceInput {
        path: path.display().to_string(),
        message: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;

    use noodles::vcf;
    use noodles_util::variant;
    use tempfile::NamedTempFile;

    use super::*;
    use crate::{
        Allele, ConsensusCaller, IndelSummary, LikelihoodVariantWriter, LikelihoodVcfSchema,
        Ploidy, SampleEvidence, SampleLikelihood,
    };

    fn sam_file(sample: &str, base: char) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:20\n\
             @RG\tID:rg\tSM:{sample}\n\
             read\t0\tchr1\t1\t60\t1M\t*\t0\t0\t{base}\tI\tRG:Z:rg"
        )
        .unwrap();
        file
    }

    fn repeated_sam_file(count: usize) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:20\n\
             @RG\tID:rg\tSM:S1"
        )
        .unwrap();
        for index in 0..count {
            writeln!(
                file,
                "read{index}\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tRG:Z:rg"
            )
            .unwrap();
        }
        file
    }

    fn indel_sam_file(sample: &str, cigar: &str, sequence: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:MX\tLN:11\n\
             @RG\tID:rg\tSM:{sample}"
        )
        .unwrap();
        let qualities = "I".repeat(sequence.len());
        for index in 1..=2 {
            writeln!(
                file,
                "r{index}\t0\tMX\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{qualities}\tRG:Z:rg"
            )
            .unwrap();
        }
        file
    }

    #[test]
    fn streams_cram_records_into_likelihood_sites() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let reference = fixtures.join("alignment-reference.fa");
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(
                3,
                fixtures.join("alignment.cram"),
                "input",
            )],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap();
        assert_eq!(
            run.samples()
                .samples()
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            ["S1"]
        );

        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 2);
        assert_eq!(
            (sites[0].position(), sites[0].reference().as_bytes()),
            (1, b"C".as_slice())
        );
        assert_eq!(
            (sites[1].position(), sites[1].reference().as_bytes()),
            (2, b"G".as_slice())
        );
        for site in sites {
            assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[0, 3, 40][..]));
            assert_eq!(site.samples()[0].evidence().depth(), 1);
        }
    }

    #[test]
    fn matches_bcftools_1_24_multisample_snp_likelihoods() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let first = sam_file("S1", 'A');
        let second = sam_file("S2", 'G');
        let run = SnpLikelihoodRun::open(
            [
                AlignmentInput::new(1, first.path(), "first"),
                AlignmentInput::new(2, second.path(), "second"),
            ],
            fixtures.join("alignment-reference.fa"),
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap();

        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 1);
        let site = &sites[0];
        assert_eq!(site.position(), 0);
        assert_eq!(site.reference().as_bytes(), b"A");
        assert_eq!(
            site.alternates()
                .iter()
                .map(|allele| allele.as_bytes())
                .collect::<Vec<_>>(),
            [b"G".as_slice(), b"<*>".as_slice()]
        );
        assert_eq!(
            site.samples()[0].phred_likelihoods(),
            Some(&[0, 3, 40, 3, 40, 40][..])
        );
        assert_eq!(
            site.samples()[1].phred_likelihoods(),
            Some(&[40, 3, 0, 40, 3, 40][..])
        );
    }

    #[test]
    fn fused_call_matches_materialized_typed_pipeline() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let reference = fixtures.join("alignment-reference.fa");
        let first = sam_file("S1", 'A');
        let second = sam_file("S2", 'G');
        let open = || {
            SnpLikelihoodRun::open(
                [
                    AlignmentInput::new(1, first.path(), "first"),
                    AlignmentInput::new(2, second.path(), "second"),
                ],
                &reference,
                SampleSelection::default(),
                PileupOptions::default(),
                SnpLikelihoodConfig::default(),
            )
            .unwrap()
        };
        let caller = crate::MultiallelicCaller::default();
        let mut materialized = Vec::new();
        open()
            .run(|site| {
                materialized.push(site);
                Ok(())
            })
            .unwrap();
        let expected = materialized
            .iter()
            .map(|site| caller.call(site).unwrap())
            .collect::<Vec<_>>();
        let mut fused = Vec::new();

        open()
            .run_called(
                |site| caller.call(site),
                |site| {
                    fused.push(site);
                    Ok(())
                },
            )
            .unwrap();

        assert_eq!(fused, expected);
    }

    #[test]
    fn matches_bcftools_1_24_per_input_depth_limit() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let input = repeated_sam_file(5);
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, input.path(), "input")],
            fixtures.join("alignment-reference.fa"),
            SampleSelection::default(),
            PileupOptions {
                maximum_depth_per_source: Some(2),
                ..PileupOptions::default()
            },
            SnpLikelihoodConfig::default(),
        )
        .unwrap();

        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].samples()[0].evidence().depth(), 2);
        assert_eq!(
            sites[0].samples()[0].phred_likelihoods(),
            Some(&[0, 6, 73][..])
        );
    }

    #[test]
    fn full_baq_matches_bcftools_1_24_likelihoods() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        let alignment = directory.path().join("alignment.sam");
        fs::write(&reference, b">MX\nCGTCTACTACG\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"MX\t11\t4\t11\t12\n").unwrap();
        fs::write(
            &alignment,
            b"@HD\tVN:1.6\tSO:coordinate\n\
              @SQ\tSN:MX\tLN:11\n\
              @RG\tID:rg\tSM:sample\n\
              M\t64\tMX\t1\t60\t11M\t*\t0\t0\tCGTCTCCTACG\tIIIIIIIIIII\tRG:Z:rg\n\
              X\t64\tMX\t1\t60\t5=1X5=\t*\t0\t0\tCGTCTCCTACG\tIIIIIIIIIII\tRG:Z:rg\n",
        )
        .unwrap();
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, alignment, "input")],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_full_baq(500, false);
        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 11);
        assert_eq!(
            sites[0].samples()[0].phred_likelihoods(),
            Some(&[0, 6, 66][..])
        );
        assert_eq!(
            sites[5].samples()[0].phred_likelihoods(),
            Some(&[73, 6, 0, 73, 6, 73][..])
        );
        assert_eq!(
            sites[10].samples()[0].phred_likelihoods(),
            Some(&[0, 6, 66][..])
        );
    }

    #[test]
    fn full_baq_and_overlap_order_matches_bcftools_1_24() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        let alignment = directory.path().join("alignment.sam");
        fs::write(&reference, b">chr1\nAAAAAAA\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"chr1\t7\t6\t7\t8\n").unwrap();
        fs::write(
            &alignment,
            b"@HD\tVN:1.6\tSO:coordinate\n\
              @SQ\tSN:chr1\tLN:7\n\
              @RG\tID:rg\tSM:sample\n\
              pair\t99\tchr1\t1\t60\t5M\t=\t3\t7\tAAAAA\tIIIII\tRG:Z:rg\n\
              pair\t147\tchr1\t3\t60\t5M\t=\t1\t-7\tAAAAA\tIIIII\tRG:Z:rg\n",
        )
        .unwrap();
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, alignment, "input")],
            reference,
            SampleSelection::default(),
            PileupOptions {
                adjust_overlaps: true,
                ..PileupOptions::default()
            },
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_full_baq(500, false);
        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 7);
        for index in [0, 1, 6] {
            assert_eq!(
                sites[index].samples()[0].phred_likelihoods(),
                Some(&[0, 3, 4][..])
            );
        }
        for site in sites.iter().take(6).skip(2) {
            assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[0, 0, 0][..]));
        }
        assert_eq!(
            sites
                .iter()
                .map(|site| site.samples()[0].evidence().depth())
                .collect::<Vec<_>>(),
            [1, 1, 2, 2, 2, 1, 1]
        );
    }

    #[test]
    fn partial_baq_matches_bcftools_1_24_indel_trigger() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        let alignment = directory.path().join("alignment.sam");
        fs::write(&reference, b">MX\nCGTCTACTACG\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"MX\t11\t4\t11\t12\n").unwrap();
        fs::write(
            &alignment,
            b"@HD\tVN:1.6\tSO:coordinate\n\
              @SQ\tSN:MX\tLN:11\n\
              @RG\tID:rg\tSM:sample\n\
              r1\t0\tMX\t1\t60\t5M1I6M\t*\t0\t0\tCGTCTCACTACG\tIIIIIIIIIIII\tRG:Z:rg\n\
              r2\t0\tMX\t1\t60\t5M1I6M\t*\t0\t0\tCGTCTCACTACG\tIIIIIIIIIIII\tRG:Z:rg\n",
        )
        .unwrap();
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, alignment, "input")],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_partial_baq(500, false);
        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 11);
        for site in [&sites[0], &sites[10]] {
            assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[0, 6, 66][..]));
        }
        for site in sites.iter().take(10).skip(1) {
            assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[0, 6, 73][..]));
        }
        assert!(
            sites
                .iter()
                .all(|site| site.samples()[0].evidence().depth() == 2)
        );
    }

    #[test]
    fn indel_likelihoods_match_bcftools_1_24() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        fs::write(&reference, b">MX\nCGTCTACTACG\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"MX\t11\t4\t11\t12\n").unwrap();
        let alignment = indel_sam_file("sample", "5M1I6M", "CGTCTCACTACG");
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, alignment.path(), "input")],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_partial_baq(500, false)
        .with_indels(IndelLikelihoodConfig::default())
        .unwrap();
        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        let site = sites
            .iter()
            .find(|site| site.alternates()[0].as_bytes() == b"TC")
            .unwrap();
        assert_eq!(site.position(), 4);
        assert_eq!(site.reference().as_bytes(), b"T");
        assert_eq!(site.alternates()[0].as_bytes(), b"TC");
        assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[56, 6, 0][..]));
        assert_eq!(site.samples()[0].evidence().depth(), 2);
        assert_eq!(site.samples()[0].evidence().allele_depths(), &[0, 2]);
        assert_eq!(site.samples()[0].evidence().allele_quality_sums(), &[0, 62]);
        assert_eq!(
            site.indel_summary(),
            Some(IndelSummary::new(2, 1.0).unwrap())
        );
    }

    #[test]
    fn deletion_likelihoods_match_bcftools_1_24() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        fs::write(&reference, b">MX\nCGTCTACTACG\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"MX\t11\t4\t11\t12\n").unwrap();
        let alignment = indel_sam_file("sample", "5M1D5M", "CGTCTCTACG");
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, alignment.path(), "input")],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_partial_baq(500, false)
        .with_indels(IndelLikelihoodConfig::default())
        .unwrap();
        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        let site = sites
            .iter()
            .find(|site| site.reference().as_bytes() == b"TA")
            .unwrap();
        assert_eq!(site.position(), 4);
        assert_eq!(site.alternates()[0].as_bytes(), b"T");
        assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[44, 6, 0][..]));
        assert_eq!(site.samples()[0].evidence().allele_depths(), &[0, 2]);
        assert_eq!(site.samples()[0].evidence().allele_quality_sums(), &[0, 48]);
    }

    #[test]
    fn multisample_indel_likelihoods_match_bcftools_1_24() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        fs::write(&reference, b">MX\nCGTCTACTACG\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"MX\t11\t4\t11\t12\n").unwrap();
        let alternate = indel_sam_file("alternate", "5M1I6M", "CGTCTCACTACG");
        let reference_reads = indel_sam_file("reference", "11M", "CGTCTACTACG");
        let trailing = indel_sam_file("trailing", "2M", "CG");
        let run = SnpLikelihoodRun::open(
            [
                AlignmentInput::new(1, alternate.path(), "alternate"),
                AlignmentInput::new(2, reference_reads.path(), "reference"),
                AlignmentInput::new(3, trailing.path(), "trailing"),
            ],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_partial_baq(500, false)
        .with_indels(IndelLikelihoodConfig::default())
        .unwrap();
        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        let site = sites
            .iter()
            .find(|site| site.alternates()[0].as_bytes() == b"TC")
            .unwrap();
        assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[56, 6, 0][..]));
        assert_eq!(site.samples()[1].phred_likelihoods(), Some(&[0, 6, 47][..]));
        assert_eq!(site.samples()[2].phred_likelihoods(), Some(&[0, 0, 0][..]));
        assert_eq!(site.samples()[0].evidence().allele_quality_sums(), &[0, 62]);
        assert_eq!(site.samples()[1].evidence().allele_quality_sums(), &[52, 0]);
        assert_eq!(site.samples()[2].evidence().depth(), 0);
    }

    #[test]
    fn streams_likelihood_calls_across_all_variant_encodings() {
        let formats = [
            VariantOutputFormat::Vcf,
            VariantOutputFormat::VcfBgzf,
            VariantOutputFormat::BcfRaw,
            VariantOutputFormat::BcfBgzf,
        ];
        let sites = [
            likelihood_site(0, Some([0, 3, 40])),
            likelihood_site(1, Some([220, 99, 0])),
        ];

        for format in formats {
            let schema = LikelihoodVcfSchema::new([(b"chr1".as_slice(), 20)], ["sample"]).unwrap();
            let mut input =
                LikelihoodVariantWriter::new(Vec::new(), schema.clone(), format).unwrap();
            for site in &sites {
                input.write_site(site).unwrap();
            }
            let input = input.finish().unwrap();
            let reader = LikelihoodVariantReader::new(&input[..]).unwrap();
            let output_schema = CalledVcfSchema::from_consensus_likelihood(reader.schema());
            let output = run_likelihood_calls(reader, Vec::new(), output_schema, format, |site| {
                ConsensusCaller::default().call(site)
            })
            .unwrap();
            let mut reader = variant::io::Reader::new(&output[..]).unwrap();
            let header = reader.read_header().unwrap();
            assert!(!header.formats().contains_key("GP"));
            let mut record = variant::Record::default();

            for position in [1, 2] {
                assert_ne!(reader.read_record(&mut record).unwrap(), 0);
                let record =
                    vcf::variant::RecordBuf::try_from_variant_record(&header, &record).unwrap();
                assert_eq!(record.variant_start().map(usize::from), Some(position));
            }
            assert_eq!(reader.read_record(&mut record).unwrap(), 0);
        }
    }

    #[test]
    fn reports_the_record_that_cannot_be_called() {
        let schema = LikelihoodVcfSchema::new([(b"chr1".as_slice(), 20)], ["sample"]).unwrap();
        let mut input =
            LikelihoodVariantWriter::new(Vec::new(), schema.clone(), VariantOutputFormat::Vcf)
                .unwrap();
        input.write_site(&likelihood_site(0, None)).unwrap();
        let input = input.finish().unwrap();
        let reader = LikelihoodVariantReader::new(&input[..]).unwrap();
        let output_schema = CalledVcfSchema::from_consensus_likelihood(reader.schema());

        assert!(matches!(
            run_likelihood_calls(
                reader,
                Vec::new(),
                output_schema,
                VariantOutputFormat::Vcf,
                |site| ConsensusCaller::default().call(site),
            ),
            Err(CallError::LikelihoodCallRecord { record: 1, .. })
        ));
    }

    fn likelihood_site(position: u64, likelihoods: Option<[u32; 3]>) -> LikelihoodSite {
        let evidence = SampleEvidence::new(1, [1, 0], [40, 0]).unwrap();
        let sample = match likelihoods {
            Some(values) => {
                SampleLikelihood::observed(Ploidy::new(2).unwrap(), values, evidence).unwrap()
            }
            None => SampleLikelihood::new(Ploidy::new(2).unwrap(), None, evidence).unwrap(),
        };
        LikelihoodSite::new(
            0,
            position,
            Allele::new(&b"A"[..]).unwrap(),
            [Allele::new(&b"G"[..]).unwrap()],
            [1.0, 1.0],
            [sample],
        )
        .unwrap()
    }
}
