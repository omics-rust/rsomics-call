use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use noodles::{
    bcf, bgzf,
    core::{Position, Region},
    vcf,
    vcf::variant::{Record as _, io::Write as _},
};
use noodles_util::variant;

use crate::{CallError, CalledSite, CalledVcfSchema, LikelihoodSite, LikelihoodVcfSchema, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantOutputFormat {
    Vcf,
    VcfBgzf,
    BcfRaw,
    BcfBgzf,
}

pub struct LikelihoodVariantReader<R>
where
    R: Read,
{
    inner: variant::io::Reader<R>,
    projection: LikelihoodProjection,
    record: variant::Record,
    record_number: u64,
    started: bool,
}

impl<R> LikelihoodVariantReader<R>
where
    R: Read,
{
    pub fn new(reader: R) -> Result<Self> {
        let mut inner = variant::io::Reader::new(reader).map_err(input_error)?;
        let header = inner.read_header().map_err(input_error)?;
        Ok(Self {
            inner,
            projection: LikelihoodProjection::new(header)?,
            record: variant::Record::default(),
            record_number: 0,
            started: false,
        })
    }

    pub fn schema(&self) -> &LikelihoodVcfSchema {
        self.projection.schema()
    }

    pub fn record_number(&self) -> u64 {
        self.record_number
    }

    pub fn select_samples(
        mut self,
        samples: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<Self> {
        self.require_unread()?;
        self.projection.select_samples(samples)?;
        Ok(self)
    }

    pub fn exclude_samples(
        mut self,
        samples: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<Self> {
        self.require_unread()?;
        self.projection.exclude_samples(samples)?;
        Ok(self)
    }

    pub fn with_prior_frequencies(
        mut self,
        total_alleles_tag: impl Into<String>,
        alternate_counts_tag: impl Into<String>,
    ) -> Result<Self> {
        self.require_unread()?;
        self.projection
            .set_prior_frequency_tags(total_alleles_tag, alternate_counts_tag)?;
        Ok(self)
    }

    pub fn read_site(&mut self) -> Result<Option<LikelihoodSite>> {
        self.started = true;
        let record_number = self.record_number + 1;
        let size = self
            .inner
            .read_record(&mut self.record)
            .map_err(|error| record_error(record_number, error))?;
        if size == 0 {
            return Ok(None);
        }
        self.record_number = record_number;
        let record = vcf::variant::RecordBuf::try_from_variant_record(
            self.projection.input_header(),
            &self.record,
        )
        .map_err(|error| record_error(record_number, error))?;
        self.projection
            .decode(&record)
            .map_err(|error| CallError::LikelihoodVariantRecord {
                record: record_number,
                message: error.to_string(),
            })
            .map(Some)
    }

    fn require_unread(&self) -> Result<()> {
        if self.started {
            return Err(CallError::LateLikelihoodReaderConfiguration);
        }
        Ok(())
    }
}

pub struct IndexedLikelihoodVariantReader {
    inner: variant::io::indexed_reader::IndexedReader<bgzf::io::Reader<File>>,
    projection: LikelihoodProjection,
}

impl IndexedLikelihoodVariantReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let mut inner = variant::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| path_input_error(path, error))?;
        let header = inner
            .read_header()
            .map_err(|error| path_input_error(path, error))?;
        Ok(Self {
            inner,
            projection: LikelihoodProjection::new(header)?,
        })
    }

    pub fn schema(&self) -> &LikelihoodVcfSchema {
        self.projection.schema()
    }

    pub fn select_samples(
        mut self,
        samples: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<Self> {
        self.projection.select_samples(samples)?;
        Ok(self)
    }

    pub fn exclude_samples(
        mut self,
        samples: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<Self> {
        self.projection.exclude_samples(samples)?;
        Ok(self)
    }

    pub fn with_prior_frequencies(
        mut self,
        total_alleles_tag: impl Into<String>,
        alternate_counts_tag: impl Into<String>,
    ) -> Result<Self> {
        self.projection
            .set_prior_frequency_tags(total_alleles_tag, alternate_counts_tag)?;
        Ok(self)
    }

    pub(crate) fn normalize_regions(
        &self,
        regions: impl IntoIterator<Item = Region>,
    ) -> Result<PreparedVariantRegions> {
        normalize_variant_regions(self.projection.input_header(), regions)
    }

    pub(crate) fn normalize_region_file(
        &self,
        regions: impl IntoIterator<Item = Region>,
    ) -> Result<PreparedVariantRegions> {
        normalize_variant_region_file(self.projection.input_header(), regions)
    }

    pub(crate) fn visit_normalized_regions(
        mut self,
        regions: PreparedVariantRegions,
        mut visit: impl FnMut(u64, LikelihoodSite) -> Result<()>,
    ) -> Result<()> {
        let header = self.projection.input_header();
        let mut record_number = 0;

        for (region_index, region) in regions.values.iter().enumerate() {
            let records = self
                .inner
                .query(header, &region.query)
                .map_err(|error| query_error(&region.query, error))?;
            for result in records {
                record_number += 1;
                let record = result.map_err(|error| record_error(record_number, error))?;
                let record =
                    vcf::variant::RecordBuf::try_from_variant_record(header, record.as_ref())
                        .map_err(|error| record_error(record_number, error))?;
                let start = record.variant_start().ok_or_else(|| {
                    record_error(
                        record_number,
                        io::Error::new(io::ErrorKind::InvalidData, "missing variant position"),
                    )
                })?;
                let end = record
                    .variant_end(header)
                    .map_err(|error| record_error(record_number, error))?;
                if regions.deduplicate
                    && regions.values[..region_index]
                        .iter()
                        .any(|previous| previous.overlaps(region.reference_id, start, end))
                {
                    continue;
                }
                let site = self.projection.decode(&record).map_err(|error| {
                    CallError::LikelihoodVariantRecord {
                        record: record_number,
                        message: error.to_string(),
                    }
                })?;
                visit(record_number, site)?;
            }
        }
        Ok(())
    }
}

pub(crate) struct VariantRegion {
    query: Region,
    reference_id: usize,
    start: Position,
    end: Position,
}

pub(crate) struct PreparedVariantRegions {
    values: Vec<VariantRegion>,
    deduplicate: bool,
}

impl VariantRegion {
    fn overlaps(&self, reference_id: usize, start: Position, end: Position) -> bool {
        self.reference_id == reference_id && start <= self.end && end >= self.start
    }
}

fn normalize_variant_regions(
    header: &vcf::Header,
    regions: impl IntoIterator<Item = Region>,
) -> Result<PreparedVariantRegions> {
    let mut values = resolve_variant_regions(header, regions)?;
    values.sort_unstable_by_key(|region| (region.reference_id, region.start, region.end));

    let mut merged: Vec<VariantRegion> = Vec::with_capacity(values.len());
    for region in values {
        if let Some(previous) = merged.last_mut()
            && previous.reference_id == region.reference_id
            && region.start <= previous.end.checked_add(1).unwrap_or(Position::MAX)
        {
            previous.end = previous.end.max(region.end);
            previous.query = Region::new(
                previous.query.name().to_vec(),
                previous.start..=previous.end,
            );
        } else {
            merged.push(region);
        }
    }
    Ok(PreparedVariantRegions {
        values: merged,
        deduplicate: true,
    })
}

fn normalize_variant_region_file(
    header: &vcf::Header,
    regions: impl IntoIterator<Item = Region>,
) -> Result<PreparedVariantRegions> {
    let mut values = resolve_variant_regions(header, regions)?;
    let mut ranks = vec![usize::MAX; header.contigs().len()];
    let mut next_rank = 0;
    for region in &values {
        if ranks[region.reference_id] == usize::MAX {
            ranks[region.reference_id] = next_rank;
            next_rank += 1;
        }
    }
    values.sort_by_key(|region| (ranks[region.reference_id], region.start, region.end));
    Ok(PreparedVariantRegions {
        values,
        deduplicate: true,
    })
}

fn resolve_variant_regions(
    header: &vcf::Header,
    regions: impl IntoIterator<Item = Region>,
) -> Result<Vec<VariantRegion>> {
    let mut values = Vec::new();
    for region in regions {
        let name =
            std::str::from_utf8(region.name().as_ref()).map_err(|_| CallError::InvalidRegion {
                region: region.to_string(),
                message: "reference sequence name is not UTF-8".to_owned(),
            })?;
        let (reference_id, _, contig) =
            header
                .contigs()
                .get_full(name)
                .ok_or_else(|| CallError::InvalidRegion {
                    region: region.to_string(),
                    message: "reference sequence is absent from the variant header".to_owned(),
                })?;
        let start = region.interval().start().unwrap_or(Position::MIN);
        let end = match contig.length().and_then(Position::new) {
            Some(length) if start > length => {
                return Err(CallError::InvalidRegion {
                    region: region.to_string(),
                    message: "interval is outside the reference sequence".to_owned(),
                });
            }
            Some(length) => region.interval().end().unwrap_or(length).min(length),
            None => region.interval().end().unwrap_or(Position::MAX),
        };
        values.push(VariantRegion {
            query: Region::new(name, start..=end),
            reference_id,
            start,
            end,
        });
    }
    if values.is_empty() {
        return Err(CallError::MissingRegions);
    }
    Ok(values)
}

pub(crate) struct LikelihoodProjection {
    input_schema: LikelihoodVcfSchema,
    schema: LikelihoodVcfSchema,
    sample_indices: Box<[usize]>,
}

impl LikelihoodProjection {
    pub(crate) fn new(header: vcf::Header) -> Result<Self> {
        let schema = LikelihoodVcfSchema::from_header(header)?;
        let sample_indices = (0..schema.header().sample_names().len()).collect();
        Ok(Self {
            input_schema: schema.clone(),
            schema,
            sample_indices,
        })
    }

    pub(crate) fn schema(&self) -> &LikelihoodVcfSchema {
        &self.schema
    }

    pub(crate) fn input_header(&self) -> &vcf::Header {
        self.input_schema.header()
    }

    pub(crate) fn select_samples(
        &mut self,
        samples: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<()> {
        let mut seen = HashSet::new();
        let mut indices = Vec::new();
        for sample in samples {
            let sample = sample.as_ref();
            if !seen.insert(sample.to_owned()) {
                return Err(CallError::DuplicateSampleSelection(sample.to_owned()));
            }
            let index = self
                .input_schema
                .header()
                .sample_names()
                .get_index_of(sample)
                .ok_or_else(|| CallError::MissingSelectedSample(sample.to_owned()))?;
            indices.push(index);
        }
        self.set_sample_indices(indices)
    }

    pub(crate) fn exclude_samples(
        &mut self,
        samples: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<()> {
        let mut excluded = HashSet::new();
        for sample in samples {
            let sample = sample.as_ref();
            if !excluded.insert(sample.to_owned()) {
                return Err(CallError::DuplicateSampleSelection(sample.to_owned()));
            }
            if !self.input_schema.header().sample_names().contains(sample) {
                return Err(CallError::MissingSelectedSample(sample.to_owned()));
            }
        }
        let indices = self
            .input_schema
            .header()
            .sample_names()
            .iter()
            .enumerate()
            .filter_map(|(index, sample)| (!excluded.contains(sample)).then_some(index))
            .collect();
        self.set_sample_indices(indices)
    }

    pub(crate) fn decode(&self, record: &vcf::variant::RecordBuf) -> Result<LikelihoodSite> {
        self.input_schema
            .decode_selected_likelihood(record, &self.sample_indices)
    }

    pub(crate) fn set_prior_frequency_tags(
        &mut self,
        total: impl Into<String>,
        alternates: impl Into<String>,
    ) -> Result<()> {
        let total = total.into();
        let alternates = alternates.into();
        self.input_schema
            .set_prior_frequency_tags(total.clone(), alternates.clone())?;
        self.schema.set_prior_frequency_tags(total, alternates)
    }

    fn set_sample_indices(&mut self, indices: Vec<usize>) -> Result<()> {
        if indices.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        let mut header = self.input_schema.header().clone();
        header.sample_names_mut().clear();
        for &index in &indices {
            header
                .sample_names_mut()
                .insert(self.input_schema.header().sample_names()[index].clone());
        }
        let prior_frequency_tags = self
            .schema
            .prior_frequency_tags()
            .map(|(total, alternates)| (total.to_owned(), alternates.to_owned()));
        self.schema = LikelihoodVcfSchema::from_header(header)?;
        if let Some((total, alternates)) = prior_frequency_tags {
            self.schema.set_prior_frequency_tags(total, alternates)?;
        }
        self.sample_indices = indices.into_boxed_slice();
        Ok(())
    }
}

pub struct LikelihoodVariantWriter<W>
where
    W: Write,
{
    inner: VariantWriter<W>,
    schema: LikelihoodVcfSchema,
}

impl<W> LikelihoodVariantWriter<W>
where
    W: Write,
{
    pub fn new(
        writer: W,
        schema: LikelihoodVcfSchema,
        format: VariantOutputFormat,
    ) -> Result<Self> {
        let mut inner = VariantWriter::new(writer, format);
        inner.write_header(schema.header())?;
        Ok(Self { inner, schema })
    }

    pub fn schema(&self) -> &LikelihoodVcfSchema {
        &self.schema
    }

    pub fn write_site(&mut self, site: &LikelihoodSite) -> Result<()> {
        let record = self.schema.encode_likelihood(site)?;
        self.inner.write_record(self.schema.header(), &record)
    }

    pub fn finish(self) -> Result<W> {
        self.inner.finish()
    }
}

pub struct CalledVariantWriter<W>
where
    W: Write,
{
    inner: VariantWriter<W>,
    schema: CalledVcfSchema,
}

impl<W> CalledVariantWriter<W>
where
    W: Write,
{
    pub fn new(writer: W, schema: CalledVcfSchema, format: VariantOutputFormat) -> Result<Self> {
        let mut inner = VariantWriter::new(writer, format);
        inner.write_header(schema.header())?;
        Ok(Self { inner, schema })
    }

    pub fn schema(&self) -> &CalledVcfSchema {
        &self.schema
    }

    pub fn write_site(&mut self, site: &CalledSite) -> Result<()> {
        let record = self.schema.encode_call(site)?;
        self.inner.write_record(self.schema.header(), &record)
    }

    pub fn finish(self) -> Result<W> {
        self.inner.finish()
    }
}

enum VariantWriter<W>
where
    W: Write,
{
    Vcf(vcf::io::Writer<BufWriter<W>>),
    VcfBgzf(vcf::io::Writer<bgzf::io::Writer<W>>),
    BcfRaw(bcf::io::Writer<BufWriter<W>>),
    BcfBgzf(bcf::io::Writer<bgzf::io::Writer<W>>),
}

impl<W> VariantWriter<W>
where
    W: Write,
{
    fn new(writer: W, format: VariantOutputFormat) -> Self {
        match format {
            VariantOutputFormat::Vcf => Self::Vcf(vcf::io::Writer::new(BufWriter::new(writer))),
            VariantOutputFormat::VcfBgzf => {
                Self::VcfBgzf(vcf::io::Writer::new(bgzf::io::Writer::new(writer)))
            }
            VariantOutputFormat::BcfRaw => {
                Self::BcfRaw(bcf::io::Writer::from(BufWriter::new(writer)))
            }
            VariantOutputFormat::BcfBgzf => Self::BcfBgzf(bcf::io::Writer::new(writer)),
        }
    }

    fn write_header(&mut self, header: &vcf::Header) -> Result<()> {
        match self {
            Self::Vcf(writer) => writer.write_header(header),
            Self::VcfBgzf(writer) => writer.write_header(header),
            Self::BcfRaw(writer) => writer.write_header(header),
            Self::BcfBgzf(writer) => writer.write_header(header),
        }
        .map_err(output_error)
    }

    fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> Result<()> {
        match self {
            Self::Vcf(writer) => writer.write_variant_record(header, record),
            Self::VcfBgzf(writer) => writer.write_variant_record(header, record),
            Self::BcfRaw(writer) => writer.write_variant_record(header, record),
            Self::BcfBgzf(writer) => writer.write_variant_record(header, record),
        }
        .map_err(output_error)
    }

    fn finish(self) -> Result<W> {
        match self {
            Self::Vcf(writer) => finish_buffer(writer.into_inner()),
            Self::VcfBgzf(writer) => writer.into_inner().finish().map_err(output_error),
            Self::BcfRaw(writer) => finish_buffer(writer.into_inner()),
            Self::BcfBgzf(writer) => writer.into_inner().finish().map_err(output_error),
        }
    }
}

fn finish_buffer<W>(writer: BufWriter<W>) -> Result<W>
where
    W: Write,
{
    writer
        .into_inner()
        .map_err(|error| output_error(error.into_error()))
}

fn input_error(error: io::Error) -> CallError {
    CallError::LikelihoodVariantInput(error.to_string())
}

fn path_input_error(path: &Path, error: io::Error) -> CallError {
    CallError::LikelihoodVariantInput(format!("{}: {error}", path.display()))
}

fn query_error(region: &Region, error: io::Error) -> CallError {
    CallError::LikelihoodVariantInput(format!("querying region {region}: {error}"))
}

fn record_error(record: u64, error: io::Error) -> CallError {
    CallError::LikelihoodVariantRecord {
        record,
        message: error.to_string(),
    }
}

fn output_error(error: io::Error) -> CallError {
    CallError::VariantOutput(error.to_string())
}

#[cfg(test)]
mod tests;
