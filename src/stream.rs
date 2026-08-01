use std::collections::HashSet;
use std::io::{self, BufWriter, Read, Write};

use noodles::{bcf, bgzf, vcf, vcf::variant::io::Write as _};
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
            return Err(CallError::LateLikelihoodSampleSelection);
        }
        Ok(())
    }
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
        self.schema = LikelihoodVcfSchema::from_header(header)?;
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
