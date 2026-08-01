use std::io::{Read, Write};
use std::path::Path;

use noodles::core::Region;

use crate::{
    CallError, CallPloidy, CalledSite, CalledVariantWriter, CalledVcfSchema, ConsensusCaller,
    ConsensusCallerConfig, GvcfBlocker, IndexedLikelihoodVariantReader, LikelihoodSite,
    LikelihoodVariantReader, LikelihoodVcfSchema, MultiallelicCaller, MultiallelicCallerConfig,
    PloidyResolver, Result, VariantOutputFormat,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CallModel {
    Multiallelic(MultiallelicCallerConfig),
    Consensus(ConsensusCallerConfig),
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CallOutputOptions {
    variants_only: bool,
    keep_masked_reference: bool,
    skip_snps: bool,
    skip_indels: bool,
}

impl CallOutputOptions {
    pub fn with_variants_only(mut self, enabled: bool) -> Self {
        self.variants_only = enabled;
        self
    }

    pub fn with_keep_masked_reference(mut self, enabled: bool) -> Self {
        self.keep_masked_reference = enabled;
        self
    }

    pub fn with_skip_snps(mut self, enabled: bool) -> Self {
        self.skip_snps = enabled;
        self
    }

    pub fn with_skip_indels(mut self, enabled: bool) -> Self {
        self.skip_indels = enabled;
        self
    }

    fn accepts(self, source: &LikelihoodSite, called: &CalledSite) -> bool {
        if !self.keep_masked_reference && source.reference().as_bytes().contains(&b'N') {
            return false;
        }
        if source.indel_summary().is_some() {
            if self.skip_indels {
                return false;
            }
        } else if self.skip_snps {
            return false;
        }
        !self.variants_only || called.is_variant()
    }
}

pub struct LikelihoodCallRun {
    caller: SiteCaller,
    ploidy: PloidyResolver,
    gvcf: Option<GvcfBlocker>,
    output: CallOutputOptions,
}

enum SiteCaller {
    Multiallelic(MultiallelicCaller),
    Consensus(ConsensusCaller),
}

impl LikelihoodCallRun {
    pub fn new(model: CallModel, ploidy: PloidyResolver) -> Self {
        let caller = match model {
            CallModel::Multiallelic(config) => {
                SiteCaller::Multiallelic(MultiallelicCaller::new(config))
            }
            CallModel::Consensus(config) => SiteCaller::Consensus(ConsensusCaller::new(config)),
        };
        Self {
            caller,
            ploidy,
            gvcf: None,
            output: CallOutputOptions::default(),
        }
    }

    pub fn with_gvcf(mut self, thresholds: impl Into<Box<[u32]>>) -> Result<Self> {
        if matches!(&self.caller, SiteCaller::Consensus(_)) {
            return Err(CallError::UnsupportedGvcfModel);
        }
        self.gvcf = Some(GvcfBlocker::new(thresholds)?);
        Ok(self)
    }

    pub fn with_output_options(mut self, options: CallOutputOptions) -> Self {
        self.output = options;
        self
    }

    pub fn run<R, W>(
        self,
        mut reader: LikelihoodVariantReader<R>,
        output: W,
        format: VariantOutputFormat,
    ) -> Result<W>
    where
        R: Read,
        W: Write,
    {
        let mut run = self.start(reader.schema(), output, format)?;
        while let Some(site) = reader.read_site()? {
            run.push(reader.record_number(), &site)?;
        }
        run.finish()
    }

    pub fn run_indexed<W>(
        self,
        reader: IndexedLikelihoodVariantReader,
        regions: impl IntoIterator<Item = Region>,
        output: W,
        format: VariantOutputFormat,
    ) -> Result<W>
    where
        W: Write,
    {
        let regions = reader.normalize_regions(regions)?;
        let mut run = self.start(reader.schema(), output, format)?;
        reader.visit_normalized_regions(regions, |record, site| run.push(record, &site))?;
        run.finish()
    }

    pub fn run_indexed_regions_file<W>(
        self,
        reader: IndexedLikelihoodVariantReader,
        regions: impl AsRef<Path>,
        output: W,
        format: VariantOutputFormat,
    ) -> Result<W>
    where
        W: Write,
    {
        let regions = crate::region_file::read_regions(regions.as_ref())?;
        let regions = reader.normalize_region_file(regions)?;
        let mut run = self.start(reader.schema(), output, format)?;
        reader.visit_normalized_regions(regions, |record, site| run.push(record, &site))?;
        run.finish()
    }

    fn start<W>(
        self,
        input_schema: &LikelihoodVcfSchema,
        output: W,
        format: VariantOutputFormat,
    ) -> Result<ActiveCallRun<W>>
    where
        W: Write,
    {
        let sample_count = self.ploidy.sample_count();
        if sample_count != input_schema.header().sample_names().len() {
            return Err(CallError::PloidySampleCountMismatch);
        }
        if self.gvcf.is_some() && !input_schema.header().formats().contains_key("DP") {
            return Err(CallError::MissingGvcfDepth);
        }
        let schema = match &self.caller {
            SiteCaller::Multiallelic(_) => CalledVcfSchema::from_likelihood(input_schema),
            SiteCaller::Consensus(_) => CalledVcfSchema::from_consensus_likelihood(input_schema),
        };
        let schema = if self.gvcf.is_some() {
            schema.with_gvcf()
        } else {
            schema
        };
        Ok(ActiveCallRun {
            caller: self.caller,
            ploidy: self.ploidy,
            gvcf: self.gvcf,
            output: self.output,
            writer: CalledVariantWriter::new(output, schema, format)?,
            reference_names: input_schema
                .header()
                .contigs()
                .keys()
                .map(|name| Box::<str>::from(name.as_str()))
                .collect(),
            ploidies: Vec::with_capacity(sample_count),
        })
    }
}

struct ActiveCallRun<W>
where
    W: Write,
{
    caller: SiteCaller,
    ploidy: PloidyResolver,
    gvcf: Option<GvcfBlocker>,
    output: CallOutputOptions,
    writer: CalledVariantWriter<W>,
    reference_names: Box<[Box<str>]>,
    ploidies: Vec<CallPloidy>,
}

impl<W> ActiveCallRun<W>
where
    W: Write,
{
    fn push(&mut self, record: u64, site: &LikelihoodSite) -> Result<()> {
        let reference = self
            .reference_names
            .get(site.reference_sequence_id())
            .expect("decoded contig ID belongs to the reader schema");
        self.ploidy
            .resolve_site_into(reference, site, &mut self.ploidies)
            .map_err(|error| call_record_error(record, error))?;
        let called = self
            .caller
            .call(site, &self.ploidies, self.ploidy.prior_chromosome_count())
            .map_err(|error| call_record_error(record, error))?;
        if !self.output.accepts(site, &called) {
            return Ok(());
        }
        if let Some(blocker) = &mut self.gvcf {
            blocker.push(called, |called| self.writer.write_site(&called))
        } else {
            self.writer.write_site(&called)
        }
    }

    fn finish(mut self) -> Result<W> {
        if let Some(blocker) = self.gvcf {
            blocker.finish(|called| self.writer.write_site(&called))?;
        }
        self.writer.finish()
    }
}

impl SiteCaller {
    fn call(
        &self,
        site: &LikelihoodSite,
        ploidies: &[CallPloidy],
        prior_chromosome_count: usize,
    ) -> Result<CalledSite> {
        match self {
            Self::Multiallelic(caller) => {
                caller.call_with_ploidies(site, ploidies, prior_chromosome_count)
            }
            Self::Consensus(caller) => caller.call_with_ploidies(site, ploidies),
        }
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
        let called = call(&site).map_err(|error| call_record_error(record, error))?;
        writer.write_site(&called)?;
    }
    writer.finish()
}

fn call_record_error(record: u64, error: CallError) -> CallError {
    CallError::LikelihoodCallRecord {
        record,
        source: Box::new(error),
    }
}

#[cfg(test)]
mod tests;
