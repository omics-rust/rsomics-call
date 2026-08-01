use std::io::{Read, Write};

use crate::{
    CallError, CallPloidy, CalledSite, CalledVariantWriter, CalledVcfSchema, ConsensusCaller,
    ConsensusCallerConfig, GvcfBlocker, LikelihoodSite, LikelihoodVariantReader,
    MultiallelicCaller, MultiallelicCallerConfig, PloidyResolver, Result, VariantOutputFormat,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CallModel {
    Multiallelic(MultiallelicCallerConfig),
    Consensus(ConsensusCallerConfig),
}

pub struct LikelihoodCallRun {
    caller: SiteCaller,
    ploidy: PloidyResolver,
    gvcf: Option<GvcfBlocker>,
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
        }
    }

    pub fn with_gvcf(mut self, thresholds: impl Into<Box<[u32]>>) -> Result<Self> {
        if matches!(&self.caller, SiteCaller::Consensus(_)) {
            return Err(CallError::UnsupportedGvcfModel);
        }
        self.gvcf = Some(GvcfBlocker::new(thresholds)?);
        Ok(self)
    }

    pub fn run<R, W>(
        mut self,
        mut reader: LikelihoodVariantReader<R>,
        output: W,
        format: VariantOutputFormat,
    ) -> Result<W>
    where
        R: Read,
        W: Write,
    {
        if self.ploidy.sample_count() != reader.schema().header().sample_names().len() {
            return Err(CallError::PloidySampleCountMismatch);
        }
        if self.gvcf.is_some() && !reader.schema().header().formats().contains_key("DP") {
            return Err(CallError::MissingGvcfDepth);
        }
        let schema = match &self.caller {
            SiteCaller::Multiallelic(_) => CalledVcfSchema::from_likelihood(reader.schema()),
            SiteCaller::Consensus(_) => CalledVcfSchema::from_consensus_likelihood(reader.schema()),
        };
        let schema = if self.gvcf.is_some() {
            schema.with_gvcf()
        } else {
            schema
        };
        let mut writer = CalledVariantWriter::new(output, schema, format)?;
        let mut ploidies = Vec::with_capacity(self.ploidy.sample_count());

        while let Some(site) = reader.read_site()? {
            let record = reader.record_number();
            let reference = reader
                .schema()
                .header()
                .contigs()
                .get_index(site.reference_sequence_id())
                .expect("decoded contig ID belongs to the reader schema")
                .0;
            self.ploidy
                .resolve_site_into(reference, &site, &mut ploidies)
                .map_err(|error| call_record_error(record, error))?;
            let called = self
                .caller
                .call(&site, &ploidies, self.ploidy.prior_chromosome_count())
                .map_err(|error| call_record_error(record, error))?;
            if let Some(blocker) = &mut self.gvcf {
                blocker.push(called, |called| writer.write_site(&called))?;
            } else {
                writer.write_site(&called)?;
            }
        }
        if let Some(blocker) = self.gvcf {
            blocker.finish(|called| writer.write_site(&called))?;
        }
        writer.finish()
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
