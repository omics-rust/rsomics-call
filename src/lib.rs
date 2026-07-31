//! Typed likelihood evidence and calling stages for `rsomics-call`.

mod alignment;
mod called;
mod caller;
mod consensus;
mod errmod;
mod error;
mod format;
mod model;
mod run;
mod samples;
mod snp;
mod stream;

pub use alignment::{AlignmentInput, AlignmentSet, ReferenceSequence};
pub use called::{CallPloidy, CalledSample, CalledSite};
pub use caller::{MultiallelicCaller, MultiallelicCallerConfig};
pub use consensus::{ConsensusCaller, ConsensusCallerConfig};
pub use errmod::{BaseObservation, ErrorModel, LikelihoodMatrix, Nucleotide};
pub use error::{CallError, Result};
pub use format::{CalledVcfSchema, LikelihoodVcfSchema};
pub use model::{Allele, LikelihoodSite, Ploidy, SampleEvidence, SampleLikelihood};
pub use run::{SnpLikelihoodRun, run_likelihood_calls};
pub use samples::{SampleMap, SampleMapBuilder, SampleSelection};
pub use snp::{SnpEvidence, SnpLikelihoodConfig, SnpSiteBuilder};
pub use stream::{
    CalledVariantWriter, LikelihoodVariantReader, LikelihoodVariantWriter, VariantOutputFormat,
};
