//! Typed likelihood evidence and calling stages for `rsomics-call`.

mod alignment;
mod caller;
mod errmod;
mod error;
mod model;
mod run;
mod samples;
mod snp;

pub use alignment::{AlignmentInput, AlignmentSet, ReferenceSequence};
pub use caller::{
    CallPloidy, CalledSample, CalledSite, MultiallelicCaller, MultiallelicCallerConfig,
};
pub use errmod::{BaseObservation, ErrorModel, LikelihoodMatrix, Nucleotide};
pub use error::{CallError, Result};
pub use model::{Allele, LikelihoodSite, Ploidy, SampleEvidence, SampleLikelihood};
pub use run::SnpLikelihoodRun;
pub use samples::{SampleMap, SampleMapBuilder, SampleSelection};
pub use snp::{SnpEvidence, SnpLikelihoodConfig, SnpSiteBuilder};
