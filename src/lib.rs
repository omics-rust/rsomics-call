//! Typed likelihood evidence and calling stages for `rsomics-call`.

mod alignment;
mod annotation;
mod call_annotation;
mod called;
mod caller;
mod consensus;
mod errmod;
mod error;
mod format;
mod glocal;
mod indel;
mod model;
mod run;
mod samples;
mod selection;
mod snp;
mod stream;
mod target_file;

pub use alignment::{AlignmentInput, AlignmentSet, ReferenceSequence};
pub use call_annotation::CalledAnnotations;
pub use called::{CallPloidy, CalledSample, CalledSite};
pub use caller::{MultiallelicCaller, MultiallelicCallerConfig};
pub use consensus::{ConsensusCaller, ConsensusCallerConfig};
pub use errmod::{BaseObservation, ErrorModel, LikelihoodMatrix, Nucleotide};
pub use error::{CallError, Result};
pub use format::{CalledVcfSchema, LikelihoodVcfSchema};
pub use indel::{IndelLikelihoodConfig, IndelSiteBuilder};
pub use model::{
    Allele, IndelSummary, LikelihoodSite, Ploidy, SampleAnnotations, SampleEvidence,
    SampleLikelihood, SiteAnnotations,
};
pub use run::{SnpLikelihoodRun, run_likelihood_calls};
pub use samples::{SampleMap, SampleMapBuilder, SampleSelection};
pub use snp::{SnpEvidence, SnpLikelihoodConfig, SnpSiteBuilder};
pub use stream::{
    CalledVariantWriter, LikelihoodVariantReader, LikelihoodVariantWriter, VariantOutputFormat,
};
