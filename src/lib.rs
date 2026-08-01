//! Typed likelihood evidence and calling stages for `rsomics-call`.

mod alignment;
mod annotation;
mod call_annotation;
mod call_run;
mod call_samples;
mod called;
mod caller;
mod consensus;
mod errmod;
mod error;
mod format;
mod glocal;
mod gvcf;
mod indel;
mod model;
mod ploidy;
mod region_file;
mod run;
mod samples;
mod selection;
mod snp;
mod stream;

pub use alignment::{AlignmentInput, AlignmentSet, ReferenceSequence};
pub use call_annotation::CalledAnnotations;
pub use call_run::{CallModel, CallOutputOptions, LikelihoodCallRun, run_likelihood_calls};
pub use call_samples::CallSampleSelection;
pub use called::{CallPloidy, CalledSample, CalledSite, GvcfSite};
pub use caller::{MultiallelicCaller, MultiallelicCallerConfig};
pub use consensus::{ConsensusCaller, ConsensusCallerConfig};
pub use errmod::{BaseObservation, ErrorModel, LikelihoodMatrix, Nucleotide};
pub use error::{CallError, Result};
pub use format::{CalledVcfSchema, LikelihoodVcfSchema};
pub use gvcf::GvcfBlocker;
pub use indel::{IndelAmbiguousReadPolicy, IndelLikelihoodConfig, IndelSiteBuilder};
pub use model::{
    Allele, IndelSummary, LikelihoodSite, Ploidy, PriorAlleleCounts, SampleAnnotations,
    SampleEvidence, SampleLikelihood, SiteAnnotations,
};
pub use ploidy::{PloidyDefinition, PloidyPreset, PloidyResolver, SamplePloidy};
pub use run::SnpLikelihoodRun;
pub use samples::{SampleMap, SampleMapBuilder, SampleSelection};
pub use snp::{SnpEvidence, SnpLikelihoodConfig, SnpSiteBuilder};
pub use stream::{
    CalledVariantWriter, IndexedLikelihoodVariantReader, LikelihoodVariantReader,
    LikelihoodVariantWriter, VariantOutputFormat,
};
