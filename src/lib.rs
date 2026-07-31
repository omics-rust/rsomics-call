//! Typed likelihood evidence and calling stages for `rsomics-call`.

mod errmod;
mod error;
mod model;
mod samples;
mod snp;

pub use errmod::{BaseObservation, ErrorModel, LikelihoodMatrix, Nucleotide};
pub use error::{CallError, Result};
pub use model::{Allele, LikelihoodSite, Ploidy, SampleLikelihood};
pub use samples::{SampleMap, SampleMapBuilder, SampleSelection};
pub use snp::{SnpEvidence, SnpLikelihoodConfig, SnpSiteBuilder};
