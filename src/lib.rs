//! Typed likelihood evidence and calling stages for `rsomics-call`.

mod errmod;
mod error;
mod model;
mod snp;

pub use errmod::{BaseObservation, ErrorModel, LikelihoodMatrix, Nucleotide};
pub use error::{CallError, Result};
pub use model::{Allele, LikelihoodSite, Ploidy, SampleLikelihood};
pub use snp::{SnpEvidence, SnpLikelihoodConfig, SnpSiteBuilder};
