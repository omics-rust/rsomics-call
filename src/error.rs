use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum CallError {
    #[error("invalid allele: {0}")]
    InvalidAllele(String),

    #[error("ploidy must be greater than zero")]
    InvalidPloidy,

    #[error("allele and genotype likelihood dimensions are inconsistent")]
    InvalidLikelihoodDimensions,

    #[error("dependency correlation must be finite and between zero and one")]
    InvalidDependencyCorrelation,

    #[error("minimum base quality cannot exceed maximum base quality")]
    InvalidSnpQualityRange,

    #[error("at least one sample is required")]
    InvalidSampleCount,

    #[error("sample index {index} is outside the configured sample count {count}")]
    InvalidSampleIndex { index: usize, count: usize },

    #[error("pileup column has an invalid reference coordinate")]
    InvalidPileupCoordinate,

    #[error("SNP evidence exceeds the supported count range")]
    SnpEvidenceOverflow,

    #[error("the MAQ error model accepts at most 255 observations per sample")]
    ErrorModelDepth,
}

pub type Result<T> = std::result::Result<T, CallError>;
