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

    #[error("invalid sample name: {0}")]
    InvalidSampleName(String),

    #[error("sample selection contains a duplicate name: {0}")]
    DuplicateSampleSelection(String),

    #[error("selected sample is absent from all inputs: {0}")]
    MissingSelectedSample(String),

    #[error("alignment source ID {0} is duplicated")]
    DuplicateAlignmentSource(u32),

    #[error("alignment source ID {0} is unknown")]
    UnknownAlignmentSource(u32),

    #[error("at least one alignment input is required")]
    MissingAlignmentInputs,

    #[error("alignment input {path}: {message}")]
    AlignmentInput { path: String, message: String },

    #[error("reference dictionary in alignment input {0} differs from the first input")]
    ReferenceDictionaryMismatch(String),

    #[error("reference input {path}: {message}")]
    ReferenceInput { path: String, message: String },

    #[error(transparent)]
    Pileup(#[from] rsomics_pileup::PileupError),

    #[error("invalid or duplicate read-group ID: {0}")]
    InvalidReadGroup(String),

    #[error("BAM RG field is not a valid Z string")]
    InvalidReadGroupField,

    #[error("pileup column has an invalid reference coordinate")]
    InvalidPileupCoordinate,

    #[error("SNP evidence exceeds the supported count range")]
    SnpEvidenceOverflow,

    #[error("multiallelic mutation rate must be finite and between zero and one")]
    InvalidMutationRate,

    #[error("the multiallelic caller currently requires diploid likelihoods")]
    UnsupportedCallerPloidy,

    #[error("called allele count exceeds the supported range")]
    CalledAlleleCountOverflow,
}

pub type Result<T> = std::result::Result<T, CallError>;
