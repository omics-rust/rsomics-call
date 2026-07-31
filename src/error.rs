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

    #[error("invalid alignment region {region}: {message}")]
    InvalidRegion { region: String, message: String },

    #[error("at least one alignment region is required")]
    MissingRegions,

    #[error("target input {path}: {message}")]
    TargetInput { path: String, message: String },

    #[error("target input {path} line {line}: {message}")]
    TargetRecord {
        path: String,
        line: u64,
        message: String,
    },

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

    #[error("invalid indel likelihood configuration")]
    InvalidIndelConfig,

    #[error("invalid indel site summary")]
    InvalidIndelSummary,

    #[error("invalid likelihood annotations")]
    InvalidLikelihoodAnnotations,

    #[error("indel evidence exceeds the supported count range")]
    IndelEvidenceOverflow,

    #[error("reference data required for indel likelihoods is invalid")]
    InvalidIndelReference,

    #[error("indel realignment failed")]
    IndelRealignment,

    #[error("multiallelic mutation rate must be finite and between zero and one")]
    InvalidMutationRate,

    #[error("consensus reference probability threshold must be finite and between zero and one")]
    InvalidConsensusThreshold,

    #[error("the consensus caller requires observed compatible likelihoods")]
    UnsupportedConsensusLikelihoods,

    #[error("the multiallelic caller requires diploid input likelihoods")]
    UnsupportedLikelihoodPloidy,

    #[error("the caller ploidy count differs from the likelihood sample count")]
    CallerPloidyCountMismatch,

    #[error("the caller group count differs from the likelihood sample count")]
    CallerGroupCountMismatch,

    #[error("caller group indices must be contiguous and nonempty")]
    InvalidCallerGroups,

    #[error("the prior chromosome count is outside the cohort ploidy range")]
    InvalidPriorChromosomeCount,

    #[error("called allele count exceeds the supported range")]
    CalledAlleleCountOverflow,

    #[error("invalid likelihood VCF/BCF: {0}")]
    InvalidLikelihoodVariant(String),

    #[error("likelihood VCF/BCF input: {0}")]
    LikelihoodVariantInput(String),

    #[error("likelihood VCF/BCF record {record}: {message}")]
    LikelihoodVariantRecord { record: u64, message: String },

    #[error("likelihood call record {record}: {source}")]
    LikelihoodCallRecord { record: u64, source: Box<CallError> },

    #[error("variant VCF/BCF output: {0}")]
    VariantOutput(String),
}

pub type Result<T> = std::result::Result<T, CallError>;
