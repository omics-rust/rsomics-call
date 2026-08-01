use crate::{Allele, CalledAnnotations, IndelSummary, SampleEvidence};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CallPloidy {
    Absent,
    Haploid,
    Diploid,
}

impl CallPloidy {
    pub(crate) fn chromosome_count(self) -> usize {
        match self {
            Self::Absent => 0,
            Self::Haploid => 1,
            Self::Diploid => 2,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CalledSample {
    pub(crate) ploidy: CallPloidy,
    pub(crate) genotype: Option<Box<[usize]>>,
    pub(crate) genotype_quality: Option<u8>,
    pub(crate) genotype_probabilities: Option<Box<[f32]>>,
    pub(crate) phred_likelihoods: Option<Box<[u32]>>,
    pub(crate) evidence: SampleEvidence,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GvcfSite {
    end_position: Option<u64>,
    minimum_depth: u32,
    collapsed: bool,
}

impl GvcfSite {
    pub(crate) fn new(end_position: Option<u64>, minimum_depth: u32, collapsed: bool) -> Self {
        Self {
            end_position,
            minimum_depth,
            collapsed,
        }
    }

    pub fn end_position(self) -> Option<u64> {
        self.end_position
    }

    pub fn minimum_depth(self) -> u32 {
        self.minimum_depth
    }

    pub fn is_collapsed(self) -> bool {
        self.collapsed
    }
}

impl CalledSample {
    pub fn ploidy(&self) -> CallPloidy {
        self.ploidy
    }

    pub fn genotype(&self) -> Option<&[usize]> {
        self.genotype.as_deref()
    }

    pub fn genotype_quality(&self) -> Option<u8> {
        self.genotype_quality
    }

    pub fn genotype_probabilities(&self) -> Option<&[f32]> {
        self.genotype_probabilities.as_deref()
    }

    pub fn phred_likelihoods(&self) -> Option<&[u32]> {
        self.phred_likelihoods.as_deref()
    }

    pub fn evidence(&self) -> &SampleEvidence {
        &self.evidence
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CalledSite {
    pub(crate) reference_sequence_id: usize,
    pub(crate) position: u64,
    pub(crate) reference: Allele,
    pub(crate) alternates: Box<[Allele]>,
    pub(crate) quality: Option<f32>,
    pub(crate) allele_counts: Box<[u32]>,
    pub(crate) samples: Box<[CalledSample]>,
    pub(crate) indel_summary: Option<IndelSummary>,
    pub(crate) annotations: Option<CalledAnnotations>,
    pub(crate) gvcf: Option<GvcfSite>,
}

impl CalledSite {
    pub fn reference_sequence_id(&self) -> usize {
        self.reference_sequence_id
    }

    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn reference(&self) -> &Allele {
        &self.reference
    }

    pub fn alternates(&self) -> &[Allele] {
        &self.alternates
    }

    pub fn quality(&self) -> Option<f32> {
        self.quality
    }

    pub fn allele_counts(&self) -> &[u32] {
        &self.allele_counts
    }

    pub fn allele_number(&self) -> u64 {
        self.allele_counts
            .iter()
            .map(|&count| u64::from(count))
            .sum()
    }

    pub fn samples(&self) -> &[CalledSample] {
        &self.samples
    }

    pub fn indel_summary(&self) -> Option<IndelSummary> {
        self.indel_summary
    }

    pub fn annotations(&self) -> Option<&CalledAnnotations> {
        self.annotations.as_ref()
    }

    pub fn gvcf(&self) -> Option<GvcfSite> {
        self.gvcf
    }

    pub fn is_variant(&self) -> bool {
        self.allele_counts[1..].iter().any(|&count| count != 0)
    }
}
