use std::num::NonZeroU8;

use crate::{CallError, Result};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Allele(Box<[u8]>);

impl Allele {
    pub fn new(value: impl Into<Box<[u8]>>) -> Result<Self> {
        let value = value.into();
        let symbolic = value.len() > 2
            && value.starts_with(b"<")
            && value.ends_with(b">")
            && value[1..value.len() - 1]
                .iter()
                .all(|byte| byte.is_ascii_alphanumeric() || b"*_:.-".contains(byte));
        let valid = !value.is_empty()
            && (value.as_ref() == b"*"
                || symbolic
                || value
                    .iter()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N')));
        if valid {
            Ok(Self(value))
        } else {
            Err(CallError::InvalidAllele(
                String::from_utf8_lossy(&value).into_owned(),
            ))
        }
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Ploidy(NonZeroU8);

impl Ploidy {
    pub fn new(value: u8) -> Result<Self> {
        NonZeroU8::new(value)
            .map(Self)
            .ok_or(CallError::InvalidPloidy)
    }

    pub fn get(self) -> u8 {
        self.0.get()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleEvidence {
    depth: u32,
    allele_depths: Box<[u32]>,
    allele_quality_sums: Box<[u32]>,
}

impl SampleEvidence {
    pub fn new(
        depth: u32,
        allele_depths: impl Into<Box<[u32]>>,
        allele_quality_sums: impl Into<Box<[u32]>>,
    ) -> Result<Self> {
        let allele_depths = allele_depths.into();
        let allele_quality_sums = allele_quality_sums.into();
        if allele_depths.is_empty() || allele_quality_sums.len() != allele_depths.len() {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        Ok(Self {
            depth,
            allele_depths,
            allele_quality_sums,
        })
    }

    pub fn empty(allele_count: usize) -> Result<Self> {
        Self::new(0, vec![0; allele_count], vec![0; allele_count])
    }

    pub fn allele_depths(&self) -> &[u32] {
        &self.allele_depths
    }

    pub fn allele_quality_sums(&self) -> &[u32] {
        &self.allele_quality_sums
    }

    pub fn depth(&self) -> u32 {
        self.depth
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleLikelihood {
    ploidy: Ploidy,
    phred_likelihoods: Option<Box<[u32]>>,
    evidence: SampleEvidence,
}

impl SampleLikelihood {
    pub fn new(
        ploidy: Ploidy,
        phred_likelihoods: Option<Box<[u32]>>,
        evidence: SampleEvidence,
    ) -> Result<Self> {
        if phred_likelihoods.as_ref().is_some_and(|likelihoods| {
            genotype_count(evidence.allele_depths.len(), ploidy.get())
                .is_none_or(|count| count != likelihoods.len())
        }) {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        Ok(Self {
            ploidy,
            phred_likelihoods,
            evidence,
        })
    }

    pub fn observed(
        ploidy: Ploidy,
        phred_likelihoods: impl Into<Box<[u32]>>,
        evidence: SampleEvidence,
    ) -> Result<Self> {
        Self::new(ploidy, Some(phred_likelihoods.into()), evidence)
    }

    pub fn missing(allele_count: usize, ploidy: Ploidy) -> Result<Self> {
        Self::new(ploidy, None, SampleEvidence::empty(allele_count)?)
    }

    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    pub fn phred_likelihoods(&self) -> Option<&[u32]> {
        self.phred_likelihoods.as_deref()
    }

    pub fn evidence(&self) -> &SampleEvidence {
        &self.evidence
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndelSummary {
    maximum_support: u32,
    maximum_fraction: f32,
}

impl IndelSummary {
    pub fn new(maximum_support: u32, maximum_fraction: f32) -> Result<Self> {
        if !maximum_fraction.is_finite() || !(0.0..=1.0).contains(&maximum_fraction) {
            return Err(CallError::InvalidIndelSummary);
        }
        Ok(Self {
            maximum_support,
            maximum_fraction,
        })
    }

    pub fn maximum_support(self) -> u32 {
        self.maximum_support
    }

    pub fn maximum_fraction(self) -> f32 {
        self.maximum_fraction
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LikelihoodSite {
    reference_sequence_id: usize,
    position: u64,
    reference: Allele,
    alternates: Box<[Allele]>,
    allele_quality_sums: Box<[f32]>,
    samples: Box<[SampleLikelihood]>,
    indel_summary: Option<IndelSummary>,
}

impl LikelihoodSite {
    pub fn new(
        reference_sequence_id: usize,
        position: u64,
        reference: Allele,
        alternates: impl Into<Box<[Allele]>>,
        allele_quality_sums: impl Into<Box<[f32]>>,
        samples: impl Into<Box<[SampleLikelihood]>>,
    ) -> Result<Self> {
        let alternates = alternates.into();
        if alternates.is_empty()
            || alternates.iter().any(|allele| allele == &reference)
            || alternates
                .iter()
                .enumerate()
                .any(|(index, allele)| alternates[..index].contains(allele))
        {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        let allele_count = alternates.len() + 1;
        let allele_quality_sums = allele_quality_sums.into();
        if allele_quality_sums.len() != allele_count
            || allele_quality_sums
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        let samples = samples.into();
        if samples.iter().any(|sample| {
            sample.evidence().allele_depths().len() != allele_count
                || sample.phred_likelihoods().is_some_and(|likelihoods| {
                    genotype_count(allele_count, sample.ploidy().get())
                        .is_none_or(|count| count != likelihoods.len())
                })
        }) {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        Ok(Self {
            reference_sequence_id,
            position,
            reference,
            alternates,
            allele_quality_sums,
            samples,
            indel_summary: None,
        })
    }

    pub fn with_indel_summary(mut self, summary: IndelSummary) -> Self {
        self.indel_summary = Some(summary);
        self
    }

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

    pub fn allele_quality_sums(&self) -> &[f32] {
        &self.allele_quality_sums
    }

    pub fn samples(&self) -> &[SampleLikelihood] {
        &self.samples
    }

    pub fn indel_summary(&self) -> Option<IndelSummary> {
        self.indel_summary
    }
}

fn genotype_count(allele_count: usize, ploidy: u8) -> Option<usize> {
    let ploidy = usize::from(ploidy);
    let n = allele_count.checked_add(ploidy)?.checked_sub(1)?;
    let k = ploidy.min(n.checked_sub(ploidy)?);
    (1..=k).try_fold(1usize, |value, divisor| {
        value
            .checked_mul(n.checked_sub(k)?.checked_add(divisor)?)
            .map(|value| value / divisor)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_dimensions_follow_vcf_genotype_order() {
        let diploid = Ploidy::new(2).unwrap();
        let evidence = SampleEvidence::new(8, [3, 4, 1], [30, 40, 10]).unwrap();
        assert!(SampleLikelihood::observed(diploid, [0, 1, 2, 3, 4, 5], evidence.clone()).is_ok());
        assert_eq!(
            SampleLikelihood::observed(diploid, [0, 1, 2], evidence),
            Err(CallError::InvalidLikelihoodDimensions)
        );
        assert_eq!(
            SampleLikelihood::missing(3, diploid)
                .unwrap()
                .phred_likelihoods(),
            None
        );
    }

    #[test]
    fn sites_reject_duplicate_alleles() {
        let reference = Allele::new(&b"A"[..]).unwrap();
        let alternate = Allele::new(&b"G"[..]).unwrap();
        let evidence = SampleEvidence::new(5, [4, 1], [40, 10]).unwrap();
        let sample =
            SampleLikelihood::observed(Ploidy::new(2).unwrap(), [0, 10, 20], evidence).unwrap();
        assert!(
            LikelihoodSite::new(0, 9, reference.clone(), [alternate], [1.0, 0.2], [sample]).is_ok()
        );
        assert_eq!(
            LikelihoodSite::new(0, 9, reference.clone(), [reference], [1.0, 0.0], []),
            Err(CallError::InvalidLikelihoodDimensions)
        );
    }

    #[test]
    fn allele_and_missing_sample_boundaries_are_checked() {
        assert!(Allele::new(&b"<NON_REF>"[..]).is_ok());
        assert_eq!(
            Allele::new(&b"<>"[..]),
            Err(CallError::InvalidAllele("<>".to_owned()))
        );
        assert_eq!(
            SampleLikelihood::missing(0, Ploidy::new(2).unwrap()),
            Err(CallError::InvalidLikelihoodDimensions)
        );
        assert_eq!(
            SampleEvidence::new(1, [1, 0], [20]),
            Err(CallError::InvalidLikelihoodDimensions)
        );
        assert_eq!(
            LikelihoodSite::new(
                0,
                9,
                Allele::new(&b"A"[..]).unwrap(),
                [Allele::new(&b"G"[..]).unwrap()],
                [1.0, f32::NAN],
                []
            ),
            Err(CallError::InvalidLikelihoodDimensions)
        );
    }
}
