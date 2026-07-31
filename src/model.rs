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
pub struct SampleLikelihood {
    ploidy: Ploidy,
    phred_likelihoods: Option<Box<[u32]>>,
    depth: u32,
    allele_depths: Box<[u32]>,
}

impl SampleLikelihood {
    pub fn observed(
        allele_count: usize,
        ploidy: Ploidy,
        phred_likelihoods: impl Into<Box<[u32]>>,
        depth: u32,
        allele_depths: impl Into<Box<[u32]>>,
    ) -> Result<Self> {
        let phred_likelihoods = phred_likelihoods.into();
        let allele_depths = allele_depths.into();
        if allele_depths.len() != allele_count
            || genotype_count(allele_count, ploidy.get())
                .is_none_or(|count| count != phred_likelihoods.len())
        {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        Ok(Self {
            ploidy,
            phred_likelihoods: Some(phred_likelihoods),
            depth,
            allele_depths,
        })
    }

    pub fn missing(allele_count: usize, ploidy: Ploidy) -> Result<Self> {
        if allele_count == 0 {
            return Err(CallError::InvalidLikelihoodDimensions);
        }
        Ok(Self {
            ploidy,
            phred_likelihoods: None,
            depth: 0,
            allele_depths: vec![0; allele_count].into_boxed_slice(),
        })
    }

    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    pub fn phred_likelihoods(&self) -> Option<&[u32]> {
        self.phred_likelihoods.as_deref()
    }

    pub fn depth(&self) -> u32 {
        self.depth
    }

    pub fn allele_depths(&self) -> &[u32] {
        &self.allele_depths
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LikelihoodSite {
    reference_sequence_id: usize,
    position: u64,
    reference: Allele,
    alternates: Box<[Allele]>,
    samples: Box<[SampleLikelihood]>,
}

impl LikelihoodSite {
    pub fn new(
        reference_sequence_id: usize,
        position: u64,
        reference: Allele,
        alternates: impl Into<Box<[Allele]>>,
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
        let samples = samples.into();
        if samples.iter().any(|sample| {
            sample.allele_depths().len() != allele_count
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
            samples,
        })
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

    pub fn samples(&self) -> &[SampleLikelihood] {
        &self.samples
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
        assert!(SampleLikelihood::observed(3, diploid, [0, 1, 2, 3, 4, 5], 8, [3, 4, 1]).is_ok());
        assert_eq!(
            SampleLikelihood::observed(3, diploid, [0, 1, 2], 8, [3, 4, 1]),
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
        let sample =
            SampleLikelihood::observed(2, Ploidy::new(2).unwrap(), [0, 10, 20], 5, [4, 1]).unwrap();
        assert!(LikelihoodSite::new(0, 9, reference.clone(), [alternate], [sample]).is_ok());
        assert_eq!(
            LikelihoodSite::new(0, 9, reference.clone(), [reference], []),
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
    }
}
