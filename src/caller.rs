use crate::{Allele, CallError, LikelihoodSite, Result, SampleEvidence, SampleLikelihood};

const PHRED_SCALE: f64 = 4.342_94;
const SITE_PHRED_SCALE: f64 = 4.343;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiallelicCallerConfig {
    mutation_rate: f64,
}

impl MultiallelicCallerConfig {
    pub fn new(mutation_rate: f64) -> Result<Self> {
        if !mutation_rate.is_finite() || mutation_rate <= 0.0 || mutation_rate >= 1.0 {
            return Err(CallError::InvalidMutationRate);
        }
        Ok(Self { mutation_rate })
    }
}

impl Default for MultiallelicCallerConfig {
    fn default() -> Self {
        Self {
            mutation_rate: 1.1e-3,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CallPloidy {
    Absent,
    Haploid,
    Diploid,
}

impl CallPloidy {
    fn chromosome_count(self) -> usize {
        match self {
            Self::Absent => 0,
            Self::Haploid => 1,
            Self::Diploid => 2,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CalledSample {
    ploidy: CallPloidy,
    genotype: Option<Box<[usize]>>,
    genotype_quality: u8,
    genotype_probabilities: Option<Box<[f32]>>,
    phred_likelihoods: Option<Box<[u32]>>,
    evidence: SampleEvidence,
}

impl CalledSample {
    pub fn ploidy(&self) -> CallPloidy {
        self.ploidy
    }

    pub fn genotype(&self) -> Option<&[usize]> {
        self.genotype.as_deref()
    }

    pub fn genotype_quality(&self) -> u8 {
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
    reference_sequence_id: usize,
    position: u64,
    reference: Allele,
    alternates: Box<[Allele]>,
    quality: Option<f32>,
    allele_counts: Box<[u32]>,
    samples: Box<[CalledSample]>,
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

    pub fn is_variant(&self) -> bool {
        self.allele_counts[1..].iter().any(|&count| count != 0)
    }
}

pub struct MultiallelicCaller {
    config: MultiallelicCallerConfig,
}

impl MultiallelicCaller {
    pub fn new(config: MultiallelicCallerConfig) -> Self {
        Self { config }
    }

    pub fn call(&self, site: &LikelihoodSite) -> Result<CalledSite> {
        self.call_inner(
            site,
            SamplePloidies::Diploid,
            site.samples().len().saturating_mul(2),
        )
    }

    /// `prior_chromosome_count` is the cohort maximum and remains fixed across sites.
    pub fn call_with_ploidies(
        &self,
        site: &LikelihoodSite,
        ploidies: &[CallPloidy],
        prior_chromosome_count: usize,
    ) -> Result<CalledSite> {
        if ploidies.len() != site.samples().len() {
            return Err(CallError::CallerPloidyCountMismatch);
        }
        self.call_inner(
            site,
            SamplePloidies::Explicit(ploidies),
            prior_chromosome_count,
        )
    }

    fn call_inner(
        &self,
        site: &LikelihoodSite,
        ploidies: SamplePloidies<'_>,
        prior_chromosome_count: usize,
    ) -> Result<CalledSite> {
        if site
            .samples()
            .iter()
            .any(|sample| sample.ploidy().get() != 2)
        {
            return Err(CallError::UnsupportedLikelihoodPloidy);
        }
        let active_chromosome_count = ploidies.chromosome_count(site.samples().len());
        let maximum_chromosome_count = site.samples().len().saturating_mul(2);
        if prior_chromosome_count < active_chromosome_count
            || prior_chromosome_count > maximum_chromosome_count
        {
            return Err(CallError::InvalidPriorChromosomeCount);
        }

        let allele_count = site.alternates().len() + 1;
        let genotype_count = allele_count * (allele_count + 1) / 2;
        let likelihoods = site
            .samples()
            .iter()
            .map(|sample| normalized_likelihoods(sample, genotype_count))
            .collect::<Vec<_>>();
        let quality_sums = normalized_quality_sums(site.allele_quality_sums());
        let log_prior = adjusted_log_prior(self.config.mutation_rate, prior_chromosome_count);
        let selection = select_alleles(&likelihoods, &quality_sums, &ploidies, log_prior);
        let unseen = site
            .alternates()
            .iter()
            .position(|allele| matches!(allele.as_bytes(), b"<*>" | b"<NON_REF>"))
            .map(|index| index + 1);
        let mut retained = selection.alleles;
        retained[0] = true;
        if let Some(index) = unseen {
            retained[index] = false;
        }
        let retained = retained
            .iter()
            .enumerate()
            .filter_map(|(index, &keep)| keep.then_some(index))
            .collect::<Vec<_>>();
        let old_to_new = allele_map(allele_count, &retained);
        let mut allele_counts = vec![0u32; retained.len()];
        let mut samples = Vec::with_capacity(site.samples().len());

        for (sample_index, sample) in site.samples().iter().enumerate() {
            let ploidy = ploidies.get(sample_index);
            if ploidy == CallPloidy::Absent {
                samples.push(call_absent_sample(sample, &retained)?);
            } else if retained.len() == 1 {
                samples.push(call_reference_sample(
                    sample,
                    &likelihoods[sample_index],
                    ploidy,
                    &mut allele_counts,
                )?);
            } else {
                samples.push(call_sample(
                    sample,
                    &likelihoods[sample_index],
                    &quality_sums,
                    &retained,
                    &old_to_new,
                    ploidy,
                    &mut allele_counts,
                )?);
            }
        }

        let alternate_count = allele_counts[1..].iter().sum::<u32>();
        let quality = if alternate_count != 0 {
            selection.maximum_quality
        } else if selection.alternate_log_sum.is_finite() {
            Some(
                (-SITE_PHRED_SCALE
                    * (selection.alternate_log_sum
                        - logsumexp(
                            selection.alternate_log_sum,
                            selection.reference_log_likelihood,
                        ))) as f32,
            )
        } else if allele_counts[0] != 0 {
            Some((-SITE_PHRED_SCALE * log_prior) as f32)
        } else {
            None
        };
        let mut alleles = std::iter::once(site.reference())
            .chain(site.alternates())
            .enumerate()
            .filter_map(|(index, allele)| old_to_new[index].map(|_| allele.clone()));
        let reference = alleles.next().unwrap();
        let alternates = alleles.collect::<Vec<_>>().into_boxed_slice();

        Ok(CalledSite {
            reference_sequence_id: site.reference_sequence_id(),
            position: site.position(),
            reference,
            alternates,
            quality,
            allele_counts: allele_counts.into_boxed_slice(),
            samples: samples.into_boxed_slice(),
        })
    }
}

impl Default for MultiallelicCaller {
    fn default() -> Self {
        Self::new(MultiallelicCallerConfig::default())
    }
}

#[derive(Clone, Copy)]
enum SamplePloidies<'a> {
    Diploid,
    Explicit(&'a [CallPloidy]),
}

impl SamplePloidies<'_> {
    fn get(self, index: usize) -> CallPloidy {
        match self {
            Self::Diploid => CallPloidy::Diploid,
            Self::Explicit(values) => values[index],
        }
    }

    fn chromosome_count(self, sample_count: usize) -> usize {
        match self {
            Self::Diploid => sample_count.saturating_mul(2),
            Self::Explicit(values) => values.iter().fold(0, |count, ploidy| {
                count.saturating_add(ploidy.chromosome_count())
            }),
        }
    }
}

struct AlleleSelection {
    alleles: Vec<bool>,
    reference_log_likelihood: f64,
    alternate_log_sum: f64,
    maximum_quality: Option<f32>,
}

fn select_alleles(
    sample_likelihoods: &[Vec<f64>],
    quality_sums: &[f64],
    ploidies: &SamplePloidies<'_>,
    log_prior: f64,
) -> AlleleSelection {
    let allele_count = quality_sums.len();
    let mut best = vec![false; allele_count];
    let mut maximum_log_likelihood = f64::NEG_INFINITY;
    let mut reference_log_likelihood = f64::NEG_INFINITY;
    let mut alternate_log_sum = f64::NEG_INFINITY;

    for allele in 0..allele_count {
        let index = genotype_index(allele, allele);
        let mut total = 0.0;
        let mut supported = false;
        for likelihoods in sample_likelihoods {
            let value = likelihoods[index];
            if value != 0.0 {
                total += value.ln();
                supported = true;
            }
        }
        if allele == 0 {
            reference_log_likelihood = total;
        } else {
            total += log_prior;
        }
        update_selection(
            &[allele],
            total,
            supported,
            allele != 0 && supported,
            &mut best,
            &mut maximum_log_likelihood,
            &mut alternate_log_sum,
        );
    }

    for first in 0..allele_count {
        if quality_sums[first] == 0.0 {
            continue;
        }
        for second in 0..first {
            if quality_sums[second] == 0.0 {
                continue;
            }
            let frequencies = normalized_frequencies(quality_sums, &[first, second]);
            let mut total = 0.0;
            let mut supported = false;
            for (sample_index, likelihoods) in sample_likelihoods.iter().enumerate() {
                let ploidy = ploidies.get(sample_index);
                let value = match ploidy {
                    CallPloidy::Absent => 0.0,
                    CallPloidy::Haploid => {
                        frequencies[0] * likelihoods[genotype_index(first, first)]
                            + frequencies[1] * likelihoods[genotype_index(second, second)]
                    }
                    CallPloidy::Diploid => {
                        frequencies[0].powi(2) * likelihoods[genotype_index(first, first)]
                            + frequencies[1].powi(2) * likelihoods[genotype_index(second, second)]
                            + 2.0
                                * frequencies[0]
                                * frequencies[1]
                                * likelihoods[genotype_index(first, second)]
                    }
                };
                if value != 0.0 {
                    total += value.ln();
                    supported = true;
                }
            }
            if first != 0 {
                total += log_prior;
            }
            if second != 0 {
                total += log_prior;
            }
            update_selection(
                &[first, second],
                total,
                supported,
                supported,
                &mut best,
                &mut maximum_log_likelihood,
                &mut alternate_log_sum,
            );
        }
    }

    for first in 0..allele_count {
        if quality_sums[first] == 0.0 {
            continue;
        }
        for second in 0..first {
            if quality_sums[second] == 0.0 {
                continue;
            }
            for third in 0..second {
                if quality_sums[third] == 0.0 {
                    continue;
                }
                let alleles = [first, second, third];
                let frequencies = normalized_frequencies(quality_sums, &alleles);
                let mut total = 0.0;
                let mut supported = false;
                for (sample_index, likelihoods) in sample_likelihoods.iter().enumerate() {
                    let ploidy = ploidies.get(sample_index);
                    let mut value = 0.0;
                    for (right, &right_allele) in alleles.iter().enumerate() {
                        match ploidy {
                            CallPloidy::Absent => {}
                            CallPloidy::Haploid => {
                                value += frequencies[right]
                                    * likelihoods[genotype_index(right_allele, right_allele)];
                            }
                            CallPloidy::Diploid => {
                                value += frequencies[right].powi(2)
                                    * likelihoods[genotype_index(right_allele, right_allele)];
                                for left in 0..right {
                                    value += 2.0
                                        * frequencies[left]
                                        * frequencies[right]
                                        * likelihoods[genotype_index(alleles[left], right_allele)];
                                }
                            }
                        }
                    }
                    if value != 0.0 {
                        total += value.ln();
                        supported = true;
                    }
                }
                for &allele in &alleles {
                    if allele != 0 {
                        total += log_prior;
                    }
                }
                update_selection(
                    &alleles,
                    total,
                    supported,
                    supported,
                    &mut best,
                    &mut maximum_log_likelihood,
                    &mut alternate_log_sum,
                );
            }
        }
    }

    let maximum_quality = maximum_log_likelihood.is_finite().then(|| {
        (-SITE_PHRED_SCALE
            * (reference_log_likelihood - logsumexp(alternate_log_sum, reference_log_likelihood)))
            as f32
    });
    AlleleSelection {
        alleles: best,
        reference_log_likelihood,
        alternate_log_sum,
        maximum_quality,
    }
}

fn update_selection(
    alleles: &[usize],
    log_likelihood: f64,
    supported: bool,
    include_in_sum: bool,
    best: &mut [bool],
    maximum: &mut f64,
    alternate_sum: &mut f64,
) {
    if supported && *maximum < log_likelihood {
        best.fill(false);
        for &allele in alleles {
            best[allele] = true;
        }
        *maximum = log_likelihood;
    }
    if include_in_sum {
        *alternate_sum = logsumexp(log_likelihood, *alternate_sum);
    }
}

fn call_sample(
    sample: &SampleLikelihood,
    likelihoods: &[f64],
    quality_sums: &[f64],
    retained: &[usize],
    old_to_new: &[Option<usize>],
    ploidy: CallPloidy,
    allele_counts: &mut [u32],
) -> Result<CalledSample> {
    let (phred_likelihoods, evidence) = trimmed_fields(sample, retained, ploidy, true)?;
    if likelihoods.iter().all(|&value| value == 0.0) {
        return Ok(CalledSample {
            ploidy,
            genotype: None,
            genotype_quality: 0,
            genotype_probabilities: None,
            phred_likelihoods,
            evidence,
        });
    }

    let genotype_count = if ploidy == CallPloidy::Diploid {
        retained.len() * (retained.len() + 1) / 2
    } else {
        retained.len()
    };
    let mut probabilities = vec![0.0f32; genotype_count];
    let mut best = 0.0;
    let mut genotype = vec![0usize; ploidy.chromosome_count()];
    for &right in retained {
        let new_right = old_to_new[right].unwrap();
        let homozygous = if ploidy == CallPloidy::Diploid {
            likelihoods[genotype_index(right, right)] * quality_sums[right].powi(2)
        } else {
            likelihoods[genotype_index(right, right)] * quality_sums[right]
        };
        let output_index = if ploidy == CallPloidy::Diploid {
            genotype_index(new_right, new_right)
        } else {
            new_right
        };
        probabilities[output_index] = homozygous as f32;
        if best < homozygous {
            best = homozygous;
            genotype.fill(new_right);
        }
        if ploidy == CallPloidy::Diploid {
            for &left in retained.iter().take_while(|&&left| left != right) {
                let new_left = old_to_new[left].unwrap();
                let heterozygous = 2.0
                    * likelihoods[genotype_index(right, left)]
                    * quality_sums[right]
                    * quality_sums[left];
                probabilities[genotype_index(new_right, new_left)] = heterozygous as f32;
                if best < heterozygous {
                    best = heterozygous;
                    genotype[0] = new_left;
                    genotype[1] = new_right;
                }
            }
        }
    }
    for &allele in &genotype {
        allele_counts[allele] = allele_counts[allele]
            .checked_add(1)
            .ok_or(CallError::CalledAlleleCountOverflow)?;
    }
    let sum = probabilities
        .iter()
        .map(|&value| f64::from(value))
        .sum::<f64>();
    let maximum = probabilities
        .iter()
        .copied()
        .map(f64::from)
        .fold(f64::NEG_INFINITY, f64::max);
    let quality = if sum == 0.0 {
        0
    } else {
        let value = -PHRED_SCALE * (1.0 - maximum / sum).ln();
        if value <= f64::from(i8::MAX) {
            value as u8
        } else {
            i8::MAX as u8
        }
    };
    let probabilities = (sum != 0.0).then(|| {
        probabilities
            .into_iter()
            .map(|value| (f64::from(value) / sum) as f32)
            .collect::<Vec<_>>()
            .into_boxed_slice()
    });
    Ok(CalledSample {
        ploidy,
        genotype: Some(genotype.into()),
        genotype_quality: quality,
        genotype_probabilities: probabilities,
        phred_likelihoods,
        evidence,
    })
}

fn call_absent_sample(sample: &SampleLikelihood, retained: &[usize]) -> Result<CalledSample> {
    let (phred_likelihoods, evidence) =
        trimmed_fields(sample, retained, CallPloidy::Absent, false)?;
    Ok(CalledSample {
        ploidy: CallPloidy::Absent,
        genotype: None,
        genotype_quality: 0,
        genotype_probabilities: None,
        phred_likelihoods,
        evidence,
    })
}

fn call_reference_sample(
    sample: &SampleLikelihood,
    likelihoods: &[f64],
    ploidy: CallPloidy,
    allele_counts: &mut [u32],
) -> Result<CalledSample> {
    let supported = likelihoods.iter().any(|&value| value != 0.0);
    if supported {
        allele_counts[0] = allele_counts[0]
            .checked_add(ploidy.chromosome_count() as u32)
            .ok_or(CallError::CalledAlleleCountOverflow)?;
    }
    let (phred_likelihoods, evidence) = trimmed_fields(sample, &[0], ploidy, false)?;
    Ok(CalledSample {
        ploidy,
        genotype: supported.then(|| vec![0; ploidy.chromosome_count()].into_boxed_slice()),
        genotype_quality: 0,
        genotype_probabilities: None,
        phred_likelihoods,
        evidence,
    })
}

fn trimmed_fields(
    sample: &SampleLikelihood,
    retained: &[usize],
    ploidy: CallPloidy,
    retain_likelihoods: bool,
) -> Result<(Option<Box<[u32]>>, SampleEvidence)> {
    let evidence = sample.evidence();
    let allele_depths = retained
        .iter()
        .map(|&index| evidence.allele_depths()[index])
        .collect::<Vec<_>>();
    let allele_quality_sums = retained
        .iter()
        .map(|&index| evidence.allele_quality_sums()[index])
        .collect::<Vec<_>>();
    let evidence = SampleEvidence::new(evidence.depth(), allele_depths, allele_quality_sums)?;
    let phred_likelihoods = retain_likelihoods
        .then(|| {
            sample.phred_likelihoods().map(|values| {
                let mut trimmed = if ploidy == CallPloidy::Diploid {
                    Vec::with_capacity(retained.len() * (retained.len() + 1) / 2)
                } else {
                    Vec::with_capacity(retained.len())
                };
                for (right_index, &right) in retained.iter().enumerate() {
                    if ploidy == CallPloidy::Diploid {
                        for &left in &retained[..=right_index] {
                            trimmed.push(values[genotype_index(right, left)]);
                        }
                    } else {
                        trimmed.push(values[genotype_index(right, right)]);
                    }
                }
                trimmed.into_boxed_slice()
            })
        })
        .flatten();
    Ok((phred_likelihoods, evidence))
}

fn normalized_likelihoods(sample: &SampleLikelihood, genotype_count: usize) -> Vec<f64> {
    let Some(values) = sample.phred_likelihoods() else {
        return vec![0.0; genotype_count];
    };
    let mut likelihoods = values
        .iter()
        .map(|&value| 10.0f64.powf(-f64::from(value) / 10.0))
        .collect::<Vec<_>>();
    let sum = likelihoods.iter().sum::<f64>();
    if sum != 0.0 {
        for value in &mut likelihoods {
            *value /= sum;
        }
    }
    likelihoods
}

fn normalized_quality_sums(values: &[f32]) -> Vec<f64> {
    let total = values.iter().sum::<f32>();
    if total == 0.0 {
        return vec![0.0; values.len()];
    }
    values
        .iter()
        .map(|&value| f64::from(value / total))
        .collect()
}

fn adjusted_log_prior(mutation_rate: f64, chromosome_count: usize) -> f64 {
    let mut factor = 1.0;
    for count in 2..chromosome_count {
        factor += 1.0 / count as f64;
    }
    (mutation_rate * factor).min(0.99).ln()
}

fn normalized_frequencies(quality_sums: &[f64], alleles: &[usize]) -> Vec<f64> {
    let total = alleles
        .iter()
        .map(|&allele| quality_sums[allele])
        .sum::<f64>();
    alleles
        .iter()
        .map(|&allele| quality_sums[allele] / total)
        .collect()
}

fn allele_map(allele_count: usize, retained: &[usize]) -> Vec<Option<usize>> {
    let mut map = vec![None; allele_count];
    for (new, &old) in retained.iter().enumerate() {
        map[old] = Some(new);
    }
    map
}

fn genotype_index(first: usize, second: usize) -> usize {
    let (left, right) = if first <= second {
        (first, second)
    } else {
        (second, first)
    };
    right * (right + 1) / 2 + left
}

fn logsumexp(first: f64, second: f64) -> f64 {
    if first > second {
        (1.0 + (second - first).exp()).ln() + first
    } else {
        (1.0 + (first - second).exp()).ln() + second
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Ploidy;

    fn sample(phred_likelihoods: [u32; 6], depths: [u32; 3], sums: [u32; 3]) -> SampleLikelihood {
        SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            phred_likelihoods,
            SampleEvidence::new(1, depths, sums).unwrap(),
        )
        .unwrap()
    }

    fn two_sample_site() -> LikelihoodSite {
        LikelihoodSite::new(
            0,
            0,
            Allele::new(&b"A"[..]).unwrap(),
            [
                Allele::new(&b"G"[..]).unwrap(),
                Allele::new(&b"<*>"[..]).unwrap(),
            ],
            [1.0, 1.0, 0.0],
            [
                sample([0, 3, 40, 3, 40, 40], [1, 0, 0], [40, 0, 0]),
                sample([40, 3, 0, 40, 3, 40], [0, 1, 0], [0, 40, 0]),
            ],
        )
        .unwrap()
    }

    #[test]
    fn matches_bcftools_1_24_two_sample_multiallelic_call() {
        let site = two_sample_site();

        let called = MultiallelicCaller::default().call(&site).unwrap();

        assert_eq!(called.reference().as_bytes(), b"A");
        assert_eq!(called.alternates()[0].as_bytes(), b"G");
        assert_eq!(called.allele_counts(), &[2, 2]);
        assert_eq!(called.allele_number(), 4);
        assert!(called.is_variant());
        assert!((called.quality().unwrap() - 7.822_08).abs() < 1e-5);
        assert_eq!(called.samples()[0].genotype(), Some(&[0, 1][..]));
        assert_eq!(called.samples()[1].genotype(), Some(&[0, 1][..]));
        assert_eq!(called.samples()[0].genotype_quality(), 3);
        assert_eq!(called.samples()[1].genotype_quality(), 3);
        assert_eq!(
            called.samples()[0].phred_likelihoods(),
            Some(&[0, 3, 40][..])
        );
        assert_eq!(
            called.samples()[1].phred_likelihoods(),
            Some(&[40, 3, 0][..])
        );
        let expected = [
            [0.499_382, 0.500_568, 4.993_82e-5],
            [4.993_82e-5, 0.500_568, 0.499_382],
        ];
        for (sample, expected) in called.samples().iter().zip(expected) {
            for (&observed, expected) in sample
                .genotype_probabilities()
                .unwrap()
                .iter()
                .zip(expected)
            {
                assert!((observed - expected).abs() < 5e-7);
            }
        }
    }

    #[test]
    fn rejects_non_diploid_likelihoods_and_invalid_prior() {
        assert_eq!(
            MultiallelicCallerConfig::new(f64::NAN),
            Err(CallError::InvalidMutationRate)
        );
        let site = LikelihoodSite::new(
            0,
            0,
            Allele::new(&b"A"[..]).unwrap(),
            [Allele::new(&b"G"[..]).unwrap()],
            [1.0, 1.0],
            [SampleLikelihood::observed(
                Ploidy::new(1).unwrap(),
                [0, 40],
                SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
            )
            .unwrap()],
        )
        .unwrap();
        assert_eq!(
            MultiallelicCaller::default().call(&site),
            Err(CallError::UnsupportedLikelihoodPloidy)
        );
    }

    #[test]
    fn matches_bcftools_1_24_haploid_call() {
        let site = two_sample_site();
        let haploid = CallPloidy::Haploid;

        let called = MultiallelicCaller::default()
            .call_with_ploidies(&site, &[haploid, haploid], 2)
            .unwrap();

        assert_eq!(called.allele_counts(), &[1, 1]);
        assert_eq!(called.allele_number(), 2);
        assert!((called.quality().unwrap() - 5.7423).abs() < 1e-4);
        assert_eq!(called.samples()[0].genotype(), Some(&[0][..]));
        assert_eq!(called.samples()[1].genotype(), Some(&[1][..]));
        assert_eq!(called.samples()[0].genotype_quality(), 40);
        assert_eq!(called.samples()[1].genotype_quality(), 40);
        assert_eq!(called.samples()[0].phred_likelihoods(), Some(&[0, 40][..]));
        assert_eq!(called.samples()[1].phred_likelihoods(), Some(&[40, 0][..]));
        let expected = [[0.9999, 9.999e-5], [9.999e-5, 0.9999]];
        for (sample, expected) in called.samples().iter().zip(expected) {
            for (&observed, expected) in sample
                .genotype_probabilities()
                .unwrap()
                .iter()
                .zip(expected)
            {
                assert!((observed - expected).abs() < 5e-7);
            }
        }
    }

    #[test]
    fn matches_bcftools_1_24_mixed_ploidy_call() {
        let site = two_sample_site();
        let haploid = CallPloidy::Haploid;
        let diploid = CallPloidy::Diploid;

        assert_eq!(
            MultiallelicCaller::default().call_with_ploidies(&site, &[haploid, diploid], 2),
            Err(CallError::InvalidPriorChromosomeCount)
        );
        assert_eq!(
            MultiallelicCaller::default().call_with_ploidies(&site, &[haploid, diploid], 5),
            Err(CallError::InvalidPriorChromosomeCount)
        );
        let called = MultiallelicCaller::default()
            .call_with_ploidies(&site, &[haploid, diploid], 4)
            .unwrap();

        assert_eq!(called.allele_counts(), &[2, 1]);
        assert_eq!(called.allele_number(), 3);
        let quality = called.quality().unwrap();
        assert!((quality - 7.817_96).abs() < 1e-5, "{quality}");
        assert_eq!(called.samples()[0].genotype(), Some(&[0][..]));
        assert_eq!(called.samples()[1].genotype(), Some(&[0, 1][..]));
        assert_eq!(called.samples()[0].genotype_quality(), 40);
        assert_eq!(called.samples()[1].genotype_quality(), 3);
    }

    #[test]
    fn matches_bcftools_1_24_absent_ploidy_call() {
        let site = two_sample_site();

        let called = MultiallelicCaller::default()
            .call_with_ploidies(&site, &[CallPloidy::Absent, CallPloidy::Diploid], 4)
            .unwrap();

        assert_eq!(called.allele_counts(), &[1, 1]);
        assert_eq!(called.allele_number(), 2);
        assert!((called.quality().unwrap() - 13.2678).abs() < 1e-4);
        assert_eq!(called.samples()[0].ploidy(), CallPloidy::Absent);
        assert_eq!(called.samples()[0].genotype(), None);
        assert_eq!(called.samples()[0].phred_likelihoods(), None);
        assert_eq!(called.samples()[0].evidence().allele_depths(), &[1, 0]);
        assert_eq!(called.samples()[1].genotype(), Some(&[0, 1][..]));
        assert_eq!(called.samples()[1].genotype_quality(), 3);
        assert_eq!(
            called.samples()[1].phred_likelihoods(),
            Some(&[40, 3, 0][..])
        );
    }

    #[test]
    fn matches_bcftools_1_24_reference_call() {
        let site = LikelihoodSite::new(
            0,
            0,
            Allele::new(&b"A"[..]).unwrap(),
            [Allele::new(&b"<*>"[..]).unwrap()],
            [1.0, 0.0],
            [SampleLikelihood::observed(
                Ploidy::new(2).unwrap(),
                [0, 3, 40],
                SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
            )
            .unwrap()],
        )
        .unwrap();

        let called = MultiallelicCaller::default().call(&site).unwrap();

        assert!(called.alternates().is_empty());
        assert_eq!(called.allele_counts(), &[2]);
        assert!(!called.is_variant());
        assert!((called.quality().unwrap() - 69.587).abs() < 1e-3);
        assert_eq!(called.samples()[0].genotype(), Some(&[0, 0][..]));
        assert_eq!(called.samples()[0].genotype_quality(), 0);
        assert_eq!(called.samples()[0].genotype_probabilities(), None);
        assert_eq!(called.samples()[0].phred_likelihoods(), None);
        assert_eq!(called.samples()[0].evidence().allele_quality_sums(), &[40]);
    }

    #[test]
    fn matches_bcftools_1_24_triallelic_call() {
        let make_sample = |likelihoods, depths, sums| {
            SampleLikelihood::observed(
                Ploidy::new(2).unwrap(),
                likelihoods,
                SampleEvidence::new(1, depths, sums).unwrap(),
            )
            .unwrap()
        };
        let site = LikelihoodSite::new(
            0,
            0,
            Allele::new(&b"A"[..]).unwrap(),
            [
                Allele::new(&b"G"[..]).unwrap(),
                Allele::new(&b"C"[..]).unwrap(),
                Allele::new(&b"<*>"[..]).unwrap(),
            ],
            [1.0, 1.0, 1.0, 0.0],
            [
                make_sample(
                    [0, 3, 40, 3, 40, 40, 3, 40, 40, 40],
                    [1, 0, 0, 0],
                    [40, 0, 0, 0],
                ),
                make_sample(
                    [40, 40, 40, 3, 3, 0, 40, 40, 3, 40],
                    [0, 0, 1, 0],
                    [0, 0, 40, 0],
                ),
                make_sample(
                    [40, 3, 0, 40, 3, 40, 40, 3, 40, 40],
                    [0, 1, 0, 0],
                    [0, 40, 0, 0],
                ),
            ],
        )
        .unwrap();

        let called = MultiallelicCaller::default().call(&site).unwrap();

        assert_eq!(
            called
                .alternates()
                .iter()
                .map(Allele::as_bytes)
                .collect::<Vec<_>>(),
            [b"G".as_slice(), b"C".as_slice()]
        );
        assert_eq!(called.allele_counts(), &[3, 2, 1]);
        assert!((called.quality().unwrap() - 15.6934).abs() < 1e-4);
        assert_eq!(called.samples()[0].genotype(), Some(&[0, 1][..]));
        assert_eq!(called.samples()[1].genotype(), Some(&[0, 2][..]));
        assert_eq!(called.samples()[2].genotype(), Some(&[0, 1][..]));
        assert_eq!(called.samples()[0].genotype_quality(), 1);
        assert_eq!(called.samples()[1].genotype_quality(), 1);
        assert_eq!(called.samples()[2].genotype_quality(), 1);
        assert_eq!(
            called.samples()[1].phred_likelihoods(),
            Some(&[40, 40, 40, 3, 3, 0][..])
        );
        let expected = [
            0.332_762,
            0.333_552,
            3.327_62e-5,
            0.333_552,
            6.655_24e-5,
            3.327_62e-5,
        ];
        for (&observed, expected) in called.samples()[0]
            .genotype_probabilities()
            .unwrap()
            .iter()
            .zip(expected)
        {
            assert!((observed - expected).abs() < 5e-7);
        }
    }
}
