use crate::{
    CallError, CallPloidy, CalledSample, CalledSite, LikelihoodSite, Result, SampleEvidence,
    SampleLikelihood,
};

const THETA: f64 = 1e-3;
const PHRED_SCALE: f64 = 4.343;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConsensusCallerConfig {
    reference_probability_threshold: f64,
}

impl ConsensusCallerConfig {
    pub fn new(reference_probability_threshold: f64) -> Result<Self> {
        if !reference_probability_threshold.is_finite()
            || !(0.0..=1.0).contains(&reference_probability_threshold)
        {
            return Err(CallError::InvalidConsensusThreshold);
        }
        Ok(Self {
            reference_probability_threshold,
        })
    }

    pub fn reference_probability_threshold(self) -> f64 {
        self.reference_probability_threshold
    }
}

impl Default for ConsensusCallerConfig {
    fn default() -> Self {
        Self {
            reference_probability_threshold: 0.5,
        }
    }
}

pub struct ConsensusCaller {
    config: ConsensusCallerConfig,
}

impl ConsensusCaller {
    pub fn new(config: ConsensusCallerConfig) -> Self {
        Self { config }
    }

    pub fn call(&self, site: &LikelihoodSite) -> Result<CalledSite> {
        self.call_with_ploidies(site, &vec![CallPloidy::Diploid; site.samples().len()])
    }

    pub fn call_with_ploidies(
        &self,
        site: &LikelihoodSite,
        ploidies: &[CallPloidy],
    ) -> Result<CalledSite> {
        if ploidies.len() != site.samples().len() {
            return Err(CallError::CallerPloidyCountMismatch);
        }
        let likelihoods = site
            .samples()
            .iter()
            .zip(ploidies)
            .map(|(sample, &ploidy)| ConsensusLikelihood::new(sample, ploidy))
            .collect::<Result<Vec<_>>>()?;
        let chromosome_count = ploidies
            .iter()
            .map(|ploidy| ploidy.chromosome_count())
            .sum::<usize>();
        if chromosome_count == 0 {
            return Err(CallError::InvalidPriorChromosomeCount);
        }
        let allele_count_likelihoods = allele_count_likelihoods(&likelihoods, chromosome_count)?;
        let posterior = posterior(&allele_count_likelihoods, chromosome_count)?;
        let reference_frequency = posterior
            .iter()
            .enumerate()
            .map(|(count, probability)| count as f64 * probability)
            .sum::<f64>()
            / chromosome_count as f64;
        let reference_probability = posterior[chromosome_count];
        let variant_probability = posterior[..chromosome_count].iter().sum::<f64>();
        let is_variant = reference_probability < self.config.reference_probability_threshold;
        let error_probability = if is_variant {
            reference_probability
        } else {
            variant_probability
        };
        let quality = Some(phred_quality(error_probability, 999.0) as f32);
        let output_allele_count = retained_allele_count(site, is_variant)?;
        let mut allele_counts = vec![0u32; output_allele_count];
        let samples = site
            .samples()
            .iter()
            .zip(likelihoods)
            .map(|(sample, likelihood)| {
                call_sample(
                    sample,
                    likelihood,
                    reference_frequency,
                    is_variant,
                    output_allele_count,
                    &mut allele_counts,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let alternates = site.alternates()[..output_allele_count - 1].to_vec();

        Ok(CalledSite {
            reference_sequence_id: site.reference_sequence_id(),
            position: site.position(),
            reference: site.reference().clone(),
            alternates: alternates.into(),
            quality,
            allele_counts: allele_counts.into(),
            samples: samples.into(),
        })
    }
}

impl Default for ConsensusCaller {
    fn default() -> Self {
        Self::new(ConsensusCallerConfig::default())
    }
}

#[derive(Clone, Copy)]
struct ConsensusLikelihood {
    ploidy: CallPloidy,
    alternate: f64,
    heterozygous: f64,
    reference: f64,
}

impl ConsensusLikelihood {
    fn new(sample: &SampleLikelihood, ploidy: CallPloidy) -> Result<Self> {
        if ploidy == CallPloidy::Absent {
            return Ok(Self {
                ploidy,
                alternate: 1.0,
                heterozygous: 0.0,
                reference: 1.0,
            });
        }
        let values = sample
            .phred_likelihoods()
            .ok_or(CallError::UnsupportedConsensusLikelihoods)?;
        let (reference, heterozygous, alternate) = match sample.ploidy().get() {
            2 if values.len() >= 3 => (
                phred_probability(values[0]),
                phred_probability(values[1]),
                phred_probability(values[2]),
            ),
            1 if ploidy == CallPloidy::Haploid && values.len() >= 2 => (
                phred_probability(values[0]),
                0.0,
                phred_probability(values[1]),
            ),
            _ => return Err(CallError::UnsupportedConsensusLikelihoods),
        };
        Ok(Self {
            ploidy,
            alternate,
            heterozygous,
            reference,
        })
    }
}

fn allele_count_likelihoods(
    samples: &[ConsensusLikelihood],
    chromosome_count: usize,
) -> Result<Vec<f64>> {
    let mut values = vec![1.0];
    let mut used_chromosomes = 0usize;
    for sample in samples {
        match sample.ploidy {
            CallPloidy::Absent => continue,
            CallPloidy::Haploid => {
                let mut next = vec![0.0; used_chromosomes + 2];
                for reference_count in 0..=used_chromosomes + 1 {
                    if reference_count <= used_chromosomes {
                        next[reference_count] += (used_chromosomes + 1 - reference_count) as f64
                            * sample.alternate
                            * values[reference_count];
                    }
                    if reference_count > 0 {
                        next[reference_count] +=
                            reference_count as f64 * sample.reference * values[reference_count - 1];
                    }
                }
                normalize(&mut next)?;
                values = next;
                used_chromosomes += 1;
            }
            CallPloidy::Diploid => {
                let mut next = vec![0.0; used_chromosomes + 3];
                for reference_count in 0..=used_chromosomes + 2 {
                    if reference_count <= used_chromosomes {
                        next[reference_count] += ((used_chromosomes - reference_count + 1)
                            * (used_chromosomes - reference_count + 2))
                            as f64
                            * sample.alternate
                            * values[reference_count];
                    }
                    if reference_count > 0 && reference_count - 1 <= used_chromosomes {
                        next[reference_count] +=
                            (reference_count * (used_chromosomes + 2 - reference_count)) as f64
                                * 2.0
                                * sample.heterozygous
                                * values[reference_count - 1];
                    }
                    if reference_count > 1 {
                        next[reference_count] += (reference_count * (reference_count - 1)) as f64
                            * sample.reference
                            * values[reference_count - 2];
                    }
                }
                normalize(&mut next)?;
                values = next;
                used_chromosomes += 2;
            }
        }
    }
    if used_chromosomes != chromosome_count {
        return Err(CallError::InvalidPriorChromosomeCount);
    }
    Ok(values)
}

fn posterior(likelihoods: &[f64], chromosome_count: usize) -> Result<Vec<f64>> {
    let mut prior = (0..chromosome_count)
        .map(|reference_count| THETA / (chromosome_count - reference_count) as f64)
        .collect::<Vec<_>>();
    prior.push(1.0 - prior.iter().sum::<f64>());
    let mut values = likelihoods
        .iter()
        .zip(prior)
        .map(|(likelihood, prior)| likelihood * prior)
        .collect::<Vec<_>>();
    normalize(&mut values)?;
    Ok(values)
}

fn normalize(values: &mut [f64]) -> Result<()> {
    let sum = values.iter().sum::<f64>();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(CallError::UnsupportedConsensusLikelihoods);
    }
    for value in values {
        *value /= sum;
    }
    Ok(())
}

fn retained_allele_count(site: &LikelihoodSite, is_variant: bool) -> Result<usize> {
    if !is_variant {
        return Ok(1);
    }
    let mut homozygous_sums = vec![0u64; site.alternates().len() + 1];
    for sample in site.samples() {
        let values = sample
            .phred_likelihoods()
            .ok_or(CallError::UnsupportedConsensusLikelihoods)?;
        for (allele, sum) in homozygous_sums.iter_mut().enumerate() {
            let index = match sample.ploidy().get() {
                2 => (allele + 1) * (allele + 2) / 2 - 1,
                1 => allele,
                _ => return Err(CallError::UnsupportedConsensusLikelihoods),
            };
            *sum = sum
                .checked_add(u64::from(values[index]))
                .ok_or(CallError::CalledAlleleCountOverflow)?;
        }
    }
    let mut ranked = homozygous_sums.into_iter().enumerate().collect::<Vec<_>>();
    ranked.sort_by_key(|&(allele, sum)| (sum, allele));
    let reference_rank = ranked.iter().position(|&(allele, _)| allele == 0).unwrap();
    Ok(if reference_rank < 2 {
        2
    } else {
        reference_rank + 1
    })
}

fn call_sample(
    sample: &SampleLikelihood,
    likelihood: ConsensusLikelihood,
    reference_frequency: f64,
    is_variant: bool,
    output_allele_count: usize,
    allele_counts: &mut [u32],
) -> Result<CalledSample> {
    let evidence = trim_evidence(sample.evidence(), output_allele_count)?;
    if likelihood.ploidy == CallPloidy::Absent {
        return Ok(CalledSample {
            ploidy: CallPloidy::Absent,
            genotype: None,
            genotype_quality: None,
            genotype_probabilities: None,
            phred_likelihoods: None,
            evidence,
        });
    }
    let (genotype, quality) = match likelihood.ploidy {
        CallPloidy::Diploid => {
            let alternate_frequency = 1.0 - reference_frequency;
            let mut probabilities = [
                likelihood.alternate * alternate_frequency * alternate_frequency,
                likelihood.heterozygous * 2.0 * reference_frequency * alternate_frequency,
                likelihood.reference * reference_frequency * reference_frequency,
            ];
            normalize(&mut probabilities)?;
            let index = if is_variant {
                maximum_index(&probabilities)
            } else {
                2
            };
            let genotype = match index {
                0 => vec![1, 1],
                1 => vec![0, 1],
                _ => vec![0, 0],
            };
            let quality = genotype_quality(1.0 - probabilities[index]);
            (genotype, Some(quality))
        }
        CallPloidy::Haploid => {
            let mut probabilities = [
                likelihood.alternate * (1.0 - reference_frequency),
                likelihood.reference * reference_frequency,
            ];
            normalize(&mut probabilities)?;
            let index = if is_variant {
                maximum_index(&probabilities)
            } else {
                1
            };
            (vec![usize::from(index == 0)], None)
        }
        CallPloidy::Absent => unreachable!(),
    };
    for &allele in &genotype {
        allele_counts[allele] = allele_counts[allele]
            .checked_add(1)
            .ok_or(CallError::CalledAlleleCountOverflow)?;
    }
    let phred_likelihoods = trim_likelihoods(sample, likelihood.ploidy, output_allele_count)?;
    Ok(CalledSample {
        ploidy: likelihood.ploidy,
        genotype: Some(genotype.into()),
        genotype_quality: quality,
        genotype_probabilities: None,
        phred_likelihoods: Some(phred_likelihoods.into()),
        evidence,
    })
}

fn trim_evidence(evidence: &SampleEvidence, allele_count: usize) -> Result<SampleEvidence> {
    SampleEvidence::new(
        evidence.depth(),
        evidence.allele_depths()[..allele_count].to_vec(),
        evidence.allele_quality_sums()[..allele_count].to_vec(),
    )
}

fn trim_likelihoods(
    sample: &SampleLikelihood,
    ploidy: CallPloidy,
    allele_count: usize,
) -> Result<Vec<u32>> {
    let values = sample
        .phred_likelihoods()
        .ok_or(CallError::UnsupportedConsensusLikelihoods)?;
    match ploidy {
        CallPloidy::Diploid => Ok(values[..allele_count * (allele_count + 1) / 2].to_vec()),
        CallPloidy::Haploid if sample.ploidy().get() == 1 => Ok(values[..allele_count].to_vec()),
        CallPloidy::Haploid => Ok((0..allele_count)
            .map(|allele| values[(allele + 1) * (allele + 2) / 2 - 1])
            .collect()),
        CallPloidy::Absent => Ok(Vec::new()),
    }
}

fn maximum_index(values: &[f64]) -> usize {
    values
        .iter()
        .enumerate()
        .fold((0, f64::NEG_INFINITY), |best, (index, &value)| {
            if value > best.1 { (index, value) } else { best }
        })
        .0
}

fn phred_probability(value: u32) -> f64 {
    10.0f64.powf(-f64::from(value) / 10.0)
}

fn phred_quality(error_probability: f64, maximum: f64) -> f64 {
    if error_probability < 1e-100 {
        maximum
    } else {
        (-PHRED_SCALE * error_probability.ln()).min(maximum)
    }
}

fn genotype_quality(error_probability: f64) -> u8 {
    if error_probability < 1e-308 {
        99
    } else {
        (-PHRED_SCALE * error_probability.ln() + 0.499).min(99.0) as u8
    }
}

#[cfg(test)]
mod tests;
