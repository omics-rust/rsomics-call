use crate::{
    CallError, CallPloidy, CalledAnnotations, CalledSample, CalledSite, LikelihoodSite, Result,
    SampleEvidence, SampleLikelihood,
};

const PHRED_SCALE: f64 = 4.342_94;
const SITE_PHRED_SCALE: f64 = 4.343;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiallelicCallerConfig {
    mutation_rate: f64,
    keep_alternates: bool,
}

impl MultiallelicCallerConfig {
    pub fn new(mutation_rate: f64) -> Result<Self> {
        if !mutation_rate.is_finite() || mutation_rate <= 0.0 || mutation_rate >= 1.0 {
            return Err(CallError::InvalidMutationRate);
        }
        Ok(Self {
            mutation_rate,
            keep_alternates: false,
        })
    }

    pub fn with_keep_alternates(mut self, enabled: bool) -> Self {
        self.keep_alternates = enabled;
        self
    }
}

impl Default for MultiallelicCallerConfig {
    fn default() -> Self {
        Self {
            mutation_rate: 1.1e-3,
            keep_alternates: false,
        }
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
            None,
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
            None,
        )
    }

    /// Groups select alleles independently while sharing the output allele record.
    pub fn call_with_groups(
        &self,
        site: &LikelihoodSite,
        ploidies: &[CallPloidy],
        prior_chromosome_count: usize,
        sample_groups: &[usize],
    ) -> Result<CalledSite> {
        if ploidies.len() != site.samples().len() {
            return Err(CallError::CallerPloidyCountMismatch);
        }
        if sample_groups.len() != site.samples().len() {
            return Err(CallError::CallerGroupCountMismatch);
        }
        self.call_inner(
            site,
            SamplePloidies::Explicit(ploidies),
            prior_chromosome_count,
            Some(sample_groups),
        )
    }

    fn call_inner(
        &self,
        site: &LikelihoodSite,
        ploidies: SamplePloidies<'_>,
        prior_chromosome_count: usize,
        sample_groups: Option<&[usize]>,
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
        let log_prior = adjusted_log_prior(self.config.mutation_rate, prior_chromosome_count);
        let groups = caller_groups(site, sample_groups, &likelihoods, &ploidies, log_prior)?;
        let mut retained = vec![false; allele_count];
        for group in &groups {
            for (retained, &selected) in retained.iter_mut().zip(&group.selection.alleles) {
                *retained |= selected;
            }
        }
        retained[0] = true;
        let has_variant = site.alternates().iter().enumerate().any(|(index, allele)| {
            retained[index + 1] && !matches!(allele.as_bytes(), b"<*>" | b"<NON_REF>")
        });
        if self.config.keep_alternates && has_variant {
            retained.fill(true);
        }
        for (index, allele) in site.alternates().iter().enumerate() {
            let remove = allele.as_bytes() == b"<*>"
                || !self.config.keep_alternates && allele.as_bytes() == b"<NON_REF>";
            if remove {
                retained[index + 1] = false;
            }
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
            let group_index = sample_groups.map_or(0, |groups| groups[sample_index]);
            let group = &groups[group_index];
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
                let context = SampleCallContext {
                    likelihoods: &likelihoods[sample_index],
                    quality_sums: &group.quality_sums,
                    retained: &retained,
                    allowed: &group.selection.alleles,
                    old_to_new: &old_to_new,
                    ploidy,
                };
                samples.push(call_sample(sample, context, &mut allele_counts)?);
            }
        }

        let alternate_count = allele_counts[1..].iter().sum::<u32>();
        let quality_group = groups.iter().fold(None, |best, group| {
            let Some(quality) = group.selection.maximum_quality else {
                return best;
            };
            match best {
                Some((_, best_quality)) if best_quality >= quality => best,
                _ => Some((group, quality)),
            }
        });
        let quality = if alternate_count != 0 {
            quality_group.map(|(_, quality)| quality)
        } else if let Some((group, _)) = quality_group
            && group.selection.alternate_log_sum.is_finite()
        {
            Some(
                (-SITE_PHRED_SCALE
                    * (group.selection.alternate_log_sum
                        - logsumexp(
                            group.selection.alternate_log_sum,
                            group.selection.reference_log_likelihood,
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
        let prior_allele_counts = site
            .prior_allele_counts()
            .map(|counts| counts.select(&retained));

        Ok(CalledSite {
            reference_sequence_id: site.reference_sequence_id(),
            position: site.position(),
            reference,
            alternates,
            quality,
            allele_counts: allele_counts.into_boxed_slice(),
            samples: samples.into_boxed_slice(),
            indel_summary: site.indel_summary(),
            annotations: site.annotations().map(CalledAnnotations::multiallelic),
            gvcf: None,
            prior_allele_counts,
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

#[derive(Clone, Copy)]
enum GroupSamples<'a> {
    All(usize),
    Explicit(&'a [usize]),
}

impl GroupSamples<'_> {
    fn for_each(self, visit: impl FnMut(usize)) {
        match self {
            Self::All(count) => (0..count).for_each(visit),
            Self::Explicit(indices) => indices.iter().copied().for_each(visit),
        }
    }
}

struct AlleleSelection {
    alleles: Vec<bool>,
    reference_log_likelihood: f64,
    alternate_log_sum: f64,
    maximum_quality: Option<f32>,
}

struct GroupCall {
    quality_sums: Vec<f64>,
    selection: AlleleSelection,
}

fn caller_groups(
    site: &LikelihoodSite,
    sample_groups: Option<&[usize]>,
    sample_likelihoods: &[Vec<f64>],
    ploidies: &SamplePloidies<'_>,
    log_prior: f64,
) -> Result<Vec<GroupCall>> {
    let Some(sample_groups) = sample_groups else {
        let quality_sums = normalized_quality_sums(
            site.allele_quality_sums()
                .iter()
                .map(|&value| f64::from(value)),
            site.prior_allele_counts(),
        );
        let selection = select_alleles(
            sample_likelihoods,
            &quality_sums,
            ploidies,
            GroupSamples::All(site.samples().len()),
            log_prior,
        );
        return Ok(vec![GroupCall {
            quality_sums,
            selection,
        }]);
    };
    let group_count = validate_sample_groups(sample_groups)?;
    let allele_count = site.alternates().len() + 1;
    let mut members = vec![Vec::new(); group_count];
    let mut quality_sums = vec![vec![0.0f32; allele_count]; group_count];
    for (sample_index, (&group, sample)) in sample_groups.iter().zip(site.samples()).enumerate() {
        members[group].push(sample_index);
        let values = sample.evidence().allele_quality_sums();
        let total = values.iter().map(|&value| value as f32).sum::<f32>();
        if total != 0.0 {
            for (sum, &value) in quality_sums[group].iter_mut().zip(values) {
                *sum += value as f32 / total;
            }
        }
    }
    Ok(quality_sums
        .into_iter()
        .zip(&members)
        .map(|(quality_sums, members)| {
            let quality_sums = normalized_quality_sums(
                quality_sums.into_iter().map(f64::from),
                site.prior_allele_counts(),
            );
            let selection = select_alleles(
                sample_likelihoods,
                &quality_sums,
                ploidies,
                GroupSamples::Explicit(members),
                log_prior,
            );
            GroupCall {
                quality_sums,
                selection,
            }
        })
        .collect())
}

pub(crate) fn validate_sample_groups(sample_groups: &[usize]) -> Result<usize> {
    let maximum_group = sample_groups
        .iter()
        .copied()
        .max()
        .ok_or(CallError::InvalidCallerGroups)?;
    if maximum_group >= sample_groups.len() {
        return Err(CallError::InvalidCallerGroups);
    }
    let group_count = maximum_group + 1;
    let mut seen = vec![false; group_count];
    for &group in sample_groups {
        seen[group] = true;
    }
    if seen.contains(&false) {
        return Err(CallError::InvalidCallerGroups);
    }
    Ok(group_count)
}

fn select_alleles(
    sample_likelihoods: &[Vec<f64>],
    quality_sums: &[f64],
    ploidies: &SamplePloidies<'_>,
    samples: GroupSamples<'_>,
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
        samples.for_each(|sample_index| {
            let likelihoods = &sample_likelihoods[sample_index];
            let value = likelihoods[index];
            if value != 0.0 {
                total += value.ln();
                supported = true;
            }
        });
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
            samples.for_each(|sample_index| {
                let likelihoods = &sample_likelihoods[sample_index];
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
            });
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
                samples.for_each(|sample_index| {
                    let likelihoods = &sample_likelihoods[sample_index];
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
                });
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

struct SampleCallContext<'a> {
    likelihoods: &'a [f64],
    quality_sums: &'a [f64],
    retained: &'a [usize],
    allowed: &'a [bool],
    old_to_new: &'a [Option<usize>],
    ploidy: CallPloidy,
}

fn call_sample(
    sample: &SampleLikelihood,
    context: SampleCallContext<'_>,
    allele_counts: &mut [u32],
) -> Result<CalledSample> {
    let (phred_likelihoods, evidence) =
        trimmed_fields(sample, context.retained, context.ploidy, true)?;
    if context.likelihoods.iter().all(|&value| value == 0.0) {
        return Ok(CalledSample {
            ploidy: context.ploidy,
            genotype: None,
            genotype_quality: None,
            genotype_probabilities: None,
            phred_likelihoods,
            evidence,
        });
    }

    let genotype_count = if context.ploidy == CallPloidy::Diploid {
        context.retained.len() * (context.retained.len() + 1) / 2
    } else {
        context.retained.len()
    };
    let mut probabilities = vec![0.0f32; genotype_count];
    let mut best = 0.0;
    let mut genotype = vec![0usize; context.ploidy.chromosome_count()];
    for &right in context.retained {
        if !context.allowed[right] {
            continue;
        }
        let new_right = context.old_to_new[right].unwrap();
        let homozygous = if context.ploidy == CallPloidy::Diploid {
            context.likelihoods[genotype_index(right, right)] * context.quality_sums[right].powi(2)
        } else {
            context.likelihoods[genotype_index(right, right)] * context.quality_sums[right]
        };
        let output_index = if context.ploidy == CallPloidy::Diploid {
            genotype_index(new_right, new_right)
        } else {
            new_right
        };
        probabilities[output_index] = homozygous as f32;
        if best < homozygous {
            best = homozygous;
            genotype.fill(new_right);
        }
        if context.ploidy == CallPloidy::Diploid {
            for &left in context
                .retained
                .iter()
                .take_while(|&&left| left != right)
                .filter(|&&left| context.allowed[left])
            {
                let new_left = context.old_to_new[left].unwrap();
                let heterozygous = 2.0
                    * context.likelihoods[genotype_index(right, left)]
                    * context.quality_sums[right]
                    * context.quality_sums[left];
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
        ploidy: context.ploidy,
        genotype: Some(genotype.into()),
        genotype_quality: Some(quality),
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
        genotype_quality: None,
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
        genotype_quality: Some(0),
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
    let evidence = sample.evidence().select(retained)?;
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

fn normalized_quality_sums(
    values: impl IntoIterator<Item = f64>,
    prior: Option<&crate::PriorAlleleCounts>,
) -> Vec<f64> {
    let mut values = values.into_iter().collect::<Vec<_>>();
    if let Some(prior) = prior
        && prior.total() != 0
    {
        let alternate_total = prior
            .alternates()
            .iter()
            .map(|&count| u64::from(count))
            .sum::<u64>();
        values[0] += 0.5 * (u64::from(prior.total()) - alternate_total) as f64;
        for (value, &count) in values[1..].iter_mut().zip(prior.alternates()) {
            *value += 0.5 * f64::from(count);
        }
    }
    let total = values.iter().sum::<f64>();
    if total == 0.0 {
        return vec![0.0; values.len()];
    }
    values.into_iter().map(|value| value / total).collect()
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
mod tests;
