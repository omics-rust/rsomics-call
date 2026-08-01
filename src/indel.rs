use std::ops::Range;

use rsomics_bamio::raw::RawRecord;
use rsomics_pileup::{Column, PileupRead};

use crate::{
    Allele, BaseObservation, CallError, ErrorModel, IndelSummary, LikelihoodMatrix, LikelihoodSite,
    Nucleotide, Ploidy, Result, SampleEvidence, SampleLikelihood,
    annotation::{AnnotationEvidence, AnnotationObservation, CigarMetrics, site_annotations},
    glocal,
};

const MAX_TYPES: usize = 64;
const MAX_SELECTED_TYPES: usize = 4;
const REVERSE: u16 = 0x10;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndelLikelihoodConfig {
    minimum_support: usize,
    minimum_fraction: f64,
    maximum_depth: usize,
    window_size: usize,
    gap_open_quality: i32,
    gap_extension_quality: i32,
    tandem_quality: i32,
    minimum_base_quality: u8,
    mapping_quality_cap: u8,
    indel_bias: f64,
    random_seed: i32,
    per_sample_support: bool,
    ambiguous_reads: IndelAmbiguousReadPolicy,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IndelAmbiguousReadPolicy {
    #[default]
    Drop,
    DistributeAlleleDepth,
    AddToReferenceAlleleDepth,
}

impl Default for IndelLikelihoodConfig {
    fn default() -> Self {
        Self {
            minimum_support: 2,
            minimum_fraction: 0.05,
            maximum_depth: 250,
            window_size: 110,
            gap_open_quality: 40,
            gap_extension_quality: 20,
            tandem_quality: 500,
            minimum_base_quality: 1,
            mapping_quality_cap: 60,
            indel_bias: 1.0,
            random_seed: 0,
            per_sample_support: false,
            ambiguous_reads: IndelAmbiguousReadPolicy::Drop,
        }
    }
}

impl IndelLikelihoodConfig {
    pub fn with_minimum_support(mut self, support: usize) -> Self {
        self.minimum_support = support;
        self
    }

    pub fn with_minimum_fraction(mut self, fraction: f64) -> Self {
        self.minimum_fraction = fraction;
        self
    }

    pub fn with_maximum_depth(mut self, depth: usize) -> Self {
        self.maximum_depth = depth;
        self
    }

    pub fn with_window_size(mut self, size: usize) -> Self {
        self.window_size = size;
        self
    }

    pub fn with_gap_open_quality(mut self, quality: i32) -> Self {
        self.gap_open_quality = quality;
        self
    }

    pub fn with_gap_extension_quality(mut self, quality: i32) -> Self {
        self.gap_extension_quality = quality;
        self
    }

    pub fn with_tandem_quality(mut self, quality: i32) -> Self {
        self.tandem_quality = quality;
        self
    }

    pub fn with_minimum_base_quality(mut self, quality: u8) -> Self {
        self.minimum_base_quality = quality;
        self
    }

    pub fn with_mapping_quality_cap(mut self, quality: u8) -> Self {
        self.mapping_quality_cap = quality;
        self
    }

    pub fn with_indel_bias(mut self, bias: f64) -> Self {
        self.indel_bias = bias;
        self
    }

    pub fn with_random_seed(mut self, seed: i32) -> Self {
        self.random_seed = seed;
        self
    }

    pub fn with_per_sample_support(mut self, enabled: bool) -> Self {
        self.per_sample_support = enabled;
        self
    }

    pub fn with_ambiguous_read_policy(mut self, policy: IndelAmbiguousReadPolicy) -> Self {
        self.ambiguous_reads = policy;
        self
    }
}

pub struct IndelSiteBuilder {
    sample_count: usize,
    config: IndelLikelihoodConfig,
    model: ErrorModel,
    observations: Vec<BaseObservation>,
}

impl IndelSiteBuilder {
    pub fn new(sample_count: usize, config: IndelLikelihoodConfig) -> Result<Self> {
        if sample_count == 0 {
            return Err(CallError::InvalidSampleCount);
        }
        if config.minimum_support == 0
            || !config.minimum_fraction.is_finite()
            || !(0.0..=1.0).contains(&config.minimum_fraction)
            || config.maximum_depth == 0
            || config.window_size == 0
            || config.gap_open_quality < 0
            || config.gap_extension_quality < 0
            || config.tandem_quality < 0
            || !config.indel_bias.is_finite()
            || config.indel_bias < 0.0
        {
            return Err(CallError::InvalidIndelConfig);
        }
        Ok(Self {
            sample_count,
            config,
            model: ErrorModel::with_random_seed(config.random_seed),
            observations: Vec::new(),
        })
    }

    pub fn build<F>(
        &mut self,
        column: &Column<'_>,
        reference_length: usize,
        mut sample_index: impl FnMut(u32, &RawRecord) -> Result<Option<usize>>,
        mut fetch_reference: F,
    ) -> Result<Option<LikelihoodSite>>
    where
        F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
    {
        let position =
            usize::try_from(column.position()).map_err(|_| CallError::InvalidPileupCoordinate)?;
        if position >= reference_length {
            return Ok(None);
        }
        if !column.entries().any(|entry| entry.projection().indel != 0) {
            return Ok(None);
        }

        let mut reads = Vec::with_capacity(column.len());
        for entry in column.entries() {
            let Some(sample) = sample_index(entry.source_id(), entry.record())? else {
                continue;
            };
            if sample >= self.sample_count {
                return Err(CallError::InvalidSampleIndex {
                    index: sample,
                    count: self.sample_count,
                });
            }
            let cigar_metrics = CigarMetrics::new(entry.cigar());
            let cigar = entry
                .cigar()
                .map(|(kind, length)| (kind, length as usize))
                .collect::<Vec<_>>();
            reads.push(Read {
                record: entry.record(),
                projection: entry.projection(),
                sample,
                query_len: query_len(&cigar)?,
                cigar,
                cigar_metrics,
                assigned_type: 0,
                sequence_quality: 0,
                indel_quality: 0,
            });
        }
        if reads.is_empty() || reads.len() >= self.config.maximum_depth {
            return Ok(None);
        }

        let Some(mut candidates) = Candidates::collect(&reads, self.sample_count, self.config)?
        else {
            return Ok(None);
        };
        let minimum_type = candidates.types[0];
        let maximum_insertion = candidates.types.last().copied().unwrap_or(0).max(0) as usize;
        let maximum_deletion = indel_length(minimum_type)?;
        let left = position.saturating_sub(self.config.window_size);
        let right = position
            .checked_add(self.config.window_size)
            .and_then(|value| value.checked_add(maximum_deletion))
            .unwrap_or(reference_length)
            .min(reference_length);
        let n_end = position
            .checked_add((self.config.window_size * 2).min(candidates.maximum_read_len))
            .unwrap_or(reference_length)
            .min(reference_length);
        let initial_end = right.max(n_end).max((position + 1).min(reference_length));
        let mut reference =
            ReferenceWindow::new(left, initial_end, reference_length, &mut fetch_reference)?;
        if too_many_ambiguous_reference_bases(
            &mut reference,
            position,
            n_end,
            &mut fetch_reference,
        )? {
            return Ok(None);
        }

        let sample_references =
            sample_references(&reads, self.sample_count, left, right, &mut reference)?;
        candidates.build_insertions(&reads);
        if candidates.types.len() < 2 {
            return Ok(None);
        }
        let reference_type = candidates
            .types
            .iter()
            .position(|&value| value == 0)
            .ok_or(CallError::IndelRealignment)?;
        let homopolymer_length =
            homopolymer_length(&mut reference, position, &mut fetch_reference)?;
        let max_context = maximum_insertion.max(maximum_deletion);
        let haplotype_capacity = right
            .saturating_sub(left)
            .saturating_add(2)
            .saturating_add(max_context.saturating_mul(2));
        let mut scores = vec![vec![i32::MAX; candidates.types.len()]; reads.len()];
        let mut indel_region = 0usize;

        for (type_index, &indel_type) in candidates.types.iter().enumerate() {
            indel_region = indel_region.max(indel_region_length(
                &mut reference,
                position,
                indel_type,
                candidates.insertion(type_index),
                &mut fetch_reference,
            )?);
            for (sample, sample_reference) in sample_references.iter().enumerate() {
                let haplotype = build_haplotype(
                    sample_reference,
                    left,
                    right,
                    position,
                    indel_type,
                    candidates.insertion(type_index),
                    haplotype_capacity,
                )?;
                for (read_index, read) in reads.iter().enumerate() {
                    if read.sample != sample || read.cigar.iter().any(|&(kind, _)| kind == 3) {
                        continue;
                    }
                    scores[read_index][type_index] = alignment_score(
                        read,
                        &haplotype,
                        left,
                        right,
                        position,
                        indel_type,
                        self.config,
                    )?;
                }
            }
        }

        let selected = assign_types(
            &mut reads,
            &scores,
            &candidates.types,
            reference_type,
            homopolymer_length,
            self.config,
        );
        if selected.len() < 2
            || reads
                .iter()
                .all(|read| read.assigned_type == 0 || read.assigned_type >= selected.len())
        {
            return Ok(None);
        }

        let alleles = build_alleles(
            &selected,
            &candidates,
            &mut reference,
            position,
            indel_region,
            &mut fetch_reference,
        )?;
        let (samples, annotation_evidence) =
            self.build_samples(&reads, self.sample_count, selected.len())?;
        let allele_quality_sums = normalized_quality_sums(&samples);
        let reference_sequence_id = usize::try_from(column.reference_id())
            .map_err(|_| CallError::InvalidPileupCoordinate)?;
        let position = u64::try_from(position).map_err(|_| CallError::InvalidPileupCoordinate)?;
        let mut alleles = alleles.into_iter();
        let reference = alleles.next().ok_or(CallError::IndelRealignment)?;
        let annotations = site_annotations(annotation_evidence.iter(), true, true)?;
        let site = LikelihoodSite::new(
            reference_sequence_id,
            position,
            reference,
            alleles.collect::<Vec<_>>(),
            allele_quality_sums,
            samples,
        )?
        .with_indel_summary(IndelSummary::new(
            candidates.maximum_support,
            candidates.maximum_fraction,
        )?)
        .with_annotations(annotations);
        Ok(Some(site))
    }

    fn build_samples(
        &mut self,
        reads: &[Read<'_>],
        sample_count: usize,
        allele_count: usize,
    ) -> Result<(Vec<SampleLikelihood>, Vec<AnnotationEvidence>)> {
        let diploid = Ploidy::new(2).unwrap();
        let mut samples = Vec::with_capacity(sample_count);
        let mut annotation_evidence = Vec::with_capacity(sample_count);
        for sample in 0..sample_count {
            self.observations.clear();
            let mut depth = 0u64;
            let mut allele_depths = vec![0u64; allele_count];
            let mut quality_sums = vec![0u64; allele_count];
            let mut forward_depths = vec![0u64; allele_count];
            let mut reverse_depths = vec![0u64; allele_count];
            let mut missed_forward = vec![0u64; allele_count];
            let mut missed_reverse = vec![0u64; allele_count];
            let mut annotations = AnnotationEvidence::default();
            let sample_depth = reads
                .iter()
                .filter(|read| read.sample == sample && !read.projection.is_reference_skip)
                .count();
            for read in reads.iter().filter(|read| read.sample == sample) {
                annotations.begin_read(read.cigar_metrics);
                annotations.observe_indel_candidate(
                    read.record,
                    read.projection,
                    read.cigar_metrics,
                );
                if read.projection.is_reference_skip {
                    continue;
                }
                annotations.add_raw_depth();
                let mut allele = read.assigned_type;
                let mut quality = read.indel_quality;
                let original_sequence_quality = read.sequence_quality;
                let mut sequence_quality = original_sequence_quality;
                if read.projection.indel == 0
                    && (usize::from(quality) < sample_depth / 2 || sample_depth > 20)
                {
                    allele = 0;
                    let base_quality = read
                        .record
                        .quality_scores()
                        .get(read.projection.qpos)
                        .copied()
                        .unwrap_or(255);
                    sequence_quality =
                        ((3 * u16::from(sequence_quality) + 2 * u16::from(base_quality)) / 8) as u8;
                }
                if sample_depth > 20 {
                    sequence_quality = sequence_quality.min(40);
                }
                if quality < self.config.minimum_base_quality {
                    if read.projection.indel == 0 && allele < allele_count {
                        if read.record.flags() & REVERSE == 0 {
                            missed_forward[allele] += 1;
                        } else {
                            missed_reverse[allele] += 1;
                        }
                    }
                    continue;
                }
                if allele >= allele_count {
                    continue;
                }
                depth += 1;
                quality = quality.min(sequence_quality);
                let raw_mapping_quality = match read.record.mapping_quality() {
                    255 => 20,
                    value => value,
                };
                let mapping_quality = raw_mapping_quality.min(self.config.mapping_quality_cap);
                quality = quality.min(mapping_quality).clamp(4, 63);
                allele_depths[allele] += 1;
                if read.record.flags() & REVERSE == 0 {
                    forward_depths[allele] += 1;
                } else {
                    reverse_depths[allele] += 1;
                }
                quality_sums[allele] += u64::from(quality);
                let observation = AnnotationObservation {
                    allele,
                    is_reference: allele == 0,
                    base_quality: original_sequence_quality,
                    raw_mapping_quality,
                    mapping_quality,
                    effective_quality: quality,
                };
                annotations.observe(read.record, read.projection, observation);
                annotations.observe_detailed(
                    read.record,
                    read.projection,
                    read.cigar_metrics,
                    observation,
                );
                self.observations.push(BaseObservation::new(
                    allele_nucleotide(allele),
                    quality,
                    read.record.flags() & REVERSE != 0,
                ));
            }
            compensate_ambiguous_depths(
                self.config.ambiguous_reads,
                &mut allele_depths,
                StrandDepths {
                    forward: &forward_depths,
                    reverse: &reverse_depths,
                },
                StrandDepths {
                    forward: &missed_forward,
                    reverse: &missed_reverse,
                },
                &mut annotations,
            );
            let matrix = self.model.calculate(&mut self.observations)?;
            let phred_likelihoods = likelihoods(&matrix, allele_count);
            let selected = (0..allele_count).collect::<Vec<_>>();
            let evidence = SampleEvidence::new(
                u32::try_from(depth).map_err(|_| CallError::IndelEvidenceOverflow)?,
                checked_counts(&allele_depths)?,
                checked_counts(&quality_sums)?,
            )?
            .with_annotations(annotations.sample_annotations(&selected)?)?;
            samples.push(SampleLikelihood::observed(
                diploid,
                phred_likelihoods,
                evidence,
            )?);
            annotation_evidence.push(annotations);
        }
        Ok((samples, annotation_evidence))
    }
}

struct StrandDepths<'a> {
    forward: &'a [u64],
    reverse: &'a [u64],
}

fn compensate_ambiguous_depths(
    policy: IndelAmbiguousReadPolicy,
    allele_depths: &mut [u64],
    observed: StrandDepths<'_>,
    missed: StrandDepths<'_>,
    annotations: &mut AnnotationEvidence,
) {
    if policy == IndelAmbiguousReadPolicy::Drop {
        return;
    }
    let missed_forward = missed.forward.iter().sum::<u64>();
    let missed_reverse = missed.reverse.iter().sum::<u64>();
    let forward_total = observed.forward.iter().sum::<u64>();
    let reverse_total = observed.reverse.iter().sum::<u64>();
    for (allele, depth) in allele_depths.iter_mut().enumerate() {
        let (forward, reverse) = match policy {
            IndelAmbiguousReadPolicy::DistributeAlleleDepth => (
                distribute_depth(observed.forward[allele], forward_total, missed_forward),
                distribute_depth(observed.reverse[allele], reverse_total, missed_reverse),
            ),
            IndelAmbiguousReadPolicy::AddToReferenceAlleleDepth if allele == 0 => {
                (missed_forward, missed_reverse)
            }
            IndelAmbiguousReadPolicy::AddToReferenceAlleleDepth => (0, 0),
            IndelAmbiguousReadPolicy::Drop => unreachable!(),
        };
        *depth += forward + reverse;
        if forward != 0 || reverse != 0 {
            annotations.add_allele_depth(allele, forward, reverse);
        }
    }
}

fn distribute_depth(observed: u64, observed_total: u64, missed_total: u64) -> u64 {
    if observed_total == 0 || missed_total == 0 {
        return 0;
    }
    (missed_total as f32 * observed as f32 / observed_total as f32).round() as u64
}

struct Read<'a> {
    record: &'a RawRecord,
    projection: &'a PileupRead,
    sample: usize,
    cigar: Vec<(u8, usize)>,
    cigar_metrics: CigarMetrics,
    query_len: usize,
    assigned_type: usize,
    sequence_quality: u8,
    indel_quality: u8,
}

struct Candidates {
    types: Vec<i64>,
    insertions: Vec<Vec<u8>>,
    maximum_read_len: usize,
    maximum_support: u32,
    maximum_fraction: f32,
}

impl Candidates {
    fn collect(
        reads: &[Read<'_>],
        sample_count: usize,
        config: IndelLikelihoodConfig,
    ) -> Result<Option<Self>> {
        let mut types = vec![0i64];
        let mut sample_depths = vec![0usize; sample_count];
        let mut sample_support = vec![0usize; sample_count];
        let mut maximum_read_len = 0usize;
        for read in reads {
            sample_depths[read.sample] += 1;
            maximum_read_len = maximum_read_len.max(read.query_len);
            if read.projection.indel != 0 {
                sample_support[read.sample] += 1;
                types.push(read.projection.indel);
            }
        }
        let total_depth = sample_depths.iter().sum::<usize>();
        let total_support = sample_support.iter().sum::<usize>();
        let support_passes = if config.per_sample_support {
            sample_support
                .iter()
                .zip(&sample_depths)
                .any(|(&support, &depth)| {
                    depth != 0
                        && support >= config.minimum_support
                        && support as f64 / depth as f64 >= config.minimum_fraction
                })
        } else {
            total_depth != 0
                && total_support >= config.minimum_support
                && total_support as f64 / total_depth as f64 >= config.minimum_fraction
        };
        if !support_passes {
            return Ok(None);
        }
        types.sort_unstable();
        types.dedup();
        if types.len() < 2 || types.len() >= MAX_TYPES {
            return Ok(None);
        }
        let insertions = vec![Vec::new(); types.len()];
        let mut maximum_support = 0usize;
        let mut maximum_fraction = 0.0f64;
        for (&support, &depth) in sample_support.iter().zip(&sample_depths) {
            let fraction = if depth == 0 {
                0.0
            } else {
                support as f64 / depth as f64
            };
            if support > maximum_support && fraction > 0.0 {
                maximum_support = support;
                maximum_fraction = fraction;
            }
        }
        Ok(Some(Self {
            types,
            insertions,
            maximum_read_len,
            maximum_support: u32::try_from(maximum_support)
                .map_err(|_| CallError::IndelEvidenceOverflow)?,
            maximum_fraction: maximum_fraction as f32,
        }))
    }

    fn build_insertions(&mut self, reads: &[Read<'_>]) {
        let mut retained_types = Vec::with_capacity(self.types.len());
        let mut retained_insertions = Vec::with_capacity(self.types.len());
        for &indel_type in &self.types {
            if indel_type <= 0 {
                retained_types.push(indel_type);
                retained_insertions.push(Vec::new());
                continue;
            }
            let length = indel_type as usize;
            let mut consensus = Vec::with_capacity(length);
            for offset in 1..=length {
                let mut counts = [0usize; 5];
                for read in reads
                    .iter()
                    .filter(|read| read.projection.indel == indel_type)
                {
                    let index = read.projection.qpos + offset;
                    let base = read
                        .record
                        .sequence_len()
                        .checked_sub(index + 1)
                        .map(|_| decode_base(read.record.seq_nibble(index)))
                        .unwrap_or(4);
                    counts[usize::from(base)] += 1;
                }
                let base = counts
                    .iter()
                    .enumerate()
                    .max_by_key(|&(index, count)| (*count, std::cmp::Reverse(index)))
                    .map(|(index, _)| index as u8)
                    .unwrap_or(4);
                consensus.push(base);
            }
            if consensus.iter().all(|&base| base < 4) {
                retained_types.push(indel_type);
                retained_insertions.push(consensus);
            }
        }
        self.types = retained_types;
        self.insertions = retained_insertions;
    }

    fn insertion(&self, type_index: usize) -> &[u8] {
        &self.insertions[type_index]
    }
}

struct ReferenceWindow {
    start: usize,
    length: usize,
    bases: Vec<u8>,
}

impl ReferenceWindow {
    fn new<F>(start: usize, end: usize, length: usize, fetch: &mut F) -> Result<Self>
    where
        F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
    {
        let mut bases = Vec::with_capacity(end.saturating_sub(start));
        fetch(start..end, &mut bases)?;
        if bases.len() != end.saturating_sub(start) {
            return Err(CallError::InvalidIndelReference);
        }
        Ok(Self {
            start,
            length,
            bases,
        })
    }

    fn end(&self) -> usize {
        self.start + self.bases.len()
    }

    fn ensure<F>(&mut self, range: Range<usize>, fetch: &mut F) -> Result<()>
    where
        F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
    {
        if range.end > self.length || range.start > range.end {
            return Err(CallError::InvalidIndelReference);
        }
        while range.start < self.start {
            let new_start = range.start.max(self.start.saturating_sub(64 * 1024));
            let mut prefix = Vec::with_capacity(self.start - new_start + self.bases.len());
            fetch(new_start..self.start, &mut prefix)?;
            if prefix.len() != self.start - new_start {
                return Err(CallError::InvalidIndelReference);
            }
            prefix.extend_from_slice(&self.bases);
            self.bases = prefix;
            self.start = new_start;
        }
        while range.end > self.end() {
            let old_end = self.end();
            let new_end = range.end.min(old_end.saturating_add(64 * 1024));
            fetch(old_end..new_end, &mut self.bases)?;
            if self.end() != new_end {
                return Err(CallError::InvalidIndelReference);
            }
        }
        Ok(())
    }

    fn base<F>(&mut self, position: usize, fetch: &mut F) -> Result<u8>
    where
        F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
    {
        let end = position
            .checked_add(1)
            .ok_or(CallError::InvalidIndelReference)?;
        self.ensure(position..end, fetch)?;
        Ok(self.bases[position - self.start])
    }

    fn slice<F>(&mut self, range: Range<usize>, fetch: &mut F) -> Result<&[u8]>
    where
        F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
    {
        self.ensure(range.clone(), fetch)?;
        Ok(&self.bases[range.start - self.start..range.end - self.start])
    }
}

fn too_many_ambiguous_reference_bases<F>(
    reference: &mut ReferenceWindow,
    position: usize,
    end: usize,
    fetch: &mut F,
) -> Result<bool>
where
    F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
{
    if position >= end {
        return Ok(false);
    }
    let sequence = reference.slice(position..end, fetch)?;
    Ok(sequence.iter().filter(|&&base| base == b'N').count() * 2 > sequence.len())
}

fn sample_references(
    reads: &[Read<'_>],
    sample_count: usize,
    left: usize,
    right: usize,
    reference: &mut ReferenceWindow,
) -> Result<Vec<Vec<u8>>> {
    let template = reference
        .bases
        .get(left - reference.start..right - reference.start)
        .ok_or(CallError::InvalidIndelReference)?
        .iter()
        .copied()
        .map(decode_reference_base)
        .collect::<Vec<_>>();
    (0..sample_count)
        .map(|sample| {
            let mut counts = vec![[0u32; 2]; right.saturating_sub(left)];
            for read in reads.iter().filter(|read| read.sample == sample) {
                let mut reference_position =
                    usize::try_from(read.record.alignment_start()).unwrap_or(usize::MAX);
                let mut query_position = 0usize;
                for &(kind, length) in &read.cigar {
                    match kind {
                        0 | 7 | 8 => {
                            if reference_position.saturating_add(length) >= left {
                                let begin = left.saturating_sub(reference_position);
                                let end = length.min(right.saturating_sub(reference_position));
                                for offset in begin..end {
                                    let index = reference_position + offset - left;
                                    let query = query_position + offset;
                                    if query >= read.record.sequence_len() {
                                        break;
                                    }
                                    let base = decode_base(read.record.seq_nibble(query));
                                    if base == template[index] {
                                        counts[index][0] += 1;
                                    } else {
                                        counts[index][1] += 1;
                                    }
                                }
                            }
                            reference_position = reference_position.saturating_add(length);
                            query_position = query_position.saturating_add(length);
                        }
                        2 | 3 => reference_position = reference_position.saturating_add(length),
                        1 | 4 => query_position = query_position.saturating_add(length),
                        _ => {}
                    }
                    if reference_position > right {
                        break;
                    }
                }
            }
            let mut consensus = template.clone();
            let mut deepest = (0u32, 0u32, -1isize);
            let mut second = (0u32, 0u32, -1isize);
            for (index, &[reference_count, alternate_count]) in counts.iter().enumerate() {
                if alternate_count >= deepest.0 {
                    second = deepest;
                    deepest = (alternate_count, reference_count, index as isize);
                } else if alternate_count >= second.0 {
                    second = (alternate_count, reference_count, index as isize);
                }
            }
            for (alternate_count, reference_count, index) in [deepest, second] {
                if index < 0 {
                    continue;
                }
                let total = reference_count + alternate_count;
                if total == 0 || f64::from(reference_count) / f64::from(total) < 0.7 {
                    consensus[index as usize] = 4;
                }
            }
            consensus.push(0);
            Ok(consensus)
        })
        .collect()
}

fn homopolymer_length<F>(
    reference: &mut ReferenceWindow,
    position: usize,
    fetch: &mut F,
) -> Result<usize>
where
    F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
{
    let Some(next) = position
        .checked_add(1)
        .filter(|&value| value < reference.length)
    else {
        return Ok(1);
    };
    let base = decode_reference_base(reference.base(next, fetch)?);
    if base == 4 {
        return Ok(1);
    }
    let mut right = next + 1;
    while right < reference.length && decode_reference_base(reference.base(right, fetch)?) == base {
        right += 1;
    }
    let mut left = position;
    while decode_reference_base(reference.base(left, fetch)?) == base {
        if left == 0 {
            return Ok(right);
        }
        left -= 1;
    }
    Ok(right - left - 1)
}

fn build_haplotype(
    sample_reference: &[u8],
    left: usize,
    right: usize,
    position: usize,
    indel_type: i64,
    insertion: &[u8],
    capacity: usize,
) -> Result<Vec<u8>> {
    let anchor_end = position - left + 1;
    let mut haplotype = Vec::with_capacity(capacity);
    haplotype.extend_from_slice(&sample_reference[..anchor_end]);
    let mut reference_position = position + 1;
    if indel_type < 0 {
        reference_position = reference_position.saturating_add(indel_length(indel_type)?);
    } else if indel_type > 0 {
        haplotype.extend_from_slice(insertion);
    }
    if reference_position < right {
        haplotype.extend_from_slice(&sample_reference[reference_position - left..right - left]);
    }
    haplotype.resize(capacity, 4);
    Ok(haplotype)
}

fn alignment_score(
    read: &Read<'_>,
    haplotype: &[u8],
    left: usize,
    right: usize,
    position: usize,
    indel_type: i64,
    config: IndelLikelihoodConfig,
) -> Result<i32> {
    let mut left_context = left;
    let mut right_context = right;
    if read.record.sequence_len() > 1000 {
        if position.saturating_sub(left) >= config.window_size {
            left_context = left_context.saturating_add(config.window_size / 2);
        }
        if right.saturating_sub(position) >= config.window_size {
            right_context = right_context.saturating_sub(config.window_size / 2);
        }
    }
    let alignment_start =
        usize::try_from(read.record.alignment_start()).map_err(|_| CallError::IndelRealignment)?;
    let alignment_span = reference_span(&read.cigar)?;
    let alignment_end = alignment_start
        .checked_add(alignment_span.saturating_sub(1))
        .ok_or(CallError::IndelRealignment)?;
    let (query_begin, mut target_begin) =
        target_to_query(alignment_start, &read.cigar, left_context, false)?;
    let (query_at_position, _) = target_to_query(alignment_start, &read.cigar, position, false)?;
    let (query_end, target_end) =
        target_to_query(alignment_start, &read.cigar, right_context, true)?;
    if indel_type < 0 {
        target_begin = target_begin
            .saturating_sub(indel_length(indel_type)?)
            .max(left_context);
    }
    if target_end <= target_begin
        || query_begin >= query_end
        || query_end > read.record.sequence_len()
    {
        return Ok(0x00ff_ffff);
    }
    let reference_begin = target_begin.saturating_sub(left);
    let reference_len = target_end
        .saturating_sub(target_begin)
        .checked_add(indel_length(indel_type)?)
        .ok_or(CallError::IndelRealignment)?;
    let reference_end = reference_begin
        .checked_add(reference_len)
        .ok_or(CallError::IndelRealignment)?;
    let reference = haplotype
        .get(reference_begin..reference_end)
        .ok_or(CallError::IndelRealignment)?;
    let query = (query_begin..query_end)
        .map(|index| decode_base(read.record.seq_nibble(index)))
        .collect::<Vec<_>>();
    let qualities = original_qualities(read.record, query_begin..query_end)?;
    let long_read = read.record.sequence_len() > 1000;
    let gap_open = if long_read { 1e-3 } else { 1e-4 };
    let gap_extension = if long_read { 1e-1 } else { 1e-2 };
    let bandwidth = indel_length(indel_type)?.saturating_add(3);
    let Some(raw_score) = glocal::score(
        reference,
        &query,
        &qualities,
        gap_open,
        gap_extension,
        bandwidth,
    ) else {
        return Ok(0x00ff_ffff);
    };
    let mut normalized = (((100.0 * f64::from(raw_score) / query.len() as f64 + 0.499) as i32)
        as f64
        * config.indel_bias) as i32;
    normalized = normalized.max(0);
    let query_position = query_at_position.saturating_sub(query_begin);
    let repeat_score = repeat_score(
        reference,
        query_position,
        alignment_start,
        alignment_end,
        target_begin,
    );
    normalized = (f64::from(normalized) * 0.8 + f64::from(repeat_score * 2)) as i32;
    Ok(raw_score
        .saturating_mul(256)
        .saturating_add(normalized.min(255)))
}

fn assign_types(
    reads: &mut [Read<'_>],
    scores: &[Vec<i32>],
    types: &[i64],
    reference_type: usize,
    homopolymer_length: usize,
    config: IndelLikelihoodConfig,
) -> Vec<usize> {
    let mut quality_sums = vec![0u64; types.len()];
    for (read, read_scores) in reads.iter_mut().zip(scores) {
        let mut order = (0..types.len()).collect::<Vec<_>>();
        order.sort_unstable_by_key(|&index| {
            (i64::from(read_scores[index]) << 6).saturating_add(index as i64)
        });
        let best = order[0];
        let comparison = if best == reference_type {
            order.get(1).copied().unwrap_or(reference_type)
        } else {
            reference_type
        };
        let best_score = read_scores[best] >> 8;
        let comparison_score = read_scores[comparison] >> 8;
        let mut indel_quality = comparison_score.saturating_sub(best_score);
        let sequence_type = if best == reference_type {
            types[comparison]
        } else {
            types[best]
        };
        let sequence_quality = sequence_quality(sequence_type, homopolymer_length, config);
        let normalized = read_scores[best] & 0xff;
        indel_quality = if normalized > 111 {
            0
        } else {
            ((1.0 - f64::from(normalized) / 111.0) * f64::from(indel_quality) + 0.499) as i32
        };
        indel_quality = indel_quality.min(sequence_quality).clamp(0, 255);
        let sequence_quality = sequence_quality.clamp(0, 255);
        read.assigned_type = best;
        read.indel_quality = indel_quality as u8;
        read.sequence_quality = sequence_quality as u8;
        quality_sums[best] =
            quality_sums[best].saturating_add(indel_quality.min(sequence_quality) as u64);
    }

    let mut selected = (0..types.len()).collect::<Vec<_>>();
    selected.sort_unstable_by(|&left, &right| {
        quality_sums[right]
            .cmp(&quality_sums[left])
            .then_with(|| left.cmp(&right))
    });
    selected.retain(|&index| index != reference_type);
    selected.insert(0, reference_type);
    selected.truncate(MAX_SELECTED_TYPES);
    for read in reads {
        read.assigned_type = selected
            .iter()
            .position(|&index| index == read.assigned_type)
            .unwrap_or(MAX_SELECTED_TYPES);
    }
    selected
}

fn sequence_quality(
    indel_type: i64,
    homopolymer_length: usize,
    config: IndelLikelihoodConfig,
) -> i32 {
    let length = indel_type.unsigned_abs() as i32;
    let gap = config.gap_open_quality
        + config
            .gap_extension_quality
            .saturating_mul(length.saturating_sub(1));
    let tandem = if homopolymer_length >= 3 {
        (f64::from(config.tandem_quality) * f64::from(length) / homopolymer_length as f64 + 0.499)
            as i32
    } else {
        1000
    };
    gap.min(tandem)
}

fn build_alleles<F>(
    selected: &[usize],
    candidates: &Candidates,
    reference: &mut ReferenceWindow,
    position: usize,
    indel_region: usize,
    fetch: &mut F,
) -> Result<Vec<Allele>>
where
    F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
{
    let anchor = reference.base(position, fetch)?.to_ascii_uppercase();
    let suffix = reference
        .slice(position + 1..position + 1 + indel_region, fetch)?
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect::<Vec<_>>();
    selected
        .iter()
        .map(|&index| {
            let indel_type = candidates.types[index];
            let mut allele = vec![anchor];
            if indel_type < 0 {
                let deletion = indel_length(indel_type)?;
                allele.extend_from_slice(
                    suffix
                        .get(deletion..)
                        .ok_or(CallError::InvalidIndelReference)?,
                );
            } else if indel_type > 0 {
                allele.extend(
                    candidates
                        .insertion(index)
                        .iter()
                        .map(|&base| b"ACGTN"[usize::from(base)]),
                );
                allele.extend_from_slice(&suffix);
            } else {
                allele.extend_from_slice(&suffix);
            }
            Allele::new(allele)
        })
        .collect()
}

fn indel_region_length<F>(
    reference: &mut ReferenceWindow,
    position: usize,
    indel_type: i64,
    insertion: &[u8],
    fetch: &mut F,
) -> Result<usize>
where
    F: FnMut(Range<usize>, &mut Vec<u8>) -> Result<()>,
{
    if indel_type == 0 {
        return Ok(0);
    }
    let length = indel_length(indel_type)?;
    let mut score = 0i32;
    let mut maximum = 0i32;
    let mut maximum_position = position;
    let mut cursor = position + 1;
    while cursor < reference.length {
        let observed = reference.base(cursor, fetch)?.to_ascii_uppercase();
        let expected = if indel_type > 0 {
            b"ACGTN"[usize::from(insertion[(cursor - position - 1) % length])]
        } else {
            reference
                .base(position + 1 + (cursor - position - 1) % length, fetch)?
                .to_ascii_uppercase()
        };
        score += if observed == expected { 1 } else { -10 };
        if score < 0 {
            break;
        }
        if score > maximum {
            maximum = score;
            maximum_position = cursor;
        }
        cursor += 1;
    }
    Ok(maximum_position - position)
}

fn target_to_query(
    alignment_start: usize,
    cigar: &[(u8, usize)],
    target: usize,
    left_edge: bool,
) -> Result<(usize, usize)> {
    let mut reference = alignment_start;
    let mut query = 0usize;
    let mut last_query = 0usize;
    for &(kind, length) in cigar {
        match kind {
            0 | 7 | 8 => {
                if alignment_start > target {
                    return Ok((query, alignment_start));
                }
                if reference.saturating_add(length) > target {
                    return Ok((query + target - reference, target));
                }
                reference = reference
                    .checked_add(length)
                    .ok_or(CallError::IndelRealignment)?;
                query = query
                    .checked_add(length)
                    .ok_or(CallError::IndelRealignment)?;
                last_query = query;
            }
            1 | 4 => {
                query = query
                    .checked_add(length)
                    .ok_or(CallError::IndelRealignment)?;
            }
            2 | 3 => {
                if reference.saturating_add(length) > target {
                    return Ok((
                        query,
                        if left_edge {
                            reference
                        } else {
                            reference + length
                        },
                    ));
                }
                reference = reference
                    .checked_add(length)
                    .ok_or(CallError::IndelRealignment)?;
            }
            _ => {}
        }
    }
    Ok((last_query, reference))
}

fn original_qualities(record: &RawRecord, range: Range<usize>) -> Result<Vec<u8>> {
    let qualities = record.quality_scores();
    let adjustment = (record.aux_type(*b"ZQ") == Some(b'Z'))
        .then(|| record.aux_value(*b"ZQ"))
        .flatten();
    range
        .map(|index| {
            let quality = qualities.get(index).copied().unwrap_or(255);
            let quality = if let Some(adjustment) = adjustment {
                let encoded = adjustment
                    .get(index)
                    .copied()
                    .ok_or(CallError::IndelRealignment)?;
                i16::from(quality) + i16::from(encoded) - 64
            } else {
                i16::from(quality)
            };
            Ok(quality.clamp(7, 30) as u8)
        })
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Repeat {
    start: usize,
    end: usize,
    length: usize,
}

fn repeat_score(
    sequence: &[u8],
    position: usize,
    alignment_start: usize,
    alignment_end: usize,
    target_begin: usize,
) -> i32 {
    repeats(sequence)
        .into_iter()
        .filter(|repeat| repeat.start <= position && repeat.end >= position)
        .map(|repeat| {
            let span = repeat.end - repeat.start;
            let mut score = (span / repeat.length) as i32;
            if repeat.start + target_begin <= alignment_start
                || repeat.end + target_begin >= alignment_end
            {
                score += (span * 2) as i32;
            }
            score
        })
        .sum()
}

fn repeats(sequence: &[u8]) -> Vec<Repeat> {
    let mut repeats = Vec::new();
    let mut word = 0u32;
    let split = sequence.len().min(15);
    for (position, &base) in sequence[..split].iter().enumerate() {
        word = word.wrapping_shl(2) | u32::from(base);
        for length in 1..=7 {
            if position >= length * 2 - 1 && repeated_word(word, length) {
                add_repeat(&mut repeats, sequence, position, length);
            }
        }
    }
    for (position, &base) in sequence.iter().enumerate().skip(split) {
        word = word.wrapping_shl(2) | u32::from(base);
        if let Some(length) = (1..=8).rev().find(|&length| repeated_word(word, length)) {
            add_repeat(&mut repeats, sequence, position, length);
        }
    }
    repeats
}

fn repeated_word(word: u32, length: usize) -> bool {
    let bits = length * 2;
    let mask = (1u32 << bits) - 1;
    word & mask == word >> bits & mask
}

fn add_repeat(repeats: &mut Vec<Repeat>, sequence: &[u8], position: usize, length: usize) {
    let start = position + 1 - length * 2;
    if repeats
        .last()
        .is_some_and(|repeat| repeat.start <= start && repeat.end >= position)
    {
        return;
    }
    let mut left = position + 1 - length;
    let mut right = position + 1;
    while right < sequence.len() && sequence[left] == sequence[right] {
        left += 1;
        right += 1;
    }
    let repeat = Repeat {
        start,
        end: position + right - (position + 1),
        length,
    };
    repeats.retain(|existing| existing.end < repeat.start || existing.start < repeat.start);
    repeats.push(repeat);
}

fn query_len(cigar: &[(u8, usize)]) -> Result<usize> {
    cigar
        .iter()
        .filter(|&&(kind, _)| matches!(kind, 0 | 1 | 4 | 7 | 8))
        .try_fold(0usize, |value, &(_, length)| {
            value.checked_add(length).ok_or(CallError::IndelRealignment)
        })
}

fn indel_length(indel_type: i64) -> Result<usize> {
    usize::try_from(indel_type.unsigned_abs()).map_err(|_| CallError::IndelRealignment)
}

fn reference_span(cigar: &[(u8, usize)]) -> Result<usize> {
    cigar
        .iter()
        .filter(|&&(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .try_fold(0usize, |value, &(_, length)| {
            value.checked_add(length).ok_or(CallError::IndelRealignment)
        })
}

fn likelihoods(matrix: &LikelihoodMatrix, allele_count: usize) -> Vec<u32> {
    let mut values = Vec::with_capacity(allele_count * (allele_count + 1) / 2);
    for second in 0..allele_count {
        for first in 0..=second {
            values.push(matrix.get(allele_nucleotide(first), allele_nucleotide(second)));
        }
    }
    let minimum = values.iter().copied().fold(f32::INFINITY, f32::min);
    values
        .into_iter()
        .map(|value| ((f64::from(value - minimum) + 0.499) as u32).min(255))
        .collect()
}

fn normalized_quality_sums(samples: &[SampleLikelihood]) -> Vec<f32> {
    let allele_count = samples
        .first()
        .map(|sample| sample.evidence().allele_quality_sums().len())
        .unwrap_or(0);
    let mut sums = vec![0.0f32; allele_count];
    for sample in samples {
        let qualities = sample.evidence().allele_quality_sums();
        let total = qualities.iter().map(|&value| u64::from(value)).sum::<u64>();
        if total != 0 {
            for (sum, &quality) in sums.iter_mut().zip(qualities) {
                *sum += quality as f32 / total as f32;
            }
        }
    }
    sums
}

fn checked_counts(values: &[u64]) -> Result<Vec<u32>> {
    values
        .iter()
        .map(|&value| u32::try_from(value).map_err(|_| CallError::IndelEvidenceOverflow))
        .collect()
}

fn allele_nucleotide(index: usize) -> Nucleotide {
    match index {
        0 => Nucleotide::A,
        1 => Nucleotide::C,
        2 => Nucleotide::G,
        3 => Nucleotide::T,
        _ => Nucleotide::N,
    }
}

fn decode_base(base: u8) -> u8 {
    match base {
        1 => 0,
        2 => 1,
        4 => 2,
        8 => 3,
        _ => 4,
    }
}

fn decode_reference_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeats_match_bcftools_1_24_str_finder() {
        type Case<'a> = (&'a [u8], &'a [(usize, usize, usize)]);
        let cases: &[Case<'_>] = &[
            (b"AAAAAA", &[(0, 5, 1)]),
            (b"ACACACAC", &[(0, 7, 2)]),
            (
                b"AGGGGAGGAGAAGAC",
                &[(1, 4, 1), (3, 9, 3), (7, 10, 2), (8, 13, 3)],
            ),
            (b"ACGTACGTACGT", &[(0, 11, 4)]),
            (
                b"TATATACCCCGGGGTTTT",
                &[(0, 5, 2), (6, 9, 1), (10, 13, 1), (14, 17, 1)],
            ),
            (b"ACGTACGTACGTACGTACGT", &[(0, 19, 4)]),
            (b"NNNN", &[]),
        ];
        for &(sequence, expected) in cases {
            let sequence = sequence
                .iter()
                .map(|base| decode_reference_base(*base))
                .collect::<Vec<_>>();
            let expected = expected
                .iter()
                .map(|&(start, end, length)| Repeat { start, end, length })
                .collect::<Vec<_>>();
            assert_eq!(repeats(&sequence), expected);
        }
    }

    #[test]
    fn validates_indel_configuration_and_sample_count() {
        assert!(matches!(
            IndelSiteBuilder::new(0, IndelLikelihoodConfig::default()),
            Err(CallError::InvalidSampleCount)
        ));
        for config in [
            IndelLikelihoodConfig::default().with_minimum_support(0),
            IndelLikelihoodConfig::default().with_minimum_fraction(f64::NAN),
            IndelLikelihoodConfig::default().with_minimum_fraction(1.1),
            IndelLikelihoodConfig::default().with_maximum_depth(0),
            IndelLikelihoodConfig::default().with_window_size(0),
            IndelLikelihoodConfig::default().with_gap_open_quality(-1),
            IndelLikelihoodConfig::default().with_gap_extension_quality(-1),
            IndelLikelihoodConfig::default().with_tandem_quality(-1),
            IndelLikelihoodConfig::default().with_indel_bias(-1.0),
        ] {
            assert!(matches!(
                IndelSiteBuilder::new(1, config),
                Err(CallError::InvalidIndelConfig)
            ));
        }
        assert!(
            IndelSiteBuilder::new(
                1,
                IndelLikelihoodConfig::default()
                    .with_minimum_base_quality(2)
                    .with_mapping_quality_cap(50)
                    .with_random_seed(7),
            )
            .is_ok()
        );
    }
}
