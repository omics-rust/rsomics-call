use std::{cell::Cell, sync::LazyLock};

use rsomics_bamio::raw::RawRecord;
use rsomics_pileup::PileupRead;
use smallvec::SmallVec;

use crate::{CallError, Result, SampleAnnotations, SiteAnnotations, model::SiteAnnotationValues};

const ALLELE_COUNT: usize = 5;
const POSITION_BINS: usize = 100;
const QUALITY_BINS: usize = 60;
const MISMATCH_BINS: usize = 32;
const EMPTY_POSITION_BINS: [u64; POSITION_BINS] = [0; POSITION_BINS];
const EMPTY_QUALITY_BINS: [u64; QUALITY_BINS] = [0; QUALITY_BINS];
static ERROR_PROBABILITIES: LazyLock<[f64; 64]> =
    LazyLock::new(|| std::array::from_fn(|quality| 10.0f64.powf(-0.1 * quality as f64)));

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CigarMetrics {
    left_soft_clip: usize,
    right_soft_clip: usize,
    has_soft_clip: bool,
    mismatch_adjustment: i64,
}

impl CigarMetrics {
    pub(crate) fn new(cigar: impl IntoIterator<Item = (u8, u32)>) -> Self {
        let mut left_soft_clip = 0;
        let mut right_soft_clip = 0;
        let mut has_core = false;
        let mut has_soft_clip = false;
        let mut mismatch_adjustment = 0;
        for (kind, length) in cigar {
            let length_usize = length as usize;
            match kind {
                4 => {
                    has_soft_clip = true;
                    mismatch_adjustment += i64::from(length);
                    if !has_core {
                        left_soft_clip += length_usize;
                    }
                    right_soft_clip += length_usize;
                }
                5 => {}
                1 | 2 => {
                    has_core = true;
                    right_soft_clip = 0;
                    mismatch_adjustment -= i64::from(length.saturating_sub(1));
                }
                _ => {
                    has_core = true;
                    right_soft_clip = 0;
                }
            }
        }
        Self {
            left_soft_clip,
            right_soft_clip,
            has_soft_clip,
            mismatch_adjustment,
        }
    }
}

#[derive(Default)]
pub(crate) struct PileupRecordState {
    cigar: Cell<Option<CigarMetrics>>,
}

impl PileupRecordState {
    pub(crate) fn cigar(&self, values: impl IntoIterator<Item = (u8, u32)>) -> CigarMetrics {
        self.cigar.get().unwrap_or_else(|| {
            let cigar = CigarMetrics::new(values);
            self.cigar.set(Some(cigar));
            cigar
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) struct AnnotationObservation {
    pub(crate) allele: usize,
    pub(crate) is_reference: bool,
    pub(crate) reverse: bool,
    pub(crate) base_quality: u8,
    pub(crate) raw_mapping_quality: u8,
    pub(crate) mapping_quality: u8,
    pub(crate) effective_quality: u8,
    pub(crate) tail_distance: u8,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AnnotationEvidence {
    forward: [u64; ALLELE_COUNT],
    reverse: [u64; ALLELE_COUNT],
    error_sums: [f64; ALLELE_COUNT],
    auxiliary: [u64; 16],
    detailed: Option<Box<DetailedAnnotationEvidence>>,
    raw_depth: u64,
    zero_mapping_quality: u64,
    soft_clipped_reads: u64,
}

#[derive(Clone, Debug, PartialEq)]
struct DetailedAnnotationEvidence {
    reference_positions: [u64; POSITION_BINS],
    alternate_positions: [u64; POSITION_BINS],
    reference_mapping_qualities: [u64; QUALITY_BINS],
    alternate_mapping_qualities: [u64; QUALITY_BINS],
    reference_base_qualities: [u64; QUALITY_BINS],
    alternate_base_qualities: [u64; QUALITY_BINS],
    forward_mapping_qualities: [u64; QUALITY_BINS],
    reverse_mapping_qualities: [u64; QUALITY_BINS],
    reference_soft_clips: [u64; POSITION_BINS],
    alternate_soft_clips: [u64; POSITION_BINS],
    indel_reference_positions: [u64; POSITION_BINS],
    indel_alternate_positions: [u64; POSITION_BINS],
    indel_reference_mapping_qualities: [u64; QUALITY_BINS],
    indel_alternate_mapping_qualities: [u64; QUALITY_BINS],
    indel_reference_soft_clips: [u64; POSITION_BINS],
    indel_alternate_soft_clips: [u64; POSITION_BINS],
    reference_mismatches: [u64; MISMATCH_BINS],
    alternate_mismatches: [u64; MISMATCH_BINS],
    mismatch_sums: [u64; 2],
    mismatch_counts: [u64; 2],
}

impl Default for DetailedAnnotationEvidence {
    fn default() -> Self {
        Self {
            reference_positions: [0; POSITION_BINS],
            alternate_positions: [0; POSITION_BINS],
            reference_mapping_qualities: [0; QUALITY_BINS],
            alternate_mapping_qualities: [0; QUALITY_BINS],
            reference_base_qualities: [0; QUALITY_BINS],
            alternate_base_qualities: [0; QUALITY_BINS],
            forward_mapping_qualities: [0; QUALITY_BINS],
            reverse_mapping_qualities: [0; QUALITY_BINS],
            reference_soft_clips: [0; POSITION_BINS],
            alternate_soft_clips: [0; POSITION_BINS],
            indel_reference_positions: [0; POSITION_BINS],
            indel_alternate_positions: [0; POSITION_BINS],
            indel_reference_mapping_qualities: [0; QUALITY_BINS],
            indel_alternate_mapping_qualities: [0; QUALITY_BINS],
            indel_reference_soft_clips: [0; POSITION_BINS],
            indel_alternate_soft_clips: [0; POSITION_BINS],
            reference_mismatches: [0; MISMATCH_BINS],
            alternate_mismatches: [0; MISMATCH_BINS],
            mismatch_sums: [0; 2],
            mismatch_counts: [0; 2],
        }
    }
}

impl Default for AnnotationEvidence {
    fn default() -> Self {
        Self {
            forward: [0; ALLELE_COUNT],
            reverse: [0; ALLELE_COUNT],
            error_sums: [0.0; ALLELE_COUNT],
            auxiliary: [0; 16],
            detailed: None,
            raw_depth: 0,
            zero_mapping_quality: 0,
            soft_clipped_reads: 0,
        }
    }
}

impl AnnotationEvidence {
    pub(crate) fn clear(&mut self) {
        self.forward.fill(0);
        self.reverse.fill(0);
        self.error_sums.fill(0.0);
        self.auxiliary.fill(0);
        self.raw_depth = 0;
        self.zero_mapping_quality = 0;
        self.soft_clipped_reads = 0;
        self.detailed = None;
    }

    pub(crate) fn begin_read(&mut self, cigar: CigarMetrics) {
        if cigar.has_soft_clip {
            self.soft_clipped_reads += 1;
        }
    }

    pub(crate) fn add_raw_depth(&mut self) {
        self.raw_depth += 1;
    }

    pub(crate) fn add_allele_depth(&mut self, allele: usize, forward: u64, reverse: u64) {
        self.forward[allele] += forward;
        self.reverse[allele] += reverse;
    }

    pub(crate) fn observe_indel_candidate(
        &mut self,
        record: &RawRecord,
        projection: &PileupRead,
        cigar: CigarMetrics,
    ) {
        let detailed = self
            .detailed
            .get_or_insert_with(|| Box::new(DetailedAnnotationEvidence::default()));
        let geometry = ReadGeometry::new(record.sequence_len(), projection.qpos, cigar);
        let position = geometry.position_bin();
        let soft_clip = geometry.soft_clip_bin();
        let mapping_quality = usize::from(record.mapping_quality().min(59));
        if projection.indel == 0 {
            detailed.indel_reference_positions[position] += 1;
            detailed.indel_reference_mapping_qualities[mapping_quality] += 1;
            detailed.indel_reference_soft_clips[soft_clip] += 1;
        } else {
            detailed.indel_alternate_positions[position] += 1;
            detailed.indel_alternate_mapping_qualities[mapping_quality] += 1;
            detailed.indel_alternate_soft_clips[soft_clip] += 1;
        }
    }

    pub(crate) fn observe(&mut self, observation: AnnotationObservation) {
        let AnnotationObservation {
            allele,
            is_reference,
            reverse,
            base_quality,
            raw_mapping_quality,
            mapping_quality,
            effective_quality,
            tail_distance,
        } = observation;
        if allele < 4 {
            if reverse {
                self.reverse[allele] += 1;
            } else {
                self.forward[allele] += 1;
            }
            self.error_sums[allele] += ERROR_PROBABILITIES[usize::from(effective_quality)];
        }

        let difference = usize::from(!is_reference);
        let strand = usize::from(reverse);
        self.auxiliary[difference * 2 + strand] += 1;
        self.auxiliary[4 + difference * 2] += u64::from(base_quality);
        self.auxiliary[5 + difference * 2] += u64::from(base_quality).pow(2);
        self.auxiliary[8 + difference * 2] += u64::from(mapping_quality);
        self.auxiliary[9 + difference * 2] += u64::from(mapping_quality).pow(2);

        let tail_distance = u64::from(tail_distance);
        self.auxiliary[12 + difference * 2] += tail_distance;
        self.auxiliary[13 + difference * 2] += tail_distance * tail_distance;

        if raw_mapping_quality == 0 {
            self.zero_mapping_quality += 1;
        }
    }

    pub(crate) fn observe_detailed(
        &mut self,
        record: &RawRecord,
        projection: &PileupRead,
        cigar: CigarMetrics,
        observation: AnnotationObservation,
    ) {
        let detailed = self
            .detailed
            .get_or_insert_with(|| Box::new(DetailedAnnotationEvidence::default()));
        let AnnotationObservation {
            is_reference,
            reverse,
            base_quality,
            mapping_quality,
            ..
        } = observation;
        let difference = usize::from(!is_reference);
        let geometry = ReadGeometry::new(record.sequence_len(), projection.qpos, cigar);
        let position = geometry.position_bin();
        let soft_clip = geometry.soft_clip_bin();
        let base_quality_bin = usize::from(base_quality.min(59));
        let mapping_quality_bin = usize::from(mapping_quality.min(59));
        if is_reference {
            detailed.reference_positions[position] += 1;
            detailed.reference_mapping_qualities[mapping_quality_bin] += 1;
            detailed.reference_base_qualities[base_quality_bin] += 1;
            detailed.reference_soft_clips[soft_clip] += 1;
        } else {
            detailed.alternate_positions[position] += 1;
            detailed.alternate_mapping_qualities[mapping_quality_bin] += 1;
            detailed.alternate_base_qualities[base_quality_bin] += 1;
            detailed.alternate_soft_clips[soft_clip] += 1;
        }
        if reverse {
            detailed.reverse_mapping_qualities[mapping_quality_bin] += 1;
        } else {
            detailed.forward_mapping_qualities[mapping_quality_bin] += 1;
        }

        if let Some(mismatches) = mismatch_count(record, cigar, is_reference) {
            detailed.mismatch_sums[difference] += mismatches as u64;
            detailed.mismatch_counts[difference] += 1;
            if is_reference {
                detailed.reference_mismatches[mismatches] += 1;
            } else {
                detailed.alternate_mismatches[mismatches] += 1;
            }
        }
    }

    pub(crate) fn sample_annotations(&self, alleles: &[usize]) -> Result<SampleAnnotations> {
        let forward = selected_counts(&self.forward, alleles)?;
        let reverse = selected_counts(&self.reverse, alleles)?;
        let quality_means = alleles
            .iter()
            .map(|&allele| {
                let depth = self.forward[allele] + self.reverse[allele];
                if depth == 0 || self.error_sums[allele] == 0.0 {
                    0
                } else {
                    (-4.3429 * (self.error_sums[allele] / depth as f64).ln() + 0.499)
                        .min(f64::from(i32::MAX)) as u32
                }
            })
            .collect::<SmallVec<[u32; 5]>>();
        let [
            forward_reference,
            reverse_reference,
            forward_alternate,
            reverse_alternate,
        ] = self.strand_counts();
        let strand_bias = if forward_reference + reverse_reference < 2
            || forward_alternate + reverse_alternate < 2
            || forward_reference + forward_alternate < 2
            || reverse_reference + reverse_alternate < 2
        {
            0
        } else {
            let (_, two_tail) = fisher_exact(
                forward_reference,
                reverse_reference,
                forward_alternate,
                reverse_alternate,
            );
            (-4.343 * two_tail.ln() + 0.499).min(255.0) as u32
        };
        Ok(SampleAnnotations::from_generated(
            forward,
            reverse,
            quality_means,
            strand_bias,
            checked_u32(self.soft_clipped_reads)?,
        ))
    }

    fn strand_counts(&self) -> [u64; 4] {
        [
            self.auxiliary[0],
            self.auxiliary[1],
            self.auxiliary[2],
            self.auxiliary[3],
        ]
    }
}

pub(crate) fn site_annotations<'a>(
    evidence: impl IntoIterator<Item = &'a AnnotationEvidence>,
    indel: bool,
    has_alternate: bool,
) -> Result<SiteAnnotations> {
    let evidence = evidence.into_iter().collect::<Vec<_>>();
    let mut total = AnnotationEvidence::default();
    let mut detailed = has_alternate.then(DetailedAnnotationEvidence::default);
    for sample in &evidence {
        add_arrays(&mut total.forward, &sample.forward);
        add_arrays(&mut total.reverse, &sample.reverse);
        add_arrays(&mut total.auxiliary, &sample.auxiliary);
        if let (Some(total), Some(sample)) = (&mut detailed, sample.detailed.as_deref()) {
            add_arrays(&mut total.reference_positions, &sample.reference_positions);
            add_arrays(&mut total.alternate_positions, &sample.alternate_positions);
            add_arrays(
                &mut total.reference_mapping_qualities,
                &sample.reference_mapping_qualities,
            );
            add_arrays(
                &mut total.alternate_mapping_qualities,
                &sample.alternate_mapping_qualities,
            );
            add_arrays(
                &mut total.reference_base_qualities,
                &sample.reference_base_qualities,
            );
            add_arrays(
                &mut total.alternate_base_qualities,
                &sample.alternate_base_qualities,
            );
            add_arrays(
                &mut total.forward_mapping_qualities,
                &sample.forward_mapping_qualities,
            );
            add_arrays(
                &mut total.reverse_mapping_qualities,
                &sample.reverse_mapping_qualities,
            );
            add_arrays(
                &mut total.reference_soft_clips,
                &sample.reference_soft_clips,
            );
            add_arrays(
                &mut total.alternate_soft_clips,
                &sample.alternate_soft_clips,
            );
            add_arrays(
                &mut total.indel_reference_positions,
                &sample.indel_reference_positions,
            );
            add_arrays(
                &mut total.indel_alternate_positions,
                &sample.indel_alternate_positions,
            );
            add_arrays(
                &mut total.indel_reference_mapping_qualities,
                &sample.indel_reference_mapping_qualities,
            );
            add_arrays(
                &mut total.indel_alternate_mapping_qualities,
                &sample.indel_alternate_mapping_qualities,
            );
            add_arrays(
                &mut total.indel_reference_soft_clips,
                &sample.indel_reference_soft_clips,
            );
            add_arrays(
                &mut total.indel_alternate_soft_clips,
                &sample.indel_alternate_soft_clips,
            );
            add_arrays(
                &mut total.reference_mismatches,
                &sample.reference_mismatches,
            );
            add_arrays(
                &mut total.alternate_mismatches,
                &sample.alternate_mismatches,
            );
            add_arrays(&mut total.mismatch_sums, &sample.mismatch_sums);
            add_arrays(&mut total.mismatch_counts, &sample.mismatch_counts);
        }
        total.raw_depth += sample.raw_depth;
        total.zero_mapping_quality += sample.zero_mapping_quality;
        total.soft_clipped_reads += sample.soft_clipped_reads;
    }

    let [
        forward_reference,
        reverse_reference,
        forward_alternate,
        reverse_alternate,
    ] = total.strand_counts();
    let strand_bias = has_alternate.then(|| {
        fisher_exact(
            forward_reference,
            reverse_reference,
            forward_alternate,
            reverse_alternate,
        )
        .0 as f32
    });
    let detailed = detailed.as_ref();
    let average_mismatches = detailed
        .filter(|values| values.mismatch_counts[0] != 0 || values.mismatch_counts[1] != 0)
        .map(|values| {
            [
                if values.mismatch_counts[0] == 0 {
                    0.0
                } else {
                    values.mismatch_sums[0] as f32 / values.mismatch_counts[0] as f32
                },
                if values.mismatch_counts[1] == 0 {
                    0.0
                } else {
                    values.mismatch_sums[1] as f32 / values.mismatch_counts[1] as f32
                },
            ]
        });
    let (reference_positions, alternate_positions) = if indel {
        detailed.map_or((&EMPTY_POSITION_BINS, &EMPTY_POSITION_BINS), |detailed| {
            (
                &detailed.indel_reference_positions,
                &detailed.indel_alternate_positions,
            )
        })
    } else {
        detailed.map_or((&EMPTY_POSITION_BINS, &EMPTY_POSITION_BINS), |detailed| {
            (&detailed.reference_positions, &detailed.alternate_positions)
        })
    };
    let (reference_mapping_qualities, alternate_mapping_qualities) = if indel {
        detailed.map_or((&EMPTY_QUALITY_BINS, &EMPTY_QUALITY_BINS), |detailed| {
            (
                &detailed.indel_reference_mapping_qualities,
                &detailed.indel_alternate_mapping_qualities,
            )
        })
    } else {
        detailed.map_or((&EMPTY_QUALITY_BINS, &EMPTY_QUALITY_BINS), |detailed| {
            (
                &detailed.reference_mapping_qualities,
                &detailed.alternate_mapping_qualities,
            )
        })
    };
    let (reference_soft_clips, alternate_soft_clips) = if indel {
        detailed.map_or((&EMPTY_POSITION_BINS, &EMPTY_POSITION_BINS), |detailed| {
            (
                &detailed.indel_reference_soft_clips,
                &detailed.indel_alternate_soft_clips,
            )
        })
    } else {
        detailed.map_or((&EMPTY_POSITION_BINS, &EMPTY_POSITION_BINS), |detailed| {
            (
                &detailed.reference_soft_clips,
                &detailed.alternate_soft_clips,
            )
        })
    };
    SiteAnnotations::new(SiteAnnotationValues {
        raw_depth: checked_u32(total.raw_depth)?,
        auxiliary: total.auxiliary.map(|value| value as f32),
        variant_distance_bias: has_alternate
            .then(|| detailed.and_then(|values| variant_distance_bias(&values.alternate_positions)))
            .flatten()
            .map(|value| value as f32),
        read_position_bias: has_alternate
            .then(|| mann_whitney_z(reference_positions, alternate_positions))
            .flatten()
            .map(|value| value as f32),
        mapping_quality_bias: has_alternate
            .then(|| mann_whitney_z(reference_mapping_qualities, alternate_mapping_qualities))
            .flatten()
            .map(|value| value as f32),
        base_quality_bias: if !has_alternate {
            None
        } else if indel {
            Some(0.0)
        } else {
            detailed
                .and_then(|values| {
                    mann_whitney_z(
                        &values.reference_base_qualities,
                        &values.alternate_base_qualities,
                    )
                })
                .map(|value| value as f32)
        },
        mapping_quality_strand_bias: if !has_alternate {
            None
        } else if indel {
            Some(0.0)
        } else {
            detailed
                .and_then(|values| {
                    mann_whitney_z(
                        &values.forward_mapping_qualities,
                        &values.reverse_mapping_qualities,
                    )
                })
                .map(|value| value as f32)
        },
        mismatch_bias: has_alternate
            .then(|| {
                detailed.and_then(|values| {
                    mann_whitney_z(&values.reference_mismatches, &values.alternate_mismatches)
                })
            })
            .flatten()
            .map(|value| value as f32),
        soft_clip_bias: has_alternate
            .then(|| mann_whitney_z(reference_soft_clips, alternate_soft_clips))
            .flatten()
            .map(|value| value as f32),
        strand_bias,
        segregation_bias: has_alternate
            .then(|| segregation_bias(&evidence))
            .flatten()
            .map(|value| value as f32),
        zero_mapping_quality_fraction: if total.raw_depth == 0 {
            0.0
        } else {
            total.zero_mapping_quality as f32 / total.raw_depth as f32
        },
        average_mismatches,
    })
}

pub(crate) fn selected_counts(
    values: &[u64; ALLELE_COUNT],
    alleles: &[usize],
) -> Result<SmallVec<[u32; 5]>> {
    alleles
        .iter()
        .map(|&allele| checked_u32(values[allele]))
        .collect()
}

fn checked_u32(value: u64) -> Result<u32> {
    u32::try_from(value).map_err(|_| CallError::SnpEvidenceOverflow)
}

fn add_arrays<const N: usize>(target: &mut [u64; N], source: &[u64; N]) {
    for (target, source) in target.iter_mut().zip(source) {
        *target += source;
    }
}

struct ReadGeometry {
    aligned_length: usize,
    position: usize,
    soft_clip_length: usize,
    soft_clip_distance: usize,
}

impl ReadGeometry {
    fn new(sequence_length: usize, qpos: usize, cigar: CigarMetrics) -> Self {
        let left = cigar.left_soft_clip;
        let right = cigar.right_soft_clip;
        let left_distance = qpos + 1 - left;
        let right_distance = sequence_length - right - qpos;
        let clip = match (left != 0, right != 0) {
            (true, true) if left_distance <= right_distance => (left, left_distance),
            (true, true) => (right, right_distance),
            (true, false) => (left, left_distance),
            (false, true) => (right, right_distance),
            (false, false) => (0, 0),
        };
        Self {
            aligned_length: sequence_length - left - right,
            position: qpos + 1 - left,
            soft_clip_length: clip.0,
            soft_clip_distance: clip.1,
        }
    }

    fn position_bin(&self) -> usize {
        self.position * (POSITION_BINS - 1) / (self.aligned_length + 1)
    }

    fn soft_clip_bin(&self) -> usize {
        (15 * self.soft_clip_length / (self.soft_clip_distance + 1)).min(POSITION_BINS - 1)
    }
}

fn mismatch_count(record: &RawRecord, cigar: CigarMetrics, is_reference: bool) -> Option<usize> {
    let value = record.aux_value(*b"NM")?;
    let mut mismatches = match record.aux_type(*b"NM")? {
        b'c' => i64::from(i8::from_le_bytes(value.try_into().ok()?)),
        b'C' => i64::from(u8::from_le_bytes(value.try_into().ok()?)),
        b's' => i64::from(i16::from_le_bytes(value.try_into().ok()?)),
        b'S' => i64::from(u16::from_le_bytes(value.try_into().ok()?)),
        b'i' => i64::from(i32::from_le_bytes(value.try_into().ok()?)),
        b'I' => i64::from(u32::from_le_bytes(value.try_into().ok()?)),
        _ => return None,
    };
    mismatches += cigar.mismatch_adjustment;
    mismatches -= if is_reference { 1 } else { 2 };
    Some(mismatches.clamp(0, (MISMATCH_BINS - 1) as i64) as usize)
}

fn mann_whitney_z<const N: usize>(reference: &[u64; N], alternate: &[u64; N]) -> Option<f64> {
    let mut equal = 0.0;
    let mut less = 0.0;
    let mut reference_count = 0.0;
    let mut alternate_count = 0.0;
    let mut ties = 0.0;
    for (&reference, &alternate) in reference.iter().zip(alternate).rev() {
        let reference = reference as f64;
        let alternate = alternate as f64;
        equal += reference * alternate;
        less += reference * alternate_count;
        reference_count += reference;
        alternate_count += alternate;
        let combined = reference + alternate;
        ties += (combined * combined - 1.0) * combined;
    }
    if reference_count == 0.0 || alternate_count == 0.0 {
        return None;
    }
    let score = less + equal * 0.5;
    let mean = reference_count * alternate_count * 0.5;
    let combined = reference_count + alternate_count;
    let variance = reference_count * alternate_count / 12.0
        * (combined + 1.0 - ties / (combined * (combined - 1.0)));
    Some(if variance <= 0.0 {
        0.0
    } else {
        (score - mean) / variance.sqrt()
    })
}

fn segregation_bias(evidence: &[&AnnotationEvidence]) -> Option<f64> {
    let alternate_count = evidence
        .iter()
        .map(|sample| sample.auxiliary[2] + sample.auxiliary[3])
        .sum::<u64>();
    if alternate_count == 0 {
        return None;
    }
    let sample_count = evidence.len() as f64;
    let depth = evidence
        .iter()
        .flat_map(|sample| sample.auxiliary[..4].iter())
        .sum::<u64>();
    let average_depth = depth / evidence.len() as u64;
    let mut variant_samples = (alternate_count as f64 / average_depth as f64 + 0.5).floor();
    variant_samples = variant_samples.clamp(1.0, sample_count);
    let frequency = variant_samples / (2.0 * sample_count);
    let error_mean = alternate_count as f64 / sample_count;
    let variant_mean = alternate_count as f64 / variant_samples;
    let mut sum = 0.0;
    for sample in evidence {
        let observed = sample.auxiliary[2] + sample.auxiliary[3];
        if observed == 0 {
            sum += (2.0 * frequency * (1.0 - frequency) * (-variant_mean).exp()
                + frequency * frequency * (-2.0 * variant_mean).exp()
                + (1.0 - frequency) * (1.0 - frequency))
                .ln()
                + error_mean;
        } else {
            let observed = observed as f64;
            let mixture = (2.0 * (1.0 - frequency))
                .ln()
                .max(frequency.ln() + observed * 2.0f64.ln() - variant_mean);
            let mixture = (1.0
                + ((2.0 * (1.0 - frequency))
                    .ln()
                    .min(frequency.ln() + observed * 2.0f64.ln() - variant_mean)
                    - mixture)
                    .exp())
            .ln()
                + mixture;
            sum += frequency.ln() + observed * (variant_mean / error_mean).ln() - variant_mean
                + mixture
                + error_mean;
        }
    }
    Some(sum)
}

fn variant_distance_bias(positions: &[u64; POSITION_BINS]) -> Option<f64> {
    const PARAMETERS: [(u64, f32, f32); 15] = [
        (3, 0.079, 18.0),
        (4, 0.09, 19.8),
        (5, 0.1, 20.5),
        (6, 0.11, 21.5),
        (7, 0.125, 21.6),
        (8, 0.135, 22.0),
        (9, 0.14, 22.2),
        (10, 0.153, 22.3),
        (15, 0.19, 22.8),
        (20, 0.22, 23.2),
        (30, 0.26, 23.4),
        (40, 0.29, 23.5),
        (50, 0.35, 23.65),
        (100, 0.5, 23.7),
        (200, 0.7, 23.7),
    ];
    let depth = positions.iter().sum::<u64>();
    if depth < 2 {
        return None;
    }
    let mut mean = positions
        .iter()
        .enumerate()
        .map(|(position, &count)| position as f32 * count as f32)
        .sum::<f32>();
    mean /= depth as f32;
    let mut mean_difference = positions
        .iter()
        .enumerate()
        .map(|(position, &count)| (position as f32 - mean).abs() * count as f32)
        .sum::<f32>();
    mean_difference /= depth as f32;
    let integer_difference = mean_difference as u64;
    if depth == 2 {
        let numerator = (200 - 2 * (integer_difference + 1) - 1) * (integer_difference + 1);
        return Some((numerator / 99) as f64 / 50.0);
    }
    let upper = PARAMETERS
        .iter()
        .position(|parameter| parameter.0 >= depth)
        .unwrap_or(PARAMETERS.len());
    let (scale, shift) = if upper == PARAMETERS.len() {
        (PARAMETERS[14].1, PARAMETERS[14].2)
    } else if upper > 0 && PARAMETERS[upper].0 != depth {
        (
            (PARAMETERS[upper - 1].1 + PARAMETERS[upper].1) * 0.5,
            (PARAMETERS[upper - 1].2 + PARAMETERS[upper].2) * 0.5,
        )
    } else {
        (PARAMETERS[upper].1, PARAMETERS[upper].2)
    };
    let argument = -((mean_difference - shift) * scale);
    Some(0.5 * complementary_error(f64::from(argument)))
}

fn complementary_error(value: f64) -> f64 {
    const SQRT_2: f64 = std::f64::consts::SQRT_2;
    let z = value.abs() * SQRT_2;
    if z > 37.0 {
        return if value > 0.0 { 0.0 } else { 2.0 };
    }
    let exponential = (-0.5 * z * z).exp();
    let probability = if z < 10.0 / SQRT_2 {
        exponential
            * ((((((0.035_262_496_599_891_09 * z + 0.700_383_064_443_688_1) * z
                + 6.373_962_203_531_65)
                * z
                + 33.912_866_078_383)
                * z
                + 112.079_291_497_870_9)
                * z
                + 221.213_596_169_931_1)
                * z
                + 220.206_867_912_376_1)
            / (((((((0.088_388_347_648_318_44 * z + 1.755_667_163_182_642) * z
                + 16.064_177_579_206_95)
                * z
                + 86.780_732_202_946_08)
                * z
                + 296.564_248_779_673_7)
                * z
                + 637.333_633_378_831_1)
                * z
                + 793.826_512_519_948_4)
                * z
                + 440.413_735_824_752_2)
    } else {
        exponential
            / 2.506_628_274_631_001
            / (z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 0.65)))))
    };
    if value > 0.0 {
        2.0 * probability
    } else {
        2.0 * (1.0 - probability)
    }
}

pub(crate) fn fisher_exact(n11: u64, n12: u64, n21: u64, n22: u64) -> (f64, f64) {
    let row = n11 + n12;
    let column = n11 + n21;
    let total = row + n21 + n22;
    let maximum = row.min(column);
    let minimum = (row + column).saturating_sub(total);
    if minimum == maximum {
        return (1.0, 1.0);
    }
    let observed = hypergeometric(n11, row, column, total);
    if observed == 0.0 {
        return (0.0, 0.0);
    }
    let mut left = 0.0;
    let mut left_boundary = minimum;
    for value in minimum..=maximum {
        let probability = hypergeometric(value, row, column, total);
        if probability < 0.999_999_99 * observed {
            left += probability;
            left_boundary = value;
        } else {
            if probability < 1.000_000_01 * observed {
                left += probability;
                left_boundary = value;
            }
            break;
        }
    }
    let mut right = 0.0;
    let mut right_boundary = maximum;
    for value in (minimum..=maximum).rev() {
        let probability = hypergeometric(value, row, column, total);
        if probability < 0.999_999_99 * observed {
            right += probability;
            right_boundary = value;
        } else {
            if probability < 1.000_000_01 * observed {
                right += probability;
                right_boundary = value;
            }
            break;
        }
    }
    let two_tail = (left + right).min(1.0);
    if left_boundary.abs_diff(n11) < right_boundary.abs_diff(n11) {
        right = 1.0 - left + observed;
    } else {
        left = 1.0 - right + observed;
    }
    let _ = (left, right);
    (observed, two_tail)
}

fn hypergeometric(value: u64, row: u64, column: u64, total: u64) -> f64 {
    (ln_binomial(row, value) + ln_binomial(total - row, column - value)
        - ln_binomial(total, column))
    .exp()
}

fn ln_binomial(n: u64, k: u64) -> f64 {
    if k == 0 || n == k {
        0.0
    } else {
        ln_gamma((n + 1) as f64) - ln_gamma((k + 1) as f64) - ln_gamma((n - k + 1) as f64)
    }
}

pub(crate) fn ln_gamma(value: f64) -> f64 {
    let mut sum = 0.165_947_018_740_846_2e-6 / (value + 7.0);
    sum += 0.993_493_711_393_074_8e-5 / (value + 6.0);
    sum -= 0.138_571_033_129_652_6 / (value + 5.0);
    sum += 12.507_343_240_090_56 / (value + 4.0);
    sum -= 176.615_029_149_838_6 / (value + 3.0);
    sum += 771.323_428_775_767_4 / (value + 2.0);
    sum -= 1_259.139_216_722_289 / (value + 1.0);
    sum += 676.520_368_121_883_5 / value;
    sum += 0.999_999_999_999_518_3;
    sum.ln() - 5.581_061_466_795_328 - value + (value - 0.5) * (value + 6.5).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pileup_record_state_computes_cigar_metrics_once() {
        let state = PileupRecordState::default();
        let first = state.cigar([(4, 2), (0, 10)]);
        let second = state.cigar([(0, 3)]);

        assert_eq!(first, second);
        assert_ne!(first, CigarMetrics::new([(0, 3)]));
    }

    #[test]
    fn fisher_exact_matches_htslib_cases() {
        for (table, expected_probability, expected_two_tail) in [
            ([2, 1, 0, 31], 0.005_347_593_583, 0.005_347_593_583),
            ([2, 1, 0, 1], 0.5, 1.0),
            ([3, 1, 0, 0], 1.0, 1.0),
            ([3, 15, 37, 45], 0.017_138_952_733, 0.033_161_943_699),
            ([12, 5, 29, 2], 0.039_079_943_857, 0.080_268_552_074),
        ] {
            let (probability, two_tail) = fisher_exact(table[0], table[1], table[2], table[3]);
            assert!((probability - expected_probability).abs() < 1e-10);
            assert!((two_tail - expected_two_tail).abs() < 1e-10);
        }
    }

    #[test]
    fn mann_whitney_ties_match_bcftools_zero_case() {
        let reference = [0, 2, 0, 0];
        let alternate = [0, 2, 0, 0];
        assert_eq!(mann_whitney_z(&reference, &alternate), Some(0.0));
    }
}
