use rsomics_bamio::raw::RawRecord;
use rsomics_pileup::Column;

use crate::{
    Allele, BaseObservation, CallError, ErrorModel, LikelihoodMatrix, LikelihoodSite, Nucleotide,
    Ploidy, Result, SampleLikelihood,
};

const NUCLEOTIDE: [Nucleotide; 16] = [
    Nucleotide::N,
    Nucleotide::A,
    Nucleotide::C,
    Nucleotide::N,
    Nucleotide::G,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::T,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::N,
    Nucleotide::N,
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SnpLikelihoodConfig {
    minimum_base_quality: u8,
    maximum_base_quality: u8,
    neighboring_quality_delta: u8,
    mapping_quality_cap: u8,
    random_seed: i32,
}

impl SnpLikelihoodConfig {
    pub fn new(
        minimum_base_quality: u8,
        maximum_base_quality: u8,
        neighboring_quality_delta: u8,
        mapping_quality_cap: u8,
    ) -> Result<Self> {
        if minimum_base_quality > maximum_base_quality {
            return Err(CallError::InvalidSnpQualityRange);
        }
        Ok(Self {
            minimum_base_quality,
            maximum_base_quality,
            neighboring_quality_delta,
            mapping_quality_cap,
            random_seed: 0,
        })
    }

    pub fn with_random_seed(mut self, seed: i32) -> Self {
        self.random_seed = seed;
        self
    }
}

impl Default for SnpLikelihoodConfig {
    fn default() -> Self {
        Self {
            minimum_base_quality: 1,
            maximum_base_quality: 60,
            neighboring_quality_delta: 30,
            mapping_quality_cap: 60,
            random_seed: 0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SnpEvidence {
    likelihoods: LikelihoodMatrix,
    depth: u32,
    allele_depths: [u32; 5],
    quality_sums: [u32; 5],
}

impl SnpEvidence {
    pub fn from_column(
        column: &Column<'_>,
        reference_base: Nucleotide,
        config: SnpLikelihoodConfig,
        model: &mut ErrorModel,
        observations: &mut Vec<BaseObservation>,
    ) -> Result<Self> {
        let mut accumulator = SnpAccumulator::with_observations(std::mem::take(observations));
        for entry in column.entries() {
            accumulator.push(entry.record(), entry.projection(), reference_base, config);
        }
        let evidence = accumulator.finish(model)?;
        *observations = accumulator.observations;
        Ok(evidence)
    }

    pub fn likelihoods(&self) -> &LikelihoodMatrix {
        &self.likelihoods
    }

    pub fn depth(&self) -> u32 {
        self.depth
    }

    pub fn allele_depths(&self) -> &[u32; 5] {
        &self.allele_depths
    }

    pub fn quality_sums(&self) -> &[u32; 5] {
        &self.quality_sums
    }
}

#[derive(Default)]
struct SnpAccumulator {
    observations: Vec<BaseObservation>,
    depth: u64,
    allele_depths: [u64; 5],
    quality_sums: [u64; 5],
}

impl SnpAccumulator {
    fn with_observations(mut observations: Vec<BaseObservation>) -> Self {
        observations.clear();
        Self {
            observations,
            ..Self::default()
        }
    }

    fn clear(&mut self) {
        self.observations.clear();
        self.depth = 0;
        self.allele_depths.fill(0);
        self.quality_sums.fill(0);
    }

    fn push(
        &mut self,
        record: &RawRecord,
        projection: &rsomics_pileup::PileupRead,
        reference_base: Nucleotide,
        config: SnpLikelihoodConfig,
    ) {
        if projection.is_deletion || projection.is_reference_skip {
            return;
        }
        self.depth += 1;

        let qualities = record.quality_scores();
        let qpos = projection.qpos;
        let quality = |index| {
            if qualities.is_empty() {
                u8::MAX
            } else {
                qualities[index]
            }
        };
        let mut base_quality = u16::from(quality(qpos));
        if qpos > 0 {
            base_quality = base_quality
                .min(u16::from(quality(qpos - 1)) + u16::from(config.neighboring_quality_delta));
        }
        if qpos + 1 < record.sequence_len() {
            base_quality = base_quality
                .min(u16::from(quality(qpos + 1)) + u16::from(config.neighboring_quality_delta));
        }
        if base_quality < u16::from(config.minimum_base_quality) {
            return;
        }
        let base_quality =
            u8::try_from(base_quality.min(u16::from(config.maximum_base_quality))).unwrap();
        let mapping_quality = match record.mapping_quality() {
            255 => 20,
            quality => quality,
        }
        .min(config.mapping_quality_cap);
        let effective_quality = base_quality.min(mapping_quality).clamp(4, 63);
        let nibble = record.seq_nibble(qpos);
        let base = if nibble == 0 {
            reference_base
        } else {
            NUCLEOTIDE[usize::from(nibble)]
        };
        self.allele_depths[base as usize] += 1;
        self.quality_sums[base as usize] += u64::from(effective_quality);
        self.observations.push(BaseObservation::new(
            base,
            effective_quality,
            record.flags() & 0x10 != 0,
        ));
    }

    fn finish(&mut self, model: &mut ErrorModel) -> Result<SnpEvidence> {
        Ok(SnpEvidence {
            likelihoods: model.calculate(&mut self.observations)?,
            depth: u32::try_from(self.depth).map_err(|_| CallError::SnpEvidenceOverflow)?,
            allele_depths: checked_counts(self.allele_depths)?,
            quality_sums: checked_counts(self.quality_sums)?,
        })
    }
}

pub struct SnpSiteBuilder {
    config: SnpLikelihoodConfig,
    model: ErrorModel,
    samples: Vec<SnpAccumulator>,
}

impl SnpSiteBuilder {
    pub fn new(sample_count: usize, config: SnpLikelihoodConfig) -> Result<Self> {
        if sample_count == 0 {
            return Err(CallError::InvalidSampleCount);
        }
        Ok(Self {
            config,
            model: ErrorModel::with_random_seed(config.random_seed),
            samples: std::iter::repeat_with(SnpAccumulator::default)
                .take(sample_count)
                .collect(),
        })
    }

    pub fn build(
        &mut self,
        column: &Column<'_>,
        reference_base: Nucleotide,
        mut sample_index: impl FnMut(u32, &RawRecord) -> Result<Option<usize>>,
    ) -> Result<LikelihoodSite> {
        for sample in &mut self.samples {
            sample.clear();
        }
        for entry in column.entries() {
            let Some(index) = sample_index(entry.source_id(), entry.record())? else {
                continue;
            };
            let count = self.samples.len();
            let sample = self
                .samples
                .get_mut(index)
                .ok_or(CallError::InvalidSampleIndex { index, count })?;
            sample.push(
                entry.record(),
                entry.projection(),
                reference_base,
                self.config,
            );
        }

        let model = &mut self.model;
        let evidence = self
            .samples
            .iter_mut()
            .map(|sample| sample.finish(model))
            .collect::<Result<Vec<_>>>()?;
        let alleles = selected_alleles(reference_base, &evidence)?;
        let allele_count = alleles.len();
        let diploid = Ploidy::new(2).unwrap();
        let samples = evidence
            .iter()
            .map(|sample| {
                let matrix = sample.likelihoods();
                let mut raw = Vec::with_capacity(allele_count * (allele_count + 1) / 2);
                for second in 0..allele_count {
                    for first in 0..=second {
                        raw.push(
                            matrix.get(alleles[first].matrix_base, alleles[second].matrix_base),
                        );
                    }
                }
                let minimum = raw.iter().copied().fold(f32::INFINITY, f32::min);
                let phred_likelihoods = raw
                    .into_iter()
                    .map(|value| ((f64::from(value - minimum) + 0.499) as u32).min(255))
                    .collect::<Vec<_>>();
                let allele_depths = alleles
                    .iter()
                    .map(|allele| sample.allele_depths()[allele.matrix_base as usize])
                    .collect::<Vec<_>>();
                SampleLikelihood::observed(
                    allele_count,
                    diploid,
                    phred_likelihoods,
                    sample.depth(),
                    allele_depths,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let reference_sequence_id = usize::try_from(column.reference_id())
            .map_err(|_| CallError::InvalidPileupCoordinate)?;
        let position =
            u64::try_from(column.position()).map_err(|_| CallError::InvalidPileupCoordinate)?;
        let mut alleles = alleles.into_iter();
        let reference = alleles.next().unwrap().allele;
        let alternates = alleles.map(|allele| allele.allele).collect::<Vec<_>>();
        LikelihoodSite::new(
            reference_sequence_id,
            position,
            reference,
            alternates,
            samples,
        )
    }
}

struct SelectedAllele {
    allele: Allele,
    matrix_base: Nucleotide,
}

fn selected_alleles(reference: Nucleotide, samples: &[SnpEvidence]) -> Result<Vec<SelectedAllele>> {
    let mut quality_sums = [0.0f32; 4];
    for sample in samples {
        let total = sample.quality_sums()[..4]
            .iter()
            .map(|&value| u64::from(value))
            .sum::<u64>();
        if total != 0 {
            for (sum, &value) in quality_sums.iter_mut().zip(&sample.quality_sums()[..4]) {
                *sum += value as f32 / total as f32;
            }
        }
    }

    let mut order = [0usize, 1, 2, 3];
    for index in 1..order.len() {
        let mut cursor = index;
        while cursor > 0 && quality_sums[order[cursor]] < quality_sums[order[cursor - 1]] {
            order.swap(cursor, cursor - 1);
            cursor -= 1;
        }
    }

    let reference_index = reference as usize;
    let mut selected = vec![SelectedAllele {
        allele: Allele::new([nucleotide_byte(reference)])?,
        matrix_base: reference,
    }];
    let mut unseen = None;
    for &index in order.iter().rev() {
        if index == reference_index {
            continue;
        }
        if quality_sums[index] == 0.0 {
            unseen = Some(index);
            break;
        }
        let base = NUCLEOTIDE[index_to_nibble(index)];
        selected.push(SelectedAllele {
            allele: Allele::new([nucleotide_byte(base)])?,
            matrix_base: base,
        });
    }
    let can_add_unseen = if reference == Nucleotide::N {
        selected.len() < 5
    } else {
        selected.len() < 4
    };
    if can_add_unseen && let Some(index) = unseen {
        selected.push(SelectedAllele {
            allele: Allele::new(&b"<*>"[..])?,
            matrix_base: NUCLEOTIDE[index_to_nibble(index)],
        });
    }
    Ok(selected)
}

fn index_to_nibble(index: usize) -> usize {
    [1, 2, 4, 8][index]
}

fn nucleotide_byte(base: Nucleotide) -> u8 {
    [b'A', b'C', b'G', b'T', b'N'][base as usize]
}

fn checked_counts(values: [u64; 5]) -> Result<[u32; 5]> {
    let mut output = [0; 5];
    for (output, value) in output.iter_mut().zip(values) {
        *output = u32::try_from(value).map_err(|_| CallError::SnpEvidenceOverflow)?;
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use rsomics_bamio::raw::RawRecord;
    use rsomics_pileup::{PileupEngine, PileupOptions};

    use super::*;
    use crate::{SampleMapBuilder, SampleSelection};

    struct RecordSpec<'a> {
        name: &'a [u8],
        flags: u16,
        cigar: &'a [(u8, u32)],
        bases: &'a [u8],
        qualities: &'a [u8],
        mapping_quality: u8,
    }

    fn record(spec: RecordSpec<'_>) -> RawRecord {
        let mut payload = Vec::new();
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.push(u8::try_from(spec.name.len() + 1).unwrap());
        payload.push(spec.mapping_quality);
        payload.extend_from_slice(&0u16.to_le_bytes());
        payload.extend_from_slice(&(spec.cigar.len() as u16).to_le_bytes());
        payload.extend_from_slice(&spec.flags.to_le_bytes());
        payload.extend_from_slice(&(spec.bases.len() as u32).to_le_bytes());
        payload.extend_from_slice(&(-1i32).to_le_bytes());
        payload.extend_from_slice(&(-1i32).to_le_bytes());
        payload.extend_from_slice(&0i32.to_le_bytes());
        payload.extend_from_slice(spec.name);
        payload.push(0);
        for &(kind, length) in spec.cigar {
            payload.extend_from_slice(&((length << 4) | u32::from(kind)).to_le_bytes());
        }
        for pair in spec.bases.chunks(2) {
            let high = base_code(pair[0]);
            let low = pair.get(1).copied().map_or(0, base_code);
            payload.push(high << 4 | low);
        }
        payload.extend_from_slice(spec.qualities);
        RawRecord::try_from(payload).unwrap()
    }

    fn base_code(base: u8) -> u8 {
        match base {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        }
    }

    #[test]
    fn column_contract_applies_neighbor_and_mapping_quality_caps() {
        let records = [
            record(RecordSpec {
                name: b"reference",
                flags: 0,
                cigar: &[(0, 3)],
                bases: b"AAA",
                qualities: &[5, 60, 5],
                mapping_quality: 60,
            }),
            record(RecordSpec {
                name: b"alternate",
                flags: 0x10,
                cigar: &[(0, 3)],
                bases: b"AGA",
                qualities: &[40, 40, 40],
                mapping_quality: 255,
            }),
            record(RecordSpec {
                name: b"deletion",
                flags: 0,
                cigar: &[(0, 1), (2, 1), (0, 1)],
                bases: b"AA",
                qualities: &[30, 30],
                mapping_quality: 60,
            }),
        ];
        let mut engine = PileupEngine::new([10], PileupOptions::default());
        for record in records {
            engine.push(record).unwrap();
        }
        engine.finish().unwrap();

        let mut model = ErrorModel::default();
        let mut observations = Vec::new();
        let mut evidence = None;
        engine
            .drain(|column| {
                if column.position() == 1 {
                    evidence = Some(
                        SnpEvidence::from_column(
                            column,
                            Nucleotide::A,
                            SnpLikelihoodConfig::default(),
                            &mut model,
                            &mut observations,
                        )
                        .unwrap(),
                    );
                }
                Ok::<_, ()>(())
            })
            .unwrap();

        let evidence = evidence.unwrap();
        assert_eq!(evidence.depth(), 2);
        assert_eq!(evidence.allele_depths(), &[1, 0, 1, 0, 0]);
        assert_eq!(evidence.quality_sums(), &[35, 0, 20, 0, 0]);
        let mut expected = [
            BaseObservation::new(Nucleotide::A, 35, false),
            BaseObservation::new(Nucleotide::G, 20, true),
        ];
        assert_eq!(
            evidence.likelihoods(),
            &model.calculate(&mut expected).unwrap()
        );
    }

    #[test]
    fn quality_range_is_checked() {
        assert_eq!(
            SnpLikelihoodConfig::new(30, 20, 10, 60),
            Err(CallError::InvalidSnpQualityRange)
        );
    }

    #[test]
    fn site_builder_matches_bcftools_1_24_reference_only_likelihood() {
        let record = record(RecordSpec {
            name: b"reference",
            flags: 0,
            cigar: &[(0, 1)],
            bases: b"A",
            qualities: &[40],
            mapping_quality: 60,
        });
        let mut engine = PileupEngine::new([10], PileupOptions::default());
        engine.push_with_source(0, record).unwrap();
        engine.finish().unwrap();

        let mut builder = SnpSiteBuilder::new(1, SnpLikelihoodConfig::default()).unwrap();
        let mut site = None;
        engine
            .drain(|column| {
                site = Some(
                    builder
                        .build(column, Nucleotide::A, |source, _| Ok(Some(source as usize)))
                        .unwrap(),
                );
                Ok::<_, ()>(())
            })
            .unwrap();
        let site = site.unwrap();

        assert_eq!(site.reference_sequence_id(), 0);
        assert_eq!(site.position(), 0);
        assert_eq!(site.reference().as_bytes(), b"A");
        assert_eq!(site.alternates()[0].as_bytes(), b"<*>");
        assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[0, 3, 40][..]));
        assert_eq!(site.samples()[0].depth(), 1);
        assert_eq!(site.samples()[0].allele_depths(), &[1, 0]);
    }

    #[test]
    fn site_builder_matches_bcftools_1_24_multisample_likelihoods() {
        let reference = record(RecordSpec {
            name: b"reference",
            flags: 0,
            cigar: &[(0, 1)],
            bases: b"A",
            qualities: &[40],
            mapping_quality: 60,
        });
        let alternate = record(RecordSpec {
            name: b"alternate",
            flags: 0,
            cigar: &[(0, 1)],
            bases: b"G",
            qualities: &[40],
            mapping_quality: 60,
        });
        let mut engine = PileupEngine::new([10], PileupOptions::default());
        engine.push_with_source(7, reference).unwrap();
        engine.push_with_source(8, alternate).unwrap();
        engine.finish().unwrap();

        let mut sample_builder = SampleMapBuilder::new(SampleSelection::default());
        sample_builder
            .add_source(7, "S1", true, Vec::<(&[u8], &str)>::new())
            .unwrap();
        sample_builder
            .add_source(8, "S2", true, Vec::<(&[u8], &str)>::new())
            .unwrap();
        let sample_map = sample_builder.finish().unwrap();
        let mut builder =
            SnpSiteBuilder::new(sample_map.samples().len(), SnpLikelihoodConfig::default())
                .unwrap();
        let mut site = None;
        engine
            .drain(|column| {
                site = Some(
                    builder
                        .build(column, Nucleotide::A, |source, record| {
                            sample_map.sample_index(source, record)
                        })
                        .unwrap(),
                );
                Ok::<_, ()>(())
            })
            .unwrap();
        let site = site.unwrap();

        assert_eq!(
            site.alternates()
                .iter()
                .map(Allele::as_bytes)
                .collect::<Vec<_>>(),
            [b"G".as_slice(), b"<*>".as_slice()]
        );
        assert_eq!(
            site.samples()[0].phred_likelihoods(),
            Some(&[0, 3, 40, 3, 40, 40][..])
        );
        assert_eq!(
            site.samples()[1].phred_likelihoods(),
            Some(&[40, 3, 0, 40, 3, 40][..])
        );
        assert_eq!(site.samples()[0].allele_depths(), &[1, 0, 0]);
        assert_eq!(site.samples()[1].allele_depths(), &[0, 1, 0]);
    }

    #[test]
    fn site_builder_checks_sample_count_and_missing_qualities() {
        assert!(matches!(
            SnpSiteBuilder::new(0, SnpLikelihoodConfig::default()),
            Err(CallError::InvalidSampleCount)
        ));

        let record = record(RecordSpec {
            name: b"missing-quality",
            flags: 0,
            cigar: &[(0, 1)],
            bases: b"A",
            qualities: &[u8::MAX],
            mapping_quality: 60,
        });
        let mut engine = PileupEngine::new([10], PileupOptions::default());
        engine.push_with_source(3, record).unwrap();
        engine.finish().unwrap();
        let mut builder = SnpSiteBuilder::new(1, SnpLikelihoodConfig::default()).unwrap();
        assert_eq!(
            engine.drain(|column| builder
                .build(column, Nucleotide::A, |_, _| Ok(Some(1)))
                .map(drop)),
            Err(CallError::InvalidSampleIndex { index: 1, count: 1 })
        );
        assert!(
            engine
                .drain(|column| builder
                    .build(column, Nucleotide::A, |_, _| Ok(Some(0)))
                    .map(drop))
                .is_ok()
        );
    }
}
