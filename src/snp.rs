use rsomics_pileup::Column;

use crate::{BaseObservation, CallError, ErrorModel, LikelihoodMatrix, Nucleotide, Result};

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
        })
    }
}

impl Default for SnpLikelihoodConfig {
    fn default() -> Self {
        Self {
            minimum_base_quality: 1,
            maximum_base_quality: 60,
            neighboring_quality_delta: 30,
            mapping_quality_cap: 60,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SnpEvidence {
    likelihoods: LikelihoodMatrix,
    allele_depths: [u32; 5],
    quality_sums: [u32; 5],
}

impl SnpEvidence {
    pub fn from_column(
        column: &Column<'_>,
        reference_base: Nucleotide,
        config: SnpLikelihoodConfig,
        model: &ErrorModel,
        observations: &mut Vec<BaseObservation>,
    ) -> Result<Self> {
        observations.clear();
        let mut allele_depths = [0; 5];
        let mut quality_sums = [0; 5];

        for entry in column.entries() {
            let projection = entry.projection();
            if projection.is_deletion || projection.is_reference_skip {
                continue;
            }
            let record = entry.record();
            let qualities = record.quality_scores();
            let qpos = projection.qpos;
            let mut base_quality = u16::from(qualities[qpos]);
            if qpos > 0 {
                base_quality = base_quality.min(
                    u16::from(qualities[qpos - 1]) + u16::from(config.neighboring_quality_delta),
                );
            }
            if qpos + 1 < qualities.len() {
                base_quality = base_quality.min(
                    u16::from(qualities[qpos + 1]) + u16::from(config.neighboring_quality_delta),
                );
            }
            if base_quality < u16::from(config.minimum_base_quality) {
                continue;
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
            allele_depths[base as usize] += 1;
            quality_sums[base as usize] += u32::from(effective_quality);
            observations.push(BaseObservation::new(
                base,
                effective_quality,
                record.flags() & 0x10 != 0,
            ));
        }

        Ok(Self {
            likelihoods: model.calculate(observations)?,
            allele_depths,
            quality_sums,
        })
    }

    pub fn likelihoods(&self) -> &LikelihoodMatrix {
        &self.likelihoods
    }

    pub fn allele_depths(&self) -> &[u32; 5] {
        &self.allele_depths
    }

    pub fn quality_sums(&self) -> &[u32; 5] {
        &self.quality_sums
    }
}

#[cfg(test)]
mod tests {
    use rsomics_bamio::raw::RawRecord;
    use rsomics_pileup::{PileupEngine, PileupOptions};

    use super::*;

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

        let model = ErrorModel::default();
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
                            &model,
                            &mut observations,
                        )
                        .unwrap(),
                    );
                }
                Ok::<_, ()>(())
            })
            .unwrap();

        let evidence = evidence.unwrap();
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
}
