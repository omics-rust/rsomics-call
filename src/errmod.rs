use std::sync::OnceLock;

use crate::{CallError, Result};

const MAX_DEPTH: usize = 255;
const QUALITY_LEVELS: usize = 64;
const BASES: usize = 5;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    N,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BaseObservation(u16);

impl BaseObservation {
    pub fn new(base: Nucleotide, quality: u8, reverse: bool) -> Self {
        Self(u16::from(quality.clamp(4, 63)) << 5 | u16::from(reverse) << 4 | base as u16)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LikelihoodMatrix([f32; BASES * BASES]);

impl LikelihoodMatrix {
    pub fn get(&self, first: Nucleotide, second: Nucleotide) -> f32 {
        self.0[first as usize * BASES + second as usize]
    }

    pub fn as_slice(&self) -> &[f32] {
        &self.0
    }
}

pub struct ErrorModel {
    correlation: [f64; MAX_DEPTH + 1],
    log_binomial: Box<[f64]>,
    heterozygous: Box<[f64]>,
    beta: Box<[OnceLock<Box<[f64]>>]>,
    random: HtsRandom,
}

impl Default for ErrorModel {
    fn default() -> Self {
        Self::new()
    }
}

impl ErrorModel {
    pub fn new() -> Self {
        Self::with_dependency_correlation(0.17).unwrap()
    }

    pub fn with_dependency_correlation(dependency_correlation: f64) -> Result<Self> {
        Self::with_dependency_correlation_and_seed(dependency_correlation, 0)
    }

    pub fn with_random_seed(seed: i32) -> Self {
        Self::with_dependency_correlation_and_seed(0.17, seed).unwrap()
    }

    pub fn with_dependency_correlation_and_seed(
        dependency_correlation: f64,
        seed: i32,
    ) -> Result<Self> {
        if !dependency_correlation.is_finite() || !(0.0..=1.0).contains(&dependency_correlation) {
            return Err(CallError::InvalidDependencyCorrelation);
        }
        let mut correlation = [0.0; MAX_DEPTH + 1];
        correlation[0] = 1.0;
        for (depth, value) in correlation.iter_mut().enumerate().skip(1) {
            *value = (1.0 - dependency_correlation).powi(depth as i32) * 0.97 + 0.03;
        }

        let mut log_binomial = vec![0.0; 256 * 256];
        for trials in 1..=MAX_DEPTH {
            let total = log_gamma(trials as f64 + 1.0);
            for successes in 1..=trials {
                log_binomial[trials << 8 | successes] = total
                    - log_gamma(successes as f64 + 1.0)
                    - log_gamma((trials - successes) as f64 + 1.0);
            }
        }
        let mut heterozygous = vec![0.0; 256 * 256];
        for trials in 0..=MAX_DEPTH {
            for successes in 0..=MAX_DEPTH {
                heterozygous[trials << 8 | successes] =
                    log_binomial[trials << 8 | successes] - std::f64::consts::LN_2 * trials as f64;
            }
        }

        Ok(Self {
            correlation,
            log_binomial: log_binomial.into_boxed_slice(),
            heterozygous: heterozygous.into_boxed_slice(),
            beta: std::iter::repeat_with(OnceLock::new)
                .take(MAX_DEPTH + 1)
                .collect(),
            random: HtsRandom::new(seed),
        })
    }

    pub fn calculate(&mut self, observations: &mut [BaseObservation]) -> Result<LikelihoodMatrix> {
        if observations.len() > MAX_DEPTH {
            self.random.shuffle(observations);
        }
        let depth = observations.len().min(MAX_DEPTH);
        let observations = &mut observations[..depth];
        let depth = observations.len();
        let mut output = [0.0; BASES * BASES];
        if depth == 0 {
            return Ok(LikelihoodMatrix(output));
        }

        observations.sort_unstable_by_key(|observation| observation.0);
        let beta = self.beta[depth].get_or_init(|| self.build_beta(depth));
        let mut base_sums = [0.0; 16];
        let mut base_counts = [0usize; 16];
        let mut strand_counts = [0usize; 32];

        for observation in observations.iter().rev() {
            let quality = usize::from(observation.0 >> 5).clamp(4, 63);
            let base_strand = usize::from(observation.0 & 0x1f);
            let base = usize::from(observation.0 & 0x0f);
            base_sums[base] += self.correlation[strand_counts[base_strand]]
                * beta[quality << 8 | base_counts[base]];
            base_counts[base] += 1;
            strand_counts[base_strand] += 1;
        }

        for first in 0..BASES {
            let mut homozygous = 0.0f32;
            let mut outside_count = 0;
            for base in 0..BASES {
                if base != first {
                    homozygous += base_sums[base] as f32;
                    outside_count += base_counts[base];
                }
            }
            if outside_count > 0 {
                output[first * BASES + first] = homozygous;
            }

            for second in first + 1..BASES {
                let selected_count = base_counts[first] + base_counts[second];
                let mut outside_sum = 0.0f32;
                for (base, sum) in base_sums.iter().enumerate().take(BASES) {
                    if base != first && base != second {
                        outside_sum += *sum as f32;
                    }
                }
                let value = (-4.343 * self.heterozygous[selected_count << 8 | base_counts[second]]
                    + f64::from(outside_sum)) as f32;
                output[first * BASES + second] = value.max(0.0);
                output[second * BASES + first] = value.max(0.0);
            }
        }

        Ok(LikelihoodMatrix(output))
    }

    fn build_beta(&self, depth: usize) -> Box<[f64]> {
        let mut beta = vec![0.0; QUALITY_LEVELS * 256];
        for quality in 1..QUALITY_LEVELS {
            let error = 10f64.powf(-(quality as f64) / 10.0);
            let log_error = error.ln();
            let log_correct = (1.0 - error).ln();
            let mut cumulative = self.log_binomial[depth << 8 | depth] + depth as f64 * log_error;
            beta[quality << 8 | depth] = f64::INFINITY;
            for matches in (0..depth).rev() {
                let previous = cumulative;
                cumulative = previous
                    + (self.log_binomial[depth << 8 | matches]
                        + matches as f64 * log_error
                        + (depth - matches) as f64 * log_correct
                        - previous)
                        .exp()
                        .ln_1p();
                beta[quality << 8 | matches] =
                    -10.0 / std::f64::consts::LN_10 * (previous - cumulative);
            }
        }
        beta.into_boxed_slice()
    }
}

struct HtsRandom {
    state: u64,
}

impl HtsRandom {
    fn new(seed: i32) -> Self {
        Self {
            state: u64::from(seed as u32) << 16 | 0x330e,
        }
    }

    fn next_f64(&mut self) -> f64 {
        self.state =
            (0x5deece66du64.wrapping_mul(self.state).wrapping_add(0xb)) & ((1u64 << 48) - 1);
        self.state as f64 / (1u64 << 48) as f64
    }

    fn shuffle<T>(&mut self, values: &mut [T]) {
        for end in (2..=values.len()).rev() {
            let index = (self.next_f64() * end as f64) as usize;
            values.swap(index, end - 1);
        }
    }
}

fn log_gamma(value: f64) -> f64 {
    const COEFFICIENTS: [f64; 6] = [
        76.180_091_729_471_46,
        -86.505_320_329_416_77,
        24.014_098_240_830_91,
        -1.231_739_572_450_155,
        1.208_650_973_866_179e-3,
        -5.395_239_384_953e-6,
    ];
    let shifted = value - 1.0;
    let series = COEFFICIENTS
        .iter()
        .enumerate()
        .fold(1.000_000_000_190_015, |sum, (index, coefficient)| {
            sum + coefficient / (shifted + index as f64 + 1.0)
        });
    let scale = shifted + 5.5;
    (shifted + 0.5) * scale.ln() - scale + (2.506_628_274_631_001 * series).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_matrix(actual: &LikelihoodMatrix, expected: [f32; 25]) {
        for (index, (&actual, expected)) in actual.as_slice().iter().zip(expected).enumerate() {
            assert!(
                (actual - expected).abs() < 0.0005,
                "matrix index {index}: {actual} != {expected}"
            );
        }
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn matches_htslib_1_24_balanced_observations() {
        let mut observations = Vec::new();
        observations.extend((0..10).map(|_| BaseObservation::new(Nucleotide::A, 40, false)));
        observations.extend((0..10).map(|_| BaseObservation::new(Nucleotide::G, 40, true)));
        let matrix = ErrorModel::default().calculate(&mut observations).unwrap();
        assert_matrix(
            &matrix,
            [
                168.366653, 198.470032, 7.54010963, 198.470032, 198.470032, 198.470032, 336.733307,
                198.470032, 336.733307, 336.733307, 7.54010963, 198.470032, 168.366653, 198.470032,
                198.470032, 198.470032, 336.733307, 198.470032, 336.733307, 336.733307, 198.470032,
                336.733307, 198.470032, 336.733307, 336.733307,
            ],
        );
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn matches_htslib_1_24_mixed_qualities() {
        let mut observations = [
            BaseObservation::new(Nucleotide::A, 60, false),
            BaseObservation::new(Nucleotide::A, 50, false),
            BaseObservation::new(Nucleotide::A, 40, false),
            BaseObservation::new(Nucleotide::A, 30, false),
            BaseObservation::new(Nucleotide::A, 20, false),
            BaseObservation::new(Nucleotide::A, 10, false),
            BaseObservation::new(Nucleotide::A, 8, false),
            BaseObservation::new(Nucleotide::A, 4, false),
            BaseObservation::new(Nucleotide::C, 35, true),
            BaseObservation::new(Nucleotide::C, 25, true),
            BaseObservation::new(Nucleotide::C, 15, true),
            BaseObservation::new(Nucleotide::G, 45, false),
        ];
        let matrix = ErrorModel::default().calculate(&mut observations).unwrap();
        assert_matrix(
            &matrix,
            [
                79.9680862, 45.147541, 63.3096428, 104.050789, 104.050789, 45.147541, 175.383667,
                147.195404, 184.414688, 184.414688, 63.3096428, 147.195404, 186.933868, 189.944214,
                189.944214, 104.050789, 184.414688, 189.944214, 221.142807, 221.142807, 104.050789,
                184.414688, 189.944214, 221.142807, 221.142807,
            ],
        );
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn matches_htslib_1_24_strand_correlation() {
        let mut observations = [BaseObservation::new(Nucleotide::A, 40, false); 8];
        let matrix = ErrorModel::default().calculate(&mut observations).unwrap();
        assert_matrix(
            &matrix,
            [
                0.0, 24.0827065, 24.0827065, 24.0827065, 24.0827065, 24.0827065, 176.440491,
                176.440491, 176.440491, 176.440491, 24.0827065, 176.440491, 176.440491, 176.440491,
                176.440491, 24.0827065, 176.440491, 176.440491, 176.440491, 176.440491, 24.0827065,
                176.440491, 176.440491, 176.440491, 176.440491,
            ],
        );
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn matches_htslib_1_24_deep_sampling_stream() {
        fn observations() -> Vec<BaseObservation> {
            (0..300)
                .map(|index| {
                    BaseObservation::new(
                        if index % 3 == 0 {
                            Nucleotide::G
                        } else {
                            Nucleotide::A
                        },
                        20 + (index % 41) as u8,
                        index % 5 == 0,
                    )
                })
                .collect()
        }

        let mut model = ErrorModel::with_random_seed(0);
        let mut first = observations();
        assert_matrix(
            &model.calculate(&mut first).unwrap(),
            [
                504.332642, 1028.13147, 88.087204, 1028.13147, 1028.13147, 1028.13147, 1161.66272,
                901.167419, 1161.66272, 1161.66272, 88.087204, 901.167419, 657.330017, 901.167419,
                901.167419, 1028.13147, 1161.66272, 901.167419, 1161.66272, 1161.66272, 1028.13147,
                1161.66272, 901.167419, 1161.66272, 1161.66272,
            ],
        );
        let mut second = observations();
        assert_matrix(
            &model.calculate(&mut second).unwrap(),
            [
                506.537384, 1039.36719, 98.3605499, 1039.36719, 1039.36719, 1039.36719, 1178.02258,
                906.291626, 1178.02258, 1178.02258, 98.3605499, 906.291626, 671.485229, 906.291626,
                906.291626, 1039.36719, 1178.02258, 906.291626, 1178.02258, 1178.02258, 1039.36719,
                1178.02258, 906.291626, 1178.02258, 1178.02258,
            ],
        );
    }

    #[test]
    fn correlation_parameter_is_checked() {
        assert!(ErrorModel::with_dependency_correlation(0.0).is_ok());
        assert!(ErrorModel::with_dependency_correlation(1.0).is_ok());
        assert!(matches!(
            ErrorModel::with_dependency_correlation(f64::NAN),
            Err(CallError::InvalidDependencyCorrelation)
        ));
        assert!(matches!(
            ErrorModel::with_dependency_correlation(1.1),
            Err(CallError::InvalidDependencyCorrelation)
        ));
    }
}
