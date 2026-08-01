use crate::{
    SiteAnnotations,
    annotation::{fisher_exact, ln_gamma},
};

#[derive(Clone, Debug, PartialEq)]
pub struct CalledAnnotations {
    pileup: SiteAnnotations,
    strand_depths: [u32; 4],
    mapping_quality: u32,
    bias_probabilities: Option<[f32; 4]>,
}

impl CalledAnnotations {
    pub(crate) fn consensus(pileup: &SiteAnnotations) -> Self {
        let auxiliary = pileup.auxiliary();
        let depth = auxiliary[..4].iter().sum::<f32>();
        let mapping_quality = if depth == 0.0 {
            0
        } else {
            (((auxiliary[9] + auxiliary[11]) / depth).sqrt() + 0.499) as u32
        };
        Self::new(pileup, mapping_quality)
    }

    pub(crate) fn multiallelic(pileup: &SiteAnnotations) -> Self {
        let auxiliary = pileup.auxiliary();
        let depth = auxiliary[..4].iter().sum::<f32>();
        let mapping_quality = if depth == 0.0 {
            0
        } else {
            ((auxiliary[8] + auxiliary[10]) / depth) as u32
        };
        Self::new(pileup, mapping_quality)
    }

    fn new(pileup: &SiteAnnotations, mapping_quality: u32) -> Self {
        let auxiliary = pileup.auxiliary();
        Self {
            pileup: pileup.clone(),
            strand_depths: std::array::from_fn(|index| auxiliary[index] as u32),
            mapping_quality,
            bias_probabilities: bias_probabilities(auxiliary),
        }
    }

    pub fn pileup(&self) -> &SiteAnnotations {
        &self.pileup
    }

    pub fn strand_depths(&self) -> &[u32; 4] {
        &self.strand_depths
    }

    pub fn mapping_quality(&self) -> u32 {
        self.mapping_quality
    }

    pub fn bias_probabilities(&self) -> Option<&[f32; 4]> {
        self.bias_probabilities.as_ref()
    }
}

fn bias_probabilities(auxiliary: &[f32; 16]) -> Option<[f32; 4]> {
    let reference_count = auxiliary[0] + auxiliary[1];
    let alternate_count = auxiliary[2] + auxiliary[3];
    if reference_count <= 0.0 || alternate_count <= 0.0 {
        return None;
    }
    let (_, strand) = fisher_exact(
        auxiliary[0] as u64,
        auxiliary[1] as u64,
        auxiliary[2] as u64,
        auxiliary[3] as u64,
    );
    Some([
        strand as f32,
        one_sided_t_test(reference_count, alternate_count, &auxiliary[4..8]) as f32,
        one_sided_t_test(reference_count, alternate_count, &auxiliary[8..12]) as f32,
        one_sided_t_test(reference_count, alternate_count, &auxiliary[12..16]) as f32,
    ])
}

fn one_sided_t_test(reference_count: f32, alternate_count: f32, moments: &[f32]) -> f64 {
    let n1 = f64::from(reference_count);
    let n2 = f64::from(alternate_count);
    if n1 == 0.0 || n2 == 0.0 || n1 + n2 < 3.0 {
        return 1.0;
    }
    let reference_mean = f64::from(moments[0]) / n1;
    let alternate_mean = f64::from(moments[2]) / n2;
    if reference_mean <= alternate_mean {
        return 1.0;
    }
    let variance = ((f64::from(moments[1]) - n1 * reference_mean * reference_mean)
        + (f64::from(moments[3]) - n2 * alternate_mean * alternate_mean))
        / (n1 + n2 - 2.0)
        * (1.0 / n1 + 1.0 / n2);
    let statistic = (reference_mean - alternate_mean) / variance.sqrt();
    if statistic < 0.0 {
        return 1.0;
    }
    let freedom = n1 + n2 - 2.0;
    0.5 * regularized_incomplete_beta(
        0.5 * freedom,
        0.5,
        freedom / (freedom + statistic * statistic),
    )
}

fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        incomplete_beta_fraction(a, b, x)
    } else {
        1.0 - incomplete_beta_fraction(b, a, 1.0 - x)
    }
}

fn incomplete_beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    let mut fraction = 1.0;
    let mut numerator = fraction;
    let mut denominator = 0.0;
    for index in 1..200 {
        let half = index >> 1;
        let half = f64::from(half);
        let coefficient = if index & 1 == 1 {
            -(a + half) * (a + b + half) * x / ((a + 2.0 * half) * (a + 2.0 * half + 1.0))
        } else {
            half * (b - half) * x / ((a + 2.0 * half - 1.0) * (a + 2.0 * half))
        };
        denominator = (1.0 + coefficient * denominator).max(1e-290);
        numerator = (1.0 + coefficient / numerator).max(1e-290);
        denominator = denominator.recip();
        let delta = numerator * denominator;
        fraction *= delta;
        if (delta - 1.0).abs() < 1e-14 {
            break;
        }
    }
    (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp()
        / a
        / fraction
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::SiteAnnotationValues;

    #[test]
    fn pv4_matches_bcftools_1_24() {
        let annotations = SiteAnnotations::new(SiteAnnotationValues {
            raw_depth: 3,
            auxiliary: [
                1.0, 0.0, 2.0, 0.0, 40.0, 1600.0, 80.0, 3200.0, 60.0, 3600.0, 120.0, 7200.0, 0.0,
                0.0, 0.0, 0.0,
            ],
            variant_distance_bias: None,
            read_position_bias: None,
            mapping_quality_bias: None,
            base_quality_bias: None,
            mapping_quality_strand_bias: None,
            mismatch_bias: None,
            soft_clip_bias: None,
            strand_bias: None,
            segregation_bias: None,
            zero_mapping_quality_fraction: 0.0,
            average_mismatches: None,
        })
        .unwrap();
        let consensus = CalledAnnotations::consensus(&annotations);
        assert_eq!(consensus.strand_depths(), &[1, 0, 2, 0]);
        assert_eq!(consensus.mapping_quality(), 60);
        assert_eq!(consensus.bias_probabilities(), Some(&[1.0; 4]));
        assert_eq!(CalledAnnotations::multiallelic(&annotations), consensus);
    }

    #[test]
    fn caller_models_use_their_bcftools_mapping_quality_definitions() {
        let annotations = SiteAnnotations::new(SiteAnnotationValues {
            raw_depth: 7,
            auxiliary: [
                3.0, 1.0, 1.0, 2.0, 160.0, 6400.0, 60.0, 1200.0, 200.0, 10000.0, 30.0, 300.0, 40.0,
                400.0, 6.0, 12.0,
            ],
            variant_distance_bias: None,
            read_position_bias: None,
            mapping_quality_bias: None,
            base_quality_bias: None,
            mapping_quality_strand_bias: None,
            mismatch_bias: None,
            soft_clip_bias: None,
            strand_bias: None,
            segregation_bias: None,
            zero_mapping_quality_fraction: 0.0,
            average_mismatches: None,
        })
        .unwrap();
        assert_eq!(
            CalledAnnotations::consensus(&annotations).mapping_quality(),
            38
        );
        assert_eq!(
            CalledAnnotations::multiallelic(&annotations).mapping_quality(),
            32
        );
    }
}
