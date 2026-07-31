const INSERTION_EMISSION: f64 = 0.25;
const MISMATCH_EMISSION: f64 = 0.333_333_333_33;

pub(crate) fn score(
    reference: &[u8],
    query: &[u8],
    qualities: &[u8],
    gap_open: f64,
    gap_extension: f64,
    bandwidth: usize,
) -> Option<i32> {
    if reference.is_empty() || query.is_empty() || qualities.len() != query.len() {
        return None;
    }

    let mut bandwidth = reference.len().max(query.len()).min(bandwidth);
    bandwidth = bandwidth.max(reference.len().abs_diff(query.len()));
    let band_width = bandwidth.checked_mul(2)?.checked_add(1)?;
    let row_width = if band_width < reference.len() {
        band_width.checked_mul(3)?.checked_add(6)?
    } else {
        reference.len().checked_mul(3)?.checked_add(6)?
    };
    let rows = query.len().checked_add(1)?;
    let mut forward = vec![0.0; rows.checked_mul(row_width)?];
    let mut scales = vec![0.0; query.len().checked_add(2)?];
    let qualities = qualities
        .iter()
        .map(|&quality| 10.0f32.powf(-f32::from(quality) / 10.0))
        .collect::<Vec<_>>();

    let start_match = 1.0 / (2.0 * query.len() as f64 + 2.0);
    let start_insertion = start_match;
    let transitions = [
        (1.0 - gap_open - gap_open) * (1.0 - start_match),
        gap_open * (1.0 - start_match),
        gap_open * (1.0 - start_match),
        (1.0 - gap_extension) * (1.0 - start_insertion),
        gap_extension * (1.0 - start_insertion),
        0.0,
        1.0 - gap_extension,
        0.0,
        gap_extension,
    ];
    let begin_match = (1.0 - gap_open) / reference.len() as f64;
    let begin_insertion = gap_open / reference.len() as f64;

    let index = band_index(0, 0, bandwidth)?;
    forward[index] = 1.0;
    scales[0] = 1.0;

    let first_row = row_width;
    let end = reference.len().min(bandwidth + 1);
    let mut sum = 0.0;
    for reference_index in 1..=end {
        let emission = emission(reference[reference_index - 1], query[0], qualities[0]);
        let column = band_index(1, reference_index, bandwidth)?;
        forward[first_row + column] = emission * begin_match;
        forward[first_row + column + 1] = INSERTION_EMISSION * begin_insertion;
        sum += forward[first_row + column] + forward[first_row + column + 1];
    }
    scales[1] = sum;

    for query_index in 2..=query.len() {
        let row = query_index * row_width;
        let previous_row = (query_index - 1) * row_width;
        let begin = 1usize.max(query_index.saturating_sub(bandwidth));
        let end = reference.len().min(query_index.saturating_add(bandwidth));
        let scale = 1.0 / scales[query_index - 1];
        let mut sum = 0.0;

        for reference_index in begin..=end {
            let column = band_index(query_index, reference_index, bandwidth)?;
            let diagonal = band_index(query_index - 1, reference_index - 1, bandwidth)?;
            let vertical = band_index(query_index - 1, reference_index, bandwidth)?;
            let horizontal = band_index(query_index, reference_index - 1, bandwidth)?;
            let emission = emission(
                reference[reference_index - 1],
                query[query_index - 1],
                qualities[query_index - 1],
            );
            let matched = emission
                * (transitions[0] * scale * forward[previous_row + diagonal]
                    + transitions[3] * scale * forward[previous_row + diagonal + 1]
                    + transitions[6] * scale * forward[previous_row + diagonal + 2]);
            let inserted = INSERTION_EMISSION
                * (transitions[1] * scale * forward[previous_row + vertical]
                    + transitions[4] * scale * forward[previous_row + vertical + 1]);
            let deleted = transitions[2] * forward[row + horizontal]
                + transitions[8] * forward[row + horizontal + 2];
            forward[row + column] = matched;
            forward[row + column + 1] = inserted;
            forward[row + column + 2] = deleted;
            sum += matched + inserted + deleted;
        }
        scales[query_index] = sum;
    }

    let last_row = query.len() * row_width;
    let scale = 1.0 / scales[query.len()];
    let mut terminal = 0.0;
    for reference_index in 1..=reference.len() {
        let Some(column) = band_index(query.len(), reference_index, bandwidth) else {
            continue;
        };
        if column < 3 || column >= row_width {
            continue;
        }
        terminal += scale * forward[last_row + column] * start_match
            + scale * forward[last_row + column + 1] * start_insertion;
    }
    scales[query.len() + 1] = terminal;

    let mut product = 1.0;
    let mut phred = 0.0;
    for scale in scales {
        product *= scale;
        if product < 1e-100 {
            phred += -4.343 * product.ln();
            product = 1.0;
        }
    }
    phred += -4.343 * (product * reference.len() as f64 * query.len() as f64).ln();
    if phred.is_finite() && phred >= f64::from(i32::MIN) && phred <= f64::from(i32::MAX) {
        Some((phred + 0.499) as i32)
    } else {
        None
    }
}

fn emission(reference: u8, query: u8, quality: f32) -> f64 {
    if reference > 3 || query > 3 {
        1.0
    } else if reference == query {
        1.0 - f64::from(quality)
    } else {
        f64::from(quality) * MISMATCH_EMISSION
    }
}

fn band_index(query_index: usize, reference_index: usize, bandwidth: usize) -> Option<usize> {
    let start = query_index.saturating_sub(bandwidth);
    reference_index
        .checked_add(1)?
        .checked_sub(start)?
        .checked_mul(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_htslib_1_24_likelihood_scores() {
        let reference = [0, 1, 2, 3, 0, 1, 2, 3];
        let exact = [0, 1, 2, 3, 0, 1, 2, 3];
        let insertion = [0, 1, 2, 3, 3, 0, 1, 2, 3];
        let deletion = [0, 1, 2, 0, 1, 2, 3];
        let mismatch = [0, 1, 2, 3, 3, 1, 2, 3];
        let ambiguous = [0, 1, 4, 3, 0, 1, 2, 3];

        assert_eq!(score(&reference, &exact, &[30; 8], 1e-4, 1e-2, 10), Some(5));
        assert_eq!(
            score(&reference, &insertion, &[30; 9], 1e-4, 1e-2, 4),
            Some(48)
        );
        assert_eq!(
            score(&reference, &deletion, &[30; 7], 1e-4, 1e-2, 4),
            Some(45)
        );
        assert_eq!(
            score(&reference, &mismatch, &[30; 8], 1e-4, 1e-2, 10),
            Some(40)
        );
        assert_eq!(
            score(&reference, &ambiguous, &[30; 8], 1e-4, 1e-2, 10),
            Some(5)
        );
    }
}
