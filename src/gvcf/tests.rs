use crate::{Allele, LikelihoodSite, MultiallelicCaller, Ploidy, SampleEvidence, SampleLikelihood};

use super::*;

fn reference_site(position: u64, depth: u32) -> CalledSite {
    let site = LikelihoodSite::new(
        0,
        position,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"<*>"[..]).unwrap()],
        [1.0, 0.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            [0, 30, 90],
            SampleEvidence::new(depth, [depth, 0], [depth.saturating_mul(30), 0]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap();
    MultiallelicCaller::default().call(&site).unwrap()
}

#[test]
fn groups_contiguous_reference_sites_by_minimum_depth() {
    let mut blocker = GvcfBlocker::new([5, 15]).unwrap();
    let mut output = Vec::new();
    for site in [
        reference_site(0, 10),
        reference_site(1, 8),
        reference_site(2, 20),
        reference_site(3, 18),
    ] {
        blocker
            .push(site, |site| {
                output.push(site);
                Ok(())
            })
            .unwrap();
    }
    blocker
        .finish(|site| {
            output.push(site);
            Ok(())
        })
        .unwrap();

    assert_eq!(output.len(), 2);
    assert_eq!(output[0].position(), 0);
    assert_eq!(output[0].gvcf().unwrap().end_position(), Some(1));
    assert_eq!(output[0].gvcf().unwrap().minimum_depth(), 8);
    assert_eq!(output[0].samples()[0].evidence().depth(), 8);
    assert_eq!(output[1].position(), 2);
    assert_eq!(output[1].gvcf().unwrap().end_position(), Some(3));
    assert_eq!(output[1].gvcf().unwrap().minimum_depth(), 18);
}

#[test]
fn low_depth_and_gaps_break_blocks() {
    let mut blocker = GvcfBlocker::new([5]).unwrap();
    let mut output = Vec::new();
    for site in [
        reference_site(0, 1),
        reference_site(1, 5),
        reference_site(3, 6),
    ] {
        blocker
            .push(site, |site| {
                output.push(site);
                Ok(())
            })
            .unwrap();
    }
    blocker
        .finish(|site| {
            output.push(site);
            Ok(())
        })
        .unwrap();
    assert_eq!(output.len(), 3);
    assert!(!output[0].gvcf().unwrap().is_collapsed());
    assert_eq!(output[0].gvcf().unwrap().minimum_depth(), 1);
    assert!(output[1].gvcf().unwrap().is_collapsed());
    assert_eq!(output[1].gvcf().unwrap().end_position(), None);
    assert_eq!(output[2].position(), 3);
}

#[test]
fn rejects_invalid_thresholds_and_order() {
    for thresholds in [vec![], vec![5, 5], vec![10, 5]] {
        assert!(matches!(
            GvcfBlocker::new(thresholds),
            Err(CallError::InvalidGvcfThresholds)
        ));
    }
    let mut blocker = GvcfBlocker::new([1]).unwrap();
    blocker.push(reference_site(2, 1), |_| Ok(())).unwrap();
    assert_eq!(
        blocker.push(reference_site(1, 1), |_| Ok(())),
        Err(CallError::InvalidGvcfOrder)
    );
}

#[test]
fn retains_the_lowest_nonreference_likelihood_pair() {
    let mut sites = [
        reference_site(0, 10),
        reference_site(1, 8),
        reference_site(2, 9),
    ];
    for (site, likelihoods) in sites.iter_mut().zip([[0, 10, 30], [0, 5, 50], [0, 5, 40]]) {
        site.samples[0].phred_likelihoods = Some(likelihoods.into());
    }
    let mut blocker = GvcfBlocker::new([5]).unwrap();
    let mut output = Vec::new();
    for site in sites {
        blocker
            .push(site, |site| {
                output.push(site);
                Ok(())
            })
            .unwrap();
    }
    blocker
        .finish(|site| {
            output.push(site);
            Ok(())
        })
        .unwrap();

    assert_eq!(output.len(), 1);
    assert_eq!(
        output[0].samples()[0].phred_likelihoods(),
        Some(&[0, 5, 40][..])
    );
    assert_eq!(output[0].samples()[0].evidence().depth(), 8);
}

#[test]
fn rejects_incompatible_block_samples_and_likelihoods() {
    let mut invalid_first = reference_site(0, 1);
    invalid_first.samples[0].phred_likelihoods = Some([0, 1].into());
    let mut blocker = GvcfBlocker::new([1]).unwrap();
    assert_eq!(
        blocker.push(invalid_first, |_| Ok(())),
        Err(CallError::InvalidGvcfLikelihoods)
    );

    let mut blocker = GvcfBlocker::new([1]).unwrap();
    blocker.push(reference_site(0, 1), |_| Ok(())).unwrap();
    let mut two_samples = reference_site(1, 1);
    two_samples.samples = vec![two_samples.samples[0].clone(); 2].into_boxed_slice();
    assert_eq!(
        blocker.push(two_samples, |_| Ok(())),
        Err(CallError::GvcfSampleCountMismatch)
    );

    let mut first = reference_site(0, 1);
    first.samples[0].phred_likelihoods = Some([0, 1, 2].into());
    let mut second = reference_site(1, 1);
    second.samples[0].phred_likelihoods = Some([0, 1].into());
    let mut blocker = GvcfBlocker::new([1]).unwrap();
    blocker.push(first, |_| Ok(())).unwrap();
    assert_eq!(
        blocker.push(second, |_| Ok(())),
        Err(CallError::InvalidGvcfLikelihoods)
    );
}

#[test]
fn excludes_a_duplicate_nonreference_coordinate_from_the_preceding_block() {
    let mut variant = reference_site(1, 10);
    variant.alternates = [Allele::new(&b"G"[..]).unwrap()].into();
    variant.allele_counts = [0, 2].into();
    variant.samples[0].genotype = Some([1, 1].into());
    let mut blocker = GvcfBlocker::new([1]).unwrap();
    let mut output = Vec::new();
    for site in [reference_site(0, 10), reference_site(1, 10), variant] {
        blocker
            .push(site, |site| {
                output.push(site);
                Ok(())
            })
            .unwrap();
    }
    blocker
        .finish(|site| {
            output.push(site);
            Ok(())
        })
        .unwrap();

    assert_eq!(output.len(), 2);
    assert_eq!(output[0].gvcf().unwrap().end_position(), None);
    assert_eq!(output[1].position(), 1);
    assert!(output[1].is_variant());
}
