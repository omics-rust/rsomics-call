use super::*;
use crate::{Allele, Ploidy};

fn sample(phred_likelihoods: [u32; 6], depths: [u32; 3], sums: [u32; 3]) -> SampleLikelihood {
    SampleLikelihood::observed(
        Ploidy::new(2).unwrap(),
        phred_likelihoods,
        SampleEvidence::new(1, depths, sums).unwrap(),
    )
    .unwrap()
}

fn two_sample_site() -> LikelihoodSite {
    LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [
            Allele::new(&b"G"[..]).unwrap(),
            Allele::new(&b"<*>"[..]).unwrap(),
        ],
        [1.0, 1.0, 0.0],
        [
            sample([0, 3, 40, 3, 40, 40], [1, 0, 0], [40, 0, 0]),
            sample([40, 3, 0, 40, 3, 40], [0, 1, 0], [0, 40, 0]),
        ],
    )
    .unwrap()
}

#[test]
fn matches_bcftools_1_24_two_sample_multiallelic_call() {
    let site = two_sample_site();

    let called = MultiallelicCaller::default().call(&site).unwrap();

    assert_eq!(called.reference().as_bytes(), b"A");
    assert_eq!(called.alternates()[0].as_bytes(), b"G");
    assert_eq!(called.allele_counts(), &[2, 2]);
    assert_eq!(called.allele_number(), 4);
    assert!(called.is_variant());
    assert!((called.quality().unwrap() - 7.822_08).abs() < 1e-5);
    assert_eq!(called.samples()[0].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(3));
    assert_eq!(called.samples()[1].genotype_quality(), Some(3));
    assert_eq!(
        called.samples()[0].phred_likelihoods(),
        Some(&[0, 3, 40][..])
    );
    assert_eq!(
        called.samples()[1].phred_likelihoods(),
        Some(&[40, 3, 0][..])
    );
    let expected = [
        [0.499_382, 0.500_568, 4.993_82e-5],
        [4.993_82e-5, 0.500_568, 0.499_382],
    ];
    for (sample, expected) in called.samples().iter().zip(expected) {
        for (&observed, expected) in sample
            .genotype_probabilities()
            .unwrap()
            .iter()
            .zip(expected)
        {
            assert!((observed - expected).abs() < 5e-7);
        }
    }
}

#[test]
fn rejects_non_diploid_likelihoods_and_invalid_prior() {
    assert_eq!(
        MultiallelicCallerConfig::new(f64::NAN),
        Err(CallError::InvalidMutationRate)
    );
    let site = LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"G"[..]).unwrap()],
        [1.0, 1.0],
        [SampleLikelihood::observed(
            Ploidy::new(1).unwrap(),
            [0, 40],
            SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap();
    assert_eq!(
        MultiallelicCaller::default().call(&site),
        Err(CallError::UnsupportedLikelihoodPloidy)
    );
}

#[test]
fn matches_bcftools_1_24_haploid_call() {
    let site = two_sample_site();
    let haploid = CallPloidy::Haploid;

    let called = MultiallelicCaller::default()
        .call_with_ploidies(&site, &[haploid, haploid], 2)
        .unwrap();

    assert_eq!(called.allele_counts(), &[1, 1]);
    assert_eq!(called.allele_number(), 2);
    assert!((called.quality().unwrap() - 5.7423).abs() < 1e-4);
    assert_eq!(called.samples()[0].genotype(), Some(&[0][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(40));
    assert_eq!(called.samples()[1].genotype_quality(), Some(40));
    assert_eq!(called.samples()[0].phred_likelihoods(), Some(&[0, 40][..]));
    assert_eq!(called.samples()[1].phred_likelihoods(), Some(&[40, 0][..]));
    let expected = [[0.9999, 9.999e-5], [9.999e-5, 0.9999]];
    for (sample, expected) in called.samples().iter().zip(expected) {
        for (&observed, expected) in sample
            .genotype_probabilities()
            .unwrap()
            .iter()
            .zip(expected)
        {
            assert!((observed - expected).abs() < 5e-7);
        }
    }
}

#[test]
fn matches_bcftools_1_24_mixed_ploidy_call() {
    let site = two_sample_site();
    let haploid = CallPloidy::Haploid;
    let diploid = CallPloidy::Diploid;

    assert_eq!(
        MultiallelicCaller::default().call_with_ploidies(&site, &[haploid, diploid], 2),
        Err(CallError::InvalidPriorChromosomeCount)
    );
    assert_eq!(
        MultiallelicCaller::default().call_with_ploidies(&site, &[haploid, diploid], 5),
        Err(CallError::InvalidPriorChromosomeCount)
    );
    let called = MultiallelicCaller::default()
        .call_with_ploidies(&site, &[haploid, diploid], 4)
        .unwrap();

    assert_eq!(called.allele_counts(), &[2, 1]);
    assert_eq!(called.allele_number(), 3);
    let quality = called.quality().unwrap();
    assert!((quality - 7.817_96).abs() < 1e-5, "{quality}");
    assert_eq!(called.samples()[0].genotype(), Some(&[0][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(40));
    assert_eq!(called.samples()[1].genotype_quality(), Some(3));
}

#[test]
fn matches_bcftools_1_24_absent_ploidy_call() {
    let site = two_sample_site();

    let called = MultiallelicCaller::default()
        .call_with_ploidies(&site, &[CallPloidy::Absent, CallPloidy::Diploid], 4)
        .unwrap();

    assert_eq!(called.allele_counts(), &[1, 1]);
    assert_eq!(called.allele_number(), 2);
    assert!((called.quality().unwrap() - 13.2678).abs() < 1e-4);
    assert_eq!(called.samples()[0].ploidy(), CallPloidy::Absent);
    assert_eq!(called.samples()[0].genotype(), None);
    assert_eq!(called.samples()[0].phred_likelihoods(), None);
    assert_eq!(called.samples()[0].evidence().allele_depths(), &[1, 0]);
    assert_eq!(called.samples()[1].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[1].genotype_quality(), Some(3));
    assert_eq!(
        called.samples()[1].phred_likelihoods(),
        Some(&[40, 3, 0][..])
    );
}

#[test]
fn matches_bcftools_1_24_alt_only_selection() {
    let site = LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [
            Allele::new(&b"G"[..]).unwrap(),
            Allele::new(&b"<*>"[..]).unwrap(),
        ],
        [0.0, 1.0, 0.0],
        [sample([40, 3, 0, 40, 3, 40], [0, 1, 0], [0, 40, 0])],
    )
    .unwrap();

    let called = MultiallelicCaller::default().call(&site).unwrap();

    assert_eq!(called.allele_counts(), &[0, 2]);
    assert!((called.quality().unwrap() - 10.7923).abs() < 1e-4);
    assert_eq!(called.samples()[0].genotype(), Some(&[1, 1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(127));
    assert_eq!(
        called.samples()[0].genotype_probabilities(),
        Some(&[0.0, 0.0, 1.0][..])
    );
}

#[test]
fn matches_bcftools_1_24_per_sample_groups() {
    let site = two_sample_site();
    let ploidies = [CallPloidy::Diploid, CallPloidy::Diploid];

    assert_eq!(
        MultiallelicCaller::default().call_with_groups(&site, &ploidies, 4, &[0, 2],),
        Err(CallError::InvalidCallerGroups)
    );
    let called = MultiallelicCaller::default()
        .call_with_groups(&site, &ploidies, 4, &[0, 1])
        .unwrap();

    assert_eq!(called.allele_counts(), &[2, 2]);
    assert!((called.quality().unwrap() - 13.2571).abs() < 1e-4);
    assert_eq!(called.samples()[0].genotype(), Some(&[0, 0][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[1, 1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(127));
    assert_eq!(called.samples()[1].genotype_quality(), Some(127));
    assert_eq!(
        called.samples()[0].genotype_probabilities(),
        Some(&[1.0, 0.0, 0.0][..])
    );
    assert_eq!(
        called.samples()[1].genotype_probabilities(),
        Some(&[0.0, 0.0, 1.0][..])
    );
}

#[test]
fn matches_bcftools_1_24_reference_call() {
    let site = LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"<*>"[..]).unwrap()],
        [1.0, 0.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            [0, 3, 40],
            SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap();

    let called = MultiallelicCaller::default().call(&site).unwrap();

    assert!(called.alternates().is_empty());
    assert_eq!(called.allele_counts(), &[2]);
    assert!(!called.is_variant());
    assert!((called.quality().unwrap() - 69.587).abs() < 1e-3);
    assert_eq!(called.samples()[0].genotype(), Some(&[0, 0][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(0));
    assert_eq!(called.samples()[0].genotype_probabilities(), None);
    assert_eq!(called.samples()[0].phred_likelihoods(), None);
    assert_eq!(called.samples()[0].evidence().allele_quality_sums(), &[40]);
}

#[test]
fn matches_bcftools_1_24_triallelic_call() {
    let make_sample = |likelihoods, depths, sums| {
        SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            likelihoods,
            SampleEvidence::new(1, depths, sums).unwrap(),
        )
        .unwrap()
    };
    let site = LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [
            Allele::new(&b"G"[..]).unwrap(),
            Allele::new(&b"C"[..]).unwrap(),
            Allele::new(&b"<*>"[..]).unwrap(),
        ],
        [1.0, 1.0, 1.0, 0.0],
        [
            make_sample(
                [0, 3, 40, 3, 40, 40, 3, 40, 40, 40],
                [1, 0, 0, 0],
                [40, 0, 0, 0],
            ),
            make_sample(
                [40, 40, 40, 3, 3, 0, 40, 40, 3, 40],
                [0, 0, 1, 0],
                [0, 0, 40, 0],
            ),
            make_sample(
                [40, 3, 0, 40, 3, 40, 40, 3, 40, 40],
                [0, 1, 0, 0],
                [0, 40, 0, 0],
            ),
        ],
    )
    .unwrap();

    let called = MultiallelicCaller::default().call(&site).unwrap();

    assert_eq!(
        called
            .alternates()
            .iter()
            .map(Allele::as_bytes)
            .collect::<Vec<_>>(),
        [b"G".as_slice(), b"C".as_slice()]
    );
    assert_eq!(called.allele_counts(), &[3, 2, 1]);
    assert!((called.quality().unwrap() - 15.6934).abs() < 1e-4);
    assert_eq!(called.samples()[0].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[0, 2][..]));
    assert_eq!(called.samples()[2].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(1));
    assert_eq!(called.samples()[1].genotype_quality(), Some(1));
    assert_eq!(called.samples()[2].genotype_quality(), Some(1));
    assert_eq!(
        called.samples()[1].phred_likelihoods(),
        Some(&[40, 40, 40, 3, 3, 0][..])
    );
    let expected = [
        0.332_762,
        0.333_552,
        3.327_62e-5,
        0.333_552,
        6.655_24e-5,
        3.327_62e-5,
    ];
    for (&observed, expected) in called.samples()[0]
        .genotype_probabilities()
        .unwrap()
        .iter()
        .zip(expected)
    {
        assert!((observed - expected).abs() < 5e-7);
    }
}
