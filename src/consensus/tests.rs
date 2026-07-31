use super::*;
use crate::{Allele, CalledVcfSchema, LikelihoodVariantReader, Ploidy, SampleEvidence};

#[test]
fn matches_bcftools_1_24_triallelic_consensus_call() {
    let data = include_bytes!("../../tests/golden/bcftools-1.24-likelihood.vcf");
    let mut reader = LikelihoodVariantReader::new(&data[..]).unwrap();
    let input_schema = reader.schema().clone();
    let site = reader.read_site().unwrap().unwrap();

    let called = ConsensusCaller::default().call(&site).unwrap();

    assert_eq!(
        called
            .alternates()
            .iter()
            .map(Allele::as_bytes)
            .collect::<Vec<_>>(),
        [b"G".as_slice()]
    );
    assert!((called.quality().unwrap() - 7.799_77).abs() < 1e-5);
    assert_eq!(called.allele_counts(), &[5, 1]);
    assert_eq!(called.samples()[0].genotype(), Some(&[0, 0][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[0, 0][..]));
    assert_eq!(called.samples()[2].genotype(), Some(&[0, 1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), Some(5));
    assert_eq!(called.samples()[1].genotype_quality(), Some(3));
    assert_eq!(called.samples()[2].genotype_quality(), Some(5));
    assert_eq!(
        called.samples()[1].phred_likelihoods(),
        Some(&[40, 40, 40][..])
    );
    assert_eq!(called.samples()[1].evidence().allele_depths(), &[0, 0]);
    let output_schema = CalledVcfSchema::from_consensus_likelihood(&input_schema);
    assert!(!output_schema.header().formats().contains_key("GP"));
    let record = output_schema.encode_call(&called).unwrap();
    assert_eq!(
        record
            .samples()
            .keys()
            .as_ref()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        ["GT", "PL", "DP", "AD", "QS", "GQ"]
    );
}

#[test]
fn matches_bcftools_1_24_haploid_consensus_call() {
    let data = include_bytes!("../../tests/golden/bcftools-1.24-likelihood.vcf");
    let mut reader = LikelihoodVariantReader::new(&data[..]).unwrap();
    let site = reader.read_site().unwrap().unwrap();
    let ploidies = vec![CallPloidy::Haploid; site.samples().len()];

    let called = ConsensusCaller::default()
        .call_with_ploidies(&site, &ploidies)
        .unwrap();

    assert!((called.quality().unwrap() - 7.788_86).abs() < 1e-5);
    assert_eq!(called.allele_counts(), &[2, 1]);
    assert_eq!(called.samples()[0].genotype(), Some(&[0][..]));
    assert_eq!(called.samples()[1].genotype(), Some(&[0][..]));
    assert_eq!(called.samples()[2].genotype(), Some(&[1][..]));
    assert_eq!(called.samples()[0].genotype_quality(), None);
    assert_eq!(called.samples()[1].genotype_quality(), None);
    assert_eq!(called.samples()[2].genotype_quality(), None);
    assert_eq!(called.samples()[0].phred_likelihoods(), Some(&[0, 40][..]));
}

#[test]
fn matches_bcftools_1_24_reference_heterozygous_and_alternate_calls() {
    let cases = [
        ([0, 3, 40], 32.9956, [0, 0], 36),
        ([182, 0, 172], 152.008, [0, 1], 99),
        ([220, 99, 0], 186.999, [1, 1], 99),
    ];
    for (phred_likelihoods, quality, genotype, genotype_quality) in cases {
        let site = biallelic_site(phred_likelihoods);
        let called = ConsensusCaller::default().call(&site).unwrap();

        assert!(
            (called.quality().unwrap() - quality).abs() < 1e-3,
            "{} != {quality}",
            called.quality().unwrap()
        );
        assert_eq!(called.samples()[0].genotype(), Some(&genotype[..]));
        assert_eq!(
            called.samples()[0].genotype_quality(),
            Some(genotype_quality)
        );
    }
}

#[test]
fn validates_threshold_ploidy_and_likelihood_inputs() {
    assert_eq!(
        ConsensusCallerConfig::new(f64::NAN),
        Err(CallError::InvalidConsensusThreshold)
    );
    let site = biallelic_site([0, 3, 40]);
    assert_eq!(
        ConsensusCaller::default().call_with_ploidies(&site, &[]),
        Err(CallError::CallerPloidyCountMismatch)
    );
    let missing = LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"G"[..]).unwrap()],
        [1.0, 1.0],
        [SampleLikelihood::missing(2, Ploidy::new(2).unwrap()).unwrap()],
    )
    .unwrap();
    assert_eq!(
        ConsensusCaller::default().call(&missing),
        Err(CallError::UnsupportedConsensusLikelihoods)
    );
}

fn biallelic_site(phred_likelihoods: [u32; 3]) -> LikelihoodSite {
    LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"G"[..]).unwrap()],
        [1.0, 1.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            phred_likelihoods,
            SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap()
}
