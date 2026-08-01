use super::*;

#[test]
fn grch37_resolves_sex_chromosomes_and_prior() {
    let definition = PloidyDefinition::preset(PloidyPreset::Grch37);
    let resolver = definition
        .resolver([SamplePloidy::sex("M"), SamplePloidy::sex("F")])
        .unwrap();
    assert_eq!(resolver.prior_chromosome_count(), 4);
    assert_eq!(
        resolver.resolve("X", 59_999).unwrap(),
        [CallPloidy::Haploid, CallPloidy::Diploid]
    );
    assert_eq!(
        resolver.resolve("X", 60_000).unwrap(),
        [CallPloidy::Diploid, CallPloidy::Diploid]
    );
    assert_eq!(
        resolver.resolve("Y", 0).unwrap(),
        [CallPloidy::Haploid, CallPloidy::Absent]
    );
    assert_eq!(
        resolver.resolve("chrM", 0).unwrap(),
        [CallPloidy::Haploid, CallPloidy::Haploid]
    );
    assert_eq!(
        definition
            .default_resolver(1)
            .unwrap()
            .resolve("Y", 0)
            .unwrap(),
        [CallPloidy::Absent]
    );
}

#[test]
fn grch38_uses_assembly_specific_par_boundaries() {
    let resolver = PloidyDefinition::preset(PloidyPreset::Grch38)
        .resolver([SamplePloidy::sex("M")])
        .unwrap();
    assert_eq!(
        resolver.resolve("chrX", 9_998).unwrap(),
        [CallPloidy::Haploid]
    );
    assert_eq!(
        resolver.resolve("chrX", 9_999).unwrap(),
        [CallPloidy::Diploid]
    );
    assert_eq!(
        resolver.resolve("chrX", 2_781_479).unwrap(),
        [CallPloidy::Haploid]
    );
}

#[test]
fn custom_definition_uses_defaults_and_fixed_assignments() {
    let definition =
        PloidyDefinition::parse("X 1 10 alpha 1\nY 1 10 alpha 0\n* * * alpha 2\n* * * beta 1\n")
            .unwrap();
    let resolver = definition
        .resolver([
            SamplePloidy::sex("alpha"),
            SamplePloidy::sex("beta"),
            SamplePloidy::Fixed(CallPloidy::Haploid),
        ])
        .unwrap();
    assert_eq!(resolver.prior_chromosome_count(), 4);
    assert_eq!(
        resolver.resolve("X", 0).unwrap(),
        [
            CallPloidy::Haploid,
            CallPloidy::Haploid,
            CallPloidy::Haploid
        ]
    );
    assert_eq!(
        resolver.resolve("Y", 0).unwrap(),
        [CallPloidy::Absent, CallPloidy::Haploid, CallPloidy::Haploid]
    );
}

#[test]
fn checked_custom_definition_rejects_ambiguous_or_malformed_records() {
    for data in [
        "X 1 10 M 1\nX 10 20 M 2\n",
        "* 1 * M 2\n",
        "X 0 10 M 1\n",
        "X 1 10 M 3\n",
        "X 1 10 M\n",
        "# empty\n",
    ] {
        assert!(PloidyDefinition::parse(data).is_err(), "{data}");
    }
    let definition = PloidyDefinition::parse("* * * M 2\n* * * F 2\n").unwrap();
    assert_eq!(
        definition
            .resolver([SamplePloidy::sex("unknown")])
            .unwrap_err(),
        CallError::UnknownPloidySex("unknown".to_owned())
    );
}

#[test]
fn constant_presets_set_the_prior_to_their_real_maximum() {
    let haploid = PloidyDefinition::preset(PloidyPreset::Haploid)
        .default_resolver(2)
        .unwrap();
    assert_eq!(haploid.prior_chromosome_count(), 2);
    assert_eq!(
        haploid.resolve("any", u64::MAX).unwrap(),
        [CallPloidy::Haploid, CallPloidy::Haploid]
    );
    let diploid = PloidyDefinition::preset(PloidyPreset::Diploid)
        .default_resolver(2)
        .unwrap();
    assert_eq!(diploid.prior_chromosome_count(), 4);
}
