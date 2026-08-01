use super::*;
use crate::PloidyPreset;

const INPUT: &[u8] = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=10>\n\
##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfirst\tsecond\tthird\n\
chr1\t1\t.\tA\t<*>\t.\t.\tQS=1,0\tPL\t0,100,200\t0,100,200\t0,100,200\n";

#[test]
fn binds_file_order_and_fixed_ploidy() {
    let selection =
        CallSampleSelection::parse("# selected cohort\nthird 0\nfirst 1\n\nsecond 2\n").unwrap();
    let definition = PloidyDefinition::parse("* * * F 2\n").unwrap();
    let (reader, resolver) = selection
        .bind(LikelihoodVariantReader::new(INPUT).unwrap(), &definition)
        .unwrap();
    assert_eq!(
        reader
            .schema()
            .header()
            .sample_names()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        ["third", "first", "second"]
    );
    assert_eq!(
        resolver.resolve("chr1", 0).unwrap(),
        [CallPloidy::Absent, CallPloidy::Haploid, CallPloidy::Diploid]
    );
}

#[test]
fn exclusion_keeps_header_order_and_uses_the_definition_default() {
    let selection = CallSampleSelection::parse_excluding("second ignored\n").unwrap();
    let definition = PloidyDefinition::preset(PloidyPreset::Haploid);
    let (reader, resolver) = selection
        .bind(LikelihoodVariantReader::new(INPUT).unwrap(), &definition)
        .unwrap();
    assert_eq!(
        reader
            .schema()
            .header()
            .sample_names()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        ["first", "third"]
    );
    assert_eq!(
        resolver.resolve("chr1", 0).unwrap(),
        [CallPloidy::Haploid, CallPloidy::Haploid]
    );
}

#[test]
fn omitted_assignment_uses_the_definition_default_sex() {
    let selection = CallSampleSelection::parse("first\n").unwrap();
    let definition = PloidyDefinition::parse("Y 1 10 F 0\n* * * F 2\n").unwrap();
    let (_, resolver) = selection
        .bind(LikelihoodVariantReader::new(INPUT).unwrap(), &definition)
        .unwrap();
    assert_eq!(resolver.resolve("Y", 0).unwrap(), [CallPloidy::Absent]);
}

#[test]
fn malformed_duplicate_unknown_and_undeclared_samples_fail() {
    assert!(matches!(
        CallSampleSelection::parse("first F extra\n"),
        Err(CallError::CallSampleRecord { line: 1, .. })
    ));
    assert_eq!(
        CallSampleSelection::parse("first F\nfirst M\n").unwrap_err(),
        CallError::DuplicateSampleSelection("first".to_owned())
    );
    let error = match CallSampleSelection::parse("missing F\n").unwrap().bind(
        LikelihoodVariantReader::new(INPUT).unwrap(),
        &PloidyDefinition::preset(PloidyPreset::Diploid),
    ) {
        Ok(_) => panic!("missing sample was accepted"),
        Err(error) => error,
    };
    assert_eq!(
        error,
        CallError::MissingSelectedSample("missing".to_owned())
    );
    let error = match CallSampleSelection::parse("first unknown\n").unwrap().bind(
        LikelihoodVariantReader::new(INPUT).unwrap(),
        &PloidyDefinition::parse("* * * F 2\n").unwrap(),
    ) {
        Ok(_) => panic!("undeclared sex was accepted"),
        Err(error) => error,
    };
    assert_eq!(error, CallError::UnknownPloidySex("unknown".to_owned()));
    assert!(matches!(
        CallSampleSelection::parse("\n# none\n"),
        Err(CallError::InvalidSampleCount)
    ));
}
