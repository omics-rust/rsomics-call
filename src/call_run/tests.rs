use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::{
    IndexedLikelihoodVariantReader, PloidyDefinition, PloidyPreset, SamplePloidy,
    VariantOutputFormat,
};

use super::*;

const INPUT: &[u8] = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=10>\n\
##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n\
chr1\t1\t.\tA\t<*>\t.\t.\tQS=1,0\tPL:DP\t0,100,200:10\n\
chr1\t2\t.\tA\t<*>\t.\t.\tQS=1,0\tPL:DP\t0,100,200:8\n";

const INPUT_WITHOUT_DP: &[u8] = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=10>\n\
##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n\
chr1\t1\t.\tA\t<*>\t.\t.\tQS=1,0\tPL\t0,100,200\n";

#[test]
fn composes_ploidy_calling_gvcf_and_output() {
    let reader = LikelihoodVariantReader::new(INPUT).unwrap();
    let ploidy = PloidyDefinition::preset(PloidyPreset::Diploid)
        .default_resolver(1)
        .unwrap();
    let output = LikelihoodCallRun::new(
        CallModel::Multiallelic(MultiallelicCallerConfig::default()),
        ploidy,
    )
    .with_gvcf([5])
    .unwrap()
    .run(reader, Vec::new(), VariantOutputFormat::Vcf)
    .unwrap();
    let record = String::from_utf8(output)
        .unwrap()
        .lines()
        .find(|line| !line.starts_with('#'))
        .unwrap()
        .to_owned();
    assert_eq!(
        record,
        "chr1\t1\t.\tA\t.\t.\t.\tEND=2;MIN_DP=8\tGT:DP\t0/0:8"
    );
}

#[test]
fn rejects_incompatible_workflow_configuration_before_output() {
    let ploidy = PloidyDefinition::preset(PloidyPreset::Diploid)
        .default_resolver(1)
        .unwrap();
    assert!(matches!(
        LikelihoodCallRun::new(
            CallModel::Consensus(ConsensusCallerConfig::default()),
            ploidy.clone(),
        )
        .with_gvcf([5]),
        Err(CallError::UnsupportedGvcfModel)
    ));

    let reader = LikelihoodVariantReader::new(INPUT).unwrap();
    let mismatched = PloidyDefinition::preset(PloidyPreset::Diploid)
        .resolver(vec![SamplePloidy::Fixed(crate::CallPloidy::Diploid); 2])
        .unwrap();
    assert!(matches!(
        LikelihoodCallRun::new(
            CallModel::Multiallelic(MultiallelicCallerConfig::default()),
            mismatched,
        )
        .run(reader, Vec::new(), VariantOutputFormat::Vcf),
        Err(CallError::PloidySampleCountMismatch)
    ));

    let reader = LikelihoodVariantReader::new(INPUT_WITHOUT_DP).unwrap();
    assert!(matches!(
        LikelihoodCallRun::new(
            CallModel::Multiallelic(MultiallelicCallerConfig::default()),
            ploidy,
        )
        .with_gvcf([5])
        .unwrap()
        .run(reader, Vec::new(), VariantOutputFormat::Vcf),
        Err(CallError::MissingGvcfDepth)
    ));
}

#[test]
fn runs_the_consensus_model_with_its_output_schema() {
    let reader = LikelihoodVariantReader::new(INPUT).unwrap();
    let ploidy = PloidyDefinition::preset(PloidyPreset::Diploid)
        .default_resolver(1)
        .unwrap();
    let output = LikelihoodCallRun::new(
        CallModel::Consensus(ConsensusCallerConfig::default()),
        ploidy,
    )
    .run(reader, Vec::new(), VariantOutputFormat::Vcf)
    .unwrap();
    let output = String::from_utf8(output).unwrap();
    assert!(!output.contains("##FORMAT=<ID=GP"));
    assert_eq!(
        output.lines().filter(|line| !line.starts_with('#')).count(),
        2
    );
}

#[test]
fn runs_calls_from_an_indexed_region() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_indexed_vcf(directory.path(), INPUT);

    let reader = IndexedLikelihoodVariantReader::open(input).unwrap();
    let ploidy = PloidyDefinition::preset(PloidyPreset::Diploid)
        .default_resolver(1)
        .unwrap();
    let output = LikelihoodCallRun::new(
        CallModel::Multiallelic(MultiallelicCallerConfig::default()),
        ploidy,
    )
    .run_indexed(
        reader,
        ["chr1:2-2".parse().unwrap()],
        Vec::new(),
        VariantOutputFormat::Vcf,
    )
    .unwrap();
    let records = String::from_utf8(output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(records.len(), 1);
    assert!(records[0].starts_with("chr1\t2\t"));
}

#[test]
fn regions_file_controls_reference_order_and_deduplicates_overlaps() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_indexed_vcf(
        directory.path(),
        b"##fileformat=VCFv4.2\n\
          ##contig=<ID=chr1,length=10>\n\
          ##contig=<ID=chr2,length=10>\n\
          ##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n\
          ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n\
          chr1\t1\t.\tA\tG\t.\t.\tQS=1,1\tPL\t40,3,0\n\
          chr1\t2\t.\tA\tG\t.\t.\tQS=1,1\tPL\t40,3,0\n\
          chr2\t2\t.\tA\tG\t.\t.\tQS=1,1\tPL\t40,3,0\n",
    );
    let regions = directory.path().join("regions.txt");
    fs::write(&regions, b"chr2\t2\t2\nchr1\t1\t2\nchr1\t2\t2\n").unwrap();
    let reader = IndexedLikelihoodVariantReader::open(input).unwrap();
    let ploidy = PloidyDefinition::preset(PloidyPreset::Diploid)
        .default_resolver(1)
        .unwrap();
    let output = LikelihoodCallRun::new(
        CallModel::Multiallelic(MultiallelicCallerConfig::default()),
        ploidy,
    )
    .run_indexed_regions_file(reader, regions, Vec::new(), VariantOutputFormat::Vcf)
    .unwrap();
    assert_eq!(
        String::from_utf8(output)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').take(2).collect::<Vec<_>>().join(":"))
            .collect::<Vec<_>>(),
        ["chr2:2", "chr1:1", "chr1:2"]
    );
}

fn write_indexed_vcf(directory: &Path, data: &[u8]) -> PathBuf {
    let input = directory.join("likelihoods.vcf.gz");
    let mut writer = noodles::bgzf::io::Writer::new(File::create(&input).unwrap());
    writer.write_all(data).unwrap();
    writer.finish().unwrap();
    let index = noodles::vcf::fs::index(&input).unwrap();
    noodles::tabix::fs::write(format!("{}.tbi", input.display()), &index).unwrap();
    input
}
