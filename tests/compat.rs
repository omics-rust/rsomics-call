use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::vcf;
use rsomics_call::{
    AlignmentInput, IndelLikelihoodConfig, LikelihoodSite, LikelihoodVcfSchema, SampleSelection,
    SnpLikelihoodConfig, SnpLikelihoodRun,
};
use rsomics_pileup::PileupOptions;

fn bcftools() -> String {
    std::env::var("BCFTOOLS").unwrap_or_else(|_| "bcftools".to_owned())
}

fn assert_bcftools_1_24() {
    let output = Command::new(bcftools()).arg("--version").output().unwrap();
    assert!(output.status.success());
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .next()
            .unwrap()
            .starts_with("bcftools 1.24")
    );
}

fn write_reference(directory: &Path, sequence: &str) -> PathBuf {
    let reference = directory.join("reference.fa");
    fs::write(&reference, format!(">MX\n{sequence}\n")).unwrap();
    fs::write(
        reference.with_extension("fa.fai"),
        format!(
            "MX\t{}\t4\t{}\t{}\n",
            sequence.len(),
            sequence.len(),
            sequence.len() + 1
        ),
    )
    .unwrap();
    reference
}

fn write_alignment(
    directory: &Path,
    name: &str,
    sample: &str,
    reference_length: usize,
    cigar: &str,
    sequence: &str,
) -> PathBuf {
    let path = directory.join(format!("{name}.sam"));
    let qualities = "I".repeat(sequence.len());
    let mut data = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:MX\tLN:{reference_length}\n@RG\tID:rg\tSM:{sample}\n"
    );
    for index in 1..=2 {
        data.push_str(&format!(
            "{name}{index}\t0\tMX\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{qualities}\tRG:Z:rg\n"
        ));
    }
    fs::write(&path, data).unwrap();
    path
}

fn bcftools_indel(reference: &Path, alignments: &[PathBuf]) -> LikelihoodSite {
    let mut command = Command::new(bcftools());
    command.args(["mpileup", "-f"]).arg(reference).args([
        "-a",
        "FORMAT/DP,FORMAT/AD,FORMAT/QS",
        "-Ov",
    ]);
    command.args(alignments);
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut reader = vcf::io::Reader::new(&output.stdout[..]);
    let header = reader.read_header().unwrap();
    let schema = LikelihoodVcfSchema::from_header(header.clone()).unwrap();
    let mut record = vcf::variant::RecordBuf::default();
    loop {
        assert_ne!(reader.read_record_buf(&header, &mut record).unwrap(), 0);
        let site = schema.decode_likelihood(&record).unwrap();
        if site.indel_summary().is_some() {
            return site;
        }
    }
}

fn rsomics_indel(reference: &Path, alignments: &[PathBuf]) -> LikelihoodSite {
    let inputs = alignments
        .iter()
        .enumerate()
        .map(|(index, path)| AlignmentInput::new(index as u32, path, path.display().to_string()));
    let run = SnpLikelihoodRun::open(
        inputs,
        reference,
        SampleSelection::default(),
        PileupOptions::default(),
        SnpLikelihoodConfig::default(),
    )
    .unwrap()
    .with_partial_baq(500, false)
    .with_indels(IndelLikelihoodConfig::default())
    .unwrap();
    let mut indel = None;
    run.run(|site| {
        if site.indel_summary().is_some() {
            assert!(indel.replace(site).is_none());
        }
        Ok(())
    })
    .unwrap();
    indel.unwrap()
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn indel_likelihoods_match_bcftools_1_24() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let reference = write_reference(directory.path(), "CGTCTACTACG");
    let insertion = write_alignment(
        directory.path(),
        "insertion",
        "alternate",
        11,
        "5M1I6M",
        "CGTCTCACTACG",
    );
    let deletion = write_alignment(
        directory.path(),
        "deletion",
        "alternate",
        11,
        "5M1D5M",
        "CGTCTCTACG",
    );
    let reference_reads = write_alignment(
        directory.path(),
        "reference",
        "reference",
        11,
        "11M",
        "CGTCTACTACG",
    );
    let trailing = write_alignment(directory.path(), "trailing", "trailing", 11, "2M", "CG");

    for alignments in [
        vec![insertion.clone()],
        vec![deletion],
        vec![insertion, reference_reads, trailing],
    ] {
        assert_eq!(
            rsomics_indel(&reference, &alignments),
            bcftools_indel(&reference, &alignments)
        );
    }

    let repeat_directory = directory.path().join("repeat");
    fs::create_dir(&repeat_directory).unwrap();
    let repeat_reference = write_reference(&repeat_directory, "CGTCAAAAAACTG");
    let repeat_insertion = write_alignment(
        &repeat_directory,
        "repeat-insertion",
        "alternate",
        13,
        "7M1I6M",
        "CGTCAAAAAAACTG",
    );
    assert_eq!(
        rsomics_indel(&repeat_reference, std::slice::from_ref(&repeat_insertion)),
        bcftools_indel(&repeat_reference, &[repeat_insertion])
    );
}
