use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, sam, vcf};
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

fn write_indexed_alignment(
    directory: &Path,
    name: &str,
    sample: &str,
    reference_length: usize,
    records: &str,
) -> PathBuf {
    let path = directory.join(format!("{name}.bam"));
    let source = format!(
        "@HD\tVN:1.6\tSO:coordinate\n\
         @SQ\tSN:MX\tLN:{reference_length}\n\
         @RG\tID:rg\tSM:{sample}\n\
         {records}"
    );
    let mut reader = sam::io::Reader::new(source.as_bytes());
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap();
    let mut writer = bam::io::Writer::new(File::create(&path).unwrap());
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    let index = bam::fs::index(&path).unwrap();
    bam::bai::fs::write(path.with_extension("bai"), &index).unwrap();
    path
}

fn decode_likelihoods(output: &[u8]) -> Vec<LikelihoodSite> {
    let mut reader = vcf::io::Reader::new(output);
    let header = reader.read_header().unwrap();
    let schema = LikelihoodVcfSchema::from_header(header.clone()).unwrap();
    let mut record = vcf::variant::RecordBuf::default();
    let mut sites = Vec::new();
    while reader.read_record_buf(&header, &mut record).unwrap() != 0 {
        sites.push(schema.decode_likelihood(&record).unwrap());
    }
    sites
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

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn indexed_regions_and_streaming_targets_match_bcftools_1_24() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let reference = write_reference(directory.path(), "ACGTACGTACG");
    let first = write_indexed_alignment(
        directory.path(),
        "first",
        "S1",
        11,
        "first\t0\tMX\t2\t60\t5M\t*\t0\t0\tCGTAC\tIIIII\tRG:Z:rg\n",
    );
    let second = write_indexed_alignment(
        directory.path(),
        "second",
        "S2",
        11,
        "second\t0\tMX\t6\t60\t3M\t*\t0\t0\tCAT\tIII\tRG:Z:rg\n",
    );
    let regions = "MX:3-5,MX:7-8";
    let targets = "MX:4-4,MX:7-8";

    let output = Command::new(bcftools())
        .args(["mpileup", "-B", "-f"])
        .arg(&reference)
        .args([
            "-a",
            "FORMAT/DP,FORMAT/AD,FORMAT/QS",
            "-r",
            regions,
            "-t",
            targets,
            "-Ov",
        ])
        .args([&first, &second])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let expected = decode_likelihoods(&output.stdout);
    let run = SnpLikelihoodRun::open_regions(
        [
            AlignmentInput::new(1, first, "first"),
            AlignmentInput::new(2, second, "second"),
        ],
        reference,
        SampleSelection::default(),
        ["MX:7-8", "MX:4-5", "MX:3-4", "MX:3-4"].map(|region| region.parse().unwrap()),
        PileupOptions::default(),
        SnpLikelihoodConfig::default(),
    )
    .unwrap()
    .with_targets(["MX:8-8", "MX:4-4", "MX:7-7"].map(|region| region.parse().unwrap()));
    let mut actual = Vec::new();
    run.run(|site| {
        actual.push(site);
        Ok(())
    })
    .unwrap();

    assert_eq!(actual, expected);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn unindexed_streaming_targets_match_bcftools_1_24() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let reference = write_reference(directory.path(), "ACGTACGTACG");
    let input = write_alignment(directory.path(), "input", "S1", 11, "11M", "ACGTACGTACG");
    let bed = directory.path().join("targets.bed");
    fs::write(&bed, b"MX\t1\t3\nMX\t6\t8\n").unwrap();
    let tab = directory.path().join("targets.txt");
    fs::write(&tab, b"MX\t2\t3\nMX\t7\t8\n").unwrap();
    let vcf = directory.path().join("targets.vcf");
    fs::write(
        &vcf,
        b"##fileformat=VCFv4.3\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
          MX\t2\t.\tC\t.\t.\t.\t.\n\
          MX\t3\t.\tG\t.\t.\t.\t.\n\
          MX\t7\t.\tG\t.\t.\t.\t.\n\
          MX\t8\t.\tT\t.\t.\t.\t.\n",
    )
    .unwrap();

    for targets in [bed, tab, vcf] {
        let output = Command::new(bcftools())
            .args(["mpileup", "-B", "-f"])
            .arg(&reference)
            .args(["-a", "FORMAT/DP,FORMAT/AD,FORMAT/QS", "-T"])
            .arg(&targets)
            .arg("-Ov")
            .arg(&input)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let expected = decode_likelihoods(&output.stdout);
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, &input, "input")],
            &reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap()
        .with_target_file(&targets)
        .unwrap();
        let mut actual = Vec::new();
        run.run(|site| {
            actual.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(actual, expected, "{}", targets.display());
    }
}
