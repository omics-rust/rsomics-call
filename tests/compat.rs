use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::sam::alignment::io::Write as _;
use noodles::vcf::variant::io::Write as _;
use noodles::{bam, sam, vcf};
use rsomics_call::{
    AlignmentInput, CalledVcfSchema, IndelLikelihoodConfig, LikelihoodSite,
    LikelihoodVariantReader, LikelihoodVcfSchema, MultiallelicCaller, PloidyDefinition,
    PloidyPreset, SamplePloidy, SampleSelection, SnpLikelihoodConfig, SnpLikelihoodRun,
};
use rsomics_pileup::PileupOptions;

const ANNOTATIONS: &str = "FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/QM,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/FS,INFO/NMBZ,INFO/NM,INFO/SCR";

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
    let mut reader = LikelihoodVariantReader::new(output).unwrap();
    let mut sites = Vec::new();
    while let Some(site) = reader.read_site().unwrap() {
        sites.push(site);
    }
    sites
}

fn assert_sites_match(actual: &[LikelihoodSite], expected: &[LikelihoodSite]) {
    assert_eq!(actual.len(), expected.len());
    if actual.is_empty() {
        return;
    }
    let reference_count = actual
        .iter()
        .chain(expected)
        .map(LikelihoodSite::reference_sequence_id)
        .max()
        .unwrap()
        + 1;
    let references = (0..reference_count)
        .map(|index| (format!("ref{index}").into_bytes(), 1_000_000))
        .collect::<Vec<_>>();
    let samples = (0..actual[0].samples().len())
        .map(|index| format!("sample{index}"))
        .collect::<Vec<_>>();
    let schema = LikelihoodVcfSchema::new(references, &samples).unwrap();
    for (actual, expected) in actual.iter().zip(expected) {
        let encode = |site: &LikelihoodSite| {
            let record = schema.encode_likelihood(site).unwrap();
            let mut data = Vec::new();
            let mut writer = vcf::io::Writer::new(&mut data);
            writer
                .write_variant_record(schema.header(), &record)
                .unwrap();
            String::from_utf8(data).unwrap()
        };
        assert_eq!(
            encode(actual),
            encode(expected),
            "site {}:{}",
            actual.reference_sequence_id(),
            actual.position()
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn multiallelic_call_annotations_match_bcftools_1_24() {
    assert_bcftools_1_24();
    let input = Path::new("tests/golden/bcftools-1.24-likelihood.vcf");
    let data = fs::read(input).unwrap();
    let mut likelihood_reader = LikelihoodVariantReader::new(&data[..]).unwrap();
    let site = likelihood_reader.read_site().unwrap().unwrap();
    assert!(likelihood_reader.read_site().unwrap().is_none());
    let schema = CalledVcfSchema::from_likelihood(likelihood_reader.schema());
    let actual = schema
        .encode_call(&MultiallelicCaller::default().call(&site).unwrap())
        .unwrap();

    let output = Command::new(bcftools())
        .args(["call", "-m", "-a", "PV4,GP,GQ", "-Ov"])
        .arg(input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut reader = vcf::io::Reader::new(&output.stdout[..]);
    let header = reader.read_header().unwrap();
    let mut expected = vcf::variant::RecordBuf::default();
    assert_ne!(reader.read_record_buf(&header, &mut expected).unwrap(), 0);
    assert_eq!(
        actual.reference_sequence_name(),
        expected.reference_sequence_name()
    );
    assert_eq!(actual.variant_start(), expected.variant_start());
    assert_eq!(actual.reference_bases(), expected.reference_bases());
    assert_eq!(actual.alternate_bases(), expected.alternate_bases());
    assert!((actual.quality_score().unwrap() - expected.quality_score().unwrap()).abs() < 1e-4);
    assert_eq!(actual.info().as_ref().len(), expected.info().as_ref().len());
    for (key, value) in expected.info().as_ref() {
        assert_eq!(actual.info().as_ref().get(key), Some(value), "INFO/{key}");
    }
    assert_eq!(actual.samples().keys(), expected.samples().keys());
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn grch37_ploidy_preset_matches_bcftools_1_24() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("ploidy.vcf");
    let samples = directory.path().join("samples.txt");
    fs::write(&samples, "male M\nfemale F\n").unwrap();
    let sites = [
        ("X", 1u64),
        ("X", 60_001),
        ("X", 2_699_521),
        ("Y", 1),
        ("MT", 1),
        ("chrX", 1),
        ("chrY", 1),
        ("chrM", 1),
    ];
    let mut data = "##fileformat=VCFv4.2\n\
                    ##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n\
                    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods\">\n"
        .to_owned();
    for reference in ["X", "Y", "MT", "chrX", "chrY", "chrM"] {
        data.push_str(&format!("##contig=<ID={reference},length=200000000>\n"));
    }
    data.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale\tfemale\n");
    for &(reference, position) in &sites {
        data.push_str(&format!(
            "{reference}\t{position}\t.\tA\tG,<*>\t.\t.\tQS=1,0,0\tPL\t0,100,200,100,200,200\t0,100,200,100,200,200\n"
        ));
    }
    fs::write(&input, data).unwrap();

    let output = Command::new(bcftools())
        .args(["call", "-m", "--ploidy", "GRCh37", "-S"])
        .arg(&samples)
        .arg("-Ov")
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let observed = called_genotypes(&output.stdout);
    let resolver = PloidyDefinition::preset(PloidyPreset::Grch37)
        .resolver([SamplePloidy::sex("M"), SamplePloidy::sex("F")])
        .unwrap();
    let expected = sites
        .iter()
        .map(|&(reference, position)| {
            let ploidies = resolver.resolve(reference, position - 1).unwrap();
            (
                reference.to_owned(),
                position,
                genotype_for_ploidy(ploidies[0]).to_owned(),
                genotype_for_ploidy(ploidies[1]).to_owned(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(observed, expected);

    let custom = directory.path().join("ploidy.txt");
    fs::write(
        &custom,
        "X 1 60000 M 1\nX 2699521 154931043 M 1\nY 1 59373566 M 1\nY 1 59373566 F 0\nMT 1 16569 M 1\nMT 1 16569 F 1\nchrX 1 60000 M 1\nchrX 2699521 154931043 M 1\nchrY 1 59373566 M 1\nchrY 1 59373566 F 0\nchrM 1 16569 M 1\nchrM 1 16569 F 1\n* * * M 2\n* * * F 2\n",
    )
    .unwrap();
    let output = Command::new(bcftools())
        .args(["call", "-m", "--ploidy-file"])
        .arg(&custom)
        .arg("-S")
        .arg(&samples)
        .arg("-Ov")
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(called_genotypes(&output.stdout), expected);
}

fn called_genotypes(data: &[u8]) -> Vec<(String, u64, String, String)> {
    String::from_utf8(data.to_vec())
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            (
                fields[0].to_owned(),
                fields[1].parse::<u64>().unwrap(),
                fields[9].split(':').next().unwrap().to_owned(),
                fields[10].split(':').next().unwrap().to_owned(),
            )
        })
        .collect()
}

fn genotype_for_ploidy(ploidy: rsomics_call::CallPloidy) -> &'static str {
    match ploidy {
        rsomics_call::CallPloidy::Absent => ".",
        rsomics_call::CallPloidy::Haploid => "0",
        rsomics_call::CallPloidy::Diploid => "0/0",
    }
}

fn bcftools_indel(reference: &Path, alignments: &[PathBuf]) -> LikelihoodSite {
    let mut command = Command::new(bcftools());
    command
        .args(["mpileup", "-f"])
        .arg(reference)
        .args(["-a", ANNOTATIONS, "-Ou"]);
    command.args(alignments);
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut reader = LikelihoodVariantReader::new(&output.stdout[..]).unwrap();
    loop {
        let site = reader.read_site().unwrap().unwrap();
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
fn complete_snp_annotations_match_bcftools_1_24() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let reference = write_reference(directory.path(), "AAAAAAAAA");
    let alignment = directory.path().join("annotations.sam");
    let mut records =
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:MX\tLN:9\n@RG\tID:rg\tSM:sample\n".to_owned();
    records
        .push_str("rf1\t0\tMX\t1\t60\t2S9M\t*\t0\t0\tTTAAAAAAAAA\tIIIIIIIIIII\tRG:Z:rg\tNM:i:0\n");
    for index in 2..=4 {
        records.push_str(&format!(
            "rf{index}\t0\tMX\t1\t60\t9M\t*\t0\t0\tAAAAAAAAA\tIIIIIIIII\tRG:Z:rg\tNM:i:0\n"
        ));
    }
    records.push_str("rr\t16\tMX\t1\t20\t9M\t*\t0\t0\tAAAAAAAAA\tIIIIIIIII\tRG:Z:rg\tNM:i:3\n");
    records.push_str("af\t0\tMX\t1\t50\t9M\t*\t0\t0\tAAAACAAAA\t555555555\tRG:Z:rg\tNM:i:1\n");
    records
        .push_str("ar1\t16\tMX\t1\t30\t9M2S\t*\t0\t0\tAAAACAAAATT\t55555555555\tRG:Z:rg\tNM:i:4\n");
    for index in 2..=4 {
        records.push_str(&format!(
            "ar{index}\t16\tMX\t1\t30\t9M\t*\t0\t0\tAAAACAAAA\t555555555\tRG:Z:rg\tNM:i:4\n"
        ));
    }
    fs::write(&alignment, records).unwrap();

    let output = Command::new(bcftools())
        .args(["mpileup", "-B", "-f"])
        .arg(&reference)
        .args(["-a", ANNOTATIONS, "-Ou"])
        .arg(&alignment)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let expected = decode_likelihoods(&output.stdout);

    let run = SnpLikelihoodRun::open(
        [AlignmentInput::new(0, &alignment, "annotations")],
        &reference,
        SampleSelection::default(),
        PileupOptions::default(),
        SnpLikelihoodConfig::default(),
    )
    .unwrap();
    let mut actual = Vec::new();
    run.run(|site| {
        actual.push(site);
        Ok(())
    })
    .unwrap();
    assert_sites_match(&actual, &expected);
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
        assert_sites_match(
            &[rsomics_indel(&reference, &alignments)],
            &[bcftools_indel(&reference, &alignments)],
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
    assert_sites_match(
        &[rsomics_indel(
            &repeat_reference,
            std::slice::from_ref(&repeat_insertion),
        )],
        &[bcftools_indel(&repeat_reference, &[repeat_insertion])],
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
        .args(["-a", ANNOTATIONS, "-r", regions, "-t", targets, "-Ou"])
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

    assert_sites_match(&actual, &expected);
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
            .args(["-a", ANNOTATIONS, "-T"])
            .arg(&targets)
            .arg("-Ou")
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

        assert_sites_match(&actual, &expected);
    }
}
