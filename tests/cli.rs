use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> &'static str {
    env!("CARGO_BIN_EXE_rsomics-call")
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(arguments: &[&str]) -> Output {
    Command::new(binary()).args(arguments).output().unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn bcftools() -> String {
    std::env::var("BCFTOOLS").unwrap_or_else(|_| "bcftools".to_owned())
}

#[test]
fn top_level_and_workflow_help_share_the_product_tree() {
    let top = run(&["--help"]);
    assert_success(&top);
    let top = String::from_utf8(top.stdout).unwrap();
    for command in ["pileup", "call", "run"] {
        assert!(top.contains(command), "{top}");
    }

    for command in ["pileup", "call", "run"] {
        let output = run(&[command, "--help"]);
        assert_success(&output);
        let help = String::from_utf8(output.stdout).unwrap();
        assert!(help.contains("--output-type <TYPE>"), "{help}");
        assert!(help.contains("Global options:"), "{help}");
    }
}

#[test]
fn all_three_commands_form_one_byte_equivalent_pipeline() {
    let directory = tempfile::tempdir().unwrap();
    let likelihoods = directory.path().join("likelihoods.vcf");
    let materialized = directory.path().join("materialized.vcf");
    let fused = directory.path().join("fused.vcf");
    let reference = fixture("alignment-reference.fa");
    let alignment = fixture("alignment.cram");

    let output = Command::new(binary())
        .arg("pileup")
        .arg("--reference")
        .arg(&reference)
        .arg(&alignment)
        .arg("-Ov")
        .arg("--output")
        .arg(&likelihoods)
        .output()
        .unwrap();
    assert_success(&output);

    let output = Command::new(binary())
        .arg("call")
        .arg(&likelihoods)
        .arg("-Ov")
        .arg("--output")
        .arg(&materialized)
        .output()
        .unwrap();
    assert_success(&output);

    let output = Command::new(binary())
        .arg("run")
        .arg("--reference")
        .arg(&reference)
        .arg(&alignment)
        .arg("-Ov")
        .arg("--output")
        .arg(&fused)
        .output()
        .unwrap();
    assert_success(&output);

    assert_eq!(fs::read(materialized).unwrap(), fs::read(fused).unwrap());
}

#[test]
fn alignment_lists_feed_the_same_pileup_command() {
    let directory = tempfile::tempdir().unwrap();
    let list = directory.path().join("alignments.txt");
    let output_path = directory.path().join("likelihoods.vcf");
    fs::write(&list, format!("{}\n", fixture("alignment.cram").display())).unwrap();

    let output = Command::new(binary())
        .arg("pileup")
        .arg("--reference")
        .arg(fixture("alignment-reference.fa"))
        .arg("--alignment-list")
        .arg(&list)
        .arg("--output")
        .arg(&output_path)
        .output()
        .unwrap();
    assert_success(&output);
    let data = fs::read_to_string(output_path).unwrap();
    assert!(data.contains("#CHROM\tPOS"), "{data}");
    assert!(data.contains("chr1\t2\t"), "{data}");
}

#[test]
fn json_is_separate_from_named_variant_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("calls.vcf");
    let output = Command::new(binary())
        .arg("--json")
        .arg("call")
        .arg(fixture("bcftools-1.24-likelihood.vcf"))
        .arg("--output")
        .arg(&output_path)
        .output()
        .unwrap();
    assert_success(&output);
    let json = String::from_utf8(output.stdout).unwrap();
    assert!(json.contains("\"status\":\"ok\""), "{json}");
    assert!(
        fs::read_to_string(output_path)
            .unwrap()
            .starts_with("##fileformat=VCF")
    );

    let output = run(&[
        "--json",
        "call",
        fixture("bcftools-1.24-likelihood.vcf").to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(!String::from_utf8_lossy(&output.stdout).contains("##fileformat"));
}

#[test]
fn output_type_binding_writes_all_four_encodings() {
    let directory = tempfile::tempdir().unwrap();
    for (kind, prefix) in [
        ("v", b"##".as_slice()),
        ("z", &[0x1f, 0x8b]),
        ("u", b"BCF".as_slice()),
        ("b", &[0x1f, 0x8b]),
    ] {
        let path = directory.path().join(format!("calls.{kind}"));
        let output = Command::new(binary())
            .arg("call")
            .arg(fixture("bcftools-1.24-likelihood.vcf"))
            .arg("--output-type")
            .arg(kind)
            .arg("--output")
            .arg(&path)
            .output()
            .unwrap();
        assert_success(&output);
        assert!(fs::read(path).unwrap().starts_with(prefix));
    }
}

#[test]
fn failed_or_aliased_output_never_replaces_input_data() {
    let directory = tempfile::tempdir().unwrap();
    let invalid = directory.path().join("invalid.vcf");
    let output_path = directory.path().join("calls.vcf");
    fs::write(&invalid, b"not a variant file\n").unwrap();
    fs::write(&output_path, b"existing output\n").unwrap();
    let output = Command::new(binary())
        .arg("call")
        .arg(&invalid)
        .arg("--output")
        .arg(&output_path)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), b"existing output\n");

    let aliased = directory.path().join("likelihoods.vcf");
    fs::copy(fixture("bcftools-1.24-likelihood.vcf"), &aliased).unwrap();
    let before = fs::read(&aliased).unwrap();
    let output = Command::new(binary())
        .arg("call")
        .arg(&aliased)
        .arg("--output")
        .arg(&aliased)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(aliased).unwrap(), before);
}

#[test]
fn ambiguous_call_sample_rows_are_not_exposed() {
    let directory = tempfile::tempdir().unwrap();
    let samples = directory.path().join("samples.txt");
    fs::write(&samples, b"s1 1\n").unwrap();
    let output = Command::new(binary())
        .arg("call")
        .arg(fixture("bcftools-1.24-likelihood.vcf"))
        .arg("--samples-file")
        .arg(samples)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("exactly one sample name"), "{stderr}");
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn complete_command_defaults_match_bcftools_1_24() {
    let version = Command::new(bcftools()).arg("--version").output().unwrap();
    assert_success(&version);
    assert!(
        String::from_utf8(version.stdout)
            .unwrap()
            .starts_with("bcftools 1.24")
    );

    let directory = tempfile::tempdir().unwrap();
    let expected_likelihoods = directory.path().join("bcftools-likelihoods.vcf");
    let expected_calls = directory.path().join("bcftools-calls.vcf");
    let actual_likelihoods = directory.path().join("rsomics-likelihoods.vcf");
    let actual_calls = directory.path().join("rsomics-calls.vcf");
    let reference = fixture("alignment-reference.fa");
    let alignment = fixture("alignment.cram");
    let annotations = "FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/QM,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/FS,INFO/NMBZ,INFO/NM,INFO/SCR";

    let output = Command::new(bcftools())
        .arg("mpileup")
        .arg("--fasta-ref")
        .arg(&reference)
        .arg("--annotate")
        .arg(annotations)
        .arg("-Ov")
        .arg("--output")
        .arg(&expected_likelihoods)
        .arg(&alignment)
        .output()
        .unwrap();
    assert_success(&output);
    let output = Command::new(bcftools())
        .arg("call")
        .arg("--multiallelic-caller")
        .arg("-Ov")
        .arg("--output")
        .arg(&expected_calls)
        .arg(&expected_likelihoods)
        .output()
        .unwrap();
    assert_success(&output);

    let output = Command::new(binary())
        .arg("pileup")
        .arg("--reference")
        .arg(&reference)
        .arg(&alignment)
        .arg("-Ov")
        .arg("--output")
        .arg(&actual_likelihoods)
        .output()
        .unwrap();
    assert_success(&output);
    let output = Command::new(binary())
        .arg("run")
        .arg("--reference")
        .arg(&reference)
        .arg(&alignment)
        .arg("-Ov")
        .arg("--output")
        .arg(&actual_calls)
        .output()
        .unwrap();
    assert_success(&output);

    let likelihood_format = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/I16\t%INFO/MQ0F[\t%PL\t%DP\t%SP\t%ADF\t%ADR\t%AD\t%SCR\t%QS\t%QM]\n";
    let call_format = "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AN\t%AC[\t%GT\t%DP\t%AD]\n";
    for (expected, actual, format) in [
        (
            &expected_likelihoods,
            &actual_likelihoods,
            likelihood_format,
        ),
        (&expected_calls, &actual_calls, call_format),
    ] {
        let query = |path: &Path| {
            let output = Command::new(bcftools())
                .arg("query")
                .arg("--format")
                .arg(format)
                .arg(path)
                .output()
                .unwrap();
            assert_success(&output);
            output.stdout
        };
        assert_eq!(query(expected), query(actual));
    }
}
