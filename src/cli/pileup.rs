use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use clap::{Args, ValueEnum};
use noodles::core::Region;
use rsomics_common::{Context, Result as CommonResult, RsomicsError};
use rsomics_pileup::{FlagFilter, PileupOptions, RecordFilter};

use crate::{
    AlignmentInput, IndelAmbiguousReadPolicy, IndelLikelihoodConfig, LikelihoodVariantWriter,
    LikelihoodVcfSchema, SampleSelection, SnpLikelihoodConfig, SnpLikelihoodRun,
};

use super::output::{VariantOutputArgs, call_result};

#[cfg(test)]
const DEFAULT_SKIP_FLAGS: u16 = 0x704;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    #[command(flatten)]
    pub(crate) input: AlignmentArgs,

    #[command(flatten)]
    pub(crate) policy: LikelihoodArgs,

    #[command(flatten)]
    output: VariantOutputArgs,
}

impl Arguments {
    pub(crate) fn execute(&self, json: bool) -> CommonResult<()> {
        let prepared = self.input.prepare()?;
        let input_paths = self.input.input_paths(&prepared);
        self.output.write(json, input_paths, |output| {
            let run = self.input.open(prepared, &self.policy)?;
            let schema = LikelihoodVcfSchema::new(
                run.reference_sequences()
                    .iter()
                    .map(|reference| (reference.name(), reference.length())),
                run.samples().samples(),
            )?;
            let mut writer =
                LikelihoodVariantWriter::new(output, schema, self.output.output_type.format())?;
            run.run(|site| writer.write_site(&site))?;
            writer.finish()?;
            Ok(())
        })
    }
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Alignment input")]
pub(crate) struct AlignmentArgs {
    /// Coordinate-sorted SAM, BAM, or CRAM inputs
    #[arg(value_name = "ALIGNMENT")]
    alignments: Vec<PathBuf>,

    /// File containing one alignment path per nonempty line
    #[arg(short = 'b', long = "alignment-list", value_name = "FILE")]
    alignment_list: Option<PathBuf>,

    /// Indexed FASTA reference
    #[arg(
        short = 'f',
        long = "reference",
        value_name = "FASTA",
        conflicts_with = "no_reference",
        required_unless_present = "no_reference"
    )]
    reference: Option<PathBuf>,

    /// Generate reference-free SNP likelihoods
    #[arg(long, conflicts_with = "reference")]
    no_reference: bool,

    /// Treat each alignment input as one sample and ignore read groups
    #[arg(long)]
    ignore_read_groups: bool,

    /// Comma-separated samples to include; prefix with ^ to exclude
    #[arg(
        short = 's',
        long,
        value_name = "LIST",
        conflicts_with = "samples_file"
    )]
    samples: Option<String>,

    /// Sample names, or input and output names, one row per sample
    #[arg(short = 'S', long, value_name = "FILE")]
    samples_file: Option<PathBuf>,

    /// Indexed regions separated by commas
    #[arg(
        short = 'r',
        long,
        value_name = "REGIONS",
        conflicts_with = "regions_file"
    )]
    regions: Option<String>,

    /// Indexed regions read from a file
    #[arg(short = 'R', long, value_name = "FILE")]
    regions_file: Option<PathBuf>,

    /// Streaming target intervals separated by commas
    #[arg(
        short = 't',
        long,
        value_name = "REGIONS",
        conflicts_with = "targets_file"
    )]
    targets: Option<String>,

    /// Streaming target intervals read from a file
    #[arg(short = 'T', long, value_name = "FILE")]
    targets_file: Option<PathBuf>,
}

pub(crate) struct PreparedAlignments {
    inputs: Vec<AlignmentInput>,
    paths: Vec<PathBuf>,
    samples: SampleSelection,
}

impl AlignmentArgs {
    pub(crate) fn prepare(&self) -> CommonResult<PreparedAlignments> {
        let mut paths = self.alignments.clone();
        if let Some(list) = &self.alignment_list {
            let data = fs::read_to_string(list)
                .rs_with_context(|| format!("reading alignment list {}", list.display()))?;
            paths.extend(
                data.lines()
                    .map(str::trim_end)
                    .filter(|line| !line.is_empty())
                    .map(PathBuf::from),
            );
        }
        if paths.is_empty() {
            return Err(RsomicsError::ConfigError(
                "at least one alignment input or --alignment-list is required".to_owned(),
            ));
        }
        let inputs = paths
            .iter()
            .enumerate()
            .map(|(index, path)| {
                let source_id = u32::try_from(index + 1).map_err(|_| {
                    RsomicsError::ConfigError("too many alignment inputs".to_owned())
                })?;
                Ok(
                    AlignmentInput::new(source_id, path, fallback_sample_name(path))
                        .ignore_read_groups(self.ignore_read_groups),
                )
            })
            .collect::<CommonResult<Vec<_>>>()?;
        let samples = read_pileup_samples(self.samples.as_deref(), self.samples_file.as_deref())?;
        Ok(PreparedAlignments {
            inputs,
            paths,
            samples,
        })
    }

    pub(crate) fn input_paths(&self, prepared: &PreparedAlignments) -> Vec<PathBuf> {
        let mut paths = prepared.paths.clone();
        paths.extend(self.alignment_list.iter().cloned());
        paths.extend(self.reference.iter().cloned());
        paths.extend(self.samples_file.iter().cloned());
        paths.extend(self.regions_file.iter().cloned());
        paths.extend(self.targets_file.iter().cloned());
        paths
    }

    pub(crate) fn open(
        &self,
        prepared: PreparedAlignments,
        likelihood: &LikelihoodArgs,
    ) -> crate::Result<SnpLikelihoodRun> {
        debug_assert!(self.reference.is_some() || self.no_reference);
        if likelihood.indel_options_changed() {
            if self.reference.is_none() {
                return Err(crate::CallError::MissingLikelihoodReference(
                    "indel likelihood options",
                ));
            }
            if likelihood.skip_indels {
                return Err(crate::CallError::InvalidIndelConfiguration(
                    "indel options cannot be combined with --skip-indels",
                ));
            }
        }
        let pileup = likelihood.pileup_options();
        let snp = likelihood.snp_config()?;
        let regions = self.regions.as_deref().map(parse_regions).transpose()?;
        let mut run = match (
            self.reference.as_deref(),
            regions,
            self.regions_file.as_deref(),
        ) {
            (Some(reference), Some(regions), None) => SnpLikelihoodRun::open_regions(
                prepared.inputs,
                reference,
                prepared.samples,
                regions,
                pileup,
                snp,
            )?,
            (Some(reference), None, Some(regions)) => SnpLikelihoodRun::open_regions_file(
                prepared.inputs,
                reference,
                prepared.samples,
                regions,
                pileup,
                snp,
            )?,
            (Some(reference), None, None) => {
                SnpLikelihoodRun::open(prepared.inputs, reference, prepared.samples, pileup, snp)?
            }
            (None, Some(regions), None) => SnpLikelihoodRun::open_regions_without_reference(
                prepared.inputs,
                prepared.samples,
                regions,
                pileup,
                snp,
            )?,
            (None, None, Some(regions)) => SnpLikelihoodRun::open_regions_file_without_reference(
                prepared.inputs,
                prepared.samples,
                regions,
                pileup,
                snp,
            )?,
            (None, None, None) => SnpLikelihoodRun::open_without_reference(
                prepared.inputs,
                prepared.samples,
                pileup,
                snp,
            )?,
            (_, Some(_), Some(_)) => unreachable!("Clap rejects two region sources"),
        };

        let baq = likelihood.baq.unwrap_or(if self.reference.is_some() {
            BaqMode::Partial
        } else {
            BaqMode::Off
        });
        run = match baq {
            BaqMode::Off => {
                if likelihood.redo_baq {
                    return Err(crate::CallError::InvalidBaqConfiguration(
                        "--redo-baq requires partial or full BAQ",
                    ));
                }
                run
            }
            BaqMode::Partial => {
                if likelihood.maximum_read_length == 0 {
                    return Err(crate::CallError::InvalidBaqConfiguration(
                        "maximum read length must be greater than zero",
                    ));
                }
                run.with_partial_baq(likelihood.maximum_read_length, likelihood.redo_baq)?
            }
            BaqMode::Full => {
                if likelihood.maximum_read_length == 0 {
                    return Err(crate::CallError::InvalidBaqConfiguration(
                        "maximum read length must be greater than zero",
                    ));
                }
                run.with_full_baq(likelihood.maximum_read_length, likelihood.redo_baq)?
            }
        };
        if self.reference.is_some() && !likelihood.skip_indels {
            run = run.with_indels(likelihood.indel_config())?;
        }
        if let Some(targets) = self.targets.as_deref() {
            run = run.with_targets(parse_regions(targets)?);
        }
        if let Some(targets) = &self.targets_file {
            run = run.with_target_file(targets)?;
        }
        Ok(run)
    }
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Pileup and likelihood policy")]
pub(crate) struct LikelihoodArgs {
    /// Skip records when all of these FLAG bits are set
    #[arg(long, value_name = "FLAGS", value_parser = parse_u16, default_value = "0")]
    skip_all_set: u16,

    /// Skip records when any of these FLAG bits are set
    #[arg(long, value_name = "FLAGS", value_parser = parse_u16, default_value = "0x704")]
    skip_any_set: u16,

    /// Skip records when all of these FLAG bits are unset
    #[arg(long, value_name = "FLAGS", value_parser = parse_u16, default_value = "0")]
    skip_all_unset: u16,

    /// Skip records when any of these FLAG bits are unset
    #[arg(long, value_name = "FLAGS", value_parser = parse_u16, default_value = "0")]
    skip_any_unset: u16,

    /// Minimum read mapping quality
    #[arg(short = 'q', long = "min-mq", value_name = "INT", default_value_t = 0)]
    minimum_mapping_quality: u8,

    /// Include paired reads not marked as proper pairs
    #[arg(short = 'A', long = "count-orphans")]
    include_anomalous_pairs: bool,

    /// Disable overlapping-mate quality adjustment
    #[arg(short = 'x', long)]
    ignore_overlaps: bool,

    /// Maximum pileup depth per alignment input; 0 disables the limit
    #[arg(short = 'd', long, value_name = "INT", default_value_t = 250)]
    maximum_depth: usize,

    /// Minimum base quality
    #[arg(short = 'Q', long = "min-bq", value_name = "INT", default_value_t = 1)]
    minimum_base_quality: u8,

    /// Maximum base quality
    #[arg(long = "max-bq", value_name = "INT", default_value_t = 60)]
    maximum_base_quality: u8,

    /// Neighboring-base quality cap delta
    #[arg(long = "delta-bq", value_name = "INT", default_value_t = 30)]
    neighboring_quality_delta: u8,

    /// Mapping-quality cap used by likelihood scoring
    #[arg(long = "likelihood-mq-cap", value_name = "INT", default_value_t = 60)]
    mapping_quality_cap: u8,

    /// Deterministic likelihood sampling seed
    #[arg(long, value_name = "INT", default_value_t = 0)]
    seed: i32,

    /// BAQ mode; defaults to partial with a reference and off without one
    #[arg(long, value_enum, value_name = "MODE")]
    baq: Option<BaqMode>,

    /// Recalculate BAQ instead of using existing BQ tags
    #[arg(short = 'E', long)]
    redo_baq: bool,

    /// Maximum read length considered by BAQ
    #[arg(short = 'M', long, value_name = "INT", default_value_t = 500)]
    maximum_read_length: usize,

    /// Disable indel likelihood generation
    #[arg(short = 'I', long)]
    skip_indels: bool,

    /// Minimum reads supporting an indel candidate
    #[arg(short = 'm', long, value_name = "INT", default_value_t = 2)]
    minimum_indel_support: usize,

    /// Minimum fraction supporting an indel candidate
    #[arg(short = 'F', long, value_name = "FLOAT", default_value_t = 0.05)]
    minimum_indel_fraction: f64,

    /// Maximum depth used for indel likelihoods
    #[arg(short = 'L', long, value_name = "INT", default_value_t = 250)]
    maximum_indel_depth: usize,

    /// Indel realignment window size
    #[arg(long, value_name = "INT", default_value_t = 110)]
    indel_window: usize,

    /// Indel gap-open quality
    #[arg(long, value_name = "INT", default_value_t = 40)]
    gap_open_quality: i32,

    /// Indel gap-extension quality
    #[arg(short = 'e', long, value_name = "INT", default_value_t = 20)]
    gap_extension_quality: i32,

    /// Tandem-repeat indel quality cap
    #[arg(long, value_name = "INT", default_value_t = 500)]
    tandem_quality: i32,

    /// Relative deletion-versus-insertion likelihood bias
    #[arg(long, value_name = "FLOAT", default_value_t = 1.0)]
    indel_bias: f64,

    /// Evaluate minimum indel support independently in each sample
    #[arg(short = 'p', long)]
    per_sample_indel_support: bool,

    /// Ambiguous indel-read allele-depth policy
    #[arg(long, value_enum, value_name = "MODE", default_value = "drop")]
    ambiguous_reads: AmbiguousReads,
}

impl LikelihoodArgs {
    fn indel_options_changed(&self) -> bool {
        self.minimum_indel_support != 2
            || self.minimum_indel_fraction != 0.05
            || self.maximum_indel_depth != 250
            || self.indel_window != 110
            || self.gap_open_quality != 40
            || self.gap_extension_quality != 20
            || self.tandem_quality != 500
            || self.indel_bias != 1.0
            || self.per_sample_indel_support
            || !matches!(self.ambiguous_reads, AmbiguousReads::Drop)
    }

    fn pileup_options(&self) -> PileupOptions {
        PileupOptions {
            filter: RecordFilter {
                flags: FlagFilter {
                    skip_all_set: self.skip_all_set,
                    skip_any_set: self.skip_any_set,
                    skip_all_unset: self.skip_all_unset,
                    skip_any_unset: self.skip_any_unset,
                },
                minimum_mapping_quality: self.minimum_mapping_quality,
                include_anomalous_pairs: self.include_anomalous_pairs,
            },
            adjust_overlaps: !self.ignore_overlaps,
            maximum_depth_per_source: (self.maximum_depth != 0).then_some(self.maximum_depth),
        }
    }

    fn snp_config(&self) -> crate::Result<SnpLikelihoodConfig> {
        Ok(SnpLikelihoodConfig::new(
            self.minimum_base_quality,
            self.maximum_base_quality,
            self.neighboring_quality_delta,
            self.mapping_quality_cap,
        )?
        .with_random_seed(self.seed))
    }

    fn indel_config(&self) -> IndelLikelihoodConfig {
        IndelLikelihoodConfig::default()
            .with_minimum_support(self.minimum_indel_support)
            .with_minimum_fraction(self.minimum_indel_fraction)
            .with_maximum_depth(self.maximum_indel_depth)
            .with_window_size(self.indel_window)
            .with_gap_open_quality(self.gap_open_quality)
            .with_gap_extension_quality(self.gap_extension_quality)
            .with_tandem_quality(self.tandem_quality)
            .with_minimum_base_quality(self.minimum_base_quality)
            .with_mapping_quality_cap(self.mapping_quality_cap)
            .with_indel_bias(self.indel_bias)
            .with_random_seed(self.seed)
            .with_per_sample_support(self.per_sample_indel_support)
            .with_ambiguous_read_policy(self.ambiguous_reads.into())
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum BaqMode {
    Off,
    Partial,
    Full,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum AmbiguousReads {
    #[default]
    Drop,
    Distribute,
    Reference,
}

impl From<AmbiguousReads> for IndelAmbiguousReadPolicy {
    fn from(value: AmbiguousReads) -> Self {
        match value {
            AmbiguousReads::Drop => Self::Drop,
            AmbiguousReads::Distribute => Self::DistributeAlleleDepth,
            AmbiguousReads::Reference => Self::AddToReferenceAlleleDepth,
        }
    }
}

fn parse_u16(value: &str) -> Result<u16, String> {
    let parsed = value
        .strip_prefix("0x")
        .or_else(|| value.strip_prefix("0X"))
        .map_or_else(
            || value.parse::<u16>(),
            |value| u16::from_str_radix(value, 16),
        );
    parsed.map_err(|error| error.to_string())
}

pub(crate) fn parse_regions(value: &str) -> crate::Result<Vec<Region>> {
    value
        .split(',')
        .map(|raw| {
            Region::from_str(raw).map_err(|error| crate::CallError::InvalidRegion {
                region: raw.to_owned(),
                message: error.to_string(),
            })
        })
        .collect()
}

fn fallback_sample_name(path: &Path) -> String {
    path.file_stem()
        .and_then(|name| name.to_str())
        .filter(|name| !name.is_empty())
        .map_or_else(|| path.display().to_string(), str::to_owned)
}

fn read_pileup_samples(list: Option<&str>, file: Option<&Path>) -> CommonResult<SampleSelection> {
    let records = if let Some(list) = list {
        list.split(',').map(str::to_owned).collect::<Vec<_>>()
    } else if let Some(file) = file {
        fs::read_to_string(file)
            .rs_with_context(|| format!("reading sample list {}", file.display()))?
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(str::to_owned)
            .collect()
    } else {
        return Ok(SampleSelection::default());
    };
    if records.is_empty() || records.iter().any(String::is_empty) {
        return Err(RsomicsError::InvalidInput(
            "sample selection is empty".to_owned(),
        ));
    }
    let exclude = records[0].starts_with('^');
    if exclude {
        let names = records
            .into_iter()
            .enumerate()
            .map(|(index, name)| {
                let name = if index == 0 {
                    name.strip_prefix('^').unwrap_or(&name).to_owned()
                } else {
                    name
                };
                if name.contains('=') || name.split_whitespace().count() != 1 {
                    return Err(RsomicsError::InvalidInput(
                        "excluded samples cannot be renamed".to_owned(),
                    ));
                }
                Ok(name)
            })
            .collect::<CommonResult<Vec<_>>>()?;
        return call_result(SampleSelection::exclude(names));
    }
    let samples = records
        .into_iter()
        .map(|record| {
            let fields = record.split_whitespace().collect::<Vec<_>>();
            let (input, output) = match fields.as_slice() {
                [single] => single
                    .split_once('=')
                    .map_or((*single, None), |(input, output)| (input, Some(output))),
                [input, output] => (*input, Some(*output)),
                _ => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "invalid sample selection row: {record}"
                    )));
                }
            };
            Ok((input.to_owned(), output.map(str::to_owned)))
        })
        .collect::<CommonResult<Vec<_>>>()?;
    call_result(SampleSelection::include(samples))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_masks_accept_decimal_and_hexadecimal() {
        assert_eq!(parse_u16("1796").unwrap(), DEFAULT_SKIP_FLAGS);
        assert_eq!(parse_u16("0x704").unwrap(), DEFAULT_SKIP_FLAGS);
    }

    #[test]
    fn alignment_lists_preserve_paths_with_spaces() {
        let directory = tempfile::tempdir().unwrap();
        let list = directory.path().join("alignments.txt");
        fs::write(&list, "reads one.bam  \n\nreads-two.cram\n").unwrap();
        let arguments = AlignmentArgs {
            alignments: Vec::new(),
            alignment_list: Some(list),
            reference: None,
            no_reference: true,
            ignore_read_groups: false,
            samples: None,
            samples_file: None,
            regions: None,
            regions_file: None,
            targets: None,
            targets_file: None,
        };
        let prepared = arguments.prepare().unwrap();
        assert_eq!(
            prepared.paths,
            [
                PathBuf::from("reads one.bam"),
                PathBuf::from("reads-two.cram")
            ]
        );
    }
}
