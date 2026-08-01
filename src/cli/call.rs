use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Context, Result as CommonResult, RsomicsError};

use crate::{
    CallModel, CallOutputOptions, CallSampleSelection, ConsensusCallerConfig,
    IndexedLikelihoodVariantReader, LikelihoodCallRun, LikelihoodVariantReader,
    MultiallelicCallerConfig, PloidyDefinition, PloidyPreset, PloidyResolver,
};

use super::output::{VariantOutputArgs, call_result};
use super::pileup::parse_regions;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    #[command(flatten)]
    input: LikelihoodInputArgs,

    #[command(flatten)]
    samples: CallSampleArgs,

    #[command(flatten)]
    pub(crate) policy: CallPolicyArgs,

    #[command(flatten)]
    output: VariantOutputArgs,
}

impl Arguments {
    pub(crate) fn execute(&self, json: bool) -> CommonResult<()> {
        self.policy
            .validate(self.input.prior_frequencies.is_some())?;
        if self.input.has_regions() {
            self.execute_indexed(json)
        } else {
            self.execute_streaming(json)
        }
    }

    fn execute_streaming(&self, json: bool) -> CommonResult<()> {
        let input = open_likelihood_input(&self.input.input)?;
        let reader = call_result(LikelihoodVariantReader::new(input))?;
        let reader = call_result(self.input.configure_reader(reader))?;
        let definition = self.policy.ploidy_definition()?;
        let selection = self.samples.selection()?;
        let (reader, resolver) = call_result(selection.bind(reader, &definition))?;
        let sample_names = reader
            .schema()
            .header()
            .sample_names()
            .iter()
            .map(|name| name.as_str().to_owned())
            .collect::<Vec<_>>();
        let groups = self.policy.sample_groups(&sample_names)?;
        let calls = call_result(self.policy.build(resolver, groups))?;
        self.output.write(json, self.input_paths(), |output| {
            calls
                .run(reader, output, self.output.output_type.format())
                .map(|_| ())
        })
    }

    fn execute_indexed(&self, json: bool) -> CommonResult<()> {
        if self.input.input == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "indexed regions require a named likelihood input".to_owned(),
            ));
        }
        let reader = call_result(IndexedLikelihoodVariantReader::open(&self.input.input))?;
        let reader = call_result(self.input.configure_indexed_reader(reader))?;
        let definition = self.policy.ploidy_definition()?;
        let selection = self.samples.selection()?;
        let (reader, resolver) = call_result(selection.bind_indexed(reader, &definition))?;
        let sample_names = reader
            .schema()
            .header()
            .sample_names()
            .iter()
            .map(|name| name.as_str().to_owned())
            .collect::<Vec<_>>();
        let groups = self.policy.sample_groups(&sample_names)?;
        let calls = call_result(self.policy.build(resolver, groups))?;
        self.output.write(json, self.input_paths(), |output| {
            if let Some(regions) = self.input.regions.as_deref() {
                calls
                    .run_indexed(
                        reader,
                        parse_regions(regions)?,
                        output,
                        self.output.output_type.format(),
                    )
                    .map(|_| ())
            } else if let Some(regions) = &self.input.regions_file {
                calls
                    .run_indexed_regions_file(
                        reader,
                        regions,
                        output,
                        self.output.output_type.format(),
                    )
                    .map(|_| ())
            } else {
                unreachable!("indexed execution requires a region source")
            }
        })
    }

    fn input_paths(&self) -> Vec<PathBuf> {
        let mut paths = vec![self.input.input.clone()];
        paths.extend(self.input.regions_file.iter().cloned());
        paths.extend(self.samples.samples_file.iter().cloned());
        paths.extend(self.policy.ploidy_file.iter().cloned());
        if let Some(groups) = self
            .policy
            .group_samples
            .as_ref()
            .filter(|path| *path != Path::new("-"))
        {
            paths.push(groups.clone());
        }
        paths
    }
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Likelihood input")]
struct LikelihoodInputArgs {
    /// Likelihood VCF or BCF; omit or use - for standard input
    #[arg(value_name = "LIKELIHOODS", default_value = "-")]
    input: PathBuf,

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

    /// INFO tags containing total and alternate panel allele counts
    #[arg(short = 'F', long = "prior-freqs", value_name = "AN,AC")]
    prior_frequencies: Option<String>,
}

impl LikelihoodInputArgs {
    fn has_regions(&self) -> bool {
        self.regions.is_some() || self.regions_file.is_some()
    }

    fn configure_reader<R: Read>(
        &self,
        reader: LikelihoodVariantReader<R>,
    ) -> crate::Result<LikelihoodVariantReader<R>> {
        match self.prior_tags()? {
            Some((total, alternates)) => reader.with_prior_frequencies(total, alternates),
            None => Ok(reader),
        }
    }

    fn configure_indexed_reader(
        &self,
        reader: IndexedLikelihoodVariantReader,
    ) -> crate::Result<IndexedLikelihoodVariantReader> {
        match self.prior_tags()? {
            Some((total, alternates)) => reader.with_prior_frequencies(total, alternates),
            None => Ok(reader),
        }
    }

    fn prior_tags(&self) -> crate::Result<Option<(String, String)>> {
        let Some(value) = self.prior_frequencies.as_deref() else {
            return Ok(None);
        };
        let Some((total, alternates)) = value.split_once(',') else {
            return Err(crate::CallError::InvalidPriorAlleleCounts);
        };
        if total.is_empty() || alternates.is_empty() || alternates.contains(',') {
            return Err(crate::CallError::InvalidPriorAlleleCounts);
        }
        Ok(Some((total.to_owned(), alternates.to_owned())))
    }
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Call sample selection")]
struct CallSampleArgs {
    /// Comma-separated sample names to include; prefix with ^ to exclude
    #[arg(
        short = 's',
        long,
        value_name = "LIST",
        conflicts_with = "samples_file"
    )]
    samples: Option<String>,

    /// File containing exactly one selected sample name per nonempty line
    #[arg(short = 'S', long, value_name = "FILE")]
    samples_file: Option<PathBuf>,
}

impl CallSampleArgs {
    fn selection(&self) -> CommonResult<CallSampleSelection> {
        let records = if let Some(samples) = self.samples.as_deref() {
            samples.split(',').map(str::to_owned).collect::<Vec<_>>()
        } else if let Some(file) = &self.samples_file {
            fs::read_to_string(file)
                .rs_with_context(|| format!("reading call sample list {}", file.display()))?
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty() && !line.starts_with('#'))
                .map(str::to_owned)
                .collect()
        } else {
            return Ok(CallSampleSelection::default());
        };
        if records.is_empty() || records.iter().any(String::is_empty) {
            return Err(RsomicsError::InvalidInput(
                "call sample selection is empty".to_owned(),
            ));
        }
        if records
            .iter()
            .any(|row| row.split_whitespace().count() != 1)
        {
            return Err(RsomicsError::InvalidInput(
                "call sample files currently require exactly one sample name per row".to_owned(),
            ));
        }
        let exclude = records[0].starts_with('^');
        if exclude {
            let names = records
                .into_iter()
                .enumerate()
                .map(|(index, name)| {
                    if index == 0 {
                        name.strip_prefix('^').unwrap_or(&name).to_owned()
                    } else {
                        name
                    }
                })
                .collect::<Vec<_>>();
            call_result(CallSampleSelection::exclude(names))
        } else {
            call_result(CallSampleSelection::include(records))
        }
    }
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Calling policy")]
pub(crate) struct CallPolicyArgs {
    /// Calling model
    #[arg(long, value_enum, default_value = "multiallelic")]
    model: CallingModel,

    /// Mutation rate for the multiallelic model
    #[arg(short = 'P', long, value_name = "FLOAT")]
    mutation_rate: Option<f64>,

    /// Reference-posterior threshold for the consensus model
    #[arg(long, value_name = "FLOAT")]
    reference_probability_threshold: Option<f64>,

    /// Retain all input alternate alleles at emitted sites
    #[arg(long = "keep-alts")]
    keep_alternates: bool,

    /// Built-in ploidy definition
    #[arg(
        long,
        value_enum,
        value_name = "PRESET",
        default_value = "diploid",
        conflicts_with = "ploidy_file"
    )]
    ploidy: PloidyChoice,

    /// Custom CHROM FROM TO SEX PLOIDY definition
    #[arg(long, value_name = "FILE")]
    ploidy_file: Option<PathBuf>,

    /// Group selected samples using SAMPLE GROUP rows; - gives each sample its own group
    #[arg(short = 'G', long, value_name = "FILE")]
    group_samples: Option<PathBuf>,

    /// Emit only variant calls
    #[arg(short = 'v', long)]
    variants_only: bool,

    /// Retain sites whose reference allele contains N
    #[arg(long)]
    keep_masked_reference: bool,

    /// Skip one called variant class
    #[arg(short = 'V', long, value_enum, value_name = "TYPE")]
    skip_variants: Option<SkipVariant>,

    /// Strictly increasing minimum per-sample depths for gVCF blocks
    #[arg(short = 'g', long, value_name = "DEPTHS")]
    gvcf: Option<String>,
}

impl CallPolicyArgs {
    pub(crate) fn validate(&self, prior_frequencies: bool) -> CommonResult<()> {
        match self.model {
            CallingModel::Multiallelic => {
                if self.reference_probability_threshold.is_some() {
                    return Err(RsomicsError::ConfigError(
                        "--reference-probability-threshold requires --model consensus".to_owned(),
                    ));
                }
            }
            CallingModel::Consensus => {
                for (enabled, option) in [
                    (self.mutation_rate.is_some(), "--mutation-rate"),
                    (self.keep_alternates, "--keep-alts"),
                    (self.group_samples.is_some(), "--group-samples"),
                    (self.gvcf.is_some(), "--gvcf"),
                    (prior_frequencies, "--prior-freqs"),
                ] {
                    if enabled {
                        return Err(RsomicsError::ConfigError(format!(
                            "{option} requires --model multiallelic"
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    pub(crate) fn ploidy_definition(&self) -> CommonResult<PloidyDefinition> {
        if let Some(file) = &self.ploidy_file {
            call_result(PloidyDefinition::read(file))
        } else {
            Ok(PloidyDefinition::preset(self.ploidy.into()))
        }
    }

    pub(crate) fn build(
        &self,
        resolver: PloidyResolver,
        groups: Option<Vec<usize>>,
    ) -> crate::Result<LikelihoodCallRun> {
        let model = match self.model {
            CallingModel::Multiallelic => {
                let config = MultiallelicCallerConfig::new(self.mutation_rate.unwrap_or(1.1e-3))?
                    .with_keep_alternates(self.keep_alternates);
                CallModel::Multiallelic(config)
            }
            CallingModel::Consensus => CallModel::Consensus(ConsensusCallerConfig::new(
                self.reference_probability_threshold.unwrap_or(0.5),
            )?),
        };
        let output = CallOutputOptions::default()
            .with_variants_only(self.variants_only)
            .with_keep_masked_reference(self.keep_masked_reference)
            .with_skip_snps(self.skip_variants == Some(SkipVariant::Snps))
            .with_skip_indels(self.skip_variants == Some(SkipVariant::Indels));
        let mut run = LikelihoodCallRun::new(model, resolver).with_output_options(output);
        if let Some(thresholds) = self.gvcf_thresholds()? {
            run = run.with_gvcf(thresholds)?;
        }
        if let Some(groups) = groups {
            run = run.with_sample_groups(groups)?;
        }
        Ok(run)
    }

    pub(crate) fn sample_groups(
        &self,
        sample_names: &[String],
    ) -> CommonResult<Option<Vec<usize>>> {
        let Some(path) = &self.group_samples else {
            return Ok(None);
        };
        if path == Path::new("-") {
            return Ok(Some((0..sample_names.len()).collect()));
        }
        let data = fs::read_to_string(path)
            .rs_with_context(|| format!("reading sample groups {}", path.display()))?;
        let selected = sample_names
            .iter()
            .map(String::as_str)
            .collect::<HashSet<_>>();
        let mut assignments = HashMap::new();
        for (index, line) in data.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split_whitespace().collect::<Vec<_>>();
            let [sample, group] = fields.as_slice() else {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample group {} line {} must contain SAMPLE and GROUP",
                    path.display(),
                    index + 1
                )));
            };
            if !selected.contains(sample) {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample group contains an unselected sample: {sample}"
                )));
            }
            if assignments
                .insert((*sample).to_owned(), (*group).to_owned())
                .is_some()
            {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample group contains a duplicate sample: {sample}"
                )));
            }
        }
        if assignments.len() != sample_names.len() {
            return Err(RsomicsError::InvalidInput(
                "sample groups must assign every selected sample".to_owned(),
            ));
        }
        let mut group_indices = HashMap::new();
        let mut next = 0;
        Ok(Some(
            sample_names
                .iter()
                .map(|sample| {
                    let group = assignments
                        .get(sample)
                        .expect("every selected sample was checked");
                    *group_indices.entry(group).or_insert_with(|| {
                        let index = next;
                        next += 1;
                        index
                    })
                })
                .collect(),
        ))
    }

    fn gvcf_thresholds(&self) -> crate::Result<Option<Vec<u32>>> {
        self.gvcf
            .as_deref()
            .map(|value| {
                value
                    .split(',')
                    .map(|value| {
                        value
                            .parse::<u32>()
                            .map_err(|_| crate::CallError::InvalidGvcfThresholds)
                    })
                    .collect::<crate::Result<Vec<_>>>()
            })
            .transpose()
    }
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum CallingModel {
    #[default]
    Multiallelic,
    Consensus,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum PloidyChoice {
    Grch37,
    Grch38,
    Haploid,
    #[default]
    Diploid,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum SkipVariant {
    Snps,
    Indels,
}

impl From<PloidyChoice> for PloidyPreset {
    fn from(value: PloidyChoice) -> Self {
        match value {
            PloidyChoice::Grch37 => Self::Grch37,
            PloidyChoice::Grch38 => Self::Grch38,
            PloidyChoice::Haploid => Self::Haploid,
            PloidyChoice::Diploid => Self::Diploid,
        }
    }
}

fn open_likelihood_input(path: &Path) -> CommonResult<Box<dyn Read>> {
    if path == Path::new("-") {
        Ok(Box::new(io::stdin()))
    } else {
        File::open(path)
            .rs_with_context(|| format!("opening likelihood input {}", path.display()))
            .map(|file| Box::new(BufReader::new(file)) as Box<dyn Read>)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn group_indices_follow_selected_sample_order() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("groups.tsv");
        fs::write(&path, "S2\tB\nS1\tA\nS3\tB\n").unwrap();
        let policy = CallPolicyArgs {
            model: CallingModel::Multiallelic,
            mutation_rate: None,
            reference_probability_threshold: None,
            keep_alternates: false,
            ploidy: PloidyChoice::Diploid,
            ploidy_file: None,
            group_samples: Some(path),
            variants_only: false,
            keep_masked_reference: false,
            skip_variants: None,
            gvcf: None,
        };
        assert_eq!(
            policy
                .sample_groups(&["S1".into(), "S2".into(), "S3".into()])
                .unwrap(),
            Some(vec![0, 1, 1])
        );
    }
}
