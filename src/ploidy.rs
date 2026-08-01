use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::{CallError, CallPloidy, LikelihoodSite, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PloidyPreset {
    Grch37,
    Grch38,
    Haploid,
    Diploid,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SamplePloidy {
    Sex(String),
    Fixed(CallPloidy),
}

impl SamplePloidy {
    pub fn sex(value: impl Into<String>) -> Self {
        Self::Sex(value.into())
    }
}

#[derive(Clone, Debug)]
pub struct PloidyDefinition {
    defaults: BTreeMap<String, CallPloidy>,
    intervals: BTreeMap<String, BTreeMap<String, Vec<PloidyInterval>>>,
    default_sex: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct PloidyInterval {
    start: u64,
    end: u64,
    ploidy: CallPloidy,
}

impl PloidyDefinition {
    pub fn preset(preset: PloidyPreset) -> Self {
        match preset {
            PloidyPreset::Grch37 => Self::from_records(grch37_records(), Some("F".to_owned())),
            PloidyPreset::Grch38 => Self::from_records(grch38_records(), Some("F".to_owned())),
            PloidyPreset::Haploid => Self::constant(CallPloidy::Haploid),
            PloidyPreset::Diploid => Self::constant(CallPloidy::Diploid),
        }
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let display = path.display().to_string();
        let file = File::open(path).map_err(|error| CallError::PloidyInput {
            path: display.clone(),
            message: error.to_string(),
        })?;
        Self::from_reader(BufReader::new(file), display)
    }

    pub fn parse(data: &str) -> Result<Self> {
        Self::from_reader(data.as_bytes(), "<ploidy>".to_owned())
    }

    fn from_reader(reader: impl BufRead, path: String) -> Result<Self> {
        let mut records = Vec::new();
        for (index, line) in reader.lines().enumerate() {
            let line_number = index as u64 + 1;
            let line = line.map_err(|error| CallError::PloidyRecord {
                path: path.clone(),
                line: line_number,
                message: error.to_string(),
            })?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            records.push(parse_record(line, &path, line_number)?);
        }
        if records.is_empty() {
            return Err(CallError::PloidyInput {
                path,
                message: "no ploidy records".to_owned(),
            });
        }
        Self::build(records, None, &path)
    }

    fn constant(ploidy: CallPloidy) -> Self {
        Self::from_records(
            [ParsedRecord::Default {
                sex: "*".to_owned(),
                ploidy,
            }],
            Some("*".to_owned()),
        )
    }

    fn from_records(
        records: impl IntoIterator<Item = ParsedRecord>,
        default_sex: Option<String>,
    ) -> Self {
        Self::build(records, default_sex, "built-in preset").unwrap()
    }

    fn build(
        records: impl IntoIterator<Item = ParsedRecord>,
        requested_default_sex: Option<String>,
        path: &str,
    ) -> Result<Self> {
        let mut explicit_defaults = BTreeMap::new();
        let mut intervals = BTreeMap::<String, BTreeMap<String, Vec<PloidyInterval>>>::new();
        let mut sexes = BTreeSet::new();
        for record in records {
            match record {
                ParsedRecord::Default { sex, ploidy } => {
                    if explicit_defaults.insert(sex.clone(), ploidy).is_some() {
                        return Err(CallError::PloidyInput {
                            path: path.to_owned(),
                            message: format!("duplicate default for sex {sex}"),
                        });
                    }
                    sexes.insert(sex);
                }
                ParsedRecord::Interval {
                    reference,
                    start,
                    end,
                    sex,
                    ploidy,
                } => {
                    sexes.insert(sex.clone());
                    intervals
                        .entry(reference)
                        .or_default()
                        .entry(sex)
                        .or_default()
                        .push(PloidyInterval { start, end, ploidy });
                }
            }
        }
        let global_default = explicit_defaults
            .get("*")
            .copied()
            .unwrap_or(CallPloidy::Diploid);
        let defaults = sexes
            .iter()
            .map(|sex| {
                (
                    sex.clone(),
                    explicit_defaults
                        .get(sex)
                        .copied()
                        .unwrap_or(global_default),
                )
            })
            .collect::<BTreeMap<_, _>>();
        for (reference, by_sex) in &mut intervals {
            for (sex, values) in by_sex {
                values.sort_unstable_by_key(|interval| (interval.start, interval.end));
                if values.windows(2).any(|pair| pair[0].end >= pair[1].start) {
                    return Err(CallError::PloidyInput {
                        path: path.to_owned(),
                        message: format!("overlapping intervals for {reference} and sex {sex}"),
                    });
                }
            }
        }
        let default_sex = requested_default_sex.or_else(|| {
            if sexes.contains("F") {
                Some("F".to_owned())
            } else if sexes.len() == 1 {
                sexes.first().cloned()
            } else {
                None
            }
        });
        Ok(Self {
            defaults,
            intervals,
            default_sex,
        })
    }

    pub fn resolver(
        &self,
        assignments: impl IntoIterator<Item = SamplePloidy>,
    ) -> Result<PloidyResolver> {
        PloidyResolver::new(self.clone(), assignments.into_iter().collect())
    }

    pub fn default_resolver(&self, sample_count: usize) -> Result<PloidyResolver> {
        let assignment = self.default_assignment()?;
        self.resolver((0..sample_count).map(|_| assignment.clone()))
    }

    pub(crate) fn default_assignment(&self) -> Result<SamplePloidy> {
        self.default_sex
            .as_ref()
            .map(|sex| SamplePloidy::Sex(sex.clone()))
            .ok_or_else(|| CallError::UnknownPloidySex("default".to_owned()))
    }

    fn ploidy(&self, reference: &str, position: u64, sex: &str) -> Result<CallPloidy> {
        let default = self
            .defaults
            .get(sex)
            .copied()
            .ok_or_else(|| CallError::UnknownPloidySex(sex.to_owned()))?;
        let Some(intervals) = self
            .intervals
            .get(reference)
            .and_then(|by_sex| by_sex.get(sex))
        else {
            return Ok(default);
        };
        let index = intervals.partition_point(|interval| interval.start <= position);
        Ok(index
            .checked_sub(1)
            .and_then(|index| intervals.get(index))
            .filter(|interval| position <= interval.end)
            .map(|interval| interval.ploidy)
            .unwrap_or(default))
    }

    fn maximum(&self, sex: &str) -> Result<usize> {
        let default = self
            .defaults
            .get(sex)
            .copied()
            .ok_or_else(|| CallError::UnknownPloidySex(sex.to_owned()))?;
        Ok(self
            .intervals
            .values()
            .filter_map(|by_sex| by_sex.get(sex))
            .flat_map(|intervals| intervals.iter().map(|interval| interval.ploidy))
            .fold(default.chromosome_count(), |maximum, ploidy| {
                maximum.max(ploidy.chromosome_count())
            }))
    }
}

#[derive(Clone, Debug)]
pub struct PloidyResolver {
    definition: PloidyDefinition,
    assignments: Box<[SamplePloidy]>,
    prior_chromosome_count: usize,
}

impl PloidyResolver {
    fn new(definition: PloidyDefinition, assignments: Vec<SamplePloidy>) -> Result<Self> {
        let mut prior_chromosome_count = 0usize;
        for assignment in &assignments {
            prior_chromosome_count = prior_chromosome_count
                .checked_add(match assignment {
                    SamplePloidy::Sex(sex) => definition.maximum(sex)?,
                    SamplePloidy::Fixed(ploidy) => ploidy.chromosome_count(),
                })
                .ok_or(CallError::InvalidPriorChromosomeCount)?;
        }
        if assignments.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        if prior_chromosome_count == 0 {
            return Err(CallError::InvalidPriorChromosomeCount);
        }
        Ok(Self {
            definition,
            assignments: assignments.into(),
            prior_chromosome_count,
        })
    }

    pub fn sample_count(&self) -> usize {
        self.assignments.len()
    }

    pub fn prior_chromosome_count(&self) -> usize {
        self.prior_chromosome_count
    }

    pub fn resolve(&self, reference: &str, position: u64) -> Result<Vec<CallPloidy>> {
        let mut ploidies = Vec::with_capacity(self.assignments.len());
        self.resolve_into(reference, position, &mut ploidies)?;
        Ok(ploidies)
    }

    pub fn resolve_into(
        &self,
        reference: &str,
        position: u64,
        ploidies: &mut Vec<CallPloidy>,
    ) -> Result<()> {
        ploidies.clear();
        ploidies.reserve(self.assignments.len());
        for assignment in &self.assignments {
            ploidies.push(match assignment {
                SamplePloidy::Sex(sex) => self.definition.ploidy(reference, position, sex),
                SamplePloidy::Fixed(ploidy) => Ok(*ploidy),
            }?);
        }
        Ok(())
    }

    pub fn resolve_site(&self, reference: &str, site: &LikelihoodSite) -> Result<Vec<CallPloidy>> {
        let mut ploidies = Vec::with_capacity(self.assignments.len());
        self.resolve_site_into(reference, site, &mut ploidies)?;
        Ok(ploidies)
    }

    pub fn resolve_site_into(
        &self,
        reference: &str,
        site: &LikelihoodSite,
        ploidies: &mut Vec<CallPloidy>,
    ) -> Result<()> {
        if site.samples().len() != self.assignments.len() {
            return Err(CallError::PloidySampleCountMismatch);
        }
        self.resolve_into(reference, site.position(), ploidies)
    }
}

enum ParsedRecord {
    Default {
        sex: String,
        ploidy: CallPloidy,
    },
    Interval {
        reference: String,
        start: u64,
        end: u64,
        sex: String,
        ploidy: CallPloidy,
    },
}

fn parse_record(line: &str, path: &str, line_number: u64) -> Result<ParsedRecord> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let invalid = |message: String| CallError::PloidyRecord {
        path: path.to_owned(),
        line: line_number,
        message,
    };
    if fields.len() != 5 {
        return Err(invalid("expected CHROM FROM TO SEX PLOIDY".to_owned()));
    }
    let ploidy = parse_ploidy(fields[4])
        .ok_or_else(|| invalid(format!("ploidy must be 0, 1, or 2, found {}", fields[4])))?;
    if fields[0] == "*" {
        if fields[1] != "*" || fields[2] != "*" {
            return Err(invalid(
                "a default record must use '* * *' coordinates".to_owned(),
            ));
        }
        return Ok(ParsedRecord::Default {
            sex: fields[3].to_owned(),
            ploidy,
        });
    }
    let start = fields[1]
        .parse::<u64>()
        .ok()
        .filter(|value| *value != 0)
        .ok_or_else(|| invalid(format!("invalid one-based start: {}", fields[1])))?;
    let end = fields[2]
        .parse::<u64>()
        .ok()
        .filter(|value| *value >= start)
        .ok_or_else(|| invalid(format!("invalid inclusive end: {}", fields[2])))?;
    Ok(ParsedRecord::Interval {
        reference: fields[0].to_owned(),
        start: start - 1,
        end: end - 1,
        sex: fields[3].to_owned(),
        ploidy,
    })
}

fn parse_ploidy(value: &str) -> Option<CallPloidy> {
    match value {
        "0" => Some(CallPloidy::Absent),
        "1" => Some(CallPloidy::Haploid),
        "2" => Some(CallPloidy::Diploid),
        _ => None,
    }
}

fn grch37_records() -> Vec<ParsedRecord> {
    human_records(60_000, 2_699_521, 154_931_043, 59_373_566)
}

fn grch38_records() -> Vec<ParsedRecord> {
    human_records(9_999, 2_781_480, 155_701_381, 57_227_415)
}

fn human_records(
    first_x_end: u64,
    second_x_start: u64,
    second_x_end: u64,
    y_end: u64,
) -> Vec<ParsedRecord> {
    let mut records = Vec::new();
    for prefix in ["", "chr"] {
        records.extend([
            interval(prefix, "X", 1, first_x_end, "M", CallPloidy::Haploid),
            interval(
                prefix,
                "X",
                second_x_start,
                second_x_end,
                "M",
                CallPloidy::Haploid,
            ),
            interval(prefix, "Y", 1, y_end, "M", CallPloidy::Haploid),
            interval(prefix, "Y", 1, y_end, "F", CallPloidy::Absent),
        ]);
    }
    records.extend([
        interval("", "MT", 1, 16_569, "M", CallPloidy::Haploid),
        interval("", "MT", 1, 16_569, "F", CallPloidy::Haploid),
        interval("chr", "M", 1, 16_569, "M", CallPloidy::Haploid),
        interval("chr", "M", 1, 16_569, "F", CallPloidy::Haploid),
        ParsedRecord::Default {
            sex: "M".to_owned(),
            ploidy: CallPloidy::Diploid,
        },
        ParsedRecord::Default {
            sex: "F".to_owned(),
            ploidy: CallPloidy::Diploid,
        },
    ]);
    records
}

fn interval(
    prefix: &str,
    reference: &str,
    start: u64,
    end: u64,
    sex: &str,
    ploidy: CallPloidy,
) -> ParsedRecord {
    ParsedRecord::Interval {
        reference: format!("{prefix}{reference}"),
        start: start - 1,
        end: end - 1,
        sex: sex.to_owned(),
        ploidy,
    }
}

#[cfg(test)]
mod tests;
