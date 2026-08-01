use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::{
    CallError, CallPloidy, IndexedLikelihoodVariantReader, LikelihoodVariantReader,
    LikelihoodVcfSchema, PloidyDefinition, PloidyResolver, Result, SamplePloidy,
    samples::validate_sample_name,
};

#[derive(Clone, Debug, Default)]
pub struct CallSampleSelection {
    mode: SelectionMode,
}

#[derive(Clone, Debug, Default)]
enum SelectionMode {
    #[default]
    All,
    Include(Box<[SampleAssignment]>),
    Exclude(Box<[Box<str>]>),
}

#[derive(Clone, Debug)]
struct SampleAssignment {
    name: Box<str>,
    ploidy: Option<SamplePloidy>,
}

impl CallSampleSelection {
    pub fn include(samples: impl IntoIterator<Item = impl Into<Box<str>>>) -> Result<Self> {
        Self::include_with_ploidy(samples.into_iter().map(|name| (name, None)))
    }

    pub fn include_with_ploidy<N>(
        samples: impl IntoIterator<Item = (N, Option<SamplePloidy>)>,
    ) -> Result<Self>
    where
        N: Into<Box<str>>,
    {
        let mut seen = HashSet::new();
        let mut assignments = Vec::new();
        for (name, ploidy) in samples {
            let name = name.into();
            validate_sample_name(&name)?;
            if !seen.insert(name.clone()) {
                return Err(CallError::DuplicateSampleSelection(name.into()));
            }
            assignments.push(SampleAssignment { name, ploidy });
        }
        if assignments.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        Ok(Self {
            mode: SelectionMode::Include(assignments.into_boxed_slice()),
        })
    }

    pub fn exclude(samples: impl IntoIterator<Item = impl Into<Box<str>>>) -> Result<Self> {
        let mut seen = HashSet::new();
        let mut names = Vec::new();
        for name in samples {
            let name = name.into();
            validate_sample_name(&name)?;
            if !seen.insert(name.clone()) {
                return Err(CallError::DuplicateSampleSelection(name.into()));
            }
            names.push(name);
        }
        if names.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        Ok(Self {
            mode: SelectionMode::Exclude(names.into_boxed_slice()),
        })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        Self::read_mode(path.as_ref(), false)
    }

    pub fn read_excluding(path: impl AsRef<Path>) -> Result<Self> {
        Self::read_mode(path.as_ref(), true)
    }

    pub fn parse(data: &str) -> Result<Self> {
        Self::from_reader(data.as_bytes(), "<samples>".to_owned(), false)
    }

    pub fn parse_excluding(data: &str) -> Result<Self> {
        Self::from_reader(data.as_bytes(), "<samples>".to_owned(), true)
    }

    fn read_mode(path: &Path, exclude: bool) -> Result<Self> {
        let display = path.display().to_string();
        let file = File::open(path).map_err(|error| CallError::CallSampleInput {
            path: display.clone(),
            message: error.to_string(),
        })?;
        Self::from_reader(BufReader::new(file), display, exclude)
    }

    pub fn bind<R>(
        self,
        reader: LikelihoodVariantReader<R>,
        definition: &PloidyDefinition,
    ) -> Result<(LikelihoodVariantReader<R>, PloidyResolver)>
    where
        R: Read,
    {
        self.bind_reader(reader, definition)
    }

    pub fn bind_indexed(
        self,
        reader: IndexedLikelihoodVariantReader,
        definition: &PloidyDefinition,
    ) -> Result<(IndexedLikelihoodVariantReader, PloidyResolver)> {
        self.bind_reader(reader, definition)
    }

    fn bind_reader<R>(self, reader: R, definition: &PloidyDefinition) -> Result<(R, PloidyResolver)>
    where
        R: ProjectLikelihoodSamples,
    {
        match self.mode {
            SelectionMode::All => {
                let resolver =
                    definition.default_resolver(reader.schema().header().sample_names().len())?;
                Ok((reader, resolver))
            }
            SelectionMode::Include(assignments) => {
                let reader = reader.include(&assignments)?;
                let mut default: Option<SamplePloidy> = None;
                let mut ploidies = Vec::with_capacity(assignments.len());
                for assignment in assignments {
                    let ploidy = match assignment.ploidy {
                        Some(ploidy) => ploidy,
                        None => match &default {
                            Some(ploidy) => ploidy.clone(),
                            None => {
                                let ploidy = definition.default_assignment()?;
                                default = Some(ploidy.clone());
                                ploidy
                            }
                        },
                    };
                    ploidies.push(ploidy);
                }
                let resolver = definition.resolver(ploidies)?;
                Ok((reader, resolver))
            }
            SelectionMode::Exclude(names) => {
                let reader = reader.exclude(&names)?;
                let resolver =
                    definition.default_resolver(reader.schema().header().sample_names().len())?;
                Ok((reader, resolver))
            }
        }
    }

    fn from_reader(reader: impl BufRead, path: String, exclude: bool) -> Result<Self> {
        let mut records = Vec::new();
        for (index, line) in reader.lines().enumerate() {
            let line_number = index as u64 + 1;
            let line = line.map_err(|error| CallError::CallSampleRecord {
                path: path.clone(),
                line: line_number,
                message: error.to_string(),
            })?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() > 2 {
                return Err(CallError::CallSampleRecord {
                    path,
                    line: line_number,
                    message: "expected SAMPLE and optional PLOIDY_OR_SEX".to_owned(),
                });
            }
            records.push((
                Box::<str>::from(fields[0]),
                fields.get(1).map(|value| parse_assignment(value)),
            ));
        }
        if exclude {
            Self::exclude(records.into_iter().map(|(name, _)| name))
        } else {
            Self::include_with_ploidy(records)
        }
    }
}

trait ProjectLikelihoodSamples: Sized {
    fn schema(&self) -> &LikelihoodVcfSchema;
    fn include(self, assignments: &[SampleAssignment]) -> Result<Self>;
    fn exclude(self, names: &[Box<str>]) -> Result<Self>;
}

impl<R> ProjectLikelihoodSamples for LikelihoodVariantReader<R>
where
    R: Read,
{
    fn schema(&self) -> &LikelihoodVcfSchema {
        self.schema()
    }

    fn include(self, assignments: &[SampleAssignment]) -> Result<Self> {
        self.select_samples(assignments.iter().map(|item| &item.name))
    }

    fn exclude(self, names: &[Box<str>]) -> Result<Self> {
        self.exclude_samples(names)
    }
}

impl ProjectLikelihoodSamples for IndexedLikelihoodVariantReader {
    fn schema(&self) -> &LikelihoodVcfSchema {
        self.schema()
    }

    fn include(self, assignments: &[SampleAssignment]) -> Result<Self> {
        self.select_samples(assignments.iter().map(|item| &item.name))
    }

    fn exclude(self, names: &[Box<str>]) -> Result<Self> {
        self.exclude_samples(names)
    }
}

fn parse_assignment(value: &str) -> SamplePloidy {
    match value {
        "0" => SamplePloidy::Fixed(CallPloidy::Absent),
        "1" => SamplePloidy::Fixed(CallPloidy::Haploid),
        "2" => SamplePloidy::Fixed(CallPloidy::Diploid),
        _ => SamplePloidy::sex(value),
    }
}

#[cfg(test)]
mod tests;
