use std::collections::{HashMap, HashSet};

use rsomics_bamio::raw::RawRecord;

use crate::{CallError, Result};

#[derive(Clone, Debug, Default)]
pub struct SampleSelection {
    mode: SelectionMode,
}

#[derive(Clone, Debug, Default)]
enum SelectionMode {
    #[default]
    All,
    Include(HashMap<Box<str>, Box<str>>),
    Exclude(HashSet<Box<str>>),
}

impl SampleSelection {
    pub fn include<N, O>(samples: impl IntoIterator<Item = (N, Option<O>)>) -> Result<Self>
    where
        N: Into<Box<str>>,
        O: Into<Box<str>>,
    {
        let mut selected = HashMap::new();
        for (name, output) in samples {
            let name = name.into();
            validate_sample_name(&name)?;
            let output = output.map_or_else(|| name.clone(), Into::into);
            validate_sample_name(&output)?;
            if selected.insert(name.clone(), output).is_some() {
                return Err(CallError::DuplicateSampleSelection(name.into()));
            }
        }
        if selected.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        Ok(Self {
            mode: SelectionMode::Include(selected),
        })
    }

    pub fn exclude(samples: impl IntoIterator<Item = impl Into<Box<str>>>) -> Result<Self> {
        let mut selected = HashSet::new();
        for name in samples {
            let name = name.into();
            validate_sample_name(&name)?;
            if !selected.insert(name.clone()) {
                return Err(CallError::DuplicateSampleSelection(name.into()));
            }
        }
        if selected.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        Ok(Self {
            mode: SelectionMode::Exclude(selected),
        })
    }
}

#[derive(Debug)]
pub struct SampleMapBuilder {
    selection: SampleSelection,
    seen: HashSet<Box<str>>,
    samples: Vec<Box<str>>,
    sample_indices: HashMap<Box<str>, usize>,
    sources: HashMap<u32, SourceSamples>,
}

impl SampleMapBuilder {
    pub fn new(selection: SampleSelection) -> Self {
        Self {
            selection,
            seen: HashSet::new(),
            samples: Vec::new(),
            sample_indices: HashMap::new(),
            sources: HashMap::new(),
        }
    }

    pub fn add_source<I, R, S>(
        &mut self,
        source_id: u32,
        name: impl Into<Box<str>>,
        ignore_read_groups: bool,
        read_groups: I,
    ) -> Result<()>
    where
        I: IntoIterator<Item = (R, S)>,
        R: Into<Box<[u8]>>,
        S: Into<Box<str>>,
    {
        if self.sources.contains_key(&source_id) {
            return Err(CallError::DuplicateAlignmentSource(source_id));
        }
        let name = name.into();
        validate_sample_name(&name)?;
        if ignore_read_groups {
            let default = self.select_sample(name);
            self.sources.insert(
                source_id,
                SourceSamples {
                    default,
                    fallback: default,
                    read_groups: HashMap::new(),
                },
            );
            return Ok(());
        }

        let mut groups = HashMap::new();
        let mut fallback = None;
        let mut group_count = 0;
        let mut all_same = true;
        let mut first_mapping = None;
        for (read_group, sample) in read_groups {
            group_count += 1;
            let read_group = read_group.into();
            validate_read_group(&read_group)?;
            let sample = sample.into();
            validate_sample_name(&sample)?;
            let mapping = self.select_sample(sample);
            if fallback.is_none() {
                fallback = mapping;
            }
            match first_mapping {
                None => first_mapping = Some(mapping),
                Some(first) if first != mapping => all_same = false,
                Some(_) => {}
            }
            if groups.insert(read_group.clone(), mapping).is_some() {
                return Err(CallError::InvalidReadGroup(
                    String::from_utf8_lossy(&read_group).into_owned(),
                ));
            }
        }

        if group_count == 0 {
            let default = self.select_sample(name);
            self.sources.insert(
                source_id,
                SourceSamples {
                    default,
                    fallback: default,
                    read_groups: groups,
                },
            );
        } else {
            let default = if all_same {
                first_mapping.flatten()
            } else {
                None
            };
            self.sources.insert(
                source_id,
                SourceSamples {
                    default,
                    fallback,
                    read_groups: groups,
                },
            );
        }
        Ok(())
    }

    pub fn finish(self) -> Result<SampleMap> {
        if let SelectionMode::Include(selected) = &self.selection.mode
            && let Some(missing) = selected
                .keys()
                .filter(|name| !self.seen.contains(*name))
                .min()
        {
            return Err(CallError::MissingSelectedSample(missing.to_string()));
        }
        if self.samples.is_empty() {
            return Err(CallError::InvalidSampleCount);
        }
        Ok(SampleMap {
            samples: self.samples.into_boxed_slice(),
            sources: self.sources,
        })
    }

    fn select_sample(&mut self, name: Box<str>) -> Option<usize> {
        self.seen.insert(name.clone());
        let output = match &self.selection.mode {
            SelectionMode::All => Some(name),
            SelectionMode::Include(selected) => selected.get(&name).cloned(),
            SelectionMode::Exclude(excluded) => (!excluded.contains(&name)).then_some(name),
        };
        output.map(|name| self.intern_sample(name))
    }

    fn intern_sample(&mut self, name: Box<str>) -> usize {
        if let Some(&index) = self.sample_indices.get(&name) {
            return index;
        }
        let index = self.samples.len();
        self.samples.push(name.clone());
        self.sample_indices.insert(name, index);
        index
    }
}

#[derive(Debug)]
pub struct SampleMap {
    samples: Box<[Box<str>]>,
    sources: HashMap<u32, SourceSamples>,
}

impl SampleMap {
    pub fn samples(&self) -> &[Box<str>] {
        &self.samples
    }

    pub fn sample_index(&self, source_id: u32, record: &RawRecord) -> Result<Option<usize>> {
        let source = if self.sources.len() == 1 {
            let (&expected, source) = self.sources.iter().next().unwrap();
            if source_id != expected {
                return Err(CallError::UnknownAlignmentSource(source_id));
            }
            source
        } else {
            self.sources
                .get(&source_id)
                .ok_or(CallError::UnknownAlignmentSource(source_id))?
        };
        if source.default.is_some() {
            return Ok(source.default);
        }
        let Some(raw) = record.aux_value(*b"RG") else {
            return Ok(source.fallback);
        };
        if record.aux_type(*b"RG") != Some(b'Z') {
            return Err(CallError::InvalidReadGroupField);
        }
        let read_group = raw
            .strip_suffix(&[0])
            .filter(|value| !value.is_empty())
            .ok_or(CallError::InvalidReadGroupField)?;
        Ok(source
            .read_groups
            .get(read_group)
            .copied()
            .unwrap_or(source.fallback))
    }
}

#[derive(Debug)]
struct SourceSamples {
    default: Option<usize>,
    fallback: Option<usize>,
    read_groups: HashMap<Box<[u8]>, Option<usize>>,
}

pub(crate) fn validate_sample_name(name: &str) -> Result<()> {
    if name.is_empty()
        || name
            .bytes()
            .any(|byte| matches!(byte, 0 | b'\t' | b'\n' | b'\r'))
    {
        Err(CallError::InvalidSampleName(name.to_owned()))
    } else {
        Ok(())
    }
}

fn validate_read_group(read_group: &[u8]) -> Result<()> {
    if read_group.is_empty()
        || matches!(read_group, b"*" | b"?")
        || read_group
            .iter()
            .any(|byte| matches!(byte, 0 | b'\t' | b'\n' | b'\r'))
    {
        Err(CallError::InvalidReadGroup(
            String::from_utf8_lossy(read_group).into_owned(),
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(read_group: Option<(&[u8], u8)>) -> RawRecord {
        let mut record = RawRecord::default();
        if let Some((read_group, kind)) = read_group {
            let mut value = read_group.to_vec();
            if kind == b'Z' {
                value.push(0);
            }
            record.append_aux(*b"RG", kind, &value).unwrap();
        }
        record
    }

    #[test]
    fn read_groups_merge_samples_across_sources() {
        let mut builder = SampleMapBuilder::new(SampleSelection::default());
        builder
            .add_source(
                4,
                "first.bam",
                false,
                [(b"rg1".as_slice(), "S1"), (b"rg2".as_slice(), "S2")],
            )
            .unwrap();
        builder
            .add_source(9, "second.bam", false, [(b"rg3".as_slice(), "S1")])
            .unwrap();
        let samples = builder.finish().unwrap();

        assert_eq!(
            samples
                .samples()
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            ["S1", "S2"]
        );
        assert_eq!(
            samples.sample_index(4, &record(Some((b"rg2", b'Z')))),
            Ok(Some(1))
        );
        assert_eq!(
            samples.sample_index(9, &record(Some((b"rg3", b'Z')))),
            Ok(Some(0))
        );
        assert_eq!(samples.sample_index(4, &record(None)), Ok(Some(0)));
    }

    #[test]
    fn ignored_read_groups_use_one_sample_per_input() {
        let mut builder = SampleMapBuilder::new(SampleSelection::default());
        builder
            .add_source(3, "input.bam", true, [(b"rg".as_slice(), "header-sample")])
            .unwrap();
        let samples = builder.finish().unwrap();

        assert_eq!(samples.samples()[0].as_ref(), "input.bam");
        assert_eq!(
            samples.sample_index(3, &record(Some((b"rg", b'Z')))),
            Ok(Some(0))
        );
    }

    #[test]
    fn inclusion_renames_and_excludes_read_groups() {
        let selection = SampleSelection::include([("S2", Some("renamed")), ("S3", None)]).unwrap();
        let mut builder = SampleMapBuilder::new(selection);
        builder
            .add_source(
                0,
                "input.bam",
                false,
                [
                    (b"one".as_slice(), "S1"),
                    (b"two".as_slice(), "S2"),
                    (b"three".as_slice(), "S3"),
                ],
            )
            .unwrap();
        let samples = builder.finish().unwrap();

        assert_eq!(
            samples
                .samples()
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            ["renamed", "S3"]
        );
        assert_eq!(
            samples.sample_index(0, &record(Some((b"one", b'Z')))),
            Ok(None)
        );
        assert_eq!(
            samples.sample_index(0, &record(Some((b"two", b'Z')))),
            Ok(Some(0))
        );
        assert_eq!(
            samples.sample_index(0, &record(Some((b"unknown", b'Z')))),
            Ok(Some(0))
        );
    }

    #[test]
    fn exclusion_and_duplicate_boundaries_are_checked() {
        let selection = SampleSelection::exclude(["S2"]).unwrap();
        let mut builder = SampleMapBuilder::new(selection);
        builder
            .add_source(
                0,
                "input.bam",
                false,
                [(b"one".as_slice(), "S1"), (b"two".as_slice(), "S2")],
            )
            .unwrap();
        assert_eq!(
            builder.add_source(0, "duplicate.bam", false, Vec::<(&[u8], &str)>::new(),),
            Err(CallError::DuplicateAlignmentSource(0))
        );
        let samples = builder.finish().unwrap();
        assert_eq!(samples.samples()[0].as_ref(), "S1");
        assert_eq!(
            samples.sample_index(0, &record(Some((b"two", b'Z')))),
            Ok(None)
        );

        let mut duplicate = SampleMapBuilder::new(SampleSelection::default());
        assert!(matches!(
            duplicate.add_source(
                1,
                "input.bam",
                false,
                [(b"rg".as_slice(), "S1"), (b"rg".as_slice(), "S2")],
            ),
            Err(CallError::InvalidReadGroup(_))
        ));
    }

    #[test]
    fn invalid_source_selection_and_rg_fields_fail() {
        let selection = SampleSelection::include([("absent", None::<&str>)]).unwrap();
        let mut missing = SampleMapBuilder::new(selection);
        missing
            .add_source(0, "input.bam", false, [(b"rg".as_slice(), "S1")])
            .unwrap();
        assert_eq!(
            missing.finish().unwrap_err(),
            CallError::MissingSelectedSample("absent".to_owned())
        );

        let mut builder = SampleMapBuilder::new(SampleSelection::default());
        builder
            .add_source(
                0,
                "input.bam",
                false,
                [(b"rg".as_slice(), "S1"), (b"rg2".as_slice(), "S2")],
            )
            .unwrap();
        let samples = builder.finish().unwrap();
        assert_eq!(
            samples.sample_index(1, &record(None)),
            Err(CallError::UnknownAlignmentSource(1))
        );
        assert_eq!(
            samples.sample_index(0, &record(Some((b"x", b'A')))),
            Err(CallError::InvalidReadGroupField)
        );
    }
}
