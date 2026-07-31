use std::ops::Range;

use noodles::core::{Position, Region};

use crate::{CallError, LikelihoodSite, ReferenceSequence, Result};

pub(crate) struct RegionSelection {
    pub(crate) query: Region,
    pub(crate) bounds: ReferenceRange,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ReferenceRange {
    reference_id: usize,
    start: u64,
    end: u64,
}

impl ReferenceRange {
    pub(crate) fn contains(self, reference_id: usize, position: u64) -> bool {
        reference_id == self.reference_id && position >= self.start && position < self.end
    }

    pub(crate) fn contains_site(self, site: &LikelihoodSite) -> bool {
        self.contains(site.reference_sequence_id(), site.position())
    }
}

pub(crate) struct TargetSet {
    intervals: Box<[Box<[Range<u64>]>]>,
}

impl TargetSet {
    pub(crate) fn from_regions(
        references: &[ReferenceSequence],
        regions: impl IntoIterator<Item = Region>,
    ) -> Self {
        let mut intervals = vec![Vec::new(); references.len()];
        for region in regions {
            let Some(reference_id) = references
                .iter()
                .position(|reference| reference.name() == region.name())
            else {
                continue;
            };
            if let Some(bounds) =
                clipped_bounds(reference_id, references[reference_id].length(), &region)
            {
                intervals[reference_id].push(bounds.start..bounds.end);
            }
        }

        Self {
            intervals: intervals
                .into_iter()
                .map(merge_intervals)
                .collect::<Box<[_]>>(),
        }
    }

    pub(crate) fn contains(&self, reference_id: usize, position: u64) -> bool {
        let Some(intervals) = self.intervals.get(reference_id) else {
            return false;
        };
        let index = intervals.partition_point(|interval| interval.start <= position);
        index != 0 && position < intervals[index - 1].end
    }

    pub(crate) fn contains_site(&self, site: &LikelihoodSite) -> bool {
        self.contains(site.reference_sequence_id(), site.position())
    }
}

pub(crate) fn normalize_regions(
    references: &[ReferenceSequence],
    regions: impl IntoIterator<Item = Region>,
) -> Result<Box<[RegionSelection]>> {
    let mut bounds = regions
        .into_iter()
        .map(|region| resolve_region(references, &region))
        .collect::<Result<Vec<_>>>()?;
    if bounds.is_empty() {
        return Err(CallError::MissingRegions);
    }
    bounds.sort_unstable_by_key(|region| (region.reference_id, region.start, region.end));

    let mut merged: Vec<ReferenceRange> = Vec::with_capacity(bounds.len());
    for region in bounds {
        if let Some(previous) = merged.last_mut()
            && previous.reference_id == region.reference_id
            && region.start <= previous.end
        {
            previous.end = previous.end.max(region.end);
        } else {
            merged.push(region);
        }
    }

    merged
        .into_iter()
        .map(|bounds| {
            let name = references[bounds.reference_id].name();
            let invalid_coordinate = || CallError::InvalidRegion {
                region: String::from_utf8_lossy(name).into_owned(),
                message: "coordinate exceeds the supported index range".to_owned(),
            };
            let start = usize::try_from(bounds.start)
                .ok()
                .and_then(|position| position.checked_add(1))
                .and_then(Position::new)
                .ok_or_else(&invalid_coordinate)?;
            let end = usize::try_from(bounds.end)
                .ok()
                .and_then(Position::new)
                .ok_or_else(invalid_coordinate)?;
            Ok(RegionSelection {
                query: Region::new(name.to_vec(), start..=end),
                bounds,
            })
        })
        .collect::<Result<Box<[_]>>>()
}

fn resolve_region(references: &[ReferenceSequence], region: &Region) -> Result<ReferenceRange> {
    let reference_id = references
        .iter()
        .position(|reference| reference.name() == region.name())
        .ok_or_else(|| CallError::InvalidRegion {
            region: region.to_string(),
            message: "reference sequence is absent from the alignment header".to_owned(),
        })?;
    clipped_bounds(reference_id, references[reference_id].length(), region).ok_or_else(|| {
        CallError::InvalidRegion {
            region: region.to_string(),
            message: "interval is outside the reference sequence".to_owned(),
        }
    })
}

fn clipped_bounds(
    reference_id: usize,
    reference_length: u64,
    region: &Region,
) -> Option<ReferenceRange> {
    let interval = region.interval();
    let start = interval
        .start()
        .map(|position| u64::try_from(usize::from(position) - 1).unwrap())
        .unwrap_or(0);
    let end = interval
        .end()
        .map(|position| u64::try_from(usize::from(position)).unwrap())
        .unwrap_or(reference_length)
        .min(reference_length);
    (start < end).then_some(ReferenceRange {
        reference_id,
        start,
        end,
    })
}

fn merge_intervals(mut intervals: Vec<Range<u64>>) -> Box<[Range<u64>]> {
    intervals.sort_unstable_by_key(|interval| (interval.start, interval.end));
    let mut merged: Vec<Range<u64>> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        if let Some(previous) = merged.last_mut()
            && interval.start <= previous.end
        {
            previous.end = previous.end.max(interval.end);
        } else {
            merged.push(interval);
        }
    }
    merged.into_boxed_slice()
}
