use std::path::Path;

use noodles::core::{Position, Region};
use rsomics_intervals::{GenomicRegion, RegionFileError, RegionFileMode, read_region_file};

use crate::{CallError, Result};

#[derive(Clone, Copy)]
enum Purpose {
    Regions,
    Targets,
}

pub(crate) fn read_regions(path: &Path) -> Result<Vec<Region>> {
    read(path, Purpose::Regions)
}

pub(crate) fn read_targets(path: &Path) -> Result<Vec<Region>> {
    read(path, Purpose::Targets)
}

fn read(path: &Path, purpose: Purpose) -> Result<Vec<Region>> {
    let mode = match purpose {
        Purpose::Regions => RegionFileMode::Intervals,
        Purpose::Targets => RegionFileMode::Targets,
    };
    read_region_file(path, mode)
        .map_err(|error| map_error(path, error, purpose))?
        .into_iter()
        .map(|region| convert_region(path, region, purpose))
        .collect()
}

fn convert_region(path: &Path, region: GenomicRegion, purpose: Purpose) -> Result<Region> {
    let name = region.chrom().as_bytes().to_vec();
    if region.is_whole_reference() {
        return Ok(Region::new(name, ..));
    }
    let start = region.start().expect("bounded region has a start");
    let end = region.end().expect("bounded region has an end");
    let start = usize::try_from(start)
        .ok()
        .and_then(|position| position.checked_add(1))
        .and_then(Position::new)
        .ok_or_else(|| range_error(path, purpose, "interval start"))?;
    let end = usize::try_from(end)
        .ok()
        .and_then(Position::new)
        .ok_or_else(|| range_error(path, purpose, "interval end"))?;
    Ok(Region::new(name, start..=end))
}

fn range_error(path: &Path, purpose: Purpose, field: &str) -> CallError {
    input_error(
        path,
        purpose,
        format!("{field} exceeds the supported coordinate range"),
    )
}

fn map_error(input_path: &Path, error: RegionFileError, purpose: Purpose) -> CallError {
    match error {
        RegionFileError::Input { path, source } => input_error(&path, purpose, source.to_string()),
        RegionFileError::Record {
            path,
            line,
            message,
        } => record_error(&path, line, purpose, message),
        error => input_error(input_path, purpose, error.to_string()),
    }
}

fn input_error(path: &Path, purpose: Purpose, message: String) -> CallError {
    match purpose {
        Purpose::Regions => CallError::RegionInput {
            path: path.display().to_string(),
            message,
        },
        Purpose::Targets => CallError::TargetInput {
            path: path.display().to_string(),
            message,
        },
    }
}

fn record_error(path: &Path, line: u64, purpose: Purpose, message: String) -> CallError {
    match purpose {
        Purpose::Regions => CallError::RegionRecord {
            path: path.display().to_string(),
            line,
            message,
        },
        Purpose::Targets => CallError::TargetRecord {
            path: path.display().to_string(),
            line,
            message,
        },
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;

    use flate2::Compression;
    use flate2::write::GzEncoder;

    use super::*;

    #[test]
    fn reads_bed_tab_and_vcf_coordinates() {
        let directory = tempfile::tempdir().unwrap();
        let bed = directory.path().join("targets.BED");
        fs::write(&bed, b"# comment\nchr1\t0\t2\tname\nchr2\n").unwrap();
        assert_eq!(
            read_targets(&bed)
                .unwrap()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            ["chr1:1-2", "chr2"]
        );

        let tab = directory.path().join("targets.txt");
        fs::write(&tab, b"chr1\t2\nchr1\t4\t5\n").unwrap();
        assert_eq!(
            read_targets(&tab)
                .unwrap()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            ["chr1:2-2", "chr1:4-5"]
        );

        let vcf = directory.path().join("targets.vcf");
        fs::write(
            &vcf,
            b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t7\t.\tA\tC\n",
        )
        .unwrap();
        assert_eq!(read_targets(&vcf).unwrap()[0].to_string(), "chr1:7-7");
    }

    #[test]
    fn detects_gzip_content_and_bed_suffix() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("targets.bed.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"chr1\t1\t3\n").unwrap();
        fs::write(&path, encoder.finish().unwrap()).unwrap();

        assert_eq!(read_targets(&path).unwrap()[0].to_string(), "chr1:2-3");
        assert_eq!(read_regions(&path).unwrap()[0].to_string(), "chr1:2-3");
    }

    #[test]
    fn reports_line_context_and_compression_failures() {
        let directory = tempfile::tempdir().unwrap();
        let invalid = directory.path().join("targets.bed");
        fs::write(&invalid, b"chr1\t4\t2\n").unwrap();
        assert!(matches!(
            read_targets(&invalid),
            Err(CallError::TargetRecord { line: 1, .. })
        ));

        let truncated = directory.path().join("targets.txt.gz");
        fs::write(&truncated, [0x1f, 0x8b, 0x08]).unwrap();
        assert!(matches!(
            read_targets(&truncated),
            Err(CallError::TargetInput { .. })
        ));
    }

    #[test]
    fn region_files_require_a_consistent_coordinate_shape() {
        let directory = tempfile::tempdir().unwrap();
        let mixed = directory.path().join("regions.txt");
        fs::write(&mixed, b"chr1\t2\nchr1\t4\t5\n").unwrap();
        assert!(matches!(
            read_regions(&mixed),
            Err(CallError::RegionRecord { line: 2, .. })
        ));

        let bed = directory.path().join("regions.bed");
        fs::write(&bed, b"chr1\n").unwrap();
        assert!(matches!(
            read_regions(&bed),
            Err(CallError::RegionRecord { line: 1, .. })
        ));
    }
}
