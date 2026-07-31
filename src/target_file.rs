use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use noodles::core::{Position, Region};

use crate::{CallError, Result};

enum Format {
    Bed,
    Tab,
    Vcf,
}

pub(crate) fn read(path: &Path) -> Result<Vec<Region>> {
    let format = format(path);
    let mut reader = open(path)?;
    let mut line = String::new();
    let mut line_number = 0;
    let mut regions = Vec::new();

    loop {
        line.clear();
        let length = reader
            .read_line(&mut line)
            .map_err(|error| input_error(path, error))?;
        if length == 0 {
            break;
        }
        line_number += 1;
        let record = line.trim();
        if record.is_empty() || record.starts_with('#') {
            continue;
        }
        let region = match format {
            Format::Bed => parse_bed(record),
            Format::Tab => parse_tab(record),
            Format::Vcf => parse_vcf(record),
        }
        .map_err(|message| CallError::TargetRecord {
            path: path.display().to_string(),
            line: line_number,
            message,
        })?;
        regions.push(region);
    }

    Ok(regions)
}

fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|error| input_error(path, error))?;
    let mut source = BufReader::new(file);
    let gzip = source
        .fill_buf()
        .map_err(|error| input_error(path, error))?
        .starts_with(&[0x1f, 0x8b]);
    if gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(source))))
    } else {
        Ok(Box::new(source))
    }
}

fn format(path: &Path) -> Format {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if name.ends_with(".bed") || name.ends_with(".bed.gz") || name.ends_with(".bed.bgz") {
        Format::Bed
    } else if name.ends_with(".vcf") || name.ends_with(".vcf.gz") {
        Format::Vcf
    } else {
        Format::Tab
    }
}

fn parse_bed(record: &str) -> std::result::Result<Region, String> {
    let fields = record.split_whitespace().collect::<Vec<_>>();
    match fields.as_slice() {
        [name] => Ok(Region::new(name.as_bytes().to_vec(), ..)),
        [name, start, end, ..] => {
            let start = coordinate(start, "BED start")?;
            let end = coordinate(end, "BED end")?;
            half_open_region(name, start, end)
        }
        _ => Err("expected CHROM or at least three BED columns".to_owned()),
    }
}

fn parse_tab(record: &str) -> std::result::Result<Region, String> {
    let fields = record.split_whitespace().collect::<Vec<_>>();
    match fields.as_slice() {
        [name] => Ok(Region::new(name.as_bytes().to_vec(), ..)),
        [name, position] => {
            let position = one_based_coordinate(position, "position")?;
            half_open_region(name, position - 1, position)
        }
        [name, start, end, ..] => {
            let start = one_based_coordinate(start, "interval start")?;
            let end = one_based_coordinate(end, "interval end")?;
            half_open_region(name, start - 1, end)
        }
        _ => Err("expected CHROM, CHROM POS, or CHROM START END".to_owned()),
    }
}

fn parse_vcf(record: &str) -> std::result::Result<Region, String> {
    let fields = record.split_whitespace().collect::<Vec<_>>();
    let [name, position, ..] = fields.as_slice() else {
        return Err("expected at least CHROM and POS VCF columns".to_owned());
    };
    let position = one_based_coordinate(position, "VCF position")?;
    half_open_region(name, position - 1, position)
}

fn half_open_region(name: &str, start: u64, end: u64) -> std::result::Result<Region, String> {
    if start >= end {
        return Err(format!(
            "interval must contain at least one coordinate: start={start}, end={end}"
        ));
    }
    let start = usize::try_from(start)
        .ok()
        .and_then(|position| position.checked_add(1))
        .and_then(Position::new)
        .ok_or_else(|| "interval start exceeds the supported coordinate range".to_owned())?;
    let end = usize::try_from(end)
        .ok()
        .and_then(Position::new)
        .ok_or_else(|| "interval end exceeds the supported coordinate range".to_owned())?;
    Ok(Region::new(name.as_bytes().to_vec(), start..=end))
}

fn coordinate(value: &str, field: &str) -> std::result::Result<u64, String> {
    value
        .parse()
        .map_err(|error| format!("invalid {field} {value}: {error}"))
}

fn one_based_coordinate(value: &str, field: &str) -> std::result::Result<u64, String> {
    let coordinate = coordinate(value, field)?;
    if coordinate == 0 {
        Err(format!("{field} must be one-based"))
    } else {
        Ok(coordinate)
    }
}

fn input_error(path: &Path, error: io::Error) -> CallError {
    CallError::TargetInput {
        path: path.display().to_string(),
        message: error.to_string(),
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
            read(&bed)
                .unwrap()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            ["chr1:1-2", "chr2"]
        );

        let tab = directory.path().join("targets.txt");
        fs::write(&tab, b"chr1\t2\nchr1\t4\t5\n").unwrap();
        assert_eq!(
            read(&tab)
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
        assert_eq!(read(&vcf).unwrap()[0].to_string(), "chr1:7-7");
    }

    #[test]
    fn detects_gzip_content_and_bed_suffix() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("targets.bed.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"chr1\t1\t3\n").unwrap();
        fs::write(&path, encoder.finish().unwrap()).unwrap();

        assert_eq!(read(&path).unwrap()[0].to_string(), "chr1:2-3");
    }

    #[test]
    fn reports_line_context_and_compression_failures() {
        let directory = tempfile::tempdir().unwrap();
        let invalid = directory.path().join("targets.bed");
        fs::write(&invalid, b"chr1\t4\t2\n").unwrap();
        assert!(matches!(
            read(&invalid),
            Err(CallError::TargetRecord { line: 1, .. })
        ));

        let truncated = directory.path().join("targets.txt.gz");
        fs::write(&truncated, [0x1f, 0x8b, 0x08]).unwrap();
        assert!(matches!(
            read(&truncated),
            Err(CallError::TargetInput { .. })
        ));
    }
}
