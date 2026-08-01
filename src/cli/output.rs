use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError, write_output};

use crate::{CallError, VariantOutputFormat};

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub(crate) enum OutputType {
    #[default]
    V,
    Z,
    U,
    B,
}

impl OutputType {
    pub(crate) fn format(self) -> VariantOutputFormat {
        match self {
            Self::V => VariantOutputFormat::Vcf,
            Self::Z => VariantOutputFormat::VcfBgzf,
            Self::U => VariantOutputFormat::BcfRaw,
            Self::B => VariantOutputFormat::BcfBgzf,
        }
    }
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Output")]
pub(crate) struct VariantOutputArgs {
    /// Output VCF or BCF; omit or use - for standard output
    #[arg(
        short,
        long,
        value_name = "VARIANT",
        default_value = "-",
        hide_default_value = true
    )]
    pub(crate) output: PathBuf,

    /// Output type: v VCF, z compressed VCF, u raw BCF, b compressed BCF
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "TYPE",
        default_value = "v"
    )]
    pub(crate) output_type: OutputType,
}

impl VariantOutputArgs {
    pub(crate) fn write<T>(
        &self,
        json: bool,
        inputs: impl IntoIterator<Item = PathBuf>,
        operation: impl FnOnce(&mut dyn Write) -> crate::Result<T>,
    ) -> Result<T> {
        if json && self.output == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "--json requires a named --output so JSON cannot mix with variant stdout"
                    .to_owned(),
            ));
        }
        reject_output_alias(&self.output, inputs)?;
        write_output(Some(&self.output), |output| {
            operation(output).map_err(map_call_error)
        })
    }
}

pub(crate) fn map_call_error(error: CallError) -> RsomicsError {
    match error {
        CallError::VariantOutput(message) => {
            RsomicsError::Io(std::io::Error::other(format!("variant output: {message}")))
        }
        error => RsomicsError::InvalidInput(error.to_string()),
    }
}

pub(crate) fn call_result<T>(result: crate::Result<T>) -> Result<T> {
    result.map_err(map_call_error)
}

fn reject_output_alias(output: &Path, inputs: impl IntoIterator<Item = PathBuf>) -> Result<()> {
    if output == Path::new("-") {
        return Ok(());
    }
    if inputs
        .into_iter()
        .filter(|input| input != Path::new("-"))
        .any(|input| paths_alias(&input, output))
    {
        return Err(RsomicsError::ConfigError(format!(
            "output {} is also an input path",
            output.display()
        )));
    }
    Ok(())
}

fn paths_alias(left: &Path, right: &Path) -> bool {
    left == right
        || same_file::is_same_file(left, right).unwrap_or(false)
        || matches!(
            (fs::canonicalize(left), fs::canonicalize(right)),
            (Ok(left), Ok(right)) if left == right
        )
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;

    #[test]
    fn named_output_cannot_alias_an_input() {
        let directory = tempdir().unwrap();
        let input = directory.path().join("input.bcf");
        fs::write(&input, b"input").unwrap();
        let output = VariantOutputArgs {
            output: input.clone(),
            output_type: OutputType::B,
        };
        let error = output.write(false, [input], |_| Ok(())).unwrap_err();
        assert!(error.to_string().contains("also an input"), "{error}");
    }

    #[test]
    fn json_requires_a_named_variant_output() {
        let output = VariantOutputArgs {
            output: PathBuf::from("-"),
            output_type: OutputType::V,
        };
        let error = output.write(true, [], |_| Ok(())).unwrap_err();
        assert!(error.to_string().contains("requires a named"), "{error}");
    }
}
