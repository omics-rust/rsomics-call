mod call;
mod output;
mod pileup;

use std::process;

use clap::{Args, Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, ToolMeta, run as run_tool};

use self::output::{VariantOutputArgs, call_result};

const META: ToolMeta = ToolMeta {
    name: "rsomics-call",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-call",
    version,
    about = "Alignment likelihood and small-variant calling workflows",
    arg_required_else_help = true,
    subcommand_required = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Generate genotype likelihoods from coordinate-sorted alignments
    Pileup(pileup::Arguments),
    /// Call variants from a likelihood VCF or BCF
    Call(call::Arguments),
    /// Generate likelihoods and call variants without an intermediate file
    Run(RunArgs),
}

#[derive(Debug, Args)]
struct RunArgs {
    #[command(flatten)]
    input: pileup::AlignmentArgs,

    #[command(flatten)]
    likelihood: pileup::LikelihoodArgs,

    #[command(flatten)]
    call: call::CallPolicyArgs,

    #[command(flatten)]
    output: VariantOutputArgs,
}

impl RunArgs {
    fn execute(&self, json: bool) -> Result<()> {
        self.call.validate(false)?;
        let prepared = self.input.prepare()?;
        let input_paths = self.input.input_paths(&prepared);
        let run = call_result(self.input.open(prepared, &self.likelihood))?;
        let sample_names = run
            .samples()
            .samples()
            .iter()
            .map(|name| name.to_string())
            .collect::<Vec<_>>();
        let definition = self.call.ploidy_definition()?;
        let resolver = call_result(definition.default_resolver(sample_names.len()))?;
        let groups = self.call.sample_groups(&sample_names)?;
        let calls = call_result(self.call.build(resolver, groups))?;
        self.output.write(json, input_paths, |output| {
            run.run_calls(calls, output, self.output.output_type.format())
                .map(|_| ())
        })
    }
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<()> {
    let json = cli.output.json;
    match cli.command {
        Command::Pileup(arguments) => arguments.execute(json),
        Command::Call(arguments) => arguments.execute(json),
        Command::Run(arguments) => arguments.execute(json),
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        rsomics_help::command::<Cli>().debug_assert();
    }

    #[test]
    fn all_three_workflows_share_one_help_tree() {
        let top = Cli::command().render_long_help().to_string();
        for command in ["pileup", "call", "run"] {
            assert!(top.contains(command), "{top}");
        }
        let error = Cli::try_parse_from(["rsomics-call", "run", "--help"]).unwrap_err();
        let help = error.to_string();
        assert!(help.contains("--reference <FASTA>"), "{help}");
        assert!(help.contains("--ploidy <PRESET>"), "{help}");
        assert!(help.contains("--output-type <TYPE>"), "{help}");
    }
}
