# AssemblyBrute documentation

AssemblyBrute creates multiple de novo genome assemblies using different tools
and parameter sets, then evaluates them with shared quality-control analyses.

## Current documentation

Start with setup and a test-data example:

1. [Installation](installation.md): controller and named tool environments.
2. [Preparing inputs](inputs.md): directory layout, formats, and naming rules.
3. [Test-data examples](examples/README.md): datasets, complete configurations,
   commands, expected outputs, and validation status.

Configure and operate a project:

- [Configuration](configuration.md): file ownership, precedence, presets,
  resource tables, and analysis settings.
- [Stages](stages.md): ordering, prerequisites, branch labels, and QC levels.
- [Running](running.md): local commands, explicit targets, Slurm settings,
  and resource overrides.

Interpret results and choose a workflow:

- [Outputs](outputs.md): assembly paths, reports, and summary tables.
- [Recipes](recipes.md): read QC, ploidy models, primary assembly, preset
  comparison, draft QC, and Hi-C scaffolding.
- [Troubleshooting](troubleshooting.md): pipeline-specific errors and limitations.

[Developer records](development/README.md) contain the implementation audit and
validation history. Previous manuals are preserved in [old/](old/README.md).
