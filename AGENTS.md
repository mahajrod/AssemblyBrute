# Repository Guidelines

## Project Structure & Module Organization
AssemblyBrute builds and evaluates de novo genome assemblies using Snakemake. `Snakefile` orchestrates the workflow; `workflow/rules/` groups rules by stage and tool. Shared Python helpers live in `workflow/include/`, standalone Python/Bash/AWK utilities in `workflow/scripts/`, and Conda definitions in `workflow/envs/`. Treat `workflow/external_tools/` as bundled upstream software.

Standalone examples live in `config/templates/<dataset>/`; `docs/examples/README.md` indexes datasets and test commands. Shared defaults and resource tables are in `config/core.yaml` and `config/resources/`; scheduler profiles are in `profile/`. Test inputs use `test_data/<dataset>/<datatype>/<format>/`, with outputs in `results/<dataset>/`. Only `DESCRIPTION` markers under `test_data/` are trackable; data are ignored. Documentation and illustrations live in `docs/` and `images/`.

## Build, Test, and Development Commands
Use Linux and Snakemake 7.32.4. Templates select existing named Conda environments with `use_existing_envs: true`; `false` selects YAML specifications relative to the main `Snakefile`. Run from the repository root; there is no separate build step. Fresh Conda solves remain unverified.

Choose a complete configuration from the example catalog and save a copy as `config/data/local.yaml`; adjust sample metadata, inputs, outputs, stages, and resources.

- `snakemake --configfile config/data/local.yaml --use-conda -c 8 --dry-run --printshellcmds` — inspect the project DAG.
- Remove `--dry-run` to execute after preparing inputs and checking the graph.
- `snakemake --profile profile/slurm/ --configfile config/data/local.yaml --use-conda` — submit jobs after adapting cluster settings.

## Coding Style & Naming Conventions
Use four-space indentation for Python/Snakemake and descriptive `snake_case` names. Reuse shared helpers and preserve configuration keys, wildcards, and output naming. No repository-wide formatter or linter is configured. Prefer configuration/documentation changes; keep workflow edits minimal and explain their effects.

## Testing Guidelines
Default to DAG checks; run scientific jobs or download sequencing data only when explicitly requested. No `tests/` directory, first-party test framework, or coverage threshold is configured. Check scripts with small fixtures. Synthetic inputs establish graph construction, not real-data execution. Dry runs write initialization files; use an isolated output root. `--cores` does not limit RAM.

Set `input_dir`, `out_dir`, and `tmp_dir` explicitly per dataset. These project/CLI paths override core defaults; other colliding keys can still be replaced by core. Changing `out_dir` does not derive `tmp_dir`.

## Commit & Pull Request Guidelines
Use concise subjects naming the affected tool or stage; history has no enforced prefix convention. PRs should describe behavior/configuration changes, validation and its limits, and dependency/output changes. Link relevant issues and include images for plot changes. Exclude sequencing data, results, credentials, and machine-specific paths.
