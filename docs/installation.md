# Installation

Use Linux, Conda or Mamba, and Snakemake **7.32.4**. Run commands from the
repository root.

## Get the repository

```bash
git clone https://github.com/mahajrod/AssemblyBrute.git
cd AssemblyBrute
```

## Create the workflow controller

```bash
conda create -n assemblybrute -c conda-forge -c bioconda python=3.11 snakemake=7.32.4
conda activate assemblybrute
snakemake --version
```

The controller also needs pandas, NumPy, and PyYAML for workflow initialization.
System dependencies include Bash, `bc`, `pigz`, `gzip`, coreutils, `curl`, and Git.

## Create tool environments

Tool environment specifications are in `workflow/envs/conda.*.yaml`. To create
all supplied environments, run the following from the repository root:

```bash
for env_file in workflow/envs/conda.*.yaml; do
  conda env create --file "$env_file" || break
done
```

To use these named environments, set the following in your project configuration
(the example templates already use this setting):

```yaml
use_existing_envs: true
```

Launch the workflow with `--use-conda`, keeping the controller environment active:

```bash
snakemake --configfile config/data/local.yaml --use-conda --cores 8
```

Snakemake selects names from `conda.<tool>.name` in the core configuration.

Alternatively, set `use_existing_envs: false` and keep `--use-conda` to let
Snakemake create only the required environments from `conda.<tool>.yaml` and
reuse them under `.snakemake/conda/`. Relative YAML paths resolve from the main
`Snakefile` directory; absolute paths are also supported.

## Dependencies for additional analyses

A bypassed FCS job still requires the environment assigned to
`conda.singularity`, even when it does not execute a container.

Other selections require additional environments from `workflow/envs/` and may
need database/container paths from
[`config/databases_and_containers.yaml`](../config/databases_and_containers.yaml).
For example, FCS needs its image/database, Kraken2 needs a database, and MitoHiFi
uses a container and reference acquisition. Replace the shipped site-specific
paths before enabling those analyses. See the
[tool support matrix](stages.md#tool-support) for integration limitations.
