# Yeast: stage tests with HiFi and Hi-C

The six standalone configurations in
[`config/templates/yeast/`](../../config/templates/yeast/) use the existing yeast
inputs under **`test_data/yeast/`** and write to **`results/yeast/`**. Run commands from the repository root
with the Snakemake 7 controller active. Complete [installation](../installation.md)
and prepare the named Conda environments needed by the selected tools.

**Validation: all six DAGs were checked. The namespaced workflows have not been
executed.** See [validation](../development/validation.md) for the checks and their limits.

## Inputs and shared settings

These tests select `test_data/yeast/` as their input root and preserve the datatype
and format subdirectories:

```text
test_data/yeast/
├── hifi/fastq/SRR13577846.half.fastq.gz
├── hic/fastq/SRR7126301.200x_1.fastq.gz
│             SRR7126301.200x_2.fastq.gz
└── draft/fasta/                  # Two FASTAs needed only for draft_qc.yaml
```

The filenames of these local subsets do not specify a reproducible downsampling
recipe. Existing
external FASTQ files can be linked into these input directories; do not commit
sequencing data or generated results.

The HiFi/Hi-C files in this checkout are symlinks to the original files under
`input/hifi/` and `input/hic/`. The two draft FASTAs are ordinary file copies.
The dataset was moved into `test_data/yeast/` without changing file contents;
relative FASTQ links were adjusted for the new location.

Each project explicitly selects:

```yaml
main_config_file: config/core.yaml
input_dir: test_data/yeast/
out_dir: results/yeast/
tmp_dir: results/yeast/tmp/
```

To use another input directory without editing the YAML, override it on the CLI:

```bash
snakemake --configfile config/templates/yeast/contig.yaml --use-conda -c 8 \
  --config input_dir=/data/yeast/ --dry-run
```

Replace `/data/yeast/` with your prepared input root. Relative paths are resolved
from the working directory, normally the repository root; absolute paths are
also accepted. Changing `input_dir` does not change the output directory.

Explicit project/CLI `input_dir`, `out_dir`, and `tmp_dir` values override core
defaults. Other colliding core settings still follow the
[documented precedence](../configuration.md#precedence-core-values-can-replace-project-values).
Initialization creates the output directory and writes logs/config snapshots even
in a dry run.

The templates retain `resource_profile: high`, `parameter_set: micro`, existing
named environments, genome prefix `testSacMyc1`, and model `ploidy: 2` from the
working local example. This is a test configuration, not an inferred biological
ploidy. Hi-C uses the current `Arima` enzyme setting; check it against your library
before interpreting scaffolding results. BUSCO and FCS/adaptor screening are skipped.

## Check each DAG separately

Each command below only builds the requested graph; it does not run scientific
tools. `-c 8` specifies the scheduling CPU budget, not a RAM limit.

```bash
snakemake --configfile config/templates/yeast/read_qc.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/yeast/all_kmer_counters.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/yeast/ploidy_check.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/yeast/contig.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/yeast/draft_qc.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/yeast/hic_scaffolding.yaml --use-conda -c 8 --dry-run
```

For shell commands and exact paths, add `--printshellcmds`.

| Configuration | Requested analyses | Main output areas under `results/yeast/` |
| --- | --- | --- |
| `read_qc.yaml` | Raw/filtered read QC, filtering, and default k-mer QC. | `qc/`, `data/`, `kmer/` |
| `all_kmer_counters.yaml` | Meryl, Yak, Jellyfish, and FastK with genome-size reports. | `kmer/hifi/final/` |
| `ploidy_check.yaml` | GenomeScope models 1–4 and Smudgeplot. | `kmer/`, including `genomescope/` and `smudgeplot/` |
| `contig.yaml` | Read preparation, Hifiasm assembly, and shared assembly QC. | `contig/`, `error_correction/` |
| `draft_qc.yaml` | Import and shared QC of two existing draft FASTAs. | `draft_qc/` |
| `hic_scaffolding.yaml` | Assembly, Purge_dups, Arima alignment, and phased YaHS scaffolding with QC. | `contig/`, `dedup/`, `hic_alignment/`, `hic_scaffolding/` |

The assembly stages use QC level 5, with BUSCO skipped. They can still request
Merqury, repeats, telomere analysis, and whole-genome alignments.

## Prepare draft-QC inputs

`draft_qc.yaml` requires exactly two uncompressed `.fasta` inputs ending in
`hap1.fasta` and `hap2.fasta`. Existing copies from the earlier completed assembly
remain valid inputs; they are not automatically refreshed from `results/yeast/`.

To evaluate a newly completed namespaced contig run, first verify its two FASTAs,
then prepare copies in an otherwise empty draft FASTA directory:

```bash
test -s results/yeast/contig/hifiasm_l1_no_ul@p2/testSacMyc1.contig.hap1.fasta
test -s results/yeast/contig/hifiasm_l1_no_ul@p2/testSacMyc1.contig.hap2.fasta
mkdir -p test_data/yeast/draft/fasta
cp -n results/yeast/contig/hifiasm_l1_no_ul@p2/testSacMyc1.contig.hap1.fasta test_data/yeast/draft/fasta/testSacMyc1.hap1.fasta
cp -n results/yeast/contig/hifiasm_l1_no_ul@p2/testSacMyc1.contig.hap2.fasta test_data/yeast/draft/fasta/testSacMyc1.hap2.fasta
```

`cp -n` preserves existing files. Review which assembly should be evaluated before
intentionally replacing those copies. Keep the HiFi/Hi-C inputs when running draft QC: included rules currently
require HiFi during initialization.
The draft-QC dry run needs these inputs too; a contig dry run does not create FASTAs.

## Execute a selected test

Remove `--dry-run` from the chosen command when ready to spend the resources.
For example:

```bash
snakemake --configfile config/templates/yeast/contig.yaml --use-conda -c 8 \
  --printshellcmds --rerun-incomplete
```

Run one configuration at a time in this checkout. Tests for the same dataset
reuse compatible intermediates in `results/yeast/`; they are not independent
clean-room runs. Changing inputs or biological settings requires reviewing which
outputs can be reused. Different configurations also replace the dataset's
`config/config.final.yaml` and `requested_files.tab` with their current selection.

Expected contig files after execution include:

```text
results/yeast/
├── config/
├── log/
├── requested_files.tab
├── data/
├── qc/
├── kmer/
└── contig/
    ├── testSacMyc1.contig.stage_stats
    └── hifiasm_l1_no_ul@p2/
        ├── testSacMyc1.contig.hap1.fasta
        └── testSacMyc1.contig.hap2.fasta
```

Existing outputs directly under `results/` are not moved or adopted automatically.
A run into the new namespace may therefore schedule substantial recomputation.

## Add another dataset

Create its own project configurations and set all three paths explicitly, for
example `input_dir: test_data/scam5/`, `out_dir: results/scam5/`, and
`tmp_dir: results/scam5/tmp/`. Adjust species, prefix, ploidy, datatypes, and library
settings for that sample. Inspect each new DAG before execution.
No symlink or work-directory trick is needed for nested output directories:
Hifiasm correction links and QC-statistics gathering now use the configured paths.
Dataset separation organizes files; it does not reserve memory or validate every
optional tool integration.
