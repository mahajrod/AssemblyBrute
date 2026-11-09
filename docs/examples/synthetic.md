# Synthetic: diploid stage tests

The six configurations in
[`config/templates/synthetic/`](../../config/templates/synthetic/) read from
`test_data/synthetic/` and write to `results/synthetic/`. Run commands from the
repository root with Snakemake 7 and the tool environments from
[installation](../installation.md).

**Validation: all six DAGs checked in named and YAML Conda modes. No scientific
jobs were executed.**

## Source and generated inputs

The [generator](../../workflow/scripts/generate_synthetic_reads.py) creates
`chr1` (750,000 bp) and `chr2` (250,000 bp), totaling **1 Mb**, with target GC 40%.
It needs no input files. Numeric seed 42 determines the reads and variants;
reference generation uses seed 43 (`--seed + 1`).
Haplotype 1 matches the generated reference; haplotype 2
contains 10,000 fixed SNPs (1% divergence), introduced before sequencing errors:

| Input | Records | Length | Coverage |
| --- | --- | --- | --- |
| HiFi-like | 20,023 reads | 1,000–3,000 bp | 40.00× |
| Illumina WGS | 166,668 pairs | 2×150 bp; insert approximately 350 bp | 50.00× |
| Artificial Hi-C | 120,000 pairs | 2×150 bp | 36.00× |

The five compressed FASTQ files total approximately 39.4 MiB.
Coverage is relative to the 1 Mb single-haplotype reference span: about 20×
HiFi and 25× WGS per haplotype. Read pairs originate from the same haplotype.
All reads use a substitution probability of 0.001 and constant Q30, without
indels or adapters. Hi-C read starts use GATC anchors, with a distance-biased cis
model and a 10% probability of trans contacts. These are technical fixtures;
the contacts do not represent experimental chromosome structure.

```text
test_data/synthetic/
├── synthetic.fasta
├── hifi/fastq/synthetic.fastq.gz
├── illumina/fastq/synthetic_1.fastq.gz
│                  synthetic_2.fastq.gz
├── hic/fastq/synthetic_1.fastq.gz
│             synthetic_2.fastq.gz
├── draft/fasta/synthetic.hap1.fasta
│               synthetic.hap2.fasta
├── variants.tsv
├── manifest.json
└── truth.tsv.gz
```

The two draft FASTAs contain the simulated haplotypes, so draft QC can run
without assembling the reads. `variants.tsv` lists SNP positions and ref/alt
alleles using zero-based source coordinates. `manifest.json` records counts,
coverage per haplotype, models, origins, and SHA-256 hashes. `truth.tsv.gz` records
each read's haplotype-specific sequence ID, interval, strand, and sequencing
substitution count using zero-based, half-open coordinates.

Generate the reference and inputs from the script and numeric parameters:

```bash
python workflow/scripts/generate_synthetic_reads.py \
  test_data/synthetic/ --chromosome-lengths 750000 250000 \
  --gc-content 0.4 --seed 42 --heterozygosity 0.01
```

Existing `DESCRIPTION` markers are preserved. The generator refuses to overwrite
data or metadata; choose a different output directory to regenerate an existing
dataset. Diploid simulation is the default; `--heterozygosity 0` selects the
single-haplotype model. Only `DESCRIPTION` markers are allowed in Git;
generated references, reads, draft FASTAs, truth
records, manifests, and results are ignored.

## Shared settings

Each configuration selects:

```yaml
main_config_file: config/core.yaml
input_dir: test_data/synthetic/
out_dir: results/synthetic/
tmp_dir: results/synthetic/tmp/
```

The templates use `resource_profile: high`, `parameter_set: micro`, existing
named environments, genome prefix `synthetic`, and `ploidy: 2`.
Contig assembly uses HiFi plus Hi-C to request both haplotypes. Hi-C uses
`MboI` and YaHS `phased`. Contig, draft, and scaffolding
QC use level 2 for QUAST and Merqury; BUSCO, FCS/adaptor screening,
microchromosome detection, and mitochondrial analyses are disabled.

The 1 Mb size gives an allelic frequency bin above the fixed threshold of 1000
used by `smudgeplot all` 0.5.3; the preceding 8.3 kb reference did not reach it.
Resource settings retain repository defaults; for a smaller Hifiasm
Bloom filter, see [preset overrides](../configuration.md#override-a-preset-in-the-private-core).

## Check each DAG separately

```bash
snakemake --configfile config/templates/synthetic/read_qc.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/synthetic/all_kmer_counters.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/synthetic/ploidy_check.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/synthetic/contig.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/synthetic/draft_qc.yaml --use-conda -c 8 --dry-run
```

```bash
snakemake --configfile config/templates/synthetic/hic_scaffolding.yaml --use-conda -c 8 --dry-run
```

| Configuration | Selected inputs | Requested analyses |
| --- | --- | --- |
| `read_qc.yaml` | HiFi, Illumina, Hi-C | Raw/filtered read QC, filtering, and default k-mer QC. |
| `all_kmer_counters.yaml` | HiFi, Illumina | Meryl, Yak, Jellyfish, and FastK with genome-size reports. |
| `ploidy_check.yaml` | HiFi, Illumina | GenomeScope models 1–4 and Smudgeplot. |
| `contig.yaml` | HiFi, Hi-C | Hifiasm haplotype assemblies and shared QC. |
| `draft_qc.yaml` | HiFi, draft | Import and QUAST/Merqury evaluation of both simulated haplotypes. |
| `hic_scaffolding.yaml` | HiFi, Hi-C | Assembly, Arima alignment, and phased YaHS scaffolding. |

## Execute a selected test

Remove `--dry-run` from a selected command, for example:

```bash
snakemake --configfile config/templates/synthetic/contig.yaml --use-conda -c 8 \
  --printshellcmds --rerun-incomplete
```

Tests share compatible intermediates under `results/synthetic/`; run one
configuration at a time. Each selection rewrites its processed configuration and
`requested_files.tab`.

Expected assembly outputs after execution include:

```text
results/synthetic/
├── config/
├── log/
├── requested_files.tab
├── data/
├── qc/
├── kmer/
├── contig/hifiasm_l1_no_ul@p2/
│   ├── synthetic.contig.hap1.fasta
│   └── synthetic.contig.hap2.fasta
├── draft_qc/draft_qc_def/
│   ├── synthetic.draft_qc.hap1.fasta
│   └── synthetic.draft_qc.hap2.fasta
└── hic_scaffolding/hifiasm_l1_no_ul@p2..arima_def..yahs_phased/
    ├── synthetic.hic_scaffolding.hap1.fasta
    └── synthetic.hic_scaffolding.hap2.fasta
```

Stage summaries use `results/synthetic/<stage>/synthetic.<stage>.stage_stats`.
The [example catalog](README.md) describes validation levels; detailed checks
are recorded in [validation history](../development/validation.md#self-contained-synthetic-reference).
