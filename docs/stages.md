# Stages and assembly branches

`stage_list` declares outputs to request and the ordered assembly chain. Snakemake
then schedules their dependencies. Selecting a stage does not isolate it from
prerequisite tools, and a stage with no requested targets can do no work.

## Start from reads or assemblies

Choose exactly one assembly starting point when building an assembly chain:

- `contig`: assemble reads using selected core tools/presets.
- `draft_qc`: import existing draft FASTAs; see [input naming](inputs.md#existing-draft-assemblies).

These two stages are mutually exclusive. Read-analysis-only runs need neither.
Put read analyses before the assembly starting point:

```yaml
stage_list: [raw_read_qc, filter_reads, kmer_qc, contig]
```

Each stage receives the **immediately preceding** stage as `prev_stage`.
Preprocessing entries are not skipped when constructing that chain. Interleaving
them between assembly transformations can fail or leave no branches; for
example, `[contig, raw_read_qc, polishing]` fails initialization. Keep each stage
name unique.

## Tool support

Available paths and integration limitations are listed below. Validation status
for runnable configurations is recorded in the [example catalog](examples/README.md).

| Stage | Tools/path | Requirements and limitations |
| --- | --- | --- |
| `raw_read_qc`, `filtered_read_qc` | FastQC, MultiQC, NanoPlot, NanoQC, TADbit | Datatype-dependent; TADbit uses Hi-C enzyme settings. |
| `raw_kmer_qc`, `kmer_qc` | Meryl, Yak, Jellyfish, FastK; GenomeScope, KrATER | Reads with the genome-size feature enabled. |
| `filter_reads` | Datatype-specific trimming/filtering | Controlled by `data_filtering` and datatype capabilities. |
| `ploidy_check` | GenomeScope; optional Smudgeplot | Genome-size-enabled reads and ploidy/k-mer selections. Standalone `smudgeplot` stage is not implemented. |
| `read_contamination_scan` | Kraken2 | Requires configured databases. |
| `mtdna` | MitoHiFi, MitoZ | MitoHiFi requires HiFi; MitoZ uses selected Hi-C/Illumina paths. Requires references, taxonomy settings, and containers. |
| `contig` | Hifiasm | Main reads and permitted phasing data determine assembly ploidy. |
| `contig` | Flye | Reads matching the preset; assembly ploidy is 1. |
| `contig` | Verkko | Phased output links have a known defect: both haplotypes link to the same source FASTA. |
| `contig` | NextDenovo | Stage integration is disabled/incomplete. |
| `draft_qc` | Existing assemblies and shared QC | Draft FASTA count must match `ploidy`; parsing also requires HiFi. Mutually exclusive with `contig`. |
| `polishing` | NextPolish2 | Previous assembly, HiFi, and Illumina-derived Yak databases at k=21/31. |
| `dedup` | purge_dups, HapSolo, combined purge | Previous assembly and read evidence; HapSolo needs BUSCO lineages. HapSolo/combined paths have source inspection only. |
| `hic_alignment` | Arima, `bwa_only`, Pairtools, Juicer | Previous assembly and Hi-C. Only Arima has recorded DAG checks; alignment needs a downstream consumer or explicit BAM target. |
| `hic_scaffolding` | YaHS, `threeddna` (3D-DNA) | Immediately follows compatible Hi-C alignment. Arima → YaHS has DAG checks; Juicer → 3D-DNA remains unchecked. |
| `ref_scaffolding` | RagTag stage | Stage include is disabled; shared RagTag QC rules are a separate path. |
| `gap_closing` | SAMBA | Previous assembly and `gap_closing_datatype`; presets can request phased reads. |
| `read_phasing` | Haplotype-specific k-mers/read extraction | Multiple haplotypes at `phasing_stage` and eligible reads. |

## Build an assembly chain incrementally

For a Hi-C workflow, a cumulative selection can grow from:

```yaml
stage_list: [contig]
```

to:

```yaml
stage_list: [contig, dedup, hic_alignment, hic_scaffolding]
stage_coretools:
  contig: [hifiasm]
  dedup: [purge_dups]
  hic_alignment: [arima]
  hic_scaffolding: [yahs]
coretool_option_sets:
  hifiasm: [l1_no_ul]
  purge_dups: [hifi_keep_hicov]
  arima: [def]
  yahs: [phased]
```

This fragment also needs Hi-C inputs, suitable ploidy/phasing settings, QC levels,
and the remaining project fields. Retain earlier stages and their preset labels
when extending a run so existing branches can be reused. To begin a different
workflow from saved FASTAs, prepare `draft` inputs and use `draft_qc` instead.

**Alignment target detail:** `request_hic_alignment_files()` currently returns
no explicit BAM targets. `[contig, hic_alignment]` alone therefore does not
request alignment BAMs through `rule all`. A downstream consumer such as
`hic_scaffolding`, or an explicit output-file target, must request them.

Branch generation pairs Juicer with `threeddna` (3D-DNA), while YaHS branches
exclude Juicer-derived predecessors. Incompatible combinations can leave no
scaffolding branches. Check nonempty requested outputs, not just a zero exit code.

Put explicitly requested `read_phasing` after the assembly chain: it uses
`phasing_stage` rather than providing the next transformed assembly. Other
consumers can request phased reads as dependencies without listing `read_phasing`.

## Understand branch labels and haplotypes

Selecting multiple presets/tools expands branches. For example:

```text
contig/hifiasm_l1_no_ul@p2/
dedup/hifiasm_l1_no_ul@p2..purge_dups_hifi_keep_hicov/
```

Each later tool/preset is appended with `..`; HapSolo lineages and references can
expand labels further. Keep the full label when targeting a file. Hifiasm groups
compatible correction options to reuse correction work, but different groups
still require separate jobs. Expensive QC can run separately for each branch.

Hifiasm/Verkko assembly ploidy depends on available and permitted phasing data and
preset/project values. `ploidy: 2` alone does not ensure two haplotypes. Flye is
assigned assembly ploidy 1. Single-assembly branches use `hap0`; multi-haplotype
branches use `hap1`, `hap2`, and so on. A combined FASTA is a concatenated
representation and does not add a separately inferred haplotype.

## Quality-control levels

`assembly_qc_level` is independent of `stage_list`; levels add target groups
subject to datatype availability, presets, and skip flags:

| Level | Added shared QC targets |
| --- | --- |
| 0 | None. Assembly outputs and their dependencies remain. |
| 1 | Stage statistics/QUAST; microchromosome candidates when enabled. |
| 2 | Merqury results through statistics dependencies. |
| 3 | Assembly tracks, repeats, and telomere analyses. |
| 4 | Whole-genome comparison and reference-based QC. |
| 5 | BUSCO results/intersections. |
| 6 | Coverage tracks. |
| 7 | Applicable Hi-C/contact-map visualization outputs. |

For example, `assembly_qc_level: {contig: 0}` disables shared QC for that stage.
Levels require the corresponding inputs and dependencies; a skip flag does not
bypass a level gate. Set a level for each selected assembly stage.
