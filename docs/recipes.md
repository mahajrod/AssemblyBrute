# Focused workflow recipes

Merge these fragments into a complete project from the
[example catalog](examples/README.md). Paths below assume `out_dir: results/`,
`tmp_dir: results/tmp/`, `genome_prefix: sample`, and `ploidy: 1`, except where
overridden. Tool-parameter overrides belong in a private core selected by
`main_config_file`; see [configuration](configuration.md#which-file-to-edit).

## 1. HiFi read QC and filtering

Prepare HiFi FASTQ and the named `nanopore` environment in addition to
`common`/`kmer`; NanoPlot and NanoQC use it even for HiFi reads.

```yaml
input_datatypes: hifi
data_filtering: [hifi]
stage_list: [raw_read_qc, filter_reads, filtered_read_qc, kmer_qc]
```

Expect raw/filtered read-QC outputs under `qc/`, filtered/final read paths under
`data/hifi/`, and final-read k-mer models under `kmer/hifi/final/`. HiFi is not
FastQC-eligible in the shipped datatype configuration; NanoPlot/NanoQC provide
its read reports.

## 2. GenomeScope ploidy comparison

Keep the named `common`/`kmer` environments and prepared HiFi reads:

```yaml
input_datatypes: hifi
stage_list: [ploidy_check]
ploidy: 1
ploidy_test_list: [1, 2]
kmer_test_list: [21]
skip_kmer_smudgeplot: true
```

The parser unions these requests with the active model ploidy and configured
counter k values. Expect `.p1.genomescope.parameters` and
`.p2.genomescope.parameters` plus model reports/plots. These requests do not
rewrite the active `ploidy`.
The older template key `kmer_ploidy_test_list` is not read by this method.
Optional Smudgeplot belongs to `ploidy_check`; it is not a standalone stage.

## 3. Hifiasm contig assembly

Select `stage_list: [contig]`, `stage_coretools.contig: [hifiasm]`, and a preset
such as `coretool_option_sets.hifiasm: [l1_no_ul]` in the complete project. Set
input datatypes and ploidy to match your reads and assembly model.

Canonical FASTAs appear under `results/contig/<branch>/`; use
`requested_files.tab` to find the actual branch and haplotype names.

## 4. Compare two Hifiasm presets with QUAST

Use HiFi reads and the `common`, `kmer`, and `hifiasm` environments:

```yaml
input_datatypes: hifi
ploidy: 1
stage_list: [contig]
stage_coretools:
  contig: [hifiasm]
coretool_option_sets:
  hifiasm: [l0, l1_no_ul]
assembly_qc_level:
  contig: 1
```

For a small genome, set `bloom_filter_bits: 0` for both `l0` and `l1_no_ul`
under `parameters.micro.tool_options.hifiasm` in the private core; use nested YAML
mappings, not a literal dotted key. Both presets retain their remaining table
defaults. See [preset overrides](configuration.md#override-a-preset-in-the-private-core).

Expect two branch FASTAs, separate QUAST reports/parameter tables, and
`results/contig/sample.contig.stage_stats`. Shared QC level 1 does not add BUSCO or
Merqury jobs. See [summary-table limitations](outputs.md#summary-tables).

## 5. QUAST evaluation of an existing draft

Prepare one uncompressed draft FASTA under `<input_dir>/draft/fasta/`. Retain HiFi
inputs: included MitoHiFi rules access them during parsing even with mtDNA disabled.

```yaml
input_datatypes: hifi,draft
ploidy: 1
phasing_stage: draft_qc
stage_list: [draft_qc]
stage_coretools:
  draft_qc: [draft_qc]
coretool_option_sets:
  draft_qc: [def]
assembly_qc_level:
  draft_qc: 1
```

The QC tools use `common`. Expect import/length files under
`draft_qc/draft_qc_def/`, a QUAST report under its `assembly_qc/quast/` directory, and
`results/draft_qc/sample.draft_qc.stage_stats`. For multiple haplotypes, adjust
`ploidy`, filenames, and evaluation assumptions as described in [inputs](inputs.md).

## 6. Hi-C alignment and YaHS scaffolding

Prepare HiFi and paired Hi-C reads, and configure the library's enzyme set:

```yaml
input_datatypes: hifi,hic
data_filtering: [hic]
ploidy: 2
phasing_stage: contig
stage_list: [contig, hic_alignment, hic_scaffolding]
stage_coretools:
  contig: [hifiasm]
  hic_alignment: [arima]
  hic_scaffolding: [yahs]
coretool_option_sets:
  hifiasm: [l1_no_ul]
  arima: [def]
  yahs: [phased]
assembly_qc_level:
  contig: 0
  hic_alignment: 0
  hic_scaffolding: 0
```

Retain `p2: null` in the project's
`tool_manually_adjusted_features.genome_assemblers.lambda` mapping. The
Arima/YaHS tools use `common`.

Expect haplotype scaffold FASTAs under a branch such as
`hic_scaffolding/hifiasm_l1_no_ul@p2..arima_def..yahs_phased/`, with AGP/scaffolding
intermediates below it. `hic_alignment` alone requests no BAM through `rule all`;
YaHS supplies the alignment consumer. Contact-map products require additional QC
settings and are not requested at level 0.
