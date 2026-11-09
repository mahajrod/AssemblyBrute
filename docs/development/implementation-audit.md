# Implementation audit

Inspected 2026-11-07 at commit `22920516`, with existing local changes left in
place. The load-order/input descriptions and initial support checks come from
Phase 1; subsequent execution evidence and follow-up findings are identified
below and in [validation](validation.md). This audit is not an assurance that
every configured tool works end to end.

## Implementation map

| Source | Responsibility |
| --- | --- |
| [`Snakefile`](../../Snakefile) | Loads helpers/configuration, discovers inputs, selects parameters, initializes stages, requests outputs, and includes rules. |
| [`workflow/include/read_configs.py`](../../workflow/include/read_configs.py) | Loads YAML files and tool/resource tables. |
| [`workflow/include/general.py`](../../workflow/include/general.py) | Recursive default filling and draft haplotype recognition. |
| [`workflow/include/general_parsing.py`](../../workflow/include/general_parsing.py) | Input discovery and stage-chain parsing. |
| [`workflow/include/stage.py`](../../workflow/include/stage.py) | Expands tools/presets into assembly branches and requests stage/QC files. |
| [`workflow/rules/`](../../workflow/rules/) | Executes preprocessing, assembly stages, tools, and QC. |
| [`config/templates/`](../../config/templates/) | Project examples; their filenames do not imply that only one stage is enabled. |
| [`workflow/envs/`](../../workflow/envs/), [`profile/`](../../profile/) | Tool environments and scheduler settings. |

## Configuration loading and precedence

Configuration precedence is implemented by pipeline code, so it cannot be
described simply as “command-line settings override YAML.” The actual sequence is:

1. Snakemake supplies the project configuration and command-line configuration.
2. `main_config_file` selects a core YAML file. Explicit project/CLI `input_dir`,
   `out_dir`, and `tmp_dir` are preserved across the core update. For other keys,
   `config.update(core_yaml)` replaces existing **top-level values** without
   recursively merging them.
3. Secondary-tool and database/container YAML files recursively fill absent
   entries through `copy_absent_entries`. Existing scalar/list values win.
4. Core-tool TSVs populate parameter-set tool dictionaries. For an existing
   parameter-set directory, each configured core-tool dictionary is reset before
   its table is loaded; subsequent default fallback supplies missing entries.
5. Cluster and datatype YAML files recursively fill missing entries. Skip flags
   are loaded only when the corresponding top-level key is absent.
6. Resource TSVs assign `threads`, `memory_mb`, and `time` dictionaries for their
   parameter-set columns. `resource_profile` chooses the table directory.
7. `Snakefile` fills absent entries in the selected `parameter_set` from
   `parameters.default`, discards other parameter sets, and selects secondary-tool
   options using `other_tool_option_sets`.

The loader and parameter fallback were exercised directly with assertions:

| Setting | Observed result |
| --- | --- |
| Project `input_dir`, `out_dir`, `tmp_dir` | Preserved when explicitly set; otherwise inherited from core. This exception was added and checked during dataset-output isolation. |
| Project `stage_list` | Retained with the repository's core YAML, which does not define it. |
| Explicit `skip_fastqc` | Retained when loading skip defaults. |
| Project `other_tool_option_sets.fastqc` | Retained; missing selections such as BUSCO filled from secondary-tool defaults. |
| Core-tool/resource tables and default fallback | Loaded and populated the selected parameter set. |

**Documentation consequence:** dataset paths can be set directly in the project
or CLI. Set `out_dir` and `tmp_dir` together; changing one does not derive the
other. A private core is still useful for tool/environment overrides. Other
colliding project/CLI keys remain subject to core precedence.
Keep paths relative to the repository root unless explicitly made absolute.

## Input recognition

Read discovery uses `input_dir/<datatype>/<format>/` and sorts matching files.
These are declarations in [`config/datatypes.yaml`](../../config/datatypes.yaml),
not proof that every conversion/tool combination has been executed:

| Datatypes | Declared formats/extensions | Recognition details |
| --- | --- | --- |
| `illumina`, `hic` | `fastq`: `.fq.gz`/`.fastq.gz`; `fasta`: `.fa.gz`/`.fasta.gz`; `mapped_bam`/`unmapped_bam`: `.bam` | Marked paired. Consecutive sorted files are paired; an odd count is rejected, and paired filenames must differ at no more than one character. |
| `hifi`, `lqccs`, `nanopore`, `simplex`, `duplex`, `ultralongnano`, `adaptivenano` | Same declared read formats/extensions | Marked single-end; datatype feature flags control applicable tools. |
| `draft` | `draft/fasta/*.fasta` | File count must equal `ploidy`; with multiple haplotypes, names must end in `hap1.fasta`, `hap2.fasta`, etc. A single assembly becomes `hap0`. |
| `reference` | `reference/<name>/*` | Parser expects one `.fasta`, `.whitelist`, and `.orderlist`; provide a `.syn` file too. Missing-file handling needs review. |

For each selected read datatype, discovery rejects missing files, multiple
recognized formats, or mixed extensions within a format. For example, `.fq.gz`
and `.fastq.gz` cannot coexist for the same datatype/format. The dry runs used
compressed FASTQ; BAM/FASTA read conversion and the paired-BAM behavior remain
unverified.

`ext_data` is handled separately in `Snakefile` for additional visualization
tracks. References to `parental` in assembler phasing logic do not establish a
working input interface: the shipped datatype configuration has no `parental`
entry. Do not advertise parental/trio input without further verification.

## Stages and tool support

**Dry-run passed** below means only the tested selection built a job graph with
synthetic inputs under Snakemake 7.32.4. **Code present** means rules/dispatch were
found but that selection was not dry-run. Phase 1 executed no tools; the Hifiasm
and draft-QC rows also identify the later executed examples.

| Stage | Tools/path | Required inputs or ordering | Evidence/status |
| --- | --- | --- | --- |
| `raw_read_qc` | FastQC, MultiQC, NanoPlot, NanoQC, TADbit as applicable | Selected read datatypes; Hi-C enzyme settings for TADbit | Dry-run passed with HiFi + Hi-C. |
| `raw_kmer_qc` | Configured counters and genome-size estimators | Datatypes with genome-size feature enabled | Dry-run passed with default counter selection. |
| `filter_reads` | Datatype-specific filtering/trimming | Reads and `data_filtering` | Dry-run passed. |
| `filtered_read_qc` | Same QC dispatch on filtered reads | Filtering outputs for selected datatypes | Dry-run passed. |
| `kmer_qc` | Meryl by default; configured Yak, Jellyfish, FastK; GenomeScope/KrATER | Final reads and selected counters/presets | Default Meryl path dry-run passed; alternative counters not checked. |
| `ploidy_check` | GenomeScope ploidy/k-mer combinations; optional Smudgeplot | Genome-size-enabled reads, ploidy/k-mer selections | Dry-run passed with Smudgeplot disabled by the default skip flag. |
| `read_contamination_scan` | Kraken2 | Reads and configured Kraken database | Reached database dependency; blocked by absent site-specific database. |
| `mtdna` | MitoHiFi and MitoZ | HiFi for MitoHiFi; Hi-C/Illumina for selected MitoZ paths; reference/taxonomy settings | HiFi + Hi-C dry-run passed; reference retrieval, containers, and execution unverified. |
| `contig` | Hifiasm | Main reads; permitted phasing data affects inferred assembly ploidy | Phase 1 diploid HiFi + Hi-C dry-run passed; Phase 2 HiFi-only primary assembly executed. Other modes remain unverified. |
| `contig` | Flye | Reads matching the chosen preset | `hifi` preset dry-run passed; stage assigns assembly ploidy 1. |
| `contig` | Verkko | Configured HiFi/Nanopore inputs and phasing choice | `default` preset dry-run passed with HiFi + Hi-C, but phased output links have a correctness defect; see below. |
| `contig` | NextDenovo | No usable stage implementation established | Rule include is disabled and rule body is a quoted TODO; test also failed on absent tool options. |
| `draft_qc` | Existing assemblies and shared QC | `draft` input; mutually exclusive with `contig`; current parsing also accesses HiFi | Phase 1 diploid import/level-5 QC dry-runs passed; Phase 4 single-draft level-1 QUAST/statistics executed. |
| `polishing` | NextPolish2 | Previous assembly, HiFi reads, and Illumina-derived Yak databases at k=21/31 | HiFi + Hi-C alone failed; adding Illumina passed. |
| `dedup` | purge_dups, HapSolo, combined purge | Previous assembly; selected read evidence; BUSCO lineage settings for HapSolo branches | purge_dups dry-run passed; HapSolo/combined selections are code present only. |
| `hic_alignment` | Arima, `bwa_only`, Pairtools, Juicer | Previous assembly, Hi-C reads, preset/phasing settings | Arima dry-run passed; other selections are code present only. |
| `hic_scaffolding` | YaHS, `threeddna` (3D-DNA) | Previous Hi-C alignment; branch generation pairs Juicer with 3D-DNA and excludes Juicer from YaHS branches | Arima → YaHS dry-run passed; Juicer → 3D-DNA not checked. |
| `ref_scaffolding` | RagTag stage | Previous assembly and reference | Stage include disabled. Test with explicit RagTag preset reached missing stage FASTA; not ready for a recipe. Shared RagTag QC rules are a separate path. |
| `gap_closing` | SAMBA | Previous assembly, `gap_closing_datatype`, optional phased reads | Dry-run passed with the template's HiFi/phased preset. |
| `read_phasing` | Haplotype-specific k-mers/read extraction | An available `phasing_stage`, its haplotypes, and eligible reads | Hifiasm → read phasing dry-run passed. |
| `smudgeplot` as a standalone stage | No matching stage dispatch | Listed in allowed stages but lacks normal wiring | Failed during stage initialization. Optional Smudgeplot belongs to `ploidy_check`; that path is still unverified. |

## Stage chains and outputs

- `stage_list` is processed in its supplied order. Each stage receives the
  immediately preceding stage as `prev_stage`; preprocessing stages are not
  automatically removed from that chain. Do not insert a preprocessing stage
  between assembly transformations without checking the resulting graph.
- `contig` and `draft_qc` are mutually exclusive assembly starting points.
  Preprocessing-only runs are valid. Assembly transformations need an upstream
  assembly chain; `hic_scaffolding` expects Hi-C alignment outputs.
- Stage expansion combines selected tools and option sets. Labels such as
  `hifiasm_l1_no_ul@p2` identify an initial assembly branch; later steps append
  `..tool_preset`. Multiple HapSolo lineages/references can expand branches further.
- Hifiasm/Verkko ploidy is inferred from available, permitted phasing inputs and
  preset/project settings. `ploidy: 2` alone does not guarantee two assembled
  haplotypes. Single-haplotype stages use `hap0`; phased stages use `hap1`, etc.
- `rule all` consumes files requested by the initialized stages. Dependencies can
  schedule analyses even when their stage is not explicitly selected. Conversely,
  requesting a stage whose enabled analyses produce no targets is not a useful
  validation.
- Initialization writes logs, `config/config.final.yaml`, `config/final_data.yaml`,
  and `requested_files.tab` under `out_dir`, even during a dry run. The saved config
  is a processed snapshot and can contain Python-specific YAML values; it has not
  been verified as a directly reusable project configuration.

## QC selection

`assembly_qc_level` is a per-stage mapping, separate from `stage_list` and the
individual `skip_*` flags. In `Stage.request_assembly_qc_files`, levels add these
target groups, subject to available datatypes, presets, and skip flags:

| Level | Added target groups |
| --- | --- |
| 0 | No shared assembly-QC targets. |
| 1 | Stage statistics/QUAST; microchromosome candidates when enabled. |
| 2 | Merqury results through stage-statistics dependencies. |
| 3 | Assembly tracks, repeats, and telomere analyses. |
| 4 | Whole-genome comparison and reference-based QC. |
| 5 | BUSCO results and intersections. |
| 6 | Coverage tracks. |
| 7 | Hi-C/contact-map visualization outputs when applicable. |

This describes target selection, not all dependencies a tool can introduce.
The contig template sets `draft_qc: 0` and `gap_closing: 0`, while the draft-QC
template sets `draft_qc: 7`. Therefore, changing only `stage_list` in the contig
template does not enable comprehensive draft QC. The audit separately checked
draft import at level 0 and draft QC at level 5; levels 6–7 remain unverified.

## Limitations to carry into the guide

1. **Output initialization (resolved locally).** Missing `out_dir` caused
   `Unable to configure handler 'debug_file'` in the initial tests. The logger now
   creates that directory before opening handlers; fresh nested roots were checked.
2. **Keep phased Verkko out of the verified quickstart.** In
   [`Verkko.smk`](../../workflow/rules/Stages/contig/Verkko.smk), both `hap1_fasta` and
   `hap2_fasta` link to `assembly.haplotype1.fasta`. Also, `--haplo-divergence` reads
   the `polishing_iterations` option rather than `haplo_divergence`. A successful
   DAG check does not detect these shell-command defects.
3. **Do not promise reference scaffolding, NextDenovo, or a standalone
   `smudgeplot` recipe.** Their incomplete wiring is distinct from optional
   analyses requiring external resources.
4. **Review reference side effects before examples.** If `.syn` is absent,
   initialization attempts to create one through a shell command; missing required
   reference files are indexed without an explicit empty-list check. Audit runs
   supplied all reference sidecars and used synthetic copies only.
5. **Replace site-specific database/container paths.** The shipped database config
   points at external locations. A dependency declaration or local container
   filename is not proof that the asset exists or that execution will work.
6. **Treat scheduler support separately.** `profile/slurm/` uses the legacy
   `cluster: sbatch ...` interface. `profile/slurm8.30/` uses `executor: slurm` and
   a different configuration shape. Neither profile was submitted to a scheduler;
   the latter's directory name is not a compatibility guarantee.

## Handoff to Phase 2

Use Hifiasm as the first assembly recipe candidate and clearly list every enabled
analysis and dependency. At the original audit, the now-removed contig template
also enabled read QC, filtering, deduplication, alignment, and scaffolding, and requested HiFi,
Hi-C, simplex, and reference inputs. It was not a minimal assembly-only example.
The old read-QC/ploidy templates also requested simplex/reference inputs; the old draft-QC
template enabled several downstream stages. The current configurations are in
the [example catalog](../examples/README.md). Tailor examples explicitly rather
than presenting an unmodified template as a minimal run.

Before publishing the quickstart: obtain a small real dataset, validate its
biological metadata, prepare portable project/core files, verify environment
creation as well as existing-environment use, and run the documented workflow.
Report code defects separately; this audit changes documentation only.

### Phase 2 follow-up

The historical HiFi-only execution used a private core with Hifiasm's initial
Bloom filter disabled and named environments; its results remain in
[validation](validation.md#executed-quickstart). The current runnable examples are
listed in the [example catalog](../examples/README.md). Additional checks established that
the `low` resource table lacks `combine_long_reads` entries required during
parsing. The previously observed Conda YAML path failure was fixed by resolving
specifications from the main `Snakefile` directory. All ten then-current example
configurations passed DAG/environment-list checks in both selection modes;
fresh Conda solves remain unverified. See
[installation](../installation.md#create-tool-environments) and the
[validation record](validation.md) for their practical limits.

### Phase 4 follow-up

Report/recipe checks identified additional limits, documented in
[outputs](../outputs.md) and [troubleshooting](../troubleshooting.md):

- The original `gather_stats_per_stage_parameter` used a literal `results/` QC
  path; it now uses `out_dir`. Duplicate headers and bare `-m` for an empty
  Merqury datatype selection remain separate limitations.
- Included `MitoHiFi.smk` definitions access HiFi during parsing, so the
  draft-only check fails even when mtDNA is not selected.
- `Stage.request_ploidy_check_files` reads `ploidy_test_list` and `kmer_test_list`,
  not the older template's `kmer_ploidy_test_list`.
- The Merqury rule supplies a single assembly or a two-haplotype pair; additional
  polyploid haplotypes are not passed to that invocation.

An isolated level-1 draft-QC run on an existing assembly produced QUAST
reports and consolidated statistics successfully with `out_dir: results/` and
HiFi retained for initialization. No workflow-code fixes were made.
