# Documentation validation

## Phase 1 scope

Checked on 2026-11-07 against commit `22920516` with existing local changes left
untouched. These are **job-graph checks**, not completed assembly runs. The local
`workflow/envs/conda.common.yaml` already had changes; no environments were built
or updated during this audit.

### Setup

- Snakemake **7.32.4**, Conda available on `PATH`, existing named tool environments
  selected through the template's `use_existing_envs: True` setting.
- A temporary working directory under `/tmp/assemblybrute-phase1/work/`, with
  symlinks to the repository's `Snakefile`, `workflow/`, `config/`, and `resources/`.
  Snakemake metadata and its redirected cache were kept under the temporary tree.
- Synthetic compressed FASTQ: one 200-base read in each HiFi file and each mate
  of the Hi-C/Illumina pairs (`reads_1.fastq.gz`, `reads_2.fastq.gz`). These fixtures
  establish file discovery only; they cannot produce meaningful assemblies.
- Two draft FASTAs named `example.hap1.fasta` and `example.hap2.fasta`. Reference
  tests supplied an isolated `reference/example/` directory with `.fasta`,
  `.whitelist`, `.orderlist`, and `.syn` files.
- Each case copied `config/core.yaml`, changed its `input_dir`, `out_dir`, and
  `tmp_dir` to temporary paths, and precreated `out_dir`. Project YAML was copied
  from the **tracked** `config/templates/Example01.yeasts.contig.yaml`, then given
  the case's `main_config_file`, `stage_list`, and input/tool overrides below.
  Remaining template values and skip flags were retained. The separate untracked
  `Example01.yeasts.contig_.yaml` was not used.

### Commands and evidence

From the temporary working directory, the command pattern was:

```bash
snakemake --snakefile /tmp/assemblybrute-phase1/work/Snakefile \
  --dry-run --configfile /tmp/assemblybrute-phase1/CASE/project.yaml \
  --cores 8 --nolock --quiet --use-conda
```

Replace `CASE` with a row's case name. The command assumes the temporary fixtures
and generated configurations still exist; it is an audit command, not a user
quickstart. Temporary case folders retain `core.yaml`, `project.yaml`,
`dry-run.log`, and processed configuration/target files under `results/`.
Those local artifacts are not committed and may disappear when `/tmp` is cleaned.

All cases use `input_datatypes: hifi,hic` unless noted. Core-tool selections and
option sets below are overrides within their existing dictionaries; other
template selections remain in place.

| Case | `stage_list` in order | Additional changes | Result |
| --- | --- | --- | --- |
| `read-qc` | `raw_read_qc, filter_reads, filtered_read_qc, kmer_qc` | None | Passed; 37 jobs. |
| `raw-kmer` | `raw_kmer_qc` | None | Passed; 15 jobs. |
| `ploidy` | `ploidy_check` | Default `skip_kmer_smudgeplot: True` retained | Passed; 15 jobs. |
| `contig` | `contig` | Hifiasm `l1_no_ul` from template | Passed; 110 jobs. |
| `draft` | `draft_qc` | Inputs `hifi,hic,draft`; `phasing_stage: draft_qc`; inherited QC level 0 | Passed; 8 jobs; import/length outputs only. |
| `draft-with-qc` | `draft_qc` | Same inputs/phasing as `draft`; `assembly_qc_level.draft_qc: 5` | Passed; 91 jobs. |
| `scaffolding` | `contig, dedup, hic_alignment, hic_scaffolding` | Template purge_dups/Arima/YaHS selections | Passed; 330 jobs. |
| `read-phasing` | `contig, read_phasing` | None | Passed; 120 jobs. |
| `flye` | `contig` | `stage_coretools.contig: [flye]`; template `hifi` preset | Passed; 54 jobs. |
| `gap-closing` | `contig, gap_closing` | Template SAMBA `phased` preset | Passed; 127 jobs. |
| `mitochondrial` | `mtdna` | Template MitoHiFi/MitoZ settings | Passed; 15 jobs. |
| `polishing` | `contig, polishing` | HiFi + Hi-C only | Failed: `yak_reads` input function requires missing `illumina` data. |
| `polishing-with-illumina` | `contig, polishing` | Inputs `hifi,hic,illumina` | Passed; 215 jobs. |
| `verkko` | `contig` | `stage_coretools.contig: [verkko]`; `coretool_option_sets.verkko: [default]` | Passed; 98 jobs; source defects remain. |
| `ref-scaffolding` | `contig, ref_scaffolding` | Inputs `hifi,hic,reference` | Failed during initialization: template lacks `coretool_option_sets.ragtag`. |
| `ref-scaffolding-configured` | `contig, ref_scaffolding` | Inputs `hifi,hic,reference`; `coretool_option_sets.ragtag: [default]` | Failed: missing reference-scaffolding stage FASTA, with no active stage producer. |
| `nextdenovo` | `contig` | `stage_coretools.contig: [nextdenovo]`; `coretool_option_sets.nextdenovo: [def]` | Failed: missing tool options. Source also has no active NextDenovo rule. |
| `standalone-smudgeplot` | `smudgeplot` | None | Failed: missing `stage_coretools.smudgeplot` during initialization. |
| `contamination` | `read_contamination_scan` | Template database configuration | Blocked: configured Kraken `nt` database path is absent. |

Job counts include the aggregate `all` job and dependencies and apply only to
these fixtures/configurations. Thirteen selections passed; six exposed missing
inputs, configuration gaps, or incomplete wiring. The first attempt also exposed
the missing-output-directory logging failure; the listed runs precreated it.

### Configuration and archive checks

The actual `read_configs.py` loader was executed with its real helper functions,
YAML files, pandas table loading, and the subsequent default-parameter fallback.
Assertions confirmed core-file precedence for project paths, retention of project
stages and explicit skip flags, recursive filling of secondary-tool selections,
and population of Hifiasm options/resource tables. This was a focused loader
check; individual command-line override combinations were not separately run.

All eight original documentation files, including assets, were moved beneath
`docs/old/`. SHA-256 comparisons verified unchanged contents and retained relative
paths within the archive.

### Limits of Phase 1 checks

- Tool execution, output correctness, scientific interpretation, performance,
  full genome datasets, and any ploidy/read mode not listed above.
- Conda environment solving/creation and the availability of every declared
  executable, container, database, and remote reference service.
- Slurm submission and compatibility of `profile/slurm8.30/`. Snakemake 9.27.1
  was detected locally but was not used for these workflow checks.
- BAM/FASTA read conversions, external visualization-track inputs, all alternate
  assemblers/presets, HapSolo/combined purge, non-Arima alignment, 3D-DNA, and the
  optional Smudgeplot path within `ploidy_check`.

Successful dry runs do not clear the Verkko shell-command defects or the other
limitations recorded in the [implementation audit](implementation-audit.md).

## Phase 2: first-run guide

This section records the historical HiFi-only assembly check. Its tutorial and
example YAML have been removed; current runnable examples are the
[six yeast stage tests](../examples/yeast.md). Checks ran under
`/tmp/assemblybrute-phase2/work/`, with independent inputs, outputs, Snakemake
metadata, and cache files. Repository workflow/config/resource files were read
through symlinks; the shared core and existing user inputs/results were not
modified.

### Dataset and configuration

The downloaded ENA FASTQ for
[SRR13577846](https://www.ebi.ac.uk/ena/browser/view/SRR13577846) matched the
published metadata:

| Property | Checked value |
| --- | --- |
| Compressed size | 618,644,816 bytes |
| MD5 | `d57a6e681f2e6f3ea205ef6c4845d2da` |
| Reads | 117,525 |
| Total read bases | 1,103,860,593 |

`gzip -t` passed. A full FASTQ scan also confirmed record structure and matching
sequence/quality lengths. Experiment metadata identifies the material as HiFi
sequencing of the yeast INSC1020; no provenance was inferred from the preexisting
local file named `SRR13577846.half.fastq.gz`.

The selected configuration uses HiFi alone, one Hifiasm `l1_no_ul` branch,
`ploidy: 1`, no trimming, shared QC level 0, and skipped FCS/adaptor filtering.
The private core sets the input/output/temp directories, maps the bypassed FCS
environment to `common`, and overrides `bloom_filter_bits: 0` for the micro
parameter set. Remaining preset fields are inherited from the shipped tool table.

During the recorded execution, the full archive was retained in
`quickstart-downloads/`, outside input discovery. The workflow input is the first
40,000 complete records copied to `quickstart-input/hifi/fastq/yeast.40k.fastq.gz`
by the preparation block used at that time. These historical paths identify the
original evidence. This ordered subset is not a random sample or a prescription
for production coverage.

The completed subset passed record/CRC checks and contained **40,000 reads and
375,169,184 bases**. Its copied FASTQ records matched the source prefix with SHA-256
`956ebcb3bdd505912ea31944c6ff3191913f9f9a5697b4a366ea0f70e675643c`.

### Documentation checks

- The exact Python setup block was executed in a temporary directory. Assertions
  checked its paths, environment mapping, and Bloom-filter override. A second run
  raised `FileExistsError` and preserved the existing core bytes.
- The subset block's final-file publication was checked with compact complete
  and short FASTQ fixtures: only a complete 40,000-record file became discoverable,
  short input left only `.part`, and rerunning preserved an existing result.
- Input-discovery checks exercised the real parser: one HiFi FASTQ extension was
  accepted, mixed extensions were rejected, and simultaneous FASTA/FASTQ formats
  were rejected.
- Bash syntax was checked for the current documentation's command blocks. Local
  Markdown links, YAML parsing, and whitespace were checked separately.

### Environment limits

Named environments were used, with these installed versions:

| Environment | Selected package versions |
| --- | --- |
| Controller | Snakemake 7.32.4, Python 3.11.18, PuLP 2.7.0, pandas 2.3.3, NumPy 2.4.6, PyYAML 6.0.3. |
| `common` | GenomeScope2 2.0.1, gfatools 0.5, samtools 1.19.2, KrATER 2.7.13. |
| `kmer` | Meryl 1.4.2. |
| `hifiasm` | Hifiasm 0.25.0-r726. |

Automatic environment resolution was checked using the same project with
`--config use_existing_envs=False --list-conda-envs`. It failed on an incorrectly
relative `conda.kmer.yaml` path before creating an environment. This historical
path failure is now fixed; see [Conda YAML path resolution](#conda-yaml-path-resolution).
The `low` resource
profile also failed at parse time on a missing `combine_long_reads` table entry;
the published example uses `high`/`micro` instead. Fresh Conda solves/creation from
the manual installation commands remain unverified.

### Executed quickstart

The accepted historical run used the completed, validated 40,000-read subset,
the small-genome core override, four scheduling cores, and named Conda environments.
A dry run, workflow execution, and final up-to-date dry run were performed in
the isolated working directory.

All three commands returned 0. The execution completed 17 scheduled jobs in
approximately 343 seconds; log directories from preparation attempts already
existed. The last command reported `Nothing to be done` with all requested files
present and up to date. Correction, Meryl counting, GenomeScope estimation,
assembly, GFA conversion, FCS bypass, indexing, and length extraction ran with the
named environments. No container was required.

| Observed primary assembly property | Value |
| --- | --- |
| Contigs | 26 |
| Total bases | 12,264,312 |
| Longest contig | 1,505,958 bases |
| Contig N50 | 808,003 bases |

An independent FASTA scan checked unique IDs, nonempty sequences, and valid
ACGTN characters. Its sequence lengths matched every row of both
`yeast.contig.hap0.fasta.fai` and `yeast.contig.hap0.len`. These are observations
from this dataset, tool version, and configuration, not guarantees for other
inputs or evidence of completeness/consensus quality.

Earlier full-depth attempts were interrupted, and an attempt against a still
being-written subset was discarded after failing. Only the subsequent run against
the closed, validated subset is reported as successful. The preparation block now
keeps unfinished output under `.part` until it can be published as pipeline input.

Logs and numeric check results are retained locally under
`/tmp/assemblybrute-phase2/`, including `execution.log`, `up-to-date.log`,
`execution-results.json`, `subset-stats.json`, and `assembly-stats.json`. These
temporary artifacts are not committed. The run does not verify a fresh Conda
installation, full-depth assembly, optional QC, phased assembly, or Slurm support.

## Phase 3: configuration, stages, and execution

The [configuration](../configuration.md), [stage](../stages.md), and
[execution](../running.md) guides were checked against the actual loader, stage
class, resource helpers, shipped profiles, and installed Snakemake 7.32.4 CLI.
No additional scientific jobs or Slurm submissions were run in this phase.

### Job-graph and CLI checks

Fourteen cases used synthetic inputs and isolated output directories under
`/tmp/assemblybrute-phase3/`. Projects were derived from the yeast example or the
tracked contig template; private cores redirected inputs/outputs and disabled
shared QC for the stage-chain checks. All graph checks used `--dry-run` with
named environments. Twelve cases returned 0; the two invalid configurations
failed as expected.

| Case | Confirmed behavior |
| --- | --- |
| CLI precedence | `genome_prefix=cli_yeast` survives; project/CLI `out_dir` values are replaced by the selected core path. Neither discarded output path is created. |
| Two presets | Both `hifiasm_l0@p1` and `hifiasm_l1_no_ul@p1` appear in requested outputs. |
| Two tools | Hifiasm and `flye_hifi@p1` branches appear. |
| Documented comparison fragment | The guide's exact YAML fragment produces three branches: two Hifiasm presets and Flye. |
| Hi-C chain | `contig → dedup → hic_alignment → hic_scaffolding` builds the purge_dups/Arima/YaHS branch with HiFi + Hi-C. |
| Alignment-only selection | `contig, hic_alignment` initializes successfully but adds no explicit BAM/alignment target through `rule all`. |
| Bad order | `contig, raw_read_qc, polishing` fails because the preceding stage does not provide an assembly chain. |
| Both assembly starts | Selecting both `contig` and `draft_qc` raises the mutual-exclusion error. |
| Legacy profile | The shipped Slurm profile parses in cluster-mode dry runs. |
| Adapted profile | A private profile, dummy account, copied cluster YAML, and altered queue name parse with `--jobs 20 --local-cores 2`; the effective queue matches the private file. |
| Summary | Output status-table columns are produced. |
| Explicit file target | An exact haplotype length-table path resolves and builds its dependency graph. |
| Resource overrides | `hifiasm_correct` prints two tool threads, `cpus=2`, and the overridden `mem` value. |
| Forced length job | `--forcerun get_seq_len` produces a valid dry-run plan. |

### Focused loader/resource checks

The real loader and default fallback were executed with private copies of the
three resource tables. Assertions confirmed that:

- The existing `l1_no_ul` preset inherits its remaining fields while retaining a
  micro-specific `bloom_filter_bits: 0` override when no micro tool directory exists.
- Custom table cells replace a competing core resource dictionary: Hifiasm gets
  two threads, 2,500 MB, and `01:23:00` from the selected private profile.
- Adding a micro tool directory resets its Hifiasm dictionary before fallback;
  a table's changed purge level takes effect and the earlier Bloom-filter override
  is removed.
- Recursive defaults preserve explicit `null`; `parse_option` omits literal
  `default` and retains numeric zero.
- Queue CPU capping and linear/exponential memory retry calculations match the
  helpers described in the guide.

The installed `--help` and official Snakemake/Slurm references were consulted for
generic CLI/submission semantics. Only parsing and graph construction were
checked for cluster examples: account/partition availability, scheduler job
execution/cancellation, compute-node access, and the executor-style profile remain
unverified. The Phase 2 narrow yeast run remains the end-to-end execution evidence.

Temporary artifacts include `results.json`, `results-extra.json`, per-case
`check.log` files, effective configuration snapshots/target lists, and
`loader/check.txt` under `/tmp/assemblybrute-phase3/`. They are not committed.

## Phase 4: outputs, recipes, and troubleshooting

The [output guide](../outputs.md), [recipes](../recipes.md), and
[troubleshooting page](../troubleshooting.md) were checked against requested target
paths, active rule definitions, report scripts, the completed yeast run, and
primary QUAST/Merqury/BUSCO documentation. Recipe fragments were read directly
from `recipes.md` and merged into the complete yeast project for testing.

### Recipe and diagnostic checks

Ten graph/initialization cases ran in separate directories under
`/tmp/assemblybrute-phase4/`. The five new recipe fragments built their graphs;
the primary-assembly recipe links to the already executed Phase 2 example.
Nine cases returned 0, including diagnostic graphs that expose limitations;
the draft-only case failed as expected.

| Case | Checked result |
| --- | --- |
| Read QC/filtering | Raw/filtered NanoPlot statistics and filtering targets requested with `data_filtering: [hifi]`. |
| Ploidy models | Explicit `ploidy_test_list: [1, 2]` requests both model parameter files. |
| Two-preset comparison | Both Hifiasm branches and the level-1 stage summary target appear. |
| Existing-draft QC | Single-draft QUAST/statistics graph builds with HiFi retained. |
| Hi-C scaffolding | Hifiasm/Arima/YaHS branch builds with paired Hi-C, ploidy 2, and the documented `lambda.p2` entry. |
| Old ploidy key | `kmer_ploidy_test_list: [3]` adds no p3 target. |
| Empty filter selection | Returns 0 but requests no files. |
| Custom statistics root | Returns 0 but prints `-q results/...` while writing targets under `custom-results/...`; dry-run success does not clear the path defect. |
| Empty Merqury selection | A graph with that selection prints bare `-m` in the statistics command. |
| Draft-only inputs | Fails at the included MitoHiFi definition with `KeyError: 'hifi'`, even with no mtDNA stage and an empty final-k-mer datatype list. |

A separate invocation of the actual statistics script with bare `-m` returned
exit code 2 and `expected one argument`, confirming that diagnostic. The documented
draft recipe retains the valid HiFi selection. Temporary synthetic fixtures test
graph construction only; they were not assembled or interpreted biologically.

### Executed draft-QC check

The completed Phase 2 yeast FASTA was supplied as one draft, with the real
40,000-read HiFi file available for initialization. The isolated core used
`out_dir: results/`; the project used the exact draft recipe fragment, QC level 1,
and the existing named environments. No assembly jobs were run.

The workflow completed all **8 scheduled jobs**, including draft import, QUAST
5.3.0, parameter/stage statistics gathering, indexing, and length extraction.
The following dry run returned 0 and reported all requested files up to date.

The imported FASTA's SHA-256 matched the source:
`284fc88784dfca6e2be770d1be4955812b141e31467ade7100c0c70151ce30bf`.
Independent checks compared the stage table with the source QUAST report for N50,
L50, largest contig, and the thresholded total length. Observed N50 was 808,003,
L50 was 6, and total length at or above 10,000 bp was 12,264,312 bases. These
describe the same existing assembly; they are not new evidence of completeness
or reference-based structural correctness.

The duplicated summary headers were reproduced, and the output guide's table-view
code ran successfully with its documented draft-QC path substitution. Actual
reports, `comparison-view.txt`, `qc-checks.json`, `execution.log`, and
`up-to-date.log` are retained in `/tmp/assemblybrute-phase4/draft-qc/`.
Other graph records are in `results.json` and per-case `dry-run.log` files;
the bare-argument probe is in `empty-merqury-argument.log`. Temporary artifacts
are not committed.

No new candidate assemblies, BUSCO/Merqury evaluations, Hi-C scaffolding jobs,
or scheduler submissions were executed in Phase 4. The additional limitations
were documented, not patched in workflow code.

## Phase 5: documentation release review

The active Markdown guide and root README were reviewed together against the
recorded evidence. The root README now links to installation, inputs,
configuration, stages, execution, outputs, recipes, and troubleshooting. Its
unfinished configuration section, missing sample-config command, and blanket
tool-status list were replaced with the current guide and explicit validation
scope. The overview and pipeline illustrations remain.

### Final checks

- Current quickstart subset/setup blocks were executed in a fresh isolated work
  directory with compact synthetic FASTQ fixtures. Checks covered final-file
  publication, 40,000-record selection, existing-file protection, private core
  paths, the Bloom-filter override, and environment mapping.
- The exact root README dry-run command built the example graph from that setup.
  Its canonical FASTA/length targets matched the guide, both Hifiasm commands
  contained `-f 0`, and no shared-QC summary was requested at level 0. Synthetic
  fixtures were not assembled or interpreted biologically.
- Active Markdown links and local section anchors, root image paths, code fences,
  shell syntax, Python syntax, YAML fragments, the example YAML, and whitespace
  were checked. Archived manuals were excluded from current-guide link checks.
- Ten public reference links, the dataset download URL, and the repository URL
  returned HTTP 200 to header requests following redirects. The linked BUSCO
  section anchor was checked in its HTML. No dataset was downloaded again.
- All eight archived files matched the original SHA-256 manifest, including
  their relative paths. Existing user changes and `AGENTS.md` were preserved.

The review distinguishes executed examples from graph/source-only checks:

| Coverage | Evidence retained |
| --- | --- |
| HiFi-only primary assembly | Phase 2 executed Hifiasm/Meryl/GenomeScope workflow; independent FASTA/index/length checks and up-to-date rerun. |
| Existing-draft QUAST/statistics | Phase 4 executed eight-job workflow; unchanged imported FASTA, metrics checked against source reports, and up-to-date rerun. |
| Additional recipe/configuration paths | Earlier phase-specific graph checks and diagnostic failures; the release pass does not promote them to executed scientific workflows. |
| New-user example preparation | Fresh-directory execution of the current preparation blocks and README dry-run command. |
| Clean installation and scheduler operation | Still unverified; existing named environments were used, and no jobs were submitted to Slurm. |

The release review did not repeat completed scientific runs, modify workflow
code, create environments, publish a site/release, or change the contributor
guide. The installation page documents the current named-environment requirement;
older automatic-environment advice is not validated by this review.

Final check artifacts are retained locally under `/tmp/assemblybrute-phase5/`,
including `fresh-quickstart-dry-run.log`, `release-checks.json`,
`external-link-results.json`, and HTTP response headers. They are not committed.
Known pipeline limitations remain documented and are outside this documentation
change.

## Test-data example catalog

The [example catalog](../examples/README.md) is the entry point for the
[yeast stage tests](../examples/yeast.md) and [Scam5 dataset](../examples/scam5.md).
The earlier `quickstart.md` URL remains a pointer. The preparation checks below
record an earlier revision; that HiFi-only tutorial is no longer distributed.

Checks for this documentation change used isolated files under
`/tmp/assemblybrute-examples-review/`:

- Executed the revised yeast subset/core preparation blocks with a synthetic
  40,000-record FASTQ outside the checkout. Confirmed the resolved external input
  path, complete subset publication, refusal to overwrite existing subset/core
  files, and refusal to publish a short input as a completed subset.
- Ran the root README's exact Snakemake dry-run command against that fixture with
  existing named environments. It returned zero, selected the primary yeast
  branch, and retained both small-genome `-f 0` options. This checks preparation
  and the job graph; the biological assembly was not repeated.
- Checked all 16 active Markdown pages, local links/anchors, shell/Python/YAML
  block syntax, and the eight archived files' original hashes.
- Matched the then-recorded Scam5 download URLs, MD5 values, byte counts,
  experiment, and sample identity against the ENA file report. Confirmed exactly
  five files from the three requested DNA runs, consistent mate filenames, and
  total compressed size 7,381,143,207 bytes. RNA-seq is excluded.

Scam5 data have not been downloaded or checked by computing file checksums;
the metadata check does not verify downloaded file contents. At that stage, its
scenarios had no DAG or execution evidence; later placeholder DAG checks are
recorded below. The correct source project for the specified
DOI and SRR accessions is `PRJNA1280708`, not `PRJNA1268985`.

Local artifacts are `check.py`, `checks.json`, `external-cache-dry-run.log`, and
the downloaded ENA metadata tables in the temporary directory above. They are
not committed. Workflow code, existing user results, and `AGENTS.md` were preserved.

## Dataset output directories: yeast

The output-layout check used the six standalone projects under
`config/templates/yeast/`, then configured with `input_dir: input/`,
`out_dir: results/yeast/`, and `tmp_dir: results/yeast/tmp/`. The current local
HiFi/Hi-C files and the two existing draft-input links were retained.

Only DAG checks and small path tests were performed for this output-layout change:

- All six projects returned zero from Snakemake 7.32.4 with `--dry-run`,
  `--use-conda`, and eight scheduling cores, in an isolated work directory.
  Every flattened entry in `requested_files.tab` stayed under `results/yeast/`.
- Printed Hifiasm commands use `ln -srf` rather than a fixed number of parent
  traversals. Printed QC-statistics commands read the selected output tree.
- An additional contig DAG used CLI paths for `results/other/` and an absolute
  input directory, confirming that dataset paths are not tied to the yeast name.
- The earlier `Example01.yeasts.contig_.yaml` (subsequently removed), with no explicit project paths,
  still built its DAG using the default `results/` and `results/tmp/` roots.
- Unit tests checked core path defaults versus explicit overrides, retained core
  precedence for other settings, and actual creation/resolution of the three
  Hifiasm correction links for default, nested, and absolute output directories.
  No assembler or QC tool was invoked by those link tests.

The logger creates missing output parents before writing initialization logs.
Project/CLI path overrides now survive loading the shared core. These are local
changes relative to the earlier audit; its historical execution paths remain
part of the original evidence, not namespaced-run results.

Local graph logs and `checks.json` are under
`.snakemake/yeast-layout-check/`. The lightweight regression checks described
above were run before the `tests/` directory was removed at the user's request.

No scientific jobs were run for this change. Runtime, peak memory, and successful
end-to-end execution in the new namespaces remain unverified. Existing outputs
directly under `results/` were not migrated. The earlier HiFi-only assembly
execution remains historical evidence, not a current runnable tutorial.

## Yeast input-root selection

The input-root-selection check used `input_dir: input/yeast/`, retaining
`hifi/fastq/`, `hic/fastq/`, and `draft/fasta/` below that root. Five file symlinks
in this checkout expose the existing HiFi, paired Hi-C, and paired draft inputs
without moving or duplicating their contents.

All six configurations passed dry runs using this layout. An additional contig
dry run set `input_dir` through the CLI to the absolute path of the earlier
`input/` root. The effective configuration and discovered read paths followed
the selected input root; outputs remained under `results/yeast/`.

Only DAG construction was checked. Workflow source files were unchanged, and no
scientific jobs were executed. Local logs and `checks.json` are under
`.snakemake/yeast-input-check/`.

## Test-data root and Scam5 configurations

The yeast tree was subsequently moved to `input/test_data/yeast/`. Its three
relative FASTQ symlinks were adjusted for the additional directory level. The
two draft FASTAs remain regular files and match their pre-move SHA-256 values.
All six yeast configurations passed dry runs with the relocated inputs.

Four complete projects in `config/templates/scam5/` select
`input_dir: input/test_data/scam5/`, `out_dir: results/scam5/`, and
`tmp_dir: results/scam5/tmp/`: `contig.yaml`, `illumina_qc.yaml`,
`hic_scaffolding.yaml`, and `combined.yaml`. Their DAGs passed using five tiny
synthetic FASTQ placeholders with the documented filenames in an isolated work
directory. These fixtures were not placed in the real Scam5 input tree.

Checks confirmed primary `hap0` assembly branches under model `ploidy: 1`,
MboI configuration, unphased YaHS targets for the Hi-C cases, Illumina Merqury
dependencies for the WGS cases, and dataset-specific output paths. The real
`input/test_data/scam5/{hifi,illumina,hic}/fastq/` directories remain empty.

No data were downloaded and no scientific jobs were executed. Workflow source
hashes were unchanged. These checks establish graph construction, not execution,
resource requirements, or biological results for Scam5. Local logs and
`checks.json` are under `.snakemake/test-data-layout-check/`.

## Repository-root test data

On 2026-11-09, `input/test_data/` was moved to `test_data/` at the repository
root. The ten project configurations now use `test_data/yeast/` or
`test_data/scam5/`, while output paths remain `results/yeast/` and
`results/scam5/`. The earlier paths above describe previous checks.

Relative yeast FASTQ links were adjusted and both regular draft FASTAs preserved.
Data files are excluded from Git; `DESCRIPTION` markers are now allowed to
preserve the directory structure without committing sequencing data.
DAG checks use existing yeast inputs and isolated synthetic Scam5 placeholders;
no sequencing data were downloaded and no scientific jobs were executed.

## Conda YAML path resolution

On 2026-11-09, environment selection was centralized in `get_conda_env(tool)`.
Named mode still returns `conda.<tool>.name`; YAML mode returns an absolute path
based on the main `Snakefile` directory. This replaces the fixed `../../../`
prefix in rules, including the shared Singularity environment selector.

All six yeast and four Scam5 configurations passed with Snakemake 7.32.4 in
both modes: 20 DAG/environment-list checks using `--use-conda --dry-run
--list-conda-envs --config use_existing_envs=True` or `False`. These checks read
the selected YAML files and calculate environment hashes without installing
packages. Requested output targets matched between modes for every configuration.

Checks used existing yeast inputs and synthetic Scam5 FASTQ placeholders in an
isolated work directory. Selector assertions also checked all core environment
names/YAML paths and preservation of absolute YAML paths. Rule edits were compared
with pre-change copies to verify that only environment selection changed.
Local checks and logs are under `.snakemake/conda-path-check/`.

No packages were installed, no sequencing data were downloaded, and no scientific
jobs were executed. Dependency solving, fresh environment creation, and tool
execution in automatically created environments remain unverified.

## Unified stage-test configurations

On 2026-11-09, four Scam5 configurations were added so both datasets provide
`read_qc.yaml`, `all_kmer_counters.yaml`, `ploidy_check.yaml`, `contig.yaml`,
`draft_qc.yaml`, and `hic_scaffolding.yaml`. Scam5 retains the additional
`illumina_qc.yaml` and `combined.yaml` scenarios.

The new read-QC test selects HiFi, WGS, and Hi-C; the counter/ploidy tests select
HiFi and WGS. The ploidy test requests models 1–4 and Smudgeplot. Draft QC selects
HiFi plus one primary draft FASTA, with QC level 2 for QUAST and Merqury.

All 14 configurations passed DAG/environment-list checks with Snakemake 7.32.4
in both named and YAML environment modes: 28 checks using `--use-conda --dry-run
--list-conda-envs --config use_existing_envs=True` or `False`. Requested targets
matched between modes and stayed within each dataset's output namespace.
Assertions checked the new read reports, all four k-mer counters, ploidy models,
Smudgeplot targets, and draft-QC outputs.

Checks used existing local reads and yeast drafts. Scam5 draft QC used a tiny
synthetic FASTA in an isolated work directory; no assembly was created under the
real input root. Two `DESCRIPTION` markers preserve its new draft input layout.
No packages were installed, no sequencing data were downloaded, and no scientific
jobs were executed. Workflow source hashes were unchanged. Local checks and logs
are under `.snakemake/scam5-unified-check/`.
