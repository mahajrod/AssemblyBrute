# Outputs

`out_dir` in the effective configuration defines the output root. Inspect
`requested_files.tab` to see what the current configuration requests; optional
files depend on stages, QC levels, presets, datatypes, and skip flags.

## Directory map

| Location under the output root | Contents |
| --- | --- |
| `config/config.final.yaml`, `config/final_data.yaml` | Processed configuration and discovered-input snapshots. Preserve the original project/core too; snapshots are not verified replacement project files. |
| `requested_files.tab` | Requested output paths, not a completion manifest. Initialization writes it even during dry runs. |
| `initialization.*.log`, `log/` | Initialization messages, tool logs, scheduler stdout/stderr, and benchmark records. Some rules place logs inside branch directories. |
| `data/<datatype>/{raw,filtered,final}/` | Normalized read links and applicable filtering outputs. `final` can link to raw reads when filtering is disabled. |
| `qc/` | Read-QC reports, grouped by tool, datatype, and read state. |
| `kmer/<datatype>/<state>/` | Counter databases/histograms, genome-size models, and ploidy-check outputs. |
| `error_correction/` | Hifiasm correction/overlap caches shared by compatible option groups. |
| `<assembly_stage>/<branch>/` | Assemblies, indexes, length tables, tool-specific intermediates, and enabled assembly QC. |
| `<assembly_stage>/<prefix>.<stage>.stage_stats` | Consolidated statistics when shared QC level is at least 1 and statistics gathering succeeds. |
| `download/`, `mtDNA/`, `ext_data/`, `ext_qc/` | Applicable downloaded resources, mitochondrial results, and external-track processing. These are not present in every run. |

## Identify the assembly to use

Canonical assembly paths follow:

```text
<out_dir>/<stage>/<branch>/<genome_prefix>.<stage>.<haplotype>.fasta
```

For example, a two-haplotype contig branch with `out_dir: results/sample/` and
`genome_prefix: sample` has this layout:

```text
results/sample/contig/hifiasm_l1_no_ul@p2/
├── sample.contig.hap1.fasta
├── sample.contig.hap1.fasta.fai
├── sample.contig.hap1.len
├── sample.contig.hap2.fasta
├── sample.contig.hap2.fasta.fai
├── sample.contig.hap2.len
├── assembly_qc/
└── log/
```

`hap0` is the single primary representation; `hap1`, `hap2`, etc. identify
multiple assembled haplotypes. `combined`/`reordered` files are additional
representations, not independently inferred haplotypes. Branch labels preserve
the tool/preset chain; see [stages](stages.md#understand-branch-labels-and-haplotypes).

The `.fai` indexes the FASTA. The `.len` file is a tab-separated sequence-ID/length
table sorted by descending length. GFA files describe assembly graphs; `.lencov`
and `.cov` supply graph-derived length/coverage information.

The canonical `.fasta` name does not certify contamination removal: the skipped
FCS path copies the unfiltered FASTA and creates an empty contaminant file.
Check `skip_fcs`, executed commands, and reports before interpreting that file.

## Read reports and model outputs

- FastQC ZIPs and applicable MultiQC HTML reports live under
  `qc/fastqc/<datatype>/<state>/` and `qc/multiqc/<datatype>/<state>/`.
- Long-read NanoPlot outputs include per-file `*.NanoStats.txt` and aggregated
  `<datatype>.<state>.NanoStats.tsv` under `qc/nanoplot/`. The HiFi recipe does
  not imply that FastQC/MultiQC are requested for HiFi.
- GenomeScope results include model/summary files and plots under `genomescope/`.
  The pipeline also extracts `*.genomescope.parameters` with `Genome size` and
  `Lambda` fields. Read the model and log alongside those extracted values.
- Ploidy-check filenames include `.p1.`, `.p2.`, etc. Additional model requests
  do not change `config.ploidy` or select a biological ploidy automatically.

## Assembly-QC reports

| Report family | Where to look and what to read |
| --- | --- |
| QUAST | `<branch>/assembly_qc/quast/<assembly_prefix>/report.tsv` and companion reports. Read metric labels and minimum-length thresholds. The pipeline's rule supplies no reference genome to QUAST, so reference-based misassembly/alignment metrics are not established by that invocation. |
| Merqury | `<branch>/assembly_qc/merqury/<datatype>/*.qv`, `*.completeness.stats`, and plots. |
| BUSCO | `<branch>/assembly_qc/busco5/`: lineage-specific summaries, full tables, archives, and intersection tables. The current shared-QC target builder requests the `busco5` paths even though separate BUSCO6 rules also exist. |
| Coverage, telomeres, repeats | Track files and plots under branch-specific QC/assembly directories. |
| Hi-C scaffolding/visualization | Scaffold FASTAs and AGP files; Pretext, `.hic`, and other contact-map files depend on QC settings. |

For definitions and interpretation, consult the primary
[QUAST manual](https://github.com/ablab/quast/blob/master/manual.html),
[Merqury documentation](https://github.com/marbl/merqury), and
[BUSCO guide](https://busco.ezlab.org/busco_userguide.html#interpreting-the-results).

The current Merqury integration passes either `hap0` or the first haplotype plus
`hap2`; it does not pass additional polyploid haplotypes. Inspect actual inputs
before treating its report as an evaluation of a complete polyploid assembly.

## Summary tables

Each branch can produce
`<branch>/assembly_qc/<prefix>.<stage>.parameter_stats`; the stage table combines
those rows and identifies `stage`, `parameters`, and `assembly_prefix`.

`gather_stats_per_stage_parameter` reads
`<out_dir>/<stage>/<branch>/assembly_qc/`, including nested output roots.

The gathering script selects `Total length (>= 10000 bp)` and repeats some
column labels, including `Largest contig`. Pandas can read duplicates with `.1`
suffixes. Select columns by their exact names and inspect the source `report.tsv`;
do not assume that the reported thresholded total includes every sequence.
At QC level 0, these summary tables are not requested.
