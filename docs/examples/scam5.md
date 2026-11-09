# Scam5: HiFi, Illumina WGS, and Hi-C

The eight standalone configurations in
[`config/templates/scam5/`](../../config/templates/scam5/) read from
**`test_data/scam5/`** and write to **`results/scam5/`**.
Run commands from the repository root with Snakemake 7 and the named tool
environments prepared as described in [installation](../installation.md).

**Status: all eight DAGs checked using local reads in an isolated work directory;
draft QC used a synthetic FASTA.** No scientific jobs were executed during these
checks. Prepare the files listed below before running a selected test.

## Source and identity

Use *Sinodiscula camellicola* strain Scam5 from
[Shi et al., Scientific Data (2026)](https://doi.org/10.1038/s41597-026-07115-4),
BioProject [PRJNA1280708](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1280708),
BioSample `SAMN50197380`. ENA currently labels the organism
`Sinodiscula camellicola (nom. inval.)`, taxon `3114500`.

BioProject [PRJNA1268985](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1268985)
belongs to a different dataset, *S. theae-sinensis*
strain Sthe1; do not mix its reads into this sample. The accessions below identify
SRA **runs**; their experiments are listed separately.

| Datatype | Run | Experiment | Layout | Compressed download |
| --- | --- | --- | --- | --- |
| `hifi` | SRR35179663 | SRX30284538 | One HiFi FASTQ | 1.54 GB |
| `illumina` | SRR35179662 | SRX30284539 | Paired FASTQ | 3.38 GB |
| `hic` | SRR35179660 | SRX30284541 | Paired FASTQ | 2.47 GB |

Download URLs, exact byte counts, and MD5 checksums for the five FASTQ files are
available from ENA. Metadata was checked on 2026-11-07 using the
[ENA file report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1280708&result=read_run&fields=run_accession,experiment_accession,sample_accession,scientific_name,tax_id,instrument_platform,instrument_model,library_strategy,library_layout,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv).
RNA-seq run `SRR35179661` is deliberately excluded; select explicit runs rather
than downloading the whole BioProject.

All five files total **7,381,143,207 bytes**, about 7.38 GB / 6.87 GiB compressed.
Allow additional space for temporary files, indexes, and results; peak disk/RAM
requirements have not been measured. The paper's Gb figures describe sequence
bases, not download sizes; its processed yields need not equal archived raw yields.

## Storage and preparation

Prepare a dataset-specific input directory; see the
[input storage guidance](README.md#keep-sequencing-data-outside-git):

```text
test_data/scam5/
├── hifi/fastq/SRR35179663_subreads.fastq.gz
├── illumina/fastq/SRR35179662_1.fastq.gz
│                 SRR35179662_2.fastq.gz
└── hic/fastq/SRR35179660_1.fastq.gz
              SRR35179660_2.fastq.gz
```

Download to temporary filenames, verify byte counts, MD5 and `gzip -t`, then
publish completed files. Reads stored elsewhere can be linked into these paths.
All projects select these paths:

```yaml
main_config_file: config/core.yaml
input_dir: test_data/scam5/
out_dir: results/scam5/
tmp_dir: results/scam5/tmp/
```

Only the datatypes selected by a project are discovered. The HiFi-only case can
run once its single HiFi file is available; the other cases require their WGS
and/or paired Hi-C files. To inspect another data location, override `input_dir`
with `--config input_dir=/path/to/scam5/`, retaining the datatype/format layout.

The PacBio archive name contains `subreads`; the paper describes HiFi sequencing.
Keep it under `hifi/`, following that provenance. Confirm paired FASTQ record
counts/identifiers as well as filename pairing before execution.

## Configurations and modeling settings

| Configuration | `input_datatypes` | Requested analyses |
| --- | --- | --- |
| `read_qc.yaml` | `hifi,illumina,hic` | Raw/filtered read QC, filtering, and default k-mer QC without assembly. |
| `all_kmer_counters.yaml` | `hifi,illumina` | Meryl, Yak, Jellyfish, and FastK with genome-size reports. |
| `ploidy_check.yaml` | `hifi,illumina` | GenomeScope models 1–4 and Smudgeplot. |
| `contig.yaml` | `hifi` | HiFi read QC/preparation, primary Hifiasm assembly, QUAST and Merqury. |
| `draft_qc.yaml` | `hifi,draft` | Import and QUAST/Merqury evaluation of one existing draft FASTA. |
| `illumina_qc.yaml` | `hifi,illumina` | Raw/filtered read QC, k-mer models, assembly and evaluation using both HiFi and WGS evidence. |
| `hic_scaffolding.yaml` | `hifi,hic` | HiFi assembly, Arima alignment and unphased YaHS scaffolding, with assembly QC. |
| `combined.yaml` | `hifi,illumina,hic` | Read/k-mer QC, assembly, Hi-C scaffolding and HiFi/WGS-based evaluation. |

The projects use `resource_profile: high`, `parameter_set: micro`, genome prefix
`scam5`, and **`ploidy: 1` as a primary-assembly model**. This setting does not
establish the biological ploidy of the sample. Both scaffolding scenarios keep that
model and select YaHS `unphased`; Hi-C is used for scaffolding, not diploid phasing.
The paper's **MboI** setting is configured as `hic_enzyme_set: MboI`.

The six common stage tests match the names used for the other dataset;
`illumina_qc.yaml` and `combined.yaml` add WGS-specific assembly scenarios.

Contig, draft, and scaffolding QC level 2 requests QUAST and Merqury. The WGS scenarios add
Illumina k-mer evidence automatically; assembly coverage estimation remains based
on HiFi reads. BUSCO, FCS/adaptor screening and mitochondrial analysis are skipped;
microchromosome detection is disabled. No polishing stage is selected.
See [recipes](../recipes.md) for other stage combinations and their limitations.

## Prepare draft-QC input

After executing `contig.yaml`, copy its primary assembly:

```bash
cp -n results/scam5/contig/hifiasm_l1_no_ul@p1/scam5.contig.hap0.fasta test_data/scam5/draft/fasta/scam5.hap0.fasta
```

Keep exactly one uncompressed FASTA in `test_data/scam5/draft/fasta/`, matching
`ploidy: 1`, and retain the HiFi reads. A contig dry run does not create this FASTA.

## Check a DAG

After preparing the selected input files, inspect the graph before execution:

```bash
snakemake --configfile config/templates/scam5/contig.yaml --use-conda -c 8 \
  --dry-run --printshellcmds
```

Use the same flags with any configuration. An empty input directory
will fail initialization; the pipeline does not download missing data.

## Check individual stages

Read QC:

```bash
snakemake --configfile config/templates/scam5/read_qc.yaml --use-conda -c 8 --dry-run
```

All k-mer counters:

```bash
snakemake --configfile config/templates/scam5/all_kmer_counters.yaml --use-conda -c 8 --dry-run
```

Ploidy models and Smudgeplot:

```bash
snakemake --configfile config/templates/scam5/ploidy_check.yaml --use-conda -c 8 --dry-run
```

Existing draft QC:

```bash
snakemake --configfile config/templates/scam5/draft_qc.yaml --use-conda -c 8 --dry-run
```

Remove `--dry-run` to execute a selected stage test. These read-analysis scenarios
write under `qc/`, `data/`, and `kmer/`; draft QC writes under `draft_qc/`, all
within `results/scam5/`.

## Run one scenario

The following commands execute jobs once their input files and environments are
available. Run one scenario at a time; completed compatible intermediates are
shared under `results/scam5/`.

HiFi primary assembly:

```bash
snakemake --configfile config/templates/scam5/contig.yaml --use-conda -c 8 --rerun-incomplete
```

HiFi plus Illumina WGS evaluation:

```bash
snakemake --configfile config/templates/scam5/illumina_qc.yaml --use-conda -c 8 --rerun-incomplete
```

HiFi plus Hi-C scaffolding:

```bash
snakemake --configfile config/templates/scam5/hic_scaffolding.yaml --use-conda -c 8 --rerun-incomplete
```

Combined HiFi, WGS and Hi-C workflow:

```bash
snakemake --configfile config/templates/scam5/combined.yaml --use-conda -c 8 --rerun-incomplete
```

## Expected outputs

The primary FASTA is
`results/scam5/contig/hifiasm_l1_no_ul@p1/scam5.contig.hap0.fasta`.
Both scaffolding scenarios request
`results/scam5/hic_scaffolding/hifiasm_l1_no_ul@p1..arima_def..yahs_unphased/scam5.hic_scaffolding.hap0.fasta`.
Stage statistics are `results/scam5/contig/scam5.contig.stage_stats` and, when
selected, `results/scam5/hic_scaffolding/scam5.hic_scaffolding.stage_stats` or
`results/scam5/draft_qc/scam5.draft_qc.stage_stats`.
Read reports, k-mer data, effective configuration and logs are under `qc/`,
`kmer/`, `config/`, and `log/` in the same dataset output root. Configuration
snapshots and summary tables reflect the most recent scenario selection.

Run the HiFi case first, then add WGS evaluation, then Hi-C scaffolding, and finally
the combined case. For each, retain the effective config, input identities, tool
versions, logs, output checks, observed runtime/resources, and an up-to-date rerun.
Use the full dataset for biological evaluation; small paired subsets are useful
for technical checks but require their own reproducible selection recipe.

The paper reports a 41.79 Mb assembly on nine chromosomes and deposits reference
`GCA_053477585.1`. Use these as comparison context. Different software, filtering,
and manual curation mean those values are not exact regression thresholds.
