# Test-data examples

Choose a dataset for source reads, preparation, complete configurations, commands,
and expected outputs. Start with [installation](../installation.md).

| Example | Inputs and purpose | Status |
| --- | --- | --- |
| [Yeast stage tests](yeast.md) | Existing HiFi + Hi-C fixtures; six separate configurations, outputs under `results/yeast/`. | All six DAGs checked; namespaced workflows not executed. |
| [Scam5](scam5.md) | Six common stage tests plus two WGS assembly scenarios, outputs under `results/scam5/`. About 7.38 GB compressed for all five read files. | All eight DAGs checked with local reads; draft QC used a synthetic FASTA. No scientific jobs executed during validation. |

Both datasets provide `read_qc.yaml`, `all_kmer_counters.yaml`,
`ploidy_check.yaml`, `contig.yaml`, `draft_qc.yaml`, and `hic_scaffolding.yaml`.
Scam5 also provides `illumina_qc.yaml` and `combined.yaml`.

The six yeast stage tests and Scam5 are opt-in integration datasets;
neither is downloaded by ordinary repository checks. Other input technologies
will need their own datasets and validation records.

## Keep sequencing data outside Git

Commit descriptions, small project configurations, and measured
validation summaries. Keep FASTQ/BAM files, subsets, and generated assemblies out
of commits. Public archives host the source data; prepare each dataset once and
reuse its reads across scenarios.
Each test-data directory contains a tracked `DESCRIPTION` file to preserve its
structure on GitHub. Other files and symlinks under `test_data/` remain ignored.
No Git LFS or additional data-versioning service is needed for these public
archive-backed examples: record run accessions and file identities in the dataset
description.

Set paths directly in the project configuration. A new dataset can use:

```yaml
input_dir: test_data/scam5/
out_dir: results/scam5/
tmp_dir: results/scam5/tmp/
```

The current yeast templates use `test_data/yeast/` and write to
`results/yeast/`. For other datasets, select a separate input directory or point
`input_dir` to reads stored elsewhere. Explicit project/CLI path settings override
core defaults. Set `out_dir` and `tmp_dir` together for each dataset.

The test-data root contains `yeast/` and `scam5/`, each retaining the usual
datatype/format layout. Supply the DNA FASTQ files listed in each guide;
draft-QC scenarios also require prepared FASTAs under `draft/fasta/`.

Keep original downloads separate from any selected input tree or derived subset.
Record run accessions, exact URLs, byte sizes, checksums, and subset method/seed.
Validate complete downloads before making them visible to input discovery.
For paired reads, preserve matching mates when subsetting.

## Separate scenarios and validation levels

Reuse prepared reads across scenarios. Select a dataset output root such as
`results/yeast/` or `results/scam5/`, with its temporary directory inside that
root. Hifiasm links and QC-statistics collection support nested output paths.
Tests for one dataset can share compatible intermediates. Use separate output
roots for independent runs; do not run multiple configurations concurrently into
the same output tree.

**Metadata checked** covers source descriptions; **dry-run checked** confirms DAG
construction, not tool execution, resource requirements, or biological quality.
**Executed** covers only the recorded workflow and settings. Fresh Conda solves
and scheduler execution remain unverified. Published metrics provide comparison
context, not exact regression thresholds.
See the [validation record](../development/validation.md) for completed checks.
