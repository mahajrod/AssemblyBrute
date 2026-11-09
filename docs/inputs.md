# Preparing inputs

Run the pipeline from the repository root. The effective configuration defines
`input_dir`; within that directory, separate reads by datatype and format.
Each project represents one assembly sample: files discovered under a selected
datatype are treated as inputs to that sample, not as separate samples.
Use paths and filenames without spaces or shell metacharacters; several workflow
commands interpolate them directly into shell commands.

Select the input root in the project YAML, for example `input_dir: input/sample/`,
or override it with `--config input_dir=/data/sample/`. All datatype/format paths
below are relative to that root. Both relative and absolute input roots are
accepted; preserve the internal directory structure when changing the root.
Output paths are configured separately.

## Read directories and formats

```text
<input_dir>/
├── hifi/fastq/sample.fastq.gz
├── hic/fastq/library_1.fastq.gz
│             library_2.fastq.gz
├── illumina/fastq/library_1.fastq.gz
│                  library_2.fastq.gz
└── nanopore/fastq/sample.fastq.gz
```

Select datatypes as a comma-separated string without spaces, for example
`input_datatypes: hifi,hic`. Include only datatypes you have prepared; selecting
one with no recognized files stops initialization. Unselected read datatypes are
excluded from discovery.

The declarations in [`config/datatypes.yaml`](../config/datatypes.yaml) are:

| Datatypes | Layout |
| --- | --- |
| `hifi`, `lqccs` | Single-end PacBio reads. |
| `nanopore`, `simplex`, `duplex`, `ultralongnano`, `adaptivenano` | Single-end Nanopore read categories. |
| `hic`, `illumina` | Paired reads. |

Each read datatype declares these format directories:

| Directory | Recognized extensions |
| --- | --- |
| `fastq/` | `.fastq.gz` or `.fq.gz` |
| `fasta/` | `.fasta.gz` or `.fa.gz` |
| `mapped_bam/`, `unmapped_bam/` | `.bam` |

Use **one format and one extension per selected datatype**. Mixing `.fq.gz` and
`.fastq.gz`, or supplying both FASTA and FASTQ for the same datatype, is rejected.
Discovery searches directly inside these directories; it does not recursively
collect nested library directories.

## Paired-read names

Use a common prefix and a single-character mate identifier:

```text
libraryA_1.fastq.gz
libraryA_2.fastq.gz
libraryB_1.fastq.gz
libraryB_2.fastq.gz
```

The parser sorts filenames, pairs consecutive entries, rejects an odd number of
files, and requires mate filenames to differ at no more than one character.
Names such as `libraryA_R1.fastq.gz`/`libraryA_R2.fastq.gz` satisfy that naming
check. Consistent suffixes must also be shared across libraries in a datatype.

## Existing draft assemblies

To start from assemblies, select `draft` among the inputs and `draft_qc` as the
assembly starting stage. `draft_qc` and `contig` are mutually exclusive.
Current included MitoHiFi definitions still access HiFi during parsing, so retain
HiFi input for the [tested draft-QC recipe](recipes.md#5-quast-evaluation-of-an-existing-draft).
A draft-only selection currently fails even when mtDNA is not selected.

```text
<input_dir>/draft/fasta/
├── sample.hap1.fasta
└── sample.hap2.fasta
```

Set `ploidy: 2` for this two-file layout. The parser requires exactly `ploidy`
files ending in `.fasta`; for multiple haplotypes, suffixes must be `hap1.fasta`,
`hap2.fasta`, and so on. A single `.fasta` with `ploidy: 1` is assigned `hap0`.
Compressed draft FASTAs are not recognized by this interface. Set
`assembly_qc_level.draft_qc` explicitly: level 0 imports the assembly but requests
no shared assembly-QC analyses.

## Reference genomes

References use their own directory layout:

```text
<input_dir>/reference/reference_name/
├── reference_name.fasta
├── reference_name.whitelist
├── reference_name.orderlist
└── reference_name.syn
```

Provide one file of each type. The whitelist identifies selected sequence IDs;
the order list specifies their ordering; the synonym file supplies tab-separated
ID mappings. Check that IDs agree with the FASTA headers. Provide `.syn` explicitly
because the current fallback that tries to create it needs review. Empty or
incomplete reference directories can fail during initialization.

Reference-based QC and the separate `ref_scaffolding` stage have different
implementation status. The latter is incomplete; see the
[tool support matrix](stages.md#tool-support).
