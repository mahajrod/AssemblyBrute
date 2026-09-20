# AssemblyBrute - Pipeline to "brute force" and evaluate de novo genome assemblies

It allows to generate multiple assemblies (using different tools and different parameters for them) from a single command
It is based on VGP-pipeline and Rapid curation but with multiple additions for QC, evaluation and curation.


# Dependencies

If you wish to run it using conda via snakemake, then you will need:
  - conda or mamba
  - snakemake
  - FCS database and FCS_GX singularity container      # optional
  - FCS_adapter singularity container                  # optional
  - Kraken databases                                   # optional
  - RapidCuration singularity containers               # this dependency will be excluded soon

# Stages of the pipeline

**Preprocessing stages**:
  - _raw_read_qc_      # comment this stage if you wish to skip quality control of the raw data
  - _raw_kmer_qc_      # runs kmer counting and genome size estimation from raw reads. Usually you don't need it
  - _filter_reads_     
  - _filtered_read_qc_
  - _kmer_qc_
  - _ploidy_check_
  - _mtdna_
  - _read_contamination_scan_

**Main stages**:
  - _draft_qc_          # SELECT EITHER 'contig' or 'draft_qc'. This stages actually initiate assembly process
  - _contig_            # SELECT EITHER 'contig' or 'draft_qc'. This stages actually initiate assembly process
  - _polishing_         # OPTIONAL. MUST NOT PRECEDE 'contig' or 'draft_qc' stages
  - _dedup_             # OPTIONAL. MUST NOT PRECEDE 'contig' or 'draft_qc' stages # TODO: do more testing for hapsolo
  - _hic_alignment_     # OPTIONAL. MUST NOT PRECEDE 'contig' or 'draft_qc' stages
  - _hic_scaffolding_   # OPTIONAL. MUST FOLLOW hic_alignment stage
  - _ref_scaffolding_   # OPTIONAL. MUST NOT PRECEDE 'contig' or 'draft_qc' stages # TODO: do more testing , at moment it requires telomere detection
  - _gap_closing_       # OPTIONAL. MUST NOT PRECEDE 'contig' or 'draft_qc' stages

Attached stages:
  - _read_phasing_ # this is an attached stage, it is attached to the stage set by "phasing_stage" parameter. However, if phased reads are required for any other stage, the corresponding rules will be called automatically. Uncomment it only if you need to get phased reads for all suitable datatypes.


# Implemented tools
**Contig assembly**:

| Assembler |                        Datatypes                        | Status |
|:---------:|:-------------------------------------------------------:|:------:|
|  hifiasm  | hifi/nanopore (+ ultra long reads)<sup>*</sup> (+ Hi-C) |   v    |
|   flye2   |                      hifi/nanopore                      |   v    | 
|  verkko   |           hifi/nanopore<sup>**</sup> (+ Hi-C)           |   v    |

<sup>*</sup> - brackets indicate an optional datatype  
<sup>**</sup> - can use both hifi and nanopore simultaneously 

**Polishing**:
 - NextPolish2

**Deduplication**:
 - purge_dups
 - hapsolo
 - combination of purge_dups and hapsolo

**Hi-C Aligners**:
 - Arima mapping pipeline
 - Pairtools
 - Juicer

**HiC scaffolding**:
 - YaHS
 - 3D-DNA

**Gap closing**:
 - SAMBA

**Reference-based scaffolding**:
 - RagTag

# Usage
I. Clone this repository
```commandline
git clone https://github.com/mahajrod/AssemblyBrute 

```

II. Place you fastq or fasta file in corresponding folders in the input directory:
```commandline
AssemblyBrute/
    input/
        hic/
            fastq/
            fasta/
        hifi/
            fastq/
            fasta/
        nanopore/
            fastq/
            fasta/
        simplex/
            fastq/
            fasta/
        duplex/
            fastq/
            fasta/
        ultralongnano/
            fastq/
            fasta/
        adaptivenano/
            fastq/
            fasta/     
        illumina/
            fastq/
            fasta/
        draft/
            fasta/
        reference/
            &lt; reference name &gt;/
            
```

III. Adjust config files.
```commandline
**To be written**
```
Some of the options (all nonested options from default.yaml) could also be set via command line. See examples before

IV. Run pipeline directly

# Examples

```commandline
snakemake --profile profile/slurm/ \
          --configfile config/data/fSynTyp1.2.yaml --printshellcmds \
          --latency-wait 120 --use-conda  \
          --res fcx=1 fcs_adaptor=1 \
          --rerun-incomplete -k
```

Running main stages one by one is recommended.

