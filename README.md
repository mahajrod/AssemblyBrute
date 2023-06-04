# AssemblyBrute - Pipeline to "bruteforce" de novo genome assembly and evaluation

It allows to generate multiple assemblies (with different parameters and different tools) from one command
It is based on VGP-pipeline (bionano scaffolding was not included) and Rapid curation but with multiple additions for QC and evaluation

# Dependencies

If you wish to run it using conda via snakemake, then you will need:
    - conda
    - mamba
    - snakemake
    - FCS database and FCS_GX singularity container      # optional
    - FCS_adapter singularity container                  # optional
    - Kraken databases                                   # optional
    - RapidCuration singularity containers               # this dependency will be excluded soon

# Usage
I. Clone this repository
```commandline
git clone https://github.com/mahajrod/AssemblyBrute 

```

II. Place you fastqs in corresponding folders in the input directory:
```commandline
AssemblyBrute/
    input/
        hic/
            fastq/
        hifi/
            fastq/
        nanopore/
            fastq/
        illumina/
            fastq/
```

III. Modify config files. I recommend to copy *default.yaml* and do all modifications in this copy.
```commandline
config/
    default.yaml   <----- modify this file, add paths to databases, set tax_id, ploidy, etc
    core.yaml      <----- modify this file only if you know what you are doing. In most case you don't need it
```
Some of the options (all nonested options from default.yaml) could also be set via command line. See examples before
IV.

# Examples

```commandline
snakemake --cores 60  --configfile config/default.yaml --printshellcmds --latency-wait 30   --config mode="assembly" "assembly_mode"="hic_scaffolding" "parameter_set"="normal" "busco_lineage_list"='["vertebrata_odb10","actinopterygii_odb10"]' "data_types"="hifi,hic" "tax_id"=206126 "use_existing_envs"=False --latency-wait 30 --use-conda --rerun-incomplete --res fcs=1 hifiasm=1
```

```commandline
snakemake --cores 70  --configfile config/default.yaml --printshellcmds --latency-wait 60   --config mode="purge_dups" "parameter_set"="big"  "data_types"="hifi,hic "use_existing_envs"=False "skip_busco"=True   --use-conda --rerun-incomplete --res hifiasm=1  fcs=1 fcs_adaptor=1 mem=800000
```