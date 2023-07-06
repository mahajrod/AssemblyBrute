# AssemblyBrute - Pipeline to "brute force" and evaluate de novo genome assemblies

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

# Input datatypes and implemented stages of the pipeline

| Datatypes | read QC | read filtration | contamination check | ec | contig | purge_dups | hic qc | hic scaffolding | curation | gap_closing | qc |
|:---------:|:-------:|:---------------:|:-------------------:|:--:|:------:|:----------:|:------:|:---------------:|:--------:|:-----------:|:--:|
| hifi + hic | v | v | v | v | v | v | v | v | v |  | v |
| hifi + hic (no phasing) | v | v | v | v | v | v | v | v | v | v |
| hifi | v | v | v | v | v | v | NA | NA |  |  | v |
| clr | v | v | v |  |  |  | NA | NA |  |  | |
| nanopore | v | v | v |  |  |  | NA | NA |  |  |  |
| illumina | v | v | v |  |  |  | NA | NA |  |  |  |
| nanopore + illumina | v | v | v |  |  |  | NA | NA |  |  |  |
| nanopore + illumina + hic | v | v | v |  |  |  |  |  |  |  |  |
| clr + illumina | v | v | v |  |  |  | NA | NA |  |  |  |
| clr + illumina + hic | v | v | v |  |  |  |  |  |  |  |  |
| assembly | NA | NA | v | NA | NA | NA | NA | NA | NA | v | v |
| assembly + hic |  |  |  | |  |  |  |  |  |  |
| assembly + nanopore + hic |  |  |  |  | NA |  |  |  |  |  |
| assembly + hifi + hic |  |  |  |  | NA |  |  |  |  | v | v |

# Implemented assemblers
|  Assembler  |       Datatypes       | Status |
|:-----------:|:---------------------:|:------:|
|   hifiasm   |      hifi + hic       |   v    |
|   hifiasm   |         hifi          |   v    | 
|   hifiasm   | hifi + nanopore + hic |        | 
|   hifiasm   |    hifi + nanopore    |   ?    | 
|   hicanu    |       nanopore        |        |
|    flye2    |       nanopore        |        | 
|     IPA     |       nanopore        |        | 
| SmartDenovo |       nanopore        |        | 
|   MaSurCa   |  nanopore + Illumina  |        | 

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

IV. Run pipeline directly or via wrapper script (**Not written yet**). See examples below

# Examples

```commandline
snakemake --cores 60  --configfile config/default.yaml --printshellcmds --latency-wait 30   --config mode="assembly" "assembly_mode"="hic_scaffolding" "parameter_set"="normal" "busco_lineage_list"='["vertebrata_odb10","actinopterygii_odb10"]' "data_types"="hifi,hic" "tax_id"=206126 "use_existing_envs"=False --latency-wait 30 --use-conda --rerun-incomplete --res fcs=1 fcs_adaptor=1 mem=800000 kmer_counter=1  telosif=1
```

```commandline
snakemake --cores 70  --configfile config/default.yaml --printshellcmds --latency-wait 60   --config mode="purge_dups" "parameter_set"="big"  "data_types"="hifi,hic "use_existing_envs"=False "skip_busco"=True   --use-conda --rerun-incomplete --res fcs=1 fcs_adaptor=1 mem=800000 kmer_counter=1  telosif=1
```

```commandline
snakemake --cores 140  --configfile config/test1.yaml --printshellcmds --latency-wait 60   --config mode="assembly" assembly_mode="curation" "parameter_set"="micro"  "data_types"="hifi,hic" "tax_id"=4932  "use_existing_envs"=False "skip_busco"=True  "final_kmer_datatype"="hifi" "busco_lineage_list"='["saccharomycetes_odb10"]' "skip_gcp"=True  --use-conda --rerun-incomplete --res   fcs=1 fcs_adaptor=1 mem=800000 kmer_counter=1  telosif=1
```

```commandline
snakemake --cores 60 --configfile config/default.yaml --printshellcmds --latency-wait 60   --config mode="qc" "parameter_set"="large"  "data_types"="illumina" "tax_id"=30532  "use_existing_envs"=False "skip_busco"=True "skip_kraken"=True "skip_gcp"=True  "final_kmer_datatype"="illumina"   --use-conda --rerun-incomplete --res hifiasm=1  fcs=1 fcs_adaptor=1 mem=800000 kmer_counter=1
```

```commandline
#shark
git pull; snakemake --cores 80  --configfile config/somniosus.yaml --printshellcmds --latency-wait 60   --config mode="assembly" assembly_mode="contig" "parameter_set"="giant"  "data_types"="hifi,hic" "tax_id"=191813  "use_existing_envs"=False "skip_busco"=True  "final_kmer_datatype"="hifi" "busco_lineage_list"='["vertebrata_odb10"]' "skip_gcp"=True  --use-conda --rerun-incomplete --res   fcs=1 fcs_adaptor=1 mem=800000 kmer_counter=1  telosif=1
```

```commandline
cd /maps/projects/codon_0000/people/svc-rapunzl-smrtanl/yggdrasil/assembly/fish/spinachia_spinachia_2_finalization
git pull; snakemake --cores 35  --configfile config/finalization.yaml --printshellcmds --latency-wait 60   --config mode="finalization" finalization_mode="curation" "parameter_set"="small"  "data_types"="hifi,hic" "use_existing_envs"=False "skip_fcs"=True "skip_fcs_adaptor"=True "skip_read_qc"=True "skip_filter_reads"=True  tax_id=206126 "busco_lineage_list"='["vertebrata_odb10", "actinopterygii_odb10"]' "skip_busco"=False "skip_higlass"=False   --use-conda --rerun-incomplete --res hifiasm=1  fcs=1 fcs_adaptor=1 telosif=1 mem=700000
```

```commandline
git pull; snakemake --cores 35  --configfile config/data/syngnathus_typhle.yaml --printshellcmds --latency-wait 60   --config mode="assembly" assembly_mode="curation" "parameter_set"="small"  "data_types"="hifi,hic" "use_existing_envs"=False "skip_fcs"=False "skip_fcs_adaptor"=False  tax_id=161592 "busco_lineage_list"='["vertebrata_odb10"]' "skip_busco"=False "skip_higlass"=True   --use-conda --rerun-incomplete --res fcs=1 fcs_adaptor=1 telosif=1 mem=700000 

```


```commandline
#Stage by stage (recommended)
#qc and filtering

#contig assembly

#purge dups

#purge_dups with calcuts thresholds adjusted, skip if the thresholds correspond to plots

# hic scaffolding

# curation

# finalization




```