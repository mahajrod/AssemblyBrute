

# Examples

```commandline
snakemake --cores 60  --configfile config/default.yaml --printshellcmds --latency-wait 30   --config mode="assembly" "assembly_mode"="hic_scaffolding" "parameter_set"="normal" "busco_lineage_list"='["vertebrata_odb10","actinopterygii_odb10"]' "data_types"="hifi,hic" "tax_id"=206126 "use_existing_envs"=False --latency-wait 30 --use-conda --rerun-incomplete --res fcs=1 hifiasm=1
```

```commandline
snakemake --cores 70  --configfile config/default.yaml --printshellcmds --latency-wait 60   --config mode="purge_dups" "parameter_set"="big"  "data_types"="hifi,hic "use_existing_envs"=False "skip_busco"=True   --use-conda --rerun-incomplete --res hifiasm=1  fcs=1 fcs_adaptor=1 mem=800000
```