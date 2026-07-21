rule hapsolo: # hapsolo renames scaffolds
    input:
        input_fasta="{fasta_dir}/{fasta_prefix}.fasta",
        busco_db_dir=config["out_dir"] / "download/busco5/lineages/{busco_lineage}",
        log_dir=ancient("{fasta_dir}/log"),
    output:
        purged_fasta="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.hapsolo.purged.fasta",
        duplication_fasta="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.hapsolo.dups.fasta",
        name_mappings="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/contigs/name_mapping.tsv",
        syn_file="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.hapsolo.syn"
    params:
        iterations_per_thread=parameters["tool_options"]["assembly_qc"]["hapsolo"]["iterations_per_thread"]
    log:
        log="{fasta_dir}/log/hapsolo.{fasta_prefix}.{busco_lineage}.log",
        cluster_log="{fasta_dir}/log/hapsolo.{fasta_prefix}.{busco_lineage}.cluster.log",
        cluster_err="{fasta_dir}/log/hapsolo.{fasta_prefix}.{busco_lineage}.cluster.err"
    benchmark:
        "{fasta_dir}/log/hapsolo.{fasta_prefix}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["hapsolo"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hapsolo"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("hapsolo"),
        cpus=get_threads(parameters["threads"]["hapsolo"], "cpu"),
        time=parameters["time"]["hapsolo"],
        mem=parameters["memory_mb"]["hapsolo"]
    threads: parameters["threads"]["hapsolo"]

    shell:
        " BUSCO_DIR=`realpath {input.busco_db_dir}`; "
        " LOG=`realpath {log.log}`; "
        " WORKDIR=`dirname {output.purged_fasta}`; "
        " mkdir -p ${{WORKDIR}}; "
        " cp -f {input.input_fasta} ${{WORKDIR}}; "
        " INPUT_FASTA=`basename {input.input_fasta}`; "
        " INPUT_FASTA_PREFIX=${{INPUT_FASTA%.fasta}}; "
        " OUTPUT_FASTA=`basename {output.purged_fasta}`; "
        " OUTPUT_DUPLICATION_FASTA=`basename {output.duplication_fasta}`; "
        " cd ${{WORKDIR}}; "
        " hapsolo_cli.py preprocess -i ${{INPUT_FASTA}} > ${{LOG}} 2>&1;  "
        " hapsolo_cli.py align -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta >> ${{LOG}} 2>&1 ; "
        " hapsolo_cli.py search -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta -l ${{BUSCO_DIR}} -o ortholog_output >> ${{LOG}} 2>&1; "
        " hapsolo_cli.py train -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta --paf ${{INPUT_FASTA_PREFIX}}_new_self_align.paf.gz -b ortholog_output -n {params.iterations_per_thread} >> ${{LOG}} 2>&1; "
        " tail -n +2 contigs/`basename {output.name_mappings}` > `basename {output.syn_file}` 2>>${{LOG}}; "
        " cp -f  asms/*_primary.fasta ${{OUTPUT_FASTA}}; "
        " cp -f  asms/*_secondary.fasta ${{OUTPUT_DUPLICATION_FASTA}}; "
        #" cat asms/*_secondary.fasta | grep -P  '^>' | sed 's/>//;s/[ \t].*//' > ${{OUTPUT_DUPLICATION_FASTA%.fasta}}.ids ; "

rule restore_orig_scaf_ids_after_hapsolo:
    input:
        purged_fasta=rules.hapsolo.output.purged_fasta,
        duplication_fasta=rules.hapsolo.output.duplication_fasta,
        syn_file=rules.hapsolo.output.syn_file,
        log_dir=ancient(rules.hapsolo.input.log_dir),
    output:
        purged_fasta="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.purged.fasta",
        purged_fasta_ids="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.purged.fasta.ids",
        duplication_fasta="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.dups.fasta",
        duplication_fasta_ids="{fasta_dir}/{fasta_prefix}/hapsolo/{busco_lineage}/{fasta_prefix}.dups.ids",
    log:
        log="{fasta_dir}/log/restore_orig_scaf_ids_after_hapsolo.{fasta_prefix}.{busco_lineage}.log",
        cluster_log="{fasta_dir}/log/restore_orig_scaf_ids_after_hapsolo.{fasta_prefix}.{busco_lineage}.cluster.log",
        cluster_err="{fasta_dir}/log/restore_orig_scaf_ids_after_hapsolo.{fasta_prefix}.{busco_lineage}.cluster.err"
    benchmark:
        "{fasta_dir}/log/restore_orig_scaf_ids_after_hapsolo.{fasta_prefix}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("gfa2fasta"),
        cpus=get_threads(parameters["threads"]["gfa2fasta"], "cpu"),
        time=parameters["time"]["gfa2fasta"],
        mem=parameters["memory_mb"]["gfa2fasta"]
    threads: parameters["threads"]["gfa2fasta"]

    shell:
        " rename_sequence_ids.py -k 1 -v 0 -i {input.purged_fasta} -s {input.syn_file} "
        "       -o {output.purged_fasta} >> {log.log} 2>&1; "
        " rename_sequence_ids.py -k 1 -v 0 -i {input.duplication_fasta} -s {input.syn_file} "
        "       -o {output.duplication_fasta} >> {log.log} 2>&1; "
        " get_sequence_ids.py -i {output.duplication_fasta} -o {output.duplication_fasta_ids} 2>>{log.log} ; "
        " get_sequence_ids.py -i {output.purged_fasta} -o {output.purged_fasta_ids} 2>>{log.log} ; "