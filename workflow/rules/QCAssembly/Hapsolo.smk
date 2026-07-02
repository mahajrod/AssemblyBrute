rule hapsolo_qc:
    input:
        input_fasta="{fasta_dir}/{fasta_prefix}.fasta",
        busco_db_dir=out_dir_path / "download/busco5/lineages/{busco_lineage}",
        log_dir="{fasta_dir}/log",
    output:
        purged_fasta="{fasta_dir}/{fasta_prefix}/hapsolo_qc/{busco_lineage}/{fasta_prefix}.purged.fasta",
        duplication_fasta="{fasta_dir}/{fasta_prefix}/hapsolo_qc/{busco_lineage}/{fasta_prefix}.dups.fasta",
        duplication_ids="{fasta_dir}/{fasta_prefix}/hapsolo_qc/{busco_lineage}/{fasta_prefix}.dups.ids",
    params:
        iterations_per_thread=parameters["tool_options"]["assembly_qc"]["hapsolo"]["iterations_per_thread"]
    log:
        preprocess=("{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.preprocess.log").resolve(),
        align=("{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.align.log").resolve(),
        search=("{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.search.log").resolve(),
        train=("{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.train.log").resolve(),
        cluster_log="{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.cluster.log",
        cluster_err="{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.cluster.err"
    benchmark:
        "{fasta_dir}/log/hapsolo_qc.{fasta_prefix}.{busco_lineage}.benchmark.txt"
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
        " WORKDIR=`dirname {output.purged_fasta}; "
        " mkdir -p ${{WORKDIR}}; "
        " cp -f {input.input_fasta} ${{WORKDIR}}; "
        " INPUT_FASTA=`basename {input.input_fasta}`; "
        " INPUT_FASTA_PREFIX=${{INPUT_FASTA%.fasta}}; "
        " OUTPUT_FASTA=`basename {output.purged_fasta}`; "
        " OUTPUT_DUPLICATION_FASTA=`basename {output.duplication_fasta}`; "
        " cd ${{WORKDIR}}; "
        " hapsolo_cli.py preprocess -i ${{INPUT_FASTA}} > {log.preprocess} 2>&1;  "
        " hapsolo_cli.py align -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta > {log.align} 2>&1 ; "
        " hapsolo_cli.py search -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta -l ${{BUSCO_DIR}} -o ortholog_output > {log.search} 2>&1; "
        " hapsolo_cli.py train -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta --paf ${{INPUT_FASTA_PREFIX}}_new_self_align.paf.gz -b ortholog_output -n {params.iterations_per_thread} > {log.train} 2>&1; "
        " cp -f  asms/*_primary.fasta ${{OUTPUT_FASTA}}; "
        " cp -f  asms/*_secondary.fasta ${{OUTPUT_DUPLICATION_FASTA}}; "
        " cat asms/*_secondary.fasta | grep -P  '^>' | sed 's/>//;s/[ \t].*//' > ${{OUTPUT_DUPLICATION_FASTA%.fasta}}.ids ; "