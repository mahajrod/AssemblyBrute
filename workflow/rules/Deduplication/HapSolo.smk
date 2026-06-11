

rule hapsolo:
    input:
        input_fasta=lambda wildcards: out_dir_path  / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["dedup"]["parameters"][wildcards.prev_stage_parameters + "..hapsolo_" + wildcards.dedup_parameters]["prev_stage"],
                                                                                                                           wildcards.prev_stage_parameters,
                                                                                                                           wildcards.genome_prefix,
                                                                                                                           wildcards.haplotype),
        busco_db_dir=(out_dir_path / "download/busco5/lineages/{busco_lineage}").resolve(),
        log_dir=out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log",
    output:
        purged_fasta=out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype}.fasta",
    params:
        iterations_per_thread=lambda wildcards: stage_dict["dedup"]["parameters"][wildcards.prev_stage_parameters + "..hapsolo_" + wildcards.dedup_parameters]["option_set"]["iterations_per_thread"],
    log:
        preprocess=(out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.preprocess.log").resolve(),
        align=(out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.{genome_prefix}.{haplotype}.align.log").resolve(),
        search=(out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.search.log").resolve(),
        train=(out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.train.log").resolve(),
        cluster_log=out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        out_dir_path  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.benchmark.txt"
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
        " WORKDIR=`dirname {output.purged_fasta}`/{wildcards.haplotype}/{wildcards.busco_lineage}; "
        " cp -f {input.input_fasta} ${{WORKDIR}}; "
        " INPUT_FASTA=`basename {input.input_fasta}`; "
        " INPUT_FASTA_PREFIX=${{INPUT_FASTA%.fasta}}; "
        " OUTPUT_FASTA=`basename {input.input_fasta}`; "
        " cd ${{WORKDIR}}; "
        " hapsolo_cli.py preprocess -i ${{INPUT_FASTA}} > {log.preprocess} 2>&1;  "
        " hapsolo_cli.py align -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta > {log.align} 2>&1 ; "
        " hapsolo_cli.py search -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta -l {input.busco_db_dir} -o ortholog_output > {log.search} 2>&1; "
        " hapsolo_cli.py train -t {threads} -i ${{INPUT_FASTA_PREFIX}}_new.fasta --paf ${{INPUT_FASTA_PREFIX}}_new_self_align.paf.gz -b ortholog_output -n {params.iterations_per_thread} > {log.train} 2>&1; "
        " cp asms/*_primary.fasta ../../${{OUTPUT_FASTA}}; "

