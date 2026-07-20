localrules: create_final_hapsolo_links

use rule create_local_links as create_final_hapsolo_links with: # abstract rule to use via redefining
    input:
        purged_fasta=config["out_dir"] / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.input.{haplotype}/{genome_prefix}.input.{haplotype}/hapsolo_qc/{busco_lineage}/{genome_prefix}.input.{haplotype}.purged.fasta",
        duplication_fasta=config["out_dir"] / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.input.{haplotype}/{genome_prefix}.input.{haplotype}/hapsolo_qc/{busco_lineage}/{genome_prefix}.input.{haplotype}.dups.fasta",
        duplication_ids=config["out_dir"] / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.input.{haplotype}/{genome_prefix}.input.{haplotype}/hapsolo_qc/{busco_lineage}/{genome_prefix}.input.{haplotype}.dups.ids",
        log_dir=ancient(config["out_dir"] / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/")
    output:
        purged_fasta=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype, hap[^_./@]+}.fasta",
        duplication_fasta=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype, hap[^_./@]+}.dups.fasta",
        duplication_ids=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype, hap[^_./@]+}.dups.ids",
    log:
        ln=config["out_dir"] / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/create_hapsolo_input_links.{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}.{genome_prefix}.dedup.{haplotype}.ln.log",


"""
rule hapsolo_assembly:
    input:
        input_fasta=lambda wildcards: config["out_dir"]  / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["dedup"].prev_stage,
                                                                                                    wildcards.prev_stage_parameters,
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.haplotype),
        busco_db_dir=config["out_dir"] / "download/busco5/lineages/{busco_lineage}",
        log_dir=ancient(config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log"),
    output:
        purged_fasta=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype, hap.*}.fasta",
        duplication_fasta=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype, hap.*}.dups.fasta",
        duplication_ids=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/{genome_prefix}.dedup.{haplotype, hap.*}.dups.ids",
    params:
        iterations_per_thread=lambda wildcards: stage_dict["dedup"].parameters["%s..hapsolo_%s@%s" % (wildcards.prev_stage_parameters,
                                                                                                         wildcards.dedup_parameters,
                                                                                                         wildcards.busco_lineage)]["option_set"]["iterations_per_thread"],
    log:
        preprocess=(config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.preprocess.log").resolve(),
        align=(config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.{genome_prefix}.{haplotype}.align.log").resolve(),
        search=(config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.search.log").resolve(),
        train=(config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.train.log").resolve(),
        cluster_log=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"]  / "dedup/{prev_stage_parameters}..hapsolo_{dedup_parameters}@{busco_lineage}/log/hapsolo.{busco_lineage}.{genome_prefix}.{haplotype}.benchmark.txt"
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
        " WORKDIR=`dirname {output.purged_fasta}`/{wildcards.genome_prefix}.dedup.{wildcards.haplotype}/hapsolo/{wildcards.busco_lineage}; "
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
        " cp -f  asms/*_primary.fasta ../../../${{OUTPUT_FASTA}}; "
        " cp -f  asms/*_secondary.fasta ../../../${{OUTPUT_DUPLICATION_FASTA}}; "
        " cat asms/*_secondary.fasta | grep -P  '^>' | sed 's/>//;s/[ \t].*//' > ../../../${{OUTPUT_DUPLICATION_FASTA%.fasta}}.ids ; "
"""