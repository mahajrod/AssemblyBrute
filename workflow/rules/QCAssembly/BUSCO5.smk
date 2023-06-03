#localrules: handle_busco5_output

rule busco5:
    priority: 500
    input:
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        tar_gz=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.tar.gz",
        summary=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.summary",
        summary_json=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.summary.json",
        busco_table=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
        missing_busco_ids=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.missing.ids",
    log:
        std=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.log",
        pigz=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.pigz.log",
        cp=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cp.log",
        cluster_log=output_dict["cluster_log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5"],
        time=parameters["time"]["busco5"],
        mem=parameters["memory_mb"]["busco5"],
    threads:
        parameters["threads"]["busco5"]
    shell:
         " BUSCO_DIR={output.tar_gz}; "
         " BUSCO_DIR=${{BUSCO_DIR%.tar.gz}}; "
         " busco -m genome -l {wildcards.busco_lineage} -c {threads} -i {input.assembly} "
         " -o `basename ${{BUSCO_DIR}}` --out_path `dirname ${{BUSCO_DIR}}` > {log.std} 2>&1;"
         " cp ${{BUSCO_DIR}}/short_summary.specific.{wildcards.busco_lineage}.{wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.haplotype}.busco5.{wildcards.busco_lineage}.txt {output.summary} ; "
         " cp ${{BUSCO_DIR}}/short_summary.specific.{wildcards.busco_lineage}.{wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.haplotype}.busco5.{wildcards.busco_lineage}.json {output.summary_json} ; "
         " cp ${{BUSCO_DIR}}/run_{wildcards.busco_lineage}/full_table.tsv {output.busco_table} ; "
         " cp ${{BUSCO_DIR}}/run_{wildcards.busco_lineage}/missing_busco_list.tsv {output.missing_busco_ids} ; "
         " tar cf - ${{BUSCO_DIR}} | pigz -p {threads} > {output.tar_gz} 2>{log.pigz} ;"
         " rm -r ${{BUSCO_DIR}}; "

"""
rule compress_busco5:
    priority: 500
    input:
        dir=rules.busco5.output.dir,
        #summary=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/short_summary.specific.{busco_lineage}.{genome_prefix}.{assembly_stage}.{haplotype}.txt",
        #summary_json=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/short_summary.specific.{busco_lineage}.{genome_prefix}.{assembly_stage}.{haplotype}.json",
        #busco_table=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/run_{busco_lineage}/full_table.tsv",
        #missing_busco_ids=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}/run_{busco_lineage}/missing_busco_list.tsv",
    output:
        tar_gz=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.tar.gz",

    log:
        tar=output_dict["log"] / "compress_busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.tar.log",
        cluster_log=output_dict["cluster_log"] / "compress_busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "compress_busco5{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "compress_busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["compress_busco5"],
        time=parameters["time"]["compress_busco5"],
        mem=parameters["memory_mb"]["compress_busco5"],
    threads:
        parameters["threads"]["compress_busco5"]
    shell:
         " tar czf {output.tar_gz} {input.dir} > {log.tar}2>&1;"
"""