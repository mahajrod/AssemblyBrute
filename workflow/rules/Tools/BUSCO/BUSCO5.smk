

rule busco5_download: #
    priority: 500
    output:
        lineage_dir=directory(config["out_dir"] / "download/busco5/lineages/{busco5_lineage}"),
    params:
        busco_download_dir=config["out_dir"] / "download/busco5/"
    log:
        std=config["out_dir"] / "log/busco5_download.{busco5_lineage}.log",
        cluster_log=config["out_dir"] / "log/busco5_download.{busco5_lineage}.cluster.log",
        cluster_err=config["out_dir"] / "log/busco5_download.{busco5_lineage}.cluster.err"
    benchmark:
        config["out_dir"] / "log/busco5_download.{busco5_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco5.8"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco5.8"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco_download"),
        cpus=parameters["threads"]["busco_download"],
        time=parameters["time"]["busco_download"],
        mem=parameters["memory_mb"]["busco_download"],
    threads:
        parameters["threads"]["busco_download"]
    shell:
         " busco --download_path {params.busco_download_dir} --download {wildcards.busco5_lineage} > {log.std} 2>&1; "

rule busco5_8_odb12: # Downloading of busco datasets is performed by a different rule to avoid conflict between different instances of busco5
    priority: 500
    input:
        busco_lineage=config["out_dir"] / "download/busco5/lineages/{busco5_lineage}",
        assembly="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        tar_gz="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb12.*}.busco5.tar.gz",
        summary="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb12.*}.busco5.summary",
        summary_json="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb12.*}.busco5.summary.json",
        busco_table="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb12.*}.busco5.full_table.tsv",
        missing_busco_ids="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb12.*}.busco5.missing.ids",
    params:
        gene_predictor="--{0}".format(parameters["tool_options"]["busco"]["gene_predictor"])
    log:
        std="{fasta_dir}/log/busco5.{fasta_prefix}.{busco5_lineage}.log",
        pigz="{fasta_dir}/log/busco5.{fasta_prefix}.{busco5_lineage}.pigz.log",
        cp="{fasta_dir}/log/busco5.{fasta_prefix}.{busco5_lineage}.cp.log",
        cluster_log="{fasta_dir}/log/busco5.{fasta_prefix}.{busco5_lineage}.cluster.log",
        cluster_err="{fasta_dir}/log/busco5.{fasta_prefix}.{busco5_lineage}.cluster.err"
    benchmark:
        "{fasta_dir}/log/busco5.{fasta_prefix}.{busco5_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco5.8"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco5.8"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco5"),
        cpus=parameters["threads"]["busco5"],
        time=parameters["time"]["busco5"],
        mem=parameters["memory_mb"]["busco5"],
    threads:
        parameters["threads"]["busco5"]
    shell:
         " BUSCO_DIR={output.tar_gz}; "
         " BUSCO_DIR=${{BUSCO_DIR%.tar.gz}}; "
         " busco --offline -f -m genome -l {input.busco_lineage} -c {threads} -i {input.assembly} {params.gene_predictor}"
         "     -o `basename ${{BUSCO_DIR}}` --out_path `dirname ${{BUSCO_DIR}}` > {log.std} 2>&1;"
         " cp ${{BUSCO_DIR}}/short_summary.specific.{wildcards.busco5_lineage}.{wildcards.fasta_prefix}.{wildcards.busco5_lineage}.busco5.txt {output.summary} ; "
         " cp ${{BUSCO_DIR}}/short_summary.specific.{wildcards.busco5_lineage}.{wildcards.fasta_prefix}.{wildcards.busco5_lineage}.busco5.json {output.summary_json} ; "
         " cp ${{BUSCO_DIR}}/run_{wildcards.busco5_lineage}/full_table.tsv {output.busco_table} ; "
         " cp ${{BUSCO_DIR}}/run_{wildcards.busco5_lineage}/missing_busco_list.tsv {output.missing_busco_ids} ; "
         " tar cf - ${{BUSCO_DIR}} | pigz -p {threads} > {output.tar_gz} 2>{log.pigz} ;"
         " rm -r ${{BUSCO_DIR}}; "

use rule busco5_8_odb12 as rule busco5_6_odb10 with: #BUSCO 5.8.3 returns strange results for vertebrata_odb10
    output:
        tar_gz="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb10.*}.busco5.tar.gz",
        summary="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb10.*}.busco5.summary",
        summary_json="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb10.*}.busco5.summary.json",
        busco_table="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb10.*}.busco5.full_table.tsv",
        missing_busco_ids="{fasta_dir}/assembly_qc/busco5/{fasta_prefix}.{busco5_lineage, .*odb10.*}.busco5.missing.ids",
    conda:
        config["conda"]["busco5.6"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco5.6"]["yaml"])