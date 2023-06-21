localrules: create_fastq_links, create_links_for_draft, create_fasta_links
rule create_fastq_links:
    priority: 1000
    input:
        input_dir_path.resolve() / ("{datatype}/fastq/{fileprefix}%s" %  config["fastq_extension"])
    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        output_dict["data"] / ("fastq/{datatype, [^/]+}/raw/{fileprefix, [^/]+}%s" % config["fastq_extension"])
    log:
        std=output_dict["log"] / "create_fastq_links.{datatype}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "create_fastq_links.{datatype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_fastq_links.{datatype}.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "create_fastq_links.{datatype}.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " ln -sf {input} {output} 2>{log.std}"

rule create_fasta_links:
    priority: 1000
    input:
        input_dir_path.resolve() / ("{datatype}/fasta/{fileprefix}%s" %  config["fasta_extension"])
    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        output_dict["data"] / ("fasta/{datatype, [^/]+}/raw/{fileprefix, [^/]+}%s" % config["fasta_extension"])
    log:
        std=output_dict["log"] / "create_fasta_links.{datatype}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "create_fasta_links.{datatype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_fasta_links.{datatype}.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "create_fasta_links.{datatype}.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " ln -sf {input} {output} 2>{log.std}"

rule create_links_for_draft:
    input:
        lambda wildcards: input_dir_path.resolve() / "draft/fasta/{0}".format(draft_file_dict[wildcards.haplotype])
    output:
        out_dir_path / "draft_qc/{parameters}/{genome_prefix}.draft_qc.{haplotype}.fasta"
        #paf=out_dir_path  / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    log:
        ln=output_dict["log"]  / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_links_for_draft"],
        time=parameters["time"]["create_links_for_draft"],
        mem=parameters["memory_mb"]["create_links_for_draft"]
    threads: parameters["threads"]["create_links_for_draft"]

    shell:
        " ln -sf {input} {output} 2>{log.ln}; "