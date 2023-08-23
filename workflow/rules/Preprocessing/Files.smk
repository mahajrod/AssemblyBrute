localrules: create_fastq_links, create_links_for_draft, create_fasta_links, create_links_for_reference
ruleorder: preprocess_hic_fastq > create_fastq_links
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
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " ln -sf {input} {output} 2>{log.std}"

rule preprocess_hic_fastq:
    priority: 2000
    input:
        input_dir_path.resolve() / ("hic/fastq/{fileprefix}%s" %  config["fastq_extension"])
    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        output_dict["data"] / ("fastq/hic/raw/{fileprefix, [^/]+}%s" % config["fastq_extension"])
    params:
        hic_type=config["hic_enzyme_set"],
        skip_trimming='skip' if config["skip_filter_reads"] else 'trim'
    log:
        std=output_dict["log"] / "preprocess_hic_fastq.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "preprocess_hic_fastq.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "preprocess_hic_fastq.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "preprocess_hic_fastq.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["preprocess_hic_fastq"],
        time=parameters["time"]["preprocess_hic_fastq"],
        mem=parameters["memory_mb"]["preprocess_hic_fastq"],
    threads:
        parameters["threads"]["preprocess_hic_fastq"]
    shell:
         " if [ '{params.hic_type}' = 'Arima' -a '{params.skip_trimming}' = 'trim' ]; "
         " then "
         "      zcat {input} | fastx_trimmer -f 8 | pigz -p {threads} > {output} 2>{log.std}; "
         " else "
         "      ln -sf {input} {output} 2>{log.std}; "
         " fi; "

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
        queue=config["queue"]["cpu"],
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
    log:
        ln=output_dict["log"]  / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["create_links_for_draft"],
        time=parameters["time"]["create_links_for_draft"],
        mem=parameters["memory_mb"]["create_links_for_draft"]
    threads: parameters["threads"]["create_links_for_draft"]

    shell:
        " ln -sf {input} {output} 2>{log.ln}; "

rule create_links_for_reference:
    input:
        fasta=lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        syn=lambda wildcards: input_reference_filedict[wildcards.ref_name]["syn"].resolve(),
        whitelist=lambda wildcards: input_reference_filedict[wildcards.ref_name]["whitelist"].resolve(),
        orderlist=lambda wildcards: input_reference_filedict[wildcards.ref_name]["orderlist"].resolve(),
    output:
        fasta=out_dir_path / "data/reference/{ref_name}/{ref_name}.softmasked.fasta",
        syn=out_dir_path / "data/reference/{ref_name}/{ref_name}.syn",
        whitelist=out_dir_path / "data/reference/{ref_name}/{ref_name}.whitelist",
        orderlist=out_dir_path / "data/reference/{ref_name}/{ref_name}.orderlist",
    log:
        ln=output_dict["log"]  / "create_links_for_reference.{ref_name}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_links_for_reference.{ref_name}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_links_for_reference.{ref_name}.err"
    benchmark:
        output_dict["benchmark"]  / "create_links_for_reference.{ref_name}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["create_links_for_draft"],
        time=parameters["time"]["create_links_for_draft"],
        mem=parameters["memory_mb"]["create_links_for_draft"]
    threads: parameters["threads"]["create_links_for_draft"]

    shell:
        " ln -sf {input.fasta} {output.fasta} 2>{log.ln}; "
        " ln -sf {input.syn} {output.syn} 2>>{log.ln}; "
        " ln -sf {input.whitelist} {output.whitelist} 2>>{log.ln}; "
        " ln -sf {input.orderlist} {output.orderlist} 2>>{log.ln}; "
