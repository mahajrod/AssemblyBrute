localrules: create_fastq_links, create_links_for_draft, create_fasta_links
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
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " ln -sf {input} {output} 2>{log.std}"

checkpoint preprocess_hic_fastq:
    priority: 1000
    input:
        forward_reads=expand(input_dir_path.resolve() / ("hic/fastq/{pairprefix}%s%s" %  (input_forward_suffix_dict["hic"],
                                                                                         config["fastq_extension"])),
               pairprefix=input_pairprefix_dict["hic"]),
        reverse_reads=expand(input_dir_path.resolve() / ("hic/fastq/{pairprefix}%s%s" %  (input_forward_suffix_dict["hic"],
                                                                                         config["fastq_extension"])),
               pairprefix=input_pairprefix_dict["hic"])
    params:
        lines_per_chunk=parameters["tool_options"]["preprocess_hic_fastq"]["reads_per_chunk"],
        chunk_prefix=parameters["tool_options"]["preprocess_hic_fastq"]["chunk_prefix"],
        arima_trim="fastx_trimmer -f {0} | ".format(parameters["tool_options"]["preprocess_hic_fastq"]["arima_trim"] + 1) if (config["hic_enzyme_set"] == "Arima") and (not config["skip_filter_reads"]) else ""
    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        dir=directory(output_dict["data"] / "fastq/hic/raw/") # {fileprefix, [^/]+}%s" % config["fastq_extension"])
    log:
        zcat_forward=output_dict["log"] / "preprocess_hic_fastq.zcat_forward.log",
        split_forward=output_dict["log"] / "preprocess_hic_fastq.split_forward.log",
        zcat_reverse=output_dict["log"] / "preprocess_hic_fastq.zcat_reverse.log",
        split_reverse=output_dict["log"] / "preprocess_hic_fastq.split_reverse.log",
        cluster_log=output_dict["cluster_log"] / "preprocess_hic_fastq.cluster.log",
        cluster_err=output_dict["cluster_error"] / "preprocess_hic_fastq.cluster.err",
    benchmark:
        output_dict["benchmark"] / "preprocess_hic_fastq.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["preprocess_hic_fastq"],
        time=parameters["time"]["preprocess_hic_fastq"],
        mem=parameters["memory_mb"]["preprocess_hic_fastq"],
    threads:
        parameters["threads"]["preprocess_hic_fastq"]
    shell:
         " mkdir -p {output.dir}; "
         " zcat {input.forward_reads} 2>{log.zcat_forward} | {params.arima_trim} "
         " split --filter 'pigz -p 30 > ${{FILE}}.gz' --additional-suffix=_1.fastq -l 80000000 - {params.chunk_prefix} > {log.split_forward} 2>&1; "
         " zcat {input.reverse_reads} 2>{log.zcat_reverse} | {params.arima_trim} "
         " split --filter 'pigz -p 30 > ${{FILE}}.gz' --additional-suffix=_2.fastq -l 80000000 - {params.chunk_prefix} > {log.split_reverse} 2>&1; "

def get_hic_chunk_pairprefix_list():
    checkpoint_output = checkpoints.preprocess_hic_fastq.output.dir
    return glob_wildcards(checkpoint_output + "/{pairprefix}_1.fastq.gz").pairprefix

def get_hic_chunk_fileprefix_list():
    checkpoint_output = checkpoints.preprocess_hic_fastq.output.dir
    return glob_wildcards(checkpoint_output + "/{fileprefix}.fastq.gz").fileprefix

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