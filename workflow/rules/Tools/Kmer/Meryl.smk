

rule meryl_assembly:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        db_dir=directory("{fasta_dir}/{fasta_prefix}/kmer/meryl/{fasta_prefix}.{kmer_length}.meryl/")
    log:
        std="{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.log",
        cluster_log="{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.cluster.log",
        cluster_err="{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.cluster.err"
    benchmark:
        "{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_assembly"),
        cpus=get_threads(parameters["threads"]["meryl_assembly"], "cpu"),
        time=parameters["time"]["meryl_assembly"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl_assembly"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl_assembly"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input.fasta} > {log.std} 2>&1;"

rule meryl_get_repetitive_kmers:
    input:
        meryl_db="{directory}/{meryl_db_prefix}.meryl/",
        log_dir=ancient("{directory}/log/"),
    output:
        repetitive_kmers="{directory}/{meryl_db_prefix}.meryl.repetitive"
    params:
        distinct=0.9998
    log:
        std="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.log",
        cluster_log="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.cluster.log",
        cluster_err="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.cluster.err"
    benchmark:
        "{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_extract"),
        cpus=parameters["threads"]["meryl_extract"],
        time=parameters["time"]["meryl_extract"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl_extract"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl_extract"]
    shell:
         " meryl print greater-than distinct={params.distinct} "
         " {input.meryl_db} > {output.repetitive_kmers} 2>{log.std}; "


rule meryl_se:
    input:
        lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}".format(wildcards.se_datatype,
                                                                               wildcards.stage,
                                                                               wildcards.fileprefix,
                                                                               config["data"][wildcards.se_datatype]["conv_ext"])
    output:
        db_dir=directory(config["out_dir"] / "kmer/{se_datatype}/{stage}/{se_datatype}.{stage}.{kmer_length}.meryl.{fileprefix}"),

    log:
        std=config["out_dir"] / "log/meryl.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.log",
        cluster_log=config["out_dir"] / "log/meryl.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/meryl.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "log/meryl.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl"),
        cpus=get_threads(parameters["threads"]["meryl"], "cpu"),
        time=parameters["time"]["meryl"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         "       output {output.db_dir} {input} 1>{log.std} 2>&1;"


rule meryl_pe:
    input:
        forward_fastq=lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                wildcards.stage,
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.pe_datatype]["conv_fwd_sfx"],
                                                                                                config["data"][wildcards.pe_datatype]["conv_ext"]),
        reverse_fastq=lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                wildcards.stage,
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.pe_datatype]["conv_rev_sfx"],
                                                                                                config["data"][wildcards.pe_datatype]["conv_ext"]),
    output:
        db_dir=directory(config["out_dir"] / "kmer/{pe_datatype}/{stage}/{pe_datatype}.{stage}.{kmer_length}.meryl.{pairprefix}")
    log:
        std=config["out_dir"] / "log/meryl.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.log",
        cluster_log=config["out_dir"] / "log/meryl.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/meryl.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "log/meryl.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_pe"),
        cpus=get_threads(parameters["threads"]["meryl"], "cpu"),
        time=parameters["time"]["meryl"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         "       output {output.db_dir} {input} 1>{log.std} 2>&1;"

def get_meryl_dbs_for_merging(wildcards):
    db_list = []
    for datatype in wildcards.datatype.split("_"):
        if datatype in config["data_feature_dict"]["paired"]:
            db_list += expand(rules.meryl_pe.output,
                              pairprefix=config["data"][datatype]["pair_prefix_list"],
                              pe_datatype=[datatype,],
                              allow_missing=True)
        else:
            db_list += expand(rules.meryl_se.output,
                              fileprefix=config["data"][datatype]["conv_file_prefix_list"], # for se_reads "conv_file_prefix_list" and "file_prefix_list" are the same
                              se_datatype=[datatype,],
                              allow_missing=True)

    return db_list

rule merge_meryl:
    input: get_meryl_dbs_for_merging
    output:
        db_dir=directory(config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"),
    log:
        count_log=config["out_dir"] / "log/merge_meryl.{datatype}.{stage}.{kmer_length}.count.log",
        cluster_log=config["out_dir"] / "log/merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "log/merge_meryl.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("merge_meryl"),
        cpus=get_threads(parameters["threads"]["meryl"], "cpu"),
        time=parameters["time"]["meryl"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl"],
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl threads={threads} memory={resources.mem}m union-sum output {output} {input} >{log.count_log} 2>&1;"

rule get_meryl_histo:
    input:
        db="{meryl_dir}/{meryl_db}/",
        log_dir=ancient("{meryl_dir}/log/")
    output:
        histo="{meryl_dir}/{meryl_db}.histo"

    log:
        histo_log="{meryl_dir}/log/get_meryl_histo.{meryl_db}.log",
        cluster_log="{meryl_dir}/log/get_meryl_histo.{meryl_db}.cluster.log",
        cluster_err="{meryl_dir}/log/get_meryl_histo.{meryl_db}.cluster.err"
    benchmark:
        "{meryl_dir}/log/get_meryl_histo.{meryl_db}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_meryl_histo"),
        cpus=parameters["threads"]["get_meryl_histo"],
        time=parameters["time"]["get_meryl_histo"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["get_meryl_histo"],
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl threads={threads} memory={resources.mem}m histogram {input.db} > {output.histo} 2>{log.histo_log}"


rule meryl_extract:
    input:
        db=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"
    output:
        kmer=temp(config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.L{min_lower_boundary}.U{max_upper_boundary}.extracted.kmer")
    log:
        meryl=config["out_dir"] / "log/meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.meryl.log",
        sort=config["out_dir"] / "log/meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.sort.log",
        pigz=config["out_dir"] / "log/meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.pigz.log",
        cluster_log=config["out_dir"] / "log/meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.cluster.log",
        cluster_err=config["out_dir"] / "log/meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.cluster.err"
    benchmark:
        config["out_dir"] / "log/meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_extract"),
        cpus=parameters["threads"]["meryl_extract"],
        time=parameters["time"]["meryl_extract"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl_extract"],
    threads:
        parameters["threads"]["meryl_extract"]
    shell:
         " meryl threads={threads} memory={resources.mem}m "
         "     print less-than {wildcards.max_upper_boundary} greater-than {wildcards.min_lower_boundary} {input.db} 2>{log.meryl} | "
         " sort > {output.kmer} 2>{log.sort};"

rule subset_extracted_kmers:
    input:
        kmer=lambda wildcards: config["out_dir"] / "kmer/{0}/{1}/{0}.{1}.{2}.meryl.L{3}.U{4}.extracted.kmer".format(wildcards.datatype,
                                                                                                                 wildcards.stage,
                                                                                                                 wildcards.kmer_length,
                                                                                                                 min(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["lower_boundary"]),
                                                                                                                 max(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["upper_boundary"]))
    output:
        kmer=temp(config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer")
    log:
        cp=config["out_dir"] / "log/subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cp.log",
        awk=config["out_dir"] / "log/subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.awk.log",
        cluster_log=config["out_dir"] / "log/subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=config["out_dir"] / "log/subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    params:
        min_lower_boundary=lambda wildcards: min(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["lower_boundary"]),
        max_upper_boundary=lambda wildcards: max(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["upper_boundary"])
    benchmark:
        config["out_dir"] / "log/subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("subset_extracted_kmers"),
        cpus=parameters["threads"]["subset_extracted_kmers"],
        time=parameters["time"]["subset_extracted_kmers"],
        mem=parameters["memory_mb"]["subset_extracted_kmers"],
    threads:
        parameters["threads"]["subset_extracted_kmers"]
    shell:
         " if [ {wildcards.lower_boundary} -eq {params.min_lower_boundary} ] && [ {wildcards.upper_boundary} -eq {params.max_upper_boundary} ] ;"
         " then "
         "     cp {input.kmer} {output.kmer} 2>{log.cp}; "
         " else "
         "     awk '{{if (($2 >= {wildcards.lower_boundary}) && ($2 <= {wildcards.upper_boundary})) print $0}}' {input.kmer} > {output.kmer} 2>{log.awk}; "
         " fi "
