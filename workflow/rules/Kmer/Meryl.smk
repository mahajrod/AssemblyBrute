#ruleorder: meryl_pe > create_fastq_links

rule meryl_assembly_new: #TODO: in future use this rule for all kmer counts on assemblies -  other rule (meryl_assembly) exists in ReadPhasing.smk
    input:
        fasta="{directory}/{fasta_prefix}.fasta",
        log_dir="{directory}/log/",
        benchmark_dir="{directory}/benchmark/"
    output:
        db_dir=directory("{directory}/kmer/meryl/{fasta_prefix, [^/]+}.{kmer_length, [^/]+}.meryl/")
    log:
        std="{directory}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.log",
        cluster_log="{directory}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.cluster.log",
        cluster_err="{directory}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.cluster.err"
    benchmark:
        "{directory}/benchmark/meryl_assembly.{fasta_prefix}.{kmer_length}.benchmark.txt"
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
         " workflow/external_tools/meryl-1.4/bin/meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input.fasta} > {log.std} 2>&1;"

rule meryl_get_repetitive_kmers:
    input:
        meryl_db="{directory}/{meryl_db_prefix}.meryl/",
        log_dir="{directory}/log/",
        benchmark_dir="{directory}/benchmark/"
    output:
        repetitive_kmers="{directory}/{meryl_db_prefix}.meryl.repetitive"
    params:
        distinct=0.9998
    log:
        std="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.log",
        cluster_log="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.cluster.log",
        cluster_err="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.cluster.err"
    benchmark:
        "{directory}/benchmark/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.benchmark.txt"
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
         " workflow/external_tools/meryl-1.4/bin/meryl print greater-than distinct={params.distinct} "
         " {input.meryl_db} > {output.repetitive_kmers} 2>{log.std}; "


rule meryl:
    input:
        lambda wildcards: output_dict["data"] / "{0}/{1}/{2}/{3}{4}".format(datatype_format_dict[wildcards.datatype],
                                                                            wildcards.datatype,
                                                                            wildcards.stage,
                                                                            wildcards.fileprefix,
                                                                            config[datatype_format_dict[wildcards.datatype] + "_extension"])
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.{fileprefix}"),
    #wildcard_constraints:
    #    fileprefix="^(?!histo).*",
    #    kmer_length="[^./]+"#, (?!^histo$)
    log:
        std=output_dict["log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.benchmark.txt"
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
         " workflow/external_tools/meryl-1.4/bin/meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"


rule meryl_pe:
    input:
        forward_fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{pairprefix}_2%s" % config["fastq_extension"]),
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.{pairprefix}") # , (?!^histo$)
    log:
        std=output_dict["log"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.benchmark.txt"
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
         " workflow/external_tools/meryl-1.4/bin/meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"

def get_meryl_dbs_for_merging(wildcards):
    db_list = []
    for datatype in wildcards.datatype.split("_"):
        db_list += expand(output_dict["kmer"] / ("%s/%s/%s.%s.%s.meryl.{fileprefix}" % (datatype,
                                                                                              wildcards.stage,
                                                                                              datatype,
                                                                                              wildcards.stage,
                                                                                              wildcards.kmer_length,)),
            fileprefix=input_file_prefix_dict[datatype] if datatype_format_dict[datatype] == "fastq" else
                       input_fasta_file_prefix_dict[datatype], allow_missing=True,) if datatype not in config["paired_fastq_based_data"] else \
            expand(rules.meryl_pe.output,
                pairprefix=input_pairprefix_dict[datatype],
                allow_missing=True,)

    return db_list
"""
rule merge_meryl:
    input: get_meryl_dbs_for_merging
        #lambda wildcards:
        #    expand(output_dict["kmer"] / ("%s/%s/%s.%s.%s.meryl.{fileprefix}" % (wildcards.datatype,
        #                                                                         wildcards.stage,
        #                                                                         wildcards.datatype,
        #                                                                         wildcards.stage,
        #                                                                         wildcards.kmer_length,)),
        #           fileprefix=input_file_prefix_dict[wildcards.datatype] if datatype_format_dict[wildcards.datatype] == "fastq" else input_fasta_file_prefix_dict[wildcards.datatype],
        #           allow_missing=True,)  if wildcards.datatype not in config["paired_fastq_based_data"] else \
        #    expand(rules.meryl_pe.output,
        #           pairprefix=input_pairprefix_dict[wildcards.datatype],
        #           allow_missing=True,)
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"),
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.histo"

    log:
        count_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.count.log",
        histo_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.histo.log",
        cluster_log=output_dict["cluster_log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("merge_meryl"),
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl"],
    threads:
        parameters["threads"]["meryl"]
    shell:
         " workflow/external_tools/meryl-1.4/bin/meryl threads={threads} memory={resources.mem}m"
         " union-sum output {output.db_dir} {input} 1>{log.count_log} 2>&1;"
         " workflow/external_tools/meryl-1.4/bin/meryl threads={threads} memory={resources.mem}m "
         " histogram {output.db_dir} > {output.histo} 2>{log.histo_log}"
"""


rule merge_meryl:
    input: get_meryl_dbs_for_merging
        #lambda wildcards:
        #    expand(output_dict["kmer"] / ("%s/%s/%s.%s.%s.meryl.{fileprefix}" % (wildcards.datatype,
        #                                                                         wildcards.stage,
        #                                                                         wildcards.datatype,
        #                                                                         wildcards.stage,
        #                                                                         wildcards.kmer_length,)),
        #           fileprefix=input_file_prefix_dict[wildcards.datatype] if datatype_format_dict[wildcards.datatype] == "fastq" else input_fasta_file_prefix_dict[wildcards.datatype],
        #           allow_missing=True,)  if wildcards.datatype not in config["paired_fastq_based_data"] else \
        #    expand(rules.meryl_pe.output,
        #           pairprefix=input_pairprefix_dict[wildcards.datatype],
        #           allow_missing=True,)
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"),
        #histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.histo"

    log:
        count_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.count.log",
        histo_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.histo.log",
        cluster_log=output_dict["cluster_log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.benchmark.txt"
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
         " workflow/external_tools/meryl-1.4/bin/meryl threads={threads} memory={resources.mem}m"
         " union-sum output {output.db_dir} {input} 1>{log.count_log} 2>&1;"
         #" workflow/external_tools/meryl-1.4/bin/meryl threads={threads} memory={resources.mem}m "
         #" histogram {output.db_dir} > {output.histo} 2>{log.histo_log}"


rule get_meryl_histo:
    input:
        db="{meryl_dir}/{meryl_db}/",
        log_dir="{meryl_dir}/log/"
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
         " workflow/external_tools/meryl-1.4/bin/meryl threads={threads} memory={resources.mem}m "
         " histogram {input.db} > {output.histo} 2>{log.histo_log}"


rule meryl_extract:
    input:
        db=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"
    output:
        kmer=temp(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.L{min_lower_boundary}.U{max_upper_boundary}.extracted.kmer")
    log:
        meryl=output_dict["log"] / "meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.meryl.log",
        sort=output_dict["log"] / "meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.sort.log",
        pigz=output_dict["log"] / "meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.pigz.log",
        cluster_log=output_dict["cluster_log"] / "meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl_extract.{datatype}.{stage}.{kmer_length}.L{min_lower_boundary}.U{max_upper_boundary}.benchmark.txt"
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
         " workflow/external_tools/meryl-1.4/bin/meryl threads={threads} memory={resources.mem}m "
         " print less-than {wildcards.max_upper_boundary} greater-than {wildcards.min_lower_boundary}  {input.db} 2>{log.meryl} | "
         " sort > {output.kmer} 2>{log.sort};"

rule subset_extracted_kmers:
    input:
        kmer=lambda wildcards: output_dict["kmer"] / "{0}/{1}/{0}.{1}.{2}.meryl.L{3}.U{4}.extracted.kmer".format(wildcards.datatype,
                                                                                                                 wildcards.stage,
                                                                                                                 wildcards.kmer_length,
                                                                                                                 min(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["lower_boundary"]),
                                                                                                                 max(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["upper_boundary"]))
    output:
        kmer=temp(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer")
    log:
        cp=output_dict["log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cp.log",
        awk=output_dict["log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.awk.log",
        cluster_log=output_dict["cluster_log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    params:
        min_lower_boundary=lambda wildcards: min(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["lower_boundary"]),
        max_upper_boundary=lambda wildcards: max(parameters["tool_options"]["smudgeplot"][wildcards.datatype]["upper_boundary"])
    benchmark:
        output_dict["benchmark"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
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
         "  then "
         " cp {input.kmer} {output.kmer} 2>{log.cp}; "
         " else "
         " awk '{{if (($2 >= {wildcards.lower_boundary}) && ($2 <= {wildcards.upper_boundary})) print $0}}' {input.kmer} > {output.kmer} 2>{log.awk}; "
         " fi "
