#ruleorder: meryl_pe > create_fastq_links
rule meryl:
    input:
        output_dict["data"] / ("fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.{fileprefix}") #, (?!^histo$)
    log:
        std=output_dict["log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=parameters["memory_mb"]["meryl"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
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
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=parameters["memory_mb"]["meryl"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"

rule merge_meryl:
    input:
        lambda wildcards:
            expand(output_dict["kmer"] / ("%s/%s/%s.%s.%s.meryl.{fileprefix}" % (wildcards.datatype,
                                                                                 wildcards.stage,
                                                                                 wildcards.datatype,
                                                                                 wildcards.stage,
                                                                                 wildcards.kmer_length,)),
                   fileprefix=input_file_prefix_dict[wildcards.datatype],
                   allow_missing=True,)  if wildcards.datatype not in config["paired_fastq_based_data"] else \
            expand(rules.meryl_pe.output,
                   pairprefix=input_pairprefix_dict[wildcards.datatype],
                   allow_missing=True,)#peoutput_dict["kmer"] / ("%s/%s/%s.%s.%s.meryl.{pairprefix}" % (wildcards.datatype,
                   #                                                                 wildcards.stage,
                   #                                                                  wildcards.datatype,
                   #                                                                  wildcards.stage,
                   #                                                                  wildcards.kmer_length,)),
                   #pairprefix=input_pairprefix_dict[wildcards.datatype],
                   #    allow_missing=True,)
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
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=parameters["memory_mb"]["meryl"],
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl threads={threads} memory={resources.mem}m"
         " union-sum output {output.db_dir} {input} 1>{log.count_log} 2>&1;"
         " meryl threads={threads} memory={resources.mem}m "
         " histogram {output.db_dir} > {output.histo} 2>{log.histo_log}"

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
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl_extract"],
        time=parameters["time"]["meryl_extract"],
        mem=parameters["memory_mb"]["meryl_extract"],
    threads:
        parameters["threads"]["meryl_extract"]
    shell:
         #" OUTPUT={output.kmer}; "
         " meryl threads={threads} memory={resources.mem}m "
         " print less-than {wildcards.max_upper_boundary} greater-than {wildcards.min_lower_boundary}  {input.db} 2>{log.meryl} | "
         " sort > {output.kmer} 2>{log.sort};"
         #" pigz -p {threads} ${{OUTPUT%.gz}} 2>{log.pigz}; "

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
        #gunzip=output_dict["log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.gunzip.log",
        awk=output_dict["log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.awk.log",
        #gzip=output_dict["log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.gzip.log",
        cluster_log=output_dict["cluster_log"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    params:
        min_lower_boundary=min(parameters["tool_options"]["smudgeplot"]["lower_boundary"]),
        max_upper_boundary=max(parameters["tool_options"]["smudgeplot"]["upper_boundary"])
    benchmark:
        output_dict["benchmark"] / "subset_extracted_kmers.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
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
         #" gzip -c > {output.kmer} 2>{log.gzip}; "
