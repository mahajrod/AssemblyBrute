#ruleorder: pretextmap > pretextsnapshot
rule pretextmap: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        bam="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.bam",
        len="{fasta_dir}/{fasta_prefix}.len",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        map="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{subset, all}.rmdup.mapq{mapq}.{pretext_res}.pretext",
        filtered_out="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{subset, all}.rmdup.mapq{mapq}.{pretext_res}.filtered_out.ids"
    params:
        resolution=lambda wildcards: " --highRes" if wildcards.pretext_res == "high_res" else "",
        max_len=lambda wildcards: parameters["tool_options"]["pretextmap"]["subsets"][wildcards.subset]["max_len"],
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
    log:
        view="{fasta_dir}/log/pretextmap.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.view.log",
        awk="{fasta_dir}/log/pretextmap.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.awk.log",
        map="{fasta_dir}/log/pretextmap.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.map.log",
        cluster_log="{fasta_dir}/log/pretextmap.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.cluster.log",
        cluster_err="{fasta_dir}/log/pretextmap.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.cluster.err"
    benchmark:
        "{fasta_dir}/log/pretextmap.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pretextmap"),
        cpus=parameters["threads"]["pretextmap"] ,
        time=parameters["time"]["pretextmap"],
        mem=parameters["memory_mb"]["pretextmap"]
    threads: parameters["threads"]["pretextmap"]

    shell:
        " MAP_LOG=`realpath -s -m {log.map}` ; "
        " VIEW_LOG=`realpath -s -m {log.view}` ; "
        " if [ '{params.max_len}' == 'None' ];"
        " then "
        "     > {output.filtered_out}; "
        " else "
        "     awk '{{if ($2 > {params.max_len}) print $1}}' {input.len} > {output.filtered_out} 2>{log.awk}; "
        " fi; "
        " if [[ -s {output.filtered_out} ]]; "
        " then "
        "     FILTER_OUT=' --filterExclude '; "
        "     FILTER_OUT=\"${{FILTER_OUT}} `cat {output.filtered_out} | tr '\\n' ',' | sed 's/,\+$//'` \"; "
        " else "
        "     FILTER_OUT=''; "
        " fi; " 
        " cd `dirname {input.bam}`; "
        " samtools view -@4 -F0x400 -h `basename {input.bam}` 2>${{VIEW_LOG}} | "
        " PretextMap -o `basename {output.map}` {params.sortby} {params.sortorder} "
        "     --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1"

rule pretextsnapshot: #Pretext-snapshot doesn't support long file names!!!!!!!!!!!
    input:
        map=expand(rules.pretextmap.output.map, res=["default"], allow_missing=True),
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        dir=directory("{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{subset}.mapq{mapq}.default.{resolution}.{ext}"),
    params:
        sequences=parameters["tool_options"]["pretextsnapshot"]["sequences"],
    log:
        std="{fasta_dir}/log/pretextsnapshot.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{ext}.{resolution}.log",
        cluster_log="{fasta_dir}/log/pretextsnapshot.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{ext}.{resolution}.cluster.log",
        cluster_err="{fasta_dir}/log/pretextsnapshot.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{ext}.{resolution}.cluster.err"
    benchmark:
        "{fasta_dir}/log/pretextsnapshot.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{ext}.{resolution}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pretextsnapshot"),
        cpus=parameters["threads"]["pretextsnapshot"] ,
        time=parameters["time"]["pretextsnapshot"],
        mem=parameters["memory_mb"]["pretextsnapshot"]
    threads: parameters["threads"]["pretextsnapshot"]
    shell:
        " LOG=`realpath -s -m {log.std}`; "
        " cd `dirname {input.map}`; "
        " PretextSnapshot --sequences {params.sequences} -r {wildcards.resolution} -f {wildcards.ext} "
        "     -m `basename {input.map}` -o `basename {output.dir}`  > ${{LOG}} 2>&1"

