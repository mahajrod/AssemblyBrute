#ruleorder: pretextmap > pretextsnapshot
rule pretextmap: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        #bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        len=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len"
    output:
        #map=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.map.pretext"  % config["genome_name"]),
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.filtered_out.ids"
    params:
        #min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        max_len=lambda wildcards: parameters["tool_options"]["pretextmap"]["subsets"][wildcards.subset]["max_len"],
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], "--sortby"),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], "--sortorder"),
    log:
        view=output_dict["log"]  / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.view.log",
        awk=output_dict["log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.awk.log",
        map=output_dict["log"]  / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.map.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pretextmap"),
        cpus=parameters["threads"]["pretextmap"] ,
        time=parameters["time"]["pretextmap"],
        mem=parameters["memory_mb"]["pretextmap"]
    threads: parameters["threads"]["pretextmap"]

    shell:
        " MAP_LOG=`realpath -s -m {log.map}` ; "
        " VIEW_LOG=`realpath -s -m {log.view}` ; "
        " if [ '{params.max_len}' == 'None' ];"
        "   then "
        "       > {output.filtered_out}; "
        "   else "
        "       awk '{{if ($2 > {params.max_len}) print $1}}' {input.len} > {output.filtered_out} 2>{log.awk}; "
        "   fi; "
        " if [[ -s {output.filtered_out} ]]; "
        "   then "
        "       FILTER_OUT=' --filterExclude '; "
        "       FILTER_OUT=\"${{FILTER_OUT}} `cat {output.filtered_out} | tr '\\n' ',' | sed 's/,\+$//'` | sed 's/,/, /g' \"; "
        "   else "
        "       FILTER_OUT=''; "
        "   fi; " 
        " cd `dirname {input.bam}`; "
        " samtools view -@4 -F0x400 -h `basename {input.bam}` 2>${{VIEW_LOG}} | "
        " PretextMap -o `basename {output.map}` {params.sortby} {params.sortorder} "
        "            --mapq {wildcards.mapq} ${{FILTER_OUT}} --highRes > ${{MAP_LOG}} 2>&1"

rule pretextsnapshot: #Pretext-snapshot doesn't support long file names!!!!!!!!!!!
    input:
        map=rules.pretextmap.output.map
        #map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.filtered.rmdup.map.pretext"
    output:
        dir=directory(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.mapq{mapq}.{resolution, [0-9]+}.{ext}"),
    params:
        sequences=parameters["tool_options"]["pretextsnapshot"]["sequences"],
    log:
        std=output_dict["log"]  / "pretextsnapshot.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{ext}.{resolution}.log",
        cluster_log=output_dict["cluster_log"] / "pretextsnapshot.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{ext}.{resolution}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextsnapshot.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{ext}.{resolution}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextsnapshot.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{ext}.{resolution}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pretextsnapshot"),
        cpus=parameters["threads"]["pretextmap"] ,
        time=parameters["time"]["pretextmap"],
        mem=parameters["memory_mb"]["pretextmap"]
    threads: parameters["threads"]["pretextmap"]
    shell:
        " LOG=`realpath -s -m {log.std}`; "
        " cd `dirname {input.map}`; "
        " PretextSnapshot --sequences {params.sequences} -r {wildcards.resolution} -f {wildcards.ext} "
        " -m `basename {input.map}` -o `basename {output.dir}`  > ${{LOG}} 2>&1"

