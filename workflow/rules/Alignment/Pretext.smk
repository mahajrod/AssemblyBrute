#ruleorder: pretextmap > pretextsnapshot
rule pretextmap: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        #bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        len=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len"
    output:
        #map=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.map.pretext"  % config["genome_name"]),
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.filtered_out.ids"
    params:
        #min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        resolution=lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
        max_len=lambda wildcards: parameters["tool_options"]["pretextmap"]["subsets"][wildcards.subset]["max_len"],
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], "--sortby"),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], "--sortorder"),
    log:
        view=output_dict["log"]  / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.view.log",
        awk=output_dict["log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.awk.log",
        map=output_dict["log"]  / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.map.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.benchmark.txt"
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
        "            --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1"

rule pretextsnapshot: #Pretext-snapshot doesn't support long file names!!!!!!!!!!!
    input:
        map=expand(rules.pretextmap.output.map, res=["default"], allow_missing=True)
        #map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.filtered.rmdup.map.pretext"
    output:
        dir=directory(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.mapq{mapq}.default.{resolution, [0-9]+}.{ext}"),
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
        cpus=parameters["threads"]["pretextsnapshot"] ,
        time=parameters["time"]["pretextsnapshot"],
        mem=parameters["memory_mb"]["pretextsnapshot"]
    threads: parameters["threads"]["pretextsnapshot"]
    shell:
        " LOG=`realpath -s -m {log.std}`; "
        " cd `dirname {input.map}`; "
        " PretextSnapshot --sequences {params.sequences} -r {wildcards.resolution} -f {wildcards.ext} "
        " -m `basename {input.map}` -o `basename {output.dir}`  > ${{LOG}} 2>&1"


rule pretext_inject_tracks:
    input:
        #bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.filtered_out.ids",
        gap_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.gap.track.bedgraph",
        canonical_telomere_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
        non_canonical_telomere_track=out_dir_path/ "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.bedgraph",
        gc_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph",
        trf_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph",
        windowmasker_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph",
        hifi_coverage_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph"

    output:
        updated_map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.tracks.win_{window, [0-9]+}.{step, [0-9]+}.pretext"
    params:
        min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        #resolution=lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
        #max_len=lambda wildcards: parameters["tool_options"]["pretextmap"]["subsets"][wildcards.subset]["max_len"],
        #sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], "--sortby"),
        #sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], "--sortorder"),
    log:
        gap=output_dict["log"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.gap.log",
        can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.can_tel.log",
        non_can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.non_can_tel.log",
        gc=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.gc.log",
        trf=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.trf.log",
        windowmasker=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.windowmasker.log",
        coverage=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.coverage.log",
        awk=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.awk.log",
        rm=output_dict["log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.rm.log",
        cluster_log=output_dict["cluster_log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pretext_inject_tracks"),
        cpus=parameters["threads"]["pretext_inject_tracks"] ,
        time=parameters["time"]["pretext_inject_tracks"],
        mem=parameters["memory_mb"]["pretext_inject_tracks"]
    threads: parameters["threads"]["pretext_inject_tracks"]

    shell:
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gap_track}  | "
        "   PretextGraph -i {input.map} -n gap -o {output.updated_map}.tmp1 > {log.gap} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.canonical_telomere_track} | "
        "   PretextGraph -i {output.updated_map}.tmp1 -n can_telomere -o {output.updated_map}.tmp2 > {log.can_tel} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.non_canonical_telomere_track} | "
        "   PretextGraph -i {output.updated_map}.tmp2 -n noncan_telomere -o {output.updated_map}.tmp3 > {log.non_can_tel} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gc_track} | "
        "   PretextGraph -i {output.updated_map}.tmp3 -n GC -o {output.updated_map}.tmp4 > {log.gc} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.trf_track} | "
        "   PretextGraph -i {output.updated_map}.tmp4  -n TRF -o {output.updated_map}.tmp5  > {log.trf} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.windowmasker_track} | "
        "   PretextGraph -i {output.updated_map}.tmp5  -n windowmasker -o {output.updated_map}.tmp6  > {log.windowmasker} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.hifi_coverage_track} | "
        "   PretextGraph -i {output.updated_map}.tmp6  -n coverage -o {output.updated_map}  > {log.coverage} 2>&1; "
        " rm -r {output.updated_map}.tmp* > {log.rm} 2>&1; "