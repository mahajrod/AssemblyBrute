#ruleorder: pretextmap > pretextsnapshot
rule pretextmap: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        len=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len"
    output:
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.filtered_out.ids"
    params:
        resolution=lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
        max_len=lambda wildcards: parameters["tool_options"]["pretextmap"]["subsets"][wildcards.subset]["max_len"],
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
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
        "       FILTER_OUT=\"${{FILTER_OUT}} `cat {output.filtered_out} | tr '\\n' ',' | sed 's/,\+$//'` \"; "
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
    output:
        dir=directory(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.mapq{mapq, [0-9]+}.default.{resolution, [0-9]+}.{ext}"),
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
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.filtered_out.ids",
        gap_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gap.track.bedgraph",
        canonical_telomere_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical.telomere.pretext.bedgraph",
        non_canonical_telomere_track=out_dir_path/ "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical.telomere.pretext.bedgraph",
        gc_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win10000.step1000.track.bedgraph" if not config["skip_pretext_10k_1k_tracks"] else [],
        gc_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win100000.step10000.track.bedgraph",
        trf_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win10000.step1000.track.bedgraph" if (not config["skip_trf"]) and (not config["skip_pretext_10k_1k_tracks"]) else [],
        trf_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win100000.step10000.track.bedgraph" if not config["skip_trf"] else [],
        windowmasker_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win10000.step1000.track.bedgraph"  if not config["skip_pretext_10k_1k_tracks"] else [],
        windowmasker_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win100000.step10000.track.bedgraph",
        all_hifi_coverage_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win10000.step1000.track.bedgraph"  if (not config["skip_pretext_10k_1k_tracks"]) and (not config["skip_pretext_coverage_tracks"]) else [],
        all_hifi_coverage_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win100000.step10000.track.bedgraph" if not config["skip_pretext_coverage_tracks"] else [],
        all_hifi_coverage_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win1000000.step100000.track.bedgraph" if (not config["skip_pretext_1000k_100k_tracks"]) and (not config["skip_pretext_coverage_tracks"]) else [],
        gc_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
        trf_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win1000000.step100000.track.bedgraph" if (not config["skip_trf"]) and (not config["skip_pretext_1000k_100k_tracks"]) else [],
        windowmasker_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
    output:
        updated_map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{subset}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.tracks.pretext",
    params:
        min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        skip_trf=config["skip_trf"],
        skip_10k_1k_tracks=config["skip_pretext_10k_1k_tracks"],
        skip_1000k_100k_tracks=config["skip_pretext_1000k_100k_tracks"],
        skip_pretext_coverage_tracks=config["skip_pretext_coverage_tracks"]
    log:
        gap=output_dict["log"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.gap.log",
        can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.can_tel.log",
        non_can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.non_can_tel.log",
        gc=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}..gc.log",
        trf=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.trf.log",
        windowmasker=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.windowmasker.log",
        coverage=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.coverage.log",
        awk=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.awk.log",
        rm=output_dict["log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.rm.log",
        cluster_log=output_dict["cluster_log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.benchmark.txt"
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
        " cp -f {input.map} {output.updated_map}; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gap_track}  | "
        "   PretextGraph -i {output.updated_map} -n gap > {log.gap} 2>&1; "
        " if [[ -s {input.canonical_telomere_track} ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.canonical_telomere_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map} -n canonical.telomere > {log.can_tel} 2>&1;"
        "   fi; "
        " if [[ -s {input.non_canonical_telomere_track} ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.non_canonical_telomere_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map} -n noncanonical.telomere > {log.non_can_tel} 2>&1;"
        "   fi;  "
        " if [[ '{params.skip_10k_1k_tracks}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.windowmasker_10k_1k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map}  -n windowmasker_10k_1k.repeat_density > {log.windowmasker} 2>&1; "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gc_10k_1k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map} -n GC_10k_1k.repeat_density  > {log.gc} 2>&1; "
        "   if [[ '{params.skip_pretext_coverage_tracks}' != 'True' ]]; "
        "       then "
        "       workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.all_hifi_coverage_10k_1k_track} | "
        "           awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "           PretextGraph -i {output.updated_map}  -n hifi_all_10k_1k.coverage  > {log.coverage} 2>&1; "
        "       fi; "
        "   if [[ '{params.skip_trf}' != 'True' ]]; "
        "       then "
        "       workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.trf_10k_1k_track} | "
        "           awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "           PretextGraph -i {output.updated_map}  -n TRF_10k_1k.repeat_density > {log.trf} 2>&1; "
        "       fi; "
        "   fi; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gc_100k_10k_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map} -n GC_100k_10k.repeat_density  > {log.gc} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.windowmasker_100k_10k_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map}  -n windowmasker_100k_10k.repeat_density > {log.windowmasker} 2>&1; "
        " if [[ '{params.skip_pretext_coverage_tracks}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.all_hifi_coverage_100k_10k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map}  -n hifi_all_100k_10k.coverage  > {log.coverage} 2>&1; "
        "   fi; "
        " if [[ '{params.skip_trf}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.trf_100k_10k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map}  -n TRF_100k_10k.repeat_density > {log.trf} 2>&1; "
        "   fi; "
        " if [[ '{params.skip_1000k_100k_tracks}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.windowmasker_1000k_100k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map}  -n windowmasker_1000k_100k.repeat_density > {log.windowmasker} 2>&1; "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gc_1000k_100k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map} -n GC_1000k_100k.repeat_density  > {log.gc} 2>&1; "
        "   if [[ '{params.skip_pretext_coverage_tracks}' != 'True' ]]; "
        "       then "
        "       workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.all_hifi_coverage_1000k_100k_track} | "
        "           awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "           PretextGraph -i {output.updated_map}  -n hifi_all_1000k_100k.coverage  > {log.coverage} 2>&1; "
        "   fi; "
        "   if [[ '{params.skip_trf}' != 'True' ]]; "
        "       then "
        "       workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.trf_1000k_100k_track} | "
        "           awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "           PretextGraph -i {output.updated_map}  -n TRF_1000k_100k.repeat_density > {log.trf} 2>&1; "
        "       fi; "
        "   fi; "