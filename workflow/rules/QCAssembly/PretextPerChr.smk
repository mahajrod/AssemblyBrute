#ruleorder: pretextmap > pretextsnapshot
#localrules: get_candidate_chr_from_painted_agp

"""
checkpoint get_candidate_chr_from_painted_agp: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        painted_agp=candidate_agp_filename
    output:
        out_dir=directory(output_dict["data"] / "candidate_chr/")
    log:
        rm=output_dict["log"]  / "get_candidate_chr_from_painted_agp.view.log",
        mkdir=output_dict["log"] / "get_candidate_chr_from_painted_agp.log",
        map=output_dict["log"]  / "get_candidate_chr_from_painted_agp.log",
        cluster_log=output_dict["cluster_log"] / "get_candidate_chr_from_painted_agp.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_candidate_chr_from_painted_agp.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_candidate_chr_from_painted_agp.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_candidate_chr_from_painted_agp"),
        cpus=parameters["threads"]["get_candidate_chr_from_painted_agp"] ,
        time=parameters["time"]["get_candidate_chr_from_painted_agp"],
        mem=parameters["memory_mb"]["get_candidate_chr_from_painted_agp"]
    threads: parameters["threads"]["get_candidate_chr_from_painted_agp"]

    shell:
        " rm -rf {output.out_dir} > {log.rm} 2>&1; "
        " mkdir -p {output.out_dir} > {log.mkdir} 2>&1; "
        " workflow/scripts/curation/extract_components_of_painted_scaffolds_from_agp.py -i {input.painted_agp} "
        "       -p {output.out_dir}/candidate "
"""
"""
rule pretextmap_chr: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        len=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len"
    output:
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.filtered_out.ids",
    params:
        resolution=lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
        #max_len=lambda wildcards: parameters["tool_options"]["pretextmap"]["subsets"][wildcards.subset]["max_len"],
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
    log:
        view=output_dict["log"]  / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{candidate_chr_id}.view.log",
        awk=output_dict["log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{candidate_chr_id}.awk.log",
        map=output_dict["log"]  / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{candidate_chr_id}.map.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{candidate_chr_id}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{candidate_chr_id}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{candidate_chr_id}.benchmark.txt"
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
        " pwd "


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

"""

rule pretextmap_chr:
    input:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        candidate_chr_black_list=output_dict["data"] / "candidate_chr/candidate.{candidate_chr_id}.pretext.blacklist"
    output:
        map=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^./]+}/per_chr/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.pretext",

        #filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/per_chr/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.{candidate_chr_id}.filtered_out.ids",
    params:
        resolution = lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
    log:
        view=output_dict["log"]  / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.view.log",
        cat=output_dict["log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cat.log",
        tr=output_dict["log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.tr.log",
        sed=output_dict["log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.sed.log",
        cd=output_dict["log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cd.log",
        map=output_dict["log"]  / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.map.log",
        echo=output_dict["log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.echo.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pretextmap_chr"),
        cpus=parameters["threads"]["pretextmap_chr"] ,
        time=parameters["time"]["pretextmap_chr"],
        mem=parameters["memory_mb"]["pretextmap_chr"]
    threads: parameters["threads"]["pretextmap_chr"]

    shell:
        " MAP_LOG=`realpath -s -m {log.map}` ; "
        " VIEW_LOG=`realpath -s -m {log.view}` ; "
        " ECHO_LOG=`realpath -s -m {log.echo}` ; "
        " CD_LOG=`realpath -s -m {log.cd}` ; "
        " if [[ -s {input.candidate_chr_black_list} ]]; "
        "   then "
        "       echo 'Blacklist is not empty...' > ${{ECHO_LOG}} 2>&1; "
        "       FILTER_OUT=`cat {input.candidate_chr_black_list} | tr '\\n' ',' | sed 's/,\+$//' `; echo ${{FILTER_OUT}} >> ${{ECHO_LOG}} 2>&1; "
        "       FILTER_OUT=\" --filterExclude ${{FILTER_OUT}}\"; echo $? >> ${{ECHO_LOG}} 2>&1; "
        "   else "
        "       echo 'Blacklist is empty...' > ${{ECHO_LOG}} 2>&1; "
        "       FILTER_OUT=''; echo $? >> ${{ECHO_LOG}} 2>&1;"
        "   fi; "
        " echo 'Entering workdir...' >> ${{ECHO_LOG}} 2>&1; "
        " cd `dirname {input.bam}` > ${{CD_LOG}} 2>&1; echo $? >> ${{ECHO_LOG}} 2>&1; "
        " echo 'Creating pretext map...' >> ${{ECHO_LOG}} 2>&1; "
        " echo \"samtools view -@ 4 -F0x400 -h `basename {input.bam}` 2>${{VIEW_LOG}} | PretextMap -o per_chr/`basename {output.map}` {params.sortby} {params.sortorder} --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1 \" >> ${{ECHO_LOG}} 2>&1; "
        " samtools view -@ 4 -F0x400 -h `basename {input.bam}` 2>${{VIEW_LOG}} | "
        " PretextMap -o per_chr/`basename {output.map}` {params.sortby} {params.sortorder} "
        "            --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1; "
        " echo 'Creating pretext map finished...' >> ${{ECHO_LOG}} 2>&1; "
        " echo $? >> ${{ECHO_LOG}} 2>&1; "

rule pretext_inject_tracks_per_chr:
    input:
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.rmdup.mapq{mapq}.{res}.pretext",
        candidate_chr_black_list=output_dict["data"] / "candidate_chr/candidate.{candidate_chr_id}.pretext.blacklist",
        gap_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gap.track.bedgraph",
        canonical_telomere_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical.telomere.pretext.bedgraph",
        non_canonical_telomere_track=out_dir_path/ "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical.telomere.pretext.bedgraph",
        gc_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win10000.step1000.track.bedgraph" if not config["skip_pretext_10k_1k_tracks"] else [],
        gc_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win100000.step10000.track.bedgraph",
        trf_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win10000.step1000.track.bedgraph" if (not config["skip_trf"]) and (not config["skip_pretext_10k_1k_tracks"]) else [],
        trf_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win100000.step10000.track.bedgraph" if not config["skip_trf"] else [],
        windowmasker_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win10000.step1000.track.bedgraph" if not config["skip_pretext_10k_1k_tracks"] else [],
        windowmasker_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win100000.step10000.track.bedgraph",
        all_hifi_coverage_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win10000.step1000.track.bedgraph" if not config["skip_pretext_10k_1k_tracks"] else [],
        all_hifi_coverage_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win100000.step10000.track.bedgraph",
        all_hifi_coverage_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
        gc_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
        trf_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win1000000.step100000.track.bedgraph" if (not config["skip_trf"]) and (not config["skip_pretext_1000k_100k_tracks"]) else [],
        windowmasker_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
    output:
        updated_map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.tracks.pretext"
    params:
        min_mapq=parameters["tool_options"]["pretextmap"]["mapq"],
        skip_trf=config["skip_trf"],
        skip_10k_1k_tracks=config["skip_pretext_10k_1k_tracks"],
        skip_1000k_100k_tracks=config["skip_pretext_1000k_100k_tracks"]
    log:
        gap=output_dict["log"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.gap.log",
        can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.can_tel.log",
        non_can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.non_can_tel.log",
        gc=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}..gc.log",
        trf=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.trf.log",
        windowmasker=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.windowmasker.log",
        coverage=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.coverage.log",
        awk=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.awk.log",
        rm=output_dict["log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.rm.log",
        cluster_log=output_dict["cluster_log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.benchmark.txt"
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
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.gap_track}  | "
        "   PretextGraph -i {output.updated_map} -n gap > {log.gap} 2>&1; "
        " if [[ -s {input.canonical_telomere_track} ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.canonical_telomere_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map} -n canonical.telomere > {log.can_tel} 2>&1;"
        "   fi; "
        " if [[ -s {input.non_canonical_telomere_track} ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.non_canonical_telomere_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map} -n noncanonical.telomere > {log.non_can_tel} 2>&1;"
        "   fi;  "
        " if [[ '{params.skip_10k_1k_tracks}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.gc_10k_1k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map} -n GC_10k_1k.repeat_density  > {log.gc} 2>&1; "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.windowmasker_10k_1k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map}  -n windowmasker_10k_1k.repeat_density > {log.windowmasker} 2>&1; "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.all_hifi_coverage_10k_1k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map}  -n hifi_all_10k_1k.coverage  > {log.coverage} 2>&1;"
        "   if [[ '{params.skip_trf}' != 'True' ]]; "
        "       then "
        "       workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.trf_10k_1k_track} | "
        "           awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "           PretextGraph -i {output.updated_map}  -n TRF_10k_1k.repeat_density > {log.trf} 2>&1; "
        "       fi; "
        " fi; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.gc_100k_10k_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map} -n GC_100k_10k.repeat_density  > {log.gc} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.windowmasker_100k_10k_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map}  -n windowmasker_100k_10k.repeat_density > {log.windowmasker} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.all_hifi_coverage_100k_10k_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "   PretextGraph -i {output.updated_map}  -n hifi_all_100k_10k.coverage  > {log.coverage} 2>&1; "
        " if [[ '{params.skip_trf}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.trf_100k_10k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map}  -n TRF_100k_10k.repeat_density > {log.trf} 2>&1; "
        " fi; "
        " if [[ '{params.skip_1000k_100k_tracks}' != 'True' ]]; "
        "   then "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.gc_1000k_100k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map} -n GC_1000k_100k.repeat_density  > {log.gc} 2>&1; "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.windowmasker_1000k_100k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "       PretextGraph -i {output.updated_map}  -n windowmasker_1000k_100k.repeat_density > {log.windowmasker} 2>&1; "
        "   workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.all_hifi_coverage_1000k_100k_track} | "
        "       awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "       PretextGraph -i {output.updated_map}  -n hifi_all_1000k_100k.coverage  > {log.coverage} 2>&1;"
        "   if [[ '{params.skip_trf}' != 'True' ]]; "
        "       then "
        "       workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.candidate_chr_black_list} -i {input.trf_1000k_100k_track} | "
        "           awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "           PretextGraph -i {output.updated_map}  -n TRF_1000k_100k.repeat_density > {log.trf} 2>&1; "
        "       fi; "
        " fi; "
