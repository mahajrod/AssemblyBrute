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

rule pretextmap_per_chr: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        candidate_chr_black_list=output_dict["data"] / "candidate_chr/candidate.{candidate_chr_id}.pretext.blacklist"
    output:
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.{candidate_chr_id}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.{candidate_chr_id}.filtered_out.ids",
    params:
        resolution = lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
        sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
        sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
    log:
        view=output_dict["log"]  / "pretextmap_per_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.view.log",
        awk=output_dict["log"] / "pretextmap_per_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.awk.log",
        map=output_dict["log"]  / "pretextmap_per_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}..map.log",
        cluster_log=output_dict["cluster_log"] / "pretextmap_per_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretextmap_per_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretextmap_per_chr.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.benchmark.txt"
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
        " if [[ -s {input.candidate_chr_black_list} ]]; "
        "   then "
        "       FILTER_OUT=' --filterExclude '; "
        "       FILTER_OUT=\"${{FILTER_OUT}} `cat {input.candidate_chr_black_list} | tr '\\n' ',' | sed 's/,\+$//'` | sed 's/,/, /g' \"; "
        "   else "
        "       FILTER_OUT=''; "
        "   fi; " 
        " cd `dirname {input.bam}`; "
        " samtools view -@4 -F0x400 -h `basename {input.bam}` 2>${{VIEW_LOG}} | "
        " PretextMap -o per_chr/`basename {output.map}` {params.sortby} {params.sortorder} "
        "            --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1"

rule pretext_inject_tracks_per_chr:
    input:
        #bam=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bam"  % config["genome_name"]),
        map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq}.{res}.{candidate_chr_id}.pretext",
        filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq}.{res}.{candidate_chr_id}.filtered_out.ids",
        gap_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gap.track.bedgraph",
        canonical_telomere_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical.telomere.pretext.bedgraph",
        non_canonical_telomere_track=out_dir_path/ "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical.telomere.pretext.bedgraph",
        gc_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win{window}.step{step}.track.bedgraph",
        trf_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win{window}.step{step}.track.bedgraph",
        windowmasker_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win{window}.step{step}.track.bedgraph",
        all_hifi_coverage_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win{window}.step{step}.track.bedgraph",
        default_hifi_coverage_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_default_mean_coverage.win{window}.step{step}.track.bedgraph"

    output:
        updated_map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.tracks.win_{window, [0-9]+}.{step, [0-9]+}.{candidate_chr_id}.pretext"
    params:
        min_mapq=parameters["tool_options"]["pretextmap"]["mapq"]
    log:
        gap=output_dict["log"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.gap.log",
        can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.can_tel.log",
        non_can_tel=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.non_can_tel.log",
        gc=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.gc.log",
        trf=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.trf.log",
        windowmasker=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.windowmasker.log",
        coverage=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.coverage.log",
        awk=output_dict["log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.awk.log",
        rm=output_dict["log"] / "pretextmap.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.rm.log",
        cluster_log=output_dict["cluster_log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{candidate_chr_id}.{mapq}.{res}.{window}.{step}.benchmark.txt"
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
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.gc_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map} -n GC.repeat_density  > {log.gc} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.trf_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map}  -n TRF.repeat_density > {log.trf} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.windowmasker_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4 }}' | "
        "   PretextGraph -i {output.updated_map}  -n windowmasker.repeat_density > {log.windowmasker} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.all_hifi_coverage_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "   PretextGraph -i {output.updated_map}  -n hifi_all.coverage  > {log.coverage} 2>&1; "
        " workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} -i {input.default_hifi_coverage_track} | "
        "   awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "   PretextGraph -i {output.updated_map}  -n hifi_default.coverage  > {log.coverage} 2>&1; "

"""
def aggregate_per_chr_maps_input(wildcards):
    checkpoint_output = checkpoints.get_candidate_chr_from_painted_agp.get(**wildcards).output[0]
    print(checkpoint_output)
    print(os.path.join(checkpoint_output,"/candidate.{candidate_chr_id}.pretext.blacklist"))
    print(glob_wildcards(os.path.join(checkpoint_output,"/candidate.{candidate_chr_id}.pretext.blacklist")).candidate_chr_id)
    return expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/per_chr/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq}.{res}.tracks.win_{window}.{step}.{candidate_chr_id}.pretext",
                  candidate_chr=glob_wildcards(os.path.join(checkpoint_output,"/candidate.{candidate_chr_id}.pretext.blacklist")).candidate_chr_id,
                  allow_missing=True)

rule aggregate_per_chr_maps: # #Pretext-map probably doesn't support long file names!!!!!!!!!!!
    input:
        aggregate_per_chr_maps_input
    output:
        out_dir=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^./]+}/per_chr/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq, [0-9]+}.{res, [^/]+}.tracks.win_{window, [0-9]+}.{step, [0-9]+}.maps.list"
    log:
        log=output_dict["log"]  / "aggregate_per_chr_maps.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{window}.{step}.log",
        cluster_log=output_dict["cluster_log"] / "aggregate_per_chr_maps.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "aggregate_per_chr_maps.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_candidate_chr_from_painted_agp.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{mapq}.{res}.{window}.{step}.benchmark.txt"
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
        " cat {input} > {output} 2>{log.log}; "
"""