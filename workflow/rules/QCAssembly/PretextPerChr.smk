#ruleorder: pretextmap > pretextsnapshot
#localrules: get_candidate_chr_from_painted_agp
localrules: get_precurated_scaffold_order_from_agp
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
        queue=config["queue"]["cpu"]["name"],
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
        queue=config["queue"]["cpu"]["name"],
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
if candidate_agp_filename:
    rule get_precurated_scaffold_order_from_agp:
        input:
            #bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
            fasta="{fasta_dir}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.fasta",
            candidate_agp=candidate_agp_filename,
            log_dir="{fasta_dir}/{merged_haplotype}/log/"
        output:
            precurated_scaffold_orderlist="{fasta_dir}/{merged_haplotype, combined|reordered}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.precurated.orderlist",
        log:
            grep1="{fasta_dir}/{merged_haplotype}/log/get_precurated_scaffold_order_from_agp.{genome_prefix}.{assembly_stage}.{merged_haplotype}.grep1.log",
            grep2="{fasta_dir}/{merged_haplotype}/log/get_precurated_scaffold_order_from_agp.{genome_prefix}.{assembly_stage}.{merged_haplotype}.grep2.log",
            cut="{fasta_dir}/{merged_haplotype}/log/get_precurated_scaffold_order_from_agp.{genome_prefix}.{assembly_stage}.{merged_haplotype}.cut.log",
            cluster_log="{fasta_dir}/{merged_haplotype}/log/get_precurated_scaffold_order_from_agp.{genome_prefix}.{assembly_stage}.{merged_haplotype}.cluster.log",
            cluster_err="{fasta_dir}/{merged_haplotype}/log/get_precurated_scaffold_order_from_agp.{genome_prefix}.{assembly_stage}.{merged_haplotype}.cluster.err"
        benchmark:
            "{fasta_dir}/{merged_haplotype}/log/get_precurated_scaffold_order_from_agp.{genome_prefix}.{assembly_stage}.{merged_haplotype}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"]["name"],
            node_options=parse_node_list("get_precurated_scaffold_order_from_agp"),
            cpus=parameters["threads"]["get_precurated_scaffold_order_from_agp"] ,
            time=parameters["time"]["get_precurated_scaffold_order_from_agp"],
            mem=parameters["memory_mb"]["get_precurated_scaffold_order_from_agp"]
        threads: parameters["threads"]["get_precurated_scaffold_order_from_agp"]

        shell:
            " grep -vP \"^#\" {input.candidate_agp} 2>{log.grep1} | "
            " grep -P \"\\tW\\t\" 2>{log.grep2} | "
            " cut -f 6 > {output.precurated_scaffold_orderlist} 2>{log.cut}; "

    rule get_precurated_fasta:
        input:
            #bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
            fasta="{fasta_dir}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.fasta",
            precurated_scaffold_orderlist="{fasta_dir}/{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.precurated.orderlist",
            log_dir="{fasta_dir}/{merged_haplotype}/log/"
        output:
            precurated_fasta="{fasta_dir}/{merged_haplotype, combined|reordered}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.precurated.fasta"
        log:
            log="{fasta_dir}/{merged_haplotype}/log/get_precurated_fasta.{genome_prefix}.{assembly_stage}.{merged_haplotype}.log",
            cluster_log="{fasta_dir}/{merged_haplotype}/log/get_precurated_fasta.{genome_prefix}.{assembly_stage}.{merged_haplotype}.cluster.log",
            cluster_err="{fasta_dir}/{merged_haplotype}/log/get_precurated_fasta.{genome_prefix}.{assembly_stage}.{merged_haplotype}.cluster.err"
        benchmark:
            "{fasta_dir}/{merged_haplotype}/log/get_precurated_fasta.{genome_prefix}.{assembly_stage}.{merged_haplotype}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"]["name"],
            node_options=parse_node_list("get_precurated_fasta"),
            cpus=parameters["threads"]["get_precurated_fasta"] ,
            time=parameters["time"]["get_precurated_fasta"],
            mem=parameters["memory_mb"]["get_precurated_fasta"]
        threads: parameters["threads"]["get_precurated_fasta"]

        shell:
            " workflow/scripts/sequence/reorder_sequences.py -i {input.fasta} -b orderlist "
            "         -r {input.precurated_scaffold_orderlist} -o {output.precurated_fasta} > {log.log} 2>&1; "

    rule get_precurated_bam:
        input:
            bam="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{merged_haplotype}.rmdup.bam",
            bam_index="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{merged_haplotype}.rmdup.bam.csi",
            precurated_scaffold_orderlist="{fasta_dir}/{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.precurated.orderlist",
            #precurated_fasta="{fasta_dir}/{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.precurated.fasta",
            #precurated_fasta_dict="{fasta_dir}/{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.precurated.dict",
            log_dir="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/"
        output:
            intermediate_bam=temp("{fasta_dir}/{merged_haplotype, combined|reordered}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.rmdup.intermediate.bam"),
            precurated_bam="{fasta_dir}/{merged_haplotype, combined|reordered}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.rmdup.precurated.bam"
        params:
            sort_threads=parameters["threads"]["samtools_sort"],
            sort_per_thread=parameters["memory_mb"]["samtools_sort_per_thread"]
        log:
            view_header="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.view_header.log",
            reorder="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.reorder.log",
            view_records="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.view_records.log",
            compress="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.compress.log",
            sort="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.sort.log",
            cluster_log="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.cluster.log",
            cluster_err="{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_bam.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.cluster.err"
        benchmark:
            "{fasta_dir}/{merged_haplotype}/alignment/{phasing_kmer_length}/log/get_precurated_fasta.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{phasing_kmer_length}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"]["name"],
            node_options=parse_node_list("get_precurated_bam"),
            cpus=parameters["threads"]["samtools_sort"] + parameters["threads"]["get_precurated_bam"],
            time=parameters["time"]["get_precurated_bam"],
            mem=parameters["memory_mb"]["get_precurated_bam"] + parameters["memory_mb"]["samtools_sort_per_thread"] * parameters["threads"]["samtools_sort"]

        threads: parameters["threads"]["samtools_sort"] + parameters["threads"]["get_precurated_bam"]

        shell:
            " TMP_DIR=`dirname {output.precurated_bam}`; "
            " (samtools view -H {input.bam} 2>{log.view_header} | "
            "      workflow/scripts/alignment/reorder_scaffolds_in_sam_header.py -r {input.precurated_scaffold_orderlist} 2>{log.reorder}; "
            "      samtools view {input.bam} 2>{log.view_records}) | "
            "    samtools view -b -@ {params.sort_threads} > {output.intermediate_bam} 2>{log.compress}; "
            "  samtools sort -T ${{TMP_DIR}}/tmp.precurated.sort -@ {params.sort_threads} -m {params.sort_per_thread}M "
            "  -o {output.precurated_bam} {output.intermediate_bam} > {log.sort} 2>&1; "


    rule pretextmap_chr:
        input:
            precurated_bam="{bam_dir}/{bam_prefix}.rmdup.precurated.bam",
            precurated_bam_index="{bam_dir}/{bam_prefix}.rmdup.precurated.bam.csi",
            #bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
            candidate_chr_black_list=output_dict["data"] / "candidate_chr/candidate.{candidate_chr_id}.pretext.blacklist",
            log_dir="{bam_dir}/per_chr/log/"
        output:
            map="{bam_dir}/per_chr/{bam_prefix}.{candidate_chr_id}.rmdup.precurated.mapq{mapq, [0-9]+}.{res, default|high_res}.pretext"

            #filtered_out=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/per_chr/{phasing_kmer_length, [^.]+}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.{candidate_chr_id}.filtered_out.ids",
        params:
            resolution = lambda wildcards: " --highRes" if wildcards.res == "high_res" else "",
            sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
            sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
        log:
            view="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.view.log",
            cat="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.cat.log",
            tr="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.tr.log",
            sed="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.sed.log",
            cd="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.cd.log",
            map="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.map.log",
            echo="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.echo.log",
            cluster_log="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.cluster.log",
            cluster_err="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.cluster.err"
        benchmark:
            "{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{res}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"]["name"],
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
            " cd `dirname {input.precurated_bam}` > ${{CD_LOG}} 2>&1; echo $? >> ${{ECHO_LOG}} 2>&1; "
            " echo 'Creating pretext map...' >> ${{ECHO_LOG}} 2>&1; "
            " echo \"samtools view -@ 4 -F0x400 -h `basename {input.precurated_bam}` 2>${{VIEW_LOG}} | PretextMap -o per_chr/`basename {output.map}` {params.sortby} {params.sortorder} --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1 \" >> ${{ECHO_LOG}} 2>&1; "
            " samtools view -@ 4 -F0x400 -h `basename {input.precurated_bam}` 2>${{VIEW_LOG}} | "
            " PretextMap -o per_chr/`basename {output.map}` {params.sortby} {params.sortorder} "
            "            --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1; "
            " echo 'Creating pretext map finished...' >> ${{ECHO_LOG}} 2>&1; "
            " echo $? >> ${{ECHO_LOG}} 2>&1; "

    use rule pretext_inject_tracks as pretext_inject_tracks_per_chr with:
        input:
            map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.{phasing_kmer_length}.{subset}.rmdup.precurated.mapq{mapq}.{res}.pretext",
            filtered_out=output_dict["data"] / "candidate_chr/candidate.{subset}.pretext.blacklist",
            gap_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gap.track.bedgraph",
            canonical_telomere_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical.telomere.pretext.bedgraph",
            non_canonical_telomere_track=out_dir_path/ "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical.telomere.pretext.bedgraph",
            canonical_telomere_tidk_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_tidk.telomere.pretext.bedgraph",
            non_canonical_telomere_tidk_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_tidk.telomere.pretext.bedgraph",
            gc_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win10000.step1000.track.bedgraph" if not config["skip_pretext_10k_1k_tracks"] else [],
            gc_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win100000.step10000.track.bedgraph",
            trf_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win10000.step1000.track.bedgraph" if (not config["skip_trf"]) and (not config["skip_pretext_10k_1k_tracks"]) else [],
            trf_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win100000.step10000.track.bedgraph" if not config["skip_trf"] else [],
            windowmasker_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win10000.step1000.track.bedgraph"  if not config["skip_pretext_10k_1k_tracks"] else [],
            windowmasker_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win100000.step10000.track.bedgraph",
            all_hifi_coverage_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win10000.step1000.track.bedgraph"  if (not config["skip_pretext_10k_1k_tracks"]) and (not config["skip_pretext_coverage_tracks"]) and ("hifi" in data_types) else [],
            all_hifi_coverage_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win100000.step10000.track.bedgraph" if (not config["skip_pretext_coverage_tracks"]) and ("hifi" in data_types) else [],
            all_hifi_coverage_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.hifi_all_nodup_reads_mean_coverage.win1000000.step100000.track.bedgraph" if (not config["skip_pretext_1000k_100k_tracks"]) and (not config["skip_pretext_coverage_tracks"]) and ("hifi" in data_types) else [],
            all_illumina_coverage_10k_1k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.illumina_all_nodup_reads_mean_coverage.win10000.step1000.track.bedgraph"  if (not config["skip_pretext_10k_1k_tracks"]) and (not config["skip_pretext_coverage_tracks"]) and ("illumina" in data_types) else [],
            all_illumina_coverage_100k_10k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.illumina_all_nodup_reads_mean_coverage.win100000.step10000.track.bedgraph" if (not config["skip_pretext_coverage_tracks"]) and ("illumina" in data_types) else [],
            all_illumina_coverage_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.illumina_all_nodup_reads_mean_coverage.win1000000.step100000.track.bedgraph" if (not config["skip_pretext_1000k_100k_tracks"]) and (not config["skip_pretext_coverage_tracks"]) and ("illumina" in data_types) else [],
            gc_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gc.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
            trf_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.trf.win1000000.step100000.track.bedgraph" if (not config["skip_trf"]) and (not config["skip_pretext_1000k_100k_tracks"]) else [],
            windowmasker_1000k_100k_track=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.windowmasker.win1000000.step100000.track.bedgraph" if not config["skip_pretext_1000k_100k_tracks"] else [],
            busco_tracks=expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.{busco_type}.track.bedgraph",
                                busco_lineage=config["busco_lineage_list"],
                                busco_type=["single_copy", "duplicated", "fragmented",],
                                allow_missing=True) if not config["skip_busco"] else[],
            purge_dups_tracks=lambda wildcards: expand(out_dir_path / ("%s/%s/assembly_qc/tracks/%s.%s.%s/%s.%s.%s.purge_dups.{datatype}.{artefact_type}.track.bedgraph" % (wildcards.assembly_stage,
                                                                                                                                                                               wildcards.parameters,
                                                                                                                                                                               wildcards.genome_prefix,
                                                                                                                                                                               wildcards.assembly_stage,
                                                                                                                                                                               wildcards.haplotype,
                                                                                                                                                                               wildcards.genome_prefix,
                                                                                                                                                                               wildcards.assembly_stage,
                                                                                                                                                                               wildcards.haplotype)),
                                    datatype=set(stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"]["purge_dups_qc_datatypes"]) & set(data_types),
                                    artefact_type=["junk", "ovlp", "haplotig", "repeat", "highcov"],
                                    allow_missing=True) if not config["skip_purge_dups_qc"] else [],
        output:
            #tmp_gap_track = temp(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, combined|reordered}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.{phasing_kmer_length}.{candidate_chr_id}.precurated.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.gap.track"),
            #gap_track_tmp="{bam_dir}/per_chr/{fasta_prefix, [^/]+combined|[^/]+reordered}.{phasing_kmer_length, [^.]+}.{candidate_chr_id}.precurated.rmdup.mapq{mapq, [0-9]+}.{res, default|high_res}.gap.track",
            updated_map=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^./]+}/alignment/{phasing_kmer_length, [^.]+}/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.{phasing_kmer_length}.{subset}.rmdup.precurated.mapq{mapq, [0-9]+}.{res, default|high_res}.tracks.pretext",
            #updated_map="{bam_dir}/per_chr/{fasta_prefix, [^/]+combined|[^/]+reordered}.{phasing_kmer_length}.{candidate_chr_id}.rmdup.precurated.mapq{mapq, [0-9]+}.{res, default|high_res}.tracks.pretext",
        log:
            preprocessing=output_dict["log"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.preprocessing.log",
            injection=output_dict["log"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.injection.log",
            cluster_log=output_dict["cluster_log"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "pretext_inject_tracks.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{subset}.{mapq}.{res}.benchmark.txt"
