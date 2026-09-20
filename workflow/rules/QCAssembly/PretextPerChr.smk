

if candidate_agp_filename:
    localrules: get_precurated_scaffold_order_from_agp
    rule get_precurated_scaffold_order_from_agp:
        input:
            fasta="{fasta_dir}/{fasta_prefix}.fasta",
            candidate_agp=candidate_agp_filename,
            log_dir=ancient("{fasta_dir}/log/")
        output:
            precurated_scaffold_orderlist="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.precurated.orderlist",
        log:
            grep1="{fasta_dir}/log/get_precurated_scaffold_order_from_agp.{fasta_prefix}.grep1.log",
            grep2="{fasta_dir}/log/get_precurated_scaffold_order_from_agp.{fasta_prefix}.grep2.log",
            cut="{fasta_dir}/log/get_precurated_scaffold_order_from_agp.{fasta_prefix}.cut.log",
            cluster_log="{fasta_dir}/log/get_precurated_scaffold_order_from_agp.{fasta_prefix}.cluster.log",
            cluster_err="{fasta_dir}/log/get_precurated_scaffold_order_from_agp.{fasta_prefix}.cluster.err"
        benchmark:
            "{fasta_dir}/log/get_precurated_scaffold_order_from_agp.{fasta_prefix}.benchmark.txt"
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
            fasta="{fasta_dir}/{fasta_prefix}.fasta",
            precurated_scaffold_orderlist="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.precurated.orderlist",
            log_dir=ancient("{fasta_dir}/log/")
        output:
            precurated_fasta="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.precurated.fasta"
        log:
            log="{fasta_dir}/log/get_precurated_fasta.{fasta_prefix}.log",
            cluster_log="{fasta_dir}/log/get_precurated_fasta.{fasta_prefix}.cluster.log",
            cluster_err="{fasta_dir}/log/get_precurated_fasta.{fasta_prefix}.cluster.err"
        benchmark:
            "{fasta_dir}/log/get_precurated_fasta.{fasta_prefix}.benchmark.txt"
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
            bam="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.bam",
            bam_index="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.bam.csi",
            precurated_scaffold_orderlist="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.precurated.orderlist",
            log_dir=ancient("{fasta_dir}/log/")
        output:
            intermediate_bam=temp("{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.intermediate.bam"),
            precurated_bam="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.precurated.bam"
        params:
            sort_threads=parameters["threads"]["samtools_sort"],
            sort_per_thread=parameters["memory_mb"]["samtools_sort_per_thread"]
        log:
            view_header="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.view_header.log",
            reorder="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.reorder.log",
            view_records="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.view_records.log",
            compress="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.compress.log",
            sort="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.sort.log",
            cluster_log="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.cluster.log",
            cluster_err="{fasta_dir}/log/get_precurated_bam.{fasta_prefix}.{phasing_kmer_length}.cluster.err"
        benchmark:
            "{fasta_dir}/log/get_precurated_fasta.{fasta_prefix}.{phasing_kmer_length}.benchmark.txt"
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
            "     workflow/scripts/alignment/reorder_scaffolds_in_sam_header.py -r {input.precurated_scaffold_orderlist} 2>{log.reorder}; "
            "     samtools view {input.bam} 2>{log.view_records}) | "
            "     samtools view -b -@ {params.sort_threads} > {output.intermediate_bam} 2>{log.compress}; "
            "  samtools sort -T ${{TMP_DIR}}/tmp.precurated.sort -@ {params.sort_threads} -m {params.sort_per_thread}M "
            "  -o {output.precurated_bam} {output.intermediate_bam} > {log.sort} 2>&1; "


    rule pretextmap_chr:
        input:
            precurated_bam="{bam_dir}/{bam_prefix}.rmdup.precurated.bam",
            precurated_bam_index="{bam_dir}/{bam_prefix}.rmdup.precurated.bam.csi",
            candidate_chr_black_list=config["out_dir"] / "data/candidate_chr/candidate.{candidate_chr_id}.pretext.blacklist",
            log_dir=ancient("{bam_dir}/per_chr/log/")
        output:
            map="{bam_dir}/per_chr/{bam_prefix}.{candidate_chr_id}.rmdup.precurated.mapq{mapq}.{pretext_res}.pretext"

        params:
            resolution = lambda wildcards: " --ultraRes" if wildcards.pretext_res == "ultra_res" else " --highRes" if wildcards.pretext_res == "high_res" else "",
            sortby=parse_option("sortby", parameters["tool_options"]["pretextmap"], " --sortby "),
            sortorder=parse_option("sortorder", parameters["tool_options"]["pretextmap"], " --sortorder "),
        log:
            view="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.view.log",
            cat="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.cat.log",
            tr="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.tr.log",
            sed="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.sed.log",
            cd="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.cd.log",
            map="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.map.log",
            echo="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.echo.log",
            cluster_log="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.cluster.log",
            cluster_err="{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.cluster.err"
        benchmark:
            "{bam_dir}/per_chr/log/pretextmap_chr.{bam_prefix}.{candidate_chr_id}.{mapq}.{pretext_res}.benchmark.txt"
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
            "     --mapq {wildcards.mapq} ${{FILTER_OUT}} {params.resolution} > ${{MAP_LOG}} 2>&1; "
            " echo 'Creating pretext map finished...' >> ${{ECHO_LOG}} 2>&1; "
            " echo $? >> ${{ECHO_LOG}} 2>&1; "

    use rule pretext_inject_tracks as pretext_inject_tracks_per_chr with:
        input:
            map="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/per_chr/{fasta_prefix}.{phasing_kmer_length}.{subset}.rmdup.precurated.mapq{mapq}.{pretext_res}.pretext",
            filtered_out=config["out_dir"] / "data/candidate_chr/candidate.{subset}.pretext.blacklist",
            track_info="{fasta_dir}/{fasta_prefix}.pretext.track.info",
            log_dir="{fasta_dir}/log/",
        output:
            updated_map="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/per_chr/{fasta_prefix}.{phasing_kmer_length}.{subset}.rmdup.precurated.mapq{mapq}.{pretext_res}.tracks.pretext",
