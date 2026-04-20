import pandas as pd

if (sum(list(pd.Series(["hic_scaffolding",
                        "gap_closing",
                        "draft_qc",
                        "contig"],).isin(config["stage_list"]))) > 0) and ("hic" in data_types) : # TODO: remove this bypass in future
    if "purge_dups" in config["stage_list"]:
        ruleorder: combine_haplotypes > create_final_links_purge_dups
        ruleorder: combine_haplotypes > create_assembly_links_if_skipping_purge_dups

    if "hic_scaffolding" in config["stage_list"]:
        ruleorder: combine_haplotypes > yahs
    else:
        pass

    if "contig" in config["stage_list"]:
        ruleorder: combine_haplotypes > gfa2fasta
    #ruleorder: bam_merge_files_for_hic_map > bam_merge_files
    if config["other_tool_option_sets"]["mapping_pipeline"] != "juicer":
        pass
        #ruleorder: rmdup_for_hic_map > rmdup
        #ruleorder: bwa_map_for_hic_map > bam_merge_pairs

    rule combine_haplotypes:
        input:
            fasta_list=lambda wildcards: expand(out_dir_path / ("%s/%s/%s.%s.{haplotype}.fasta" % (wildcards.assembly_stage,
                                                                                                         wildcards.parameters,
                                                                                                         wildcards.genome_prefix,
                                                                                                         wildcards.assembly_stage)),
                                                 haplotype=stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"] ,
                                                 allow_missing=True)
        output:
            combined_fasta=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.combined.fasta",
        params:
            haplotype_list=lambda wildcards: stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"],
            out_dir=str(out_dir_path)
        log:
            log=output_dict["log"]  / "combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.log",
            cluster_log=output_dict["cluster_log"] / "combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"],
            node_options=parse_node_list("combine_haplotypes"),
            cpus=parameters["threads"]["combine_haplotypes"],
            time=parameters["time"]["combine_haplotypes"],
            mem=parameters["memory_mb"]["combine_haplotypes"]
        threads: parameters["threads"]["combine_haplotypes"]

        shell:
            " > {output.combined_fasta}; "
            " > {log.log}; "
            " for HAP in {params.haplotype_list};"
            " do "
            "       label_sequences.py -i {params.out_dir}/{wildcards.assembly_stage}/{wildcards.parameters}/{wildcards.genome_prefix}.{wildcards.assembly_stage}.${{HAP}}.fasta "
            "                          -l ${{HAP}} -s '.' >> {output.combined_fasta} 2>{log.log} ; "
            " done "

    hic_curation_haplotype = "reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"

    rule bam2bed_for_hic_map:
        input:
            bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.rmdup.bam",

        output:
            bed=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.nodup.bed"
        log:
            samtools=output_dict["log"]  / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.samtools.log",
            bam2bed=output_dict["log"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.bam2bed.log",
            sort=output_dict["log"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.sort.log",
            cluster_log=output_dict["cluster_log"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"],
            node_options=parse_node_list("bam2bed_for_hic_map"),
            cpus=parameters["threads"]["bam2bed_for_hic_map"],
            time=parameters["time"]["bam2bed_for_hic_map"],
            mem=parameters["memory_mb"]["bam2bed_for_hic_map"]
        threads: parameters["threads"]["bam2bed_for_hic_map"]

        shell:
            " TMP_DIR={output.bed}.tmp; "
            " mkdir -p ${{TMP_DIR}}; "
            " samtools view -@ 4 -u -F0x400 {input.bam} 2>{log.samtools} | bamToBed 2>{log.bam2bed} | "
            " sort -k4 --parallel=20 -S100G -T ${{TMP_DIR}} > {output.bed} 2>{log.sort}; "
            " rm -r ${{TMP_DIR}}; "

    rule bed2pairs_for_hic_map:
        input:
            bed=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.nodup.bed",
        output:
            pairs=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.nodup.pairs",
        log:
            paste=output_dict["log"]  / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.paste.log",
            awk1=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.awk1.log",
            awk2=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.awk2.log",
            tr=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.tr.log",
            sort=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.sort.log",
            cluster_log=output_dict["cluster_log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"],
            node_options=parse_node_list("bed2pairs"),
            cpus=parameters["threads"]["bed2pairs"],
            time=parameters["time"]["bed2pairs"],
            mem=parameters["memory_mb"]["bed2pairs"]
        threads: parameters["threads"]["bed2pairs"]

        shell:
            " TMP_DIR={output.pairs}.tmp; "
            " mkdir -p ${{TMP_DIR}}; "
            " paste -d '\t' - - < {input.bed} 2>{log.paste} | "
            " awk 'BEGIN {{FS=\"\t\"; OFS=\"\t\"}} "
            "      {{if ($1 > $7) {{print substr($4,1,length($4)-2),$12,$7,$8,\"16\",$6,$1,$2,\"8\",$11,$5}} "
            "      else {{print substr($4,1,length($4)-2),$6,$1,$2,\"8\",$12,$7,$8,\"16\",$5,$11}} }}' 2>{log.awk1} | "
            " tr '\-+' '01' 2>{log.tr} | sort -k3,3d -k7,7d --parallel=20 -S100G -T ${{TMP_DIR}} 2>{log.sort} | awk 'NF==11' > {output.pairs} 2>{log.awk2};"
            " rm -r ${{TMP_DIR}}; "

    rule create_higlass_track_from_bed: #
        input:
            fai=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.fasta.fai",
            pairs=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.hic_scaffolding.NA.{haplotype}.nodup.pairs"
        output:
            genome_higlass=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.higlass.genome",
            higlass_cool=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.hic_scaffolding.NA.{haplotype}.nodup.higlass.cool",
            higlass_mcool=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.hic_scaffolding.NA.{haplotype}.nodup.higlass.mcool",
        log:
            cut=output_dict["log"]  / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cut.log",
            sed1=output_dict["log"] / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.sed1.log",
            sort=output_dict["log"] / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.sort.log",
            cload=output_dict["log"]  / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cload.log",
            zoomify=output_dict["log"]  / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.zoomify.log",
            cluster_log=output_dict["cluster_log"] / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "create_higlass_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"],
            node_options=parse_node_list("create_higlass_track_from_bed"),
            cpus=parameters["threads"]["create_higlass_track"] ,
            time=parameters["time"]["create_higlass_track"],
            mem=parameters["memory_mb"]["create_higlass_track"]
        threads: parameters["threads"]["create_higlass_track"]
        shell:
            " cut -f1,2 {input.fai} 2>{log.cut} | sed 's/-/_/g' 2>{log.sed1} | sort -k2,2 -nr > {output.genome_higlass} 2>{log.sort}; "
            " cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 {output.genome_higlass}:1000 {input.pairs} {output.higlass_cool} 2>{log.cload}; "
            " cooler zoomify --resolutions 5000,10000,25000,50000,100000,150000,200000,300000,400000,500000,1000000,2500000 "
            " -o {output.higlass_mcool} {output.higlass_cool} 2>{log.zoomify}; "

    rule bwa_map_for_hic_map: #
        input:
            index=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.ann",
            reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
            forward_fastq=lambda wildcards: ancient(output_dict["data"] / "fastq/hic/{0}/{1}{2}{3}".format("filtered" if "hic" in config["filtered_data"] else "raw",
                                                                                                  wildcards.pairprefix,
                                                                                                  "_1" if "hic" in config["filtered_data"] else input_forward_suffix_dict["hic"],
                                                                                                  config["fastq_extension"])),
            reverse_fastq=lambda wildcards: ancient(output_dict["data"] / "fastq/hic/{0}/{1}{2}{3}".format("filtered" if "hic" in config["filtered_data"] else "raw",
                                                                                                   wildcards.pairprefix,
                                                                                                   "_2" if "hic" in config["filtered_data"] else input_reverse_suffix_dict["hic"],
                                                                                                   config["fastq_extension"])) ,
        output:
            #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
            bam=temp(out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{haplotype, combined|reordered}/alignment/{phasing_kmer_length, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix, [^/]+}.bwa.bam")
        params:
            id="{0}_hic".format(config["genome_prefix"]),
            bwa_tool=config["bwa_tool"]
        log:
            map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.map.log",
            sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.sort.log",
            #filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.filter.log",
            cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"],
            node_options=parse_node_list("bwa_map"),
            cpus=parameters["threads"]["bwa_map"] ,
            time=parameters["time"]["bwa_map"],
            mem=parameters["memory_mb"]["bwa_map"]
        threads: parameters["threads"]["bwa_map"]
        shell:
            " {params.bwa_tool} mem  -T 0 -SP5 -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
            " {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "
"""
    rule rmdup_for_hic_map:
        input:
            bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.bam"
        output:
            bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, combined|reordered}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
            #bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
        params:
            sort_threads=parameters["threads"]["samtools_sort"],
            collate_threads=parameters["threads"]["samtools_collate"],
            fixmate_threads=parameters["threads"]["samtools_fixmate"],
            markdup_threads=parameters["threads"]["samtools_markdup"],
            sort_per_thread=parameters["memory_mb"]["samtools_sort_per_thread"]
        log:
            collate=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.collate.log",
            fixmate=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.fixmate.log",
            sort=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort.log",
            markdup=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.markdup.log",
            cluster_log=output_dict["cluster_log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
            cluster_err=output_dict["cluster_error"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
        benchmark:
            output_dict["benchmark"]  / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"],
            node_options=parse_node_list("rmdup"),
            cpus=parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"],
            time=parameters["time"]["rmdup"],
            mem=50000 + parameters["memory_mb"]["samtools_collate"] + parameters["memory_mb"]["samtools_fixmate"] + parameters["memory_mb"]["samtools_markdup"] + parameters["memory_mb"]["samtools_sort_per_thread"] * parameters["threads"]["samtools_sort"]
        threads: parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"]
        shell:
            " TMP_DIR=`dirname {output.bam}`; "
            " samtools collate -T ${{TMP_DIR}}/tmp.collate  -@ {params.collate_threads}  -O {input.bam} 2>{log.collate} | "
            " samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
            " samtools sort -T ${{TMP_DIR}}/tmp.sort -@ {params.sort_threads} -m {params.sort_per_thread}M 2>{log.sort} | "
            " samtools markdup -@ {params.markdup_threads} - {output.bam} > {log.markdup} 2>&1; "
"""