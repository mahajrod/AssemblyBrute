

rule combine_haplotypes:
    input:
        fasta_list=lambda wildcards: expand(out_dir_path / ("%s/%s/%s.%s.{haplotype}.fasta" % (wildcards.assembly_stage,
                                                                                                     wildcards.parameters,
                                                                                                     wildcards.genome_prefix,
                                                                                                     wildcards.assembly_stage)),
                         haplotype=stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"] ,
                         allow_missing=True)
    output:
        combined_fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.combined.fasta",
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

rule bam2bed_for_hic_map:
    input:
        bam=out_dir_path / "{assembly_stage}/{parameters}/combined/alignment/NA/{genome_prefix}.{assembly_stage}.NA.combined.rmdup.bam",

    output:
        bed=out_dir_path / "{assembly_stage}/{parameters}/combined/alignment/NA/{genome_prefix}.{assembly_stage}.NA.combined.nodup.bed"
    log:
        samtools=output_dict["log"]  / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.samtools.log",
        bam2bed=output_dict["log"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.bam2bed.log",
        sort=output_dict["log"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam2bed_for_hic_map.{assembly_stage}..{parameters}.{genome_prefix}.benchmark.txt"
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
        " samtools view -@ 4 -u -F0x400 {input.bam} 2>{log.samtools} | bamToBed 2>{log.bam2bed} | "
        " sort -k4 --parallel=10 -S50G -T . > {output.bed} 2>{log.sort}; "


rule bed2pairs_for_hic_map:
    input:
        bed=out_dir_path / "{assembly_stage}/{parameters}/combined/alignment/NA/{genome_prefix}.{assembly_stage}.NA.combined.nodup.bed",
    output:
        pairs=out_dir_path / "{assembly_stage}/{parameters}/combined/alignment/NA/{genome_prefix}.{assembly_stage}.NA.combined.nodup.pairs",
    log:
        paste=output_dict["log"]  / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.paste.log",
        awk1=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.awk1.log",
        awk2=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.awk2.log",
        tr=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.tr.log",
        sort=output_dict["log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bed2pairs.{assembly_stage}..{parameters}.{genome_prefix}.benchmark.txt"
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
        " paste -d '\t' - - < {input.bed} 2>{log.paste} | "
        " awk 'BEGIN {{FS=\"\t\"; OFS=\"\t\"}} "
        "      {{if ($1 > $7) {{print substr($4,1,length($4)-2),$12,$7,$8,\"16\",$6,$1,$2,\"8\",$11,$5}} "
        "      else {{print substr($4,1,length($4)-2),$6,$1,$2,\"8\",$12,$7,$8,\"16\",$5,$11}} }}' 2>{log.awk1} | "
        " tr '\-+' '01' 2>{log.tr} | sort -k3,3d -k7,7d 2>{log.sort} | awk 'NF==11' > {output.pairs} 2>{log.awk2} "

        #" cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 {output.genome_higlass}:1000 {output.higlass_bed} {output.higlass_cool} 2>{log.cload}; "
        #" cooler zoomify --resolutions 5000,10000,20000,40000,60000,80000,100000,120000,150000,200000,300000,400000,500000,1000000,2500000 "
        #" -o {output.higlass_mcool} {output.higlass_cool} 2>{log.zoomify}; "

rule bwa_map_for_hic_map: #
    input:
        index=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.combined.fasta.ann",
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.combined.fasta",
        #fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}".format(wildcards.fileprefix, config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
        #                        out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
        #                                                                                 detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
        #                                                                                 wildcards.haplotype,
        #                                                                                 wildcards.phasing_kmer_length,
        #                                                                                 wildcards.fileprefix,
        #                                                                                 config["fastq_extension"]
        #                                                                                 ),
        forward_fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1",
                                                                                               config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}{6}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                            detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                            wildcards.haplotype,
                                                                                            wildcards.phasing_kmer_length,
                                                                                            wildcards.pairprefix,
                                                                                            input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1",
                                                                                            config["fastq_extension"]),
        reverse_fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2",
                                                                                               config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}{6}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                            detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                            wildcards.haplotype,
                                                                                            wildcards.phasing_kmer_length,
                                                                                            wildcards.pairprefix,
                                                                                            input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2",
                                                                                            config["fastq_extension"]),
    output:
        #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
        bam=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/combined/alignment/{phasing_kmer_length, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.combined.{pairprefix, [^/]+}.bwa.bam"
    params:
        id="{0}_hic".format(config["genome_prefix"]),
        bwa_tool=config["bwa_tool"]
    log:
        map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.combined.{phasing_kmer_length}.{pairprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.combined.{phasing_kmer_length}.{pairprefix}.sort.log",
        #filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.combined.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
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
        " {params.bwa_tool} mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "

rule bam_merge_files_for_hic_map:
    input:
        bams=expand(rules.bwa_map.output.bam, #out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.filtered.{pairprefix}.bam",
                    allow_missing=True,
                    pairprefix=input_pairprefix_dict["hic"]), #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        reference_fai=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.combined.fasta.fai",
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.combined.fasta"
    output:
        bam=out_dir_path / "{assembly_stage}/{parameters}/combined/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.combined.bwa.bam" # TODO: make temp
    params:
        sort_threads=parameters["threads"]["samtools_sort"]
    log:
        std=output_dict["log"] / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("bwa_merge_files"),
        cpus=parameters["threads"]["samtools_sort"] ,
        time=parameters["time"]["samtools_sort"],
        mem=parameters["memory_mb"]["samtools_sort"]
    threads: parameters["threads"]["samtools_sort"]
    shell:
        " samtools merge -@ {params.sort_threads} -o {output.bam} {input.bams} 1>{log.std} 2>&1"

rule rmdup_for_hic_map:
    input:
        bam=rules.bam_merge_files.output.bam
    output:
        bam=out_dir_path / "{assembly_stage}/{parameters}/combined/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.combined.rmdup.bam",
        #bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        collate_threads=parameters["threads"]["samtools_collate"],
        fixmate_threads=parameters["threads"]["samtools_fixmate"],
        markdup_threads=parameters["threads"]["samtools_markdup"],
        sort_per_thread=parameters["memory_mb"]["samtools_sort_per_thread"]
    log:
        collate=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.collate.log",
        fixmate=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.fixmate.log",
        sort=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.sort.log",
        markdup=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.markdup.log",
        cluster_log=output_dict["cluster_log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.cluster.log",
        cluster_err=output_dict["cluster_error"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.combined.benchmark.txt"
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