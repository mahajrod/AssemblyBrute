"""

rule bwa_map: #
    input:
        index=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.ann",
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}".format(wildcards.fileprefix,
                                                                                    config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                         detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                         wildcards.haplotype,
                                                                                         wildcards.phasing_kmer_length,
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]
                                                                                         )
    output:
        bam=temp(out_dir_path  / "{assembly_stage}/{parameters}/{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam")
    params:
        id="{0}_hic".format(config["genome_prefix"]),
        bwa_tool=config["bwa_tool"],
    log:
        fastx=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.fastx.log",
        map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.sort.log",
        filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=parameters["threads"]["bwa_map_arima"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"],
    threads: parameters["threads"]["bwa_map_arima"]
    shell:
        " {params.bwa_tool} mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} <(zcat {input.fastq}  2>{log.fastx}) 2>{log.map} |"
        " workflow/external_tools/arima_mapping_pipeline/filter_five_end.pl 2>{log.filter} | samtools view -Sb - > {output.bam} 2>{log.sort} "
"""
rule arima_bwa_map: #
    input:
        index=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.ann",
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}".format(wildcards.fileprefix,
                                                                                    config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                         detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                         wildcards.haplotype,
                                                                                         wildcards.phasing_kmer_length,
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]
                                                                                         )
    output:
        raw_bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.raw.bam"
    params:
        id="{0}_hic".format(config["genome_prefix"]),
        bwa_tool=config["bwa_tool"],
    log:
        fastx=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.fastx.log",
        map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.sort.log",
        filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=get_threads(parameters["threads"]["bwa_map_arima"], "cpu"),
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"],
    threads: parameters["threads"]["bwa_map_arima"]
    shell:
        " {params.bwa_tool} mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} <(zcat {input.fastq}  2>{log.fastx}) 2>{log.map} | samtools view -Sb - > {output.raw_bam} 2>{log.sort} "

rule arima_filter_five_end: #
    input:
        raw_bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.raw.bam",
        stats=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.raw.bam.general_stats"
    output:
        bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
    params:
        id="{0}_hic".format(config["genome_prefix"]),
        bwa_tool=config["bwa_tool"],
    log:
        view=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.view.log",
        filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.filter.log",
        view2=output_dict["log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.view2.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=parameters["threads"]["arima_filter_five_end"] ,
        time=parameters["time"]["arima_filter_five_end"],
        mem=parameters["memory_mb"]["arima_filter_five_end"],
    threads: parameters["threads"]["arima_filter_five_end"]
    shell:
        " samtools view -h {input.raw_bam} 2>{log.view} | "
        " workflow/external_tools/arima_mapping_pipeline/filter_five_end.pl 2>{log.filter} | "
        " samtools view -Sb - > {output.bam} 2>{log.view2} "

rule arima_two_read_bam_combiner:
    input:
        forward_bam=lambda wildcards: out_dir_path / ("{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.{5}{6}.bwa.bam".format(wildcards.assembly_stage,
                                                                                                                                  wildcards.parameters,
                                                                                                                                  wildcards.haplotype,
                                                                                                                                  wildcards.phasing_kmer_length,
                                                                                                                                  wildcards.genome_prefix,
                                                                                                                                  wildcards.pairprefix,
                                                                                                                                  input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1")),
        forward_stats=lambda wildcards: out_dir_path / ("{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.{5}{6}.bwa.bam.general_stats".format(wildcards.assembly_stage,
                                                                                                                                        wildcards.parameters,
                                                                                                                                        wildcards.haplotype,
                                                                                                                                        wildcards.phasing_kmer_length,
                                                                                                                                        wildcards.genome_prefix,
                                                                                                                                        wildcards.pairprefix,
                                                                                                                                        input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1")),
        reverse_bam=lambda wildcards: out_dir_path / ("{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.{5}{6}.bwa.bam".format(wildcards.assembly_stage,
                                                                                                                                  wildcards.parameters,
                                                                                                                                  wildcards.haplotype,
                                                                                                                                  wildcards.phasing_kmer_length,
                                                                                                                                  wildcards.genome_prefix,
                                                                                                                                  wildcards.pairprefix,
                                                                                                                                  input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2")),
        reverse_stats=lambda wildcards: out_dir_path / ("{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.{5}{6}.bwa.bam.general_stats".format(wildcards.assembly_stage,
                                                                                                                                        wildcards.parameters,
                                                                                                                                        wildcards.haplotype,
                                                                                                                                        wildcards.phasing_kmer_length,
                                                                                                                                        wildcards.genome_prefix,
                                                                                                                                        wildcards.pairprefix,
                                                                                                                                        input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2")),
        reference_fai=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai"
    output:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix}.bwa.bam", # TODO: make_tem
    params:
        min_mapq=parameters["tool_options"]["two_read_bam_combiner"]["mapq"],
        sort_threads=parameters["threads"]["samtools_sort"],
        sort_memory=parameters["memory_mb"]["samtools_sort"],
    log:
        merge=output_dict["log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.merge.log",
        view=output_dict["log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.view.log",
        sort=output_dict["log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{phasing_kmer_length}.{genome_prefix}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{phasing_kmer_length}.{genome_prefix}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_pairs.{assembly_stage}.{parameters}.{phasing_kmer_length}.{genome_prefix}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["arima_two_read_bam_combiner"] ,
        node_options=parse_node_list("bwa_merge_pairs"),
        time=parameters["time"]["arima_two_read_bam_combiner"],
        mem=parameters["memory_mb"]["arima_two_read_bam_combiner"] + parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"] + 30000,
        queue=config["queue"]["cpu"]["name"]
    threads: parameters["threads"]["arima_two_read_bam_combiner"] + parameters["threads"]["samtools_sort"]
    shell:
        " TMP_PREFIX=`dirname {output.bam}`/{wildcards.pairprefix}; "
        " workflow/external_tools/arima_mapping_pipeline/two_read_bam_combiner.pl {input.forward_bam} "
        " {input.reverse_bam} samtools {params.min_mapq} 2>{log.merge} | "
        " samtools view -bS -t {input.reference_fai} - 2>{log.view} | "
        " samtools sort -T ${{TMP_PREFIX}} -m {params.sort_memory}M -@ {params.sort_threads} -o {output.bam} 2>{log.sort}"
