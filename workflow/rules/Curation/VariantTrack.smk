if "purge_dups" in config["stage_list"]:
    ruleorder: minimap2_cov > minimap2_purge_dups_reads

rule add_basequalities_to_bam: #adds basequalities to bam generated from fasta (for deepvariant)
    input:
        bam=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.bam",
        bai=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.bam.bai",
        #reference=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta"
    output:
        bam=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.with_qual.bam",
        bai=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.with_qual.bam.bai",
    params:
        dataformat=lambda wildcards: datatype_format_dict[wildcards.datatype]
    log:
        view1=output_dict["log"]  / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.view1.log",
        add=output_dict["log"]  / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.add.log",
        view2=output_dict["log"]  / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.view2.log",
        index=output_dict["log"]  / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.index.log",
        ln=output_dict["log"]  / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "add_basequalities_to_bam.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("add_basequalities_to_bam"),
        cpus=parameters["threads"]["add_basequalities_to_bam"],
        time=parameters["time"]["add_basequalities_to_bam"],
        mem=parameters["memory_mb"]["add_basequalities_to_bam"]
    threads: parameters["threads"]["add_basequalities_to_bam"]

    shell:
        " if [ '{params.dataformat}' == 'fasta' ] ;"
        "   then "
        "       samtools view -@ 4 -h {input.bam} 2>{log.view1} | "
        "       workflow/scripts/curation/add_basequalities_to_bam.py 2>{log.add} |  "
        "       samtools view -b -@ 4 > {output.bam} 2>{log.view2}; "
        "       samtools index {output.bam} 2>{log.index}; "
        "   else "
        "       LN_LOG=`realpath {log.ln}`; "
        "       cd `dirname {input.bam}`; "
        "       ln -sf `basename {input.bam}` `basename {output.bam}` 2>${{LN_LOG}}; "
        "       ln -sf `basename {input.bai}` `basename {output.bai}` 2>>${{LN_LOG}}; "
        "   fi"

rule deepvariant: #
    input:
        bam=rules.add_basequalities_to_bam.output.bam,
        bai=rules.add_basequalities_to_bam.output.bai,
        reference=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        vcf=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.vcf.gz",
        gvcf=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.g.vcf.gz",
        #paf=out_dir_path  / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        sif=config["tool_containers"]["deepvariant"]["gpu"] if config["queue"]["gpu"] and config["queue"]["gpu"] and config["tool_containers"]["deepvariant"]["gpu"] else config["tool_containers"]["deepvariant"]["cpu"],
        model=lambda wildcards: parameters["tool_options"]["deepvariant"][wildcards.datatype]["model"],
        gpu_options=" --nv " if config["queue"]["gpu"] and config["queue"]["gpu"] and config["tool_containers"]["deepvariant"]["gpu"] else " " # enables nVidia support
        #tmp_dir=config["alternative_tmp_dir"],
    log:
        deepvariant=output_dict["log"]  / "deepvariant.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.deepvariant.log",
        cluster_log=output_dict["cluster_log"] / "deepvariant.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "deepvariant.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "deepvariant.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{seq_type}.{genome_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        queue=config["queue"]["gpu"] if config["queue"]["gpu"] and config["queue"]["gpu"] and config["tool_containers"]["deepvariant"]["gpu"] else config["queue"]["cpu"],
        node_options=parse_node_list("deepvariant"),
        cpus=parameters["threads"]["deepvariant"],
        time=parameters["time"]["deepvariant"],
        mem=parameters["memory_mb"]["deepvariant"]
    threads: parameters["threads"]["deepvariant"]

    shell:
        " WORKDIR=`dirname {output.vcf}`; "
        " SIF=`realpath {params.sif}`; "
        " LOG=`realpath {log.deepvariant}`; "
        " cd ${{WORKDIR}}; "
        " singularity run {params.gpu_options} -B /usr/lib/locale/:/usr/lib/locale/ ${{SIF}} /opt/deepvariant/bin/run_deepvariant "
        " --model_type={params.model} --ref=`basename {input.reference}` --reads=`basename {input.bam}` "
        " --output_vcf=`basename {output.vcf}` --output_gvcf=`basename {output.gvcf}` "
        " --intermediate_results_dir ./deepvariant --num_shards={threads} > ${{LOG}} 2>&1; "

"""
        " mkdir -p {params.tmp_dir}; "
        " TMPDIR=`realpath {params.tmp_dir}`;"
        " WORKDIR=`dirname {output.vcf}`; "
        " SIF=`realpath {params.sif}`; "
        " LOG=`realpath {log.deepvariant}`; "
        " cd ${{WORKDIR}}; "
        " export SINGULARITYENV_TMPDIR=/mnt/tmp; "
        " export SINGULARITYENV_SQLITE_TMPDIR=/mnt/tmp; "
        " TMPDIR=${{TMPDIR}} SINGULARITYENV_TMPDIR=/mnt/tmp SINGULARITYENV_SQLITE_TMPDIR=/mnt/tmp "
        " singularity run -B /usr/lib/locale/:/usr/lib/locale/ -B `realpath {params.tmp_dir}`:/mnt/tmp/  ${{SIF}} /opt/deepvariant/bin/run_deepvariant "
        " --model_type={params.model} --ref=`basename {input.reference}` --reads=`basename {input.bam}` "
        " --output_vcf=`basename {output.vcf}` --output_gvcf=`basename {output.gvcf}` "
        " --intermediate_results_dir ./deepvariant --num_shards={threads} > ${{LOG}} 2>&1; "
        #" --call_variants_extra_args='config_string=\"device_count {{key: \'cpu\' value: {threads}}} intra_op_parallelism_threads:{threads} inter_op_parallelism_threads:{threads}\"' "
"""
"""
rule deepvariant_filter: #
    input:
        vcf=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.vcf.gz",
    output:
        vcf=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.filtered.vcf.gz",
    log:
        filter=output_dict["log"]  / "deepvariant_filter.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.filter.log",
        cluster_log=output_dict["cluster_log"] / "deepvariant_filter.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "deepvariant_filter.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "deepvariant_filter.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["deepvariant_filter"],
        time=parameters["time"]["deepvariant_filter"],
        mem=parameters["memory_mb"]["deepvariant_filter"]
    threads: parameters["threads"]["deepvariant_filter"]

    shell:
        " bcftools filter --exclude 'FILTER~\"\"' > {log.filter} 2>&1; "

"""