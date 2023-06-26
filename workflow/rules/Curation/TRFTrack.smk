
rule trf: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta"
    output:
        bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.trf.simple.bed",
        merged_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.trf.track.bed",
    params:
        #lambda wildcards: parse_option("window_size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -w ")
        matching_weight=parse_option("matching_weight", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -m "),
        mismatching_penalty=parse_option("mismatching_penalty", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -s "),
        indel_penalty=parse_option("indel_penalty", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -l "),
        match_probability=parse_option("match_probability", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -a "),
        indel_probability=parse_option("indel_probability", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -d "),
        min_alignment_score=parse_option("min_alignment_score", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -c "),
        max_period=parse_option("max_period", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -e "),
        max_repeat_length=parse_option("max_repeat_length", parameters["tool_options"]["trf"][config["other_tool_option_sets"]["trf"]], " -g "),
    log:
        trf=output_dict["log"]  / "trf.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.trf.log",
        sort=output_dict["log"]  / "trf.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort.log",
        merge=output_dict["log"]  / "trf.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.merge.log",
        cluster_log=output_dict["cluster_log"] / "trf.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "trf.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "trf.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["trf"] ,
        time=parameters["time"]["trf"],
        mem=parameters["memory_mb"]["trf"]
    threads: parameters["threads"]["trf"]
    shell:
        " WORK_DIR=`dirname {output.bed}`; "
        " INPUT_FASTA=`basename  {input.fasta}`; "
        " OUTPUT_PREFIX=`basename {output.bed}`; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.simple.bed}};"
        " cd ${{WORK_DIR}}; "
        " tandem_repeat_masking.py -t {threads} {params.matching_weight} {params.mismatching_penalty} "
        " {params.indel_penalty} {params.match_probability} {params.indel_probability} {params.min_alignment_score} "
        " {params.max_period} {params.max_repeat_length} -i ${{INPUT_FASTA}} -o ${{OUTPUT_PREFIX}} > {log.trf} 2>&1; "
        " sort -k1,1V -k2,2n -k3,3n {output.bed} 2>{log.sort} | bedtools merge -i stdin > {output.merged_bed} 2>{log.merge}; "