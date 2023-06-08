rule create_higlass_track: #
    input:
        fai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta.fai",
        repeat_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.repeat.bed",
    output:
        genome_higlass=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.higlass.genome",
        higlass_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.higlass.bed",
        higlass_cool=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.higlass.cool",
        higlass_mcool=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.higlass.mcool",
    log:
        paste=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.paste.log",
        cut=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cut.log",
        sed1=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sed1.log",
        sort=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort.log",
        sed2=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sed2.log",
        awk=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.awk.log",
        tr=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.tr.log",
        sort2=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort.log",
        cload=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cload.log",
        zoomify=output_dict["log"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.zoomify.log",
        cluster_log=output_dict["cluster_log"] / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_higlass_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_higlass_track"] ,
        time=parameters["time"]["create_higlass_track"],
        mem=parameters["memory_mb"]["create_higlass_track"]
    threads: parameters["threads"]["create_higlass_track"]
    shell:
        " cut -f1,2 {input.fai} 2>{log.cut} | sed 's/-/_/g' 2>{log.sed1} | sort -k2,2 -nr > {output.genome_higlass} 2>{log.sort}; "
        " paste -d '\\t' - - < $bed 2>{log.paste} | sed 's/-/_/g' 2>{log.sed2} | "
        " awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} {{if ($1 > $7) {{print substr($4,1,length($4)-2),$12,$7,$8,\"16\",$6,$1,$2,\"8\",$11,$5}} "
        "                                          else {{ print substr($4,1,length($4)-2),$6,$1,$2,\"8\",$12,$7,$8,\"16\",$5,$11}} }}' 2>{log.awk} | "
        " tr '\\-+' '01' 2>{log.tr} | sort --parallel={threads} -S{resources.mem}M -k3,3d -k7,7d > {output.higlass_bed} 2>{log.sort2}; "
        " cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 {output.genome_higlass}:1000 {output.higlass_bed} {output.higlass_cool} 2>{log.cload}; "
        " cooler zoomify -o {output.higlass_mcool} {output.higlass_cool} 2>{log.zoomify}; "
