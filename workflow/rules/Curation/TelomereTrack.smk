rule telo_finder: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta",
    output:
        canonical=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.canonical.txt",
        canonical_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.canonical.kmer",
        canonical_top_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.canonical.top.kmer",
        non_canonical=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.non_canonical.txt",
        non_canonical_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.non_canonical.kmer",
        non_canonical_top_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.non_canonical.top.kmer",
    params:
        size=parse_option("size", parameters["tool_options"]["telo_finder"],  "--size", default_value="default"),
        min_kmer=parse_option("min_kmer", parameters["tool_options"]["telo_finder"], "--klo", default_value="default"),
        max_kmer=parse_option("max_kmer", parameters["tool_options"]["telo_finder"], "--khi", default_value="default"),
        ends=parse_option("ends", parameters["tool_options"]["telo_finder"], "--ends", default_value="default")
    log:
        std=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.log",
        cp=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cp.log",
        grep=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.grep.log",
        sed=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sed.log",
        tee=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.tee.log",
        head=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.head.log",
        cp1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cp1.log",
        grep1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.grep1.log",
        sed1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sed1.log",
        tee1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.tee1.log",
        head1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.head1.log",
        cluster_log=output_dict["cluster_log"] / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"]
    threads: parameters["threads"]["telo_finder"]

    shell:
        " STD_LOG=`realpath {log.std}`; "
        " CP_LOG=`realpath {log.cp}`; "
        " GREP_LOG=`realpath {log.grep}`; "
        " SED_LOG=`realpath {log.sed}`; "
        " TEE_LOG=`realpath {log.tee}`; "
        " HEAD_LOG=`realpath {log.head}`; "
        " CP1_LOG=`realpath {log.cp1}`; "
        " GREP1_LOG=`realpath {log.grep1}`; "
        " SED1_LOG=`realpath {log.sed1}`; "
        " TEE1_LOG=`realpath {log.tee1}`; "
        " HEAD1_LOG=`realpath {log.head1}`; "
        " WORKDIR=`dirname {output.canonical}`; "
        " SCRIPT=`realpath workflow/scripts/rapid_curation/telo_finder.py `; "
        " FASTA=`realpath {input.fasta}`;"
        " cd ${{WORKDIR}}; "
        " ${{SCRIPT}} {params.size} {params.max_kmer} {params.max_kmer} "
        " {params.ends} ${{FASTA}} > ${{STD_LOG}} 2>&1; "
        " cp canonical.txt `basename {output.canonical}` > ${{CP_LOG}} 2>&1; "
        " grep 'candidate_telo_kmer:' `basename {output.canonical}` 2>${{GREP_LOG}} | sed 's/.*\t//' 2>${{SED_LOG}} | "
        " tee `basename {output.canonical_kmer}` 2>${{TEE_LOG}} | head -n 1 > `basename {output.canonical_top_kmer}` 2>${{HEAD_LOG}}; "
        " cp non_canonical.txt `basename {output.non_canonical}` > ${{CP1_LOG}} 2>&1; "
        " grep 'candidate_telo_kmer:' `basename {output.non_canonical}` 2>${{GREP1_LOG}} | sed 's/.*\t//' 2>${{SED1_LOG}} | "
        " tee `basename {output.canonical_kmer}` 2>${{TEE1_LOG}} | head -n 1 > `basename {output.canonical_top_kmer}` 2>${{HEAD1_LOG}}; "