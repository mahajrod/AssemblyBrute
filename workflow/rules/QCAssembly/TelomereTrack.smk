ruleorder: telo_container > create_bedgraph_track
ruleorder: get_telomere_warning > create_bedgraph_track

rule telo_finder:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta"
    output:
        canonical="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.canonical.txt",
        canonical_kmer="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.canonical.kmer",
        canonical_top_kmer="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.canonical.top.kmer",
        non_canonical="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_canonical.txt",
        non_canonical_kmer="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_canonical.kmer",
        non_canonical_top_kmer="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_canonical.top.kmer",
    params:
        size=parse_option("size", parameters["tool_options"]["telo_finder"],  "--size", default_value="default"),
        min_kmer=parse_option("min_kmer", parameters["tool_options"]["telo_finder"], "--klo", default_value="default"),
        max_kmer=parse_option("max_kmer", parameters["tool_options"]["telo_finder"], "--khi", default_value="default"),
        ends=parse_option("ends", parameters["tool_options"]["telo_finder"], "--ends", default_value="default")
    log:
        std="{fasta_dir}/telomere/{fasta_prefix}.log",
        cluster_log="{fasta_dir}/telomere/{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/telomere/{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/telomere/{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_finder"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"],
        telosif=1
    threads: parameters["threads"]["telo_finder"]

    shell:
        " STD_LOG=`realpath -s {log.std}`; "
        " > ${{STD_LOG}}; "
        " WORKDIR=`dirname {output.canonical}`; "
        " SCRIPT=`realpath -s workflow/external_tools/rapid_curation/telo_finder.py `; "
        " FASTA=`realpath -s {input.fasta}`;"
        " cd ${{WORKDIR}}; "
        " ${{SCRIPT}} {params.size} {params.min_kmer} {params.max_kmer} {params.ends} ${{FASTA}} >> ${{STD_LOG}} 2>&1; "
        " cp canonical.txt `basename {output.canonical}` >> ${{STD_LOG}} 2>&1; "
        " if [ -s `basename {output.canonical}` ]; "
        " then "
        "       grep 'telo_kmer:' `basename {output.canonical}` 2>>${{STD_LOG}} | sed 's/.*\\t//' 2>>${{STD_LOG}} | "
        "       tee `basename {output.canonical_kmer}` 2>>$${{STD_LOG}} | head -n 1 > `basename {output.canonical_top_kmer}` 2>>${{STD_LOG}}; "
        " else "
        "       touch `basename {output.canonical_kmer}`; "
        "       touch `basename {output.canonical_top_kmer}`;  "
        " fi; "
        " cp non_canonical.txt `basename {output.non_canonical}` >> ${{STD_LOG}} 2>&1; "
        " if [ -s `basename {output.non_canonical}` ]; "
        " then "
        "       grep 'candidate_telo_kmer:' `basename {output.non_canonical}` 2>>${{STD_LOG}} | sed 's/.*\\t//' 2>>${{STD_LOG}} | "
        "       tee `basename {output.non_canonical_kmer}` 2>>${{STD_LOG}}| head -n 1 > `basename {output.non_canonical_top_kmer}` 2>>${{STD_LOG}}; "
        " else "
        "       touch `basename {output.non_canonical_kmer}`; "
        "       touch `basename {output.non_canonical_top_kmer}`;  "
        " fi; "

rule telo_container: #TODO: add possibility to use custom telomere c
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        non_cannonical_top_kmer=rules.telo_finder.output.non_canonical_top_kmer,
        cannonical_top_kmer=rules.telo_finder.output.canonical_top_kmer
    output:
        cannonical_telo_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.cannonical_telomere.win1000.step200.track.bedgraph",
        cannonical_telo_bed="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.cannonical.telomere.bed",
        cannonical_telo="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.cannonical.telomere",
        cannonical_telo_win="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.cannonical.telomere.windows",
        non_cannonical_telo_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_cannonical_telomere.win1000.step200.track.bedgraph",
        non_cannonical_telo_bed="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_cannonical.telomere.bed",
        non_cannonical_telo="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_cannonical.telomere",
        non_cannonical_telo_win="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.non_cannonical.telomere.windows",
    params:
        container=config["tool_containers"]["rapid_telomere"]
    log:
        std="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.log",
        cluster_log="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.cluster.log",
        cluster_err="{fasta_dir}/telomere/{fasta_prefix, [^/]+}.cluster.err"
    benchmark:
        "{fasta_dir}/telomere/{fasta_prefix, [^/]+}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_container"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"],
        telosif=1
    threads: parameters["threads"]["telo_finder"]

    shell: # TODO pack code below as script. Issue with container - it return 111 inside of 0 and snakemake breaks. Added ' || true'
        " LOG=`realpath -s {log.std}`; "
        " > ${{LOG}}; "
        " FASTA=`realpath -s {input.fasta}`; "
        " INPUT_CANNONICAL_KMER=`realpath -s {input.cannonical_top_kmer}`; "
        " INPUT_NONCANNONICAL_KMER=`realpath -s {input.non_cannonical_top_kmer}`; "
        " FINALDIR=`dirname {output.cannonical_telo_track}`; "
        " cd ${{FINAL_DIR}}; "
        " WORKDIR=./telo_tmp/; "
        " DESTDIR=${{WORKDIR}}/results/; "
        " TEMPDIR=${{WORKDIR}}/tmp/; "
        " HICDIR=${{WORKDIR}}/hic/; "
        " mkdir -p ${{DESTDIR}} ${{TEMPDIR}} ${{HICDIR}} > ${{LOG}} 2>&1; "
        " cp ${{FASTA}} ${{WORKDIR}}/ref.fa >> ${{LOG}} 2>&1; "
        " export SINGULARITY_BIND=${{WORKDIR}}:/data,${{HICDIR}}:/hic,${{DESTDIR}}:/output,${{TEMPDIR}}:/tmp; "
        " if [ -s ${{INPUT_CANNONICAL_KMER}} ]; "
        " then "
        "       CANNONICAL_TEL_KMER=`head -n 1 ${{INPUT_CANNONICAL_KMER}}`; "
        "       CANNONICAL_OUTPUT_PREFIX=`basename {output.cannonical_telo_bed}`'.tmp'; "
        "       singularity run {params.container} -t ${{CANNONICAL_OUTPUT_PREFIX}} "
        "       -s ${{CANNONICAL_TEL_KMER}} >> ${{LOG}} 2>&1 || true ; "
        "       sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{CANNONICAL_OUTPUT_PREFIX}}_telomere.bedgraph > {output.cannonical_telo_track} 2>>${{LOG}}; "
        "       sort -k1,1V -k2,2n -k3,3n  ${{DESTDIR}}/${{CANNONICAL_OUTPUT_PREFIX}}_telomere.bed > {output.cannonical_telo_bed} 2>>${{LOG}}; "
        "       cp ${{DESTDIR}}/ref.telomere {output.cannonical_telo} >> ${{LOG}} 2>&1; "
        "       cp ${{DESTDIR}}/ref.windows {output.cannonical_telo_win} >> ${{LOG}} 2>&1; "
        "       rm -r ${{DESTDIR}}/* >> ${{LOG}} 2>&1; "
        " else "
        "       touch {output.cannonical_telo_track} >> ${{LOG}} 2>&1; "
        "       touch {output.cannonical_telo_bed} >> ${{LOG}} 2>&1; "
        "       touch {output.cannonical_telo} >> ${{LOG}} 2>&1; "
        "       touch {output.cannonical_telo_win} >> ${{LOG}}  2>&1; "
        " fi; "
        " if [ -s ${{INPUT_NONCANNONICAL_KMER}} ]; "
        " then "
        "       NON_CANNONICAL_TEL_KMER=`head -n 1 ${{INPUT_NONCANNONICAL_KMER}}`; "
        "       NON_CANNONICAL_OUTPUT_PREFIX=`basename {output.non_cannonical_telo_bed}`'.tmp'; "   
        "       singularity run {params.container} -t ${{NON_CANNONICAL_OUTPUT_PREFIX}} "
        "       -s ${{NON_CANNONICAL_TEL_KMER}} >> ${{LOG}} 2>&1 || true; "
        "       if [ -s '${{DESTDIR}}/ref.telomere' ]; "
        "       then"
        "           sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{NON_CANNONICAL_OUTPUT_PREFIX}}_telomere.bedgraph > {output.non_cannonical_telo_track} 2>>${{LOG}}; "
        "           sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{NON_CANNONICAL_OUTPUT_PREFIX}}_telomere.bed > {output.non_cannonical_telo_bed} 2>>${{LOG}}; "
        "           cp ${{DESTDIR}}/ref.telomere {output.non_cannonical_telo} >> ${{LOG}}} 2>&1; "
        "           cp ${{DESTDIR}}/ref.windows {output.non_cannonical_telo_win} >> ${{LOG}} 2>&1; "
        "       else"
        "           touch  {output.non_cannonical_telo_track}  >> ${{LOG}} 2>&1; "
        "           touch  {output.non_cannonical_telo_bed} >> ${{LOG}} 2>&1; "
        "           touch  {output.non_cannonical_telo}  >> ${{LOG}} 2>&1; "
        "           touch  {output.non_cannonical_telo_win} >> ${{LOG}} 2>&1; "
        "       fi; "
        "       rm -r ${{DESTDIR}}/* >> ${{LOG}} 2>&1; "
        " else "
        "       touch {output.non_cannonical_telo_track} >> ${{LOG}} 2>&1; "
        "       touch {output.non_cannonical_telo_bed} >> ${{LOG}} 2>&1; "
        "       touch {output.non_cannonical_telo} >> ${{LOG}} 2>&1; "
        "       touch {output.non_cannonical_telo_win} >> ${{LOG}}  2>&1; "
        " fi; "
        " rm -r ${{WORKDIR}} >> ${{LOG}} 2>&1; "

rule get_telomere_warning:
    input:
        cannonical_telo_track=rules.telo_container.output.cannonical_telo_track, #out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.cannonical.telomere.bedgraph",
        non_cannonical_telo_track=rules.telo_container.output.non_cannonical_telo_track, # out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.non_cannonical.telomere.bedgraph",
        fai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta.fai",
    output:
        cannonical_telo_warning_track=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical_telomere_warning.win1000.step200.track.bedgraph",
        non_cannonical_telo_warning_track=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.non_cannonical_telomere_warning.win1000.step200.track.bedgraph",
    log:
        cannonical=output_dict["log"]  / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cannonical.log",
        non_cannonical=output_dict["log"]  / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.non_cannonical.log",
        cannonical_touch=output_dict["log"]  / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cannonical.touch.log",
        non_cannonical_touch=output_dict["log"]  / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.non_cannonical.touch.log",
        cluster_log=output_dict["cluster_log"] / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_telomere_warning.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_telomere_warning"),
        cpus=parameters["threads"]["get_telomere_warning"] ,
        time=parameters["time"]["get_telomere_warning"],
        mem=parameters["memory_mb"]["get_telomere_warning"]
    threads: parameters["threads"]["get_telomere_warning"]

    shell:
        " if [ -s {input.cannonical_telo_track} ]; "
        " then "
        "       workflow/scripts/curation/find_internal_telomere.py  -i {input.cannonical_telo_track}  "
        "                   -f {input.fai} > {output.cannonical_telo_warning_track} 2>{log.cannonical}; "
        " else"
        "       touch {output.cannonical_telo_warning_track} > {log.cannonical_touch} 2>&1; "
        " fi;"
        " if [ -s {input.non_cannonical_telo_track} ]; "
        " then "
        "       workflow/scripts/curation/find_internal_telomere.py  -i {input.non_cannonical_telo_track}  "
        "                   -f {input.fai} > {output.non_cannonical_telo_warning_track} 2>{log.non_cannonical}; "
        " else"
        "       touch {output.non_cannonical_telo_warning_track} > {log.non_cannonical_touch} 2>&1; "
        " fi; "