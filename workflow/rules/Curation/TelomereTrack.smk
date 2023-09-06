ruleorder: telo_container > create_bedgraph_track
ruleorder: get_telomere_warning > create_bedgraph_track

rule telo_finder:
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta",
    output:
        canonical=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.canonical.txt",
        canonical_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.canonical.kmer",
        canonical_top_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.canonical.top.kmer",
        non_canonical=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.non_canonical.txt",
        non_canonical_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.non_canonical.kmer",
        non_canonical_top_kmer=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.non_canonical.top.kmer",
    params:
        size=parse_option("size", parameters["tool_options"]["telo_finder"],  "--size", default_value="default"),
        min_kmer=parse_option("min_kmer", parameters["tool_options"]["telo_finder"], "--klo", default_value="default"),
        max_kmer=parse_option("max_kmer", parameters["tool_options"]["telo_finder"], "--khi", default_value="default"),
        ends=parse_option("ends", parameters["tool_options"]["telo_finder"], "--ends", default_value="default")
    log:
        std=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.log",
        cp=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cp.log",
        grep=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.grep.log",
        sed=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.sed.log",
        tee=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.tee.log",
        head=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.head.log",
        cp1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cp1.log",
        grep1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.grep1.log",
        sed1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.sed1.log",
        tee1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.tee1.log",
        head1=output_dict["log"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.head1.log",
        cluster_log=output_dict["cluster_log"] / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "telo_finder.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_finder"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"]
    threads: parameters["threads"]["telo_finder"]

    shell:
        " STD_LOG=`realpath -s {log.std}`; "
        " CP_LOG=`realpath -s {log.cp}`; "
        " GREP_LOG=`realpath -s {log.grep}`; "
        " SED_LOG=`realpath -s {log.sed}`; "
        " TEE_LOG=`realpath -s {log.tee}`; "
        " HEAD_LOG=`realpath -s {log.head}`; "
        " CP1_LOG=`realpath -s {log.cp1}`; "
        " GREP1_LOG=`realpath -s {log.grep1}`; "
        " SED1_LOG=`realpath -s {log.sed1}`; "
        " TEE1_LOG=`realpath -s {log.tee1}`; "
        " HEAD1_LOG=`realpath -s {log.head1}`; "
        " WORKDIR=`dirname {output.canonical}`; "
        " SCRIPT=`realpath -s workflow/external_tools/rapid_curation/telo_finder.py `; "
        " FASTA=`realpath -s {input.fasta}`;"
        " cd ${{WORKDIR}}; "
        " ${{SCRIPT}} {params.size} {params.max_kmer} {params.max_kmer} "
        " {params.ends} ${{FASTA}} > ${{STD_LOG}} 2>&1; "
        " cp canonical.txt `basename {output.canonical}` > ${{CP_LOG}} 2>&1; "
        " if [ -s `basename {output.canonical}` ]; "
        " then "
        "       grep 'telo_kmer:' `basename {output.canonical}` 2>${{GREP_LOG}} | sed 's/.*\\t//' 2>${{SED_LOG}} | "
        "       tee `basename {output.canonical_kmer}` 2>${{TEE_LOG}} | head -n 1 > `basename {output.canonical_top_kmer}` 2>${{HEAD_LOG}}; "
        " else "
        "       touch `basename {output.canonical_kmer}`; "
        "       touch `basename {output.canonical_top_kmer}`;  "
        " fi; "
        " cp non_canonical.txt `basename {output.non_canonical}` > ${{CP1_LOG}} 2>&1; "
        " if [ -s `basename {output.non_canonical}` ]; "
        " then "
        "       grep 'candidate_telo_kmer:' `basename {output.non_canonical}` 2>${{GREP1_LOG}} | sed 's/.*\\t//' 2>${{SED1_LOG}} | "
        "       tee `basename {output.non_canonical_kmer}` 2>${{TEE1_LOG}} | head -n 1 > `basename {output.non_canonical_top_kmer}` 2>${{HEAD1_LOG}}; "
        " else "
        "       touch `basename {output.non_canonical_kmer}`; "
        "       touch `basename {output.non_canonical_top_kmer}`;  "
        " fi; "

rule telo_container: #TODO: add possibility to use custom telomere c
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta",
        non_cannonicaL_top_kmer=rules.telo_finder.output.non_canonical_top_kmer,
        cannonicaL_top_kmer=rules.telo_finder.output.canonical_top_kmer
    output:
        cannonical_telo_track=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical_telomere.win1000.step200.track.bedgraph",
        cannonical_telo_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical.telomere.bed",
        cannonical_telo=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical.telomere",
        cannonical_telo_win=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical.telomere.windows",
        non_cannonical_telo_track=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.non_cannonical_telomere.win1000.step200.track.bedgraph",
        non_cannonical_telo_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.non_cannonical.telomere.bed",
        non_cannonical_telo=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.non_cannonical.telomere",
        non_cannonical_telo_win=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.non_cannonical.telomere.windows",
    params:
        container=config["tool_containers"]["rapid_telomere"]
    log:
        std=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.log",
        mkdir=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.mkdir.log",
        cp=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cp.log",
        cp_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cp_cannonical.log",
        cp_non_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cp_non_cannonical.log",
        touch_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.touch_cannonical.log",
        touch_non_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.rm_non_cannonical.log",
        rm_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.rm_cannonical.log",
        rm_non_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.touch_non_cannonical.log",
        rm=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.rm.log",
        cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cannonical.log",
        non_cannonical=output_dict["log"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.non_cannonical.log",
        cluster_log=output_dict["cluster_log"] / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "telo_container.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.benchmark.txt"
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
        " FINALDIR=`dirname {output.cannonical_telo_track}`; "
        " WORKDIR=${{FINALDIR}}/telo_tmp/; "
        " DESTDIR=${{WORKDIR}}/results/; "
        " TEMPDIR=${{WORKDIR}}/tmp/; "
        " HICDIR=${{WORKDIR}}/hic/; "
        " mkdir -p ${{DESTDIR}} ${{TEMPDIR}} ${{HICDIR}} > {log.mkdir} 2>&1; "
        " cp {input.fasta} ${{WORKDIR}}/ref.fa > {log.cp} 2>&1; "
        " export SINGULARITY_BIND=${{WORKDIR}}:/data,${{HICDIR}}:/hic,${{DESTDIR}}:/output,${{TEMPDIR}}:/tmp; "
        " if [ -s {input.cannonicaL_top_kmer} ]; "
        " then "
        "       CANNONICAL_TEL_KMER=`head -n 1 {input.cannonicaL_top_kmer}`; "
        "       singularity run {params.container} -t {wildcards.genome_prefix}_cannonical "
        "       -s ${{CANNONICAL_TEL_KMER}} > {log.cannonical} 2>&1 || true ; "
        "       sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/{wildcards.genome_prefix}_cannonical_telomere.bedgraph > {output.cannonical_telo_track} 2>{log.cp_cannonical}; "
        "       sort -k1,1V -k2,2n -k3,3n  ${{DESTDIR}}/{wildcards.genome_prefix}_cannonical_telomere.bed > {output.cannonical_telo_bed} 2>>{log.cp_cannonical}; "
        "       cp ${{DESTDIR}}/ref.telomere {output.cannonical_telo} >> {log.cp_cannonical} 2>&1; "
        "       cp ${{DESTDIR}}/ref.windows {output.cannonical_telo_win} >> {log.cp_cannonical} 2>&1; "
        "       rm -r ${{DESTDIR}}/* > {log.rm_cannonical} 2>&1; "
        " else "
        "       touch {output.cannonical_telo_track} > {log.touch_cannonical} 2>&1; "
        "       touch {output.cannonical_telo_bed} >> {log.touch_cannonical} 2>&1; "
        "       touch {output.cannonical_telo} >> {log.touch_cannonical} 2>&1; "
        "       touch {output.cannonical_telo_win} >>  {log.touch_cannonical}  2>&1; "
        " fi; "
        " > {log.touch_non_cannonical}; "
        " if [ -s {input.non_cannonicaL_top_kmer} ]; "
        " then "
         "      NON_CANNONICAL_TEL_KMER=`head -n 1 {input.non_cannonicaL_top_kmer}`; "
        "       singularity run {params.container} -t {wildcards.genome_prefix}_non_cannonical "
        "       -s ${{NON_CANNONICAL_TEL_KMER}}> {log.non_cannonical} 2>&1 || true; "
        "       if [ -s '${{DESTDIR}}/ref.telomere' ]; "
        "       then"
        "           sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/{wildcards.genome_prefix}_non_cannonical_telomere.bedgraph > {output.non_cannonical_telo_track} 2>{log.cp_non_cannonical}; "
        "           sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/{wildcards.genome_prefix}_non_cannonical_telomere.bed > {output.non_cannonical_telo_bed} 2>>{log.cp_non_cannonical}; "
        "           cp ${{DESTDIR}}/ref.telomere {output.non_cannonical_telo} >> {log.cp_non_cannonical} 2>&1; "
        "           cp ${{DESTDIR}}/ref.windows {output.non_cannonical_telo_win} >> {log.cp_non_cannonical} 2>&1; "
        "       else"
        "           touch  {output.non_cannonical_telo_track}  >> {log.touch_non_cannonical} 2>&1; "
        "           touch  {output.non_cannonical_telo_bed} >> {log.touch_non_cannonical} 2>&1; "
        "           touch  {output.non_cannonical_telo}  >> {log.touch_non_cannonical} 2>&1; "
        "           touch  {output.non_cannonical_telo_win} >> {log.touch_non_cannonical} 2>&1; "
        "       fi; "
        "       rm -r ${{DESTDIR}}/* > {log.rm_non_cannonical} 2>&1; "
        " else "
        "       touch {output.non_cannonical_telo_track} >> {log.touch_non_cannonical} 2>&1; "
        "       touch {output.non_cannonical_telo_bed} >> {log.touch_non_cannonical} 2>&1; "
        "       touch {output.non_cannonical_telo} >> {log.touch_non_cannonical} 2>&1; "
        "       touch {output.non_cannonical_telo_win} >>  {log.touch_non_cannonical}  2>&1; "
        " fi; "
        " rm -r ${{WORKDIR}} > {log.rm} 2>&1; "

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