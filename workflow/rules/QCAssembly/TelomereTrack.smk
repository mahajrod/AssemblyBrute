ruleorder: telo_container > create_bedgraph_track
ruleorder: get_telomere_warning > create_bedgraph_track

localrules: copy_telomere_files

rule telo_finder:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta"
    output:
        canonical="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.canonical.txt",
        canonical_kmer="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.canonical.kmer",
        canonical_top_kmer="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.canonical.top.kmer",
        non_canonical="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical.txt",
        non_canonical_kmer="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical.kmer",
        non_canonical_top_kmer="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical.top.kmer",
    params:
        size=parse_option("size", parameters["tool_options"]["telo_finder"],  "--size", default_value="default"),
        min_kmer=parse_option("min_kmer", parameters["tool_options"]["telo_finder"], "--klo", default_value="default"),
        max_kmer=parse_option("max_kmer", parameters["tool_options"]["telo_finder"], "--khi", default_value="default"),
        ends=parse_option("ends", parameters["tool_options"]["telo_finder"], "--ends", default_value="default")
    log:
        std="{fasta_dir}/telo_finder.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/telo_finder.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/telo_finder.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/telo_finder.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_finder"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"],
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
        non_canonical_top_kmer=rules.telo_finder.output.non_canonical_top_kmer,
        canonical_top_kmer=rules.telo_finder.output.canonical_top_kmer
    output:
        canonical_telo_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere.win1000.step200.track.bedgraph",
        canonical_telo_bed="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical.telomere.bed",
        canonical_telo="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical.telomere",
        canonical_telo_win="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical.telomere.windows",
        non_canonical_telo_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.bedgraph",
        non_canonical_telo_bed="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical.telomere.bed",
        non_canonical_telo="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical.telomere",
        non_canonical_telo_win="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical.telomere.windows",
    params:
        container=config["tool_containers"]["rapid_telomere"]
    log:
        std="{fasta_dir}/telo_container.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/telo_container.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/telo_container.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/telo_container.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_container"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"],
    threads: parameters["threads"]["telo_finder"]

    shell: # TODO pack code below as script. Issue with container - it return 111 inside of 0 and snakemake breaks. Added ' || true'
        " STARTDIR=`pwd`; "
        " LOG=`realpath -s {log.std}`; "
        " echo \"Started logging...\" > ${{LOG}}; "
        " FASTA=`realpath -s {input.fasta}`; "
        " INPUT_CANONICAL_KMER=`realpath -s {input.canonical_top_kmer}`; "
        " INPUT_NONCANONICAL_KMER=`realpath -s {input.non_canonical_top_kmer}`; "
        " FINALDIR=`dirname {output.canonical_telo_track}`; "
        " FINALDIR=`realpath -s ${{FINALDIR}}`; "
        " cd ${{FINALDIR}} >> ${{LOG}} 2>&1; "
        " WORKDIR=${{FINALDIR}}/telo_tmp/; "
        " DESTDIR=${{WORKDIR}}/results/; "
        " TEMPDIR=${{WORKDIR}}/tmp/; "
        " HICDIR=${{WORKDIR}}/hic/; "
        " mkdir -p ${{DESTDIR}} ${{TEMPDIR}} ${{HICDIR}} >> ${{LOG}} 2>&1; "
        " echo \"Current directory:\" `pwd` >> ${{LOG}}; "
        " ls ./ >> ${{LOG}}; "
        " cp ${{FASTA}} ${{WORKDIR}}/ref.fa >> ${{LOG}} 2>&1; "
        " export SINGULARITY_BIND=${{WORKDIR}}:/data,${{HICDIR}}:/hic,${{DESTDIR}}:/output,${{TEMPDIR}}:/tmp; "
        " if [ -s ${{INPUT_CANONICAL_KMER}} ]; "
        " then "
        "       echo 'Canonical kmer found. Statring analysis...' >> ${{LOG}}; "
        "       CANONICAL_TEL_KMER=`head -n 1 ${{INPUT_CANONICAL_KMER}}`; "
        "       CANONICAL_OUTPUT_PREFIX=`basename {output.canonical_telo_bed}`'.tmp'; "
        "       echo 'Running telo container for canonical kmer...' >> ${{LOG}}; "
        "       singularity run {params.container} -t ${{CANONICAL_OUTPUT_PREFIX}} "
        "       -s ${{CANONICAL_TEL_KMER}} >> ${{LOG}} 2>&1 || true ; "
        "       cd ${{STARTDIR}} >> ${{LOG}} 2>&1; "
        "       echo 'Sorting results...' >> ${{LOG}}; "
        "       sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{CANONICAL_OUTPUT_PREFIX}}_telomere.bedgraph > {output.canonical_telo_track} 2>>${{LOG}}; "
        "       sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{CANONICAL_OUTPUT_PREFIX}}_telomere.bed > {output.canonical_telo_bed} 2>>${{LOG}}; "
        "       echo 'Copying files to the final destination...'  >> ${{LOG}}; "
        "       cp ${{DESTDIR}}/ref.telomere {output.canonical_telo} >> ${{LOG}} 2>&1; "
        "       cp ${{DESTDIR}}/ref.windows {output.canonical_telo_win} >> ${{LOG}} 2>&1; "
        "       rm -r ${{DESTDIR}}/* >> ${{LOG}} 2>&1; "
        " else "
        "       cd ${{STARTDIR}} >> ${{LOG}} 2>&1; "
        "       echo 'No canonical kmer found. Creating empty output files...'  >> ${{LOG}};  "
        "       touch {output.canonical_telo_track} >> ${{LOG}} 2>&1; "
        "       touch {output.canonical_telo_bed} >> ${{LOG}} 2>&1; "
        "       touch {output.canonical_telo} >> ${{LOG}} 2>&1; "
        "       touch {output.canonical_telo_win} >> ${{LOG}}  2>&1; "
        " fi; "
        " cd ${{FINALDIR}} >> ${{LOG}} 2>&1; "
        " if [ -s ${{INPUT_NONCANONICAL_KMER}} ]; "
        " then "
        "       echo 'Non canonical kmer found. Statring analysis...'  >> ${{LOG}}; "
        "       NON_CANONICAL_TEL_KMER=`head -n 1 ${{INPUT_NONCANONICAL_KMER}}`; "
        "       NON_CANONICAL_OUTPUT_PREFIX=`basename {output.non_canonical_telo_bed}`'.tmp'; "
        "       echo 'Running telo container for non canonical kmer...'  >> ${{LOG}}; "
        "       singularity run {params.container} -t ${{NON_CANONICAL_OUTPUT_PREFIX}} "
        "       -s ${{NON_CANONICAL_TEL_KMER}} >> ${{LOG}} 2>&1 || true; "
        "       if [ -s '${{DESTDIR}}/ref.telomere' ]; "
        "       then"
        "           cd ${{STARTDIR}} >> ${{LOG}} 2>&1; "
        "           echo 'Sorting results...' >> ${{LOG}};  "
        "           sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{NON_CANONICAL_OUTPUT_PREFIX}}_telomere.bedgraph > {output.non_canonical_telo_track} 2>>${{LOG}}; "
        "           sort -k1,1V -k2,2n -k3,3n ${{DESTDIR}}/${{NON_CANONICAL_OUTPUT_PREFIX}}_telomere.bed > {output.non_canonical_telo_bed} 2>>${{LOG}}; "
        "           echo 'Copying files to the final destination...' >> ${{LOG}}; "
        "           cp ${{DESTDIR}}/ref.telomere {output.non_canonical_telo} >> ${{LOG}} 2>&1; "
        "           cp ${{DESTDIR}}/ref.windows {output.non_canonical_telo_win} >> ${{LOG}} 2>&1; "
        "       else"
        "           cd ${{STARTDIR}} >> ${{LOG}} 2>&1; "
        "           touch  {output.non_canonical_telo_track}  >> ${{LOG}} 2>&1; "
        "           touch  {output.non_canonical_telo_bed} >> ${{LOG}} 2>&1; "
        "           touch  {output.non_canonical_telo}  >> ${{LOG}} 2>&1; "
        "           touch  {output.non_canonical_telo_win} >> ${{LOG}} 2>&1; "
        "       fi; "
        "       rm -r ${{DESTDIR}}/* >> ${{LOG}} 2>&1; "
        " else "
        "       cd ${{STARTDIR}} >> ${{LOG}} 2>&1; "
        "       echo 'No non canonical kmer found. Creating empty output files...'  >> ${{LOG}};  "
        "       touch {output.non_canonical_telo_track} >> ${{LOG}} 2>&1; "
        "       touch {output.non_canonical_telo_bed} >> ${{LOG}} 2>&1; "
        "       touch {output.non_canonical_telo} >> ${{LOG}} 2>&1; "
        "       touch {output.non_canonical_telo_win} >> ${{LOG}}  2>&1; "
        " fi; "
        " echo 'Removing temporary workdir...'  >> ${{LOG}};  "
        " rm -r ${{WORKDIR}} >> ${{LOG}} 2>&1; "


rule collapse_overlapping_telomere_windows:
    input:
        canonical_telo_track = rules.telo_container.output.canonical_telo_track,
        non_canonical_telo_track = rules.telo_container.output.non_canonical_telo_track,
        #canonical_telo_track="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.bedgraph",
        #non_canonical_telo_track="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.bedgraph",
    output:
        canonical_collapsed_telo_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.bed",
        non_canonical_collapsed_telo_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.bed",
    params:
        container=config["tool_containers"]["rapid_telomere"]
    log:
        std="{fasta_dir}/collapse_telomere_track.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/collapse_telomere_track.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/collapse_telomere_track.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/collapse_telomere_track.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_container"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"],
    threads: parameters["threads"]["telo_finder"]

    shell:
        " echo -e '#scaffold\tstart\tend\tlength\tmedian\tmean\tstdev\tmode\tabsmin' > {output.canonical_collapsed_telo_track} 2>{log.std}; "
        " bedtools merge -c 4 -o median,mean,stdev,mode,absmin -i {input.canonical_telo_track} 2>>{log.std} | "
        " awk '{{print $1\"\t\"$2\"\t\"$3\"\t\"$3-$2\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$8}}' >> {output.canonical_collapsed_telo_track} 2>>{log.std}; "
        " echo -e '#scaffold\tstart\tend\tmedian\tmean\tstdev\tmode\tabsmin' > {output.non_canonical_collapsed_telo_track} 2>>{log.std}; "
        " bedtools merge -c 4 -o median,mean,stdev,mode,absmin -i {input.non_canonical_telo_track} | "
        " awk '{{print $1\"\t\"$2\"\t\"$3\"\t\"$3-$2\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$8}}' >> {output.non_canonical_collapsed_telo_track} 2>>{log.std}; "

rule classify_telomeric_regions_windows:
    input:
        canonical_collapsed_telo_track="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.bed",
        non_canonical_collapsed_telo_track="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.bed",
        fai="{fasta_dir}/{fasta_prefix}.fasta.fai"
    output:
        canonical_region_all_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.all.status",
        canonical_region_filtered_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        canonical_region_filtered_count="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        canonical_region_filtered_scaffold_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        non_canonical_region_all_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.all.status",
        non_canonical_region_filtered_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        non_canonical_region_filtered_count="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        non_canonical_region_filtered_scaffold_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
    params:
        fraction_threshold=parameters["tool_options"]["assembly_qc"]["telomere"]["fraction_threshold"]
    log:
        std="{fasta_dir}/classify_telomeric_regions_windows.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/classify_telomeric_regions_windows{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/classify_telomeric_regions_windows.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/classify_telomeric_regions_windows.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("telo_container"),
        cpus=parameters["threads"]["telo_finder"] ,
        time=parameters["time"]["telo_finder"],
        mem=parameters["memory_mb"]["telo_finder"],
    threads: parameters["threads"]["telo_finder"]

    shell:
        " CANONICAL_OUT_PREFIX={input.canonical_collapsed_telo_track}; "
        " CANONICAL_OUT_PREFIX=${{CANONICAL_OUT_PREFIX%.bed}}; "
        " NON_CANONICAL_OUT_PREFIX={input.non_canonical_collapsed_telo_track}; "
        " NON_CANONICAL_OUT_PREFIX=${{NON_CANONICAL_OUT_PREFIX%.bed}}; "
        " workflow/scripts/curation/classify_collapsed_telomere_track.py -i {input.canonical_collapsed_telo_track} "
        "          -p ${{CANONICAL_OUT_PREFIX}} -f {input.fai} -s {params.fraction_threshold} >{log.std} 2>&1; "

rule get_telomere_warning:
    input:
        canonical_telo_track=rules.telo_container.output.canonical_telo_track,
        non_canonical_telo_track=rules.telo_container.output.non_canonical_telo_track,
        fai="{fasta_dir}/{fasta_prefix}.fasta.fai",
    output:
        canonical_telo_warning_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere_warning.win1000.step200.track.bedgraph",
        non_canonical_telo_warning_track="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere_warning.win1000.step200.track.bedgraph",
        canonical_telo_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_telomere_warning.win1000.step200.track.status",
        non_canonical_telo_status="{fasta_dir}/telomere/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_telomere_warning.win1000.step200.track.status"
    params:
        fraction_threshold=parameters["tool_options"]["assembly_qc"]["telomere"]["fraction_threshold"]
    log:
        canonical="{fasta_dir}/get_telomere_warning.{fasta_prefix}.canonical.log",
        non_canonical="{fasta_dir}/get_telomere_warning.{fasta_prefix}.non_canonical.log",
        cluster_log="{fasta_dir}/get_telomere_warning.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/get_telomere_warning.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/get_telomere_warning.{fasta_prefix}.benchmark.txt"
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
        " CANONICAL_PREFIX={output.canonical_telo_warning_track}; "
        " CANONICAL_PREFIX=${{CANONICAL_PREFIX%.bedgraph}}; "
        " NON_CANONICAL_PREFIX={output.non_canonical_telo_warning_track}; "
        " NON_CANONICAL_PREFIX=${{NON_CANONICAL_PREFIX%.bedgraph}}; "
        " echo 'Checking positions of canonical telomeres...' > {log.canonical} 2>&1; "
        " if [ -s {input.canonical_telo_track} ]; "
        " then "
        "       workflow/scripts/curation/check_telomere.py  -i {input.canonical_telo_track} -f {input.fai} "
        "                  -s {params.fraction_threshold} -p ${{CANONICAL_PREFIX}} >> {log.canonical} 2>&1; "
        " else"
        "       touch {output.canonical_telo_warning_track} >> {log.canonical} 2>&1; "
        "       touch {output.canonical_telo_status} >> {log.canonical} 2>&1; "
        " fi;"
        " echo 'Checking positions of non canonical telomeres...' > {log.non_canonical} 2>&1; "
        " if [ -s {input.non_canonical_telo_track} ]; "
        " then "
        "       workflow/scripts/curation/check_telomere.py  -i {input.non_canonical_telo_track} -f {input.fai} "
        "                  -s {params.fraction_threshold} -p ${{NON_CANONICAL_PREFIX}}  >> {log.non_canonical} 2>&1; "
        " else"
        "       touch {output.non_canonical_telo_warning_track} >> {log.non_canonical} 2>&1; "
        "       touch {output.non_canonical_telo_status} >> {log.non_canonical} 2>&1; "
        " fi; "


rule copy_telomere_files:
    input:
        canonical_telo_track=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
        canonical_telo_warning_track=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere_warning.win1000.step200.track.bedgraph",
        non_canonical_telo_track=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.bedgraph",
        non_canonical_telo_warning_track=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere_warning.win1000.step200.track.bedgraph",
        canonical_region_all_status=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.all.status",
        canonical_region_filtered_status=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        canonical_region_filtered_count=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        canonical_region_filtered_scaffold_status=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        non_canonical_region_all_status=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.all.status",
        non_canonical_region_filtered_status=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        non_canonical_region_filtered_count=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        non_canonical_region_filtered_scaffold_status=out_dir_path / "{assembly_stage}/{parameters}/telomere/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
    output:
        canonical_telo_track = out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
        canonical_telo_warning_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere_warning.win1000.step200.track.bedgraph",
        non_canonical_telo_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.bedgraph",
        non_canonical_telo_warning_track=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere_warning.win1000.step200.track.bedgraph",
        canonical_region_all_status=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qctelomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.all.status",
        canonical_region_filtered_status=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        canonical_region_filtered_count=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        canonical_region_filtered_scaffold_status=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        non_canonical_region_all_status=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.all.status",
        non_canonical_region_filtered_status=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        non_canonical_region_filtered_count=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        non_canonical_region_filtered_scaffold_status=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/telomere/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",

    log:
        ln=out_dir_path / "{assembly_stage}/{parameters}/create_telomere_links.{genome_prefix}.{assembly_stage}.{haplotype}.ln.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/create_telomere_links.{genome_prefix}.{assembly_stage}.{haplotype}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/create_telomere_links.{genome_prefix}.{assembly_stage}.{haplotype}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/create_telomere_links.{genome_prefix}.{assembly_stage}.{haplotype}.benchmark.txt"
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
        " TRACK_DIR=`dirname {output.canonical_telo_track}`; "
        " TELOMERE_DIR=`dirname {output.canonical_region_all_status}`; "
        " cp {input.canonical_telo_track} {input.canonical_telo_warning_track} {input.non_canonical_telo_track} "
        "    {input.non_canonical_telo_warning_track} ${{TRACK_DIR}} > {log.ln} 2>&1; "
        " cp {input.canonical_region_all_status} {input.canonical_region_filtered_status} "
        "    {input.canonical_region_filtered_count} {input.canonical_region_filtered_scaffold_status} ${{TELOMERE_DIR}} > {log.ln} 2>&1; "
        " cp {input.non_canonical_region_all_status} {input.non_canonical_region_filtered_status} "
        "    {input.non_canonical_region_filtered_count} {input.non_canonical_region_filtered_scaffold_status} ${{TELOMERE_DIR}} > {log.ln} 2>&1; "