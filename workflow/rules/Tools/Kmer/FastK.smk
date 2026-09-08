
"""
rule fastk_assembly:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        db_dir=directory("{fasta_dir}/{fasta_prefix}/kmer/meryl/{fasta_prefix}.{kmer_length}.meryl/")
    log:
        std="{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.log",
        cluster_log="{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.cluster.log",
        cluster_err="{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.meryl.cluster.err"
    benchmark:
        "{fasta_dir}/log/meryl_assembly.{fasta_prefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_assembly"),
        cpus=get_threads(parameters["threads"]["meryl_assembly"], "cpu"),
        time=parameters["time"]["meryl_assembly"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl_assembly"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl_assembly"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input.fasta} > {log.std} 2>&1;"

rule meryl_get_repetitive_kmers:
    input:
        meryl_db="{directory}/{meryl_db_prefix}.meryl/",
        log_dir=ancient("{directory}/log/"),
    output:
        repetitive_kmers="{directory}/{meryl_db_prefix}.meryl.repetitive"
    params:
        distinct=0.9998
    log:
        std="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.log",
        cluster_log="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.cluster.log",
        cluster_err="{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.meryl.cluster.err"
    benchmark:
        "{directory}/log/meryl_get_repetitive_kmers.{meryl_db_prefix}.meryl.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_extract"),
        cpus=parameters["threads"]["meryl_extract"],
        time=parameters["time"]["meryl_extract"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["meryl_extract"],
        kmer_counter=1
    threads:
        parameters["threads"]["meryl_extract"]
    shell:
         " meryl print greater-than distinct={params.distinct} "
         " {input.meryl_db} > {output.repetitive_kmers} 2>{log.std}; "

"""

rule fastk_se:
    input:
        lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}".format(wildcards.se_datatype,
                                                                               wildcards.stage,
                                                                               wildcards.fileprefix,
                                                                               config["data"][wildcards.se_datatype]["conv_ext"])

    output:
        db=directory(config["out_dir"] / "kmer/{se_datatype}/{stage}/{se_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{fileprefix}/"),
        #ktab=config["out_dir"] / "kmer/{se_datatype}/{stage}/{se_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{fileprefix}/fastk_db.ktab",
    log:
        std=config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.log",
        cluster_log=config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk"),
        cpus=get_threads(parameters["threads"]["fastk"], "cpu"),
        time=parameters["time"]["fastk"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["fastk"],
        kmer_counter=1
    threads:
        parameters["threads"]["fastk"]
    shell:
         " MEM_GB=`echo '{resources.mem}/1024' | bc`; "
         " TMP_DIR={output.db}/tmp_`basename {output.db}`; "
         " mkdir -p {output.db} ${{TMP_DIR}}; "
         " FastK -v -t{wildcards.min_kmer_count} -k{wildcards.kmer_length} -M${{MEM_GB}} -T{threads} "
         "       -P${{TMP_DIR}} -N{output.db}/fastk_db {input}  > {log.std} 2>&1; "
         " rm -r ${{TMP_DIR}}; "



rule fastk_pe:
    input:
        forward_fastq=lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                wildcards.stage,
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.pe_datatype]["conv_fwd_sfx"],
                                                                                                config["data"][wildcards.pe_datatype]["conv_ext"]),
        reverse_fastq=lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                wildcards.stage,
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.pe_datatype]["conv_rev_sfx"],
                                                                                                config["data"][wildcards.pe_datatype]["conv_ext"]),
    output:
        db=directory(config["out_dir"] / "kmer/{pe_datatype}/{stage}/{pe_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{pairprefix}/"),
        #ktab=config["out_dir"] / "kmer/{pe_datatype}/{stage}/{pe_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{pairprefix}/fastk_db.ktab"
    log:
        std=config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.log",
        cluster_log=config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk"),
        cpus=get_threads(parameters["threads"]["fastk"], "cpu"),
        time=parameters["time"]["fastk"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["fastk"],
        kmer_counter=1
    threads:
        parameters["threads"]["fastk"]
    shell:
         " MEM_GB=`echo '{resources.mem}/1024' | bc`; "
         " TMP_DIR={output.db}/tmp_`basename {output.db}`; "
         " mkdir -p {output.db} ${{TMP_DIR}}; "
         " FastK -v -t{wildcards.min_kmer_count} -k{wildcards.kmer_length} -M${{MEM_GB}} -T{threads} "
         "       -P${{TMP_DIR}} -N{output.db}/fastk_db {input}  > {log.std} 2>&1; "
         " rm -r ${{TMP_DIR}}; "

def get_fastk_dbs_for_merging(wildcards):
    db_list = []
    for datatype in wildcards.datatype.split("_"):
        if datatype in config["data_feature_dict"]["paired"]:
            db_list += expand(rules.fastk_pe.output.db,
                              pairprefix=config["data"][datatype]["pair_prefix_list"],
                              pe_datatype=[datatype,],
                              allow_missing=True)
        else:
            db_list += expand(rules.fastk_se.output.db,
                              fileprefix=config["data"][datatype]["conv_file_prefix_list"], # for se_reads "conv_file_prefix_list" and "file_prefix_list" are the same
                              se_datatype=[datatype,],
                              allow_missing=True)

    return db_list

rule merge_fastk:
    input: get_fastk_dbs_for_merging
    output:
        db=directory(config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/"),
        #ktab=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/fastk_db.ktab",
    log:
        std=config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.std.log",
        cluster_log=config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk"),
        cpus=get_threads(parameters["threads"]["fastk"], "cpu"),
        time=parameters["time"]["fastk"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["fastk"],
        kmer_counter=1
    threads:
        parameters["threads"]["fastk"]
    shell: # Fastmerge return an error is input is a single database
         " INPUT_DB_ARRAY=({input}); "
         " if [ ${{#INPUT_DB_ARRAY[@]}} == '1' ] ;"
         " then "
         "     INPUT_DB_DIR=`dirname {input}`; "
         "     ln -sf `basename {input}` {output.db}; "
         " else "
         "      INPUT_DB_ARRAY=(\"${{INPUT_DB_ARRAY[@]/%/fastk_db.ktab}}\"); "
         "      TMP_DIR={output.db}/tmp_`basename {output.db}`; "
         "      mkdir -p {output.db} ${{TMP_DIR}}; "
         "      Fastmerge -t -T{threads} ${{INPUT_DB_ARRAY}} {output.db}/fastk_db.ktab > {log.std} 2>&1; "
         "      Fastmerge -h -T{threads} ${{INPUT_DB_ARRAY}} {output.db}/fastk_db.hist > {log.std} 2>&1; "
         "      rm -r ${{TMP_DIR}}; "
         " fi; "



rule get_fastk_histo:
    input:
        db="{fastk_dir}/{fastk_db_prefix}.fastk_min{min_kmer_count}/",
        log_dir=ancient("{fastk_dir}/log/")
    output:
        histo="{fastk_dir}/{fastk_db_prefix}.fastk_min{min_kmer_count}.histo"

    log:
        histo_log="{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.log",
        cluster_log="{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.cluster.log",
        cluster_err="{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.cluster.err"
    benchmark:
        "{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk_histo"),
        cpus=parameters["threads"]["fastk_histo"],
        time=parameters["time"]["fastk_histo"],
        mem=parameters["memory_mb"]["fastk_histo"],
    threads:
        parameters["threads"]["fastk_histo"]
    shell:
         " Histex -G -h32637 {input.db}/fastk_db.hist > {output.histo} 2>{log.histo_log}"

#ruleorder: create_final_fastk_histo_link > get_fastk_histo

use rule create_local_links as create_final_fastk_db_link with:
    input:
        input=lambda wildcards: config["out_dir"] / ("kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min%s/" % (parameters["tool_options"]["fastk"][wildcards.datatype]["min_kmer_count"])),
        log_dir=ancient(config["out_dir"] / "kmer/log/")
    output:
        input=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk"
    log:
        ln=config["out_dir"] / "kmer/create_final_fastk_db_link.{datatype}.{stage}.{stage}.{kmer_length}.fastk.ln.log",



use rule create_local_links as create_final_fastk_histo_link with:
    input:
        input=lambda wildcards: config["out_dir"] / ("kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min%s.histo" % (parameters["tool_options"]["fastk"][wildcards.datatype]["min_kmer_count"])),
        log_dir=ancient(config["out_dir"] / "kmer/log/")
    output:
        input=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk.histo"
    log:
        ln=config["out_dir"] / "kmer/create_final_fastk_histo_link.{datatype}.{stage}.{stage}.{kmer_length}.fastk.ln.log",


"""
rule fastk_se:
    input:
        lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}".format(wildcards.se_datatype,
                                                                               wildcards.stage,
                                                                               wildcards.fileprefix,
                                                                               config["data"][wildcards.se_datatype]["conv_ext"])

    output:
        hist=config["out_dir"] / "kmer/{se_datatype}/{stage}/{se_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{fileprefix}/fastk_db.hist",
        ktab=config["out_dir"] / "kmer/{se_datatype}/{stage}/{se_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{fileprefix}/fastk_db.ktab",
    log:
        std=config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.log",
        cluster_log=config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fastk_se.{se_datatype}.{stage}.{fileprefix}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk"),
        cpus=get_threads(parameters["threads"]["fastk"], "cpu"),
        time=parameters["time"]["fastk"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["fastk"],
        kmer_counter=1
    threads:
        parameters["threads"]["fastk"]
    shell:
         " MEM_GB=`echo '{resources.mem}/1024' | bc`; "
         " FILE_PREFIX=`basename {output.ktab}`; "
         " FILE_PREFIX=${{FILE_PREFIX%.ktab}}; "
         " DB_DIR=`dirname {output.ktab}`; "
         " TMP_DIR=${{DB_DIR}}/tmp_${{FILE_PREFIX}}; "
         " mkdir -p ${{TMP_DIR}}; "
         " FastK -v -t{wildcards.min_kmer_count} -k{wildcards.kmer_length} -M${{MEM_GB}} -T{threads} "
         "       -P${{TMP_DIR}} -N${{DB_DIR}}/${{FILE_PREFIX}} {input}  > {log.std} 2>&1; "
         " rm -r ${{TMP_DIR}}; "



rule fastk_pe:
    input:
        forward_fastq=lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                wildcards.stage,
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.pe_datatype]["conv_fwd_sfx"],
                                                                                                config["data"][wildcards.pe_datatype]["conv_ext"]),
        reverse_fastq=lambda wildcards: config["out_dir"] / "data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                wildcards.stage,
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.pe_datatype]["conv_rev_sfx"],
                                                                                                config["data"][wildcards.pe_datatype]["conv_ext"]),
    output:
        hist=config["out_dir"] / "kmer/{pe_datatype}/{stage}/{pe_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{pairprefix}/fastk_db.hist",
        ktab=config["out_dir"] / "kmer/{pe_datatype}/{stage}/{pe_datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.{pairprefix}/fastk_db.ktab"
    log:
        std=config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.log",
        cluster_log=config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fastk_pe.{pe_datatype}.{stage}.{pairprefix}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk"),
        cpus=get_threads(parameters["threads"]["fastk"], "cpu"),
        time=parameters["time"]["fastk"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["fastk"],
        kmer_counter=1
    threads:
        parameters["threads"]["fastk"]
    shell:
         " MEM_GB=`echo '{resources.mem}/1024' | bc`; "
         " FILE_PREFIX=`basename {output.ktab}`; "
         " FILE_PREFIX=${{FILE_PREFIX%.ktab}}; "
         " DB_DIR=`dirname {output.ktab}`; "
         " TMP_DIR=${{DB_DIR}}/tmp_${{FILE_PREFIX}}; "
         " mkdir -p ${{TMP_DIR}}; "
         " FastK -v -t{wildcards.min_kmer_count} -k{wildcards.kmer_length} -M${{MEM_GB}} -T{threads} "
         "       -P${{TMP_DIR}} -N${{DB_DIR}}/${{FILE_PREFIX}} {input}  > {log.std} 2>&1; "
         " rm -r ${{TMP_DIR}}; "

def get_fastk_dbs_for_merging(wildcards):
    db_list = []
    for datatype in wildcards.datatype.split("_"):
        if datatype in config["data_feature_dict"]["paired"]:
            db_list += expand(rules.fastk_pe.output.ktab,
                              pairprefix=config["data"][datatype]["pair_prefix_list"],
                              pe_datatype=[datatype,],
                              allow_missing=True)
        else:
            db_list += expand(rules.fastk_se.output.ktab,
                              fileprefix=config["data"][datatype]["conv_file_prefix_list"], # for se_reads "conv_file_prefix_list" and "file_prefix_list" are the same
                              se_datatype=[datatype,],
                              allow_missing=True)

    return db_list

rule merge_fastk:
    input: get_fastk_dbs_for_merging
    output:
        hist=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/fastk_db.hist",
        ktab=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/fastk_db.ktab",
    log:
        std=config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.std.log",
        cluster_log=config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/merge_fastk.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk"),
        cpus=get_threads(parameters["threads"]["fastk"], "cpu"),
        time=parameters["time"]["fastk"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["fastk"],
        kmer_counter=1
    threads:
        parameters["threads"]["fastk"]
    shell: # Fastmerge return an error is input is a single database
         " FILE_PREFIX=`basename {output.ktab}`; "
         " FILE_PREFIX=${{FILE_PREFIX%.ktab}}; "
         " DB_DIR=`dirname {output.ktab}`; "
         " TMP_DIR=${{DB_DIR}}/tmp_${{FILE_PREFIX}}; "
         " mkdir -p ${{TMP_DIR}}; "
         " INPUT_ARRAY=({input}); "
         " if [ ${{#INPUT_ARRAY[@]}} == '1' ] ;"
         " then "
         "     INPUT_DB_DIR=`dirname {input}`; "
         "     ln -sf `basename ${{INPUT_DB_DIR}}` ${{DB_DIR}}; "
         " else "
         "      Fastmerge -t -T{threads} ${{DB_DIR}}/${{FILE_PREFIX}} {input} {output.ktab} > {log.std} 2>&1; "
         "      Fastmerge -h -T{threads} ${{DB_DIR}}/${{FILE_PREFIX}} {input} {output.hist} > {log.std} 2>&1; "
         " fi; "
         " rm -r ${{TMP_DIR}}; "


rule get_fastk_histo:
    input:
        hist="{fastk_dir}/{fastk_db_prefix}.fastk_min{min_kmer_count}/fastk_db.hist",
        log_dir=ancient("{fastk_dir}/log/")
    output:
        histo="{fastk_dir}/{fastk_db_prefix}.fastk_min{min_kmer_count}.histo"

    log:
        histo_log="{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.log",
        cluster_log="{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.cluster.log",
        cluster_err="{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.cluster.err"
    benchmark:
        "{fastk_dir}/log/get_fastk_histo.{fastk_db_prefix}.fastk_min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastk_histo"),
        cpus=parameters["threads"]["fastk_histo"],
        time=parameters["time"]["fastk_histo"],
        mem=parameters["memory_mb"]["fastk_histo"],
    threads:
        parameters["threads"]["fastk_histo"]
    shell:
         " Histex -G -h32637 {input.hist} > {output.histo} 2>{log.histo_log}"

#ruleorder: create_final_fastk_histo_link > get_fastk_histo

use rule create_local_links as create_final_fastk_db_link with:
    input:
        input=lambda wildcards: config["out_dir"] / ("kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min%s/" % (parameters["tool_options"]["fastk"][wildcards.datatype]["min_kmer_count"])),
        log_dir=ancient(config["out_dir"] / "kmer/log/")
    output:
        input=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk"
    log:
        ln=config["out_dir"] / "kmer/create_final_fastk_db_link.{datatype}.{stage}.{stage}.{kmer_length}.fastk.ln.log",



use rule create_local_links as create_final_fastk_histo_link with:
    input:
        input=lambda wildcards: config["out_dir"] / ("kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min%s.histo" % (parameters["tool_options"]["fastk"][wildcards.datatype]["min_kmer_count"])),
        log_dir=ancient(config["out_dir"] / "kmer/log/")
    output:
        input=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk.histo"
    log:
        ln=config["out_dir"] / "kmer/create_final_fastk_histo_link.{datatype}.{stage}.{stage}.{kmer_length}.fastk.ln.log",

"""