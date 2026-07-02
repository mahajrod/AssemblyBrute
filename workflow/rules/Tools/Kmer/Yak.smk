
def get_files_for_yak(wildcards):
    file_list = []
    for datatype in wildcards.datatype.split("_"):
        file_list += expand(config["out_dir"] / ("data/%s/%s/{fileprefix}%s" % (datatype,
                                                                                          wildcards.stage,
                                                                                          config["data"][wildcards.datatype]["conv_ext"])),
                            fileprefix=config["data"][wildcards.datatype]["conv_file_prefix_list"],
                            allow_missing=True)

    return file_list


rule yak_reads:
    # yak reads data twice, so if input is a stream, it should be supplied twice!!!!!!!! for pe reads it also means two identical streams!!!!
    # so just combine all the files in a single list and supply it to yak as a stream twice
    input: get_files_for_yak
    output:
        db=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.yak"
    log:
        std=config["out_dir"] / "log/yak_reads.{datatype}.{stage}.{kmer_length}.log",
        cluster_log=config["out_dir"] / "log/yak_reads.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/yak_reads.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "log/yak_reads.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("yak_reads"),
        cpus=get_threads(parameters["threads"]["yak_reads"], "cpu"),
        time=parameters["time"]["yak_reads"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["yak_reads"],
        yak_kmer_counter=1
    threads:
        parameters["threads"]["yak_reads"]
    shell:
         " yak count -o {output} -k {wildcards.kmer_length} -b 37 <(zcat {input}) <(zcat {input}) > {log.std} 2>&1; "

rule yak_histo:
    input:
        db="{yak_db_dir}/{yak_db_prefix}.yak",
        log_dir=ancient("{yak_db_dir}/log/")
    output:
        original_histo="{yak_db_dir}/{yak_db_prefix}.yak.histo.original",
        histo="{yak_db_dir}/{yak_db_prefix}.yak.histo",
    log:
        std="{yak_db_dir}/log/yak_histo.{yak_db_prefix}.yak.log",
        cluster_log="{yak_db_dir}/log/yak_histo.{yak_db_prefix}.yak.cluster.log",
        cluster_err="{yak_db_dir}/log/yak_histo.{yak_db_prefix}.yak.cluster.err"
    benchmark:
        "{yak_db_dir}/log/yak_histo.{yak_db_prefix}.yak.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("yak_reads"),
        cpus=get_threads(parameters["threads"]["yak_reads"], "cpu"),
        time=parameters["time"]["yak_reads"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["yak_reads"],
    threads:
        parameters["threads"]["yak_reads"]
    shell:
         " yak inspect {input.db} > {output.original_histo} 2>{log.std}; "
         " cut -f 2,4 {output.original_histo} | tac | tail -n +3 > {output.histo} 2>>{log.std}; "
