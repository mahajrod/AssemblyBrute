def get_kraken2_compression_flag_from_extension(extension):
    if len(extension) > 2:
        if extension[-2:] == "gz":
            return " --gzip-compressed "
    elif len(extension) > 3:
        if extension[-3:] == "bz2":
            return " --bzip2-compressed "
    else:
        return " "


rule kraken2: #
    input:
        se_fastq=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{fileprefix}%s" % (wildcards.datatype,
                                                                                                           config["data"][wildcards.datatype]["conv_ext"])),
                    fileprefix=config["data"][wildcards.datatype]["conv_file_prefix_list"],
                    allow_missing=True) if wildcards.datatype not in config["data_feature_dict"]["paired"] else [],
        forward_fastq=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{pairprefix}%s%s" % (wildcards.datatype,
                                                                                                                  config["data"][wildcards.datatype]["conv_fwd_sfx"],
                                                                                                                  config["data"][wildcards.datatype]["conv_ext"])),
                                               pairprefix=config["data"][wildcards.datatype]["pair_prefix_list"],
                                               allow_missing=True) if wildcards.datatype in config["data_feature_dict"]["paired"] else [],
        reverse_fastq=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{pairprefix}%s%s" % (wildcards.datatype,
                                                                                                                  config["data"][wildcards.datatype]["conv_rev_sfx"],
                                                                                                                  config["data"][wildcards.datatype]["conv_ext"])),
                                               pairprefix=config["data"][wildcards.datatype]["pair_prefix_list"],
                                               allow_missing=True) if wildcards.datatype in config["data_feature_dict"]["paired"] else [],
        db=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["path"]
    output:
        summary=config["out_dir"] / "contamination_scan/kraken2/{datatype}/kraken2.{database}.report",
        out=config["out_dir"] / "contamination_scan/kraken2/{datatype}/kraken2.{database}.out.gz",
    params:
        forward_fastq=lambda wildcards: " <(cat {0} )".format(" ".join(expand(config["out_dir"] / ("data/%s/filtered/{pairprefix}%s%s" % (wildcards.datatype,
                                                                                                                                                    config["data"][wildcards.datatype]["conv_fwd_sfx"],
                                                                                                                                                    config["data"][wildcards.datatype]["conv_ext"])),
                                                                       pairprefix=input_pairprefix_dict[wildcards.datatype],
                                                                       allow_missing=True) )) if wildcards.datatype in config["data_feature_dict"]["paired"] else "",
        reverse_fastq=lambda wildcards: " <(cat {0} )".format(" ".join(expand(config["out_dir"] / ("data/%s/filtered/{pairprefix}%s%s" % (wildcards.datatype,
                                                                                                                                                    config["data"][wildcards.datatype]["conv_rev_sfx"],
                                                                                                                                                    config["data"][wildcards.datatype]["conv_ext"])),
                                                                       pairprefix=input_pairprefix_dict[wildcards.datatype],
                                                                       allow_missing=True) )) if wildcards.datatype in config["data_feature_dict"]["paired"] else "",
        memory_mapping=lambda wildcards: "" if config["allowed_databases"]["kraken2"][wildcards.database]["in_memory"] else  " --memory-mapping ",
        paired=lambda wildcards: " --paired " if wildcards.datatype in config["data_feature_dict"]["paired"] else "",
        compressed=lambda wildcards: get_kraken2_compression_flag_from_extension(config["data"][wildcards.datatype]["conv_ext"])
    log:
        std=config["out_dir"] / "log/kraken2.{database}.{datatype}.log",
        pigz=config["out_dir"] / "log/kraken2.{database}.{datatype}.pigz.log",
        cluster_log=config["out_dir"] / "log/kraken2.{database}.{datatype}.cluster.log",
        cluster_err=config["out_dir"] / "log/kraken2.{database}.{datatype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/kraken2.{database}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("kraken2"),
        cpus=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["threads"] ,
        time=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["time"] ,
        mem=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["memory_mb"] ,
    threads: lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["threads"] ,

    shell:
        " OUT_FILE={output.out}; "
        " kraken2 --threads {threads} {params.memory_mapping} {params.paired} {params.compressed}  --db {input.db} "
        "     --report-minimizer-data --output ${{OUT_FILE%.gz}} --report {output.summary} "
        "     {input.se_fastq} {params.forward_fastq} {params.reverse_fastq} > {log.std} 2>&1;"
        " pigz -p {threads} ${{OUT_FILE%.gz}} > {log.pigz} 2>&1"
