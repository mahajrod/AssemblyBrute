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
        se_fastq=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (wildcards.datatype,
                                                                                                    config["fastq_extension"])),
                    fileprefix=input_file_prefix_dict[wildcards.datatype],
                    allow_missing=True) if wildcards.datatype not in config["paired_fastq_based_data"] else [],
        forward_fastq=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{pairprefix}_1%s" % (wildcards.datatype,
                                                                                                              config["fastq_extension"])),
                                               pairprefix=input_pairprefix_dict[wildcards.datatype],
                                               allow_missing=True) if wildcards.datatype in config["paired_fastq_based_data"] else [],
        reverse_fastq=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{pairprefix}_2%s" % (wildcards.datatype,
                                                                                                              config["fastq_extension"])),
                                               pairprefix=input_pairprefix_dict[wildcards.datatype],
                                               allow_missing=True) if wildcards.datatype in config["paired_fastq_based_data"] else [],
        #reverse_fastq=lambda wildcards: output_dict["data"] / ("fastq/{datatype}/filtered/{pairprefix}_2%s" % config["fastq_extension"]) if wildcards.datatype in config["paired_fastq_based_data"] else [],
        db=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["path"]
    output:
        summary=out_dir_path / "contamination_scan/kraken2/{datatype}/kraken2.{database}.report",
        out=out_dir_path / "contamination_scan/kraken2/{datatype}/kraken2.{database}.out.gz",
    params:
        forward_fastq=lambda wildcards: " <(cat {0} )".format(" ".join(expand(output_dict["data"] / ("fastq/%s/filtered/{pairprefix}_1%s" % (wildcards.datatype,
                                                                                                              config["fastq_extension"])),
                                                                       pairprefix=input_pairprefix_dict[wildcards.datatype],
                                                                       allow_missing=True) )) if wildcards.datatype in config["paired_fastq_based_data"] else "",
        reverse_fastq=lambda wildcards: " <(cat {0} )".format(" ".join(expand(output_dict["data"] / ("fastq/%s/filtered/{pairprefix}_2%s" % (wildcards.datatype,
                                                                                                              config["fastq_extension"])),
                                                                       pairprefix=input_pairprefix_dict[wildcards.datatype],
                                                                       allow_missing=True) )) if wildcards.datatype in config["paired_fastq_based_data"] else "",
        memory_mapping=lambda wildcards: "" if config["allowed_databases"]["kraken2"][wildcards.database]["in_memory"] else  " --memory-mapping ",
        paired=lambda wildcards: " --paired " if wildcards.datatype in config["paired_fastq_based_data"] else "",
        compressed=get_kraken2_compression_flag_from_extension(config["fastq_extension"])
        #forward_fastq=lambda wildcards: " <(cat {0}) ".format(" ".join(input_filedict[wildcards.datatype][::2])) if wildcards.datatype in config["paired_fastq_based_data"] else "",
        #reverse_fastq=lambda wildcards: " <(cat {0}) ".format(" ".join(input_filedict[wildcards.datatype][1::2])) if wildcards.datatype in config["paired_fastq_based_data"] else "",
    log:
        std=output_dict["log"]  / "kraken2.{database}.{datatype}.log",
        pigz=output_dict["log"]  / "kraken2.{database}.{datatype}.pigz.log",
        cluster_log=output_dict["cluster_log"] / "kraken2.{database}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "kraken2.{database}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "kraken2.{database}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("kraken2"),
        cpus=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["threads"] ,
        time=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["time"] ,
        mem=lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["memory_mb"] ,
    threads: lambda wildcards: config["allowed_databases"]["kraken2"][wildcards.database]["threads"] ,

    shell:
        " OUT_FILE={output.out}; "
        " kraken2 --threads {threads} {params.memory_mapping} {params.paired} {params.compressed}  --db {input.db} "
        " --output ${{OUT_FILE%.gz}} --report {output.summary} "
        " {input.se_fastq} {params.forward_fastq} {params.reverse_fastq} > {log.std} 2>&1;"
        " pigz -p {threads} ${{OUT_FILE%.gz}} > {log.pigz} 2>&1"

