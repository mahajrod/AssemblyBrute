
rule trim_hic: # trims arima hic files from 5', but creates links to raw files if hic is of other type
    priority: 2000
    input:
        fastq=config["out_dir"] / ("data/hic/raw/{fileprefix}%s" % config["fastq_ext"])
    output:
        fastq=config["out_dir"] / ("data/hic/trimmed/{fileprefix}%s" % config["fastq_ext"])
    params:
        hic_type=config["hic_enzyme_set"],
    log:
        std=config["out_dir"] / "log/trim_hic.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/trim_hic.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/trim_hic.{fileprefix}.cluster.err",
    benchmark:
        config["out_dir"] / "log/trim_hic.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("preprocess_hic_fastq"),
        cpus=parameters["threads"]["preprocess_hic_fastq"],
        time=parameters["time"]["preprocess_hic_fastq"],
        mem=parameters["memory_mb"]["preprocess_hic_fastq"],
    threads:
        parameters["threads"]["preprocess_hic_fastq"]
    shell:
         " if [ '{params.hic_type}' == 'Arima' ]; "
         " then "
         "     echo 'Input is Arima Hi-C! Trimming first 8 bp from both forward and reverse reads...' >> {log.std}; "
         "     zcat {input} | fastx_trimmer -f 8 | pigz -p {threads} > {output.fastq} 2>>{log.std}; "
         " else "
         "     echo 'Input is not Arima Hi-C! Trimming skipped...' >> {log.std}; "
         "     ln -sf {input} {output.fastq} 2>>{log.std}; "
         " fi; "