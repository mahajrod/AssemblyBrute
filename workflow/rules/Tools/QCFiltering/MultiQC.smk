
rule multiqc:
    input:
        fastqc_reports=lambda wildcards: expand(config["out_dir"] / ("qc/fastqc/%s/%s/{fileprefix}_fastqc.zip" % (wildcards.fastqc_datatype,
                                                                                                                            wildcards.stage)),
                                 fileprefix=config["data"][wildcards.fastqc_datatype]["conv_file_prefix_list"],
                                 allow_missing=True)
    output:
        dir=directory(config["out_dir"] / "qc/multiqc/{fastqc_datatype}/{stage}/"),
        report=config["out_dir"] / "qc/multiqc/{fastqc_datatype}/{stage}/multiqc.{fastqc_datatype}.{stage}.report.html"
    params:
        # multiqc adds report filename to outdir path and even creates additional subdirectories if necessary.
        # So if you set --outdir option --filename should not contain directories.
        # Moreover, --filename is in fact not filename but prefix
        input_dir=lambda wildcards: config["out_dir"] / "qc/fastqc/{0}/{1}/".format(wildcards.fastqc_datatype,
                                                                                 wildcards.stage)
    log:
        std=config["out_dir"] / "log/multiqc.{fastqc_datatype}.{stage}.log",
        cluster_log=config["out_dir"] / "log/multiqc.{fastqc_datatype}.{stage}.cluster.log",
        cluster_err=config["out_dir"] / "log/multiqc.{fastqc_datatype}.{stage}.cluster.err"
    benchmark:
        config["out_dir"] / "log/multiqc.{fastqc_datatype}.{stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("multiqc"),
        cpus=parameters["threads"]["multiqc"],
        time=parameters["time"]["multiqc"],
        mem=parameters["memory_mb"]["multiqc"],
    threads:
        parameters["threads"]["multiqc"]
    shell:
        " REPORT_PREFIX={output.report}; "
        " REPORT_PREFIX=`basename ${{REPORT_PREFIX%.html}}`; "
        " OUTDIR=`dirname {output.report}`; "
        " multiqc --filename ${{REPORT_PREFIX}} -p --outdir ${{OUTDIR}} --comment {wildcards.fastqc_datatype} {params.input_dir} > {log.std} 2>&1; "
