
rule multiqc:
    input:
        #output_dict["qc"] / "fastqc/{datatype}/{stage}/",
        fastqc_reports=lambda wildcards: expand(output_dict["qc"] / ("fastqc/%s/%s/{fileprefix}_fastqc.zip" % (wildcards.datatype,
                                                                                                               wildcards.stage)),
                                 fileprefix=input_file_prefix_dict[wildcards.datatype],
                                 allow_missing=True)
    output:
        dir=directory(output_dict["qc"] / "multiqc/{datatype}/{stage}/"),
        report=output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html"
        #stats=merged_raw_multiqc_dir_path / "{library_id}/{library_id}.raw.multiqc.stats"
    params:
        # multiqc adds report filename to outdir path and even creates additional subdirectories if necessary.
        # So if you set --outdir option --filename should not contain directories.
        # Moreover, --filename is in fact not filename but prefix
        #report_filename=lambda wildcards: "multiqc.{0}.{1}.report".format(wildcards.datatype,
        #                                                                  wildcards.stage),
        input_dir=lambda wildcards: output_dict["qc"] / "fastqc/{0}/{1}/".format(wildcards.datatype,
                                                                                 wildcards.stage)
    log:
        std=output_dict["log"]/ "multiqc.{datatype}.{stage}.log",
        #stats=log_dir_path / "{library_id}/multiqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"]/ "multiqc.{datatype}.{stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "multiqc.{datatype}.{stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "multiqc.{datatype}.{stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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
        " multiqc --filename ${{REPORT_PREFIX}} -p --outdir ${{OUTDIR}} "
        " --comment {wildcards.datatype} {params.input_dir} > {log.std} 2>&1; "