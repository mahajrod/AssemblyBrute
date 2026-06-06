
rule porechop_abi:
    input:
        fastq=output_dict["data"] / ("fastq/{datatype}/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        trimmed_fastq=output_dict["data"] / ("fastq/{datatype, nanopore|simplex|duplex|ultralongnano|adaptivenano}/trimmed/{fileprefix, [^/]+}%s" % config["fastq_extension"]),
    params:
         ab_initio=lambda wildcards: parse_option_flag("ab_initio", parameters["tool_options"]["porechop_abi"]["nanopore"], "--ab_initio"),
         verbosity=lambda wildcards: parse_option("verbosity", parameters["tool_options"]["porechop_abi"]["nanopore"], "-v"),
    log:
        porechop_abi=output_dict["log"]/ "porechop_abi.{datatype}.trimmed.{fileprefix}.porechop_abi.log",
        cluster_log=output_dict["cluster_log"] / "porechop_abi.{datatype}.trimmed.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "porechop_abi.{datatype}.trimmed.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "porechop_abi.{datatype}.trimmed.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("porechop_abi"),
        cpus=parameters["threads"]["porechop_abi"],
        time=parameters["time"]["porechop_abi"],
        mem=parameters["memory_mb"]["porechop_abi"],
    threads:
        parameters["threads"]["porechop_abi"],
    shell:
        " porechop_abi {params.ab_initio} {params.verbosity} -t {threads} "
        " -i {input.fastq} -o {output.trimmed_fastq} 1>{log.porechop_abi} 2>&1; "


rule chopper:
    input:
        input_fastq=output_dict["data"] / ("fastq/{datatype}/%s/{fileprefix}%s" % ("trimmed" if not config["skip_porechop_abi"] else "raw",
                                                                                   config["fastq_extension"]))
    output:
        filtered_fastq=output_dict["data"] / ("fastq/{datatype, nanopore|simplex|duplex|ultralongnano|adaptivenano}/filtered/{fileprefix, [^/]+}%s" % config["fastq_extension"]),
    params:
         headcrop  = parse_option("headcrop",  parameters["tool_options"]["chopper"]["nanopore"], "--headcrop"),
         maxlength = parse_option("maxlength", parameters["tool_options"]["chopper"]["nanopore"], "--maxlength"),
         minlength = parse_option("minlength", parameters["tool_options"]["chopper"]["nanopore"], "--minlength"),
         quality   = parse_option("quality",   parameters["tool_options"]["chopper"]["nanopore"], "--quality"),
         tailcrop  = parse_option("tailcrop",  parameters["tool_options"]["chopper"]["nanopore"], "--tailcrop"),
    log:
        zcat=output_dict["log"]/ "chopper.{datatype}.filtered.{fileprefix}.zcat.log",
        chopper=output_dict["log"]/ "chopper.{datatype}.filtered.{fileprefix}.chopper.log",
        pigz_filtered=output_dict["log"]/ "chopper.{datatype}.filtered.{fileprefix}.pigz_filtered.log",
        cluster_log=output_dict["cluster_log"] / "chopper.{datatype}.filtered.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "chopper.{datatype}.filtered.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "chopper.{datatype}.filtered.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("chopper"),
        cpus=parameters["threads"]["chopper"],
        time=parameters["time"]["chopper"],
        mem=parameters["memory_mb"]["chopper"],
    threads:
        parameters["threads"]["chopper"]
    shell:
        " zcat {input} 2>{log.zcat} | "
        " chopper {params.headcrop} {params.maxlength} {params.minlength} {params.quality} {params.tailcrop} "
        " -t {threads} 2>{log.chopper} | pigz -p 4 > {output.filtered_fastq} 2>{log.pigz_filtered} ; "


