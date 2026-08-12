
rule porechop_abi:
    input:
        fastq=config["out_dir"] / ("data/{nanopore_datatype}/raw/{fileprefix}%s" % config["fastq_ext"])
    output:
        trimmed_fastq=config["out_dir"] / ("data/{nanopore_datatype}/trimmed/{fileprefix}%s" % config["fastq_ext"]),
    params:
         ab_initio=lambda wildcards: parse_option_flag("ab_initio", parameters["tool_options"]["porechop_abi"]["nanopore"], "--ab_initio"),
         verbosity=lambda wildcards: parse_option("verbosity", parameters["tool_options"]["porechop_abi"]["nanopore"], "-v"),
    log:
        porechop_abi=config["out_dir"] / "log/porechop_abi.{nanopore_datatype}.trimmed.{fileprefix}.porechop_abi.log",
        cluster_log=config["out_dir"] / "log/porechop_abi.{nanopore_datatype}.trimmed.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/porechop_abi.{nanopore_datatype}.trimmed.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/porechop_abi.{nanopore_datatype}.trimmed.{fileprefix}.benchmark.txt"
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
        "     -i {input.fastq} -o {output.trimmed_fastq} 1>{log.porechop_abi} 2>&1; "

use rule porechop_abi as porechop_abi_track_data with:
    input:
        fastq=config["out_dir"] / ("ext_data/{nanopore_datatype}/{track_name}/raw/{fileprefix}%s" % config["fastq_ext"])
    output:
        trimmed_fastq=config["out_dir"] / ("ext_data/{nanopore_datatype}/{track_name}/trimmed/{fileprefix}%s" % config["fastq_ext"]),
    log:
        porechop_abi=config["out_dir"] / "log/porechop_abi_track_data.{nanopore_datatype}.{track_name}.trimmed.{fileprefix}.porechop_abi.log",
        cluster_log=config["out_dir"] / "log/porechop_abi_track_data.{nanopore_datatype}.{track_name}.trimmed.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/porechop_abi_track_data.{nanopore_datatype}.{track_name}.trimmed.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/porechop_abi_track_data.{nanopore_datatype}.{track_name}.trimmed.{fileprefix}.benchmark.txt"

rule chopper:
    input:
        input_fastq=config["out_dir"] / ("data/{nanopore_datatype}/%s/{fileprefix}%s" % ("trimmed" if not config["skip_porechop_abi"] else "raw",
                                                                                   config["fastq_ext"]))
    output:
        filtered_fastq=config["out_dir"] / ("data/{nanopore_datatype}/filtered/{fileprefix}%s" % config["fastq_ext"]),
    params:
         headcrop  = parse_option("headcrop",  parameters["tool_options"]["chopper"]["nanopore"], "--headcrop"),
         maxlength = parse_option("maxlength", parameters["tool_options"]["chopper"]["nanopore"], "--maxlength"),
         minlength = parse_option("minlength", parameters["tool_options"]["chopper"]["nanopore"], "--minlength"),
         quality   = parse_option("quality",   parameters["tool_options"]["chopper"]["nanopore"], "--quality"),
         tailcrop  = parse_option("tailcrop",  parameters["tool_options"]["chopper"]["nanopore"], "--tailcrop"),
    log:
        zcat=config["out_dir"] / "log/chopper.{nanopore_datatype}.filtered.{fileprefix}.zcat.log",
        chopper=config["out_dir"] / "log/chopper.{nanopore_datatype}.filtered.{fileprefix}.chopper.log",
        pigz_filtered=config["out_dir"] / "log/chopper.{nanopore_datatype}.filtered.{fileprefix}.pigz_filtered.log",
        cluster_log=config["out_dir"] / "log/chopper.{nanopore_datatype}.filtered.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/chopper.{nanopore_datatype}.filtered.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/chopper.{nanopore_datatype}.filtered.{fileprefix}.benchmark.txt"
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
        "     chopper {params.headcrop} {params.maxlength} {params.minlength} {params.quality} {params.tailcrop} "
        "         -t {threads} 2>{log.chopper} | pigz -p 4 > {output.filtered_fastq} 2>{log.pigz_filtered} ; "


use rule chopper as chopper_track_data with:
    input:
        input_fastq=config["out_dir"] / ("ext_data/{nanopore_datatype}/{track_name}/%s/{fileprefix}%s" % ("trimmed" if not config["skip_porechop_abi"] else "raw",
                                                                                   config["fastq_ext"]))
    output:
        filtered_fastq=config["out_dir"] / ("ext_data/{nanopore_datatype}/{track_name}/filtered/{fileprefix}%s" % config["fastq_ext"]),
    log:
        zcat=config["out_dir"] / "log/chopper_track_data.{nanopore_datatype}.{track_name}.filtered.{fileprefix}.zcat.log",
        chopper=config["out_dir"] / "log/chopper_track_data.{nanopore_datatype}.{track_name}.filtered.{fileprefix}.chopper.log",
        pigz_filtered=config["out_dir"] / "log/chopper_track_data.{nanopore_datatype}.{track_name}.filtered.{fileprefix}.pigz_filtered.log",
        cluster_log=config["out_dir"] / "log/chopper_track_data.{nanopore_datatype}.{track_name}.filtered.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/chopper_track_data.{nanopore_datatype}.{track_name}.filtered.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/chopper_track_data.{nanopore_datatype}.{track_name}.filtered.{fileprefix}.benchmark.txt"