

rule trimmomatic_pe:
    input:
        forward_fastq=lambda wildcards: config["out_dir"] / ("data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                 "raw" if wildcards.pe_datatype != "hic" else "trimmed",
                                                                                                  wildcards.pairprefix,
                                                                                                  config["fwd_fastq_sfx"],
                                                                                                  config["fastq_ext"])),
        reverse_fastq=lambda wildcards: config["out_dir"] / ("data/{0}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                 "raw" if wildcards.pe_datatype != "hic" else "trimmed",
                                                                                                 wildcards.pairprefix,
                                                                                                 config["rev_fastq_sfx"],
                                                                                                 config["fastq_ext"])),
    output:
        forward_fastq=config["out_dir"] / ("data/{pe_datatype}/filtered/{pairprefix}_1%s" % config["fastq_ext"]),
        forward_se_fastq=config["out_dir"] / ("data/{pe_datatype}/filtered/{pairprefix}_1.se%s" % config["fastq_ext"]),
        reverse_fastq=config["out_dir"] / ("data/{pe_datatype}/filtered/{pairprefix}_2%s" % config["fastq_ext"]),
        reverse_se_fastq=config["out_dir"] / ("data/{pe_datatype}/filtered/{pairprefix}_2.se%s" % config["fastq_ext"]),
        stats=config["out_dir"] / "data/{pe_datatype}/filtered/{pairprefix}.trimmomatic.stats"
    params:
        min_read_length=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.pe_datatype]["min_read_length"],
        sliding_window_size=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.pe_datatype]["sliding_window_size"],
        sliding_window_quality=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.pe_datatype]["sliding_window_quality"],
        illumina_clip=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.pe_datatype]["illumina_clip"],
        adapter_file=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.pe_datatype]["adapter_file"],
        mem=int(0.5 * parameters["memory_mb"]["trimmomatic_pe"]) # trimmomatic (at least version 0.40) eats more memory than it is allowed by -Xmx option
    log:
        std=config["out_dir"] / "log/trimmomatic_pe.{pe_datatype}.{pairprefix}.log",
        cluster_log=config["out_dir"] / "log/trimmomatic_pe.{pe_datatype}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/trimmomatic_pe.{pe_datatype}.{pairprefix}.cluster.log"
    benchmark:
        config["out_dir"] / "log/trimmomatic_pe.{pe_datatype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("trimmomatic_pe"),
        cpus=parameters["threads"]["trimmomatic_pe"],
        time=parameters["time"]["trimmomatic_pe"],
        mem=parameters["memory_mb"]["trimmomatic_pe"],
    threads:
        parameters["threads"]["trimmomatic_pe"]
    shell:
         "trimmomatic -Xmx{params.mem}m  PE -threads {threads} -phred33 {input.forward_fastq} {input.reverse_fastq} "
         "     {output.forward_fastq} {output.forward_se_fastq} {output.reverse_fastq} {output.reverse_se_fastq} "
         "     ILLUMINACLIP:{params.adapter_file}:{params.illumina_clip} "
         "     SLIDINGWINDOW:{params.sliding_window_size}:{params.sliding_window_quality} "
         "     MINLEN:{params.min_read_length} > {output.stats} 2>{log.std}; "

use rule trimmomatic_pe as trimmomatic_pe_track_data with:
    input:
        forward_fastq=lambda wildcards: config["out_dir"] / ("ext_track_data/{0}/{5}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                       "raw" if wildcards.pe_datatype != "hic" else "trimmed",
                                                                                                        wildcards.pairprefix,
                                                                                                        config["fwd_fastq_sfx"],
                                                                                                        config["fastq_ext"],
                                                                                                        wildcards.track_name)),
        reverse_fastq=lambda wildcards: config["out_dir"] / ("ext_track_data/{0}/{5}/{1}/{2}{3}{4}".format(wildcards.pe_datatype,
                                                                                                       "raw" if wildcards.pe_datatype != "hic" else "trimmed",
                                                                                                       wildcards.pairprefix,
                                                                                                       config["rev_fastq_sfx"],
                                                                                                       config["fastq_ext"],
                                                                                                       wildcards.track_name)),
    output:
        forward_fastq=config["out_dir"] / ("ext_track_data/{pe_datatype}/{track_name}/filtered/{pairprefix}_1%s" % config["fastq_ext"]),
        forward_se_fastq=config["out_dir"] / ("ext_track_data/{pe_datatype}/{track_name}/filtered/{pairprefix}_1.se%s" % config["fastq_ext"]),
        reverse_fastq=config["out_dir"] / ("ext_track_data/{pe_datatype}/{track_name}/filtered/{pairprefix}_2%s" % config["fastq_ext"]),
        reverse_se_fastq=config["out_dir"] / ("ext_track_data/{pe_datatype}/{track_name}/filtered/{pairprefix}_2.se%s" % config["fastq_ext"]),
        stats=config["out_dir"] / "ext_track_data/{pe_datatype}/{track_name}/filtered/{pairprefix}.trimmomatic.stats"
    log:
        std=config["out_dir"] / "log/trimmomatic_pe_track_data.{pe_datatype}.{track_name}.{pairprefix}.log",
        cluster_log=config["out_dir"] / "log/trimmomatic_pe_track_data.{pe_datatype}.{track_name}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/trimmomatic_pe_track_data.{pe_datatype}.{track_name}.{pairprefix}.cluster.log"
    benchmark:
        config["out_dir"] / "log/trimmomatic_pe_track_data.{pe_datatype}.{track_name}.{pairprefix}.benchmark.txt"