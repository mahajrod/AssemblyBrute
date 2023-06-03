
"""
rule hifiasm_correct:
    priority: 2000
    input:
        hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                    fileprefix=input_file_prefix_dict["hifi"],
                    allow_missing=True),
    output:
        ec_bin=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.ec.bin",
        ovlp_reverse_bin=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.ovlp.reverse.bin",
        #primary_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.p_ctg.gfa" % config["genome_name"]),
        #alternative_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.a_ctg.gfa" % config["genome_name"]),
    params:
        dir=lambda wildcards: output_dict["contig"] / "hifiasm_{0}/".format(wildcards.contig_options),
        ec_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.ec.bin".format(wildcards.genome_prefix),
        ovlp_reverse_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.ovlp.reverse.bin".format(wildcards.genome_prefix),
        ovlp_source_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.ovlp.source.bin".format(wildcards.genome_prefix),

        output_prefix=lambda wildcards: output_dict["contig"] / "hifiasm_{0}/{1}.contig.hifi".format(wildcards.contig_options,
                                                                                                     wildcards.genome_prefix),
        window_size=lambda wildcards: parse_option("window_size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -k "),
        D=lambda wildcards: parse_option("D", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -N "), #" -N {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["N"]) if "N" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
    log:
        std=output_dict["log"] / "hifiasm_correct.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm_correct{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm_correct{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm_correct.{contig_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
        hifiasm=1
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " hifiasm -t {threads} -e --write-ec {params.window_size} {params.bloom_filter_bits} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         " {params.kmer_length} {params.D} {params.N} "
         " -o {params.output_prefix} {input.hifi}  1>{log.std} 2>&1;"

"""
rule hifiasm: # TODO: implement modes without hic data
    priority: 1000
    input:
        hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                    fileprefix=input_file_prefix_dict["hifi"],
                    allow_missing=True),
        hic_forward=input_filedict["hic"][::2] if "hic" in input_filedict else [],
        hic_reverse=input_filedict["hic"][1::2] if "hic" in input_filedict else [],
        genomescope_report=output_dict["kmer"] / ("%s/filtered/genomescope/{genome_prefix}.%s.filtered.%s.%s.genomescope.parameters" % (config["final_kmer_datatype"],
                                                                                                                                        config["final_kmer_datatype"],
                                                                                                                                        config["final_kmer_length"],
                                                                                                                                        config["final_kmer_counter"])),
    output:
        ec_bin=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.ec.bin",
        ovlp_reverse_bin=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.ovlp.reverse.bin",
        primary_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.hic.hap1.p_ctg.gfa",
        alternative_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.hic.hap2.p_ctg.gfa",
        alt_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hifi.hic.a_ctg.gfa",
        primary_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hap1.gfa",
        alternative_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hap2.gfa",
        alt_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.alt.gfa",
        #primary_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.p_ctg.gfa" % config["genome_name"]),
        #alternative_contig_graph=output_dict["contig"] / ("hifiasm/%s.contig.hifiasm.hifi.hic.a_ctg.gfa" % config["genome_name"]),
    params:
        dir=lambda wildcards: output_dict["contig"] / "hifiasm_{0}/".format(wildcards.contig_options),
        ec_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.ec.bin".format(wildcards.genome_prefix),
        ovlp_reverse_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.ovlp.reverse.bin".format(wildcards.genome_prefix),
        ovlp_source_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.ovlp.source.bin".format(wildcards.genome_prefix),
        hic_lk_bin= lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.hic.lk.bin".format(wildcards.genome_prefix),
        hic_tlb_bin=lambda wildcards: output_dict["contig"] / "hifiasm_*/{0}.contig.hifi.hic.tlb.bin".format(wildcards.genome_prefix),
        purge_level=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["purge level"],
        ploidy=config["ploidy"],
        output_prefix=lambda wildcards: output_dict["contig"] / "hifiasm_{0}/{1}.contig.hifi".format(wildcards.contig_options,
                                                                                                     wildcards.genome_prefix),
        cov_multiplicator=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["cov_multiplicator"],

        #primary=lambda wildcards: parse_option_flag("primary", parameters["tool_options"]["hifiasm"][wildcards.contig_options],  "--primary"), # "--primary" if parameters["tool_options"]["hifiasm"][wildcards.contig_options]["primary"] else "",
        window_size=lambda wildcards: parse_option("window_size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -k "),
        D=lambda wildcards: parse_option("D", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -N "),
        ignore_bin=lambda wildcards: " -i " if ("ignore_bin" in parameters["tool_options"]["hifiasm"][wildcards.contig_options]) and parameters["tool_options"]["hifiasm"][wildcards.contig_options]["ignore_bin"] else "",
        hic_forward=" --h1 " + ",".join(map(str, input_filedict["hic"][::2]) if "hic" in input_filedict else []), #in case of multiple hic libraries files in the list MUST be COMMA-separated
        hic_reverse=" --h2 " + ",".join(map(str, input_filedict["hic"][1::2]) if "hic" in input_filedict else [])
    log:
        std=output_dict["log"] / "hifiasm.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm.{contig_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
        hifiasm=1
    threads:
        parameters["threads"]["hifiasm"]
    shell: # TODO: rewrite test for presence of previous run as a python script and add additional tests and detection of previous run type (hifi + hic or hifi only)
         " # check if there was a hifiasm run\n"
         " if [ '{params.ignore_bin}' = ' -i ' ]; then"
         "      echo 'Option -i is set. Ignoring previous hifiasm runs...'; "
         " else"
         "      echo 'Option -i is not set. Seeking for previous hifiasm runs...'; "
         "      EC_BIN=(`find ./ -wholename \"*{params.ec_bin}\"`); "
         "      OVLP_REVERSE_BIN=(`find ./ -wholename \"*{params.ovlp_reverse_bin}\"`); "
         "      OVLP_SOURCE_BIN=(`find ./ -wholename \"*{params.ovlp_source_bin}\"`); "
         "      HIC_LK_BIN=(`find ./ -wholename \"*{params.hic_lk_bin}\"`); "
         "      HIC_TLB_BIN=(`find ./ -wholename \"*{params.hic_tlb_bin}\"`); "
         "      if [ ${{#EC_BIN[@]}} -eq 0 ]; then "
         "              echo 'First hifiasm run!'; "
         "      else "
         "          echo 'Previous hifiasm run detected! *.bin files from previous run will be used...'; "
         "          cp ${{EC_BIN[0]}} ${{OVLP_REVERSE_BIN[0]}} ${{OVLP_SOURCE_BIN[0]}} ${{HIC_LK_BIN[0]}} ${{HIC_TLB_BIN[0]}} {params.dir};  "
         "      fi ; "
         " fi;"
         " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "
         " hifiasm {params.D} {params.N} {params.ignore_bin} --primary -t {threads} -l {params.purge_level}  -o {params.output_prefix} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         " {params.hic_forward} {params.hic_reverse} "
         " {input.hifi}  1>{log.std} 2>&1;"
         " ln {output.primary_contig_graph} {output.primary_alias};"
         " ln {output.alternative_contig_graph} {output.alternative_alias};"
         " ln {output.alt_contig_graph} {output.alt_alias}; "
