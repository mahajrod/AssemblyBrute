
rule hifiasm_correct:
    priority: 2000
    input:
        hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                    fileprefix=input_file_prefix_dict["hifi"],
                    allow_missing=True),
    output:
        ec_bin=output_dict["error_correction"] / "hifiasm_{correction_options}/{genome_prefix}.contig.ec.bin",
        ec_fasta=output_dict["error_correction"] / "hifiasm_{correction_options}/{genome_prefix}.contig.ec.fasta.gz",
        alias_ec_fasta=out_dir_path / "data/fastq/hifi/error_corrected_hifiasm_{correction_options}/{genome_prefix}.contig.ec.fasta.gz",
        ovlp_reverse_bin=output_dict["error_correction"] / "hifiasm_{correction_options}/{genome_prefix}.contig.ovlp.reverse.bin",
        ovlp_source_bin=output_dict["error_correction"] / "hifiasm_{correction_options}/{genome_prefix}.contig.ovlp.source.bin",
    params:
        window_size=lambda wildcards: parse_option("window_size", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -k "),
        D=lambda wildcards: parse_option("D", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -N "), #" -N {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["N"]) if "N" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
    log:
        std=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.log",
        pigz=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.pigz.log",
        mv=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.mv.log",
        ln=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.ln.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm_correct{correction_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm_correct{correction_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm_correct.{correction_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.ec_bin}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.ec.bin}}; "
         " UNCOMPRESSED_FASTA={output.ec_fasta}; "
         " UNCOMPRESSED_FASTA=${{UNCOMPRESSED_FASTA%sta.gz}}; "
         " hifiasm -t {threads} -e --write-ec {params.window_size} {params.bloom_filter_bits} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         " {params.kmer_length} {params.D} {params.N} "
         " -o ${{OUTPUT_PREFIX}} {input.hifi}  1>{log.std} 2>&1;"
         " pigz -p {threads} ${{UNCOMPRESSED_FASTA}} > {log.pigz} 2>&1 ; "
         " mv ${{UNCOMPRESSED_FASTA}}.gz {output.ec_fasta} > {log.mv} 2>&1; "
         " ln {output.ec_fasta} {output.alias_ec_fasta} > {log.ln} 2>&1; "

rule hifiasm_hic: # TODO: implement modes without hic data as independent rule
    priority: 1000
    input:
        hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                    fileprefix=input_file_prefix_dict["hifi"],
                    allow_missing=True),
        hic_forward=input_filedict["hic"][::2] if "hic" in input_filedict else [],
        hic_reverse=input_filedict["hic"][1::2] if "hic" in input_filedict else [],
        nanopore=input_filedict["nanopore"] if "nanopore" in input_filedict else [], # TODO: test this option
        ec_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ec.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                          wildcards.genome_prefix),
        ovlp_reverse_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.reverse.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                              wildcards.genome_prefix),
        ovlp_source_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.source.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                            wildcards.genome_prefix),
        genomescope_report=output_dict["kmer"] / ("%s/filtered/genomescope/{genome_prefix}.%s.filtered.%s.%s.genomescope.parameters" % (config["final_kmer_datatype"],
                                                                                                                                        config["final_kmer_datatype"],
                                                                                                                                        config["final_kmer_length"],
                                                                                                                                        config["final_kmer_counter"])),
    output:
        primary_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        alternative_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        alt_contig_graph=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hic.a_ctg.gfa",
        primary_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hap1.gfa",
        alternative_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hap2.gfa",
        alt_alias=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.alt.gfa",
    params:
        purge_level=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["purge level"],
        ploidy=config["ploidy"],
        cov_multiplicator=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["cov_multiplicator"],
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
        hic_forward=(" --h1 " + ",".join(map(str, input_filedict["hic"][::2]))) if "hic" in input_filedict else "", #in case of multiple hic libraries files in the list MUST be COMMA-separated
        hic_reverse=(" --h2 " + ",".join(map(str, input_filedict["hic"][1::2]))) if "hic" in input_filedict else "",
        nanopore=(" --ul " + ",".join(map(str, input_filedict["nanopore"]))) if "nanopore" in input_filedict else "",
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
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.primary_alias}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.hap1.gfa}}; "
         " OUT_DIR=`dirname ${{OUTPUT_PREFIX}}`; "
         " ln -s ../../../{input.ec_bin} ${{OUT_DIR}}; "
         " ln -s ../../../{input.ovlp_reverse_bin} ${{OUT_DIR}}; "
         " ln -s ../../../{input.ovlp_source_bin} ${{OUT_DIR}}; "
         " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "
         " hifiasm {params.window_size} {params.bloom_filter_bits} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         " {params.kmer_length} {params.D} {params.N} {params.ignore_bin} --primary -t {threads} -l {params.purge_level}  -o ${{OUTPUT_PREFIX}} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         " {params.hic_forward} {params.hic_reverse} {params.nanopore} "
         " {input.hifi}  1>{log.std} 2>&1;"
         " ln {output.primary_contig_graph} {output.primary_alias};"
         " ln {output.alternative_contig_graph} {output.alternative_alias};"
         " ln {output.alt_contig_graph} {output.alt_alias}; "
