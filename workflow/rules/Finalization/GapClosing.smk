
def get_read_files_for_samba(wildcards):
    phasing_kmer_length = stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"]["phasing_kmer_length"]
    print("AAAAAA")
    if phasing_kmer_length == "NA":
        print("BBBBBB")
        filelist = expand(output_dict["data"] / ("%s/%s/raw/{fileprefix}%s" % (datatype_format_dict[config["gap_closing_datatype"]],
                                                                              config["gap_closing_datatype"],
                                                                              config[datatype_format_dict[config["gap_closing_datatype"]] + "_extension"])),
                          allow_missing=True,
                          fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]] if datatype_format_dict[config["gap_closing_datatype"]] == "fastq" else input_fasta_file_prefix_dict[config["gap_closing_datatype"]])
    else:
        print("CCCCCCCCCC")
        filelist = expand(out_dir_path / ("%s/%s/%s/{haplotype}/%s/%s/{fileprefix}%s" % (config["phasing_stage"],
                                                                                          detect_phasing_parameters(wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters, config["phasing_stage"], stage_separator=".."),
                                                                                          datatype_format_dict[config["gap_closing_datatype"]] ,
                                                                                          stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"]["phasing_kmer_length"],
                                                                                          config["gap_closing_datatype"],
                                                                                          config[datatype_format_dict[config["gap_closing_datatype"]] + "_extension"])),
                         fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]] if datatype_format_dict[config["gap_closing_datatype"]] == "fastq" else input_fasta_file_prefix_dict[config["gap_closing_datatype"]],
                         allow_missing=True)
    print(filelist)

    return list(map(lambda s: Path(s).resolve(), filelist))


rule samba:
    priority: 500
    input:
        #reads=lambda wildcards: list(map(lambda s: s.resolve(), expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (config["gap_closing_datatype"],
        #                                                                                           config["fastq_extension"])),
        #                                fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]]))) if not stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"][config["gap_closing_datatype"]]["use_corrected_reads"] \
        #                                                    else (out_dir_path / ("data/fastq/%s/error_corrected_hifiasm_option_set_1/%s.contig.ec.fasta.gz" % (config["gap_closing_datatype"], wildcards.genome_prefix))).resolve(),
        reads=get_read_files_for_samba,
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters, wildcards.genome_prefix, wildcards.haplotype)


    output: 
        fasta=out_dir_path / "gap_closing/{prev_stage_parameters}..samba_{gap_closing_parameters}/{genome_prefix}.gap_closing.{haplotype}.fasta",
    params:
        datatype=lambda wildcards: parse_option("datatype", parameters["tool_options"]["samba"][wildcards.gap_closing_parameters][config["gap_closing_datatype"]], " -d "),
        matching_len=lambda wildcards: parse_option("matching_len", parameters["tool_options"]["samba"][wildcards.gap_closing_parameters][config["gap_closing_datatype"]], " -m ")
    log:
        samba=output_dict["log"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.samba.log",
        ln=output_dict["log"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["masurca"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["masurca"]["yaml"])
    resources:
        cpus=parameters["threads"]["samba"],
        time=parameters["time"]["samba"],
        mem=parameters["memory_mb"]["samba"],
    threads:
        parameters["threads"]["samba"]
    shell:
         " OUTPUT_DIR=`dirname {output.fasta}`/{wildcards.haplotype}; "
         " mkdir -p ${{OUTPUT_DIR}}; "
         " INPUT_FASTA=`realpath -s {input.fasta}`; "
         " INPUT_FASTA_BASENAME=`basename {input.fasta}`; "
         " LOG_SAMBA=`realpath -s {log.samba}`; "
         " LOG_LN=`realpath -s {log.ln}`; "
         " cd ${{OUTPUT_DIR}}; "
         " close_scaffold_gaps.sh -t {threads} -q <(zcat {input.reads}) {params.datatype} -r ${{INPUT_FASTA}} "
         " {params.matching_len} -v > ${{LOG_SAMBA}} 2>&1; "
         " ln -sf {wildcards.haplotype}/${{INPUT_FASTA_BASENAME}}.split.joined.fa ../`basename {output.fasta}` > ${{LOG_LN}} 2>&1; "
