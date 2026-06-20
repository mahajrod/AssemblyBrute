
def get_read_files_for_samba(wildcards):
    phasing_kmer_length = stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"]["phasing_kmer_length"]
    #print("AAAAAA")
    if phasing_kmer_length == "NA":
        #print("BBBBBB")
        filelist = expand(output_dict["data"] / ("%s/%s/raw/{fileprefix}%s" % (datatype_format_dict[config["gap_closing_datatype"]],
                                                                              config["gap_closing_datatype"],
                                                                              config[datatype_format_dict[config["gap_closing_datatype"]] + "_extension"])),
                          allow_missing=True,
                          fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]] if datatype_format_dict[config["gap_closing_datatype"]] == "fastq" else input_fasta_file_prefix_dict[config["gap_closing_datatype"]])
    else:
        #print("CCCCCCCCCC")
        filelist = expand(out_dir_path / ("%s/%s/%s/{haplotype}/%s/%s/{fileprefix}%s" % (config["phasing_stage"],
                                                                                          detect_phasing_parameters(wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters, config["phasing_stage"], stage_separator=".."),
                                                                                          datatype_format_dict[config["gap_closing_datatype"]] ,
                                                                                          stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"]["phasing_kmer_length"],
                                                                                          config["gap_closing_datatype"],
                                                                                          config[datatype_format_dict[config["gap_closing_datatype"]] + "_extension"])),
                         fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]] if datatype_format_dict[config["gap_closing_datatype"]] == "fastq" else input_fasta_file_prefix_dict[config["gap_closing_datatype"]],
                         allow_missing=True)
    #print(filelist)

    return filelist


rule samba:
    priority: 500
    input:
        reads=get_read_files_for_samba,
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters, wildcards.genome_prefix, wildcards.haplotype)
    output:
        fasta=out_dir_path / "gap_closing/{prev_stage_parameters}..samba_{gap_closing_parameters}/{genome_prefix}.gap_closing.{haplotype, hap.*}/{genome_prefix}.gap_closing.{haplotype}.split.joined.fa" ,
    params:
        datatype=lambda wildcards: parse_option("datatype", parameters["tool_options"]["samba"][wildcards.gap_closing_parameters][config["gap_closing_datatype"]], " -d "),
        matching_len=lambda wildcards: parse_option("matching_len", parameters["tool_options"]["samba"][wildcards.gap_closing_parameters][config["gap_closing_datatype"]], " -m ")
    log:
        samba=output_dict["log"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.samba.log",
        reorder=output_dict["log"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.reorder.log",
        cluster_log=output_dict["cluster_log"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["masurca"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["masurca"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("samba"),
        cpus=parameters["threads"]["samba"],
        time=parameters["time"]["samba"],
        mem=parameters["memory_mb"]["samba"],
    threads:
        parameters["threads"]["samba"]
    shell:
         " OUTPUT_DIR=`dirname {output.fasta}`/{wildcards.genome_prefix}.gap_closing.{wildcards.haplotype}; "
         " mkdir -p ${{OUTPUT_DIR}}; "
         " INPUT_FASTA=`realpath -s {input.fasta}`; "
         " INPUT_FASTA_BASENAME=`basename {input.fasta}`; "
         " LOG_SAMBA=`realpath -s {log.samba}`; "
         " LOG_REORDER=`realpath -s {log.reorder}`; "
         " INPUT_FILES='';"
         " for FILE in {input.reads}; do INPUT_FILES=\"${{INPUT_FILES}} \"`realpath -s ${{FILE}}`; done; "
         " cd ${{OUTPUT_DIR}}; "
         " close_scaffold_gaps.sh -t {threads} -q <(zcat ${{INPUT_FILES}}) {params.datatype} -r ${{INPUT_FASTA}} "
         " {params.matching_len} -v > ${{LOG_SAMBA}} 2>&1; "

rule reorder_samba_output:
    priority: 500
    input:
        fasta=out_dir_path / "gap_closing/{prev_stage_parameters}..samba_{gap_closing_parameters}/{genome_prefix}.gap_closing.{haplotype, hap.*}/{genome_prefix}.gap_closing.{haplotype}.split.joined.fa" ,
    output:
        fasta=out_dir_path / "gap_closing/{prev_stage_parameters}..samba_{gap_closing_parameters}/{genome_prefix}.gap_closing.{haplotype, hap.*}.fasta" ,
    log:
        reorder=output_dict["log"] / "reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.reorder.log",
        cluster_log=output_dict["cluster_log"] / "reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("reorder_samba_output"),
        cpus=parameters["threads"]["reorder_samba_output"],
        time=parameters["time"]["reorder_samba_output"],
        mem=parameters["memory_mb"]["reorder_samba_output"],
    threads:
        parameters["threads"]["reorder_samba_output"]
    shell:
         " INPUT_PREFIX={input.fasta}; "
         " INPUT_PREFIX=${{INPUT_PREFIX%.fa}}; "
         " grep -P '^>' {input.fasta} | sed 's/>//;s/[ \t].*//' | sort -V > ${{INPUT_PREFIX}}.sorted.ids;  "  
         " ./workflow/scripts/sequence/reorder_sequences.py -i  {input.fasta} -r ${{INPUT_PREFIX}}.sorted.ids "
         "                                                  -o {output.fasta}  > {log.reorder} 2>&1; "