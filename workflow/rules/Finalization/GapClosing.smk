
rule samba: #TODO: FIX CASE OF ERROR CORRECTED READS,NOW IT IS HARDCODED TO error_corrected_hifiasm_option_set_1
    priority: 500
    input:
        reads=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (config["gap_closing_datatype"],
                                                                                                   config["fastq_extension"])),
                                        fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]],
                                        allow_missing=True) if not stage_dict["gap_closing"]["parameters"][wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"][config["gap_closing_datatype"]]["use_corrected_reads"] \
                                                            else out_dir_path / ("data/fastq/%s/error_corrected_hifiasm_option_set_1/%s.contig.ec.fasta.gz" % (config["gap_closing_datatype"], wildcards.genome_prefix)),
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
         " cd ${{OUTPUT_DIR}}; "
         " close_scaffold_gaps.sh -t {threads} -q <(zcat {input.reads}) {params.datatype} -r ${{INPUT_FASTA}} "
         " {params.matching_len} -v > {log.samba} 2>&1; "
         " ln -sf {wildcards.haplotype}/${{INPUT_FASTA_BASENAME}}.split.joined.fa {output.fasta} > {log.ln} 2>&1; "
