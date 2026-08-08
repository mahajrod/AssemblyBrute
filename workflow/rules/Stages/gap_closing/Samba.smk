
ruleorder: reorder_samba_output > maskfasta

def get_read_files_for_samba(wildcards):
    phasing_kmer_length = stage_dict["gap_closing"].parameters[wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"]["phasing_kmer_length"]
    datatype = config["gap_closing_datatype"]
    if phasing_kmer_length == "NA":
        filelist = expand(config["out_dir"] / ("data/%s/final/{fileprefix}%s" % (config["gap_closing_datatype"],
                                                                                           config["data"][datatype]["conv_ext"])),
                          allow_missing=True,
                          fileprefix=config["data"][datatype]["conv_file_prefix_list"])
    else:
        filelist = expand(config["out_dir"] / ("%s/%s/%s.%s.{haplotype}/reads/%s/%s/{fileprefix}%s" % (config["phasing_stage"],
                                                                                                             detect_phasing_parameters(wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters, config["phasing_stage"], stage_separator=".."),
                                                                                                             config["genome_prefix"],
                                                                                                             config["phasing_stage"],
                                                                                                             datatype,
                                                                                                             stage_dict["gap_closing"].parameters[wildcards.prev_stage_parameters + "..samba_" + wildcards.gap_closing_parameters]["option_set"]["phasing_kmer_length"],
                                                                                                             config["data"][datatype]["conv_ext"])),
                         fileprefix=config["data"][datatype]["conv_file_prefix_list"],
                         allow_missing=True)

    return filelist


rule samba:
    priority: 500
    input:
        reads=get_read_files_for_samba,
        fasta=lambda wildcards: config["out_dir"] / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["gap_closing"].prev_stage,
                                                                                       wildcards.prev_stage_parameters, wildcards.genome_prefix, wildcards.haplotype)
    output:
        fasta=config["out_dir"] / "gap_closing/{prev_stage_parameters}..samba_{gap_closing_parameters}/{genome_prefix}.gap_closing.{haplotype, hap[^/]*}/samba/{genome_prefix}.gap_closing.{haplotype}.fasta" ,
    params:
        datatype=lambda wildcards: parse_option("datatype", parameters["tool_options"]["samba"][wildcards.gap_closing_parameters][config["gap_closing_datatype"]], " -d "),
        matching_len=lambda wildcards: parse_option("matching_len", parameters["tool_options"]["samba"][wildcards.gap_closing_parameters][config["gap_closing_datatype"]], " -m ")
    log:
        samba=(config["out_dir"] / "log/samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.samba.log").resolve(),
        cluster_log=config["out_dir"] / "log/samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/samba.gap_closing.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
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
         " OUTPUT_DIR=`dirname {output.fasta}`/; "
         " mkdir -p ${{OUTPUT_DIR}}; "
         " INPUT_FASTA=`realpath -s {input.fasta}`; "
         " ln -sf ${{INPUT_FASTA}} ${{OUTPUT_DIR}}; "
         " INPUT_FILES=''; "
         " for FILE in {input.reads}; "
         "     do "
         "     ln -sf `realpath -s ${{FILE}}` ${{OUTPUT_DIR}}; "
         "     INPUT_FILES=\"${{INPUT_FILES}} \"`basename ${{FILE}}`; "
         "     done; "
         " cd ${{OUTPUT_DIR}}; "
         " close_scaffold_gaps.sh -t {threads} -q <(zcat ${{INPUT_FILES}}) {params.datatype} -r `basename ${{INPUT_FASTA}}` "
         "     {params.matching_len} -v > {log.samba} 2>&1; "
         " ln -sf `basename {input.fasta}`.split.joined.fa `basename {output.fasta}`"

rule reorder_samba_output:
    priority: 500
    input:
        fasta=rules.samba.output.fasta
    output:
        fasta=config["out_dir"] / "gap_closing/{prev_stage_parameters}..samba_{gap_closing_parameters}/{genome_prefix}.gap_closing.{haplotype, hap.*}.fasta" ,
    log:
        reorder=config["out_dir"] / "log/reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.reorder.log",
        cluster_log=config["out_dir"] / "log/reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/reorder_samba_output.{prev_stage_parameters}..samba_{gap_closing_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
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
         "     -o {output.fasta} --by orderlist  > {log.reorder} 2>&1; "