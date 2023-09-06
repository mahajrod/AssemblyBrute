localrules: parse_genomescope_output

rule genomescope:
    input:
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo"
    output:
        summary=output_dict["kmer"] / "{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}_summary.txt",
        model=output_dict["kmer"] / "{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}_model.txt",
    params:
        ploidy=config["ploidy"],
        genome_name=lambda wildcards: wildcards.genome_prefix,
        #max_coverage=lambda wildcards: parameters["tool_options"][wildcards.kmer_tool][wildcards.datatype]["max_coverage"],
        out_dir=lambda wildcards: output_dict["kmer"] / "{0}/{1}/genomescope/{0}.{1}.{2}.{3}".format(wildcards.datatype,
                                                                                                     wildcards.stage,
                                                                                                     wildcards.kmer_length,
                                                                                                     wildcards.kmer_tool)
    log:
        std=output_dict["log"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("genomescope"),
        cpus=parameters["threads"]["genomescope"],
        time=parameters["time"]["genomescope"],
        mem=parameters["memory_mb"]["genomescope"],
    threads:
        parameters["threads"]["genomescope"]
    shell:
         " genomescope.R --start_shift 2 -i {input.histo} -p {params.ploidy} -k {wildcards.kmer_length}  "
         " -n {params.genome_name} --fitted_hist  --testing  -o {params.out_dir} > {log.std} 2>&1" # -m {params.max_coverage}


rule parse_genomescope_output:
    input:
        summary=rules.genomescope.output.summary,
        model=rules.genomescope.output.model
        #summary=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/%s_summary.txt" % config["genome_name"]),
        #model=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/%s_model.txt" % config["genome_name"])
    output:
        output_dict["kmer"] / "{datatype}/{stage}/genomescope/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters"
    log:
        std=output_dict["log"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("parse_genomescope_output"),
        cpus=parameters["threads"]["parse_genomescope_output"],
        time=parameters["time"]["parse_genomescope_output"],
        mem=parameters["memory_mb"]["parse_genomescope_output"],
    threads:
        parameters["threads"]["parse_genomescope_output"]
    shell:
         " GENLEN=`grep 'Genome Haploid Length' {input.summary} | sed 's/,//g;s/ \{{2,\}}/\t/g' | cut -f 3 | sed 's/ .*//'`;   "
         " LAMBDA=`grep 'kmercov' {input.model} | tail -n  1 | awk '{{printf \"%.0f\", $2}}'`;"
         " echo -e \"Genome size\\t${{GENLEN}}\\nLambda\\t${{LAMBDA}}\\n\" > {output} 2>{log.std}"
