localrules: parse_genomescope_output, parse_genomescope_output_per_lib, parse_genomescope_output_ploidy_test

def get_starting_lambda(wildcards):
    if "genomescope" in config["tool_manually_adjusted_features"]:
        if "starting_lambda" in config["tool_manually_adjusted_features"]["genomescope"]:
            if config["tool_manually_adjusted_features"]["genomescope"]["starting_lambda"] is not None:
                return " --lambda {0}".format(config["tool_manually_adjusted_features"]["genomescope"]["starting_lambda"])
            else:
                return ""
        else:
            return ""
    else:
        return ""

def get_start_shift(wildcards):
    if "genomescope" in config["tool_manually_adjusted_features"]:
        if "start_shift" in config["tool_manually_adjusted_features"]["genomescope"]:
            if config["tool_manually_adjusted_features"]["genomescope"]["start_shift"] is not None:
                return " --start_shift {0} ".format(config["tool_manually_adjusted_features"]["genomescope"]["start_shift"])
            else:
                return " --start_shift 2 "
        else:
            return " --start_shift 2 "
    else:
        return " --start_shift 2 "

rule genomescope:
    input:
        histo="{db_dir}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo",
        log="{db_dir}/log/"
    output:
        summary="{db_dir}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/%s_summary.txt" % config["genome_prefix"],
        model="{db_dir}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}/%s_model.txt" % config["genome_prefix"],
    params:
        ploidy=config["ploidy"],
        genome_name=config["genome_prefix"],
        #max_coverage=lambda wildcards: parameters["tool_options"][wildcards.kmer_tool][wildcards.datatype]["max_coverage"],
        start_shift=get_start_shift,
        starting_lambda=get_starting_lambda
    log:
        std="{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.log",
        cluster_log="{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.log",
        cluster_err="{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.err"
    benchmark:
        "{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("genomescope"),
        cpus=parameters["threads"]["genomescope"],
        time=parameters["time"]["genomescope"],
        mem=parameters["memory_mb"]["genomescope"],
    threads:
        parameters["threads"]["genomescope"]
    shell:
         " OUT_DIR=`dirname {output.summary}`; "
         " genomescope2 {params.start_shift} {params.starting_lambda} -i {input.histo} -p {params.ploidy} -k {wildcards.kmer_length}  "
         "     -n {params.genome_name} --fitted_hist  --testing  -o ${{OUT_DIR}} > {log.std} 2>&1" # -m {params.max_coverage}

use rule genomescope as genomescope_ploidy_test with:
    output:
        summary="{db_dir}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}@p{ploidy}/%s_summary.txt" % config["genome_prefix"],
        model="{db_dir}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}@p{ploidy}/%s_model.txt" % config["genome_prefix"],
    params:
        ploidy=lambda wildcards: wildcards.ploidy,
        genome_name=config["genome_prefix"],
        start_shift=get_start_shift,
        starting_lambda=get_starting_lambda
    log:
        std="{db_dir}/log/genomescope_ploidy_test.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.log",
        cluster_log="{db_dir}/log/genomescope_ploidy_test.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.cluster.log",
        cluster_err="{db_dir}/log/genomescope_ploidy_test.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.cluster.err"
    benchmark:
        "{db_dir}/log/genomescope_ploidy_test.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.benchmark.txt"

use rule genomescope as genomescope_per_lib with:
    input:
        histo="{db_dir}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.histo",
        log="{db_dir}/log/"
    output:
        summary="{db_dir}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}/%s_summary.txt" % config["genome_prefix"],
        model="{db_dir}/genomescope/{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}/%s_model.txt" % config["genome_prefix"],
    log:
        std="{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.log",
        cluster_log="{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.cluster.log",
        cluster_err="{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.cluster.err"
    benchmark:
        "{db_dir}/log/genomescope.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.benchmark.txt"


rule parse_genomescope_output:
    input:
        summary=rules.genomescope.output.summary,
        model=rules.genomescope.output.model,
        log="{db_dir}/log/"
    output:
        "{db_dir}/genomescope/%s.{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters"  % config["genome_prefix"]
    log:
        std="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.log",
        cluster_log="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.log",
        cluster_err="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.err"
    benchmark:
        "{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
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

use rule parse_genomescope_output as parse_genomescope_output_per_lib with:
    input:
        summary=rules.genomescope_per_lib.output.summary,
        model=rules.genomescope_per_lib.output.model,
        log="{db_dir}/log/"
    output:
        "{db_dir}/genomescope/%s.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.genomescope.parameters"  % config["genome_prefix"]
    log:
        std="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.log",
        cluster_log="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.cluster.log",
        cluster_err="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.cluster.err"
    benchmark:
        "{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.benchmark.txt"

use rule parse_genomescope_output as parse_genomescope_output_ploidy_test with:
    input:
        summary=rules.genomescope_ploidy_test.output.summary,
        model=rules.genomescope_ploidy_test.output.model,
        log="{db_dir}/log/"
    output:
        "{db_dir}/genomescope/%s.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.genomescope.parameters"  % config["genome_prefix"]
    log:
        std="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.log",
        cluster_log="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.cluster.log",
        cluster_err="{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.cluster.err"
    benchmark:
        "{db_dir}/log/parse_genomescope_output.{datatype}.{stage}.{kmer_length}.{kmer_tool}.p{ploidy}.benchmark.txt"