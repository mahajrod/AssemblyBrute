localrules: smudgeplot_assess

rule smudgeplot_assess:
    input:
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo"
    output:
        boundaries=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.smudgeplot.boundaries",
    log:
        upper=output_dict["log"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.upper.log",
        lower=output_dict["log"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.lower.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot_assess.{datatype}.{stage}.{kmer_length}.{kmer_tool}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("smudgeplot_assess"),
        cpus=parameters["threads"]["smudgeplot_plot"],
        time=parameters["time"]["smudgeplot_plot"],
        mem=parameters["memory_mb"]["smudgeplot_plot"],
    threads:
        parameters["threads"]["smudgeplot_plot"]
    shell:
         " LOWER_BOUNDARY=`smudgeplot.py cutoff {input.histo} L 2>{log.lower}`; "
         " UPPER_BOUNDARY=`smudgeplot.py cutoff {input.histo} U 2>{log.upper}`; "
         " echo -e \"low_boundary\tupper_boundary\n${{LOWER_BOUNDARY}}\t${{UPPER_BOUNDARY}}\n\" > {output.boundaries}"

rule smudgeplot_hetkmers:
    input:
        kmer=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer",
    output:
        coverages=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_coverages.tsv",
        sequences=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_sequences.tsv",
        #smudgeplot=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_smudgeplot.png",
        #summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv"
    log:
        hetkmers=output_dict["log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.hetkmers.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("smudgeplot_hetkmers"),
        cpus=parameters["threads"]["smudgeplot_hetkmers"],
        time=parameters["time"]["smudgeplot_hetkmers"],
        mem=parameters["memory_mb"]["smudgeplot_hetkmers"],
        smudgeplot_hetkmers=1
    threads:
        parameters["threads"]["smudgeplot_hetkmers"]
    shell:
         " COV_OUT={output.coverages}; "
         " PREFIX=${{COV_OUT%_coverages.tsv}}; "
         #" HAPLOID_COVERAGE=`awk 'NR==2 {{print 2 * $2}}' {input.genomescope_report}`; "
         " smudgeplot.py hetkmers -o ${{PREFIX}} {input.kmer} > {log.hetkmers} 2>&1; "
         #" smudgeplot.py plot -k {wildcards.kmer_length} -n ${{HAPLOID_COVERAGE}}  -o ${{PREFIX}} {output.coverages} > {log.plot} 2>&1; "

rule smudgeplot_plot: # in some cases smudgeplot could fail in geneeration of image. It's not crucial analysis, so better to continue pipeline execution. Because of it warning file was set as output
    input:
        coverages=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_coverages.tsv",
        genomescope_report=output_dict["kmer"] / ("{datatype}/{stage}/genomescope/%s.%s.filtered.%s.%s.genomescope.parameters" % (config["genome_prefix"],
                                                                                                                                  config["final_kmer_datatype"],
                                                                                                                                  config["final_kmer_length"],
                                                                                                                                  config["final_kmer_counter"])),
    output:
        warnings=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_warnings.txt",
        warnings_no_priors=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.no_priors_warnings.txt",
        #smudgeplot=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_smudgeplot.png",
        #summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv",
        #smudgeplot_no_priors=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.no_priors_smudgeplot.png",
        #summary_no_priors=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.no_priors_summary_table.tsv"
    log:
        plot=output_dict["log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.plot.log",
        cluster_log=output_dict["cluster_log"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "smudgeplot.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("smudgeplot_plot"),
        cpus=parameters["threads"]["smudgeplot_plot"],
        time=parameters["time"]["smudgeplot_plot"],
        mem=parameters["memory_mb"]["smudgeplot_plot"],
    threads:
        parameters["threads"]["smudgeplot_plot"]
    shell:
         " WARNINGS={output.warnings}; "
         " PREFIX=${{WARNINGS%_warnings.txt}}; "
         " HAPLOID_COVERAGE=`awk 'NR==2 {{print $2}}' {input.genomescope_report}`; "
         " smudgeplot.py plot -k {wildcards.kmer_length} -n ${{HAPLOID_COVERAGE}}  -o ${{PREFIX}} {input.coverages} > {log.plot} 2>&1; "
         " WARNINGS_NO_PRIORS={output.warnings_no_priors}; "
         " PREFIX_NO_PRIORS=${{WARNINGS_NO_PRIORS%_warnings.txt}}; "
         " smudgeplot.py plot -k {wildcards.kmer_length} -o ${{PREFIX_NO_PRIORS}} {input.coverages} > {log.plot} 2>&1; "


rule compress_kmer:
    input:
        kmer=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer",
        summary=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_summary_table.tsv"
    output:
        kmer_gz=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer.gz"
    log:
        std=output_dict["log"] / "compress_kmer.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.log",
        cluster_log=output_dict["cluster_log"] / "compress_kmer{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "compress_kmer.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "compress_kmer.{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("compress_kmer"),
        cpus=parameters["threads"]["compress_kmer"],
        time=parameters["time"]["compress_kmer"],
        mem=parameters["memory_mb"]["compress_kmer"],
        smudgeplot=1
    threads:
        parameters["threads"]["compress_kmer"]
    shell:
         " pigz -p 10 {input.kmer} > {log.std} 2>&1; "
