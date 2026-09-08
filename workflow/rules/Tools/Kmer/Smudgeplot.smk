
# here a FastK based version of the smudgeplot is used

localrules: smudgeplot_assess

rule smudgeplot_assess:
    input:
        fastk_histo=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.histo",
    output:
        boundaries=config["out_dir"] / "kmer/{datatype}/{stage}/smudgeplot/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/smudgeplot.boundaries",
    log:
        upper=config["out_dir"] / "log/smudgeplot_assess.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.upper.log",
        lower=config["out_dir"] / "log/smudgeplot_assess.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.lower.log",
        cluster_log=config["out_dir"] / "log/smudgeplot_assess.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/smudgeplot_assess.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/smudgeplot_assess.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("smudgeplot_assess"),
        cpus=parameters["threads"]["smudgeplot_assess"],
        time=parameters["time"]["smudgeplot_assess"],
        mem=parameters["memory_mb"]["smudgeplot_assess"],
    threads:
        parameters["threads"]["smudgeplot_assess"]
    shell:
         " LOWER_BOUNDARY=`smudgeplot cutoff {input.fastk_histo} L 2>{log.lower}`; "
         " UPPER_BOUNDARY=`smudgeplot cutoff {input.fastk_histo} U 2>{log.upper}`; "
         " echo -e \"low_boundary\tupper_boundary\n${{LOWER_BOUNDARY}}\t${{UPPER_BOUNDARY}}\" > {output.boundaries}"

rule smudgeplot_hetmers:
    input:
        fastk_db=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/",
        boundaries=config["out_dir"] / "kmer/{datatype}/{stage}/smudgeplot/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/smudgeplot.boundaries"
    output:
        smu_file=config["out_dir"] / "kmer/{datatype}/{stage}/smudgeplot/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/smudgeplot_hetmers.smu",
    log:
        std=config["out_dir"] / "log/smudgeplot_hetmers.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.std.log",
        cluster_log=config["out_dir"] / "log/smudgeplot_hetmers.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/smudgeplot_hetmers.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}..cluster.err"
    benchmark:
        config["out_dir"] / "log/smudgeplot_hetmers.{datatype}.{stage}.{kmer_length}.min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("smudgeplot_hetmers"),
        cpus=parameters["threads"]["smudgeplot_hetmers"],
        time=parameters["time"]["smudgeplot_hetmers"],
        mem=parameters["memory_mb"]["smudgeplot_hetmers"],
        smudgeplot_hetkmers=1
    threads:
        parameters["threads"]["smudgeplot_hetmers"]
    shell:
         " OUT_DIR=`dirname {output.smu_file}`; "
         " OUT_PREFIX={output.smu_file}; "
         " OUT_PREFIX=${{OUT_PREFIX%.smu}}; "
         " TMP_DIR=${{OUT_DIR}}/tmp/; "
         " mkdir -p ${{TMP_DIR}}; "
         " smudgeplot hetmers -L `cut -f 1 {input.boundaries} | sed -n 2p` -t {threads} -o ${{OUT_PREFIX}} "
         "      --verbose -tmp ${{TMP_DIR}}  {input.fastk_db}/fastk_db.ktab > {log.std} 2>&1; "
         " rm -r ${{TMP_DIR}}; "

rule smudgeplot:
    input:
        smu_file=config["out_dir"] / "kmer/{datatype}/{stage}/smudgeplot/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/smudgeplot_hetmers.smu",
        boundaries=config["out_dir"] / "kmer/{datatype}/{stage}/smudgeplot/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/smudgeplot.boundaries"
    output:
        centralities_file=config["out_dir"] / "kmer/{datatype}/{stage}/smudgeplot/{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}/smudgeplot_hetmers_centralities.txt"
    log:
        plot=config["out_dir"] / "log/smudgeplot.{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.plot.log",
        cluster_log=config["out_dir"] / "log/smudgeplot.{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.cluster.log",
        cluster_err=config["out_dir"] / "log/smudgeplot.{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.cluster.err"
    benchmark:
        config["out_dir"] / "log/smudgeplot.{datatype}.{stage}.{kmer_length}.fastk_min{min_kmer_count}.benchmark.txt"
    conda:
        config["conda"]["smudgeplot"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["smudgeplot"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("smudgeplot"),
        cpus=parameters["threads"]["smudgeplot"],
        time=parameters["time"]["smudgeplot"],
        mem=parameters["memory_mb"]["smudgeplot"],
    threads:
        parameters["threads"]["smudgeplot"]
    shell:
         " OUT_PREFIX={output.centralities_file}; "
         " OUT_PREFIX=${{OUT_PREFIX%_centralities.txt}}; "
         " smudgeplot all -t {wildcards.datatype}.{wildcards.stage}.{wildcards.kmer_length}.fastk_min{wildcards.min_kmer_count} "
         "      -cov_min `cut -f 1 {input.boundaries} | sed -n 2p` "
         "      -cov_max `cut -f 2 {input.boundaries} | sed -n 2p` "
         "      -ylim `cut -f 2 {input.boundaries} | sed -n 2p` "
         "      -o ${{OUT_PREFIX}} {input.smu_file} > {log.plot} 2>&1; "



