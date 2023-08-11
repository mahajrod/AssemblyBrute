rule gc_count:
    input:
        db=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"
    output:
        counts=output_dict["kmer"] / "{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.counts",
    params:
        tmp_dir=config["tmp_dir"]
        #kmer_length=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["kmer_length"],
        #min_coverage=lambda wildcards: "greater-than {0}".format(parameters["tool_options"]["gcp"][wildcards.datatype]["min_coverage"]) if "min_coverage" in parameters["tool_options"]["gcp"][wildcards.datatype] else "",
        #max_coverage=lambda wildcards: "less-than {0}".format(parameters["tool_options"]["gcp"][wildcards.datatype]["max_coverage"]) if "max_coverage" in parameters["tool_options"]["gcp"][wildcards.datatype] else "",
    log:
        gc_count=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.log",
        meryl=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.meryl.log",
        sort=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.sort.log",
        uniq=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.uniq.log",
        sed=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.sed.log",
        awk=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.awk.log",
        cluster_log=output_dict["cluster_log"] / "gc_plot.{datatype}.{stage}.L{min_coverage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gc_plot{datatype}.{stage}.L{min_coverage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["gc_count"],
        time=parameters["time"]["gc_count"],
        mem=parameters["memory_mb"]["gc_count"] + 30000,
    threads:
        parameters["threads"]["gc_count"]
    shell: # output: coverage\tgc\tcount\n
         " meryl threads={threads} memory={resources.mem}m greater-than {wildcards.min_coverage} "
         " print {input.db} 2>{log.meryl} | count_kmer_gc.py 2>{log.gc_count} | "
         " sort -S30000M -T {params.tmp_dir} -k2,2n -k1,1n 2>{log.sort} | "
         " uniq -c 2>{log.uniq} |  sed 's/^\s\+//;s/ /\\t/' 2>{log.sed} | "
         " awk '{{printf \"%i\\t%i\\t%i\\n\", $3,$2,$1 }}' > {output.counts} 2>{log.awk} "

rule gc_plot:
    input:
        counts=output_dict["kmer"] / "{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.counts",
        genomescope_report=output_dict["kmer"] / ("%s/{stage}/genomescope/%s.%s.filtered.%s.%s.genomescope.parameters" % (config["final_kmer_datatype"],
                                                                                                                          config["genome_prefix"],
                                                                                                                          config["final_kmer_datatype"],
                                                                                                                          config["final_kmer_length"],
                                                                                                                          config["final_kmer_counter"])),
        pip="results/config/pip.common.requirements" # added to ensure that distinctipy package was installed
    output:
        heatmap_png=output_dict["kmer"] / "{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.heatmap.png",
    params:
        ploidy=config["ploidy"],
    log:
        gc_count=output_dict["log"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.log",
        meryl=output_dict["log"] / "gc_plot.{datatype}.{stage}.{stage}.{kmer_length}.L{min_coverage}.log",
        cluster_log=output_dict["cluster_log"] / "gc_plot.{datatype}.{datatype}.{stage}.{kmer_length}.L{min_coverage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gc_plot.{datatype}.{datatype}.{stage}.{kmer_length}.L{min_coverage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["gc_plot"],
        time=parameters["time"]["gc_plot"],
        mem=parameters["memory_mb"]["gc_plot"],
    threads:
        parameters["threads"]["gc_plot"]
    shell:
         " LAMBDA=`awk 'NR==2 {{printf \"%.0f\", $2}}' {input.genomescope_report}`; "
         " HEATMAP_PNG_NAME={output.heatmap_png}; "
         " HEATMAP_PNG_PREFIX=${{HEATMAP_PNG_NAME%.heatmap.png}}; "
         " draw_gc_plot.py -i {input.counts}  -k {wildcards.kmer_length} -l ${{LAMBDA}} "
         " -p {params.ploidy} -m 4 -o ${{HEATMAP_PNG_PREFIX}} > {log.gc_count} 2>&1; " # -g 8 TODO: implement GC fraction calculation
