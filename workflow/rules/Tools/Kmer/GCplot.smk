rule gc_count:
    input:
        db=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"
    output:
        counts=config["out_dir"] / "kmer/{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.counts",
    params:
        tmp_dir=config["tmp_dir"]
    log:
        gc_count=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.log",
        meryl=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.meryl.log",
        sort=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.sort.log",
        uniq=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.uniq.log",
        sed=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.sed.log",
        awk=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.awk.log",
        cluster_log=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.L{min_coverage}.{kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/gc_plot{datatype}.{stage}.L{min_coverage}.{kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("gc_count"),
        cpus=parameters["threads"]["gc_count"],
        time=parameters["time"]["gc_count"],
        mem=parameters["memory_mb"]["gc_count"] + 30000,
    threads:
        parameters["threads"]["gc_count"]
    shell: # output: coverage\tgc\tcount\n
         " meryl threads={threads} memory={resources.mem}m greater-than {wildcards.min_coverage} "
         "     print {input.db} 2>{log.meryl} | count_kmer_gc.py 2>{log.gc_count} | "
         "     sort -S30000M -T {params.tmp_dir} -k2,2n -k1,1n 2>{log.sort} | "
         "     uniq -c 2>{log.uniq} |  sed 's/^\s\+//;s/ /\\t/' 2>{log.sed} | "
         "     awk '{{printf \"%i\\t%i\\t%i\\n\", $3,$2,$1 }}' > {output.counts} 2>{log.awk} "

rule gc_plot:
    input:
        counts=config["out_dir"] / "kmer/{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.counts",
        genomescope_report=config["out_dir"] / ("kmer/{datatype}/{stage}/genomescope/%s.{datatype}.final.%s.%s.genomescope.parameters" % (
                                                                                                                          config["genome_prefix"],
                                                                                                                          config["final_kmer_length"],
                                                                                                                          config["final_kmer_counter"])),
    output:
        heatmap_png=config["out_dir"] / "kmer/{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.heatmap.png",
    params:
        ploidy=config["ploidy"],
    log:
        gc_count=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.log",
        meryl=config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{stage}.{kmer_length}.L{min_coverage}.log",
        cluster_log=config["out_dir"] / "log/gc_plot.{datatype}.{datatype}.{stage}.{kmer_length}.L{min_coverage}.cluster.log",
        cluster_err=config["out_dir"] / "log/gc_plot.{datatype}.{datatype}.{stage}.{kmer_length}.L{min_coverage}.cluster.err"
    benchmark:
        config["out_dir"] / "log/gc_plot.{datatype}.{stage}.{kmer_length}.L{min_coverage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("gc_plot"),
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
         "     -p {params.ploidy} -m 4 -o ${{HEATMAP_PNG_PREFIX}} > {log.gc_count} 2>&1; " # -g 8 TODO: implement GC fraction calculation
