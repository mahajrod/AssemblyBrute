

"""
rule nextdenovo: # TODO: implement
    priority: 1000
    input:
        main_reads=get_main_read_filelist,
    output:
        hap1_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_.*@p2}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        hap2_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_.*@p2}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_.*@p2}/{genome_prefix}.contig.hic.a_ctg.gfa",

    log:
        std=config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["nextdenovo"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nextdenovo"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("nextdenovo_hic"),
        cpus=get_threads(parameters["threads"]["nextdenovo"], "cpu"),
        time=parameters["time"]["nextdenovo"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["nextdenovo"], coeff=1.4, mode="exp"), #parameters["memory_mb"]["hifiasm"]
    threads:
        parameters["threads"]["nextdenovo"]
    shell:
         " "
"""