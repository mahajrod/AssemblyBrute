

rule flye:
    priority: 1000
    input:
        main_reads=get_main_read_filelist,

    output:
        hap0_graph=config["out_dir"] / "contig/{parameters, flye[^/]*}/{genome_prefix}.contig.hap0.unfiltered.gfa",
    params:
        input_type=lambda wildcards: "--" + stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["input_type"],
        polishing_iterations=lambda wildcards: parse_option("polishing_iterations", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --iterations "),
        read_error=lambda wildcards: parse_option("read_error", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --read-error "),
        keep_haplotypes=lambda wildcards: parse_option_flag("keep_haplotypes", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --keep-haplotypes "),
        no_alt_contigs=lambda wildcards: parse_option_flag("no_alt_contigs", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-alt-contigs "),
        scaffolding=lambda wildcards: parse_option_flag("scaffolding", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --scaffold "),
    log:
        std=config["out_dir"] / "log/flye.{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/flye.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/flye.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/flye.{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["flye"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["flye"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("flye"),
        cpus=get_threads(parameters["threads"]["flye"], "cpu"),
        time=parameters["time"]["flye"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["flye"], coeff=1.4, mode="exp"),
    threads:
        parameters["threads"]["flye"]
    shell:
         " OUT_DIR=`dirname {output.hap0_graph}`/assembly; "
         " mkdir -p ${{OUT_DIR}};"
         " flye -t {threads} {params.input_type} {input.main_reads} {params.polishing_iterations} {params.read_error} "
         "     {params.keep_haplotypes} {params.no_alt_contigs} {params.scaffolding} "
         "     -o ${{OUT_DIR}}  > {log.std} 2>&1; "
         " ln -sf assembly/assembly_graph.gfa {output.hap0_graph} >> {log.std} 2>&1; "
