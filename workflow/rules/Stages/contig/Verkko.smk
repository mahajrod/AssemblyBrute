

rule verkko_phasing:
    priority: 1000
    input:
        hifi_reads=get_hifi_read_filelist,
        nano_reads=get_nano_read_filelist,
        hic_forward=lambda wildcards: expand(config["out_dir"] / "data/hic/final/{pairprefix}{forward_suffix}{extension}",
                                             forward_suffix=[config["data"]["hic"]["conv_fwd_sfx"], ],
                                             pairprefix=config["data"]["hic"]["pair_prefix_list"],
                                             extension=[config["data"]["hic"]["conv_ext"]]) if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_hic"] else [],
        hic_reverse=lambda wildcards: expand(config["out_dir"] / "data/hic/final/{pairprefix}{reverse_suffix}{extension}",
                                             reverse_suffix=[config["data"]["hic"]["conv_rev_sfx"], ],
                                             pairprefix=config["data"]["hic"]["pair_prefix_list"],
                                             extension=[config["data"]["hic"]["conv_ext"]]) if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_hic"] else [],
    output:
        consensus_fasta=config["out_dir"] / "contig/{parameters, verkko[^/]*@p2}/{genome_prefix}.contig.consensus.unfiltered.fasta",
        hap1_fasta=config["out_dir"] / "contig/{parameters, verkko[^/]*@p2}/{genome_prefix}.contig.hap1.unfiltered.fasta",
        hap2_fasta=config["out_dir"] / "contig/{parameters, verkko[^/]*@p2}/{genome_prefix}.contig.hap2.unfiltered.fasta",
    params:
        hifi_reads=lambda wildcards: " --hifi %s " % " ".join(get_hifi_read_filelist(wildcards)) if get_hifi_read_filelist(wildcards) else "",
        nano_reads=lambda wildcards: " --nano %s " % " ".join(get_nano_read_filelist(wildcards)) if get_nano_read_filelist(wildcards) else "",
        haplo_divergence=lambda wildcards: parse_option("polishing_iterations", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --haplo-divergence "),
        uneven_cov=lambda wildcards: parse_option_flag("uneven_cov", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --uneven-depth "),
        no_rdna_tangle=lambda wildcards: parse_option_flag("no_rdna_tangle", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-rdna-tangle "),
        no_nano=lambda wildcards: parse_option_flag("no_nano", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-nano "),
        no_correction=lambda wildcards: parse_option_flag("no_correction", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-correction "),
        telomere_motif= lambda wildcards: parse_option("telomere_motif", config, " --telomere-motif ") if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_telomere"] else "",
    log:
        std=config["out_dir"] / "log/verkko.{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/verkko.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/verkko.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/flye.{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["verkko"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["verkko"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("verkko"),
        cpus=get_threads(parameters["threads"]["verkko"], "cpu"),
        time=parameters["time"]["verkko"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["verkko"], coeff=1.4, mode="exp"),
    threads:
        parameters["threads"]["verkko"]
    shell:
         " verkko {params.haplo_divergence} {params.uneven_cov} {params.no_rdna_tangle} {params.no_nano} "
         "        {params.no_correction} {params.telomere_motif} "
         "        -d results/verkko/ "
         "        {params.nano_reads} {params.hifi_reads} "
         "        --hic1 {input.hic_forward} --hic2 {input.hic_reverse} "
         "        --local --local-cpus 16 > {log.std} 2>&1;"
         " ln -sf assembly.fasta `basename {output.consensus_fasta}`; "
         " ln -sf assembly.haplotype1.fasta `basename {output.hap1_fasta}`; "
         " ln -sf assembly.haplotype1.fasta `basename {output.hap2_fasta}`; "

rule verkko_no_phasing:
    priority: 1000
    input:
        hifi_reads=get_hifi_read_filelist,
        nano_reads=get_nano_read_filelist,
    output:
        hap0_fasta=config["out_dir"] / "contig/{parameters, verkko[^/]*@p1}/{genome_prefix}.contig.hap0.unfiltered.fasta",
    params:
        hifi_reads=lambda wildcards: " --hifi %s " % " ".join(get_hifi_read_filelist(wildcards)) if get_hifi_read_filelist(wildcards) else "",
        nano_reads=lambda wildcards: " --nano %s " % " ".join(get_nano_read_filelist(wildcards)) if get_nano_read_filelist(wildcards) else "",
        haplo_divergence=lambda wildcards: parse_option("polishing_iterations", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --haplo-divergence "),
        uneven_cov=lambda wildcards: parse_option_flag("uneven_cov", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --uneven-depth "),
        no_rdna_tangle=lambda wildcards: parse_option_flag("no_rdna_tangle", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-rdna-tangle "),
        no_nano=lambda wildcards: parse_option_flag("no_nano", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-nano "),
        no_correction=lambda wildcards: parse_option_flag("no_correction", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --no-correction "),
        telomere_motif= lambda wildcards: parse_option("telomere_motif",config," --telomere-motif ") if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_telomere"] else "",

    log:
        std=config["out_dir"] / "log/verkko.{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/verkko.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/verkko.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/flye.{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["verkko"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["verkko"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("verkko"),
        cpus=get_threads(parameters["threads"]["verkko"], "cpu"),
        time=parameters["time"]["verkko"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["verkko"], coeff=1.4, mode="exp"),
    threads:
        parameters["threads"]["verkko"]
    shell:
         " verkko {params.haplo_divergence} {params.uneven_cov} {params.no_rdna_tangle} {params.no_nano} "
         "        {params.no_correction} {params.telomere_motif} "
         "        -d results/verkko/ "
         "        {params.hifi_reads} {params.nano_reads} "
         "        --local --local-cpus 16 > {log.std} 2>&1; "
         " ln -sf assembly.fasta `basename {output.hap0_fasta}`; "
