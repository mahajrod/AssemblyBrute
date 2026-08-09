
localrules: extract_lambda_value

rule hifiasm_correct:
    priority: 2000
    input:
        main_reads=get_main_read_filelist_for_correction,
    output:
        ec_bin=config["out_dir"] / "error_correction/hifiasm_{correction_options}@{correction_mode}_mode/{genome_prefix}.contig.ec.bin",
        #ec_fasta=config["out_dir"] / "error_correction/hifiasm_{correction_options}@hifi_mode/{genome_prefix}.hifi.ec.fasta.gz",
        #alias_ec_fasta=config["out_dir"]  / "data/hifi/error_corrected_hifiasm_{correction_options}@hifi_mode/{genome_prefix}.hifi.ec.fasta.gz",
        ovlp_reverse_bin=config["out_dir"] / "error_correction/hifiasm_{correction_options}@{correction_mode}_mode/{genome_prefix}.contig.ovlp.reverse.bin",
        ovlp_source_bin=config["out_dir"] / "error_correction/hifiasm_{correction_options}@{correction_mode}_mode/{genome_prefix}.contig.ovlp.source.bin",
    params:
        window_size=lambda wildcards: parse_option("window_size", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -k "),
        D=lambda wildcards: parse_option("D", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -D "),
        N=lambda wildcards: parse_option("N", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -N "),
        ont_assembly=lambda wildcards: parse_option_flag("ont_mode", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --ont "),
        reads=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{1}_mode/{2}.contig.ec{3}".format(wildcards.correction_options,
                                                                                                                    wildcards.correction_mode,
                                                                                                                    wildcards.genome_prefix,
                                                                                                                    ".fq" if assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options']["ont_mode"] else ".fa"),
        renamed_reads=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{1}_mode/{2}.ec.{1}_mode{3}".format(wildcards.correction_options,
                                                                                                                               wildcards.correction_mode,
                                                                                                                               "_".join(set(assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options']["main_datatypes"]) & set(config["data"].keys())),
                                                                                                                               ".fq" if assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options']["ont_mode"] else ".fa",)

    log:
        std=config["out_dir"] / "log/hifiasm_correct.{correction_options}.{correction_mode}.{genome_prefix}.log",
        pigz=config["out_dir"] / "log/hifiasm_correct.{correction_options}.{correction_mode}.{genome_prefix}.pigz.log",
        mv=config["out_dir"] / "log/hifiasm_correct.{correction_options}.{correction_mode}.{genome_prefix}.mv.log",
        ln=config["out_dir"] / "log/hifiasm_correct.{correction_options}.{correction_mode}.{genome_prefix}.ln.log",
        cluster_log=config["out_dir"] / "log/hifiasm_correct{correction_options}.{correction_mode}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/hifiasm_correct{correction_options}.{correction_mode}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/hifiasm_correct.{correction_options}.{correction_mode}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("hifiasm_correct"),
        cpus=get_threads(parameters["threads"]["hifiasm"], "cpu"),
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.ec_bin}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.ec.bin}}; "
         " hifiasm -t {threads} -e --write-ec {params.window_size} {params.bloom_filter_bits} {params.ont_assembly} "
         "     {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size} "
         "     {params.kmer_length} {params.D} {params.N} -o ${{OUTPUT_PREFIX}} {input.main_reads} > {log.std} 2>&1; "
         " mv {params.reads} {params.renamed_reads} > {log.mv} 2>&1; "
         " pigz -p {threads} {params.renamed_reads} > {log.pigz} 2>&1 ; "


rule hifiasm_hic_2p:
    priority: 1000
    input:
        main_reads=get_main_read_filelist,
        ultralong_reads=get_ultralong_read_filelist,
        hic_forward=lambda wildcards: expand(config["out_dir"] / "data/hic/final/{pairprefix}{forward_suffix}{extension}",
                                             forward_suffix=[config["data"]["hic"]["conv_fwd_sfx"], ],
                                             pairprefix=config["data"]["hic"]["pair_prefix_list"],
                                             extension=[config["data"]["hic"]["conv_ext"]]) if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_hic"] else [],
        hic_reverse=lambda wildcards: expand(config["out_dir"] / "data/hic/final/{pairprefix}{reverse_suffix}{extension}",
                                             reverse_suffix=[config["data"]["hic"]["conv_rev_sfx"], ],
                                             pairprefix=config["data"]["hic"]["pair_prefix_list"],
                                             extension=[config["data"]["hic"]["conv_ext"]]) if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_hic"] else [],
        ec_bin=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{2}_mode/{1}.contig.ec.bin".format(stage_dict["contig"].parameters[wildcards.parameters]["option_set_group"],
                                                                                                                 wildcards.genome_prefix,
                                                                                                                 "nano" if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"] else "hifi") ,
        ovlp_reverse_bin=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{2}_mode/{1}.contig.ovlp.reverse.bin".format(stage_dict["contig"].parameters[wildcards.parameters]["option_set_group"],
                                                                                                                                      wildcards.genome_prefix,
                                                                                                                                      "nano" if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"] else "hifi"),
        ovlp_source_bin=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{2}_mode/{1}.contig.ovlp.source.bin".format(stage_dict["contig"].parameters[wildcards.parameters]["option_set_group"],
                                                                                                                                   wildcards.genome_prefix,
                                                                                                                                   "nano" if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"] else "hifi"),
        lambda_file=rules.extract_lambda_value.output.lambda_file
    output:
        hap1_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p2}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        hap2_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p2}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p2}/{genome_prefix}.contig.hic.a_ctg.gfa",

    params:
        output_prefix=lambda wildcards: config["out_dir"] / "contig/{0}/{1}.contig".format(wildcards.parameters, wildcards.genome_prefix),
        purge_level=lambda wildcards: stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["purge_level"],
        ploidy=lambda wildcards: stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["assembly_ploidy"], #config["ploidy"],
        cov_multiplicator=lambda wildcards: stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["cov_multiplicator"],
        window_size=lambda wildcards: parse_option("window_size", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -k "),
        dual_scaf=lambda wildcards: parse_option_flag("dual_scaf", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --dual-scaf "),
        D=lambda wildcards: parse_option("D", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -D "),
        N=lambda wildcards: parse_option("N", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -N "),
        ignore_bin=lambda wildcards: " -i " if ("ignore_bin" in stage_dict["contig"].parameters[wildcards.parameters]["option_set"]) and stage_dict["contig"].parameters[wildcards.parameters]["option_set"] else "",
        sim_threshold_for_hapdup_reads=lambda wildcards: parse_option("sim_threshold_for_hapdup_reads", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -s "),
        hic_forward=lambda wildcards: (" --h1 " + ",".join(map(str, expand(config["out_dir"] / "data/hic/final/{pairprefix}{forward_suffix}{extension}",
                                                         forward_suffix=[config["data"]["hic"]["conv_fwd_sfx"], ],
                                                         pairprefix=config["data"]["hic"]["pair_prefix_list"],
                                                         extension=[config["data"]["hic"]["conv_ext"]] )))) if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_hic"] and ("hic" in config["data"]) else [], #in case of multiple hic libraries files in the list MUST be COMMA-separated
        hic_reverse=lambda wildcards: (" --h2 " + ",".join(map(str, expand(config["out_dir"] / "data/hic/final/{pairprefix}{reverse_suffix}{extension}",
                                                         reverse_suffix=[config["data"]["hic"]["conv_rev_sfx"], ],
                                                         pairprefix=config["data"]["hic"]["pair_prefix_list"],
                                                         extension=[config["data"]["hic"]["conv_ext"]] )))) if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_hic"] and ("hic" in config["data"]) else [],
        ultralong_reads=lambda wildcards: (" --ul " + ",".join(map(str, get_ultralong_read_filelist(wildcards)))) if get_ultralong_read_filelist(wildcards) else "",
        telomere_motif= lambda wildcards: parse_option("telomere_motif", config, " --telo-m ") if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["use_telomere"] else "",

        ul_cut=lambda wildcards: parse_option("ul-cut", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " --ul-cut "),
        ont_assembly= lambda wildcards: parse_option_flag("ont_mode", stage_dict["contig"].parameters[wildcards.parameters]["option_set"]," --ont "),
        ont_mode=lambda wildcards: stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"],
        post_join=lambda wildcards: parse_option("post_join", stage_dict["contig"].parameters[wildcards.parameters]["option_set"], " -u ")
    log:
        std=config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/hifiasm.{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("hifiasm_hic"),
        cpus=get_threads(parameters["threads"]["hifiasm"], "cpu"),
        time=parameters["time"]["hifiasm"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["hifiasm"], coeff=1.4, mode="exp"),
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={params.output_prefix}; "
         " OUT_DIR=`dirname ${{OUTPUT_PREFIX}}`; "
         " if [[ '{params.ont_mode}' != 'True' ]]; "
         "     then "
         "     ln -sf ../../../{input.ec_bin} ${{OUT_DIR}} 1>{log.std} 2>&1; "
         "     ln -sf ../../../{input.ovlp_reverse_bin} ${{OUT_DIR}} 1>>{log.std} 2>&1; "
         "     ln -sf ../../../{input.ovlp_source_bin} ${{OUT_DIR}} 1>>{log.std} 2>&1; "
         "     fi; "
         " LAMBDA=`head -n 1 {input.lambda_file}` 1>>{log.std} 2>&1;  "
         " COV_UPPER_BOUNDARY=`echo \"{params.cov_multiplicator}*${{LAMBDA}}\" | bc` >>{log.std} 2>&1;  "
         " COV_UPPER_BOUNDARY=${{COV_UPPER_BOUNDARY%.*}}; "
         " hifiasm {params.window_size} {params.bloom_filter_bits} {params.ont_assembly} "
         "     {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         "     {params.kmer_length} {params.D} {params.N} {params.ignore_bin} {params.sim_threshold_for_hapdup_reads} "
         "     --primary -t {threads} -l {params.purge_level}  -o ${{OUTPUT_PREFIX}} "
         "     --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         "     {params.hic_forward} {params.hic_reverse} {params.ultralong_reads} {params.ul_cut} {params.dual_scaf} "
         "     {params.telomere_motif} {params.post_join} "
         "     {input.main_reads}  >>{log.std} 2>&1; "

use rule hifiasm_hic_2p as hifiasm_hic_3p with:
    output:
        hap1_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p3}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        hap2_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p3}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        hap3_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p3}/{genome_prefix}.contig.hic.hap3.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p3}/{genome_prefix}.contig.hic.a_ctg.gfa",

use rule hifiasm_hic_2p as hifiasm_hic_4p with:
    output:
        hap1_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p4}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        hap2_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p4}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        hap3_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p4}/{genome_prefix}.contig.hic.hap3.p_ctg.gfa",
        hap4_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p4}/{genome_prefix}.contig.hic.hap4.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p4}/{genome_prefix}.contig.hic.a_ctg.gfa",

use rule hifiasm_hic_2p as hifiasm_hic_5p with:
    output:
        hap1_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p5}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        hap2_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p5}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        hap3_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p5}/{genome_prefix}.contig.hic.hap3.p_ctg.gfa",
        hap4_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p5}/{genome_prefix}.contig.hic.hap4.p_ctg.gfa",
        hap5_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p5}/{genome_prefix}.contig.hic.hap5.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p5}/{genome_prefix}.contig.hic.a_ctg.gfa",

use rule hifiasm_hic_2p as hifiasm_hic_6p with:
    output:
        hap1_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.hap1.p_ctg.gfa",
        hap2_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.hap2.p_ctg.gfa",
        hap3_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.hap3.p_ctg.gfa",
        hap4_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.hap4.p_ctg.gfa",
        hap5_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.hap5.p_ctg.gfa",
        hap6_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.hap6.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p6}/{genome_prefix}.contig.hic.a_ctg.gfa",

use rule hifiasm_hic_2p as hifiasm_long_reads_only with:
    input:
        main_reads=get_main_read_filelist,
        ultralong_reads=get_ultralong_read_filelist,
        hic_forward=[],
        hic_reverse=[],
        ec_bin=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{2}_mode/{1}.contig.ec.bin".format(stage_dict["contig"].parameters[wildcards.parameters]["option_set_group"],
                                                                                                                 wildcards.genome_prefix,
                                                                                                                 "nano" if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"] else "hifi") ,
        ovlp_reverse_bin=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{2}_mode/{1}.contig.ovlp.reverse.bin".format(stage_dict["contig"].parameters[wildcards.parameters]["option_set_group"],
                                                                                                                                      wildcards.genome_prefix,
                                                                                                                                      "nano" if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"] else "hifi"),
        ovlp_source_bin=lambda wildcards: config["out_dir"] / "error_correction/hifiasm_{0}@{2}_mode/{1}.contig.ovlp.source.bin".format(stage_dict["contig"].parameters[wildcards.parameters]["option_set_group"],
                                                                                                                                   wildcards.genome_prefix,
                                                                                                                                   "nano" if stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["ont_mode"] else "hifi"),
        lambda_file=rules.extract_lambda_value.output.lambda_file

    output:
        primary_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p1}/{genome_prefix}.contig.p_ctg.gfa",
        alt_contig_graph=config["out_dir"] / "contig/{parameters, hifiasm_[^/]*@p1}/{genome_prefix}.contig.a_ctg.gfa",

ruleorder: create_haploid_gfa_link > create_primary_gfa_link
localrules: create_primary_gfa_link, create_alt_gfa_link, create_haploid_gfa_link, create_haploid_alt_gfa_link

use rule create_local_links as create_primary_gfa_link with:
    input:
        gfa=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.hic.{haplotype}.p_ctg.gfa",
        log_dir=ancient(config["out_dir"] / "contig/{parameters}/log")
    output:
        gfa=config["out_dir"] / "contig/{parameters, hifiasm[^/]*}/{genome_prefix}.contig.{haplotype, hap[^/]+}.unfiltered.gfa"
    log:
        ln=config["out_dir"] / "contig/{parameters}/log/create_gfa_links.{parameters}.{genome_prefix}.contig.{haplotype}.ln.log",

use rule create_local_links as create_alt_gfa_link with:
    input:
        gfa=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.hic.a_ctg.gfa",
        log_dir=ancient(config["out_dir"] / "contig/{parameters}/log")
    output:
        gfa=config["out_dir"] / "contig/{parameters, hifiasm[^/]*}/{genome_prefix}.contig.alt.unfiltered.gfa"
    log:
        ln=config["out_dir"] / "contig/{parameters}/log/create_gfa_links.{parameters}.{genome_prefix}.contig.alt.ln.log",

use rule create_local_links as create_haploid_gfa_link with:
    input:
        gfa=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.p_ctg.gfa",
        log_dir=ancient(config["out_dir"] / "contig/{parameters}/log")
    output:
        gfa=config["out_dir"] / "contig/{parameters, hifiasm_.*@p1}/{genome_prefix}.contig.hap0.unfiltered.gfa"
    log:
        ln=config["out_dir"] / "contig/{parameters}/log/create_gfa_links.{parameters}.{genome_prefix}.contig.hap0.ln.log",

use rule create_local_links as create_haploid_alt_gfa_link with:
    input:
        gfa=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.a_ctg.gfa",
        log_dir=ancient(config["out_dir"] / "contig/{parameters}/log")
    output:
        gfa=config["out_dir"] / "contig/{parameters, hifiasm_.*@p1}/{genome_prefix}.contig.alt0.unfiltered.gfa"
    log:
        ln=config["out_dir"] / "contig/{parameters}/log/create_gfa_links.{parameters}.{genome_prefix}.contig.alt0.ln.log",


rule get_length_and_coverage_from_hifiasm_graph:
    input:
        gfa="{gfa_dir}/{gfa_prefix}.gfa",
        log_dir=ancient("{gfa_dir}/log/")
    output:
        cov="{gfa_dir, .*hifiasm.*}/{gfa_prefix}.gfa.cov",
        len_cov="{gfa_dir}/{gfa_prefix}.gfa.lencov",

    log:
        std="{gfa_dir}/log/{gfa_prefix}.get_length_and_coverage_from_hifiasm_graph.log",
        cluster_log="{gfa_dir}/log/{gfa_prefix}.get_length_and_coverage_from_hifiasm_graph.cluster.log",
        cluster_err="{gfa_dir}/log/{gfa_prefix}.get_length_and_coverage_from_hifiasm_graph.cluster.err"
    benchmark:
        "{gfa_dir}/log/{gfa_prefix}.get_length_and_coverage_from_hifiasm_graph.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_length_and_coverage_from_hifiasm_graph"),
        cpus=parameters["threads"]["get_coverage_from_hifiasm_graph"],
        time=parameters["time"]["get_coverage_from_hifiasm_graph"],
        mem=parameters["memory_mb"]["get_coverage_from_hifiasm_graph"],
    threads:
        parameters["threads"]["get_coverage_from_hifiasm_graph"]
    shell:
         " workflow/scripts/extract_length_and_coverage_from_hifiasm_gfa.bash {input.gfa} > {output.len_cov} 2>{log.std}; "
         " cut -f 1,3 {output.len_cov} > {output.cov} 2>>{log.std}; "


"""
rule get_lowcoverage_contig_ids:
    input:
        cov=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.unfiltered.gfa.cov"
    output:
        low_cov_ids=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype, [^/]+}.unfiltered.gfa.lowcov.ids",
    log:
        std=config["out_dir"] / "log/get_lowcoverage_contig_ids.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=config["out_dir"] / "log/get_lowcoverage_contig_ids.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/get_lowcoverage_contig_ids.{parameters}.{genome_prefix}.{haplotype}.cluster.err",
    params:
        min_coverage=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["min_contig_coverage"]
    benchmark:
        config["out_dir"] / "log/get_lowcoverage_contig_ids.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_lowcoverage_contig_ids"),
        cpus=parameters["threads"]["get_lowcoverage_contig_ids"],
        time=parameters["time"]["get_lowcoverage_contig_ids"],
        mem=parameters["memory_mb"]["get_lowcoverage_contig_ids"],
    threads:
        parameters["threads"]["get_lowcoverage_contig_ids"]
    shell:
         " if [ '{wildcards.haplotype}' == 'alt' ];"
         " then "
         "   > {output.low_cov_ids}; "
         " else "
         "   awk '{{if ($2 < {params.min_coverage}) print $1}}' {input.cov} > {output.low_cov_ids} 2>{log.std}; "
         " fi; "

#ruleorder: filter_contigs_by_coverage > gfa2fasta

rule filter_contigs_by_coverage:
    input:
        low_cov_ids=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.unfiltered.gfa.lowcov.ids",
        unfiltered_fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.unfiltered.fasta"
    output:
        filtered_fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype, [^/]+}.lenfiltered.fasta",
    log:
        std=config["out_dir"] / "log/filter_contigs_by_coverage.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=config["out_dir"] / "log/filter_contigs_by_coverage.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/filter_contigs_by_coverage.{parameters}.{genome_prefix}.{haplotype}.cluster.err",
    benchmark:
        config["out_dir"] / "log/filter_contigs_by_coverage.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("filter_contigs_by_coverage"),
        cpus=parameters["threads"]["filter_contigs_by_coverage"],
        time=parameters["time"]["filter_contigs_by_coverage"],
        mem=parameters["memory_mb"]["filter_contigs_by_coverage"],
    threads:
        parameters["threads"]["filter_contigs_by_coverage"]
    shell:
         " extract_sequences_by_ids.py -i {input.unfiltered_fasta} -d {input.low_cov_ids} -r "
         " -o {output.filtered_fasta} > {log.std} 2>&1; "
"""
