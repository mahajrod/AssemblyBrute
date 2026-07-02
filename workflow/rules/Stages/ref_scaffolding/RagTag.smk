localrules: create_links_ragtag_scaffolds

rule ragtag_ref_scaf: #
    input:
        fasta=config["out_dir"] / ("%s/{prev_stage_parameters}/%s.%s.{haplotype}.fasta" % (stage_dict["ref_scaffolding"].prev_stage,
                                                                                           config["genome_prefix"],
                                                                                           stage_dict["ref_scaffolding"].prev_stage)),
        contig_both_blacklist=config["out_dir"] / ("%s/{prev_stage_parameters}/assembly_qc/%s.%s.{haplotype}/telomere/%s.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.both.ids" %
                                                  (stage_dict["ref_scaffolding"].prev_stage,
                                                   config["genome_prefix"],
                                                   stage_dict["ref_scaffolding"].prev_stage,
                                                   config["genome_prefix"],
                                                   stage_dict["ref_scaffolding"].prev_stage,
                                                   "canonical" if config["use_canonical_telomere"] else "non_canonical"
                                                   )),
        reference_fasta=config["out_dir"] / "data/reference/{reference}/{reference}/masking/{reference}.softmasked.fasta",
    output:
        ragtag_fasta=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters, [^/]+}@{reference}/{haplotype}/ragtag.scaffold.fasta",
        ragtag_agp=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters, [^/]+}@{reference}/{haplotype}/ragtag.scaffold.agp",
    params:
        min_aln_len=lambda wildcards: parse_option("min_aln_len", parameters["tool_options"]["ragtag"][wildcards.ref_scaf_parameters], " -f ")

    log:
        ragtag=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{haplotype}.log",
        ln=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{haplotype}.ln.log",
        awk=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{haplotype}.awk.log",
        cluster_log=config["out_dir"] / "log/ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["ragtag"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["ragtag"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("ragtag"),
        cpus=parameters["threads"]["ragtag"],
        time=parameters["time"]["ragtag"],
        mem=parameters["memory_mb"]["ragtag"]
    threads: parameters["threads"]["ragtag"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_fasta}`; "
        " echo 'Scaffolding...' > {log.ragtag} 2>&1;"
        " ragtag.py scaffold -t {threads} -j {input.contig_both_blacklist} {params.min_aln_len} "
        "     -w {input.reference_fasta} {input.fasta} -o ${{RAGTAG_DIR}}  >> {log.ragtag} 2>&1; "


rule ragtag_filter: #
    input:
        ragtag_agp=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.scaffold.agp",
        contig_five_prime_blacklist=config["out_dir"] / ("%s/{prev_stage_parameters}/assembly_qc/telomere/{genome_prefix}.%s.{haplotype}/{genome_prefix}.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids" %
                                                  (stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   "canonical" if config["use_canonical_telomere"] else "non_canonical"
                                                   )),
        contig_three_prime_blacklist=config["out_dir"] / ("%s/{prev_stage_parameters}/assembly_qc/telomere/{genome_prefix}.%s.{haplotype}/{genome_prefix}.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids" %
                                                  (stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   "canonical" if config["use_canonical_telomere"] else "non_canonical"
                                                   )),
    output:
        ragtag_filtered_agp=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters, [^/]+}@{reference}/{haplotype}/{genome_prefix}.ref_scaffolding.{haplotype}.agp",
        ragtag_filtered_stats=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters, [^/]+}@{reference}/{haplotype}/{genome_prefix}.ref_scaffolding.{haplotype}.stats",
    log:
        ragtag=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag_filter.{genome_prefix}.ragtag.{haplotype}.log",
        ln=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag_filter.{genome_prefix}.ragtag.{haplotype}.ln.log",
        cluster_log=config["out_dir"] / "log/ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("ragtag_filter"),
        cpus=parameters["threads"]["ragtag_filter"],
        time=parameters["time"]["ragtag_filter"],
        mem=parameters["memory_mb"]["ragtag_filter"]
    threads: parameters["threads"]["ragtag_filter"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_filtered_agp}`; "
        " echo 'Spliting agp at internal telomeres...' > {log.ragtag} 2>&1; "
        " workflow/scripts/curation/split_agp_by_telomeric_merges.py -a {input.ragtag_agp} "
        "     --five_prime_blacklist {input.contig_five_prime_blacklist} "
        "     --three_prime_blacklist {input.contig_three_prime_blacklist} -p ${{RAGTAG_DIR}}/ragtag.scaffold >> {log.ragtag} 2>&1; "
        " ln -sf ragtag.scaffold.telomere.split.final.renamed.agp {output.ragtag_filtered_agp} >> {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.telomere.split.final.renamed.stats {output.ragtag_filtered_stats} >> {log.ln} 2>&1; "

rule ragtag_generate_filtered_fasta: #
    input:
        ragtag_filtered_agp=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/{genome_prefix}.ref_scaffolding.{haplotype}.agp",
        fasta=config["out_dir"] / ("%s/{prev_stage_parameters}/%s.%s.{haplotype}.fasta" % (stage_dict["ref_scaffolding"]["prev_stage"],
                                                                                      config["genome_prefix"],
                                                                                      stage_dict["ref_scaffolding"]["prev_stage"])),
    output:
        ragtag_filtered_fasta=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters, [^/]+}@{reference}/{haplotype}/{genome_prefix}.ref_scaffolding.{haplotype}.fasta",
    log:
        ragtag=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag_filter.{genome_prefix}.ragtag.{haplotype}.log",
        ln=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag_filter.{genome_prefix}.ragtag.{haplotype}.ln.log",
        cluster_log=config["out_dir"] / "log/ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["ragtag"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["ragtag"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("ragtag_filter"),
        cpus=parameters["threads"]["ragtag_filter"],
        time=parameters["time"]["ragtag_filter"],
        mem=parameters["memory_mb"]["ragtag_filter"]
    threads: parameters["threads"]["ragtag_filter"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_filtered_fasta}`; "
        " echo 'Assembling fasta from modified agp...' > {log.ragtag} 2>&1; "
        " ragtag.py agp2fa {input.ragtag_filtered_agp} {input.fasta} > ${{RAGTAG_DIR}}/ragtag.scaffold.telomere.split.final.renamed.fasta 2>>{log.ragtag}; "
        " ln -sf ragtag.scaffold.telomere.split.final.renamed.fasta {output.ragtag_filtered_fasta} > {log.ln} 2>&1; "

use rule create_local_links as create_links_ragtag_scaffolds with:
    input:
        final_fasta=rules.ragtag_generate_filtered_fasta.output.ragtag_filtered_fasta,
        final_agp=rules.ragtag_filter.output.ragtag_filtered_agp,
        final_stats=rules.ragtag_filter.output.ragtag_filtered_stats,
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        final_fasta=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{genome_prefix}.ref_scaffolding.{haplotype, [^/]+}.fasta",
        final_agp=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{genome_prefix}.ref_scaffolding.{haplotype, [^/]+}.agp",
        final_stats=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{genome_prefix}.ref_scaffolding.{haplotype, [^/]+}.stats",
    log:
        ln=config["out_dir"] / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.ln.log",
