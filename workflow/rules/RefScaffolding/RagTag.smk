localrules: create_links_ragtag_scaffolds

rule ragtag: #
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/%s.%s.{haplotype}.fasta" % (config["genome_prefix"],
                                                                                      stage_dict["ref_scaffolding"]["prev_stage"],
                                                                                      stage_dict["ref_scaffolding"]["prev_stage"])),
        contig_both_blacklist=out_dir_path / ("%s/{prev_stage_parameters}/assembly_qc/telomere/%s.%s.{haplotype}/%s.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.both.ids" %
                                                  (stage_dict["ref_scaffolding"]["prev_stage"],
                                                   config["genome_prefix"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   config["genome_prefix"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   "canonical" if config["use_canonical_telomere"] else "non_canonical"
                                                   )),
        #contig_five_prime_blacklist=out_dir_path / ("%s/{prev_stage_parameters}/assembly_qc/telomere/{genome_prefix}.%s.{haplotype}/{genome_prefix}.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids" %
        #                                          (stage_dict["ref_scaffolding"]["prev_stage"],
        #                                           stage_dict["ref_scaffolding"]["prev_stage"],
        #                                           stage_dict["ref_scaffolding"]["prev_stage"],
        #                                           "canonical" if config["use_canonical_telomere"] else "non_canonical"
        #                                           )),
        #contig_three_prime_blacklist=out_dir_path / ("%s/{prev_stage_parameters}/assembly_qc/telomere/{genome_prefix}.%s.{haplotype}/{genome_prefix}.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids" %
        #                                          (stage_dict["ref_scaffolding"]["prev_stage"],
        #                                           stage_dict["ref_scaffolding"]["prev_stage"],
        #                                           stage_dict["ref_scaffolding"]["prev_stage"],
        #                                           "canonical" if config["use_canonical_telomere"] else "non_canonical"
        #                                           )),
        reference_fasta=out_dir_path / "data/reference/{reference}/repeats/{reference}.softmasked.fasta",
    output:
        ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^./]+}/ragtag.scaffold.fasta",
        ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^./]+}/ragtag.scaffold.agp",
        #ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^./]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.stats",
    params:
        min_aln_len=lambda wildcards: parse_option("min_aln_len", parameters["tool_options"]["ragtag"][wildcards.ref_scaf_parameters], " -f ")

    log:
        ragtag=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{haplotype}.log",
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{haplotype}.ln.log",
        awk=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{haplotype}.awk.log",
        cluster_log=output_dict["cluster_log"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["ragtag"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["ragtag"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("ragtag"),
        cpus=parameters["threads"]["ragtag"],
        time=parameters["time"]["ragtag"],
        mem=parameters["memory_mb"]["ragtag"]
    threads: parameters["threads"]["ragtag"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_fasta}`; "
        " echo 'Scaffolding...' > {log.ragtag} 2>&1;"
        " ragtag.py scaffold -t {threads} -j {input.contig_both_blacklist} {params.min_aln_len} "
        " -w {input.reference_fasta} {input.fasta} -o ${{RAGTAG_DIR}}  >> {log.ragtag} 2>&1; "


rule ragtag_filter: #
    input:
        ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.scaffold.agp",
        contig_five_prime_blacklist=out_dir_path / ("%s/{prev_stage_parameters}/assembly_qc/telomere/{genome_prefix}.%s.{haplotype}/{genome_prefix}.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids" %
                                                  (stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   "canonical" if config["use_canonical_telomere"] else "non_canonical"
                                                   )),
        contig_three_prime_blacklist=out_dir_path / ("%s/{prev_stage_parameters}/assembly_qc/telomere/{genome_prefix}.%s.{haplotype}/{genome_prefix}.%s.{haplotype}.%s_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids" %
                                                  (stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   stage_dict["ref_scaffolding"]["prev_stage"],
                                                   "canonical" if config["use_canonical_telomere"] else "non_canonical"
                                                   )),
        reference_fasta=out_dir_path / "data/reference/{reference}/repeats/{reference}.softmasked.fasta",
    output:
        ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^./]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.fasta",
        ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^./]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.agp",
        ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^./]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.stats",
    params:
        min_aln_len=lambda wildcards: parse_option("min_aln_len", parameters["tool_options"]["ragtag"][wildcards.ref_scaf_parameters], " -f ")

    log:
        ragtag=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag_filter.{genome_prefix}.ragtag.{haplotype}.log",
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag_filter.{genome_prefix}.ragtag.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ragtag_filter.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("ragtag"),
        cpus=parameters["threads"]["ragtag"],
        time=parameters["time"]["ragtag"],
        mem=parameters["memory_mb"]["ragtag"]
    threads: parameters["threads"]["ragtag"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_fasta}`; "
        " echo 'Spliting agp at internal telomeres...' >> {log.ragtag} 2>&1; "
        " workflow/scripts/curation/split_agp_by_telomeric_merges.py -a {input.ragtag_agp} "
        "   --five_prime_blacklist {input.contig_five_prime_blacklist} "
        "   --three_prime_blacklist {input.contig_three_prime_blacklist} -p ${{RAGTAG_DIR}}/ragtag.scaffold >> {log.ragtag} 2>&1; "
        " echo 'Assembly fasta from modified agp...' >> {log.ragtag} 2>&1; "
        " ragtag.py agp2fa ${{RAGTAG_DIR}}/ragtag.scaffold.telomere.split.final.renamed.agp {input.reference_fasta} > ${{RAGTAG_DIR}}/ragtag.scaffold.telomere.split.final.renamed.fasta 2>>{log.ragtag}; "
        " ln -sf ragtag.scaffold.telomere.split.final.renamed.fasta {output.ragtag_fasta} > {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.telomere.split.final.renamed.agp {output.ragtag_agp} >> {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.telomere.split.final.renamed.stats {output.ragtag_stats} >> {log.ln} 2>&1; "

rule create_links_ragtag_scaffolds:
    input:
        fasta=rules.ragtag_filter.output.ragtag_fasta,
        agp=rules.ragtag_filter.output.ragtag_agp,
        stats=rules.ragtag_filter.output.ragtag_stats,
        #reference_syn=out_dir_path / "data/reference/{reference}/{reference}.syn",
    output:
        final_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{genome_prefix, [^/]+}.ref_scaffolding.{haplotype, [^/]+}.fasta",
        final_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{genome_prefix, [^/]+}.ref_scaffolding.{haplotype, [^/]+}.agp",
        final_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{genome_prefix, [^/]+}.ref_scaffolding.{haplotype, [^/]+}.stats",
        #ragtag_syn=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{haplotype, [^/]+}/{genome_prefix}.ref_scaffolding.{haplotype}.syn",
    log:
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.ln.log",
        #awk=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.awk.log",
        #rename_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.rename_fasta.log",
        #rename_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.rename_agp.log",
        cluster_log=output_dict["cluster_log"] / "rename_ragtag_scaffolds.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "rename_ragtag_scaffolds{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "rename_ragtag_scaffolds.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("rename_ragtag_scaffolds"),
        cpus=parameters["threads"]["rename_ragtag_scaffolds"],
        time=parameters["time"]["rename_ragtag_scaffolds"],
        mem=parameters["memory_mb"]["rename_ragtag_scaffolds"]
    threads: parameters["threads"]["rename_ragtag_scaffolds"]

    shell:
        " > {log.ln}; "
        " for FILENAME in {output.final_fasta} {output.final_agp} {output.final_stats}; "
        "   do"
        "   ln -sf {wildcards.haplotype}/`basename ${{FILENAME}}` ${{FILENAME}} >> {log.ln} 2>&1; "
        "   done "
        #" RAGTAG_DIR=`dirname {input.fasta}`; "
        #" awk '{{print $2\"_RagTag\tps\"$1 }}' {input.reference_syn} > {output.ragtag_syn} 2>{log.awk}; "
        #" rename_sequence_ids.py -i {input.fasta} -o {output.final_fasta} -l -s {output.ragtag_syn} -k 0 -c 1 > {log.rename_fasta} 2>&1; "
        #" ln -sf {wildcards.haplotype}/{wildcards.genome_prefix}.ref_scaffolding.{wildcards.haplotype}.stats {output.final_stats} > {log.ln} 2>&1; "
        #" replace_column_value_by_syn.py -i {input.agp} -o {output.final_agp} -s {output.ragtag_syn} -c 0 > {log.rename_agp} 2>&1; " # -k 0 -c 1