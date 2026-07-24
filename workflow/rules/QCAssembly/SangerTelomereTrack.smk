
ruleorder: copy_telomere_files > classify_telomeric_regions_windows # TODO: check if it is still necessary

localrules: copy_telomere_files
localrules: copy_telomere_track_for_pretext, create_telomere_track_for_pretext

rule create_telomere_track_for_pretext:
    input:
        canonical_telo="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical.telomere",
        non_canonical_telo="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical.telomere",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        canonical_telo_bedgraph="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical.telomere.pretext.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical.telomere.pretext.bedgraph"
    log:
        canonical="{fasta_dir}/log/create_telomere_track_for_pretext.{fasta_prefix}.canonical.log",
        non_canonical="{fasta_dir}/log/create_telomere_track_for_pretext.{fasta_prefix}.non_canonical.log",
        cluster_log="{fasta_dir}/log/create_telomere_track_for_pretext.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/create_telomere_track_for_pretext.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/create_telomere_track_for_pretext.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_telomere_track_for_pretext"),
        cpus=parameters["threads"]["create_telomere_track_for_pretext"] ,
        time=parameters["time"]["create_telomere_track_for_pretext"],
        mem=parameters["memory_mb"]["create_telomere_track_for_pretext"],
    threads: parameters["threads"]["create_telomere_track_for_pretext"]

    shell:
        " awk '{{print $1\"\\t\"$4\"\\t\"$5\"\\t\"((($5-$4)<0)?-($5-$4):($5-$4))}}' {input.canonical_telo} | "
        "     sed 's/>//g' > {output.canonical_telo_bedgraph} 2>{log.canonical}; "
        " awk '{{print $1\"\\t\"$4\"\\t\"$5\"\\t\"((($5-$4)<0)?-($5-$4):($5-$4))}}' {input.non_canonical_telo} | "
        "     sed 's/>//g' > {output.non_canonical_telo_bedgraph} 2>{log.non_canonical}; "


rule copy_telomere_track_for_pretext:
    input:
        canonical_telo_bedgraph="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical.telomere.pretext.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical.telomere.pretext.bedgraph",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        canonical_telo_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.canonical.telomere.pretext.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.non_canonical.telomere.pretext.bedgraph",
    log:
        canonical="{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.canonical.log",
        non_canonical="{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.non_canonical.log",
        cluster_log="{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("telo_container"),
        cpus=parameters["threads"]["copy_telomere_track_for_pretext"] ,
        time=parameters["time"]["copy_telomere_track_for_pretext"],
        mem=parameters["memory_mb"]["copy_telomere_track_for_pretext"],
    threads: parameters["threads"]["copy_telomere_track_for_pretext"]

    shell:
        " cp -f {input.canonical_telo_bedgraph} {output.canonical_telo_bedgraph} > {log.canonical} 2>&1; "
        " cp -f {input.non_canonical_telo_bedgraph} {output.non_canonical_telo_bedgraph} > {log.non_canonical} 2>&1; "

rule copy_telomere_files:
    input:
        fai="{fasta_dir}/{fasta_prefix}.fasta.fai",
        canonical_telo_track="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.bedgraph",
        canonical_telo_warning_track="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere_warning.win1000.step200.track.bedgraph",
        non_canonical_telo_track="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.bedgraph",
        non_canonical_telo_warning_track="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere_warning.win1000.step200.track.bedgraph",
        canonical_region_all_status="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.all.status",
        canonical_region_filtered_status="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        canonical_region_filtered_count="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        canonical_region_filtered_scaffold_status="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        non_canonical_region_all_status="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.all.status",
        non_canonical_region_filtered_status="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        non_canonical_region_filtered_count="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        non_canonical_region_filtered_scaffold_status="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        canonical_region_filtered_scaffold_both_telomeres_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.both.ids",
        canonical_region_filtered_scaffold_five_prime_only_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime_only.ids",
        canonical_region_filtered_scaffold_three_prime_only_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime_only.ids",
        non_canonical_region_filtered_scaffold_both_telomeres_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.both.ids",
        non_canonical_region_filtered_scaffold_five_prime_only_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime_only.ids",
        non_canonical_region_filtered_scaffold_three_prime_only_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime_only.ids",
        canonical_region_filtered_scaffold_five_prime_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids",
        canonical_region_filtered_scaffold_three_prime_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids",
        non_canonical_region_filtered_scaffold_five_prime_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids",
        non_canonical_region_filtered_scaffold_three_prime_id_file="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        canonical_telo_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.bedgraph",
        canonical_telo_warning_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.canonical_telomere_warning.win1000.step200.track.bedgraph",
        non_canonical_telo_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.bedgraph",
        non_canonical_telo_warning_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere_warning.win1000.step200.track.bedgraph",
        canonical_region_all_status="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.all.status",
        canonical_region_filtered_status="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        canonical_region_filtered_count="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        canonical_region_filtered_scaffold_status="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        non_canonical_region_all_status="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.all.status",
        non_canonical_region_filtered_status="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.status",
        non_canonical_region_filtered_count="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.count",
        non_canonical_region_filtered_scaffold_status="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.status",
        canonical_region_filtered_scaffold_both_telomeres_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.both.ids",
        canonical_region_filtered_scaffold_five_prime_only_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime_only.ids",
        canonical_region_filtered_scaffold_three_prime_only_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime_only.ids",
        non_canonical_region_filtered_scaffold_both_telomeres_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.both.ids",
        non_canonical_region_filtered_scaffold_five_prime_only_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime_only.ids",
        non_canonical_region_filtered_scaffold_three_prime_only_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime_only.ids",
        canonical_region_filtered_scaffold_five_prime_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids",
        canonical_region_filtered_scaffold_three_prime_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids",
        non_canonical_region_filtered_scaffold_five_prime_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.five_prime.ids",
        non_canonical_region_filtered_scaffold_three_prime_id_file="{fasta_dir}/assembly_qc/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical_telomere.win1000.step200.track.collapsed.filtered.scaffold.telomeres.three_prime.ids",
    log:
        std="{fasta_dir}/log/create_telomere_links.{fasta_prefix}.std.log",
        cluster_log="{fasta_dir}/log/create_telomere_links.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/create_telomere_links.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/create_telomere_links.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_telomere_warning"),
        cpus=parameters["threads"]["get_telomere_warning"] ,
        time=parameters["time"]["get_telomere_warning"],
        mem=parameters["memory_mb"]["get_telomere_warning"]
    threads: parameters["threads"]["get_telomere_warning"]

    shell:
        " TRACK_DIR=`dirname {output.canonical_telo_track}`; "
        " TELOMERE_DIR=`dirname {output.canonical_region_all_status}`; "
        " cp {input.canonical_telo_track} {input.canonical_telo_warning_track} {input.non_canonical_telo_track} "
        "     {input.non_canonical_telo_warning_track} ${{TRACK_DIR}} > {log.std} 2>&1; "
        " cp {input.canonical_region_all_status} {input.canonical_region_filtered_status} "
        "     {input.canonical_region_filtered_count} {input.canonical_region_filtered_scaffold_status} ${{TELOMERE_DIR}} > {log.std} 2>&1; "
        " cp {input.non_canonical_region_all_status} {input.non_canonical_region_filtered_status} "
        "     {input.non_canonical_region_filtered_count} {input.non_canonical_region_filtered_scaffold_status} ${{TELOMERE_DIR}} > {log.std} 2>&1; "
        " cp {input.canonical_region_filtered_scaffold_both_telomeres_id_file} {input.canonical_region_filtered_scaffold_five_prime_only_id_file} "
        "     {input.canonical_region_filtered_scaffold_three_prime_only_id_file} {input.canonical_region_filtered_scaffold_five_prime_id_file}"
        "     {input.canonical_region_filtered_scaffold_three_prime_id_file} ${{TELOMERE_DIR}} > {log.std} 2>&1; "
        " cp {input.non_canonical_region_filtered_scaffold_both_telomeres_id_file} {input.non_canonical_region_filtered_scaffold_five_prime_only_id_file} "
        "     {input.non_canonical_region_filtered_scaffold_three_prime_only_id_file} {input.non_canonical_region_filtered_scaffold_five_prime_id_file}"
        "     {input.non_canonical_region_filtered_scaffold_three_prime_id_file} ${{TELOMERE_DIR}} > {log.std} 2>&1; "
