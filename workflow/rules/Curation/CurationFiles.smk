localrules: gather_curation_files, gather_curation_tracks

rule gather_curation_files: #
    input:
        hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
        assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.assembly",
    output:
        hic=out_dir_path / "curation_files/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.hic",
        assembly=out_dir_path / "curation_files/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.assembly",
    log:
        cp=output_dict["log"]  / "gather_curation_files.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "gather_curation_files.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gather_curation_files.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "gather_curation_files.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("gather_curation_files"),
        cpus=parameters["threads"]["gather_curation_files"],
        time=parameters["time"]["gather_curation_files"],
        mem=parameters["memory_mb"]["gather_curation_files"]
    threads: parameters["threads"]["gather_curation_files"]

    shell:
        " cp {input.hic} `dirname {output.hic}` > {log.cp} 2>&1; "
        " cp {input.assembly} `dirname {output.assembly}` >> {log.cp} 2>&1; "

rule gather_curation_tracks: #
    input:
        contig_dir=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/contigs/",
        scaffolds_dir=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/scaffolds/",
    output:
        contig_dir=directory(out_dir_path / "curation_files/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^/]+}/contigs/"),
        scaffolds_dir=directory(out_dir_path / "curation_files/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^/]+}/scaffolds/"),
    log:
        cp=output_dict["log"]  / "gather_curation_tracks.{prev_stage_parameters}..{curation_parameters}.{haplotype}.cp.log",
        mkdir=output_dict["log"]  / "gather_curation_tracks.{prev_stage_parameters}..{curation_parameters}.{haplotype}.mkdir.log",
        cluster_log=output_dict["cluster_log"] / "gather_curation_tracks.{prev_stage_parameters}..{curation_parameters}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gather_curation_tracks.{prev_stage_parameters}..{curation_parameters}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "gather_curation_tracks.{prev_stage_parameters}..{curation_parameters}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("gather_curation_files"),
        cpus=parameters["threads"]["gather_curation_files"],
        time=parameters["time"]["gather_curation_files"],
        mem=parameters["memory_mb"]["gather_curation_files"]
    threads: parameters["threads"]["gather_curation_files"]

    shell:
        " mkdir -p {output.contig_dir} {output.scaffolds_dir} > {log.mkdir} 2>&1 ; "
        " cp {input.contig_dir}/*.bedgraph {output.contig_dir} > {log.cp} 2>&1; "
        " cp {input.scaffolds_dir}/*.png {input.scaffolds_dir}/*.svg {input.scaffolds_dir}/*.bedgraph {input.scaffolds_dir}/*.tab.gz {output.scaffolds_dir} >> {log.cp} 2>&1; "