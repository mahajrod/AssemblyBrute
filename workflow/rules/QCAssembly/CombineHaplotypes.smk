rule combine_haplotypes:
    input:
        fasta_list=lambda wildcards: expand(config["out_dir"] / ("%s/%s/%s.%s.{haplotype}.fasta" % (wildcards.assembly_stage,
                                                                                                     wildcards.parameters,
                                                                                                     wildcards.genome_prefix,
                                                                                                     wildcards.assembly_stage)),
                                             haplotype=stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"] ,
                                             allow_missing=True)
    output:
        combined_fasta=config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.combined.fasta",
    params:
        haplotype_list=lambda wildcards: stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
        out_dir=str(config["out_dir"])
    log:
        log=config["out_dir"] / "log/combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/combine_haplotypes.{assembly_stage}..{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("combine_haplotypes"),
        cpus=parameters["threads"]["combine_haplotypes"],
        time=parameters["time"]["combine_haplotypes"],
        mem=parameters["memory_mb"]["combine_haplotypes"]
    threads: parameters["threads"]["combine_haplotypes"]

    shell:
        " > {output.combined_fasta}; "
        " > {log.log}; "
        " for HAP in {params.haplotype_list};"
        " do "
        "    label_sequences.py -i {params.out_dir}/{wildcards.assembly_stage}/{wildcards.parameters}/{wildcards.genome_prefix}.{wildcards.assembly_stage}.${{HAP}}.fasta "
        "        -l ${{HAP}} -s '.' >> {output.combined_fasta} 2>{log.log} ; "
        " done "