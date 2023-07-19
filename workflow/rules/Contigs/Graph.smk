rule gfa2fasta:
    input:
        gfa="{gra_prefix}.gfa"
    output:
        fasta="{gra_prefix}.fasta"
    log:
        std="{gra_prefix}.gfa2fasta.log",
        cluster_log="{gra_prefix}.gfa2fasta.cluster.log",
        cluster_err="{gra_prefix}.gfa2fasta.cluster.err"
    benchmark:
        "{gra_prefix}.gfa2fasta.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["gfa2fasta"],
        time=parameters["time"]["gfa2fasta"],
        mem=parameters["memory_mb"]["gfa2fasta"],
    threads:
        parameters["threads"]["gfa2fasta"]
    shell:
         " gfatools gfa2fa {input.gfa} > {output.fasta} 2>{log.std};"

rule get_length_and_coverage_from_hifiasm_graph:
    input:
        gfa="{gra_prefix}.gfa"
    output:
        cov="{gra_prefix}.gfa.cov",
        len_cov="{gra_prefix}.gfa.lencov"
    log:
        std="{gra_prefix}.get_length_and_coverage_from_hifiasm_graph.log",
        cluster_log="{gra_prefix}.get_length_and_coverage_from_hifiasm_graph.cluster.log",
        cluster_err="{gra_prefix}.get_length_and_coverage_from_hifiasm_graph.cluster.err"
    benchmark:
        "{gra_prefix}.get_length_and_coverage_from_hifiasm_graph.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["get_coverage_from_hifiasm_graph"],
        time=parameters["time"]["get_coverage_from_hifiasm_graph"],
        mem=parameters["memory_mb"]["get_coverage_from_hifiasm_graph"],
    threads:
        parameters["threads"]["get_coverage_from_hifiasm_graph"]
    shell:
         " workflow/scripts/extract_length_and_coverage_from_hifiasm_gfa.bash {input.gfa} > {output.len_cov} 2>{log.std}; "
         " cut -f 1,3 {output.len_cov} > {output.cov} 2>>{log.std};"

