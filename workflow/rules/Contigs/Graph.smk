rule gfa2fasta:
    input:
        gfa=output_dict["contig"] / "{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.gfa"
    output:
        fasta=output_dict["contig"] / "{parameters}/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.fasta"
    log:
        std=output_dict["log"] / "gfa2fasta.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "gfa2fasta.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gfa2fasta.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gfa2fasta.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
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
        gfa=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.{assembly_stage}.{haplotype}.gfa"
    output:
        cov=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.gfa.cov",
        len_cov=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.gfa.lencov"
    log:
        std=output_dict["log"] / "gfa2fasta.{assembly_stage}.hifiasm_{contig_options}.{genome_prefix}.{haplotype}.log",
        cut=output_dict["log"] / "gfa2fasta.{assembly_stage}.hifiasm_{contig_options}.{genome_prefix}.{haplotype}.cut.log",
        cluster_log=output_dict["cluster_log"] / "gfa2fasta.{assembly_stage}.hifiasm_{contig_options}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gfa2fasta.{assembly_stage}.hifiasm_{contig_options}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gfa2fasta.{assembly_stage}.hifiasm_{contig_options}.{genome_prefix}.{haplotype}.benchmark.txt"
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
         " cut -f 1,3 {output.len_cov} > {output.cov} 2>{log.cut};"

