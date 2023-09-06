rule quast:
    input:
        #primary_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.hifi.hic.p_ctg.fasta" % config["genome_name"]),
        #alternative_assembly=out_dir_path / ("{assembly_stage}/{assembler}/%s.{assembly_stage}.{assembler}.hifi.hic.a_ctg.fasta" % config["genome_name"])
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
    output:
        #summary=output_dict["assembly_qc"] /("{assembly_stage}/{assembler}/{assembly_stage}.{assembler}.hifi.hic.{haplotype}_ctg.gfa" % config["genome_name"]),
        dir=directory(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/quast/{genome_prefix}.{assembly_stage}.{haplotype}"),

    params:
        large_genome_flag="--large" if parameters["tool_options"]["quast"]["large_genome"] else "",
    log:
        std=output_dict["log"] / "quast.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "quast.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "quast.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "quast.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("quast"),
        cpus=parameters["threads"]["quast"],
        time=parameters["time"]["quast"],
        mem=parameters["memory_mb"]["quast"],
    threads:
        parameters["threads"]["quast"]
    shell:
         " quast -o {output.dir} -t {threads} -l {wildcards.haplotype} {params.large_genome_flag} "
         " {input.assembly}  1>{log.std} 2>&1;"
