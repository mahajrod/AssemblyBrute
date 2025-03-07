
rule miniprot:
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        miniprot_gff=out_dir_path / "{assembly_stage}/{parameters, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.miniprot.gff",
        candidate_tsv=out_dir_path / "{assembly_stage}/{parameters, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.candidates.microchromosomes.tsv"
    params:
        microchromosome_prot_set=config["microchromosome_prot_set"]
    log:
        miniprot=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.miniprot.log",
        awk=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.awk.log",
        grep=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.grep.log",
        cut=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cut.log",
        tr=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.tr.log",
        cut2=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cut2.log",
        sed=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.sed.log",
        awk2=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.awk2.log",
        cut3=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cut3.log",
        sort=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.sort.log",
        uniq=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.uniq.log",
        sort2=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.sort2.log",
        awk3=output_dict["log"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.awk3.log",
        cluster_log=output_dict["cluster_log"] / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "miniprot.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["microchromosomes"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["microchromosomes"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("miniprot_microchromosomes"),
        cpus=parameters["threads"]["miniprot_microchromosomes"],
        time=parameters["time"]["miniprot_microchromosomes"],
        mem=parameters["memory_mb"]["miniprot_microchromosomes"]
    threads: parameters["threads"]["miniprot_microchromosomes"]
    shell:
        " miniprot -t {threads} --gff {input.fasta} {params.microchromosome_prot_set} > {output.miniprot_gff} 2>{log.miniprot}; "
        " awk '$3== \"mRNA\"' {output.miniprot_gff} 2>{log.awk} | grep -w \"Rank=1\" 2>{log.grep} | "
        " cut -f1,9 2>{log.cut} | tr \";\" \"\\t\" 2>{log.tr} | cut -f1,4,5,6,7  2>{log.cut2} | "
        " sed 's/Identity=//g;s/Positive=//g' 2>{log.sed} | awk '$2 >= 0.7' 2>{log.awk2} |  cut -f1 2>{log.cut3} | "
        " sort 2>{log.sort} | uniq -c 2>{log.uniq} | sort -k1,1nr 2>{log.sort2} | awk '{{print $2 \"\\t\" $1}}' > {output.candidate_tsv} 2>{log.awk3} "


rule place_microsomes_first:
    input:
        candidate_tsv=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.candidates.microchromosomes.tsv",
        len_file=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        assembly= lambda wildcards: (out_dir_path / "{0}/{1}/{2}.{0}.{3}.assembly".format(wildcards.assembly_stage,
                                                                                                wildcards.parameters,
                                                                                                wildcards.genome_prefix,
                                                                                                wildcards.haplotype,
                                                                                                )) if (wildcards.assembly_stage == "hic_scaffolding") and (wildcards.haplotype != "combined") else []
    output:
        filtered_tsv=out_dir_path / "{assembly_stage}/{parameters, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.max{max_length, [^/]+}.candidates.microchromosomes.filtered.tsv",
        reordered_fasta=out_dir_path / "{assembly_stage}/{parameters, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.max{max_length, [^/]+}.reordered.fasta",
    params:
        assembly_option= lambda wildcards: " -a " + str(out_dir_path / "{0}/{1}/{2}.{0}.{3}.assembly".format(wildcards.assembly_stage,
                                                                                                            wildcards.parameters,
                                                                                                            wildcards.genome_prefix,
                                                                                                            wildcards.haplotype,
                                                                                                            )) if (wildcards.assembly_stage == "hic_scaffolding") and (wildcards.haplotype != "combined") else ""
    log:
        log=output_dict["log"]  / "place_microsomes_first.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.max{max_length}.log",
        cluster_log=output_dict["cluster_log"] / "place_microsomes_first.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.max{max_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "place_microsomes_first.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.max{max_length}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "place_microsomes_first.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{max_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("place_microsomes_first"),
        cpus=parameters["threads"]["place_microsomes_first"],
        time=parameters["time"]["place_microsomes_first"],
        mem=parameters["memory_mb"]["place_microsomes_first"]
    threads: parameters["threads"]["place_microsomes_first"]
    shell:
        " OUTPUT_PREFIX={output.reordered_fasta}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.reordered.fasta}}; "
        " workflow/scripts/curation/move_microchromosomes_first.py  -i {input.candidate_tsv}  {params.assembly_option} "
        " -l {input.len_file} -f {input.fasta} -m {wildcards.max_length} -o ${{OUTPUT_PREFIX}} > {log.log} 2>&1; "
