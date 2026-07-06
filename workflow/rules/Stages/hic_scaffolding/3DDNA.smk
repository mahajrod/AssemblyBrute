localrules: threeddna_create_links


rule threeddna: #
    input:
        reference=lambda wildcards: config["out_dir"]  / ("hic_alignment/%s/%s.%s.%s.fasta" % (get_prev_stage_parameters(wildcards.parameters),
                                                                                               wildcards.genome_prefix,
                                                                                               stage_dict["hic_scaffolding"].prev_stage,
                                                                                               wildcards.haplotype)),
        merged_nodups=lambda wildcards: expand(config["out_dir"] / "hic_alignment/{prev_stage_parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.merged_nodups.txt",
                                               phasing_kmer_length=[stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"]["phasing_kmer_length"]],
                                               prev_stage_parameters=[get_prev_stage_parameters(wildcards.parameters),],
                                               allow_missing=True),
        bam_general_stats=lambda wildcards: expand(config["out_dir"] / "hic_alignment/{prev_stage_parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.rmdup.bam.general_stats",
                                               phasing_kmer_length=[stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"]["phasing_kmer_length"]],
                                               prev_stage_parameters=[get_prev_stage_parameters(wildcards.parameters),],
                                               allow_missing=True)

    params:
        restriction_seq=config["hic_enzyme_set"]  if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else "none",
        fastq_extensions=config["fastq_ext"],
        editor_repeat_coverage=lambda wildcards: parse_option("editor_repeat_coverage",
                                                              stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"], " --editor-repeat-coverage "),
        min_contig_length=lambda wildcards: parse_option("min_contig_len",
                                                         stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"], " --input "),
        min_mapping_quality=lambda wildcards: parse_option("min_mapping_quality",
                                                           stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"], " --mapq "),
    output:
        draft_fasta=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}/scaffolding/{genome_prefix}.draft.{haplotype}.fasta",
        rawchrom_hic=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}/scaffolding/{genome_prefix}.draft.{haplotype}.rawchrom.hic",
        rawchrom_assembly=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}/scaffolding/{genome_prefix}.draft.{haplotype}.rawchrom.assembly",
        hic_fasta=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}/scaffolding/{genome_prefix}.draft.{haplotype}_HiC.fasta",
    log:
        threeddna=config["out_dir"] / "log/threeddna.{parameters}.{genome_prefix}.{haplotype}.threeddna.log",
        ln=config["out_dir"] / "log/threeddna.{parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=config["out_dir"] / "log/threeddna.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/threeddna.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/threeddna.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("threeddna"),
        cpus=parameters["threads"]["threeddna"] ,
        time=parameters["time"]["threeddna"],
        mem=parameters["memory_mb"]["threeddna"]
    threads: parameters["threads"]["threeddna"]

    shell:
        " OUTPUT_DIR=`dirname {output.draft_fasta}`; "
        " SCRIPT=`realpath -s workflow/external_tools/3d-dna/run-3ddna-pipeline.sh`; "
        " INPUT_FASTA=`realpath -s {input.reference}`; "
        " MERGED_NODUPS=`realpath -s {input.merged_nodups}`; "
        " THREEDDNA_LOG=`realpath -s  {log.threeddna}`; "
        " LN_LOG=`realpath -s {log.ln}`; "
        " > ${{LN_LOG}}; "
        " ln -s ${{INPUT_FASTA}} {output.draft_fasta} >> ${{LN_LOG}} 2>&1; "
        " ln -s ${{MERGED_NODUPS}} ${{OUTPUT_DIR}}/`basename {input.merged_nodups}` >> ${{LN_LOG}} 2>&1; "
        " cd ${{OUTPUT_DIR}}; "
        " ${{SCRIPT}} {params.editor_repeat_coverage} {params.min_contig_length} {params.min_mapping_quality} `basename {output.draft_fasta}` `basename {input.merged_nodups}` > ${{THREEDDNA_LOG}} 2>&1; "

localrules: threeddna_create_links

use rule create_local_links as threeddna_create_links with: #
    input:
        rawchrom_hic=config["out_dir"] / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.draft.{haplotype}.rawchrom.hic",
        rawchrom_assembly=config["out_dir"] / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.draft.{haplotype}.rawchrom.assembly",
        hic_fasta=config["out_dir"] / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.draft.{haplotype}_HiC.fasta",
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        rawchrom_hic=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}.rawchrom.hic",
        rawchrom_assembly=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}.rawchrom.assembly",
        hic_fasta=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*threeddna[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^./]*}.fasta",
    log:
        ln=config["out_dir"] / "log/threeddna_create_links.{parameters}.{genome_prefix}.{haplotype}.ln.log",
