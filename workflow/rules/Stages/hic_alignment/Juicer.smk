localrules: create_fastq_links_for_juicer
#ruleorder: create_fastq_links_for_juicer > create_final_links

use rule create_local_links as create_fastq_links_for_juicer with:
    input:
        forward_fastq=lambda wildcards: config["out_dir"] / "data//hic/final/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               config["data"]["hic"]["conv_fwd_sfx"],
                                                                                               config["data"]["hic"]["conv_ext"]) if wildcards.phasing_kmer_length == "NA" else \
                                config["out_dir"] / "{0}/{1}/{2}.{0}.{3}/reads/hic/{4}/{5}{6}{7}".format(config["phasing_stage"],
                                                                                                    detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."),
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.haplotype,
                                                                                                    wildcards.phasing_kmer_length,
                                                                                                    wildcards.pairprefix,
                                                                                                    config["data"]["hic"]["conv_fwd_sfx"],
                                                                                                    config["data"]["hic"]["conv_ext"]),
        reverse_fastq=lambda wildcards: config["out_dir"] / "data//hic/final/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               config["data"]["hic"]["conv_rev_sfx"],
                                                                                               config["data"]["hic"]["conv_ext"]) if wildcards.phasing_kmer_length == "NA" else \
                                config["out_dir"] / "{0}/{1}/{2}.{0}.{3}/reads/hic/{4}/{5}{6}{7}".format(config["phasing_stage"],
                                                                                                    detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."),
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.haplotype,
                                                                                                    wildcards.phasing_kmer_length,
                                                                                                    wildcards.pairprefix,
                                                                                                    config["data"]["hic"]["conv_rev_sfx"],
                                                                                                    config["data"]["hic"]["conv_ext"]),
        log_dir=ancient(config["out_dir"] / "log/")

    output:
        forward_fastq=config["out_dir"] / ("hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/fastq/{pairprefix}_R1_001%s" % config["data"]["hic"]["conv_ext"]),
        reverse_fastq=config["out_dir"] / ("hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/fastq/{pairprefix}_R2_001%s" % config["data"]["hic"]["conv_ext"]),
    log:
        ln=config["out_dir"] / "log/create_fastq_links_for_juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.ln.log",

rule juicer: #
    input:
        fasta=config["out_dir"]  / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta",
        index=config["out_dir"]  / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta.ann",
        restriction_site_file=config["out_dir"]  / ("hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}_%s.txt" % config["hic_enzyme_set"]),
        forward_fastq=expand(config["out_dir"] / ("hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/fastq/{pairprefix}_R1_001%s" % config["data"]["hic"]["conv_ext"]),
                             pairprefix=config["data"]["hic"]["pair_prefix_list"],
                             allow_missing=True),
        reverse_fastq=expand(config["out_dir"] / ("hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/fastq/{pairprefix}_R2_001%s" % config["data"]["hic"]["conv_ext"]),
                             pairprefix=config["data"]["hic"]["pair_prefix_list"],
                             allow_missing=True),
    params:
        restriction_seq=config["hic_enzyme_set"]  if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else "none",
        fastq_extensions=config["fastq_ext"]
    output:
        aln_dir=directory(config["out_dir"] / "hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/aligned/"),
        merged_no_dups=config["out_dir"] / "hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.merged_nodups.txt",
        merged_dedup_bam=config["out_dir"] / "hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.rmdup.bam",
        merged_inter_30=config["out_dir"] / "hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.inter_30.txt",
        merged_inter=config["out_dir"] / "hic_alignment/{parameters, [^/]*juicer[^/]*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.inter.txt",
    log:
        juicer=config["out_dir"] / "log/juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.juicer.log",
        mkdir=config["out_dir"] / "log/juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.mkdir.log",
        ln=config["out_dir"] / "log/juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.ln.log",
        cluster_log=config["out_dir"] / "log/juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "log/juicer.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("juicer"),
        cpus=parameters["threads"]["juicer"] ,
        time=parameters["time"]["juicer"],
        mem=parameters["memory_mb"]["juicer"]
    threads: parameters["threads"]["juicer"]

    shell:
        " OUTPUT_DIR=`dirname {output.merged_no_dups}`; "
        " OUTPUT_DIR=`realpath -s ${{OUTPUT_DIR}}`; "
        " > {log.ln}; "
        " SCRIPT=`realpath -s ./workflow/external_tools/juicer/scripts/juicer.sh`;"
        " JUICER_DIR=`realpath -s ./workflow/external_tools/juicer/`; "
        " FASTA=`realpath -s {input.fasta}`; "
        " if [ \"{params.restriction_seq}\" == \"none\" ]; "
        " then "
        "     RESTRICTION_SITE_OPTION=''; "
        " else "
        "     RESTRICTION_SITE_OPTION=\" -y `realpath -s {input.restriction_site_file}` \"; "
        " fi; "
        " ${{SCRIPT}} -t {threads} -T {threads} -D ${{JUICER_DIR}} -g {wildcards.genome_prefix} -s {params.restriction_seq} "
        "     -z ${{FASTA}} ${{RESTRICTION_SITE_OPTION}} -S early --assembly -d ${{OUTPUT_DIR}} > {log.juicer} 2>&1; "
        " mv ${{OUTPUT_DIR}}/aligned/merged_nodups.txt {output.merged_no_dups}; "
        " mv ${{OUTPUT_DIR}}/aligned/merged_dedup.bam {output.merged_dedup_bam}; " 
        " mv ${{OUTPUT_DIR}}/aligned/inter_30.txt {output.merged_inter_30}; "
        " mv ${{OUTPUT_DIR}}/aligned/inter.txt {output.merged_inter}; " 
        " rm -r ${{OUTPUT_DIR}}/splits ; "
