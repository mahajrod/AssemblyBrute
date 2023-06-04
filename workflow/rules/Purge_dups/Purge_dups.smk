
localrules: create_primary_contig_link, merge_pri_hapdups_with_alt, extract_stats_from_purge_dups_file #create_link_for_purged_fasta,
localrules: merge_pri_hapdups_with_alt_for_len_files # create_primary_contig_len_file_link,
ruleorder: create_primary_contig_link > merge_pri_hapdups_with_alt
ruleorder:  minimap2_purge_dups_qc > minimap2_purge_dups_reads
ruleorder: get_purge_dups_read_stat_qc > get_purge_dups_read_stat

rule create_primary_contig_link:
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap1.fasta" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                            stage_dict["purge_dups"]["prev_stage"])),
        len=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap1.len" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                            stage_dict["purge_dups"]["prev_stage"]))
    output:
        fasta=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap1.fasta",
        len=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap1.len"
        #fasta=out_dir_path / ("purge_dups/{assembler}/input/%s.contig.{assembler}.hap1.fasta" % config["genome_name"])
    log:
        ln1=output_dict["log"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.ln1.log",
        ln2=output_dict["log"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.ln2.log",
        cluster_log=output_dict["cluster_log"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.hap1.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln -s ../../../../{input.fasta} {output.fasta} 1>{log.ln1} 2>&1; "
        " ln -s ../../../../{input.len} {output.len} 1>{log.ln2} 2>&1"

rule minimap2_purge_dups_reads:
    input:
        fastq=lambda wildcards: output_dict["data"] / "fastq/{0}/filtered/{1}{2}".format(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"],
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]),
        reference=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta"
    output:
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.{haplotype}.{fileprefix}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -x "),
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " minimap2 {params.alignment_scheme} {params.index_size} -t {threads}  {input.reference} "
        " {input.fastq} 2>{log.std} |  gzip -c - > {output.paf} "

rule get_purge_dups_read_stat: #TODO: adjust -d -m -u options for calcuts
    input:
        paf=lambda wildcards: expand(rules.minimap2_purge_dups_reads.output.paf,
                   fileprefix=input_file_prefix_dict[stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]],
                   genome_prefix=[config["genome_prefix"]],
                   allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{1}.{0}.filtered.{2}.{3}.genomescope.parameters".format(config["final_kmer_datatype"],
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"])
    output:
        pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.stat",
        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov",
        cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/cutoffs",
        len=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov.len",
        stat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov.stat",
        bed=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov.bed"
    params:
        #out_dir=lambda wildcards: out_dir_path  / "purge_dups/{0}..{1}/{2}".format(wildcards.prev_stage_parameters,
        #                                                                          wildcards.purge_dups_parameters,
        #                                                                          wildcards.haplotype),
        cov_multiplicator=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["cov_multiplicator"]
        #cov_multiplicator=parameters["tool_options"]["purge_dups"]["cov_multiplicator"]

    log:
        pbstat=output_dict["log"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.pbstat.log",
        calcuts=output_dict["log"]  / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.calcuts.log",
        convert=output_dict["log"]  / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.convert.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]

    shell:
        " OUT_DIR=`dirname {output.pbstat}`; "
        " LEN_FILE={output.len}; "
        " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`;"
        " pbcstat -O ${{OUT_DIR}} {input.paf} 1>{log.pbstat} 2>&1; "
        " calcuts -d 1 -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts}; " #check parameters for calcuts
        " convert_coverage_file_to_bed.py -i {output.pbbasecov}  -o ${{LEN_FILE%.len}} > {log.convert} 2>&1; "

rule minimap2_purge_dups_assembly:
    input:
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta"
    output:
        split_reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.purge_dups_input.{haplotype}.split.fasta",
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.purge_dups_input.{haplotype}.split.minimap2.self.paf.gz"
    params:
        index_size=parse_option("index_size", parameters["tool_options"]["minimap2"]["self"], " -I "),
        alignment_scheme=parse_option("alignment_scheme", parameters["tool_options"]["minimap2"]["self"], " -x "),
    log:
        split_fa=output_dict["log"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.split_fa.log",
        minimap2=output_dict["log"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.minimap2.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " split_fa {input.reference} > {output.split_reference} 2>{log.split_fa};"
        " minimap2 -DP {params.alignment_scheme} {params.index_size} -t {threads}  {output.split_reference} "
        " {output.split_reference} 2>{log.minimap2} |  gzip -c - > {output.paf}; "

rule purge_dups: # TODO: find what options are used in ERGA for get_seqs
    input:
        cutoffs=rules.get_purge_dups_read_stat.output.cutoffs,
        pbbasecov=rules.get_purge_dups_read_stat.output.pbbasecov,
        self_paf=rules.minimap2_purge_dups_assembly.output.paf,
        reference = out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta"
    output:
        bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.bed",
        purged=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.purge_dups.{haplotype}.purged.fasta",
        hapdups=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.purge_dups.{haplotype}.hap.fasta",
        purged_alias=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta",
    params:
        #bed_local_path=lambda wildcards: "{0}.dups.bed".format(wildcards.genome_prefix),
        #out_dir=lambda wildcards: out_dir_path  / "purge_dups/{0}..{1}/{2}".format(wildcards.prev_stage_parameters,
        #                                                                         wildcards.purge_dups_parameters,
        #                                                                         wildcards.haplotype),
        get_seq_prefix=lambda wildcards: "{0}.purge_dups.{1}".format(wildcards.genome_prefix, wildcards.haplotype)
    log:
        purge_dups=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.purge_dups.log",
        get_seqs=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.get_seqs.log",
        ln=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " OUT_DIR=`dirname {output.bed}`; "
        " purge_dups -2 -T {input.cutoffs} -c {input.pbbasecov} {input.self_paf} > {output.bed} 2>{log.purge_dups}; "
        " PURGE_DUPS_BED=`realpath {output.bed}`; "
        " REFERENCE=`realpath {input.reference}`; "
        " GET_SEQ_LOG=`realpath {log.get_seqs}`; "
        " cd ${{OUT_DIR}}; "
        " get_seqs -p {params.get_seq_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1; "
        " for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done; "
        " ln `basename {output.purged}` ../`basename {output.purged_alias}` > {log.ln} 2>&1; "

rule merge_pri_hapdups_with_alt: # TODO: add handling of polyploid cases
    input:
        alt_contig=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap2.fasta" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                                 stage_dict["purge_dups"]["prev_stage"])),
        pri_hapdups=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/hap1/{genome_prefix}.purge_dups.hap1.hap.fasta",
    output:
        alt_plus_pri_hapdup=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap2.fasta",

    log:
        std=output_dict["log"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.log",
        cluster_log=output_dict["cluster_log"] / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["merge_pri_hapdups_with_alt"] ,
        time=parameters["time"]["merge_pri_hapdups_with_alt"],
        mem=parameters["memory_mb"]["merge_pri_hapdups_with_alt"]
    threads: parameters["threads"]["merge_pri_hapdups_with_alt"]

    shell:
        " cat {input.alt_contig} {input.pri_hapdups} > {output.alt_plus_pri_hapdup} 2>{log.std}; "

rule merge_pri_hapdups_with_alt_for_len_files: # TODO: add handling of polyploid cases
    input:
        alt_len=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap2.len" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                                 stage_dict["purge_dups"]["prev_stage"])),
        pri_hapdups_len=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/hap1/{genome_prefix}.purge_dups.hap1.hap.len",
    output:
        alt_plus_pri_len=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap2.len",
    log:
        std=output_dict["log"]  / "merge_pri_hapdups_with_alt_for_len_files.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.log",
        cluster_log=output_dict["cluster_log"] / "merge_pri_hapdups_with_alt_for_len_files.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_pri_hapdups_with_alt_for_len_files.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt_for_len_files.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["merge_pri_hapdups_with_alt"] ,
        time=parameters["time"]["merge_pri_hapdups_with_alt"],
        mem=parameters["memory_mb"]["merge_pri_hapdups_with_alt"]
    threads: parameters["threads"]["merge_pri_hapdups_with_alt"]

    shell:
        " cat {input.alt_len} {input.pri_hapdups_len} > {output.alt_plus_pri_len} 2>{log.std}"

#rule create_link_for_purged_fasta: # command moved to rules.purge_dups
#    input:
#        purged=rules.purge_dups.output.purged
#    output:
#        purged=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta"
#    log:
#        std=output_dict["log"]  / "create_link_for_purged_fasta.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.log",
#        cluster_log=output_dict["cluster_log"] / "create_link_for_purged_fasta.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
#        cluster_err=output_dict["cluster_error"] / "create_link_for_purged_fasta.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
#    benchmark:
#        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
#    conda:
#        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
#    resources:
#        cpus=parameters["threads"]["create_links"] ,
#        time=parameters["time"]["create_links"],
#        mem=parameters["memory_mb"]["create_links"]
#    threads: parameters["threads"]["create_links"]
#
#    shell:
#        " ln {input.purged} {output.purged} > {log.std} 2>&1;"

#rule extract_coverage_from_purge_dups_file: # moved to rules.get_purge_dups_read_stat
#    input:
#        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.base.cov"
#    output:
#        len=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov.len",
#        stat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov.stat",
#        bed=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/PB.base.cov.bed"
#    log:
#        std=output_dict["log"]  / "extract_coverage_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.log",
#        cluster_log=output_dict["cluster_log"] / "extract_coverage_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.log",
#        cluster_err=output_dict["cluster_error"] / "extract_coverage_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.err"
#    benchmark:
#        output_dict["benchmark"]  / "extract_coverage_from_purge_dups_file..{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
#    conda:
#        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
#    resources:
#        cpus=parameters["threads"]["extract_coverage_from_purge_dups_file"] ,
#        time=parameters["time"]["extract_coverage_from_purge_dups_file"],
#        mem=parameters["memory_mb"]["extract_coverage_from_purge_dups_file"]
#    threads: parameters["threads"]["extract_coverage_from_purge_dups_file"]
#
#    shell:
#        " LEN_FILE={output.len}; "
#        " convert_coverage_file_to_bed.py -i {input.pbbasecov}  -o ${{LEN_FILE%.len}} > {log.std} 2>&1; "

rule extract_stats_from_purge_dups_file:
    input:
        stat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.base.cov.stat",
        bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.dups.bed",
        len=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.len"
    output:
        extended_bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.extended.bed",
        stat=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.stat",
        junk_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.junk.ids",
        ovlp_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.ovlp.ids",
        haplotig_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.haplotig.ids",
        repeat_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.repeat.ids",
        highcov_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.highcov.ids"
    log:
        std=output_dict["log"]  / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_stats_from_purge_dups_file"] ,
        time=parameters["time"]["extract_stats_from_purge_dups_file"],
        mem=parameters["memory_mb"]["extract_stats_from_purge_dups_file"]
    threads: parameters["threads"]["extract_stats_from_purge_dups_file"]

    shell:
        " STATS_FILE={output.stat}; "
        " ./workflow/scripts/purge_dups/calculate_purge_dups_stats.py  -b {input.bed} -s {input.stat} -l {input.len} "
        " -o ${{STATS_FILE%.stat}} > {log.std} 2>&1; "

rule extract_artefact_sequences:
    input:
        artefact_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/{genome_prefix}.dups.{artefact}.ids",
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.fasta",
        len_file=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.{haplotype}.len"
    output:
        artefact_fasta=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype, [^.]+}/{genome_prefix}.dups.{artefact}.fasta"
    log:
        std=output_dict["log"]  / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{artefact}.log",
        cluster_log=output_dict["cluster_log"] / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.{artefact}.log",
        cluster_err=output_dict["cluster_error"] / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.{artefact}.err"
    benchmark:
        output_dict["benchmark"]  / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{artefact}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_artefact_sequences"] ,
        time=parameters["time"]["extract_artefact_sequences"],
        mem=parameters["memory_mb"]["extract_artefact_sequences"]
    threads: parameters["threads"]["extract_artefact_sequences"]

    shell:
        " extract_sequences_by_ids.py -i {input.reference} -d {input.artefact_ids} "
        " -o {output.artefact_fasta} > {log.std} 2>&1 ; "


rule minimap2_purge_dups_qc:
    input:
        fastq=lambda wildcards: output_dict["data"] / "fastq/{0}/filtered/{1}{2}".format(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"],
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]),
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta"
    output:
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/{genome_prefix}.{haplotype}.{fileprefix}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -x "),
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " minimap2 {params.alignment_scheme} -I {params.index_size} -t {threads}  {input.reference} "
        " {input.fastq} 2>{log.std} |  gzip -c - > {output.paf}; "

rule get_purge_dups_read_stat_qc:
    input:
        paf=lambda wildcards: expand(rules.minimap2_purge_dups_qc.output.paf,
                           genome_prefix=[config["genome_prefix"]],
                           fileprefix=input_file_prefix_dict[stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]],
                           allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{1}.{0}.filtered.{2}.{3}.genomescope.parameters".format(config["final_kmer_datatype"],
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"]),
        before_pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.stat",
        before_pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.base.cov",
        before_cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/cutoffs"
    output:
        pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/PB.stat",
        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/PB.base.cov",
        cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/cutoffs",
        coverage_plot=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/{haplotype}.before-after.comparison.coverage.png"
    params:
        cov_multiplicator=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["cov_multiplicator"]
    log:
        pbstat=output_dict["log"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.pbstat.log",
        png=output_dict["log"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.png.log",
        calcuts=output_dict["log"]  / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]

    shell:
        " OUT_DIR=`dirname {output.pbbasecov}`;"
        " COV_PLOT={output.coverage_plot}; "
        " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`;"
        " pbcstat -O ${{OUT_DIR}} {input.paf} 1>{log.pbstat} 2>&1; "
        " calcuts -d 1 -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts}; " #check parameters for calcuts
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py -b {input.before_pbstat},{output.pbstat} "
        " -l before,after -c {input.before_cutoffs},{output.cutoffs} -e png,svg -o ${{COV_PLOT%.png}} > {log.png} 2>&1; "

rule get_purge_stat_haplotype_comparison: #TODO: adjust -d -m -u options for calcuts
    input:
        before_pbstat=expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.stat",
                             haplotype=haplotype_list,
                             allow_missing=True,),
        before_cutoffs=expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/cutoffs",
                             haplotype=haplotype_list,
                             allow_missing=True,),
        after_pbstat=expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                             haplotype=haplotype_list,
                             allow_missing=True,),
        after_cutoffs=expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/cutoffs",
                             haplotype=haplotype_list,
                             allow_missing=True,),
    output:
        before_coverage_plot=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
        after_coverage_plot=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/after.comparison.coverage.png"
    params:
        label_list=haplotype_list,
        label_string=",".join(haplotype_list)
    log:
        before=output_dict["log"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.before.log",
        after=output_dict["log"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.after.log",
        calcuts=output_dict["log"]  / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["get_purge_stat_haplotype_comparison"] ,
        time=parameters["time"]["get_purge_stat_haplotype_comparison"],
        mem=parameters["memory_mb"]["get_purge_stat_haplotype_comparison"]
    threads: parameters["threads"]["get_purge_stat_haplotype_comparison"]

    shell:
        " BEFORE_COV_PLOT={output.before_coverage_plot}; "
        " AFTER_COV_PLOT={output.after_coverage_plot}; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py "
        " -b `echo {input.before_pbstat} | tr ' ' ','` "
        " -l {params.label_string} -c `echo {input.before_cutoffs} | tr ' ' ','` "
        " -e png,svg -o ${{BEFORE_COV_PLOT%.png}} > {log.before} 2>&1; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py "
        " -b `echo {input.after_pbstat} | tr ' ' ','` "
        " -l {params.label_string} -c `echo {input.after_cutoffs} | tr ' ' ','` "
        " -e png,svg -o ${{AFTER_COV_PLOT%.png}} > {log.after} 2>&1; "