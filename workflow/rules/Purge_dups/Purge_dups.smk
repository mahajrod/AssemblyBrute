ruleorder: create_contig_links > get_seq_len
#ruleorder: create_links_if_skipping_purge_dups > create_final_links_purge_dups # rule to skip purge_dups when necessary
#ruleorder: create_assembly_qc_if_skipping_purge_dups > # rule to skip assembly qc rules when skipping purge_dups
#ruleorder: purge_dups > gfa2fasta
#ruleorder: create_final_links_purge_dups > gfa2fasta
localrules: create_contig_links, create_final_links_purge_dups, extract_stats_from_purge_dups_file #create_link_for_purged_fasta,

#localrules: create_primary_contig_link_hap1, create_primary_contig_link_hap0, merge_pri_hapdups_with_alt, extract_stats_from_purge_dups_file #create_link_for_purged_fasta,
#localrules: merge_pri_hapdups_with_alt_for_len_files # create_primary_contig_len_file_link,
#ruleorder: create_primary_contig_link_hap1 > merge_pri_hapdups_with_alt
#ruleorder: create_primary_contig_link_hap0 > merge_pri_hapdups_with_alt
"""
rule create_assembly_links_if_skipping_purge_dups:
    input:
        fasta=out_dir_path  / "{prev_stage}/{prev_stage}..{prev_stage_parameters}/{genome_prefix}.{prev_stage}.{haplotype, [^.]+}.fasta",
        len=out_dir_path  / "{prev_stage}/{prev_stage}..{prev_stage_parameters}/{genome_prefix}.{prev_stage}.{haplotype, [^.]+}.len",
        fai=out_dir_path  / "{prev_stage}/{prev_stage}..{prev_stage_parameters}/{genome_prefix}.{prev_stage}.{haplotype, [^.]+}.fasta.fai",

    output:
        fasta_alias=out_dir_path / "purge_dups/{prev_stage}_{prev_stage_parameters}..purge_dups_skipped/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta",
        len_alias=out_dir_path / "purge_dups/{prev_stage}_{prev_stage_parameters}..purge_dups_skipped/{genome_prefix}.purge_dups.{haplotype, [^.]+}.len",
        fai_alias=out_dir_path / "purge_dups/{prev_stage}_{prev_stage_parameters}..purge_dups_skipped/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta.fai"
        #assembly_qc_alias=directory(out_dir_path / "purge_dups/{prev_stage}_{prev_stage_parameters}..purge_dups_skipped/assembly_qc/")
    log:
        ln=output_dict["log"]  / "create_assembly_links_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_assembly_links_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_assembly_links_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_assembly_links_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_contig_links"),
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln -sf {input.fasta} {output.fasta_alias} > {log.ln} 2>&1; "
        " ln -sf {input.len} {output.len_alias} >> {log.ln} 2>&1; "
        " ln -sf {input.fai} {output.fai_alias} >> {log.ln} 2>&1; "

rule create_assembly_qc_if_skipping_purge_dups:
    input:
        assembly_qc=out_dir_path  / "{prev_stage}/{prev_stage}..{prev_stage_parameters}/assembly_qc",
    output:
        assembly_qc=directory(out_dir_path / "purge_dups/{prev_stage}_{prev_stage_parameters}..purge_dups_skipped/assembly_qc"),
    log:
        ln=output_dict["log"]  / "create_assembly_qc_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_assembly_qc_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_assembly_qc_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_assembly_qc_if_skipping_purge_dups.purge_dups.{prev_stage}_{prev_stage_parameters}.purge_dups_skipped.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_contig_links"),
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln -sf {input.assembly_qc} {output.assembly_qc} > {log.ln} 2>&1; "
"""
rule create_contig_links:
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                            stage_dict["purge_dups"]["prev_stage"])),
        len=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.len" % (stage_dict["purge_dups"]["prev_stage"],
                                                                                            stage_dict["purge_dups"]["prev_stage"]))
    output:
        fasta=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.fasta",
        len=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.len"
        #fasta=out_dir_path / ("purge_dups/{assembler}/input/%s.contig.{assembler}.hap1.fasta" % config["genome_name"])
    log:
        ln1=output_dict["log"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.ln1.log",
        ln2=output_dict["log"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.ln2.log",
        cluster_log=output_dict["cluster_log"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_contig_links"),
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln -sf ../../../../../{input.fasta} {output.fasta} 1>{log.ln1} 2>&1; "
        " ln -sf ../../../../../{input.len} {output.len} 1>{log.ln2} 2>&1"

rule minimap2_purge_dups_reads:
    input:
        fastq=lambda wildcards: output_dict["data"] / "fastq/{0}/filtered/{1}{2}".format(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"],
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]),
        reference=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.{haplotype}.{fileprefix}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -x "),
        #min_mapq=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]]["min_mapping_quality"]
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{purge_stage}.{genome_prefix}.{fileprefix}.log",
        #awk=output_dict["log"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}.awk.log",
        gzip=output_dict["log"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{purge_stage}.{genome_prefix}.{fileprefix}.gzip.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{purge_stage}.{genome_prefix}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{purge_stage}.{genome_prefix}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{purge_stage}.{genome_prefix}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("minimap2_purge_dups_reads"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell: # awk -F'\t' '{{if ($12 >= {params.min_mapq}) print $0 }}' 2>{log.awk} | "
        " minimap2 {params.alignment_scheme} {params.index_size} -t {threads}  {input.reference} "
        " {input.fastq} 2>{log.std} | "
        "  gzip -c - > {output.paf} 2>{log.gzip}; "

rule get_purge_dups_read_stat:
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
        pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/PB.stat",
        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/PB.base.cov",
        cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/cutoffs",
        len=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/PB.base.cov.len",
        stat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/PB.base.cov.stat",
        bed=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/PB.base.cov.bed"
    params:
        cov_multiplicator=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["cov_multiplicator"],
        calcuts_lower_threshold=lambda wildcards: parse_option("lower_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -l "),
        calcuts_haploid_diploid_threshold=lambda wildcards: parse_option("haploid_diploid_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -m "),
        calcuts_upper_threshold=str(config["tool_manually_adjusted_features"]["calcuts"]["upper_threshold"]), # None needs to be converted to "None"
    log:
        pbstat=output_dict["log"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.{purge_stage}.pbstat.log",
        calcuts=output_dict["log"]  / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.{purge_stage}.calcuts.log",
        convert=output_dict["log"]  / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.{purge_stage}.convert.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.{purge_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_dups_read_stat.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.{purge_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.{purge_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_purge_dups_read_stat"),
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]
    shell:
        " OUT_DIR=`dirname {output.pbstat}`; "
        " LEN_FILE={output.len}; "
        " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "
        " if [ '{params.calcuts_upper_threshold}' != 'None' ] ; then COV_UPPER_BOUNDARY={params.calcuts_upper_threshold}; fi; "
        " pbcstat -O ${{OUT_DIR}} {input.paf} 1>{log.pbstat} 2>&1; "
        " calcuts -d 1 {params.calcuts_lower_threshold} {params.calcuts_haploid_diploid_threshold} "
        " -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts}; " #check parameters for calcuts
        " convert_coverage_file_to_bed.py -i {output.pbbasecov}  -o ${{LEN_FILE%.len}} > {log.convert} 2>&1; "

rule minimap2_purge_dups_assembly:
    input:
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        split_reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.split.fasta",
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.split.minimap2.self.paf.gz"
    params:
        index_size=parse_option("index_size", parameters["tool_options"]["minimap2"]["self"], " -I "),
        alignment_scheme=parse_option("alignment_scheme", parameters["tool_options"]["minimap2"]["self"], " -x "),
    log:
        split_fa=output_dict["log"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.split_fa.log",
        minimap2=output_dict["log"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.minimap2.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_assembly.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("minimap2_purge_dups_assembly"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " split_fa {input.reference} > {output.split_reference} 2>{log.split_fa};"
        " minimap2 -DP {params.alignment_scheme} {params.index_size} -t {threads}  {output.split_reference} "
        " {output.split_reference} 2>{log.minimap2} |  gzip -c - > {output.paf}; "

rule purge_dups: #
    input:
        cutoffs=rules.get_purge_dups_read_stat.output.cutoffs,
        pbbasecov=rules.get_purge_dups_read_stat.output.pbbasecov,
        self_paf=rules.minimap2_purge_dups_assembly.output.paf,
        #reference = out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.raw.bed",
        #purged=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.purged.fasta",
        #hapdups=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.hap.fasta",
        #purged_alias=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta",
    #params:
        #get_seq_prefix=lambda wildcards: "{0}.input.{1}".format(wildcards.genome_prefix, wildcards.haplotype)
    log:
        purge_dups=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.purge_dups.log",
        #get_seqs=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.get_seqs.log",
        #ln=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.ln.log",
        cluster_log=output_dict["cluster_log"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("purge_dups"),
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " OUT_DIR=`dirname {output.bed}`; "
        " purge_dups -2 -T {input.cutoffs} -c {input.pbbasecov} {input.self_paf} > {output.bed} 2>{log.purge_dups}; "
        #" PURGE_DUPS_BED=`realpath -s {output.bed}`; "
        #" REFERENCE=`realpath -s {input.reference}`; "
        #" GET_SEQ_LOG=`realpath -s {log.get_seqs}`; "
        #" LN_LOG=`realpath -s {log.ln}`; "
        #" cd ${{OUT_DIR}}; "
        #" get_seqs -p {params.get_seq_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1; "
        #" for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done; "

rule get_purged_seqs: #
    input:
        #cutoffs=rules.get_purge_dups_read_stat.output.cutoffs,
        #pbbasecov=rules.get_purge_dups_read_stat.output.pbbasecov,
        #self_paf=rules.minimap2_purge_dups_assembly.output.paf,
        raw_dups_bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.raw.bed",
        reference = out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        filtered_bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.filtered.bed",
        purged=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.purged.fasta",
        hapdups=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.hap.fasta",
        #purged_alias=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta",
    params:
        blacklist_option=lambda wildcards: parse_option("purging_blacklist",
                                                        stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"],
                                                        option_prefix="-b"),
        whitelist_option=lambda wildcards: parse_option("purging_whitelist",
                                                        stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"],
                                                        option_prefix="-w"),
        get_seq_prefix=lambda wildcards: "{0}.input.{1}".format(wildcards.genome_prefix, wildcards.haplotype)
    log:
        #purge_dups=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.purge_dups.log",
        get_seqs=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.get_seqs.log",
        filter=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.filter.log",
        ln=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.ln.log",
        cluster_log=output_dict["cluster_log"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("purge_dups"),
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " workflow/scripts/purge_dups/filter_dups.bed.py -i {input.raw_dups_bed} "
        " {params.blacklist_option} {params.whitelist_option} -o {output.filtered_bed} > {log.filter} 2>&1; "
        " OUT_DIR=`dirname {output.filtered_bed}`; "
        " PURGE_DUPS_BED=`realpath -s {output.filtered_bed}`; "
        " REFERENCE=`realpath -s {input.reference}`; "
        " GET_SEQ_LOG=`realpath -s {log.get_seqs}`; "
        " LN_LOG=`realpath -s {log.ln}`; "
        " cd ${{OUT_DIR}}; "
        " get_seqs -p {params.get_seq_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1; "
        " for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done; "


        #" ln -sf {wildcards.haplotype}/`basename {output.purged}` ../`basename {output.purged_alias}` > ${{LN_LOG}} 2>&1; "
"""
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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("merge_pri_hapdups_with_alt"),
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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("merge_pri_hapdups_with_alt_for_len_files"),
        cpus=parameters["threads"]["merge_pri_hapdups_with_alt"] ,
        time=parameters["time"]["merge_pri_hapdups_with_alt"],
        mem=parameters["memory_mb"]["merge_pri_hapdups_with_alt"]
    threads: parameters["threads"]["merge_pri_hapdups_with_alt"]

    shell:
        " cat {input.alt_len} {input.pri_hapdups_len} > {output.alt_plus_pri_len} 2>{log.std}"
"""

rule filter_removed_contigs: # TODO: find what options are used in ERGA for get_seqs
    input:
        stat=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.stat",
        #junk_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.junk.ids",
        #ovlp_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.ovlp.ids",
        #highcov_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.highcov.ids",
        hapdups=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.hap.fasta",
    output:
        interhaplotype_transfer_id_blacklist=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.transfer_blacklist.ids",
        filtered_hapdups=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.hap.for_transfer.fasta",
    params:
        blacklist=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["interhaplotype_transfer_blacklist"]
    log:
        extract=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.extract.log",
        cluster_log=output_dict["cluster_log"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("purge_dups"),
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " STAT={input.stat}; "
        " > {output.interhaplotype_transfer_id_blacklist};"
        " for ENTRY in {params.blacklist};"
        "   do "
        "   cat ${{STAT%.stat}}.${{ENTRY}}.ids >> {output.interhaplotype_transfer_id_blacklist}; "
        "   done; "
        " extract_sequences_by_ids.py -i {input.hapdups} -o {output.filtered_hapdups} "
        " -d {output.interhaplotype_transfer_id_blacklist} -r  -p parse > {log.extract} 2>&1; "

rule crossmerge_hapdups_with_deduped_contigs: # TODO: add handling of polyploid cases
    input:
        purged=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype}/{genome_prefix}.input.{haplotype}.purged.fasta",
        hapdups=lambda wildcards: out_dir_path  / "purge_dups/{0}..{1}/first_stage/{2}/{3}.input.{2}.hap.for_transfer.fasta".format(wildcards.prev_stage_parameters,
                                                                                                                wildcards.purge_dups_parameters,
                                                                                                                list(set(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"]) - set([wildcards.haplotype]))[0],
                                                                                                                wildcards.genome_prefix)
        #alt_contig=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.hap2.fasta" % (stage_dict["purge_dups"]["prev_stage"],
        #                                                                                         stage_dict["purge_dups"]["prev_stage"])),
        #pri_hapdups=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/hap1/{genome_prefix}.purge_dups.hap1.hap.fasta",
    output:
        merged=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/second_stage/{haplotype, [^.]+}/{genome_prefix}.input.{haplotype}.fasta"
        #alt_plus_pri_hapdup=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/input/{genome_prefix}.purge_dups_input.hap2.fasta",

    log:
        std=output_dict["log"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.purge_dups.log",
        cluster_log=output_dict["cluster_log"] / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "merge_pri_hapdups_with_alt.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("merge_pri_hapdups_with_alt"),
        cpus=parameters["threads"]["merge_pri_hapdups_with_alt"] ,
        time=parameters["time"]["merge_pri_hapdups_with_alt"],
        mem=parameters["memory_mb"]["merge_pri_hapdups_with_alt"]
    threads: parameters["threads"]["merge_pri_hapdups_with_alt"]

    shell:
        " cat {input.purged} {input.hapdups} > {output.merged} 2>{log.std}; "

rule create_final_links_purge_dups:
    input:
        fasta=lambda wildcards: out_dir_path  / "purge_dups/{0}..{1}/{2}/{3}/{4}.input.{3}.purged.fasta".format(wildcards.prev_stage_parameters,
                                                                                                                wildcards.purge_dups_parameters,
                                                                                                                "first_stage" if wildcards.haplotype == "hap0" else "second_stage",
                                                                                                                wildcards.haplotype,
                                                                                                                wildcards.genome_prefix),
    output:
        purged_alias=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype, [^.]+}.fasta"
    log:
        ln=output_dict["log"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_contig_links.purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_contig_links"),
        cpus=parameters["threads"]["create_links"] ,
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " ln -sf second_stage/{wildcards.haplotype}/{wildcards.genome_prefix}.input.{wildcards.haplotype}.purged.fasta {output.purged_alias} > {log.ln} 2>&1; "

rule extract_stats_from_purge_dups_file:
    input:
        stat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/PB.base.cov.stat",
        bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.filtered.bed",
        len=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.len"
    output:
        extended_bed=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.extended.bed",
        stat=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.stat",
        junk_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.junk.ids",
        ovlp_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.ovlp.ids",
        haplotig_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.haplotig.ids",
        repeat_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.repeat.ids",
        highcov_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.highcov.ids"
    log:
        std=output_dict["log"]  / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.log",
        cluster_log=output_dict["cluster_log"] / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "extract_stats_from_purge_dups_file.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("extract_stats_from_purge_dups_file"),
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
        artefact_ids=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.{artefact}.ids",
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.fasta",
        len_file=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype}/{genome_prefix}.input.{haplotype}.len"
    output:
        artefact_fasta=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{purge_stage}/{haplotype, [^.]+}/{genome_prefix}.dups.{artefact}.fasta"
    log:
        std=output_dict["log"]  / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.{artefact}.log",
        cluster_log=output_dict["cluster_log"] / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.{purge_stage}.{artefact}.log",
        cluster_err=output_dict["cluster_error"] / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.cluster.{purge_stage}.{artefact}.err"
    benchmark:
        output_dict["benchmark"]  / "extract_artefact_sequences.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.{artefact}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("extract_artefact_sequences"),
        cpus=parameters["threads"]["extract_artefact_sequences"] ,
        time=parameters["time"]["extract_artefact_sequences"],
        mem=parameters["memory_mb"]["extract_artefact_sequences"]
    threads: parameters["threads"]["extract_artefact_sequences"]

    shell:
        " extract_sequences_by_ids.py -i {input.reference} -d {input.artefact_ids} "
        " -o {output.artefact_fasta} > {log.std} 2>&1 ; "
