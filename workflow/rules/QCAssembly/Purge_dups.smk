
wildcard_constraints:
    haplotype="[^.]+"


rule qc_minimap2_purge_dups_reads:
    input:
        fastq=lambda wildcards: output_dict["data"] / "fastq/{0}/filtered/{1}{2}".format(wildcards.datatype, #stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"],
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]),
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        log_dir=out_dir_path / "{assembly_stage}/{parameters}/log/"
    output:
        paf=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, [^./]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{haplotype, hap[^.]+}.{fileprefix, [^/]+}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.datatype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.datatype], " -x "),
        #min_mapq=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]]["min_mapping_quality"]
    log:
        std=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_reads.{assembly_stage}.{parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}.log",
        #awk=output_dict["log"]  / "minimap2_purge_dups_reads.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}.awk.log",
        gzip=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_reads.{assembly_stage}.{parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}.gzip.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_reads.{assembly_stage}.{parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}.{fileprefix}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_reads.{assembly_stage}.{parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}.{fileprefix}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_reads.{assembly_stage}.{parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("minimap2_purge_dups_reads"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["minimap2"], mode="linear")
    threads: parameters["threads"]["minimap2"]

    shell: # awk -F'\t' '{{if ($12 >= {params.min_mapq}) print $0 }}' 2>{log.awk} | "
        " minimap2 {params.alignment_scheme} {params.index_size} -t {threads}  {input.assembly} "
        " {input.fastq} 2>{log.std} | "
        "  gzip -c - > {output.paf} 2>{log.gzip}; "

#def get_paf_list(wildcards):
#    paf_list = []
#    for datatype in stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"]["main_datatypes"]:
#        paf_list += expand(rules.minimap2_purge_dups_reads.output.paf,
#                           datatype=[datatype],
#                           fileprefix=input_file_prefix_dict[datatype],
#                           genome_prefix=[config["genome_prefix"]],
#                           allow_missing=True)
#    return paf_list

rule qc_get_purge_dups_read_stat:
    input:
        paf=lambda wildcards: expand(rules.minimap2_purge_dups_reads.output.paf,
                                     fileprefix=input_file_prefix_dict[wildcards.datatype],
                                     genome_prefix=[config["genome_prefix"]],
                                     allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{1}.{0}.filtered.{2}.{3}.genomescope.parameters".format("_".join(config["final_kmer_datatypes"]),
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"]),
        log_dir=out_dir_path / "{assembly_stage}/{parameters}/log/"
    output:
        pbstat=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/PB.stat",
        pbbasecov=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/PB.base.cov",
        cutoffs=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/cutoffs",
        len=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/PB.base.cov.len",
        stat=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/PB.base.cov.stat",
        bed=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/PB.base.cov.bed"
    params:
        cov_multiplicator=lambda wildcards: stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"]["cov_multiplicator"],
        calcuts_lower_threshold=parse_option("lower_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -l "),
        calcuts_haploid_diploid_threshold=parse_option("haploid_diploid_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -m "),
        calcuts_upper_threshold=str(config["tool_manually_adjusted_features"]["calcuts"]["upper_threshold"]), # None needs to be converted to "None"
    log:
        pbstat=out_dir_path / "{assembly_stage}/{parameters}/log/get_purge_dups_read_stat.{haplotype}.{datatype}.pbstat.log",
        calcuts=out_dir_path / "{assembly_stage}/{parameters}/log/get_purge_dups_read_stat.{haplotype}.{datatype}.calcuts.log",
        convert=out_dir_path / "{assembly_stage}/{parameters}/log/get_purge_dups_read_stat.{haplotype}.{datatype}.convert.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/log/get_purge_dups_read_stat.{haplotype}.{datatype}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/log/get_purge_dups_read_stat.{haplotype}.{datatype}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_reads.{haplotype}.{datatype}.benchmark.txt"
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

rule qc_minimap2_purge_dups_assembly:
    input:
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        log_dir=out_dir_path / "{assembly_stage}/{parameters}/log/"
    output:
        split_assembly=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.split.fasta",
        paf=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.split.minimap2.self.paf.gz"
    params:
        index_size=parse_option("index_size", parameters["tool_options"]["minimap2"]["self"], " -I "),
        alignment_scheme=parse_option("alignment_scheme", parameters["tool_options"]["minimap2"]["self"], " -x "),
    log:
        split_fa=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_assembly.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.split_fa.log",
        minimap2=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_assembly.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.minimap2.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_assembly.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_assembly.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/log/minimap2_purge_dups_assembly.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("minimap2_purge_dups_assembly"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["minimap2"], mode="linear")
    threads: parameters["threads"]["minimap2"]

    shell:
        " split_fa {input.assembly} > {output.split_assembly} 2>{log.split_fa};"
        " minimap2 -DP {params.alignment_scheme} {params.index_size} -t {threads}  {output.split_assembly} "
        " {output.split_assembly} 2>{log.minimap2} |  gzip -c - > {output.paf}; "

rule qc_purge_dups: #
    input:
        cutoffs=rules.get_purge_dups_read_stat.output.cutoffs,
        pbbasecov=rules.get_purge_dups_read_stat.output.pbbasecov,
        self_paf=rules.minimap2_purge_dups_assembly.output.paf,
        log_dir=out_dir_path / "{assembly_stage}/{parameters}/log/"
    output:
        bed=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.raw.bed",
    log:
        purge_dups=out_dir_path / "{assembly_stage}/{parameters}/log/qc_purge_dups.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.purge_dups.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/log/qc_purge_dups.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/log/qc_purge_dups.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/log/qc_purge_dups.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.benchmark.txt"
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

rule qc_get_purged_seqs: #
    input:
        raw_dups_bed=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.raw.bed",
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        log_dir=out_dir_path / "{assembly_stage}/{parameters}/log/"
    output:
        filtered_bed=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.filtered.bed",
        purged=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.purged.fasta",
        hapdups=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{assembly_stage}.{haplotype}.hap.fasta",
    params:
        blacklist_option=lambda wildcards: parse_option("purging_blacklist",
                                                        stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"],
                                                        option_prefix="-b", expression=lambda l: ",".join(l)),
        whitelist_option=lambda wildcards: parse_option("purging_whitelist",
                                                        stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"],
                                                        option_prefix="-w", expression=lambda l: ",".join(l)),
        get_seq_prefix=lambda wildcards: "{0}.{1}.{2}".format(wildcards.genome_prefix, wildcards.assembly_stage, wildcards.haplotype)
    log:
        #purge_dups=output_dict["log"]  / "purge_dups.{prev_stage_parameters}.{purge_dups_parameters}.{genome_prefix}.purge_dups.{haplotype}.{purge_stage}.purge_dups.log",
        get_seqs=out_dir_path / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.get_seqs.log",
        filter=out_dir_path / "{assembly_stage}/{parameters}/log/qc_get_purged_seqspurge_.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.filter.log",
        ln=out_dir_path / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.ln.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.benchmark.txt"
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
        " workflow/scripts/purge_dups/filter_dups_bed.py -i {input.raw_dups_bed} "
        " {params.blacklist_option} {params.whitelist_option} -o {output.filtered_bed} > {log.filter} 2>&1; "
        " OUT_DIR=`dirname {output.filtered_bed}`; "
        " PURGE_DUPS_BED=`realpath -s {output.filtered_bed}`; "
        " REFERENCE=`realpath -s {input.assembly}`; "
        " GET_SEQ_LOG=`realpath -s {log.get_seqs}`; "
        " LN_LOG=`realpath -s {log.ln}`; "
        " cd ${{OUT_DIR}}; "
        " get_seqs -p {params.get_seq_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1; "
        " for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done; "

rule qc_extract_stats_from_purge_dups_file:
    input:
        bed=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.raw.bed",
        stat=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/PB.base.cov.stat",
        len=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
        log_dir=out_dir_path / "{assembly_stage}/{parameters}/log/",
    output:
        extended_bed=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.extended.bed",
        stat=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.stat",
        junk_ids=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.junk.ids",
        ovlp_ids=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.ovlp.ids",
        haplotig_ids=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.haplotig.ids",
        repeat_ids=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.repeat.ids",
        highcov_ids=out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype, hap[^.]+}/{datatype, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.dups.highcov.ids"
    log:
        std=out_dir_path / "{assembly_stage}/{parameters}/log/qc_extract_stats_from_purge_dups_file.{parameters}.{genome_prefix}.purge_dups.{haplotype}.{datatype}.log",
        cluster_log=out_dir_path / "{assembly_stage}/{parameters}/log/qc_extract_stats_from_purge_dups_file.{parameters}.{genome_prefix}.purge_dups.{haplotype}.{datatype}.cluster.log",
        cluster_err=out_dir_path / "{assembly_stage}/{parameters}/log/qc_extract_stats_from_purge_dups_file.{parameters}.{genome_prefix}.purge_dups.{haplotype}.{datatype}.cluster.err"
    benchmark:
        out_dir_path / "{assembly_stage}/{parameters}/log/qc_extract_stats_from_purge_dups_file.{parameters}.{genome_prefix}.purge_dups.{haplotype}.{datatype}.benchmark.txt"
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
        artefact_fasta=out_dir_path  / "purge_dups/{prev_stage_parameters, [^/]+}..{purge_dups_parameters, [^/]+}/{purge_stage, [^/]+}/{haplotype, [^.]+}/{genome_prefix, [^/]+}.dups.{artefact, [^/]+}.fasta"
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

