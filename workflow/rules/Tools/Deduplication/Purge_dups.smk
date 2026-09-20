localrules: extract_stats_from_purge_dups_file
#ruleorder: purge_dups > bam2bed
#ruleorder: purge_dups > bam2bed_for_hic_map

rule minimap2_purge_dups_reads:
    input:
        fastq=lambda wildcards: config["out_dir"] / "data/{0}/final/{1}{2}".format(wildcards.dtype,
                                                                                         wildcards.fileprefix,
                                                                                         config["data"][wildcards.dtype]["conv_ext"]),
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        paf="{fasta_dir}/{fasta_prefix}/purge_dups/{dtype}/{fasta_prefix}.{fileprefix}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.dtype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.dtype], " -x "),
    log:
        std="{fasta_dir}/log/minimap2_purge_dups_reads.{fasta_prefix}.{dtype}.{fileprefix}.log",
        gzip="{fasta_dir}/log/minimap2_purge_dups_reads.{fasta_prefix}.{dtype}.{fileprefix}.gzip.log",
        cluster_log="{fasta_dir}/log/minimap2_purge_dups_reads.{fasta_prefix}.{dtype}.{fileprefix}.cluster.log",
        cluster_err="{fasta_dir}/log/minimap2_purge_dups_reads.{fasta_prefix}.{dtype}.{fileprefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/minimap2_purge_dups_reads.{fasta_prefix}.{dtype}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("minimap2_purge_dups_reads"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["minimap2"], mode="linear")
    threads: parameters["threads"]["minimap2"]

    shell:
        " minimap2 {params.alignment_scheme} {params.index_size} -t {threads} {input.fasta} "
        "     {input.fastq} 2>{log.std} | gzip -c - > {output.paf} 2>{log.gzip}; "

def get_paf_list(wildcards):
    paf_list = []
    for dtype in wildcards.datatype.split("_"):
        paf_list += expand(rules.minimap2_purge_dups_reads.output.paf,
                                         dtype=[dtype],
                                         fileprefix=config["data"][dtype]["conv_file_prefix_list"],
                                         allow_missing=True)
    return paf_list

rule get_purge_dups_read_stat:
    input:
        paf=get_paf_list,
        genomescope_report=config["out_dir"] / "kmer/{0}/final/genomescope/{1}.{0}.final.{2}.{3}.genomescope.parameters".format("_".join(config["final_kmer_datatypes"]),
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"]),
        log_dir=ancient("{fasta_dir}/log/")
    output:
        pbstat="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/PB.stat",
        pbbasecov="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/PB.base.cov",
        cutoffs="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/cutoffs",
        len="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/PB.base.cov.len",
        stat="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/PB.base.cov.stat",
        bed="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/PB.base.cov.bed"
    params:
        cov_multiplicator=parameters["tool_options"]["assembly_qc"]["purge_dups"]["cov_multiplicator"],
        calcuts_lower_threshold=parse_option("lower_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -l "),
        calcuts_haploid_diploid_threshold=parse_option("haploid_diploid_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -m "),
        calcuts_upper_threshold=str(config["tool_manually_adjusted_features"]["calcuts"]["upper_threshold"]), # None needs to be converted to "None"
    log:
        pbstat="{fasta_dir}/log/get_purge_dups_read_stat.{fasta_prefix}.{datatype}.pbstat.log",
        calcuts="{fasta_dir}/log/get_purge_dups_read_stat.{fasta_prefix}.{datatype}.calcuts.log",
        convert="{fasta_dir}/log/get_purge_dups_read_stat.{fasta_prefix}.{datatype}.convert.log",
        cluster_log="{fasta_dir}/log/get_purge_dups_read_stat.{fasta_prefix}.{datatype}.cluster.log",
        cluster_err="{fasta_dir}/log/get_purge_dups_read_stat.{fasta_prefix}.{datatype}.cluster.err"
    benchmark:
        "{fasta_dir}/log/minimap2_purge_dups_reads.{fasta_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
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
        "     -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts}; " #check parameters for calcuts
        " convert_coverage_file_to_bed.py -i {output.pbbasecov}  -o ${{LEN_FILE%.len}} > {log.convert} 2>&1; "

rule minimap2_purge_dups_assembly:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        split_assembly="{fasta_dir}/{fasta_prefix}/purge_dups/{fasta_prefix}.split.fasta",
        paf="{fasta_dir}/{fasta_prefix}/purge_dups/{fasta_prefix}.split.minimap2.self.paf.gz"
    params:
        index_size=parse_option("index_size", parameters["tool_options"]["minimap2"]["self"], " -I "),
        alignment_scheme=parse_option("alignment_scheme", parameters["tool_options"]["minimap2"]["self"], " -x "),
    log:
        split_fa="{fasta_dir}/log/minimap2_purge_dups_assembly.{fasta_prefix}.split_fa.log",
        minimap2="{fasta_dir}/log/minimap2_purge_dups_assembly.{fasta_prefix}.minimap2.log",
        cluster_log="{fasta_dir}/log/minimap2_purge_dups_assembly.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/minimap2_purge_dups_assembly.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/minimap2_purge_dups_assembly.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("minimap2_purge_dups_assembly"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["minimap2"], mode="linear")
    threads: parameters["threads"]["minimap2"]
    shell:
        " split_fa {input.fasta} > {output.split_assembly} 2>{log.split_fa};"
        " minimap2 -DP {params.alignment_scheme} {params.index_size} -t {threads} {output.split_assembly} "
        "     {output.split_assembly} 2>{log.minimap2} |  gzip -c - > {output.paf}; "

rule purge_dups: #
    input:
        cutoffs=rules.get_purge_dups_read_stat.output.cutoffs,
        pbbasecov=rules.get_purge_dups_read_stat.output.pbbasecov,
        self_paf=rules.minimap2_purge_dups_assembly.output.paf,
        log_dir="{fasta_dir}/log/",
    output:
        bed="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.raw.bed",
    log:
        purge_dups="{fasta_dir}/log/qc_purge_dups.{fasta_prefix}.{datatype}.purge_dups.log",
        cluster_log="{fasta_dir}/log/qc_purge_dups.{fasta_prefix}.{datatype}.cluster.log",
        cluster_err="{fasta_dir}/log/qc_purge_dups.{fasta_prefix}.{datatype}.cluster.err"
    benchmark:
        "{fasta_dir}/log/qc_purge_dups.{fasta_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("purge_dups"),
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]
    shell:
        " purge_dups -2 -T {input.cutoffs} -c {input.pbbasecov} {input.self_paf} > {output.bed} 2>{log.purge_dups}; "


rule extract_stats_from_purge_dups_file:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        raw_dups_bed=rules.purge_dups.output.bed, #"{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.raw.bed",
        stat="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/PB.base.cov.stat",
        len="{fasta_dir}/{fasta_prefix}.len",
        log_dir="{fasta_dir}/log/",
    output:
        extended_bed="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.extended.bed",
        stat="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.stat",
        id_files=expand("{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.{artefact_type}.ids",
                        artefact_type=["junk", "ovlp", "haplotig", "repeat", "highcov"],
                        allow_missing=True),
        bed_files=expand("{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.{artefact_type}.extended.bed",
                        artefact_type=["junk", "ovlp", "haplotig", "repeat", "highcov"],
                        allow_missing=True)
    log:
        std="{fasta_dir}/log/extract_stats_from_purge_dups_file.{fasta_prefix}.purge_dups.{datatype}.log",
        cluster_log="{fasta_dir}/log/extract_stats_from_purge_dups_file.{fasta_prefix}.purge_dups.{datatype}.cluster.log",
        cluster_err="{fasta_dir}/log/extract_stats_from_purge_dups_file.{fasta_prefix}.purge_dups.{datatype}.cluster.err"
    benchmark:
        "{fasta_dir}/log/extract_stats_from_purge_dups_file.{fasta_prefix}.purge_dups.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("extract_stats_from_purge_dups_file"),
        cpus=parameters["threads"]["extract_stats_from_purge_dups_file"] ,
        time=parameters["time"]["extract_stats_from_purge_dups_file"],
        mem=parameters["memory_mb"]["extract_stats_from_purge_dups_file"]
    threads: parameters["threads"]["extract_stats_from_purge_dups_file"]

    shell:
        " STATS_FILE={output.stat}; "
        " ./workflow/scripts/purge_dups/calculate_purge_dups_stats.py  -b {input.raw_dups_bed} -s {input.stat} -l {input.len} "
        "     -o ${{STATS_FILE%.stat}} > {log.std} 2>&1; "

rule create_purge_dups_track:
    priority: 500
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir="{fasta_dir}/log/",
        extended_bed=rules.extract_stats_from_purge_dups_file.output.extended_bed, #"{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.dups.extended.bed",
    output:
        bedgraph="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.purge_dups.{datatype}.{artefact_type}.track.bedgraph",
        bed="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.purge_dups.{datatype}.{artefact_type}.track.bed",
    log:
        log="{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.log",
        cluster_log="{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.cluster.log",
        cluster_err="{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.cluster.err"
    params:
        artefact=lambda wildcards: f"{wildcards.artefact_type}".upper()
    benchmark:
        "{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco5_intersect_all"),
        cpus=parameters["threads"]["busco5_intersect_all"],
        time=parameters["time"]["busco5_intersect_all"],
        mem=parameters["memory_mb"]["busco5_intersect_all"],
    threads:
        parameters["threads"]["busco5_intersect_all"]
    shell: # true was added to grep as it return 1 if no matching lines were found - snakemake recognizes it as error
        " (tail -n +2 {input.extended_bed} | grep {params.artefact} || true) |  awk -F '\\t' '{{ printf \"%s\\t%i\\t%i\\t1\\n\",$1,$2,$3 }}' | "
        "     sort -k1,1V -k2,2n -k3,3n  > {output.bedgraph} 2>{log.log}; "
        " head -n 1 {input.extended_bed} > {output.bed} 2>> {log.log}; "
        " (tail -n +2 {input.extended_bed} | grep {params.artefact} || true) >> {output.bed} 2>>{log.log}; "
