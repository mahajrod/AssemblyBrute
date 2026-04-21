rule get_recommended_mtDNA_reference: # reference is inferred by NCBI taxid
    input:
        []
    output:
        mtdna_ref_fasta=out_dir_path / "data/mtDNA/recommended/recommended.fasta",
        mtdna_ref_gb=out_dir_path / "data/mtDNA/recommended/recommended.gb",
    params:
        latin_name=config["species_name"]
    log:
        log=output_dict["log"]  / "get_recommended_mtDNA_reference.log",
        cluster_log=output_dict["cluster_log"] / "get_recommended_mtDNA_reference.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_recommended_mtDNA_reference.err"
    benchmark:
        output_dict["benchmark"]  / "get_recommended_mtDNA_reference.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_recommended_mtDNA_reference"),
        cpus=parameters["threads"]["get_recommended_mtDNA_reference"],
        time=parameters["time"]["get_recommended_mtDNA_reference"],
        mem=parameters["memory_mb"]["get_recommended_mtDNA_reference"]
    threads: parameters["threads"]["get_recommended_mtDNA_reference"]
    shell:
        " OUT_DIR=`dirname {output.mtdna_ref_fasta}`; "
        " workflow/external_tools/mitohifi/src/findMitoReference.py "
        "      --species '{params.latin_name}' "
        "      --outfolder ${{OUT_DIR}} "
        "      --type mitochondrion > {log.log} 2>&1; "
        " REF_PREFIX=`grep 'output is written' {log.log} | cut -d ' ' -f 5 `;"
        " REF_PREFIX=`basename ${{REF_PREFIX}} | sed 's/\.\[gb,fasta\]//'`; "
        " cp -f ${{OUT_DIR}}/${{REF_PREFIX}}.fasta {output.mtdna_ref_fasta}; "
        " cp -f ${{OUT_DIR}}/${{REF_PREFIX}}.gb {output.mtdna_ref_fasta}; "


rule mitohifi_reads:
    input:
        mtdna_ref_fasta=out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.fasta",#lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        mtdna_ref_gb=out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.gb",#lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        #hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #            fileprefix=input_file_prefix_dict["hifi"],
        #            allow_missing=True),
        hifi_reads=output_dict["data"] / ("fastq/hifi/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        stats=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/contig_stats.tsv",
        mtDNA_gb=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.gb",
        mtDNA_fasta=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.fasta",
        coverage_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.coverage.png",
        annotation_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.annotation.png"
    params:
        sif=config["tool_containers"]["mitohifi"],
        kingdom=config["kingdom"],
        genetic_code=config["mtdna_genetic_code"],
        min_mapping_quality=parameters["tool_options"]["mitohifi"]["min_mapping_quality"]
    log:
        cp=output_dict["log"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.cp.log",
        cd=output_dict["log"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.cd.log",
        mitohifi=output_dict["log"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.mitohifi.log",
        cluster_log=output_dict["cluster_log"] / "mitohifi_reads.{stage}.{mtdna_ref}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "mitohifi_reads.{stage}.{mtdna_ref}.{fileprefix}.err"
    benchmark:
        output_dict["benchmark"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("mitohifi_reads"),
        cpus=parameters["threads"]["mitohifi_reads"],
        time=parameters["time"]["mitohifi_reads"],
        mem=parameters["memory_mb"]["mitohifi_reads"]
    threads: parameters["threads"]["mitohifi_reads"]
    shell:
        " OUT_DIR=`realpath -m {output.stats}`; "
        " OUT_DIR=`dirname ${{OUT_DIR}}`; "
        " HIFI_READS=`realpath {input.hifi_reads}`; "
        " HIFI_DIR=`dirname ${{HIFI_READS}}`; "
        " REF_FASTA=`realpath {input.mtdna_ref_fasta}`; "
        " REF_GB=`realpath {input.mtdna_ref_gb}`; "
        " REF_DIR=`dirname ${{REF_FASTA}}`; "
        " cp -f {input.mtdna_ref_fasta} {input.mtdna_ref_gb} ${{OUT_DIR}} > {log.cp} 2>&1; "
        " cd ${{OUT_DIR}} > {log.cd} 2>&1; "
        " singularity run --pid --contain --pid --contain "
        "                 --bind ${{OUT_DIR}}:${{OUT_DIR}} "
        "                 --bind ${{HIFI_DIR}}:${{HIFI_DIR}} "
        "                 --bind ${{REF_DIR}}:${{REF_DIR}} "
        "                 {params.sif} mitohifi.py "
        "                 -r ${{HIFI_READS}} -f ${{REF_FASTA}} -g ${{REF_GB}} "
        "                 -t {threads} -a {params.kingdom} -covMap {params.min_mapping_quality} "
        "                 -o {params.genetic_code} > {log.mitohifi} 2>&1; "

rule combine_long_reads:
    input:
        long_reads=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (wildcards.datatype,
                                                                                                         config["fastq_extension"])),
                          fileprefix=input_file_prefix_dict[wildcards.datatype])
    output:
        combined_long_reads=output_dict["data"] / ("fastq/{datatype, hifi|nanopore|simplex|duplex}/combined/{datatype}.combined%s" % config["fastq_extension"])
    log:
        log=output_dict["log"] / "combine_long_reads.{datatype}.log",
        cluster_log=output_dict["cluster_log"] / "combine_long_reads.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "combine_long_reads.{datatype}.err"
    benchmark:
        output_dict["benchmark"] / "combine_long_reads.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" %config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("combine_long_reads"),
        cpus=parameters["threads"]["combine_long_reads"],
        time=parameters["time"]["combine_long_reads"],
        mem=parameters["memory_mb"]["combine_long_reads"],
    threads: parameters["threads"]["combine_long_reads"]
    shell:
        "cat {input.long_reads} > {output.combined_long_reads} 2>{log.log}; "

use rule mitohifi_reads as mitohifi_combined_reads with:
    input:
        mtdna_ref_fasta=out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.fasta",
        mtdna_ref_gb = out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.gb",
        hifi_reads = output_dict["data"] / ("fastq/hifi/combined/hifi.combined%s" % config["fastq_extension"])
    output:
        stats=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/contig_stats.tsv",
        mtDNA_gb=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.gb",
        mtDNA_fasta=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.fasta",
        coverage_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.coverage.png",
        annotation_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.annotation.png"