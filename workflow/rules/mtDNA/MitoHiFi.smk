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
        " cp -f ${{OUT_DIR}}/${{REF_PREFIX}}.gb {output.mtdna_ref_gb}; "


rule mitohifi_reads:
    input:
        mtdna_ref_fasta=out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.fasta",#lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        mtdna_ref_gb=out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.gb",#lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        #hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #            fileprefix=input_file_prefix_dict["hifi"],
        #            allow_missing=True),
        hifi_reads=output_dict["data"] / ("fastq/hifi/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        finish_flag=out_dir_path / "mtDNA/mitohifi/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/FINISH_FLAG",
        #mtDNA_gb=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.gb",
        #mtDNA_fasta=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.fasta",
        #coverage_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.coverage.png",
        #annotation_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, filtered}/{fileprefix}/final_mitogenome.annotation.png"
    params:
        sif=config["tool_containers"]["mitohifi"],
        kingdom=config["kingdom"],
        genetic_code=config["mtdna_genetic_code"],
        min_mapping_quality=parameters["tool_options"]["mitohifi"]["min_mapping_quality"],
        genome_prefix=config["genome_prefix"]
    log:
        cp=output_dict["log"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.cp.log",
        cd=output_dict["log"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.cd.log",
        mitohifi=output_dict["log"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.mitohifi.log",
        cluster_log=output_dict["cluster_log"] / "mitohifi_reads.{stage}.{mtdna_ref}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "mitohifi_reads.{stage}.{mtdna_ref}.{fileprefix}.err"
    benchmark:
        output_dict["benchmark"]  / "mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("mitohifi_reads"),
        cpus=parameters["threads"]["mitohifi_reads"],
        time=parameters["time"]["mitohifi_reads"],
        mem=parameters["memory_mb"]["mitohifi_reads"]
    threads: parameters["threads"]["mitohifi_reads"]
    shell:
        " OUT_DIR=`realpath -m {output.finish_flag}`; "
        " OUT_DIR=`dirname ${{OUT_DIR}}`; "
        " HIFI_READS=`realpath {input.hifi_reads}`; "
        " HIFI_DIR=`dirname ${{HIFI_READS}}`; "
        " REF_FASTA=`realpath {input.mtdna_ref_fasta}`; "
        " REF_GB=`realpath {input.mtdna_ref_gb}`; "
        " REF_DIR=`dirname ${{REF_FASTA}}`; "
        " MITOHIFI_LOG=`realpath {log.mitohifi}`; "
        " OUTPUT_PREFIX={params.genome_prefix}.mtdna.mitohifi.ref_{wildcards.mtdna_ref}.hifi.{wildcards.fileprefix}.{wildcards.stage}; "
        " cp -f {input.mtdna_ref_fasta} {input.mtdna_ref_gb} ${{OUT_DIR}} > {log.cp} 2>&1; "
        " cd ${{OUT_DIR}} > {log.cd} 2>&1; "
        " singularity run --pid "
        "                 --bind ${{OUT_DIR}}:${{OUT_DIR}} "
        "                 --bind ${{HIFI_DIR}}:${{HIFI_DIR}} "
        "                 --bind ${{REF_DIR}}:${{REF_DIR}} "
        "                 {params.sif} mitohifi.py "
        "                 -r ${{HIFI_READS}} -f ${{REF_FASTA}} -g ${{REF_GB}} "
        "                 -t {threads} -a {params.kingdom} -covMap {params.min_mapping_quality} "
        "                 -o {params.genetic_code} > ${{MITOHIFI_LOG}} 2>&1 || true; "
        " > FINISH_FLAG; "
        " if [[ -f 'final_mitogenome.fasta' ]]; "
        " then "
        "       cp final_mitogenome.fasta ${{OUTPUT_PREFIX}}.fasta; "
        "       cp final_mitogenome.gb ${{OUTPUT_PREFIX}}.gb; "
        "       cp contigs_stats.tsv ${{OUTPUT_PREFIX}}.contigs_stats.tsv; "
        "       cp final_mitogenome.annotation.png ${{OUTPUT_PREFIX}}.annotation.png; "
        "       cp final_mitogenome.coverage.png ${{OUTPUT_PREFIX}}.coverage.png; "        
        " fi; "

use rule mitohifi_reads as mitohifi_combined_reads with:
    input:
        mtdna_ref_fasta=out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.fasta",
        mtdna_ref_gb = out_dir_path / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.gb",
        hifi_reads = output_dict["data"] / ("fastq/hifi/combined/hifi.combined%s" % config["fastq_extension"])
    output:
        stats=out_dir_path / "mtDNA/mitohifi/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/contigs_stats.tsv",
        #mtDNA_gb=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.gb",
        #mtDNA_fasta=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.fasta",
        #coverage_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.coverage.png",
        #annotation_plot=out_dir_path / "mtDNA/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/final_mitogenome.annotation.png"