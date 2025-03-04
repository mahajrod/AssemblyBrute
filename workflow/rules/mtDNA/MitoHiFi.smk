
rule mitohifi:
    input:
        mito_ref_fasta=lambda wildcards: input_reference_filedict[wildcards.ref_name]["mtdna.fasta"],#lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        mito_ref_gb=lambda wildcards: input_reference_filedict[wildcards.ref_name]["mtdna.gb"],#lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                    fileprefix=input_file_prefix_dict["hifi"],
                    allow_missing=True),
    output:
        out_dir=directory(out_dir_path / "mtDNA/{ref_name}/"),
        temp_merged_reads=temp(out_dir_path / "mtDNA/{ref_name}/all.fastq.gz"),
        mtDNA_gb=,
        mtDNA_fasta=
        final_mitogenome.fasta
        final_mitogenome.gb
        shared_genes.tsv
    params:
        sif=config["tool_containers"]["mitohifi"],
        kingdom=config["kingdom"],
        genetic_code=config["mitochondrial_genetic_code"],
        min_mapping_quality=parameters["tool_options"]["mitohifi"]["min_mapping_quality"]
    log:
        mkdir=output_dict["log"]  / "mitohifi.{ref_name}.mkdir.log",
        cat=output_dict["log"]  / "mitohifi.{ref_name}.cat.log",
        cp=output_dict["log"]  / "mitohifi.{ref_name}.cp.log",
        cd=output_dict["log"]  / "mitohifi.{ref_name}.cd.log",
        mitohifi=output_dict["log"]  / "mitohifi.{ref_name}.mitohifi.log",
        cluster_log=output_dict["cluster_log"] / "mitohifi.{ref_name}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "mitohifi.{ref_name}.err"
    benchmark:
        output_dict["benchmark"]  / "mitohifi.{ref_name}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("mitohifi"),
        cpus=parameters["threads"]["mitohifi"],
        time=parameters["time"]["mitohifi"],
        mem=parameters["memory_mb"]["mitohifi"]
    threads: parameters["threads"]["mitohifi"]
    shell:
        " OUT_DIR=`realpath {output.out_dir}`; "
        " mkdir -p {output.out_dir} > {log.mkdir} 2>&1; "
        " cat {input.hifi} > {output.temp_merged_reads} 2>{log.cat}; "
        " cp {input.mito_ref_fasta}` {output.out_dir}/reference.mtdna.fasta > {log.cp} 2>&1; "
        " cp {input.mito_ref_gb}` {output.out_dir}/reference.mtdna.gb >> {log.cp} 2>&1;; "
        " cd {output.out_dir} > {log.cd} 2>&1; "
        " singularity run --bind ${{OUT_DIR}}:${{OUT_DIR}} "
        " {params.sif} mitohifi.py -r `basename {output.temp_merged_reads}`"
        " -f reference.mtdna.fasta -g reference.mtdna.gb -t {threads} "
        " -a {params.kingdom} -covMap {params.min_mapping_quality} -o {params.genetic_code} > {log.mitohifi} 2>&1; "


REF_MTDNA_DIR="/maps/projects/tomg/people/xsg178/yggdrasil/mito/eudromia_elegans/assembly/"
    sbatch -n 30 -t "02:00:00" --mem 15000 --nodes 1 --wrap="singularity run --bind ${REF_MTDNA_DIR}:${REF_MTDNA_DIR} docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c bEudEle1.sanger.hap1.fasta -f ${REF_MTDNA_DIR}/eudromia_elegans.mtDNA.fasta  -g ${REF_MTDNA_DIR}/eudromia_elegans.mtDNA.gb  -t 30 -a animal -covMap 20 -o 2 "

    FASTQ_DIR=/projects/tomg/people/xsg178/yggdrasil/assembly/eudromia_elegans/input/hifi/fastq/
    REF_MTDNA_DIR="/maps/projects/tomg/people/xsg178/yggdrasil/mito/eudromia_elegans/assembly/"
    sbatch -n 30 -t "02:00:00" --mem 45000 --nodes 1 --wrap="singularity run --bind ${REF_MTDNA_DIR}:${REF_MTDNA_DIR} --bind ${FASTQ_DIR}:${FASTQ_DIR} docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r ${FASTQ_DIR}/*.fastq.gz -f ${REF_MTDNA_DIR}/NC_002772.2.fasta  -g ${REF_MTDNA_DIR}/NC_002772.2.gb  -t 30 -a animal -covMap 20 -o 2 "