rule get_recommended_mtDNA_reference: # reference is inferred by NCBI taxid
    input:
        []
    output:
        mtdna_ref_fasta=config["out_dir"]  / "data/mtDNA/recommended/recommended.fasta",
        mtdna_ref_gb=config["out_dir"]  / "data/mtDNA/recommended/recommended.gb",
    params:
        latin_name=config["species_name"]
    log:
        log=config["out_dir"] / "log/get_recommended_mtDNA_reference.log",
        cluster_log=config["out_dir"] / "log/get_recommended_mtDNA_reference.cluster.log",
        cluster_err=config["out_dir"] / "log/get_recommended_mtDNA_reference.err"
    benchmark:
        config["out_dir"] / "log/get_recommended_mtDNA_reference.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_recommended_mtDNA_reference"),
        cpus=parameters["threads"]["get_recommended_mtDNA_reference"],
        time=parameters["time"]["get_recommended_mtDNA_reference"],
        mem=parameters["memory_mb"]["get_recommended_mtDNA_reference"]
    threads: parameters["threads"]["get_recommended_mtDNA_reference"]
    shell:
        " OUT_DIR=`dirname {output.mtdna_ref_fasta}`; "
        " workflow/external_tools/mitohifi/src/findMitoReference.py "
        "     --species '{params.latin_name}' "
        "     --outfolder ${{OUT_DIR}} "
        "     --type mitochondrion > {log.log} 2>&1; "
        " REF_PREFIX=`grep 'output is written' {log.log} | cut -d ' ' -f 5 `;"
        " REF_PREFIX=`basename ${{REF_PREFIX}} | sed 's/\.\[gb,fasta\]//'`; "
        " cp -f ${{OUT_DIR}}/${{REF_PREFIX}}.fasta {output.mtdna_ref_fasta}; "
        " cp -f ${{OUT_DIR}}/${{REF_PREFIX}}.gb {output.mtdna_ref_gb}; "


rule mitohifi_reads:
    input:
        mtdna_ref_fasta=config["out_dir"]  / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.fasta",
        mtdna_ref_gb=config["out_dir"]  / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.gb",
        hifi_reads=lambda wildcards: config["out_dir"] / ("data/hifi/%s/%s%s" % (wildcards.stage,
                                                                                 wildcards.fileprefix,
                                                                                 config["data"]["hifi"]["conv_ext"]))
    output:
        finish_flag=config["out_dir"]  / "mtDNA/mitohifi/{mtdna_ref}/hifi/{stage, final}/{fileprefix}/FINISH_FLAG",
    params:
        singularity_load_str=(config["singularity_load_str"] + "; ") if config["singularity_load_mode"] == "cluster" else "",
        sif=config["tool_containers"]["mitohifi"],
        kingdom=config["kingdom"],
        genetic_code=config["mtdna_genetic_code"],
        min_mapping_quality=parameters["tool_options"]["mitohifi"]["min_mapping_quality"],
        genome_prefix=config["genome_prefix"]
    log:
        cp=(config["out_dir"] / "log/mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.cp.log").resolve(),
        cd=(config["out_dir"] / "log/mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.cd.log").resolve(),
        sed=(config["out_dir"] / "log/mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.sed.log").resolve(),
        mitohifi=(config["out_dir"] / "log/mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.mitohifi.log").resolve(),
        cluster_log=config["out_dir"] / "log/mitohifi_reads.{stage}.{mtdna_ref}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/mitohifi_reads.{stage}.{mtdna_ref}.{fileprefix}.err"
    benchmark:
        config["out_dir"] / "log/mitohifi_reads.{mtdna_ref}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        singularity_conda_env
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("mitohifi_reads"),
        cpus=parameters["threads"]["mitohifi_reads"],
        time=parameters["time"]["mitohifi_reads"],
        mem=parameters["memory_mb"]["mitohifi_reads"]
    threads: parameters["threads"]["mitohifi_reads"]
    shell:
        " {params.singularity_load_str} "
        " OUT_DIR=`realpath -m {output.finish_flag}`; "
        " OUT_DIR=`dirname ${{OUT_DIR}}`; "
        " HIFI_READS=`realpath {input.hifi_reads}`; "
        " HIFI_DIR=`dirname ${{HIFI_READS}}`; "
        " REF_FASTA=`realpath {input.mtdna_ref_fasta}`; "
        " REF_GB=`realpath {input.mtdna_ref_gb}`; "
        " REF_DIR=`dirname ${{REF_FASTA}}`; "
        " OUTPUT_PREFIX={params.genome_prefix}.mtdna.mitohifi.ref_{wildcards.mtdna_ref}.hifi.{wildcards.fileprefix}.{wildcards.stage}; "
        " cp -f {input.mtdna_ref_fasta} {input.mtdna_ref_gb} ${{OUT_DIR}} > {log.cp} 2>&1; "
        " cd ${{OUT_DIR}} > {log.cd} 2>&1; "
        " singularity run --pid "
        "     --bind ${{OUT_DIR}}:${{OUT_DIR}} "
        "     --bind ${{HIFI_DIR}}:${{HIFI_DIR}} "
        "     --bind ${{REF_DIR}}:${{REF_DIR}} "
        "     {params.sif} mitohifi.py "
        "     -r ${{HIFI_READS}} -f ${{REF_FASTA}} -g ${{REF_GB}} "
        "     -t {threads} -a {params.kingdom} -covMap {params.min_mapping_quality} "
        "     -o {params.genetic_code} > {log.mitohifi} 2>&1 || true; "
        " > FINISH_FLAG; "
        " if [[ -f \"final_mitogenome.fasta\" ]]; "
        " then "
        "     sed 's/^>/>{params.genome_prefix}.{wildcards.fileprefix} /' final_mitogenome.fasta > ${{OUTPUT_PREFIX}}.fasta 2>{log.sed}; "
        "     cp -f final_mitogenome.gb ${{OUTPUT_PREFIX}}.gb >> {log.cp} 2>&1; "
        "     cp -f contigs_stats.tsv ${{OUTPUT_PREFIX}}.contigs_stats.tsv >> {log.cp} 2>&1; " 
        "     cp -f final_mitogenome.annotation.png ${{OUTPUT_PREFIX}}.annotation.png >> {log.cp} 2>&1; " 
        "     cp -f final_mitogenome.coverage.png ${{OUTPUT_PREFIX}}.coverage.png >> {log.cp} 2>&1; "         
        " fi; "



use rule mitohifi_reads as mitohifi_combined_reads with:
    input:
        mtdna_ref_fasta=config["out_dir"] / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.fasta",
        mtdna_ref_gb=config["out_dir"] / "data/mtDNA/{mtdna_ref}/{mtdna_ref}.gb",
        hifi_reads=config["out_dir"]/ ("data/hifi/combined/hifi.combined%s" % config["data"]["hifi"]["conv_ext"])
    output:
        finish_flag=config["out_dir"] / "mtDNA/mitohifi/{mtdna_ref}/hifi/{stage, combined}/{fileprefix}/FINISH_FLAG",

"""        
rule downsample_se_reads:
    input:
        se_reads=config["out_dir"] / ("data/fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_ext"]),
        forward_fastqc=config["out_dir"] / "qc/fastqc/{datatype}/filtered/{pairprefix}_1_fastqc.zip",
        reverse_fastqc=config["out_dir"] / "qc/fastqc/{datatype}/filtered/{pairprefix}_2_fastqc.zip",
    output:
        forward_reads=config["out_dir"] / ("data/fastq/{datatype, hic|illumina}/downsampled_mitoz/{pairprefix}_1%s" % config["fastq_ext"]),
        reverse_reads=config["out_dir"] / ("data/fastq/{datatype, hic|illumina}/downsampled_mitoz/{pairprefix}_2%s" % config["fastq_ext"]),
    params:
        max_data_gbp=0.5, # TODO set as option for mitoz
    log:
        log=config["out_dir"] / "log/downsample_pe_reads.{datatype}.{pairprefix}.log",
        cluster_log=config["out_dir"] / "log/downsample_pe_reads.{datatype}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/downsample_pe_reads.{datatype}.{pairprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/downsample_pe_reads.{datatype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("downsample_pe_reads"),
        cpus=parameters["threads"]["downsample_pe_reads"],
        time=parameters["time"]["downsample_pe_reads"],
        mem=parameters["memory_mb"]["downsample_pe_reads"]
    threads: parameters["threads"]["downsample_pe_reads"]
    shell:
        " DATA=0; "
        " for FILE in {input.forward_fastqc} {input.reverse_fastqc}; "
        "   do "
        "   TMP=(`unzip -p ${{FILE}} | grep 'Total Bases' | cut -f 2`); "
        "   if [[ \"${{TMP[1]}}\" == 'Mbp' ]];"
        "      then "
        "      TMP=`echo \"scale=3; ${{TMP}} / 1000\" | bc | sed 's/^\./0./'`; "
        "   elif [[ \"${{TMP[1]}}\" != 'Gbp' ]];"
        "      then "
        "      echo \"Unrecognized data unit - ${{TMP[1]}} \"; "
        "      exit 1;"
        "   fi;"
        " DATA=`echo \"scale=3; ${{DATA}} + ${{TMP}}\" | bc | sed 's/^\./0./'`;  "
        " DOWNSAMPLING_FRACTION=`echo \"scale=3; {params.max_data_gbp} / ${{DATA}} \" | bc | sed 's/^\./0./'`;  "
        " if [[  $(echo \"${{DOWNSAMPLING_FRACTION}} < 0.8\" | bc -l) -eq 1 ]];"
        "    then "
        "    seqtk sample -2 -s1000 {input.forward_reads} ${{DOWNSAMPLING_FRACTION}} > {output.forward_reads} 2>{log.log}; "
        "    seqtk sample -2 -s1000 {input.reverse_reads} ${{DOWNSAMPLING_FRACTION}} > {output.reverse_reads} 2>>{log.log}; "
        "    else "
        "    cp -f {input.forward_reads} {output.forward_reads} >> {log.log} 2>&1; "
        "    cp -f {input.reverse_reads} {output.reverse_reads} >> {log.log} 2>&1; "
        " fi; "
"""