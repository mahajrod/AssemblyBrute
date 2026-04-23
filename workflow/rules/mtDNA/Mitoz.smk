rule downsample_pe_reads:
    input:
        forward_reads=output_dict["data"] / ("fastq/{datatype}/filtered/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_reads=output_dict["data"] / ("fastq/{datatype}/filtered/{pairprefix}_2%s" % config["fastq_extension"]),
        forward_fastqc=output_dict["qc"] / "fastqc/{datatype}/filtered/{pairprefix}_1_fastqc.zip",
        reverse_fastqc=output_dict["qc"] / "fastqc/{datatype}/filtered/{pairprefix}_2_fastqc.zip",
    output:
        forward_reads=output_dict["data"] / ("fastq/{datatype, hic|illumina}/downsampled_mitoz/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_reads=output_dict["data"] / ("fastq/{datatype, hic|illumina}/downsampled_mitoz/{pairprefix}_2%s" % config["fastq_extension"]),
    params:
        max_data_gbp=0.5, # TODO set as option for mitoz
    log:
        log=output_dict["log"]  / "downsample_pe_reads.{datatype}.{pairprefix}.log",
        cluster_log=output_dict["cluster_log"] / "downsample_pe_reads.{datatype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "downsample_pe_reads.{datatype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "downsample_pe_reads.{datatype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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

rule mitoz:
    input:
        forward_reads=output_dict["data"] / ("fastq/{datatype}/filtered/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_reads=output_dict["data"] / ("fastq/{datatype}/filtered/{pairprefix}_2%s" % config["fastq_extension"]),

        #forward_reads=output_dict["data"] / ("fastq/{datatype, hic|illumina}/downsampled_mitoz/{pairprefix}_1%s" % config["fastq_extension"]),
        #reverse_reads=output_dict["data"] / ("fastq/{datatype, hic|illumina}/downsampled_mitoz/{pairprefix}_2%s" % config["fastq_extension"]),
    output:
        finish_flag=out_dir_path / "mtDNA/mitoz/denovo/{datatype, hic|illumina}/{stage, filtered}/{pairprefix}/FINISH_FLAG"

    params:
        genome_prefix=config["genome_prefix"],
        clade=config["mitoz_clade"],
        species_name=config["species_name"],
        genetic_code=config["mtdna_genetic_code"],
        assembler="megahit",
        max_raw_data=3, # Gbp
        max_filtered_data=0, # Gbp
    log:
        log=(output_dict["log"]  / "mitoz.{datatype}.{stage}.{pairprefix}.log").resolve(),
        cluster_log=output_dict["cluster_log"] / "mitoz.{datatype}.{stage}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "mitoz.{datatype}.{stage}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "mitoz.{datatype}.{stage}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["mitoz"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["mitoz"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("mitoz"),
        cpus=parameters["threads"]["mitoz"],
        time=parameters["time"]["mitoz"],
        mem=parameters["memory_mb"]["mitoz"]
    threads: parameters["threads"]["mitoz"]
    shell:
        " set +e; "
        " WORKDIR=`dirname {output.finish_flag}`; "
        " TMPDIR=${{WORKDIR}}/tmp; "
        " mkdir -p ${{TMPDIR}}; "
        " OUTPUT_PREFIX={params.genome_prefix}.mtdna.{wildcards.datatype}.mitoz; "
        " FINAL_FASTA=${{WORKDIR}}/{params.genome_prefix}.mtdna.mitoz.denovo.{wildcards.datatype}.{wildcards.pairprefix}.{wildcards.stage}.fasta; "
        " mitoz all --workdir ${{WORKDIR}} --thread_number {threads} --assembler {params.assembler} "
        "                --tmp_dir ${{TMPDIR}} "
        "                --fq1 {input.forward_reads} --fq2 {input.reverse_reads} "
        "                --outprefix ${{OUTPUT_PREFIX}} "
        "                --clade {params.clade} "
        "                --requiring_taxa {params.clade} "
        "                --species_name '{params.species_name}' "
        "                --genetic_code {params.genetic_code} "
        "                --data_size_for_mt_assembly {params.max_raw_data},{params.max_filtered_data} > {log.log} 2>&1 || true; "
        " > {output.finish_flag}; "
        " if [[ -f '${{WORKDIR}}/${{OUTPUT_PREFIX}}.result/${{OUTPUT_PREFIX}}.megahit.result/${{OUTPUT_PREFIX}}.megahit.mitogenome.fa' ]]; "
        " then"
        "     cp ${{WORKDIR}}/${{OUTPUT_PREFIX}}.result/${{OUTPUT_PREFIX}}.megahit.result/${{OUTPUT_PREFIX}}.megahit.mitogenome.fa ${{FINAL_FASTA}}; "
        " fi; "
        " exit 0; "

