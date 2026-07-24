
rule miniprot:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        miniprot_gff="{fasta_dir}/{fasta_prefix}.miniprot.gff",
        candidate_tsv="{fasta_dir}/{fasta_prefix}.candidates.microchromosomes.tsv"
    params:
        microchromosome_prot_set=config["microchromosome_prot_set"]
    log:
        miniprot="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.miniprot.log",
        awk="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.awk.log",
        grep="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.grep.log",
        cut="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.cut.log",
        tr="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.tr.log",
        cut2="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.cut2.log",
        sed="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.sed.log",
        awk2="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.awk2.log",
        cut3="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.cut3.log",
        sort="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.sort.log",
        uniq="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.uniq.log",
        sort2="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.sort2.log",
        awk3="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.awk3.log",
        cluster_log="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.cluster.log",
        cluster_err="{fasta_dir}/log/miniprot.{fasta_prefix}.combined.cluster.err"
    benchmark:
        "{fasta_dir}/log/miniprot.{fasta_prefix}.combined.benchmark.txt"
    conda:
        config["conda"]["microchromosomes"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["microchromosomes"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("miniprot_microchromosomes"),
        cpus=parameters["threads"]["miniprot_microchromosomes"],
        time=parameters["time"]["miniprot_microchromosomes"],
        mem=parameters["memory_mb"]["miniprot_microchromosomes"]
    threads: parameters["threads"]["miniprot_microchromosomes"]
    shell:
        " miniprot -t {threads} --gff {input.fasta} {params.microchromosome_prot_set} > {output.miniprot_gff} 2>{log.miniprot}; "
        " awk '$3== \"mRNA\"' {output.miniprot_gff} 2>{log.awk} | grep -w \"Rank=1\" 2>{log.grep} | "
        "     cut -f1,9 2>{log.cut} | tr \";\" \"\\t\" 2>{log.tr} | cut -f1,4,5,6,7  2>{log.cut2} | "
        "     sed 's/Identity=//g;s/Positive=//g' 2>{log.sed} | awk '$2 >= 0.7' 2>{log.awk2} |  cut -f1 2>{log.cut3} | "
        "     sort 2>{log.sort} | uniq -c 2>{log.uniq} | sort -k1,1nr 2>{log.sort2} | awk '{{print $2 \"\\t\" $1}}' > {output.candidate_tsv} 2>{log.awk3} "

rule place_microchromosomes_first:
    input:
        candidate_tsv="{fasta_dir}/{fasta_prefix_part}.combined.candidates.microchromosomes.tsv",
        len_file="{fasta_dir}/{fasta_prefix_part}.combined.len",
        fasta="{fasta_dir}/{fasta_prefix_part}.combined.fasta",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        filtered_tsv="{fasta_dir}/{fasta_prefix_part}.reordered.candidates.microchromosomes.filtered.tsv",
        reordered_fasta="{fasta_dir}/{fasta_prefix_part}.reordered.fasta",
    params:
        max_length=parameters["tool_options"]["microchromosome_detection"]["max_length"]
    log:
        log="{fasta_dir}/log/place_microchromosomes_first.{fasta_prefix_part}.combined.log",
        cluster_log="{fasta_dir}/log/place_microchromosomes_first.{fasta_prefix_part}.combined.cluster.log",
        cluster_err="{fasta_dir}/log/place_microchromosomes_first.{fasta_prefix_part}.combined.cluster.err"
    benchmark:
        "{fasta_dir}/log/place_microchromosomes_first.{fasta_prefix_part}.combined.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("place_microchromosomes_first"),
        cpus=parameters["threads"]["place_microchromosomes_first"],
        time=parameters["time"]["place_microchromosomes_first"],
        mem=parameters["memory_mb"]["place_microchromosomes_first"]
    threads: parameters["threads"]["place_microchromosomes_first"]
    shell:
        " OUTPUT_PREFIX={output.reordered_fasta}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.reordered.fasta}}; "
        " workflow/scripts/curation/move_microchromosomes_first.py  -i {input.candidate_tsv}  " 
        "     -l {input.len_file} -f {input.fasta} -m {params.max_length} -o ${{OUTPUT_PREFIX}} > {log.log} 2>&1; "
