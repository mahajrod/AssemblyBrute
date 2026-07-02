rule simulate_pcr:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        nhr="{fasta_dir}/blast2.2.27/{fasta_prefix}.fasta.nhr"
    output:
        #primers="{fasta_dir}/simulate_pcr/{primers}.pair.{fasta_prefix}.fasta.amplicons",
        amplicons="{fasta_dir}/simulate_pcr/{primer_prefix}/{primer_prefix}.fasta.pair.{fasta_prefix}.fasta.amplicons",
        nhr="{fasta_dir}/simulate_pcr/{primer_prefix}/{fasta_prefix}.fasta.nhr",
    params:
        min_amplicon_len=,
        max_amplicon_len=,
        mix_primers=,
        max_mismatches=,
        three_prime_fixed_bases=
    log:
        std="{fasta_dir}/blast2.2.27.{fasta_prefix}.{primer_prefix}.log",
        cluster_log="{fasta_dir}/blast2.2.27.{fasta_prefix}.{primer_prefix}.cluster.log",
        cluster_err="{fasta_dir}/blast2.2.27.{fasta_prefix}.{primer_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/blast2.2.27.{fasta_prefix}.{primer_prefix}.benchmark.txt"
    conda:
        config["conda"]["simulate_pcr"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["simulate_pcr"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("simulate_pcr"),
        cpus=parameters["threads"]["simulate_pcr"] ,
        time=parameters["time"]["simulate_pcr"],
        mem=parameters["memory_mb"]["simulate_pcr"]
    threads: parameters["threads"]["simulate_pcr"]

    shell:
        " BLASTDB={input.nhr}; "
        " BLASTDB=${{BLASTDB%.nhr}}; "
        " WORKDIR=`basename {output.amplicons}`; "
        " workflow/external_tools/simulate_pcr/simulate_PCR -num_threads {threads} -db ${{BLASTDB}} "
        "        -primers primers.fasta -minlen {params.min_amplicon_len} -maxlen {params.max_amplicon_len} "
        "        -mux {params.mix_primers}  -extract_amp 1 -mm {params.max_mismatches} "
        "        -3prime {params.three_prime_fixed_bases} > {log.std} 2>&1; "