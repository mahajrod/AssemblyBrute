rule get_seq_len:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta.fai",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        len_file="{fasta_dir}/{fasta_prefix}.len",
    log:
        std="{fasta_dir}/log/{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_seq_len"),
        cpus=parameters["threads"]["get_seq_len"],
        time=parameters["time"]["get_seq_len"],
        mem=parameters["memory_mb"]["get_seq_len"],
    threads:
        parameters["threads"]["get_seq_len"]
    shell:
         " workflow/scripts/get_sequence_lengths_from_fai.py -i {input.fasta} -o {output.len_file} >{log.std} 2>&1;"
