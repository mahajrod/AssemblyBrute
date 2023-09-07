rule get_seq_len:
    input:
        fasta="{fasta_prefix}.fasta.fai",
    output:
        #dir=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/",
        len_file="{fasta_prefix}.len",
    #params:
    #    MAVR_dir=config["MAVR_path"]
    log:
        std="{fasta_prefix}.log",
        cluster_log="{fasta_prefix}.cluster.log",
        cluster_err="{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_seq_len"),
        cpus=parameters["threads"]["get_seq_len"],
        time=parameters["time"]["get_seq_len"],
        mem=parameters["memory_mb"]["get_seq_len"],
    threads:
        parameters["threads"]["get_seq_len"]
    shell:
         " workflow/scripts/get_sequence_lengths_from_fai.py -i {input.fasta} -o {output.len_file} 1>{log.std} 2>&1;" #{params.MAVR_dir}/scripts/sequence/
         #" OMP_NUM_THREADS={threads} merqury.sh {input.meryl_db_dir} "
         #" {input.primary_assembly} {input.alternative_assembly} {params.out_prefix}  1>{log.std} 2>&1 || true;"


