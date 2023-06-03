rule get_seq_len:
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
    output:
        #dir=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/",
        len_file=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
    #params:
    #    MAVR_dir=config["MAVR_path"]
    log:
        std=output_dict["log"].resolve() / "get_seq_len.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=(output_dict["cluster_log"]).resolve() / "get_seq_len.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=(output_dict["cluster_error"]).resolve() / "get_seq_len.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "get_seq_len.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["get_seq_len"],
        time=parameters["time"]["get_seq_len"],
        mem=parameters["memory_mb"]["get_seq_len"],
    threads:
        parameters["threads"]["get_seq_len"]
    shell:
         " get_sequence_lengths.py -i {input.fasta} -o {output.len_file} 1>{log.std} 2>&1;" #{params.MAVR_dir}/scripts/sequence/
         #" OMP_NUM_THREADS={threads} merqury.sh {input.meryl_db_dir} "
         #" {input.primary_assembly} {input.alternative_assembly} {params.out_prefix}  1>{log.std} 2>&1 || true;"


