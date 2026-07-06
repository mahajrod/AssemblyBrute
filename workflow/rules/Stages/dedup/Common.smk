localrules: create_dedup_input_links

use rule create_local_links as create_dedup_input_links with:
    input:
        input=lambda wildcards: config["out_dir"]  / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["dedup"].prev_stage,
                                                                                         get_prev_stage_parameters(wildcards.parameters),
                                                                                         wildcards.genome_prefix,
                                                                                         wildcards.haplotype),
        log_dir=ancient(config["out_dir"] / "dedup/{parameters}/log/")
    output:
        input=config["out_dir"] / "dedup/{parameters, [^/]*hapsolo[^/]*|[^/]*purge_dups[^/]*}/{genome_prefix}.input.{haplotype, hap[^_./@]+}/{genome_prefix}.input.{haplotype}.fasta"
    log:
        ln=config["out_dir"] / "dedup/{parameters}/log/create_dedup_input_links.{genome_prefix}.{parameters}.{haplotype}.ln.log",

