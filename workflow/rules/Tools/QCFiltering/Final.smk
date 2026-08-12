
localrules: create_final_fastq_links, create_final_fasta_links
localrules: create_final_track_fastq_links, create_final_track_fasta_links


use rule create_local_links as create_final_fastq_links with:
    input:
        fastq=lambda wildcards: config["out_dir"]/ ("data/%s/%s/%s%s" % (wildcards.datatype,
                                                                          "filtered" if wildcards.datatype in config["data_feature_dict"]["filter"] else "raw",
                                                                          wildcards.fileprefix,
                                                                          config["fastq_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fastq=config["out_dir"] / ("data/{datatype}/final/{fileprefix}%s" % config["fastq_ext"])
    log:
        ln=config["out_dir"] / "log/create_final_fastq_links.{datatype}.{fileprefix}.ln.log",

use rule create_local_links as create_final_fasta_links with:
    input:
        fastq=lambda wildcards: config["out_dir"]/ ("data/%s/%s/%s%s" % (wildcards.datatype,
                                                                         "filtered" if wildcards.datatype in config["data_feature_dict"]["filter"] else "raw",
                                                                         wildcards.fileprefix,
                                                                         config["fasta_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fastq=config["out_dir"] / ("data/{datatype}/final/{fileprefix}%s" % config["fasta_ext"])
    log:
        ln=config["out_dir"] / "log/create_final_fasta_links.{datatype}.{fileprefix}.ln.log",

use rule create_local_links as create_final_track_fastq_links with:
    input:
        fastq=lambda wildcards: config["out_dir"]/ ("ext_data/%s/%s/%s/%s%s" % (wildcards.datatype, wildcards.track_name,
                                                                                  "filtered" if wildcards.datatype in config["data_feature_dict"]["filter"] else "raw",
                                                                                  wildcards.fileprefix,
                                                                                  config["fastq_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fastq=config["out_dir"] / ("ext_data/{datatype}/{track_name}/final/{fileprefix}%s" % config["fastq_ext"])
    log:
        ln=config["out_dir"] / "log/create_final_track_fastq_links.{datatype}.{track_name}.{fileprefix}.ln.log",

use rule create_local_links as create_final_track_fasta_links with:
    input:
        fastq=lambda wildcards: config["out_dir"]/ ("ext_data/%s/%s/%s/%s%s" % (wildcards.datatype, wildcards.track_name,
                                                                                  "filtered" if wildcards.datatype in config["data_feature_dict"]["filter"] else "raw",
                                                                                  wildcards.fileprefix,
                                                                                  config["fasta_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fastq=config["out_dir"] / ("ext_data/{datatype}/{track_name}/final/{fileprefix}%s" % config["fasta_ext"])
    log:
        ln=config["out_dir"] / "log/create_final_track_fasta_links.{datatype}.{track_name}.{fileprefix}.ln.log",
