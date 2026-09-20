

localrules: create_se_fastq_links, create_pe_fastq_links
localrules: create_se_fasta_links, create_pe_fasta_links
localrules: create_links_for_reference, create_link_for_draft

localrules: create_track_se_fastq_links, create_track_pe_fastq_links
localrules: create_track_se_fasta_links, create_track_pe_fasta_links

use rule create_local_links as create_se_fastq_links with:
    input:
        fastq=lambda wildcards: config["data"][wildcards.se_datatype]["in_dir"] / ("{fileprefix}%s" % config["data"][wildcards.se_datatype]["in_ext"]),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fastq=config["out_dir"] / ("data/{se_datatype}/raw/{fileprefix}%s" % config["fastq_ext"])
    log:
        ln=config["out_dir"] / "log/create_se_fastq_links.{se_datatype}.{fileprefix}.ln.log",

use rule create_local_links as create_track_se_fastq_links with:
    input:
        fastq=lambda wildcards: config["ext_data"][wildcards.se_datatype][wildcards.track_name]["in_dir"] / ("{fileprefix}%s" % config["ext_data"][wildcards.se_datatype][wildcards.track_name]["in_ext"]),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fastq=config["out_dir"] / ("ext_data/{se_datatype}/{track_name}/raw/{fileprefix}%s" % config["fastq_ext"])
    log:
        ln=config["out_dir"] / "log/create_track_se_fastq_links.{se_datatype}.{track_name}.{fileprefix}.ln.log",

use rule create_local_links as create_pe_fastq_links with: # abstract rule to use via redefining
    input:
        forward_fastq=lambda wildcards: config["data"][wildcards.pe_datatype]["in_dir"] / ("{pairprefix}%s%s" % (config["data"][wildcards.pe_datatype]["in_fwd_sfx"],
                                                                                                                 config["data"][wildcards.pe_datatype]["in_ext"])),
        reverse_fastq=lambda wildcards: config["data"][wildcards.pe_datatype]["in_dir"] / ("{pairprefix}%s%s" % (config["data"][wildcards.pe_datatype]["in_rev_sfx"],
                                                                                                                 config["data"][wildcards.pe_datatype]["in_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        forward_fastq=config["out_dir"] / ("data/{pe_datatype}/raw/{pairprefix}%s%s" % (config["fwd_fastq_sfx"],
                                                                                        config["fastq_ext"])),
        reverse_fastq=config["out_dir"] / ("data/{pe_datatype}/raw/{pairprefix}%s%s" % (config["rev_fastq_sfx"],
                                                                                        config["fastq_ext"])),
    log:
        ln=config["out_dir"] / "log/create_pe_fastq_links.{pe_datatype}.{pairprefix}.ln.log",

use rule create_local_links as create_track_pe_fastq_links with: # abstract rule to use via redefining
    input:
        forward_fastq=lambda wildcards: config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_dir"] / ("{pairprefix}%s%s" % (config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_fwd_sfx"],
                                                                                                                                             config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_ext"])),
        reverse_fastq=lambda wildcards: config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_dir"] / ("{pairprefix}%s%s" % (config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_rev_sfx"],
                                                                                                                                             config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        forward_fastq=config["out_dir"] / ("ext_data/{pe_datatype}/{track_name}/raw/{pairprefix}%s%s" % (config["fwd_fastq_sfx"],
                                                                                        config["fastq_ext"])),
        reverse_fastq=config["out_dir"] / ("ext_data/{pe_datatype}/{track_name}/raw/{pairprefix}%s%s" % (config["rev_fastq_sfx"],
                                                                                        config["fastq_ext"])),
    log:
        ln=config["out_dir"] / "log/create_track_pe_fastq_links.{pe_datatype}.{track_name}.{pairprefix}.ln.log",

use rule create_local_links as create_se_fasta_links with:
    input:
        fasta=lambda wildcards: config["data"][wildcards.se_datatype]["in_dir"] / ("{fileprefix}%s" % config["data"][wildcards.se_datatype]["in_ext"]),
        #fasta=input_dir_path / ("{se_datatype}/fasta/{fileprefix}%s" %  config["fasta_ext"]),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fasta=config["out_dir"] / ("data/{se_datatype}/raw/{fileprefix}%s" % config["fasta_ext"])
    log:
        ln=config["out_dir"] / "log/create_se_fasta_links.{se_datatype}.{fileprefix}.log",

use rule create_local_links as create_track_se_fasta_links with:
    input:
        fasta=lambda wildcards: config["ext_data"][wildcards.se_datatype][wildcards.track_name]["in_dir"] / ("{fileprefix}%s" % config["ext_data"][wildcards.se_datatype][wildcards.track_name]["in_ext"]),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fasta=config["out_dir"] / ("ext_data/{se_datatype}/{track_name}/raw/{fileprefix}%s" % config["fasta_ext"])
    log:
        ln=config["out_dir"] / "log/create_track_se_fasta_links.{track_name}.{se_datatype}.{fileprefix}.log",

use rule create_local_links as create_pe_fasta_links with:
    input:
        forward_fastq=lambda wildcards: config["data"][wildcards.pe_datatype]["in_dir"] / ("{pairprefix}%s%s" % (config["data"][wildcards.pe_datatype]["in_fwd_sfx"],
                                                                                                          config["data"][wildcards.pe_datatype]["in_ext"])),
        reverse_fastq=lambda wildcards: config["data"][wildcards.pe_datatype]["in_dir"] / ("{pairprefix}%s%s" % (config["data"][wildcards.pe_datatype]["in_rev_sfx"],
                                                                                                                  config["data"][wildcards.pe_datatype]["in_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        forward_fastq=config["out_dir"] / ("data/{pe_datatype}/raw/{pairprefix}%s%s" % (config["fwd_fasta_sfx"],
                                                                                        config["fasta_ext"])),
        reverse_fastq=config["out_dir"] / ("data/{pe_datatype}/raw/{pairprefix}%s%s" % (config["rev_fasta_sfx"],
                                                                                        config["fasta_ext"])),
    log:
        ln=config["out_dir"] / "log/create_pe_fasta_links.{pe_datatype}.{pairprefix}.ln.log",

use rule create_local_links as create_track_pe_fasta_links with:
    input:
        forward_fastq=lambda wildcards: config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_dir"] / ("{pairprefix}%s%s" % (config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_fwd_sfx"],
                                                                                                                                             config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_ext"])),
        reverse_fastq=lambda wildcards: config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_dir"] / ("{pairprefix}%s%s" % (config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_rev_sfx"],
                                                                                                                                             config["ext_data"][wildcards.pe_datatype][wildcards.track_name]["in_ext"])),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        forward_fastq=config["out_dir"] / ("ext_data/{pe_datatype}/{track_name}/raw/{pairprefix}%s%s" % (config["fwd_fasta_sfx"],
                                                                                                           config["fasta_ext"])),
        reverse_fastq=config["out_dir"] / ("ext_data/{pe_datatype}/{track_name}raw/{pairprefix}%s%s" % (config["rev_fasta_sfx"],
                                                                                                          config["fasta_ext"])),
    log:
        ln=config["out_dir"] / "log/create_track_pe_fasta_links.{pe_datatype}.{track_name}.{pairprefix}.ln.log",


use rule create_local_links as create_link_for_draft with:
    input:
        draft=lambda wildcards: input_dir_path / "draft/fasta/{0}".format(config["data"]["draft"]["haplotypes"][wildcards.haplotype]),
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        draft=config["out_dir"] / "draft_qc/{parameters}/{genome_prefix}.draft_qc.{haplotype, hap[^/]*}.fasta"
    log:
        ln=config["out_dir"] / "log/create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.ln.log",

use rule create_local_links as create_links_for_reference with: # abstract rule to use via redefining
    input:
        fasta=lambda wildcards: config["data"]["reference"]["ref_dict"][wildcards.ref_name]["fasta"],
        syn=lambda wildcards: config["data"]["reference"]["ref_dict"][wildcards.ref_name]["syn"],
        whitelist=lambda wildcards: config["data"]["reference"]["ref_dict"][wildcards.ref_name]["whitelist"],
        orderlist=lambda wildcards: config["data"]["reference"]["ref_dict"][wildcards.ref_name]["orderlist"],
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        fasta=config["out_dir"] / "data/reference/{ref_name}/{ref_name}.fasta",
        syn=config["out_dir"] / "data/reference/{ref_name}/{ref_name}.syn",
        whitelist=config["out_dir"] / "data/reference/{ref_name}/{ref_name}.custom.whitelist",
        orderlist=config["out_dir"] / "data/reference/{ref_name}/{ref_name}.custom.orderlist",
    log:
        ln=config["out_dir"] / "log/create_links_for_reference.{ref_name}.ln.log",

