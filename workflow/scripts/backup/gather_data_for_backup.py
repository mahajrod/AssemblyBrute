#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse
from pathlib import Path


def backup_stage_files(stage_name, results_path, backup_path, file_pattern_list):
    print(f"Backuping data from the {stage_name} stage ...")
    stage_dir_path = results_path / stage_name    # Example: results/gap_closing/
    if stage_dir_path.exists():
        backup_stage_dir_path = backup_path / stage_name
        os.makedirs(backup_stage_dir_path, exist_ok=True)
        stat_file_path_list = list(stage_dir_path.glob("*.stage_stats"))
        for filename in stat_file_path_list:
            print(f"\tCopying {filename}/ ...")
            os.system(f"cp -rL {filename} {backup_stage_dir_path}")

        for stage_option_dir_path in stage_dir_path.glob("*"): # Example: results/gap_closing/draft_qc_def..samba_phased/
            if stage_option_dir_path.is_dir():
                print(f"\tCopying files for {stage_option_dir_path}/ ...")

                stage_option = stage_option_dir_path.name
                backup_stage_option_dir_path = backup_stage_dir_path / stage_option
                os.makedirs(backup_stage_option_dir_path, exist_ok=True)

                for pattern in file_pattern_list:
                    for filepath in stage_option_dir_path.glob(pattern):
                        print(f"\t\tCopying {filepath} ...")
                        os.system(f"cp -rL {filepath} {backup_stage_option_dir_path}")

                for common_dir in "assembly_qc", "wga":
                    common_dir_path = stage_option_dir_path / common_dir
                    if common_dir_path.exists():
                        backup_common_dir_stage_option_dir_path = backup_stage_option_dir_path / "assembly_qc/"
                        os.makedirs(backup_common_dir_stage_option_dir_path, exist_ok=True)
                        print(f"\t\tCopying {common_dir_path}/ ...")
                        os.system(f"cp -rL {common_dir_path} {backup_common_dir_stage_option_dir_path}")

                for hap_pattern in ".hap*", ".reordered", ".combined":
                    for hap_dir_path in stage_option_dir_path.glob(f"*{hap_pattern}"):
                        if hap_dir_path.is_dir():
                            hap_dir_name = hap_dir_path.name
                            backup_hap_dir_path = backup_stage_option_dir_path / hap_dir_name
                            os.makedirs(backup_hap_dir_path, exist_ok=True)
                            print(f"\t\tCopying {hap_dir_path}/ ...")
                            for analysis_dir_path in hap_dir_path.glob("*"):
                                analysis_dir_name = analysis_dir_path.name
                                print(f"\t\t\tCopying files for {analysis_dir_name}/ ...")
                                backup_analysis_dir_path = backup_hap_dir_path / analysis_dir_name
                                if analysis_dir_name == "alignment": # maybe do a selective copying in future
                                    os.system(f"cp -rL {analysis_dir_path} {backup_analysis_dir_path}")
                                else:
                                    os.system(f"cp -rL {analysis_dir_path} {backup_analysis_dir_path}")

                #for set_type in "reordered", "combined":
                #    set_type_path = stage_option_dir_path / set_type
                #    if set_type_path.exists:
                #        print(f"\t\tCopying files for {set_type} dataset...")
                #        backup_set_type_stage_option_dir_path = backup_stage_option_dir_path / set_type
                #        os.makedirs(backup_set_type_stage_option_dir_path, exist_ok=True)
                #        for pattern in "*.png", "*.svg", "per_chr", "*.pretext", ".rmdup.bam", ".rmdup.bam.csi", ".rmdup.bam.bai", ".assembly", ".agp", ".bed", ".syn":
                #            for filepath in (set_type_path / "alignment/NA/").glob(pattern):
                #                print(f"\t\t\tCopying {filepath}...")
                #                os.system(f"cp -rL {filepath} {backup_set_type_stage_option_dir_path}")
    else:
        print(f"\t{stage_name} dir was not found, skipping...")


parser = argparse.ArgumentParser()

parser.add_argument("-r", "--results_dir", action="store", dest="results_dir", required=True,
                    help="Directory containing results of the pipeline")


args = parser.parse_args()

results_dir_path = Path(args.results_dir)

backup_dir_path = results_dir_path / "backup/"
os.makedirs(backup_dir_path, exist_ok=True)

# backup QC
qc_dir_path = results_dir_path / "qc/"
print("Backuping QC data...")
if qc_dir_path.exists():
    print(f"\tCopying {qc_dir_path}...")

    os.system(f"cp -rL {qc_dir_path} {backup_dir_path}")
else:
    print("\tQC data not found, skipping...")

# backup kmer
kmer_dir_path = results_dir_path / "kmer/"
print("Backuping kmer data...")
if kmer_dir_path.exists():
    for folder_path in kmer_dir_path.glob("*"):
        kmer_datatype = folder_path.name
        filtered_kmer_dir_path = kmer_dir_path / kmer_datatype / "filtered/"

        backup_kmer_datatype_path = backup_dir_path / "kmer/" / kmer_datatype / "filtered/"
        os.makedirs(backup_kmer_datatype_path, exist_ok=True)
        print(filtered_kmer_dir_path)
        histo_file_path = list(filtered_kmer_dir_path.glob("*.histo"))
        print(histo_file_path)
        if histo_file_path:
            histo_file_path = histo_file_path[0]
            print(f"\tCopying {histo_file_path}...")
            os.system(f"cp -rL {histo_file_path} {backup_kmer_datatype_path}")

        genomescope_dir_path = list(filtered_kmer_dir_path.glob("genomescope"))
        print(genomescope_dir_path)
        if genomescope_dir_path:
            genomescope_dir_path = genomescope_dir_path[0]
            print(f"\tCopying {genomescope_dir_path}...")
            os.system(f"cp -rL {genomescope_dir_path} {backup_kmer_datatype_path}")
else:
    print("\tKmer data not found, skipping...")
"""
# backup contamination scan
print("Backuping contamination scan data...")
contamination_scan_dir_path = results_dir_path / "contamination_scan/kraken2/"
if contamination_scan_dir_path.exists():
    for folder_path in contamination_scan_dir_path.glob("*"):
        scan_datatype = folder_path.name
        kraken2_report_path = contamination_scan_dir_path  / scan_datatype / "kraken2.nt.report"

        backup_scan_datatype_path = backup_dir_path  / "contamination_scan/kraken2/" / scan_datatype
        os.makedirs(backup_scan_datatype_path, exist_ok=True)
        print(f"\tCopying {kraken2_report_path}...")
        os.system(f"cp -rL {kraken2_report_path} {backup_scan_datatype_path}")
else:
    print("\tContamination scan data not found, skipping...")
# backup error corrected hifi reads
print("Backuping corrected hifi reads data...")

error_correction_read_dir_path = results_dir_path / "error_correction/"
if error_correction_read_dir_path.exists():
    error_correction_read_backup_dir_path = backup_dir_path / "error_correction/"
    os.makedirs(error_correction_read_backup_dir_path, exist_ok=True)
    print(f"\tCopying {error_correction_read_dir_path }...")
    os.system(f"cp -rL {error_correction_read_dir_path } {error_correction_read_backup_dir_path}")

else:
    print("\tError correction read dir was not found, skipping...")

backup_stage_files("contig", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.gfa", "*.bed","*.bedgraph", "*.len", "*.cov", "*.lencov", "*.ids", "telomere"])
backup_stage_files("dedup", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.bed", "*.bedgraph", "*.len", "*.assembly", "*.agp", "*.ids", "telomere", "*.png", "*.svg"])
backup_stage_files("polishing", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.bed", "*.bedgraph", "*.len", "*.assembly", "*.agp", "*.ids", "telomere", "*.png", "*.svg"])
backup_stage_files("hic_alignment", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.bed", "*.bedgraph", "*.len", "*.assembly", "*.agp", "*.ids", "telomere", "*.png", "*.svg"])
backup_stage_files("hic_scaffolding", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.bed", "*.bedgraph", "*.len", "*.assembly", "*.agp", "*.ids", "*.hic", "telomere", "*.tab.gz", "*.png", "*.svg"])
for stage in "draft_qc", "gap_closing", "ref_scaffolding":
    backup_stage_files(stage, results_dir_path, backup_dir_path,
                    ["*.fasta", ".fai", "*.bed", "*.bedgraph", "*.len", "*.assembly", "*.agp", "*.ids", "telomere", "*.png", "*.svg", "*.tab.gz"])

"""