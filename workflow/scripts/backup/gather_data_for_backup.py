#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse
from pathlib import Path


def backup_stage_files(stage_name, results_path, backup_path, file_pattern_list):
    print(f"Backuping data from the {stage_name} stage...")
    stage_dir_path = results_path / stage_name
    if stage_dir_path.exists():
        backup_stage_dir_path = backup_path / stage_name
        os.makedirs(backup_stage_dir_path, exist_ok=True)
        stat_file_path_list = list(stage_dir_path.glob("*.stage_stats"))
        for filename in stat_file_path_list:
            print(f"\tCopying {filename}...")
            os.system(f"cp -r {filename} {backup_stage_dir_path}")

        for stage_option_dir_path in stage_dir_path.glob("*"):
            if stage_option_dir_path.is_dir():
                print(f"\tCopying files for {stage_option_dir_path}...")

                stage_option = stage_option_dir_path.name
                backup_stage_option_dir_path = backup_stage_dir_path / stage_option
                os.makedirs(backup_stage_option_dir_path, exist_ok=True)

                for pattern in file_pattern_list:
                    for filepath in stage_option_dir_path.glob(pattern):
                        print(f"\t\tCopying {filepath}...")
                        os.system(f"cp -r {filepath} {backup_stage_option_dir_path}")

                assembly_qc_dir_path = stage_option_dir_path / "assembly_qc/"
                if assembly_qc_dir_path.exists:
                    backup_assembly_qc_stage_option_dir_path = backup_stage_option_dir_path / "assembly_qc/"
                    os.makedirs(backup_assembly_qc_stage_option_dir_path, exist_ok=True)

                    for filepath in assembly_qc_dir_path.glob("*") :
                        #print(filepath.name)
                        if filepath.name == "merqury":
                            print("\t\tCopying merqury files...")
                            backup_merqury_qc_path = backup_assembly_qc_stage_option_dir_path / "merqury/"
                            os.makedirs(backup_merqury_qc_path, exist_ok=True)
                            for merqury_filepath in filepath.glob("*"):
                                if merqury_filepath.name[-6:] == ".fasta":
                                    continue
                                if merqury_filepath.name[-6:] == ".meryl":
                                    continue
                                print(f"\t\t\tCopying {merqury_filepath}...")
                                os.system(f"cp -r {merqury_filepath} {backup_merqury_qc_path}")

                        else:
                            print(f"\t\tCopying {filepath}...")
                            os.system(f"cp -r {filepath} {backup_assembly_qc_stage_option_dir_path}")
                for set_type in "reordered", "combined":
                    set_type_path = stage_option_dir_path / set_type
                    if set_type_path.exists:
                        print(f"\t\tCopying files for {set_type} dataset...")
                        print(set_type_path)
                        backup_set_type_stage_option_dir_path = backup_stage_option_dir_path / set_type
                        os.makedirs(backup_set_type_stage_option_dir_path, exist_ok=True)
                        for pattern in "*.png", "*.svg", "per_chr", "*.pretext", ".rmdup.bam", ".rmdup.bam.csi", ".rmdup.bam.bai", ".assembly", ".agp", ".bed", ".syn":
                            for filepath in (set_type_path / "alignment/NA/").glob(pattern):
                                print(f"\t\t\tCopying {filepath}...")
                                os.system(f"cp -r {filepath} {backup_set_type_stage_option_dir_path}")
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

    os.system(f"cp -r {qc_dir_path} {backup_dir_path}")
else:
    print("\tQC data not found, skipping...")

# backup kmer
kmer_dir_path = results_dir_path / "kmer/"
print("Backuping kmer data...")
if kmer_dir_path.exists():
    for folder_path in kmer_dir_path.glob("*"):
        kmer_datatype = folder_path.name
        filtered_kmer_dir_path = kmer_dir_path / kmer_datatype / "filtered/"
        histo_file_path = list(filtered_kmer_dir_path.glob("*.histo"))[0]
        genomescope_dir_path = list(filtered_kmer_dir_path.glob("genomescope"))[0]

        backup_kmer_datatype_path = backup_dir_path / "kmer/" / kmer_datatype / "filtered/"
        os.makedirs(backup_kmer_datatype_path, exist_ok=True)

        print(f"\tCopying {histo_file_path}...")
        os.system(f"cp -r {histo_file_path} {backup_kmer_datatype_path}")
        #shutil.copy(histo_file_path, backup_kmer_datatype_path)
        print(f"\tCopying {genomescope_dir_path}...")
        os.system(f"cp -r {genomescope_dir_path} {backup_kmer_datatype_path}")
else:
    print("\tKmer data not found, skipping...")

# backup contamination scan
print("Backuping contamination scan data...")
contamination_scan_dir_path = results_dir_path / "contamination_scan/kraken2/"
if contamination_scan_dir_path.exists():
    for folder_path in contamination_scan_dir_path.glob("*"):
        scan_datatype = folder_path.name
        print(scan_datatype)
        kraken2_report_path = contamination_scan_dir_path  / scan_datatype / "kraken2.nt.report"

        backup_scan_datatype_path = backup_dir_path  / "contamination_scan/kraken2/" / scan_datatype
        os.makedirs(backup_scan_datatype_path, exist_ok=True)
        print(f"\tCopying {kraken2_report_path}...")
        os.system(f"cp -r {kraken2_report_path} {backup_scan_datatype_path}")
else:
    print("\tContamination scan data not found, skipping...")
# backup error corrected hifi reads
print("Backuping corrected hifi reads data...")

error_correction_read_dir_path = results_dir_path / "error_correction/"
if error_correction_read_dir_path.exists():
    hifi_read_dir_path = backup_dir_path / "data/hifi/"
    os.makedirs(hifi_read_dir_path, exist_ok=True)
    print(f"\tCopying {error_correction_read_dir_path }...")
    os.system(f"cp -r {error_correction_read_dir_path } {hifi_read_dir_path}")

else:
    print("\tError correction read dir was not found, skipping...")

backup_stage_files("contig", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.gfa", "*.bed", "*.len", "*.cov", "*.lencov", "*.ids", "telomere"])
backup_stage_files("purge_dups", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.bed", "*.len", "*.assembly", "*.agp", "*.ids", "telomere",
                              "*.png", "*.svg"])
backup_stage_files("hic_scaffolding", results_dir_path, backup_dir_path,
                ["*.fasta", ".fai", "*.bed", "*.len", "*.assembly", "*.agp", "*.ids", "*.hic", "telomere",
                              "*.tab.gz", "*.png", "*.svg"])
for stage in "draft_qc", "gap_closing", "ref_scaffolding":
    backup_stage_files(stage, results_dir_path, backup_dir_path,
                    ["*.fasta", ".fai", "*.bed", "*.len", "*.assembly", "*.agp", "*.ids", "telomere",
                                  "*.png", "*.svg", "*.tab.gz"])

