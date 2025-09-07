#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import shutil
import argparse

from pathlib import Path
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


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
print(f"\tCopying {qc_dir_path}...")

os.system(f"cp -r {qc_dir_path} {backup_dir_path}")
#shutil.copy(qc_dir_path, backup_dir)

# backup kmer
kmer_dir_path = results_dir_path / "kmer/"
for folder_path in kmer_dir_path.glob("*"):
    kmer_datatype = folder_path.name
    filtered_kmer_dir_path = kmer_dir_path / kmer_datatype / "filtered/"
    histo_file_path = list(filtered_kmer_dir_path.glob("*.histo"))[0]
    genomescope_dir_path = list(filtered_kmer_dir_path.glob("genomescope"))[0]

    backup_kmer_datatype_path = backup_dir_path / "kmer/" / kmer_datatype / "filtered/"
    os.makedirs(backup_kmer_datatype_path, exist_ok=True)

    print("Backuping kmer data...")
    print(f"\tCopying {histo_file_path}...")
    shutil.copy(histo_file_path, backup_kmer_datatype_path)
    print(f"\tCopying {genomescope_dir_path}...")
    os.system(f"cp -r {genomescope_dir_path} {backup_kmer_datatype_path}")

