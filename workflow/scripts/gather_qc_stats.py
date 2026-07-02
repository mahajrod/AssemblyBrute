#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse

from pathlib import Path
from copy import deepcopy

import pandas as pd

from RouToolPa.Parsers.BUSCO import BUSCOtable


def read_gene_string_from_busco_summary(in_file):
    with open(in_file, "r") as in_fd:
        for line in in_fd:
            if "***** Results: *****" in line:
                break
        in_fd.readline()
        return in_fd.readline().strip()


parser = argparse.ArgumentParser()

parser.add_argument("-q", "--qc_folder", action="store", dest="qc_folder", required=True,
                    help="Folder with qc data")
parser.add_argument("-e", "--input_prefix", action="store", dest="input_prefix", required=True,
                    help="Prefix of input files.")
parser.add_argument("-b", "--busco_database_list", action="store", dest="busco_database_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of busco databases used for QC. "
                         "If absent, busco report will be ignored.")
parser.add_argument("-a", "--assembly_prefix_list", action="store", dest="assembly_prefix_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of prefixes of qc'd assemblies. "
                         "If absent, busco report will be ignored.")
parser.add_argument("-m", "--merqury_datatype_list", action="store", dest="merqury_datatype_list", required=False,
                    type=lambda s: s.split(","), default=[],
                    help="Comma-separated list of datatypes used for merqury qc. "
                         "If absent, merqury reports will be ignored.")
parser.add_argument("-s", "--stage", action="store", dest="stage", required=True,
                    help="Assembly stage")
parser.add_argument("-p", "--parameters", action="store", dest="parameters", required=True,
                    help="Assembly parameters")
parser.add_argument("-o", "--output", action="store", dest="output",
                    default=sys.stdout, help="output file - default: stdout")

args = parser.parse_args()

qc_folder_path = Path(args.qc_folder)

if not qc_folder_path.exists():
    raise ValueError("ERROR!!! QC folder '{}' doesn't exist".format(args.qc_folder))

df_dict = {}

quast_columns = ["# contigs (>= 0 bp)", "# contigs (>= 10000 bp)",
                 "# contigs (>= 10000 bp)", "Total length (>= 10000 bp)",
                 "Largest contig", "Largest contig", "GC (%)", "N50", "L50"]

for assembly_prefix in args.assembly_prefix_list:
    df_dict[assembly_prefix] = {}
    df_dict[assembly_prefix]["quast"] = pd.read_csv(qc_folder_path / f"quast/{assembly_prefix}/report.tsv", sep="\t",
                                              header=0, index_col=0).transpose()
    df_dict[assembly_prefix]["busco5"] = {}
    for busco_db in args.busco_database_list:
        df_dict[assembly_prefix]["busco5"][busco_db] = pd.DataFrame([read_gene_string_from_busco_summary(qc_folder_path / "busco5/{0}.{1}.busco5.summary".format(assembly_prefix, busco_db))],
                                                              columns=[busco_db, ], index=pd.Index([assembly_prefix, ]))

merqury_df_list = []
for merqury_datatype in args.merqury_datatype_list:
    merqury_qv_path = qc_folder_path / "merqury/{1}/{0}.{1}.qv".format(args.input_prefix, merqury_datatype)
    if merqury_qv_path.exists():
        merqury_qv_df = pd.read_csv(merqury_qv_path,
                                    sep="\t", index_col=0, header=None,
                                    names=["assembly_prefix", "unique_kmers", "read_and_assembly_kmers", "qv", "error_rate"]).iloc[0: len(args.assembly_prefix_list)]

    else:
        sys.stdout.write(f"WARNING!!! Merqury QV file for datatype {merqury_datatype} is absent! Skipping...\n")
        merqury_qv_df = None

    merqury_completeness_stats_path = qc_folder_path / "merqury/{1}/{0}.{1}.completeness.stats".format(args.input_prefix, merqury_datatype)
    if merqury_completeness_stats_path.exists():
        merqury_completeness_df = pd.read_csv(merqury_completeness_stats_path,
                                              sep="\t", index_col=0, header=None,
                                              names=["assembly_prefix", "kmer_set", "assembly_solid_kmers",
                                                     "read_solid_kmers", "completeness"]).iloc[0: len(args.assembly_prefix_list)]
    else:
        sys.stdout.write(f"WARNING!!! Merqury completeness file for datatype {merqury_datatype} is absent! Skipping...\n")
        merqury_completeness_df = None

    if (merqury_qv_df is None) and (merqury_completeness_df is None):
        continue
    else:
        if (merqury_qv_df is not None) and (merqury_completeness_df is not None):
            merqury_concatenated_df = pd.concat([merqury_qv_df, merqury_completeness_df], axis=1)
        elif merqury_qv_df is not None:
            merqury_concatenated_df = merqury_qv_df
        else:
            merqury_concatenated_df = merqury_completeness_df
        merqury_concatenated_df.rename(columns={column: merqury_datatype + "@" + column for column in list(merqury_concatenated_df.columns)},
                                       inplace=True)
        merqury_df_list.append(deepcopy(merqury_concatenated_df))

final_df = pd.DataFrame([[stage, parameters] for stage, parameters in zip([args.stage] * len(args.assembly_prefix_list),
                                                                          [args.parameters] * len(args.assembly_prefix_list))],
                        index=pd.Index([assembly_prefix for assembly_prefix in args.assembly_prefix_list], name="assembly_prefix"),
                        columns=["stage", "parameters"])

final_df = pd.concat([final_df,
                      pd.concat([df_dict[assembly_prefix]["quast"][quast_columns] for assembly_prefix in args.assembly_prefix_list]),
                      *merqury_df_list,
                      *[pd.concat([df_dict[assembly_prefix]["busco5"][busco_db] for assembly_prefix in args.assembly_prefix_list]) for busco_db in args.busco_database_list]
                      ],
                     axis=1)

for column in ["# contigs (>= 0 bp)", "# contigs (>= 10000 bp)",
               "# contigs (>= 10000 bp)", "Total length (>= 10000 bp)",
               "Largest contig", "Largest contig", "L50"]:
    final_df[column] = final_df[column].astype("Int64")

final_df.index.name = "assembly_prefix"
final_df.to_csv(args.output, sep="\t", header=True, index=True)
