#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from pathlib import Path

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
parser.add_argument("-a", "--haplotype_list", action="store", dest="haplotype_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of qc'd haplotypes. "
                         "If absent, busco report will be ignored.")
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

#if args.busco_database_list:
#    busco_path_list = list(qc_folder_path.glob("busco5/*full_table.tsv"))

df_dict = {}  # {haplotype: {} for haplotype in args.haplotype_list}

quast_columns = ["# contigs (>= 0 bp)", "# contigs (>= 10000 bp)",
                 "# contigs (>= 10000 bp)", "Total length (>= 10000 bp)",
                 "Largest contig", "Largest contig", "GC (%)", "N50", "L50"]

for haplotype in args.haplotype_list:
    df_dict[haplotype] = {}
    df_dict[haplotype]["quast"] = pd.read_csv(qc_folder_path / "quast/{0}.{1}/report.tsv".format(args.input_prefix, haplotype), sep="\t",
                                              header=0, index_col=0).transpose()
    df_dict[haplotype]["busco5"] = {}
    for busco_db in args.busco_database_list:
        df_dict[haplotype]["busco5"][busco_db] = pd.DataFrame([read_gene_string_from_busco_summary(qc_folder_path / "busco5/{0}.{1}.busco5.{2}.summary".format(args.input_prefix, haplotype, busco_db))],
                                                              columns=[busco_db, ], index=pd.Index([haplotype, ]))
        #BUSCOtable(in_file=qc_folder_path / "busco5/{0}.{1}.busco5.{2}.full_table.tsv".format(args.input_prefix,
        #           haplotype,
        #           busco_db))
merqury_qv_df = pd.read_csv(qc_folder_path / "merqury/{0}.qv".format(args.input_prefix),
                            sep="\t", index_col=0, header=None,
                            names=["haplotype", "unique_kmers", "read_and_assembly_kmers", "qv", "error_rate"])
merqury_qv_df.rename(index={"{0}.{1}".format(args.input_prefix,
                                             haplotype): haplotype for haplotype in args.haplotype_list},
                     inplace=True)

merqury_completeness_df = pd.read_csv(qc_folder_path / "merqury/{0}.completeness.stats".format(args.input_prefix),
                                      sep="\t", index_col=0, header=None,
                                      names=["haplotype", "kmer_set", "assembly_solid_kmers",
                                             "read_solid_kmers", "completeness"])
merqury_completeness_df.rename(index={"{0}.{1}".format(args.input_prefix,
                                                       haplotype): haplotype for haplotype in args.haplotype_list},
                               inplace=True)


final_df = pd.DataFrame([[stage, parameters] for stage, parameters in zip([args.stage] * len(args.haplotype_list),
                                                                          [args.parameters] * len(args.haplotype_list))],
                        index=pd.Index([haplotype for haplotype in args.haplotype_list], name="haplotype"),
                        columns=["stage", "parameters"])

final_df = pd.concat([final_df,
                      pd.concat([df_dict[haplotype]["quast"][quast_columns] for haplotype in args.haplotype_list]),
                      merqury_qv_df.loc[args.haplotype_list],
                      merqury_completeness_df[["assembly_solid_kmers", "read_solid_kmers", "completeness"]].loc[args.haplotype_list] if not merqury_completeness_df.empty else pd.DataFrame(),
                      *[pd.concat([df_dict[haplotype]["busco5"][busco_db] for haplotype in args.haplotype_list]) for busco_db in args.busco_database_list]
                      ],
                     axis=1)

for column in ["# contigs (>= 0 bp)", "# contigs (>= 10000 bp)",
               "# contigs (>= 10000 bp)", "Total length (>= 10000 bp)",
               "Largest contig", "Largest contig", "L50"]:
    final_df[column] = final_df[column].astype("Int64")

final_df.index.name = "haplotype"
final_df.to_csv(args.output, sep="\t", header=True, index=True)
