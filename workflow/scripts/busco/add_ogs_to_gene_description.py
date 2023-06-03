#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("-b", "--busco_dataset_info", action="store", dest="busco_dataset_info", required=True,
                    help="BUSCO dataset info file. "
                         "Usually it is located at busco_downloads/lineages/mammalia_odb10/info/ogs.id.info")
parser.add_argument("-g", "--gene_description", action="store", dest="gene_description", required=True,
                    help="OrthoDB gene description file."
                         "For OrthoDB it is odb10v0_genes.tab.gz")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file tab file")

args = parser.parse_args()

busco_dataset_df = pd.read_csv(args.busco_dataset_info, sep="\t", names=["busco_gene_id", "OG"], index_col=0)

genes_df = pd.read_csv(args.gene_description, sep="\t", names=["busco_gene_id", "busco_species_id", "ncbi_id",
                                                               "common_name", "description"], index_col=0)

busco_dataset_df.merge(genes_df, how="left", left_on="busco_gene_id",
                       right_on="busco_gene_id").to_csv(args.output, sep="\t", index=True, header=False)

