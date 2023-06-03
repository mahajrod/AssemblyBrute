#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path


parser = argparse.ArgumentParser()


parser.add_argument("-s", "--stat_file", action="store", dest="stat_file", required=True,
                    help="File wit statistics extracted from coverage file.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

artefact_set = {"HAPLOTIG", "JUNK", "REPEAT", "OVLP", "HIGHCOV"}
stat_cov_df = pd.read_csv(args.stat_file, sep="\t", header=None, names=["#scaffold", "length", "mean_cov", "median_cov"],
                          index_col=0)

