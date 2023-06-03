#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--pb_stat_file_list", action="store", dest="pb_stat_file_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of per base stats (usually PB.stat) "
                         "produced by purge dups for all haplotypes.")
parser.add_argument("-c", "--cutoff_file_list", action="store", dest="cutoff_file_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of cutoffs (usually cutoffs) produced by purge dups for all haplotypes.")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list labels for all haplotypes.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", default=["svg", "png"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of graphical formats supported by matplotlib. Default: 'png,svg'")

parser.add_argument("--x_label", action="store", dest="x_label", default="coverage",
                    help="Label for x axis. Default: 'coverage'")
parser.add_argument("--y_label", action="store", dest="y_label", default="counts",
                    help="Label for y axis. Default: 'counts")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")
parser.add_argument("-s", "--suptitle", action="store", dest="suptitle", default=None,
                    help="Suptitle of figure. Default: not set")

parser.add_argument("--xmin", action="store", dest="xmin", type=float, default=1,
                    help="Minimum limit for Y axis . Default: 1")
parser.add_argument("--xmax",  action="store", dest="xmax", type=float, default=None,
                    help="Maximum limit for X axis. Default: 3.3 * diploid coverage")

parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.1,
                    help="Adjust left border of subplots on the figure. Default: 0.05")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float, default=0.9,
                    help="Adjust top border of subplots on the figure. Default: 0.05")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float, default=0.9,
                    help="Adjust right border of subplots on the figure. Default: 0.95")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float, default=0.1,
                    help="Adjust bottom border of subplots on the figure. Default: 0.95")
parser.add_argument("--subplots_adjust_hspace", action="store", dest="subplots_adjust_hspace", type=float, default=0.1,
                    help="Adjust height space between subplots on the figure. Default: 0.1")
parser.add_argument("--subplots_adjust_wspace", action="store", dest="subplots_adjust_wspace", type=float, default=0.1,
                    help="Adjust width space between  subplots on the figure. Default: 0.1")

args = parser.parse_args()

number_of_haplotypes = len(args.label_list)

haplotype_color_list = ["blue", "green", "orange", "purple", "black"] + ["black"] * max(0, number_of_haplotypes - 5)
haplotype_style_list = ["-", "--", "-.", ":", ","] + [","] * max(0, number_of_haplotypes - 5)


pb_stat_df_dict = {haplotype: pd.read_csv(filename, sep="\t", header=None,
                                          names=["coverage", "counts"]) for filename, haplotype in zip(args.pb_stat_file_list, args.label_list)}

cutoff_df_dict = {haplotype: pd.read_csv(filename, sep="\t", header=None,
                                         names=["junk_cov", "param1", "inter_cov", "param2",
                                                "param3", "high_cov"]) for filename, haplotype in zip(args.cutoff_file_list, args.label_list)}

# junk cov
# param1 = ?
# inter_cov = coverage at the gap between haploid and diploid peaks
# param2 = inter_cov + 1
# param3 = diploid_coverage + (diploid_coverage - param2)
# high_cov = 3 * diploid_coverage

fig, ax = plt.subplots(1, figsize=(8, 8), dpi=args.dpi)

max_high_cov = 0

for haplotype, cuttoff_df, color, linestyle in zip(pb_stat_df_dict.keys(), cutoff_df_dict,
                                                   haplotype_color_list, haplotype_style_list):
    ax.plot(pb_stat_df_dict[haplotype]["coverage"], pb_stat_df_dict[haplotype]["counts"],
            color=color, label=haplotype, linestyle=linestyle)
    max_high_cov = max(max_high_cov, cutoff_df_dict[haplotype]["high_cov"][0])
    for line_name, line_color in zip(["junk_cov", "param2", "high_cov"], ["red", "magenta", "cyan"]):
        ax.axvline(cutoff_df_dict[haplotype][line_name][0], color=line_color, zorder=1, linestyle=linestyle)

ax.set_ylabel(args.y_label)
ax.set_xlabel(args.x_label)
ax.legend(loc='upper right')
ax.grid(visible=True, linestyle="--")
ax.set_xlim(xmin=args.xmin, xmax=args.xmax if args.xmax is not None else 1.1 * max_high_cov)
plt.suptitle(args.suptitle)
plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right,
                    top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom,
                    hspace=args.subplots_adjust_hspace, wspace=args.subplots_adjust_wspace)

for ext in args.extension_list:
    fig.savefig("{0}.{1}".format(args.output_prefix, ext))
