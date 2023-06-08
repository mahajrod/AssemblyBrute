#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
from copy import deepcopy

import pandas as pd

from RouToolPa.Parsers.BUSCO import BUSCOtable


def create_track_from_series(series,):
    prev_element = series.iloc[0]
    count_list = []
    counter = 1
    for element in series.iloc[1:]:
        if element == prev_element:
            counter += 1
        else:
            count_list.append([prev_element, counter])
            prev_element = element
            counter = 1
    count_list.append([element, counter])
    count_df = pd.DataFrame.from_records(count_list, columns=("status", "counts"))
    count_df["end"] = count_df["counts"].cumsum()
    count_df["start"] = count_df["end"].shift(1, fill_value=0)
    count_df["scaffold"] = series.name
    return count_df[["scaffold", "start", "end", "status"]]


def create_track_from_df(dataframe, color_dict, output_prefix=None):

    track_dict = {label: create_track_from_series(dataframe[label]) for label in dataframe.columns}
    track_merged_df = pd.concat(track_dict.values())
    track_merged_df["status"] = track_merged_df["status"].replace(color_dict)
    if output_prefix:
        track_merged_df.to_csv("{0}.busco.counts.bedgraph".format(output_prefix), sep="\t", header=False, index=False)
    return track_merged_df


def get_informative_buscos(dataframe):
    return dataframe[dataframe.apply(lambda s: True if s.nunique() > 1 else False, axis=1)]


def prepare_len_df(dataframe, output_prefix=None):
    n_buscos = len(dataframe)
    len_df = pd.DataFrame.from_dict({label: n_buscos for label in dataframe.columns}, orient='index', columns=["length"])
    len_df.index.name = "scaffold"
    len_df.to_csv("{0}.busco.len".format(output_prefix), sep="\t", header=False, index=True)
    return len_df


parser = argparse.ArgumentParser()

parser.add_argument("-b", "--busco_table_list", action="store", dest="busco_table_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with tables generated by BUSCO. One file per assembly. "
                         "Required.")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of labels for assemblies used in BUSCO runs. "
                         "Must have exactly the same length as --busco_table_list. Required.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files. Required.")

args = parser.parse_args()


busco_categories = ['Complete', 'Duplicated', 'Fragmented', 'Missing', ]
busco_color_dict = {'Complete':     "green",
                    'Duplicated':   "blue",
                    'Fragmented':   "orange",
                    'Missing':      "red"}

hap_file_dict = {label: filename for label, filename in zip(args.label_list, args.busco_table_list)}

busco_table_dict = {label: BUSCOtable(in_file=hap_file_dict[label]) for label in hap_file_dict}

number_of_buscos = len(busco_table_dict[args.label_list[0]].records[["OG"]].drop_duplicates().reset_index(level=1, drop=True))
print(busco_table_dict)
pd.Series(args.label_list).to_csv("{0}.busco.orderlist".format(args.output_prefix), sep="\t", index=False, header=False)


busco_legend_df = pd.DataFrame.from_dict(busco_color_dict, columns=["color"], orient="index")
busco_legend_df.reset_index(drop=False).set_index("color").to_csv("{0}.busco.legend".format(args.output_prefix),
                                                                  sep="\t", header=False, index=True)

busco_status_df_dict = {label: busco_table_dict[label].records[["OG"]].drop_duplicates().reset_index(level=1, drop=True).reset_index(drop=False).rename(columns={"status": label}).set_index("OG") for label in busco_table_dict}

busco_status_df = pd.concat(busco_status_df_dict.values(), axis=1).sort_values(by=args.label_list)
busco_status_df.to_csv("{0}.busco.merged.tsv".format(args.output_prefix), sep="\t", header=True, index=True)

count_track_merged_df = create_track_from_df(busco_status_df, busco_color_dict, output_prefix=args.output_prefix)
busco_len_df = prepare_len_df(busco_status_df, output_prefix="{0}".format(args.output_prefix))

informative_busco_status_df = get_informative_buscos(busco_status_df)
informative_len_df = prepare_len_df(informative_busco_status_df, output_prefix="{0}.informative".format(args.output_prefix))
informative_count_track_merged_df = create_track_from_df(informative_busco_status_df, busco_color_dict,
                                                         output_prefix="{0}.informative".format(args.output_prefix))
#print(informative_busco_status_df)
#print(informative_len_df)

offset = min(count_track_merged_df[count_track_merged_df["start"] == 0]["end"])

busco_no_complete_len_df = pd.DataFrame.from_dict({label: number_of_buscos - offset for label in busco_table_dict},
                                                  orient='index', columns=["length"])
busco_no_complete_len_df.index.name = "scaffold"
busco_no_complete_len_df.to_csv("{0}.no_complete.busco.len".format(args.output_prefix), sep="\t",
                                header=False, index=True)

count_track_merged_df[["start", "end"]] = count_track_merged_df[["start", "end"]] - offset
count_track_merged_df.loc[count_track_merged_df["start"] < 0, "start"] = 0
no_complete_track_merged_df = count_track_merged_df = count_track_merged_df[count_track_merged_df["end"] > 0]
no_complete_track_merged_df.to_csv("{0}.no_complete.busco.counts.bedgraph".format(args.output_prefix), sep="\t",
                                   header=False, index=False)
