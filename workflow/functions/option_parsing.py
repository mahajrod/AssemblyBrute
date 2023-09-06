#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for manipulations with options inside Snakemake rules
"""


def parse_option(option, option_dict, option_prefix, default_value="default", none_value=None, expression=None):
    if option not in option_dict:
        if none_value is None:
            return ""
        else:
            " {0} {1}".format(option_prefix, none_value)
    if option_dict[option] is None:
        if none_value is None:
            return ""
        else:
            " {0} {1}".format(option_prefix, none_value)
    if option_dict[option] == default_value:
        return ""
    return " {0} {1}".format(option_prefix, expression(option_dict[option]) if expression else option_dict[option])


def parse_option_flag(option, option_dict, option_prefix):
    if option not in option_dict:
        return ""
    if option_dict[option]:
        return " {0}".format(option_prefix)
    else:
        return ""


def group_option_sets(option_set_dict, grouping_option_list):
    option_set_group_dict = {}
    option_group_dict = {}
    for option_set in option_set_dict:
        option_group_label_list = []
        for grouping_option in grouping_option_list:
            if grouping_option in option_set_dict[option_set]:
                option_group_label_list.append("{0}_{1}".format(grouping_option, option_set_dict[option_set][grouping_option]))
            else:
                option_group_label_list.append("{0}_{1}".format(grouping_option, "default"))
        option_group_label = ".".join(option_group_label_list)
        option_set_group_dict[option_set] = option_group_label
        if option_group_label in option_group_dict:
            option_group_dict[option_group_label].append(option_set)
        else:
            option_group_dict[option_group_label] = [option_set]

    option_group_syn_dict = {"option_set_{0}".format(index): group for group, index in zip(option_group_dict,
                                                                                           range(1, len(option_group_dict) + 1))}
    final_dict = {}
    for option_group_syn in option_group_syn_dict:
        final_dict[option_group_syn] = {}
        final_dict[option_group_syn]["option_set_list"] =  option_group_dict[option_group_syn_dict[option_group_syn]]
        final_dict[option_group_syn]["grouping_options"] = {}
        for grouping_option in grouping_option_list:
            if grouping_option in option_set_dict[final_dict[option_group_syn]["option_set_list"][0]]:
                 final_dict[option_group_syn]["grouping_options"][grouping_option] = option_set_dict[final_dict[option_group_syn]["option_set_list"][0]][grouping_option]
            else:
                final_dict[option_group_syn]["grouping_options"][grouping_option] = "default"

    final_option_set_group_dict = {}

    for option_group_syn in final_dict:
        for option_set in final_dict[option_group_syn]["option_set_list"]:
            final_option_set_group_dict[option_set] = option_group_syn

    return final_dict, final_option_set_group_dict


def parse_node_list(rulename, grid_system="slurm"):
    black_list = set(config["nodes"]["blacklist"])
    white_list = set(config["nodes"]["whitelist"])
    if rulename in config["rule_nodes"]:
        if "blacklist" in config["rule_nodes"][rulename]:
            black_list = set(config["rule_nodes"][rulename]["blacklist"]) if config["rule_nodes"][rulename]["blacklist"] else black_list
        if "whitelist" in config["rule_nodes"][rulename]:
            white_list = set(config["rule_nodes"][rulename]["whitelist"]) if config["rule_nodes"][rulename]["whitelist"] else white_list
    if grid_system == "slurm":
        whitelist_option = " --nodelist={0} ".format(",".join(white_list)) if white_list else " "
        blacklist_option = " --exclude={0} ".format(",".join(black_list)) if black_list else " "

        return whitelist_option + blacklist_option
    else:
        print("White and black node lists are implemented only for slurm. "
              "Modify 'parse_node_list' function in 'workflow/functions/option_parsing_py' "
              "if you need such functionality for other grid systems")
