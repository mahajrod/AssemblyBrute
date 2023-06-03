#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for manipulations with options inside Snakemake rules
"""


def parse_option(option, option_dict, option_prefix, default_value="default"):
    if option not in option_dict:
        return ""
    if option_dict[option] is None:
        return ""
    if option_dict[option] == default_value:
        return ""
    return " {0} {1}".format(option_prefix, option_dict[option])


def parse_option_flag(option, option_dict, option_prefix):
    if option not in option_dict:
        return ""
    if option_dict[option]:
        return " {0}".format(option_prefix)
    else:
        return ""
