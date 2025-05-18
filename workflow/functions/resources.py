#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for Snakemake file
"""
import os
import yaml
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath


def get_memory(wildcards, attempt, start_mem, coeff=2, mode="linear"):
    if mode == "exp":
        if coeff <= 1:
            raise ValueError(f"ERROR!!! Coefficient for exponential resource selection should be above 1 and not {coeff}!")
        print([attempt, start_mem, coeff, mode])
        return int((coeff ** (attempt - 1)) * start_mem)
    elif mode == "linear":
        print([attempt, start_mem, coeff, mode])
        return attempt * start_mem
    else:
        raise ValueError("ERROR!!! Unknown mode for memory selection")


