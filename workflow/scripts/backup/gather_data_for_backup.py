#!/usr/bin/env python
__author__ = 'mahajrod'

import sys
import argparse

from pathlib import Path
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-r", "--results_dir", action="store", dest="results_dir", required=True,
                    help="Directory containing results of the pipeline")


args = parser.parse_args()


