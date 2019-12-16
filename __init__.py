# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

import os
from pathlib import Path

# APPROACH 1.
# from . src import cli as cli
from . src import convert as cv
from . src import ds as ds
from . src import io as io
from . src import mathfunctions as mf

# __all__ = ["cv", "ds", "io", "mf"] # "cli",
__all__ = ["*"]

# APPROACH 2.
# from . src import *

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = (current_dir).parent

# Get version info.
# version_file = open(os.path.join(parent_dir, 'VERSION'))
# version = str(version_file.read().strip())
with open('VERSION') as version_file:
    version = version_file.read().strip()
    print("*** __init__.py - Current Version: ", version, "\n")

# Get authors info.
# authors_file = open(os.path.join(parent_dir, 'AUTHORS'))
# authors = str(authors_file.read().strip())
with open('AUTHORS') as authors_file:
    authors = authors_file.read().strip()
    # print("*** __init__.py - Authors: ", authors, "\n")

# Get value form VERSION file. Avoids syncing between multiple locations problem.
__version__ = version
__author__ = authors   # 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
__name__ = 'genomedashboard'
