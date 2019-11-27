# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

import os
from pathlib import Path

from . src import cli as cli
from . src import convert as cv
from . src import ds as ds
from . src import io as io
from . src import mathfunctions as mf

__all__ = ["cli", "cv", "ds", "io", "mf"]

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = (current_dir).parent

# Get version info.
version_file = open(os.path.join(parent_dir, 'VERSION'))
version = str(version_file.read().strip())
print("*** __init__.py - Current Version: ", version, "\n")

# Get authors info.
authors_file = open(os.path.join(parent_dir, 'AUTHORS'))
authors = str(authors_file.read().strip())
print("*** __init__.py - Authors: ", authors, "\n")

# Get value form VERSION file. Avoids syncing between multiple locations problem.
__version__ = version
__author__ = authors   # 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
__name__ = 'genomedashboard'
