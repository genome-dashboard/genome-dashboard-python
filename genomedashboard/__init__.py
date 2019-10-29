# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

# from os import path
# from io import open
# from genomedashboard import genomedashboard as gd

from genomedashboard.cli import cli as cli
from genomedashboard.convert import convert as cv
from genomedashboard.ds import ds as ds
from genomedashboard.io import io as io

# from . import cli, convert, ds, io
# from . import cli
# from . import convert
# from . import ds
# from . import io

#import genomedashboard.cli as cli
#import genomedashboard.convert as cv
#import genomedashboard.ds as ds
#import genomedashboard.io as io

# __all__ = ["cli", "convert", "ds", "io"]

# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
__version__ = '0.0.41'   # This value MUST be aligned with the value in .genome-dashboard-python/setup.py!!!
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
