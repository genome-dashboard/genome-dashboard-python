# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

# from os import path
# from io import open
# from genomedashboard import genomedashboard as gd

# from . import cli as cli
# from . import convert as cv
# from . import ds as ds
# from . import io as io

# from . import cli, convert, ds, io
from . import cli
from . import convert
from . import ds
from . import io



# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
__version__ = '0.0.39'   # This value MUST be aligned with the value in .genome-dashboard-python/setup.py!!!
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'

__all__ = ["cli", "convert", "ds", "io"]
