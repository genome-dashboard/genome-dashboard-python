# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

import genomedashboard.cli as cli
from genomedashboard.convert import convert as cv
from genomedashboard.ds import ds as ds
from genomedashboard.io import io as io
from genomedashboard.mathfunction import mathfunction as mf

__all__ = ["cli", "convert", "ds", "io", "mathfunction"]

# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
__version__ = '0.0.52'   # This value MUST be aligned with the value in .genome-dashboard-python/setup.py!!!
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
__name__ = 'genomedashboard'
