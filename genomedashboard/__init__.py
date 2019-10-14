# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""


# ABSOLUTE IMPORTS. WORKS.
import genomedashboard.cli as cli
from genomedashboard.convert import *
from genomedashboard.ds import *
from genomedashboard.io import *

# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
__version__ = '0.0.45'   # This value MUST be aligned with the value in .genome-dashboard-python/setup.py!!!
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
__name__ = 'genomedashboard'
