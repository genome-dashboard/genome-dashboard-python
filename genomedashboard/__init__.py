# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

# from os import path
# from io import open
# from genomedashboard import genomedashboard as gd
# import cli as cli
from .convert import convert as cv
from .ds import ds as ds
from .io import io as io


# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
__version__ = '0.0.34'
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
