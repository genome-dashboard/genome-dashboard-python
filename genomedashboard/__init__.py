# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

# from os import path
# from io import open
from genomedashboard import *
from cli import *
from .io import io as io
from .converter import converter as converter
from .ds import ds as ds


# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
__version__ = '0.0.33'
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
