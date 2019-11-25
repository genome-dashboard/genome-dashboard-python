# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""

# import genomedashboard.cli as cli
# from genomedashboard.convert import convert as cv
# from genomedashboard.ds import ds as ds
# from genomedashboard.io import io as io
# from genomedashboard.mathfunction import mathfunction as mf
#
# __all__ = ["cli", "convert", "ds", "io", "mathfunction"]

import os
import pkgutil


__all__ = []

for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):
    __all__.append(module_name)
    _module = loader.find_module(module_name).load_module(module_name)
    globals()[module_name] = _module

# Get path for file.
current_dir = os.path.dirname(os.path.abspath(__file__))

def read_text_lines(fname):
    with io.open(os.path.join(current_dir, fname)) as fd:
        return fd.readlines()

# Get version info.
version_file = open(os.path.join(current_dir, 'VERSION'))
version = version_file.read().strip()
# print("*** Current Version: ", version, "\n")


# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
# This value MUST be aligned with the value in .genome-dashboard-python/setup.py!!!
# __version__ = '0.0.54'
# Get value form VERSION file.
__version__ = str(version)
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
__name__ = 'genomedashboard'
