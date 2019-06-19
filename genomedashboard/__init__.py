# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""


from io import open
from . import gdash as g


# Get the current version from a single source of truth.
with open('VERSION.rst') as version_file:
    current_version = version_file.read()
    print(current_version)

# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
# __version__ =  str(current_version)
__version__ = '0.0.25'
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
