# -*- coding: utf-8 -*-

"""Top-level package for Genome Dashboard"""


from os import path
from io import open
from . import gdash


# Get the current version from a single source of truth.
# with open('../VERSION') as version_file:
#     current_version = str(version_file.read())
#     print(current_version)


# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
# Get the current version from a single source of truth.
# here = path.abspath(path.dirname(__file__))
# with open(path.join(here, '../VERSION.md'), encoding='utf-8') as e:
#     current_version = str(e.read())
    # print(current_version)


# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
# __version__ =  str(current_version)

__version__ = '0.0.30'
__author__ = 'Zilong Li, Ran Sun, Thomas C. Bishop'
__email__ = 'genome.dashboard@gmail.com'
