#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
A setuptools based setup module (replaces disutils).
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

debug = True
print(debug)

if debug:
    print("\n>>> STARTING SETUP SCRIPT <<<\n")

from setuptools import setup, find_packages
from os import path
from io import open
from pathlib import Path
from sphinx_content_filter import *

# import docutils.nodes
# import docutils.parsers.rst
# import docutils.utils
# import docutils.frontend
#
# def parse_rst(text:str) -> docutils.nodes.document:
#     parser = docutils.parsers.rst.Parser()
#     components = (docutils.parsers.rst.Parser,)
#     settings = docutils.frontend.OptionParser(components=components).get_default_values()
#     document = docutils.utils.new_document('<rst-doc>', settings=settings)
#     parser.parse(text, document)
#     return document

if debug:
    print("\n>>> GETTING PATHS <<<\n")

here = path.abspath(path.dirname(__file__))

if debug:
    print(here)

# current_dir = path.abspath(path.dirname(__file__))
current_dir = Path(here)

if debug:
    print(current_dir)

# parent_dir = path.abspath(path.dirname(__file__)).parent
parent_dir = current_dir.parent

if debug:
    print(parent_dir)

if debug:
    print("\n>>> PARSING VERSION <<<\n")

with open(path.join(current_dir, 'VERSION'), encoding='utf-8') as version_file:
    current_version = version_file.read().strip()

if debug:
    print(current_version)

if debug:
    print("\n>>> PARSING AUTHORS <<<\n")

with open(path.join(current_dir, 'AUTHORS.rst'), encoding='utf-8') as authors_file:
    authors = authors_file.read().strip()

if debug:
    print(authors)

if debug:
    print("\n>>> PARSING DESCRIPTION RST <<<\n")

# ver 1.
# def read_text_lines(fname):
#     with open(path.join(current_dir, fname)) as fd:
#         return fd.readlines()
#
# description = read_text_lines('DESCRIPTION.rst')
# print(description)
# filtered_description = ''.join(yield_sphinx_only_markup(description)),
# print(filtered_description)
# print(str(filtered_description))

# ver 2.
with open(path.join(current_dir, 'DESCRIPTION.rst'), encoding='utf-8') as description_file:
    description = description_file.read()

if debug:
    print(description)

if debug:
    print("\n>>> PARSING README RST <<<\n")

with open(path.join(current_dir, 'README.rst'), encoding='utf-8') as readme_file:
    readme = readme_file.read()

if debug:
    print(readme)

"""
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
"""

if debug:
    print("\n>>> PARSING HISTORY RST <<<\n")

with open(path.join(current_dir, 'HISTORY.rst'), encoding='utf-8') as history_file:
    history = history_file.read()

if debug:
    print(history)

if debug:
    print("\n>>> PARSING LONG_DESCRIPTION <<<\n")

linebreak = str('\n')
# long_description = description + linebreak + readme + linebreak + history
long_description = str(description + linebreak + readme + linebreak + history)
# long_description = parse_rst(description + linebreak + readme + linebreak + history)
# long_description = str(description + readme + history)

if debug:
    print(long_description)

if debug:
    print("\n>>> ASSIGNING CONFIGURATION VALUES <<<\n")

module_name                             = 'genomedashboard'
module_version                          = current_version
# module_authors                          = 'Zilong Li, Ran Sun, Thomas C. Bishop'
module_authors                          = authors
module_authors_email                    = 'zli007@latech.edu, rsu007@latech.edu, bishop@latech.edu, taoteg@gmail.com'
module_license_type                     = "MIT license"
module_url                              = 'http://dna.engr.latech.edu/~gdash/GDash-landing-page/'
module_download_url                     = 'https://pypi.org/project/genomedashboard/#files'
module_keywords                         = 'python biology genomics genome dashboard'
module_python                           = '>=2.7'
module_description                      = description
module_long_description                 = description
# module_long_description                 = description + '\n\n' + readme + '\n\n' + history
# module_long_description                 = long_description
module_long_description_content_type    = 'text/x-rst'   # 'text/plain',  'text/markdown' or 'text/x-rst'.
module_data_included                    = True
module_enable_compression               = False
module_test_suite                       = 'tests'
module_setup_requires                   = [ ]
module_test_requires                    = [ ]

module_project_urls                     = {
                                            'PyPI': 'https://pypi.org/project/genomedashboard/',
                                            'Documentation': 'https://genomedashboard.readthedocs.io/en/latest/readme.html',
                                            'Source Code': 'https://github.com/genomeDashboard/genomedashboard',
                                            'Issue Tracker': 'https://github.com/genome-dashboard/genome-dashboard-python/issues',
                                            'Demo': 'http://dna.engr.latech.edu/~gdash/',
                                            # 'Funding': 'https://donate.pypi.org',
                                            # 'Say Thanks!': 'http://saythanks.io/to/example',
                                        }

# module_includes                         = [
#                                             # 'genomedashboard',
#                                             'cli',
#                                             'genomedashboard.convert',
#                                             'genomedashboard.data',
#                                             'genomedashboard.ds',
#                                             'genomedashboard.io'
#                                         ]

# module_excludes                         = [
#                                             'contrib',
#                                             'docs',
#                                             'tests'
#                                         ]

# When your source code is in a subdirectory under the project root, e.g.
# `src/`, it is necessary to specify the `package_dir` argument.
module_package_dir                      = {'': 'src'}  # Optional

# module_packages                         = find_packages()
# module_packages                         = find_packages(exclude = module_excludes)
# module_packages                         = find_packages('genomedashboard', exclude = module_excludes)
# module_packages                         = find_packages(include = module_includes, exclude = module_excludes)
module_packages                         = find_packages(where='src'),  # Required

# module_scripts                          = ['src/cli.py','src/convert.py','src/ds.py','src/io.py','src/mathfunction.py']
module_scripts                          = ['cli.py','convert.py','ds.py','io.py','mathfunction.py']

module_install_requires                 = [
                                            'docutils>=0.3',
                                            'click>=6.0',
                                            'twobitreader',
                                            'pyBigWig',
                                            'numpy',
                                            'scipy',
                                            'matplotlib'
                                        ]

module_classifiers                      = [
                                            'Development Status :: 2 - Pre-Alpha',
                                            'Intended Audience :: Developers',
                                            'Intended Audience :: Science/Research',
                                            'Intended Audience :: Healthcare Industry',
                                            'Intended Audience :: Education',
                                            'Topic :: Software Development :: Libraries',
                                            'License :: OSI Approved :: MIT License',
                                            'Natural Language :: English',
                                            "Programming Language :: Python :: 2",
                                            'Programming Language :: Python :: 2.7',
                                            'Programming Language :: Python :: 3',
                                            'Programming Language :: Python :: 3.4',
                                            'Programming Language :: Python :: 3.5',
                                            'Programming Language :: Python :: 3.6',
                                            'Programming Language :: Python :: 3.7',
                                        ]

module_entry_points                     = {
                                            'console_scripts': [
                                                # 'genomedashboard=genomedashboard.cli:main',
                                                'genomedashboard=src.cli:main',
                                            ],
                                        }

module_package_data                     = {
                                            '': ['*.txt', '*.rst'],
                                            # 'genomedashboard': ['*.msg'],
                                            'src': ['*.msg'],
                                            # 'genomedashboard': ['data/*.dat'],
                                            'src': ['data/*.dat'],
                                        }

module_extras_require                   = {
                                            'dev': ['check-manifest'],
                                            'test': ['coverage'],
                                        }

if debug:
    print("\n>>> RUNNING SETUP METHOD <<<\n")

# DO NOT EDIT BELOW THIS LINE.
setup(
    name                            = module_name,
    version                         = module_version,
    description                     = module_description,
    package_dir                     = module_package_dir,                   # new
    packages                        = module_packages,
    scripts                         = module_scripts,                       # new
    python_requires                 = module_python,
    author                          = module_authors,
    author_email                    = module_authors_email,
    long_description_content_type   = module_long_description_content_type,
    long_description                = module_long_description,
    license                         = module_license_type,
    url                             = module_url,
    download_url                    = module_download_url,                  # new
    project_urls                    = module_project_urls,
    classifiers                     = module_classifiers,
    keywords                        = module_keywords,
    install_requires                = module_install_requires,
    extras_require                  = module_extras_require,
    package_data                    = module_package_data,
    entry_points                    = module_entry_points,
    include_package_data            = module_data_included,
    zip_safe                        = module_enable_compression,
    setup_requires                  = module_setup_requires,
    test_suite                      = module_test_suite,
    tests_require                   = module_test_requires,
)

if debug:
    print("\n>>> SETUP COMPLETED. <<<\n")
    print("\n>>> PLEASE TEST AND PUBLISH. <<<\n")

"""
# USAGE
# Build the package:    > python setup.standard.py sdist bdist
# Check with package:   > twine check dist/*
# Test PyPI upload:     > twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# Upload to PyPI:       > twine upload --repository-url https://upload.pypi.org/legacy/ dist/*
# Note: uploads require a PyPI user account.
"""
