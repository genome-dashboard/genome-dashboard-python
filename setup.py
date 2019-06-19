#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always use setuptools over distutils.
from setuptools import setup, find_packages
from os import path
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import
from io import open


# here = path.abspath(path.dirname(__file__))

# with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#     long_description = f.read()

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


# Configuration for package when publishing.
# Edit these values to reflect yourpackage details.
current_version = '0.0.22'    # -->>> !!!! IMPORTANT: BUMP VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  !!!! <<<--
module_name = 'genome_dashboard'     # Using Python convention of a hyphen in package name but an underscore in build name used during installation.
module_description = "Genome Dashboard is the logic behind a web-based prototype of a genomics dashboard, specifically designed to integrate informatics and 4D material studies of chromatin. Genome Dashboard unites our Interactive Chromatin Modeling (ICM) tools with the Biodalliance genome browser and the JSMol molecular viewer to rapidly fold any DNA sequence into atomic or coarse-grained models of DNA, nucleosomes or chromatin."
module_python = '>=2.7'
module_authors = 'Zilong Li, Ran Sun, Thomas C. Bishop'
module_authors_email = 'zli007@latech.edu, rsu007@latech.edu, bishop@latech.edu'
module_long_description = readme + '\n\n' + history
# 'text/plain', 'text/x-rst', or 'text/markdown'
module_long_description_content_type = 'text/x-rst'
module_license = "MIT license"
module_url = 'https://github.com/genomeDashboard/genome-dashboard'
module_keywords = 'python biology genomics'
module_data_included = True
module_enable_compression = False
module_test_suite = 'tests'
module_includes = [
    'genome_dashboard.gdash',
    'genome_dashboard.htp',
    'genome_dashboard.ui',
]
module_excludes = ['contrib', 'docs', 'tests']
requirements = ['Click>=6.0'] # , 'peppercorn'
setup_requirements = [ ]
test_requirements = [ ]

module_classifiers = [
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

module_entry_points = {
    'console_scripts': [
        'genome_dashboard=genome_dashboard.cli:main',
    ],
}

module_package_data = {
    '': ['*.txt'],
    'genome_dashboard': ['data/*.dat'],
}

module_extras_require = {
    'dev': ['check-manifest'],
    'test': ['coverage'],
}


# Setup method to publish package.
# DO NOT EDIT BELOW THIS LINE.
setup(
    name=module_name,
    version=current_version,
    description=module_description,
    packages=find_packages(include=module_includes, exclude=module_excludes),
    python_requires=module_python,
    author=module_authors,
    author_email=module_authors_email,
    long_description=module_long_description,
    long_description_content_type=module_long_description_content_type,
    license=module_license,
    url=module_url,
    classifiers=module_classifiers,
    keywords=module_keywords,
    install_requires=requirements,
    extras_require=module_extras_require,
    package_data=module_package_data,
    entry_points=module_entry_points,
    include_package_data=module_data_included,
    zip_safe=module_enable_compression,
    setup_requires=setup_requirements,
    test_suite=module_test_suite,
    tests_require=test_requirements,
)
