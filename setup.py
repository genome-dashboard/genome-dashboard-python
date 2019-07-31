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


print("\n>>> STARTING PACKAGE SETUP <<<\n")

print("\n... PARSING CONFIG FILES ...\n")

print("\n... PARSING DESCRIPTION RST ...\n")
with open('DESCRIPTION.rst') as description_file:
    description_text = str(description_file.read())

print("\n... PARSING README RST ...\n")
with open('README.rst') as readme_file:
    readme = readme_file.read()

print("\n... PARSING HISTORY RST ...\n")
with open('HISTORY.rst') as history_file:
    history = history_file.read()

# For Markdown.
# here = path.abspath(path.dirname(__file__))

# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
# Get the current version from a single source of truth.

# print("\n... PARSING VERSION RST ...\n")
# with open('VERSION') as version_file:
#     current_version = str(version_file.read())

# print("\n... PARSING VERSION MD ...\n")
# with open(path.join(here, 'VERSION.md'), encoding='utf-8') as f:
#     current_version = str(f.read())

# print("\n... PARSING DESCRIPTION MD ...\n")
# with open(path.join(here, 'DESCRIPTION.md'), encoding='utf-8') as f:
#     description_text = f.read()

# print("\n... PARSING LONG_DESCRIPTION MD ...\n")
# with open(path.join(here, 'LONG_DESCRIPTION.md'), encoding='utf-8') as f:
#     long_description_text = f.read()

# print("\n... PARSING README MD ...\n")
# with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#     readme = f.read()

# print("\n... PARSING HISTORY MD ...\n")
# with open(path.join(here, 'HISTORY.md'), encoding='utf-8') as f:
#     history = f.read()

print("\n... ASSIGNING CONFIGURATION VALUES ...\n")

# Configuration for package when publishing.
# Edit these values to reflect your package details.
# -->>> !!!! IMPORTANT: BUMP THE VERSION WITH EVERY COMMIT USING SEMVER CONVENTIONS  <Major.minor.patch> !!!! <<<--
module_version = '0.0.34'
module_name = 'genomedashboard'
module_authors = 'Zilong Li, Ran Sun, Thomas C. Bishop'
module_authors_email = 'zli007@latech.edu, rsu007@latech.edu, bishop@latech.edu'
module_license = "MIT license"
module_url = 'https://github.com/genomeDashboard/genomedashboard'
module_keywords = 'python biology genomics genome dashboard'
module_python = '>=3.5'
module_description = description_text

# 'text/plain', 'text/x-rst', or 'text/markdown'
module_long_description_content_type = 'text/x-rst'
# module_long_description = readme + '\n\n' + history
# module_long_description = readme
module_long_description = "THIS DOESNT WORK IN BUILD"

module_data_included = True
module_enable_compression = False
module_test_suite = 'tests'

module_includes = [
    'genomedashboard',
    'genomedashboard.htp',
    'genomedashboard.ui',
]

module_excludes = [
    'contrib',
    'docs',
    'tests'

]
module_install_requires = [
    'Click>=6.0',
    # 'peppercorn',
]

module_setup_requires = [ ]

module_test_requires = [ ]

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
        'genomedashboard=genomedashboard.cli:main',
    ],
}

module_package_data = {
    '': ['*.txt'],
    'genomedashboard': ['data/*.dat'],
}

module_extras_require = {
    'dev': ['check-manifest'],
    'test': ['coverage'],
}

module_options = {
    'bdist_wheel':{'universal':True}
}


print("\n... BUILDING PACKAGE ...\n")

# Setup method to publish package.
# DO NOT EDIT BELOW THIS LINE.
setup(
    name=module_name,
    version=module_version,
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
    install_requires=module_install_requires,
    extras_require=module_extras_require,
    package_data=module_package_data,
    entry_points=module_entry_points,
    include_package_data=module_data_included,
    zip_safe=module_enable_compression,
    setup_requires=module_setup_requires,
    test_suite=module_test_suite,
    tests_require=module_test_requires,
    options=module_options,
)

print("\n>>> PACKAGE SETUP FINISHED <<<\n")
