"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import
from io import open


#######################################
############# START ADDED.

from pathlib import Path

here = path.abspath(path.dirname(__file__))
current_dir = Path(here)
parent_dir = current_dir.parent

with open('VERSION') as version_file:
    current_version = version_file.read().strip()
    # print(current_version)

# with open('AUTHORS.rst') as authors_file:
#     authors = authors_file.read()
    # print(authors)

with open(path.join(here, 'AUTHORS.md'), encoding='utf-8') as f:
    authors = f.read()
    # print(authors)

# with open('DESCRIPTION.rst') as description_file:
#     description = description_file.read()
    # print(description)

with open(path.join(here, 'DESCRIPTION.md'), encoding='utf-8') as f:
    description = f.read()
    # print(description)

# with open('README.rst') as readme_file:
#     readme = readme_file.read()
    # print(readme)

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    readme = f.read()
    # print(readme)

# with open('HISTORY.rst') as history_file:
#     history = history_file.read()
#     print(history)

with open(path.join(here, 'HISTORY.md'), encoding='utf-8') as f:
    history = f.read()
    # print(history)

# unsure of this...
linebreak = str('  ')   # linebreaks in markdown: two spaces, backslash, <br/>, &nbsp;
# long_description = str(description + linebreak + readme + linebreak + history)
long_description = str(description + readme + history)
print(long_description)

############# END ADDED.
#######################################


# Get the long description from the README file
# with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#     long_description = f.read()


# ASSIGN VALUES.
module_name                             = 'genomedashboard'
module_version                          = current_version
# module_authors                          = 'Zilong Li, Ran Sun, Thomas C. Bishop'
module_authors                          = authors
module_authors_email                    = 'zli007@latech.edu, rsu007@latech.edu, bishop@latech.edu, taoteg@gmail.com'
module_license_type                     = "MIT license"
module_url                              = 'http://dna.engr.latech.edu/~gdash/GDash-landing-page/'
module_download_url                     = 'https://pypi.org/project/genomedashboard/#files'
module_keywords                         = 'python biology genomics genome dashboard'
module_python                           = '>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4'
module_description                      = description
# module_long_description                 = description + '\n\n' + readme + '\n\n' + history
module_long_description                 = long_description
# module_long_description                 = description   # Cant get text syntax to pass properly.
module_long_description_content_type    = 'text/markdown'   # 'text/plain',  'text/markdown' or 'text/x-rst'.
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
module_packages                         = find_packages(where = 'src')

module_scripts                          = ['src/genomedashboard.py']

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
                                                'genomedashboard=src.genomedashboard:main',
                                            ],
                                        }

module_package_data                     = {
                                            '': ['*.txt', '*.rst'],
                                            '': ['data/*.dat','data/*.txt'],
                                        }

module_extras_require                   = {
                                            'dev': ['check-manifest'],
                                            'test': ['coverage'],
                                        }


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
    tests_require                   = module_test_requires
)

"""
# USAGE
# Build the package:    > python setup.standard.py sdist bdist
# Check with package:   > twine check dist/*
# Test PyPI upload:     > twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# Upload to PyPI:       > twine upload --repository-url https://upload.pypi.org/legacy/ dist/*
# Note: uploads require a PyPI user account.
"""
