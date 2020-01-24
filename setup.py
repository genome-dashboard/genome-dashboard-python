"""A setuptools based setup module.

See:
    https://packaging.python.org/guides/distributing-packages-using-setuptools/
    https://github.com/pypa/sampleproject
"""

"""
PACKAGE PUBLISHING STEPS:

    1. Update the ./HISTORY.md file with the latest release notes.
    2. Increment the version number.
        - Approach 1:
            -- Increment the version number in ./VERSION.
            ** Keeps throwing install error on VERSION file, manifest should resolve this.
        - Approach 2:
            -- Increment the version number in docs/version.md manually.
            -- Increment the version number in setup.py (below) manually.
            ** These must match!!!
    4. Build the package:
        > python setup.py sdist
    5. Check with package:
        > twine check dist/*
    6. Test PyPI upload:
        > twine upload --repository-url https://test.pypi.org/legacy/ dist/*
    7. Upload to PyPI:
        > twine upload --repository-url https://upload.pypi.org/legacy/ dist/*

    Notes:
    - Uploads require a PyPI user account.
    - Use the twine keyring feature for cli credential management locally: https://pypi.org/project/twine/
    - Use the pypi-cli package to inspect package info and status: https://pypi.org/project/pypi-cli/
        ex. > pypi info genomedashboard
"""


# Favor setuptools over distutils.
# io.open is needed for projects that support Python 2.7
# It ensures open() defaults to text mode with universal newlines,
# and accepts an argument to specify the text encoding
# Python 3 only projects can skip this import

from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'VERSION'), encoding='utf-8') as f:
    mod_version = f.read()

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

setup(
    name='genomedashboard',
    version=mod_version,
    description='Genome Dashboard is the logic behind a web-based prototype of a genomics dashboard, specifically designed to integrate informatics and 4D material studies of chromatin. Genome Dashboard unites our Interactive Chromatin Modeling (ICM) tools with the Biodalliance genome browser and the JSMol molecular viewer to rapidly fold any DNA sequence into atomic or coarse-grained models of DNA, nucleosomes or chromatin.',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    scripts=['src/genomedashboard.py'],
    python_requires='>=2.7, ==3.*, <4',
    author='Genome Dashboard Team',
    author_email='genome.dashboard@gmail.com',
    long_description=readme,
    long_description_content_type='text/markdown',
    license='Louisiana Tech University License',
    url='http://dna.engr.latech.edu/~gdash/GDash-landing-page/',
    download_url='https://pypi.org/project/genomedashboard/#files',
    project_urls={
        'PyPI': 'https://pypi.org/project/genomedashboard/',
        'Documentation': 'https://genomedashboard.readthedocs.io/en/latest/readme.html',
        'Source Code': 'https://github.com/genomeDashboard/genomedashboard',
        'Issue Tracker': 'https://github.com/genome-dashboard/genome-dashboard-python/issues',
        'Demo': 'http://dna.engr.latech.edu/~gdash/',
        # 'Funding': 'https://donate.pypi.org',
        # 'Say Thanks!': 'http://saythanks.io/to/example',
    },
    classifiers=[
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
        'Programming Language :: Python :: 3.8',
    ],
    keywords='python biology genomics genome dashboard',
    install_requires=[
        'docutils>=0.3',
        'click>=6.0',
        'twobitreader',
        'pyBigWig',
        'numpy',
        'scipy',
        'matplotlib',
    ],
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },
    package_data={
        '': ['data/*.dat', 'data/*.txt'],
    },
    entry_points={
        'console_scripts': ['genomedashboard=src.genomedashboard:main'],
    },
    include_package_data=True,
    zip_safe=False,
    setup_requires=[ ],
    test_suite='tests',
    tests_require=[
        'unittest',
        'click>=6.0',
        'numpy',
    ]
)
