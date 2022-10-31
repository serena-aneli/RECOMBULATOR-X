#! /usr/bin/env python
#
# Copyright (C) 2012-2020 Michael Waskom

if False:
    DESCRIPTION = "seaborn: statistical data visualization"
    LONG_DESCRIPTION = """\
    Seaborn is a library for making statistical graphics in Python. It is built on top of `matplotlib <https://matplotlib.org/>`_ and closely integrated with `pandas <https://pandas.pydata.org/>`_ data structures.

    Here is some of the functionality that seaborn offers:

    - A dataset-oriented API for examining relationships between multiple variables
    - Specialized support for using categorical variables to show observations or aggregate statistics
    - Options for visualizing univariate or bivariate distributions and for comparing them between subsets of data
    - Automatic estimation and plotting of linear regression models for different kinds dependent variables
    - Convenient views onto the overall structure of complex datasets
    - High-level abstractions for structuring multi-plot grids that let you easily build complex visualizations
    - Concise control over matplotlib figure styling with several built-in themes
    - Tools for choosing color palettes that faithfully reveal patterns in your data

    Seaborn aims to make visualization a central part of exploring and understanding data. Its dataset-oriented plotting functions operate on dataframes and arrays containing whole datasets and internally perform the necessary semantic mapping and statistical aggregation to produce informative plots.
    """

    MAINTAINER = 'Michael Waskom'
    MAINTAINER_EMAIL = 'mwaskom@gmail.com'
    URL = 'https://seaborn.pydata.org'
    LICENSE = 'BSD (3-clause)'
    DOWNLOAD_URL = 'https://github.com/mwaskom/seaborn/'
DISTNAME = 'recombulatorx'
VERSION = '0.0.0'
PYTHON_REQUIRES = ">=3.7"

INSTALL_REQUIRES = [
    'numpy>=1.15.4',
    'scipy>=1.2.0',
    'networkx>=2.0',
    'pandas>=0.23.4',
]
EXTRAS_REQUIRE = {
    'all': [
        'numba>=0.49.0'
        #'matplotlib>=1.0',
    ]
}


'''
python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install setuptools wheel
pip install numpy==1.15.4 scipy==1.2.0 networkx==2.0 pandas==0.23.4
python3 xstr_recomb/cli.py testsim.tsv

pip install numba==0.49.0

'''

if False:
    PACKAGES = [
        'seaborn',
        'seaborn.colors',
        'seaborn.external',
        'seaborn.tests',
    ]

    CLASSIFIERS = [
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Multimedia :: Graphics',
        'Operating System :: OS Independent',
        'Framework :: Matplotlib',
    ]


if __name__ == "__main__":

    from setuptools import setup

    import sys
    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("recombulatorx requires python >= 3.7.")

    setup(
        name=DISTNAME,
#        author=MAINTAINER,
#        author_email=MAINTAINER_EMAIL,
#        maintainer=MAINTAINER,
#        maintainer_email=MAINTAINER_EMAIL,
#        description=DESCRIPTION,
#        long_description=LONG_DESCRIPTION,
#        license=LICENSE,
#        url=URL,
        version=VERSION,
#        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        extras_require=EXTRAS_REQUIRE,
#        packages=PACKAGES,
#        classifiers=CLASSIFIERS
    )
