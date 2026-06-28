#!/usr/bin/env python3
"""
Setup script for compGWAS - Comprehensive GWAS Analysis Toolkit
"""

from setuptools import setup, find_packages
import os
import sys

# Read the contents of README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Minimum version requirements
PYTHON_REQUIRES = '>=3.7'

# List of required Python packages
INSTALL_REQUIRES = [
    'numpy>=1.16.0',
    'pandas>=0.24.0',
]

# Optional dependencies for development
EXTRAS_REQUIRE = {
    'dev': [
        'pytest>=5.0',
        'pytest-cov>=2.8',
    ],
}

setup(
    name='compGWAS',
    version='1.0.0',
    description='A python toolkit for comprehensive GWAS analysis of SNPs/Indels',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='BaoCodeLab',
    author_email='',
    url='https://github.com/BaoCodeLab/compGWAS',
    license='GPLv3',
    
    # Packages
    packages=find_packages(),
    py_modules=['compGWAS'],
    
    # Include package data
    include_package_data=True,
    package_data={
        'allGWAS': ['GWASlib/*.R', 'preGWASlib/*.py', 'GWASannolib/*.py', 
                    'GWASLDlib/*.py', 'GWASfilterlib/*.py'],
    },
    
    # Python version requirement
    python_requires=PYTHON_REQUIRES,
    
    # Dependencies
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    
    # Entry point for CLI
    entry_points={
        'console_scripts': [
            'compGWAS=compGWAS:main',
        ],
    },
    
    # Classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    
    keywords=[
        'GWAS',
        'genomics',
        'SNP',
        'Indel',
        'genome-wide association',
        'bioinformatics',
    ],
    
    zip_safe=False,
)
