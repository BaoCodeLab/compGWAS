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

setup(
    name='compGWAS',
    version='1.0.0',
    description='A python toolkit for comprehensive GWAS analysis of SNPs/Indels',
    author='BaoCodeLab',
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
    
    # Entry point for CLI
    entry_points={
        'console_scripts': [
            'compGWAS=compGWAS:main',
        ],
    },
    
    # Classifiers
    classifiers=[
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
    ],
    
    zip_safe=False,
)
