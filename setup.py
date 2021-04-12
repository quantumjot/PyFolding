#!/usr/bin/env python
from setuptools import find_packages, setup

def get_install_required():
    with open("./requirements.txt", "r") as reqs:
        requirements = reqs.readlines()
    return [r.rstrip() for r in requirements]

setup(
    name='PyFolding',
    version='1.0',
    description='PyFolding is a simple Python based framework for fitting \
    kinetic and thermodynamic models to protein folding data. It includes \
    several basic models and is extensible to enable fitting of more \
    interesting models.',
    author='Alan R. Lowe',
    author_email='a.lowe@ucl.ac.uk',
    url='https://github.com/quantumjot/PyFolding',
    packages=find_packages(),
    install_requires=get_install_required(),
    python_requires=">=3.6",
)
