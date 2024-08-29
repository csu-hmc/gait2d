#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy
from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize

exec(open('pygait2d/version.py').read())

extension = Extension(
    name="gait2de",
    sources=[
        os.path.join("algait2de", "gait2de.pyx"),
        os.path.join("algait2de", "gait2de_al.c"),
    ],
    include_dirs=[numpy.get_include()],
)

description = "An implementation of a planar human gait model."

setup(
    name='Gait2D',
    author='Jason K. Moore',
    author_email='moorepants@gmail.com',
    version=__version__,
    url="http://github.com/csu-hmc/gait2d",
    description=description,
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=[
        'sympy',
        'pydy',
        'pyyaml',
    ],
    extras_require={
        'examples': ['numpy', 'scipy', 'cython'],
        'doc': ['sphinx>=1.1.0', 'numpydoc>=0.4'],
    },
    tests_require=['pytest'],
    ext_modules=cythonize([extension]),
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
