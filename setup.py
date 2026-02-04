#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

exec(open('pygait2d/version.py').read())

extension = Extension(
    name="algait2de.gait2de",
    sources=[
        os.path.join("algait2de", "gait2de.pyx"),
        os.path.join("algait2de", "gait2de_al.c"),
    ],
    include_dirs=[numpy.get_include()],
)

description = "An implementation of a planar human gait model."

setup(
    name='gait2d',
    author='Jason K. Moore',
    author_email='moorepants@gmail.com',
    version=__version__,
    url="http://github.com/csu-hmc/gait2d",
    description=description,
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=[
        'pydy',
        'pyyaml',
        'sympy',
    ],
    extras_require={
        'examples': ['numpy', 'scipy', 'cython'],
        'doc': ['sphinx', 'numpydoc'],
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
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
