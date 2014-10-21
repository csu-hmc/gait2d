#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

from pygait2d import __version__

description = \
    """An implementation of a planar human gait model."""

setup(name='Gait2D',
      author='Jason K. Moore',
      author_email='moorepants@gmail.com',
      version=__version__,
      url="http://github.com/csu-hmc/gait2d",
      description=description,
      license='UNLICENSE',
      packages=find_packages(),
      install_requires=['sympy',
                        'pydy',
                        'pyyaml',
                        ],
      extras_require={'examples': ['numpy', 'scipy', 'cython'],
                      'tests': ['numpy', 'cython'],
                      'doc': ['sphinx>=1.1.0',
                              'numpydoc>=0.4'],
                      },
      tests_require=['nose>1.3.0'],
      test_suite='nose.collector',
      long_description=open('README.rst').read(),
      classifiers=[
                   'Development Status :: 4 - Beta',
                   'Intended Audience :: Science/Research',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Physics',
                  ],
      )
