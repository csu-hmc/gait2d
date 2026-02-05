Introduction
============

gait2d is a dynamics model that contains the essential elements for simulating
human gait in the sagittal plane. The model has been used in earlier forms by
Gerritsen et al. (1998) [Gerritsen1998]_ and Ackermann and van den Bogert
(2010) [Ackermann2010]_. The main features of the model are:

- Seven body segments
- Hunt-Crossley style foot contact model
- Sixteen musculotendon actuators
- Differentiable dynamics
- Fast execution

Model dynamics and outputs are twice differentiable with respect to all inputs
and the numerical model is designed to execute at high computational
efficiency, both of which are important for certain numerical methods for
simulation and optimal control. The model is intended for education and
research. The model can be used for studies such as Ackermann and van den
Bogert (2010) [Ackermann2010]_ and Geyer and Herr (2010) [Geyer2010]_.

The included Autolev dynamics model was originally written by Ton van den
Bogert and adapted by Jason K. Moore. The SymPy Mechanics model was created by
Jason to match the Autolev model results and then futher extended for use with
PyDy and opty. These two models are provided in two packages in this
distribution: algait2de and pygait2d.

License
=======

All files in this repository are licensed under the Apache 2.0 open source
license.

Dependencies
============

Run time dependencies:

- matplotlib
- numpy
- pydy
- python
- pyyaml
- symmeplot
- sympy

Build dependencies:

- cython [with a C compiler]
- numpy
- python
- setuptools

Development dependencies:

- pytest
- sphinx
- sphinx_gallery

Dependencies used in the examples:

- cython [with a C compiler]
- opty
- scipy

Install
=======

This package has not yet been released to PyPi, so install from the development
source::

   git clone https://github.com/csu-hmc/gait2d
   cd gait2d
   conda env create -f gait2d-dev.yml
   conda activate gait2d-dev
   python -m pip install --no-deps --no-build-isolation --editable .

Usage
=====

pygait2d
--------

See ``examples/plot_forward_sim.py`` for a basic use example.

algait2de
---------

To manually rebuild the Autolev model and make use of it in Python::

   $ cd algait2d
   $ al gait2de.al
   $ python autolevclean.py
   $ cd ..
   $ python -m pip install --no-deps --no-build-isolation --editable .

The above will install algait2de and then it is importable and used with::

   $ python
   >>> from algait2de import gait2de
   >>> gait2de.evaluate_autolev_rhs(...)

References
==========

.. [Gerritsen1998] Gerritsen, K. G. M., Bogert, A. J. van den, Hulliger, M., &
   Zernicke, R. F.  (1998). Intrinsic Muscle Properties Facilitate Locomotor
   Control—A Computer Simulation Study. Motor Control, 2(3), 206–220.
   https://doi.org/10.1123/mcj.2.3.206
.. [Ackermann2010] Ackermann, M., & van den Bogert, A. J. (2010). Optimality
   principles for model-based prediction of human gait. Journal of
   Biomechanics, 43(6), 1055–1060.
   https://doi.org/10.1016/j.jbiomech.2009.12.012
.. [Geyer2010] Geyer, H., & Herr, H. (2010). A Muscle-Reflex Model that Encodes
   Principles of Legged Mechanics Produces Human Walking Dynamics and Muscle
   Activities. Neural Systems and Rehabilitation Engineering, IEEE Transactions
   On, 18(3).  https://doi.org/10.1109/TNSRE.2010.2047592
