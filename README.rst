Introduction
============

Gait2D is a dynamic model that contains the essential elements for simulating
human gait in the sagittal plane. The model has been used in earlier forms by
Ackermann and van den Bogert (2010) and by Gerritsen et al. (1998). The main
features of the model are:

- Seven body segments
- Fast execution

Model dynamics and outputs are twice differentiable with respect to all inputs,
which is important for certain numerical methods for simulation and optimal
control. The model is intended for education and basic research. The model can
be used for studies such as Ackermann and van den Bogert (2010) and Geyer and
Herr (2010).

Speed
=====

The cythonized Autolev C code takes about 30 micro seconds per rhs eval and the
pydy cython version takes about 70 microseconds (the slow part is, of course,
the Python level solve on the full mass matrix).

Usage
=====

pygait2d
--------

See ``example/run.py``.

algait2de
---------

To manually build the Autolev model and make use of it in Python::

   $ cd algait2d
   $ al gait2de.al
   $ python autolevclean.py
   $ python setup.py build_ext --inplace
   $ python
   >>> import gait2de
   >>> gait2de.evaluate_autolev_rhs(...)
