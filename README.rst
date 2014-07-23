Gait2de is a dynamic model that contains the essential elements for simulating
human gait in the sagittal plane.  The model has been used in earlier forms by
Ackermann and van den Bogert (2010) and by Gerritsen et al. (1998).  The main
features of the model are:

- Seven body segments, sixteen muscles
- A Matlab programming interface
- Fast execution, ~0.03 ms to compute the state derivatives and other model outputs

Model dynamics and outputs are twice differentiable with respect to all inputs,
which is important for certain numerical methods for simulation and optimal
control The model is intended for education and basic research.  We envision
that the model can be used for studies such as Ackermann and van den Bogert
(2010) and Geyer and Herr (2010).

For more complex 3D models, see, for instance, the Opensim project
(www.simtk.org/home/opensim).  These more complex models will probably execute
significantly slower and may not be twice differentiable.

algait2de
---------

To manually build the model and use in Python::

  $ al gait2de.al
  $ python autolevclean.py
  $ python setup.py build_ext --inplace
  $ python
  >>> import gait2de
  >>> gait2de.evaluate_autolev_rhs(...)
