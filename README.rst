Introduction
============

gait2D is a dynamic model that contains the essential elements for simulating
human gait in the sagittal plane. The model has been used in earlier forms by
Ackermann and van den Bogert (2010) [Ackermann2010]_ and by Gerritsen et al.
(1998) [Gerritsen1998]_. The main features of the model are:

- Seven body segments
- Fast execution

Model dynamics and outputs are twice differentiable with respect to all inputs,
which is important for certain numerical methods for simulation and optimal
control. The model is intended for education and basic research. The model can
be used for studies such as Ackermann and van den Bogert (2010) and Geyer and
Herr (2010) [Geyer2010]_.

Speed
=====

The cythonized Autolev C code takes about 5 micro seconds per rhs eval and the
PyDy cython version takes about 15 microseconds (the slow part is, of course,
the Python level solve on the full mass matrix).

Usage
=====

pygait2d
--------

See ``examples/run.py``.

algait2de
---------

To manually build the Autolev model and make use of it in Python::

   $ cd algait2d
   $ al gait2de.al
   $ python autolevclean.py
   $ cd ..
   $ python setup.py build_ext --inplace
   $ python
   >>> import gait2de
   >>> gait2de.evaluate_autolev_rhs(...)

Model Description
=================

The model is planar (xy plane) and there are seven rigid bodies. All rotations
are defined to be about the z axis. When all configuration variables are set to
zero, the model is standing upright and the trunk's hip joint is at the
inertial origin with the human facing in the positive x direction.

Rigid bodies and constants
--------------------------

- Trunk (A):

  - ``ma``: mass [kg]
  - ``ia``: moment of inertia about mass center wrt to the trunk reference
    frame [kg*m^2]
  - ``xa``: local x location of mass center wrt to the hip joint [m]
  - ``ya``: local y location of mass center wrt to the hip joint [m]

- Right Thigh (B)

  - ``mb``: mass [kg]
  - ``ib``: moment of inertia about mass center wrt to right thigh reference
    frame [kg*m^2]
  - ``xb``: local x location of mass center wrt to the hip joint [m]
  - ``yb``: local y location of mass center wrt to the hip joint [m]
  - ``lb``: joint to joint segment length [m]

- Right Shank (C)

  - ``mc``: mass [kg]
  - ``ic``: moment of inertia about mass center wrt to right shank reference
    frame [kg*m^2]
  - ``xc``: x location of mass center wrt to the knee joint [m]
  - ``yc``: y location of mass center wrt to the knee joint [m]
  - ``lc``: joint to joint segment length [m]

- Right Foot (D)

  - ``md``: mass [kg]
  - ``id``: moment of inertia about mass center wrt to right foot reference
    frame [kg*m^2]
  - ``xd``: local x location of mass center wrt to the ankle joint [m]
  - ``yd``: local y location of mass center wrt to the ankle joint [m]
  - ``hxd``: local x location of heel wrt to the ankle joint [m]
  - ``txd``: local x location of toe wrt to the ankle joint [m]
  - ``fyd``: local y location of heel and toe relative to ankle joint [m]

- Left Thigh (E)

  - ``me``: mass [kg]
  - ``ie``: moment of inertia about mass center wrt to left thigh reference
    frame [kg*m^2]
  - ``xe``: local x location of mass center wrt to the hip joint [m]
  - ``ye``: local y location of mass center wrt to the hip joint [m]
  - ``le``: joint to joint segment length [m]

- Left Shank (F)

  - ``mf``: mass [kg]
  - ``if``: moment of inertia about mass center wrt to left shank reference
    frame [kg*m^2]
  - ``xf``: x location of mass center wrt to the knee joint [m]
  - ``yf``: y location of mass center wrt to the knee joint [m]
  - ``lf``: joint to joint segment length [m]

- Left Foot (G)

  - ``mg``: mass [kg]
  - ``ig``: moment of inertia about mass center wrt to left foot reference
    frame [kg*m^2]
  - ``xg``: local x location of mass center wrt to the ankle joint [m]
  - ``yg``: local y location of mass center wrt to the ankle joint [m]
  - ``hxg``: local x location of heel wrt to the ankle joint [m]
  - ``txg``: local x location of toe wrt to the ankle joint [m]
  - ``fyg``: local y location of heel and toe relative to ankle joint [m]

- Other constants

  - ``kc``: ground contact stiffness [N/m^3]
  - ``cc``: ground contact damping [s/m]
  - ``mu``: friction coefficient
  - ``vs``: velocity constant [m/s]
  - ``g``: acceleration due to gravity [m/s^2]

Generalized coordinates
-----------------------

- ``qax, qay``: location of trunk hip joint relative to inertial origin
- ``qa``: angle of trunk relative to inertial reference frame, ``qa=0`` makes
  trunk standing upright and ``qa>0`` leans trunk backwards
- ``qb``: angle of right thigh relative to trunk (hip), ``qb=0`` makes thigh
  aligned with trunk and ``qb>0`` abducts the hip
- ``qc``: angle of right shank relative to right thigh (knee), ``qc=0`` makes
  shank aligned with thigh and ``qc>0`` extends the knee
- ``qd``: angle of right foot relative to right shank (ankle), ``qd=0`` makes
  foot 90 deg to shank and ``qd>0`` dorsiflexes the foot
- ``qe``: angle of left thigh relative to trunk (hip), ``qe=0`` makes thigh
  aligned with trunk and ``qe>0`` abducts the hip
- ``qf``: angle of left shank relative to left thigh (knee), ``qf=0`` makes
  shank aligned with thigh and ``qf>0`` extends the knee
- ``qg``: angle of left foot relative to left shank (ankle), ``qg=0`` makes
  foot 90 deg to shank and ``qg>0`` dorsiflexes the foot

Specified inputs
----------------

- ``Fax, Fay``: "hand of god", forces acting on the trunk mass center relative
  to inertial origin
- ``Ta``: "hand of god", torque acting on trunk relative to inertial frame
- ``Tb``: hip joint torque, ``Tb>0`` extends the hip
- ``Tc``: knee joint torque, ``Tc>0`` abducts the knee
- ``Td``: ankle joint torque, ``Td>0`` plantarflexes the foot
- ``Te``: hip joint torque, ``Te>0`` extends the hip
- ``Tf``: knee joint torque, ``Tf>0`` abducts the knee
- ``Tg``: ankle joint torque, ``Tg>0`` plantarflexes the foot

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
