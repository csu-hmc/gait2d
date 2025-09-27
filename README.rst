Introduction
============

gait2d is a dynamic model that contains the essential elements for simulating
human gait in the sagittal plane. The model has been used in earlier forms by
Gerritsen et al.  (1998) [Gerritsen1998]_ and Ackermann and van den Bogert
(2010) [Ackermann2010]_. The main features of the model are:

- Seven body segments
- Hunt-Crossley style foot contact model
- Sixteen musculotendon actuators
- Differentiable dynamics
- Fast execution

Model dynamics and outputs are twice differentiable with respect to all inputs,
which is important for certain numerical methods for simulation and optimal
control. The model is intended for education and basic research. The model can
be used for studies such as Ackermann and van den Bogert (2010) and Geyer and
Herr (2010) [Geyer2010]_.

The included Autolev model was originally written by Ton van den Bogert and
adapted by Jason K. Moore. The SymPy Mechanics model was created by Jason to
match the Autolev model results for use with PyDy and opty.

Speed
=====

The cythonized Autolev C code takes about 5 micro seconds per rhs eval and the
PyDy Cython version takes about 15 microseconds (the PyDy version is slower
because of the Python level linear solve on the full mass matrix, but only for
forward dynamics).

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
   >>> from algait2de import gait2de
   >>> gait2de.evaluate_autolev_rhs(...)

Model Description
=================

The model is planar (xy plane) and there are seven rigid bodies. All rotations
are defined to be positive about the z axis for a child body relative to its
parent body, e.g thigh relative to trunk. When all configuration variables are
set to zero, the model is standing upright and the trunk's hip joint is at the
inertial origin with the human facing in the positive x direction.

Rigid bodies and constants
--------------------------

- Trunk (A):

  - ``ma``: mass [kg]
  - ``ia``: moment of inertia about mass center wrt to the trunk reference
    frame [kg*m^2]
  - ``xa``: local x location of mass center wrt to the hip joint [m]
  - ``ya``: local y location of mass center wrt to the hip joint [m]
  - ``glut_r_A_origin_x``: local right x location of origin wrt to the hip joint [m]
  - ``glut_r_A_origin_y``: local right y location of origin wrt to the hip joint [m]
  - ``glut_r_A_middle_x``: local right x location of middle wrt to the hip joint [m]
  - ``glut_r_A_middle_y``: local right y location of middle wrt to the hip joint [m]
  - ``hams_r_A_origin_x``: local right x location of origin wrt to the hip joint [m]
  - ``hams_r_A_origin_y``: local right y location of origin wrt to the hip joint [m]
  - ``ilio_r_A_origin_x``: local right x location of origin wrt to the hip joint [m]
  - ``ilio_r_A_origin_y``: local right y location of origin wrt to the hip joint [m]
  - ``rect_r_A_origin_x``: local right x location of origin wrt to the hip joint [m]
  - ``rect_r_A_origin_y``: local right y location of origin wrt to the hip joint [m]
  - ``glut_l_A_origin_x``: local left x location of origin wrt to the hip joint [m]
  - ``glut_l_A_origin_y``: local left y location of origin wrt to the hip joint [m]
  - ``glut_l_A_middle_x``: local left x location of middle wrt to the hip joint [m]
  - ``glut_l_A_middle_y``: local left y location of middle wrt to the hip joint [m]
  - ``hams_l_A_origin_x``: local left x location of origin wrt to the hip joint [m]
  - ``hams_l_A_origin_y``: local left y location of origin wrt to the hip joint [m]
  - ``ilio_l_A_origin_x``: local left x location of origin wrt to the hip joint [m]
  - ``ilio_l_A_origin_y``: local left y location of origin wrt to the hip joint [m]
  - ``rect_l_A_origin_x``: local left x location of origin wrt to the hip joint [m]
  - ``rect_l_A_origin_y``: local left y location of origin wrt to the hip joint [m]

- Right Thigh (B)

  - ``mb``: mass [kg]
  - ``ib``: moment of inertia about mass center wrt to right thigh reference
    frame [kg*m^2]
  - ``xb``: local x location of mass center wrt to the hip joint [m]
  - ``yb``: local y location of mass center wrt to the hip joint [m]
  - ``lb``: joint to joint segment length [m]
  - ``gast_r_B_origin_x``: local x location of origin wrt to the hip joint [m]
  - ``gast_r_B_origin_y``: local y location of origin wrt to the hip joint [m]
  - ``glut_r_B_insert_x``: local x location of insertion wrt to the hip joint [m]
  - ``glut_r_B_insert_y``: local y location of insertion wrt to the hip joint [m]
  - ``ilio_r_B_insert_x``: local x location of insertion wrt to the hip joint [m]
  - ``ilio_r_B_insert_y``: local y location of insertion wrt to the hip joint [m]
  - ``rect_r_B_middle_x``: local x location of middle wrt to the hip joint [m]
  - ``rect_r_B_middle_y``: local y location of middle wrt to the hip joint [m]
  - ``rect_r_B_middle_r``: radius of knee wrapping circle [m]
  - ``vast_r_B_middle_x``: local x location of middle wrt to the hip joint [m]
  - ``vast_r_B_middle_y``: local y location of middle wrt to the hip joint [m]
  - ``vast_r_B_middle_r``: radius of knee wrapping circle [m]
  - ``vast_r_B_origin_x``: local x location of origin wrt to the hip joint [m]
  - ``vast_r_B_origin_y``: local y location of origin wrt to the hip joint [m]

- Right Shank (C)

  - ``mc``: mass [kg]
  - ``ic``: moment of inertia about mass center wrt to right shank reference
    frame [kg*m^2]
  - ``xc``: x location of mass center wrt to the knee joint [m]
  - ``yc``: y location of mass center wrt to the knee joint [m]
  - ``lc``: joint to joint segment length [m]
  - ``gast_r_C_middle_x``: local x location of middle wrt to the knee joint [m]
  - ``gast_r_C_middle_y``: local y location of middle wrt to the knee joint [m]
  - ``hams_r_C_insert_x``: local x location of insertion wrt to the knee joint [m]
  - ``hams_r_C_insert_y``: local y location of insertion wrt to the knee joint [m]
  - ``rect_r_C_insert_x``: local x location of insertion wrt to the knee joint [m]
  - ``rect_r_C_insert_y``: local y location of insertion wrt to the knee joint [m]
  - ``sole_r_C_origin_x``: local x location of origin wrt to the knee joint [m]
  - ``sole_r_C_origin_y``: local y location of origin wrt to the knee joint [m]
  - ``tibi_r_C_origin_x``: local x location of origin wrt to the knee joint [m]
  - ``tibi_r_C_origin_y``: local y location of origin wrt to the knee joint [m]
  - ``tibi_r_C_middle_x``: local x location of middle wrt to the knee joint [m]
  - ``tibi_r_C_middle_y``: local y location of middle wrt to the knee joint [m]
  - ``vast_r_C_insert_x``: local x location of insertion wrt to the knee joint [m]
  - ``vast_r_C_insert_y``: local y location of insertion wrt to the knee joint [m]

- Right Foot (D)

  - ``md``: mass [kg]
  - ``id``: moment of inertia about mass center wrt to right foot reference
    frame [kg*m^2]
  - ``xd``: local x location of mass center wrt to the ankle joint [m]
  - ``yd``: local y location of mass center wrt to the ankle joint [m]
  - ``hxd``: local x location of heel wrt to the ankle joint [m]
  - ``txd``: local x location of toe wrt to the ankle joint [m]
  - ``fyd``: local y location of heel and toe relative to ankle joint [m]
  - ``gast_r_D_insert_x``: local x location of insertion wrt to the ankle joint [m]
  - ``gast_r_D_insert_y``: local y location of insertion wrt to the ankle joint [m]
  - ``sole_r_D_insert_x``: local x location of insertion wrt to the ankle joint [m]
  - ``sole_r_D_insert_y``: local y location of insertion wrt to the ankle joint [m]
  - ``tibi_r_D_insert_x``: local x location of insertion wrt to the ankle joint [m]
  - ``tibi_r_D_insert_y``: local y location of insertion wrt to the ankle joint [m]

- Left Thigh (E)

  - ``me``: mass [kg]
  - ``ie``: moment of inertia about mass center wrt to left thigh reference
    frame [kg*m^2]
  - ``xe``: local x location of mass center wrt to the hip joint [m]
  - ``ye``: local y location of mass center wrt to the hip joint [m]
  - ``le``: joint to joint segment length [m]
  - ``gast_l_E_origin_x``: local x location of origin wrt to the hip joint [m]
  - ``gast_l_E_origin_y``: local y location of origin wrt to the hip joint [m]
  - ``glut_l_E_insert_x``: local x location of insertion wrt to the hip joint [m]
  - ``glut_l_E_insert_y``: local y location of insertion wrt to the hip joint [m]
  - ``ilio_l_E_insert_x``: local x location of insertion wrt to the hip joint [m]
  - ``ilio_l_E_insert_y``: local y location of insertion wrt to the hip joint [m]
  - ``rect_l_E_middle_x``: local x location of middle wrt to the hip joint [m]
  - ``rect_l_E_middle_y``: local y location of middle wrt to the hip joint [m]
  - ``rect_l_E_middle_r``: radius of knee wrapping circle [m]
  - ``vast_l_E_middle_x``: local x location of middle wrt to the hip joint [m]
  - ``vast_l_E_middle_y``: local y location of middle wrt to the hip joint [m]
  - ``vast_l_E_middle_r``: radius of knee wrapping circle [m]
  - ``vast_l_E_origin_x``: local x location of origin wrt to the hip joint [m]
  - ``vast_l_E_origin_y``: local y location of origin wrt to the hip joint [m]

- Left Shank (F)

  - ``mf``: mass [kg]
  - ``if``: moment of inertia about mass center wrt to left shank reference
    frame [kg*m^2]
  - ``xf``: x location of mass center wrt to the knee joint [m]
  - ``yf``: y location of mass center wrt to the knee joint [m]
  - ``lf``: joint to joint segment length [m]
  - ``gast_l_F_middle_x``: local x location of middle wrt to the knee joint [m]
  - ``gast_l_F_middle_y``: local y location of middle wrt to the knee joint [m]
  - ``hams_l_F_insert_x``: local x location of insertion wrt to the knee joint [m]
  - ``hams_l_F_insert_y``: local y location of insertion wrt to the knee joint [m]
  - ``rect_l_F_insert_x``: local x location of insertion wrt to the knee joint [m]
  - ``rect_l_F_insert_y``: local y location of insertion wrt to the knee joint [m]
  - ``sole_l_F_origin_x``: local x location of origin wrt to the knee joint [m]
  - ``sole_l_F_origin_y``: local y location of origin wrt to the knee joint [m]
  - ``tibi_l_F_origin_x``: local x location of origin wrt to the knee joint [m]
  - ``tibi_l_F_origin_y``: local y location of origin wrt to the knee joint [m]
  - ``tibi_l_F_middle_x``: local x location of middle wrt to the knee joint [m]
  - ``tibi_l_F_middle_y``: local y location of middle wrt to the knee joint [m]
  - ``vast_l_F_insert_x``: local x location of insertion wrt to the knee joint [m]
  - ``vast_l_F_insert_y``: local y location of insertion wrt to the knee joint [m]

- Left Foot (G)

  - ``mg``: mass [kg]
  - ``ig``: moment of inertia about mass center wrt to left foot reference
    frame [kg*m^2]
  - ``xg``: local x location of mass center wrt to the ankle joint [m]
  - ``yg``: local y location of mass center wrt to the ankle joint [m]
  - ``hxg``: local x location of heel wrt to the ankle joint [m]
  - ``txg``: local x location of toe wrt to the ankle joint [m]
  - ``fyg``: local y location of heel and toe relative to ankle joint [m]
  - ``gast_l_G_insert_x``: local x location of insertion wrt to the ankle joint [m]
  - ``gast_l_G_insert_y``: local y location of insertion wrt to the ankle joint [m]
  - ``sole_l_G_insert_x``: local x location of insertion wrt to the ankle joint [m]
  - ``sole_l_G_insert_y``: local y location of insertion wrt to the ankle joint [m]
  - ``tibi_l_G_insert_x``: local x location of insertion wrt to the ankle joint [m]
  - ``tibi_l_G_insert_y``: local y location of insertion wrt to the ankle joint [m]

- Muscle Properties

  - Gastrocnemius

     - ``F_M_max_gast_r``: maximum isometric force [N]
     - ``l_M_opt_gast_r``: optimal fiber length [m]
     - ``l_T_slack_gast_r``: tendon slack length [m]
     - ``F_M_max_gast_l``: maximum isometric force [N]
     - ``l_M_opt_gast_l``: optimal fiber length [m]
     - ``l_T_slack_gast_l``: tendon slack length [m]

  - Gluteus maximus, medius, and minimus

     - ``F_M_max_glut_r``: maximum isometric force [N]
     - ``l_M_opt_glut_r``: optimal fiber length [m]
     - ``l_T_slack_glut_r``: tendon slack length [m]
     - ``F_M_max_glut_l``: maximum isometric force [N]
     - ``l_M_opt_glut_l``: optimal fiber length [m]
     - ``l_T_slack_glut_l``: tendon slack length [m]

  - Hamstrings

     - ``F_M_max_hams_r``: maximum isometric force [N]
     - ``l_M_opt_hams_r``: optimal fiber length [m]
     - ``l_T_slack_hams_r``: tendon slack length [m]
     - ``F_M_max_hams_l``: maximum isometric force [N]
     - ``l_M_opt_hams_l``: optimal fiber length [m]
     - ``l_T_slack_hams_l``: tendon slack length [m]

  - Ilioposas

     - ``F_M_max_ilio_r``: maximum isometric force [N]
     - ``l_M_opt_ilio_r``: optimal fiber length [m]
     - ``l_T_slack_ilio_r``: tendon slack length [m]
     - ``F_M_max_ilio_l``: maximum isometric force [N]
     - ``l_M_opt_ilio_l``: optimal fiber length [m]
     - ``l_T_slack_ilio_l``: tendon slack length [m]

  - Rectus femoris

     - ``F_M_max_rect_r``: maximum isometric force [N]
     - ``l_M_opt_rect_r``: optimal fiber length [m]
     - ``l_T_slack_rect_r``: tendon slack length [m]
     - ``F_M_max_rect_l``: maximum isometric force [N]
     - ``l_M_opt_rect_l``: optimal fiber length [m]
     - ``l_T_slack_rect_l``: tendon slack length [m]

  - Soleus

     - ``F_M_max_sole_r``: maximum isometric force [N]
     - ``l_M_opt_sole_r``: optimal fiber length [m]
     - ``l_T_slack_sole_r``: tendon slack length [m]
     - ``F_M_max_sole_l``: maximum isometric force [N]
     - ``l_M_opt_sole_l``: optimal fiber length [m]
     - ``l_T_slack_sole_l``: tendon slack length [m]

  - Tibialis anterior

     - ``F_M_max_tibi_r``: maximum isometric force [N]
     - ``l_M_opt_tibi_r``: optimal fiber length [m]
     - ``l_T_slack_tibi_r``: tendon slack length [m]
     - ``F_M_max_tibi_l``: maximum isometric force [N]
     - ``l_M_opt_tibi_l``: optimal fiber length [m]
     - ``l_T_slack_tibi_l``: tendon slack length [m]

  - Vastus intermedius, medialis, and lateralis

     - ``F_M_max_vast_r``: maximum isometric force [N]
     - ``l_M_opt_vast_r``: optimal fiber length [m]
     - ``l_T_slack_vast_r``: tendon slack length [m]
     - ``F_M_max_vast_l``: maximum isometric force [N]
     - ``l_M_opt_vast_l``: optimal fiber length [m]
     - ``l_T_slack_vast_l``: tendon slack length [m]

- Other constants

  - ``kc``: ground contact stiffness [N/m^3]
  - ``cc``: ground contact damping [s/m]
  - ``mu``: ground contact friction coefficient
  - ``vs``: ground contact velocity constant [m/s]
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

Muscle activation states
------------------------

- ``a_ilio_r``: muscle activation
- ``a_hams_r``: muscle activation
- ``a_glut_r``: muscle activation
- ``a_rect_r``: muscle activation
- ``a_vast_r``: muscle activation
- ``a_gast_r``: muscle activation
- ``a_sole_r``: muscle activation
- ``a_tibi_r``: muscle activation
- ``a_ilio_l``: muscle activation
- ``a_hams_l``: muscle activation
- ``a_glut_l``: muscle activation
- ``a_rect_l``: muscle activation
- ``a_vast_l``: muscle activation
- ``a_gast_l``: muscle activation
- ``a_sole_l``: muscle activation
- ``a_tibi_l``: muscle activation

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
- ``e_ilio_r``: muscle excitation
- ``e_hams_r``: muscle excitation
- ``e_glut_r``: muscle excitation
- ``e_rect_r``: muscle excitation
- ``e_vast_r``: muscle excitation
- ``e_gast_r``: muscle excitation
- ``e_sole_r``: muscle excitation
- ``e_tibi_r``: muscle excitation
- ``e_ilio_l``: muscle excitation
- ``e_hams_l``: muscle excitation
- ``e_glut_l``: muscle excitation
- ``e_rect_l``: muscle excitation
- ``e_vast_l``: muscle excitation
- ``e_gast_l``: muscle excitation
- ``e_sole_l``: muscle excitation
- ``e_tibi_l``: muscle excitation

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
