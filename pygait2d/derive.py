#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
import sympy as sy
import sympy.physics.mechanics as me

# internal libraries
from .segment import (BodySegment, TrunkSegment, FootSegment, contact_force,
                      time_varying, time_symbol)

me.dynamicsymbols._t = time_symbol


def derive_equations_of_motion(trig_simp=False, seat_force=False,
                               gait_cycle_control=False):
    """Returns the equations of motion for the walking model along with all
    of the constants, coordinates, speeds, joint torques, visualization
    frames, inertial reference frame, and origin point.

    Parameters
    ==========
    trig_simp : boolean, optional, default=False
        sympy.trigsimp() will be applied to each expression in the mass
        matrix and forcing vector. This will slow the derivation down but
        give smaller expressions. TODO: May be smarter to do this on each
        component that builds the EoMs instead of at the end.
    seat_force : boolean, optional, default=False
        If true, a contact force will be added to the hip joint to represent a
        surface higher than the ground to sit on.
    gait_cycle_control : boolean, optinal, default=False
        If true, the specified forces and torques are replaced with a full
        state feeback controller summed with the forces and torques.

    Returns
    =======
    mass_matrix : Matrix
    forcing_vector : Matrix
    kane : sympy.physics.mechanics.Kane
        A Kane object in which the EoMs have been derived.
    constants : list of sympy.core.symbol.Symbol
        The constants in the equations of motion.
    coordinates : list of sympy.core.function.Function
        The generalized coordinates of the system.
    speeds : list of sympy.core.function.Function
        The generalized speeds of the system.
    specified : list of sympy.core.function.Function, optional, default=None
        The specifed quantities of the system.
    visualization_frames : list of VizFrame
    ground : ReferenceFrame
    origin : Point
    segments : list of Segment

    """
    if trig_simp is True:
        me.Vector.simp = True

    print('Forming positions, velocities, accelerations and forces.')
    segment_descriptions = {'A': (TrunkSegment, 'Trunk', 'Hip'),
                            'B': (BodySegment, 'Right Thigh', 'Right Knee'),
                            'C': (BodySegment, 'Right Shank', 'Right Ankle'),
                            'D': (FootSegment, 'Right Foot', 'Right Heel'),
                            'E': (BodySegment, 'Left Thigh', 'Left Knee'),
                            'F': (BodySegment, 'Left Shank', 'Left Ankle'),
                            'G': (FootSegment, 'Left Foot', 'Left Heel')}

    ground = me.ReferenceFrame('N')
    origin = me.Point('O')
    origin.set_vel(ground, 0)

    segments = []
    constants = []
    coordinates = []
    speeds = []
    specified = []
    kinematic_equations = []
    external_forces_torques = []
    bodies = []
    visualization_frames = []

    for label in sorted(segment_descriptions.keys()):

        segment_class, desc, joint_desc = segment_descriptions[label]

        if label == 'A':  # trunk
            parent_reference_frame = ground
            origin_joint = origin
        elif label == 'E':  # left thigh
            # For the left thigh, set the trunk and hip as the
            # reference_frame and origin joint.
            parent_reference_frame = segments[0].reference_frame
            origin_joint = segments[0].joint
        else:  # thighs, shanks
            parent_reference_frame = segments[-1].reference_frame
            origin_joint = segments[-1].joint

        segment = segment_class(label, desc, parent_reference_frame,
                                origin_joint, joint_desc, ground)
        segments.append(segment)

        # constants, coordinates, speeds, kinematic differential equations
        if label == 'A':  # trunk
            coordinates += segment.qa
            speeds += segment.ua
            constants += segment.constants
        else:
            # skip g for all segments but the trunk
            constants += segment.constants[1:]

        coordinates.append(segment.generalized_coordinate_symbol)
        speeds.append(segment.generalized_speed_symbol)

        kinematic_equations += segment.kinematic_equations

        # gravity
        external_forces_torques.append((segment.mass_center,
                                        segment.gravity))

        # joint torques
        external_forces_torques.append((segment.reference_frame,
                                        segment.torque))
        external_forces_torques.append((segment.parent_reference_frame,
                                        -segment.torque))
        specified.append(segment.joint_torque_symbol)

        # contact force
        if label == 'D' or label == 'G':  # foot
            external_forces_torques.append((segment.heel,
                                            contact_force(segment.heel,
                                                          ground, origin)))
            external_forces_torques.append((segment.toe,
                                            contact_force(segment.toe,
                                                          ground, origin)))
        else:
            external_forces_torques.append((segment.joint,
                                            contact_force(segment.joint,
                                                          ground, origin)))

        # bodies
        bodies.append(segment.rigid_body)

        visualization_frames += segment.visualization_frames()

    if seat_force:
        seat_level = origin.locatenew(
            'seat', (segments[2].length_symbol -
                     segments[3].foot_depth)*ground.y)
        external_forces_torques.append((segments[0].joint,
                                        contact_force(segments[0].joint,
                                                      ground, seat_level)))

    # add contact force for trunk mass center.
    external_forces_torques.append((segments[0].mass_center,
                                    contact_force(segments[0].mass_center,
                                                  ground, origin)))
    # add hand of god
    # TODO : move this into segment.py
    trunk_force_x, trunk_force_y = time_varying('Fax, Fay')
    specified = [trunk_force_x, trunk_force_y] + specified
    external_forces_torques.append((segments[0].mass_center, trunk_force_x *
                                    ground.x + trunk_force_y * ground.y))

    # add contact model constants
    # TODO : these should be grabbed from the segments, not recreated.
    constants += list(sy.symbols('kc, cc, mu, vs', real=True, positive=True))

    # equations of motion
    print("Initializing Kane's Method.")
    kane = me.KanesMethod(ground, coordinates, speeds, kinematic_equations)
    print("Forming Kane's Equations.")
    kane.kanes_equations(bodies, loads=external_forces_torques)
    mass_matrix = kane.mass_matrix_full
    forcing_vector = kane.forcing_full

    if trig_simp is True:
        # If trig_simp is used, which takes a long time, it would be nice to
        # pickle the results. Seems that the standard pickle module may have
        # trouble with that, but the dill package can probably do it.
        # https://pypi.python.org/pypi/dill
        # TODO : This should be done in parallel.
        # TODO : Maybe I should enable Vector.simp == True so this happens
        # as things go along instead of all at the end.
        # TODO : Simplifying the mass matrix doesn't take too long, but the
        # forcing vector takes really long.
        for i, expression in enumerate(kane.mass_matrix_full):
            print("Simplifying matrix expression {}".format(i))
            kane.mass_matrix_full[i] = expression.trigsimp()

        for i, expression in enumerate(kane.forcing_full):
            print("Simplifying forcing expression {}".format(i))
            kane.forcing_full[i] = expression.trigsimp()

    if gait_cycle_control:
        # joint_torques(phase) = mean_joint_torque + K*(joint_state_desired -
        # joint_state)
        # r = [Fax(t), Fay(t), Ta(t), Tb(t), Tc(t), Td(t), Te(t), Tf(t), Tg(t)]
        # x = [qax(t), qay(t), qa(t), qb(t), qc(t), qd(t), qe(t), qf(t), qg(t),
        #      uax(t), uay(t), ua(t), ub(t), uc(t), ud(t), ue(t), uf(t), ug(t)]
        # commanded states
        # xc = [qax(t), qay(t), qa(t), qb(t), qc(t), qd(t), qe(t), qf(t), qg(t)]
        #       uax(t), uay(t), ua(t), ub(t), uc(t), ud(t), ue(t), uf(t), ug(t)]
        # controls
        # uc(t) = r(t) + K(t)*(xc(t) - x(t))
        # K is, in general, 9 x 18
        # the first three rows and columns will be zero if hand of god is
        # absent, which effectively makes it a 6x6
        # K = |kax_qax, kax_qay, kax_qa, kax_qb, kax_qc, kax_qd, kax_qe, kax_qf, kax_qg,
        #      kax_uax, kax_uay, kax_ua, kax_ub, kax_uc, kax_ud, kax_ue, kax_uf, kax_ug|
        #     |kay_qax, kay_qay, kay_qa, kay_qb, kay_qc, kay_qd, kay_qe, kay_qf, kay_qg|
        #     |ka_qax, ka_qay, ka_qa, ka_qb, ka_qc, ka_qd, ka_qe, ka_qf, ka_qg|
        #     ...
        #     |kg_qax, kg_qay, kg_qa, kg_qb, kg_qc, kg_qd, kg_qe, kg_qf, kg_qg|
        # We can just go through the final equations of motion and replace the
        # joint torques Tb through Tg with Tb -> Tb + kb_qb*(qb_des - qb) +
        # kb_ub*(ub_des - qb) + ...
        K = []
        for ri in specified:
            row = []
            for xi in coordinates + speeds:
                row.append(sy.Function('k_{}_{}'.format(ri.name, xi.name),
                                       real=True)(time_symbol))
            K.append(row)
        K = sy.Matrix(K)

        xc = []
        for xi in coordinates + speeds:
            xc.append(sy.Function('{}_c'.format(xi.name),
                                  real=True)(time_symbol))
        xc = sy.Matrix(xc)

        uc = (sy.Matrix(specified) + K@(xc - sy.Matrix(coordinates + speeds)))

        repl = {k: v for k, v in zip(specified, uc)}

        forcing_vector = forcing_vector.xreplace(repl)

        specified += K[:]
        specified += xc[:]

    return (mass_matrix, forcing_vector, kane, constants, coordinates, speeds,
            specified, visualization_frames, ground, origin, segments)
