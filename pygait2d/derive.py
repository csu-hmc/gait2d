#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
import sympy as sy
import sympy.physics.mechanics as me

# internal libraries
from .segment import (BodySegment, TrunkSegment, FootSegment, contact_force,
                      time_varying, time_symbol)

me.dynamicsymbols._t = time_symbol


def derive_equations_of_motion(trig_simp=False, seat_force=False):
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

    return (mass_matrix, forcing_vector, kane, constants, coordinates, speeds,
            specified, visualization_frames, ground, origin, segments)
