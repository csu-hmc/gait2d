#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import symbols
import sympy.physics.mechanics as me

# internal libraries
from segment import BodySegment, TrunkSegment, FootSegment, contact_force


def derive_equations_of_motion(trig_simp=False):
    """Returns the full mass matrix and forcing vector for the walking model
    along with all of the constants, coordinates, speeds, and joint torques.

    Parameters
    ==========
    trig_simp : boolean, optional, default=False
        trigsimp will be applied to each expression in the mass matrix and
        forcing vector. This will slow the derivation down but give smaller
        expressions. TODO: May be smarter to do this on each component that
        buils the EoMs instead of at the end.


    Returns
    ==========
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

    """
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

    for label in sorted(segment_descriptions.keys()):

        segment_class, desc, joint_desc = segment_descriptions[label]

        if label == 'A':  # trunk
            parent_reference_frame = ground
            origin_joint = origin
        elif label == 'E':  # left thigh
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

    # add contact force for trunk mass center.
    external_forces_torques.append((segments[0].mass_center,
                                    contact_force(segments[0].mass_center,
                                                  ground, origin)))
    # add hand of god
    trunk_force_x, trunk_force_y = me.dynamicsymbols('Fax, Fay')
    specified = [trunk_force_x, trunk_force_y] + specified
    external_forces_torques.append((segments[0].mass_center, trunk_force_x *
                                    ground.x + trunk_force_y * ground.y))

    # add contact model constants
    constants += list(symbols('kc, cc, mu, vs'))

    # equations of motion
    print("Initializing Kane's Method.")
    kane = me.KanesMethod(ground, coordinates, speeds, kinematic_equations)
    print("Forming Kane's Equations.")
    kane.kanes_equations(external_forces_torques, bodies)

    if trig_simp is True:
        # If trig_simp is used, which takes a long time, it would be nice to
        # pickle the results. Seems that the standard pickle module may have
        # trouble with that, but the dill package can probably do it.
        # https://pypi.python.org/pypi/dill
        for i, expression in enumerate(kane.mass_matrix_full):
            print("Simplifying matrix expression {}".format(i))
            kane.mass_matrix_full[i] = expression.trigsimp()

        for i, expression in enumerate(kane.forcing_full):
            print("Simplifying forcing expression {}".format(i))
            kane.forcing_full[i] = expression.trigsimp()

    return kane, constants, coordinates, speeds, specified


def analytic_solve(kane):
    """Computes the right hand side of the ODE's symbolically."""
    rhs = kane.mass_matrix.cholesky_solve(kane.forcing)
    return rhs


def compute_jacobians(fr_plus_frstar, coordinates, speeds):

    fr_plus_frstar.jacobian(coordinates)
    fr_plus_frstar.jacobian(speeds)
