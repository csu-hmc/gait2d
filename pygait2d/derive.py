#!/usr/bin/env python
# -*- coding: utf-8 -*-

# builtin
from dataclasses import dataclass

# external libraries
import sympy as sm
import sympy.physics.mechanics as me
import sympy.physics.biomechanics as bm

# internal libraries
from .segment import (BodySegment, TrunkSegment, FootSegment, contact_force,
                      time_varying, time_symbol)
from .utils import ExtensorPathway

me.dynamicsymbols._t = time_symbol


@dataclass
class Symbolics():
    """Storage for all of the SymPy and SymPy Mechanics objects.

    Parameters
    ==========
    kanes_method: KanesMethod
        A ``KanesMethod`` object in which the equations of motion have been
        derived.
    dyn_diff_eqs: sm.Matrix
        Kane's Fr + Fr* expression.
    constants : list of Symbol
        Constants in the equations of motion.
    specifieds : list of Function(t), optional
        Specifed variables in the equations of motion.
    inertial_frame: me.ReferenceFrame
        An inertial reference frame representing the Earth and the direction of
        the uniform gravitational field.
    origin : Point
        A point fixed in the ground reference frame used for calculating
        translational velocities.
    segments : list of Segment
        All of the segment objects that make up the human.
    viz_frames: list of VisualizationFrame, optional
    muscles: list of MusculotendonDeGroote2016, optional
        All of the musculotendon actuators in the human.
    controller_repl: dictionary

    Attributes
    ==========
    activations: list of Function(t), optional
        Muscle activation state variables.
    coordinates : list of Function(t)
        The generalized coordinates of the system.
    excitations: list of Function(t), optional
        Specified muscle excitation input variables.
    mus_diff_eqs: sm.Matrix = None
    speeds : list of Function(t)
        The generalized speeds of the system.
    states : list of Function(t)

    """
    # TODO : Add these as properties.
    # forcing_vector : Matrix, shape(18, 1)
    # Full forcing vector where: ``mass_matrix*x' = forcing vector``.
    # mass_matrix : Matrix, shape(18, 18)
    # Full mass matrix of the system to be multiplied by ``x' =
    # [coordinates', speeds']``.

    # TODO : Sort out the order of these for each model combination.
    # states = [coordinates, speeds, activations]
    # specifieds = [forces, torques, excitations]
    # eoms = [kinematical, dynamical, muscle]

    kanes_method: me.KanesMethod
    dyn_diff_eqs: sm.Matrix
    constants: list
    specifieds: list
    inertial_frame: me.ReferenceFrame
    origin: me.Point
    segments: list
    viz_frames: list = None
    muscles: list = None
    controller_repl: dict = None

    def __str__(self):
        states = self.states
        specifieds = self.specifieds
        constants = self.constants
        msg = (
            f'States ({len(states)}): {states}\n\n'
            f'Specifieds ({len(specifieds)}): {specifieds}\n\n'
            f'Constants ({len(constants)}): {constants}'
        )
        return msg

    @property
    def kin_diff_eqs(self):
        return sm.Matrix([k - v for k, v in
                          self.kanes_method.kindiffdict().items()])

    @property
    def mus_diff_eqs(self):
        if self.muscles is not None:
            return sm.Matrix([mus.a.diff() - mus.rhs()[0, 0]
                              for mus in self.muscles])
        else:
            return None

    @property
    def equations_of_motion(self):
        eoms = self.kin_diff_eqs.col_join(self.dyn_diff_eqs)
        if self.muscles is not None:
            eoms = eoms.col_join(self.mus_diff_eqs)
        if self.controller_repl is not None:
            eoms = eoms.xreplace(self.controller_repl)
        return eoms

    @property
    def coordinates(self):
        return self.kanes_method.q[:]

    @property
    def speeds(self):
        return self.kanes_method.u[:]

    @property
    def activations(self):
        if self.muscles is not None:
            return [mus.a for mus in self.muscles]
        else:
            return []

    @property
    def states(self):
        states = self.coordinates + self.speeds
        if self.activations:
            states += self.activations
        return states

    @property
    def excitations(self):
        if self.muscles is not None:
            return [mus.e for mus in self.muscles]
        else:
            return []


def generate_gait_cycle_torque_controller(coordinates, speeds, specified):
    # joint_torques(phase) = mean_joint_torque + K*(joint_state_desired -
    # joint_state)
    # r = [Fax(t), Fay(t), Ta(t), Tb(t), Tc(t), Td(t), Te(t), Tf(t), Tg(t)]
    # x = [qax(t), qay(t), qa(t), qb(t), qc(t), qd(t), qe(t), qf(t), qg(t),
    #      uax(t), uay(t), ua(t), ub(t), uc(t), ud(t), ue(t), uf(t), ug(t)]
    # commanded states
    # xc = [qax_c(t), qay_c(t), qa_c(t), qb_c(t), qc_c(t), qd_c(t), qe_c(t),
    #       qf_c(t), qg_c(t),
    #       uax_c(t), uay_c(t), ua_c(t), ub_c(t), uc_c(t), ud_c(t), ue_c(t),
    #       uf_c(t), ug_c(t)]
    # controlled joint torques
    # uc(t) = r(t) + K(t)*(xc(t) - x(t))
    # r(t) : force or torque
    # K(t) : time varying full state feedback gain matrix
    # xc(t) : commanded (desired) states
    # x(t) : states
    # K is, in general, 9 x 18
    # the first three rows and columns will be zero if hand of god is
    # absent, which effectively makes it a 6x18
    # K = |kax_qax, kax_qay, kax_qa, kax_qb, kax_qc, kax_qd, kax_qe, kax_qf,
    #      kax_qg, kax_uax, kax_uay, kax_ua, kax_ub, kax_uc, kax_ud, kax_ue,
    #      kax_uf, kax_ug|
    #     |kay_qax, kay_qay, kay_qa, kay_qb, kay_qc, kay_qd, kay_qe, kay_qf,
    #      kay_qg|
    #     |ka_qax, ka_qay, ka_qa, ka_qb, ka_qc, ka_qd, ka_qe, ka_qf, ka_qg|
    #     ...
    #     |kg_qax, kg_qay, kg_qa, kg_qb, kg_qc, kg_qd, kg_qe, kg_qf, kg_qg|
    # We can just go through the final equations of motion and replace the
    # joint torques Tb through Tg with Tb -> Tb + kb_qb*(qb_des - qb) +
    # kb_ub*(ub_des - qb) + ...
    print('Generating gait cycle torque controller.')
    K = []
    for ri in specified:
        row = []
        for xi in coordinates + speeds:
            row.append(time_varying('k_{}_{}'.format(ri.name, xi.name)))
        K.append(row)
    K = sm.Matrix(K)

    xc = []
    for xi in coordinates + speeds:
        xc.append(time_varying('{}_c'.format(xi.name))),

    r = sm.Matrix(specified)
    xc = sm.Matrix(xc)
    x = sm.Matrix(coordinates + speeds)

    uc = r + K@(xc - x)

    repl = {k: v for k, v in zip(specified, uc)}

    specified += K[:]
    specified += xc[:]

    return repl, specified


def generate_muscles(segments):
    """Returns the loads due to the musculotendon actuators and the activation
    dynamics differential equations."""

    print('Generating musculotendon pathways and activation dynamics.')

    # The Pathway type followed by the origin, (middle,) inersetion bodies
    muscle_descriptions = {
        'ilio_r': ('Linear', 'A', 'B'),
        'hams_r': ('Linear', 'A', 'C'),
        'glut_r': ('Obstacle', 'A', 'A', 'B'),
        'rect_r': ('Extensor', 'A', 'B', 'C'),
        'vast_r': ('Extensor', 'B', 'B', 'C'),
        'gast_r': ('Obstacle', 'B', 'C', 'D'),
        'sole_r': ('Linear', 'C', 'D'),
        'tibi_r': ('Obstacle', 'C', 'C', 'D'),
        'ilio_l': ('Linear', 'A', 'E'),
        'hams_l': ('Linear', 'A', 'F'),
        'glut_l': ('Obstacle', 'A', 'A', 'E'),
        'rect_l': ('Extensor', 'A', 'E', 'F'),
        'vast_l': ('Extensor', 'E', 'E', 'F'),
        'gast_l': ('Obstacle', 'E', 'F', 'G'),
        'sole_l': ('Linear', 'F', 'G'),
        'tibi_l': ('Obstacle', 'F', 'F', 'G'),
    }

    def get_segment_by_label(label):
        """Returns Segment based on label A, B, C, D, E, F, G."""
        for seg in segments:
            if seg.reference_frame.name == label:
                return seg
        raise ValueError(f'No segment with label: {label}!')

    def setup_point(muscle_label, body_label, point_name):
        """Creates and attaches a point to a segment relative to its reference
        point."""
        label = '_'.join([muscle_label, body_label, point_name])
        point = me.Point(label)
        # point will be fixed on this segment:
        seg = get_segment_by_label(body_label)
        x, y, r = sm.symbols(label + '_x, ' + label + '_y, ' + label + '_r',
                             real=True)
        # the muscle points are defined using the body fixed unit vectors for
        # the body that the point is fixed in
        if body_label == 'A':
            origin_point = seg.joint
        else:
            origin_point = seg.origin_joint
        point.set_pos(origin_point,
                      x*seg.reference_frame.x + y*seg.reference_frame.y)
        point.v2pt_theory(origin_point,
                          seg.inertial_frame,
                          seg.reference_frame)
        return point, x, y, r

    muscles = []
    mus_loads = []
    mus_excit = []
    mus_const = []
    for mus_label, pathway_data in muscle_descriptions.items():
        pathway_type = pathway_data[0]
        body_labels = pathway_data[1:]

        origin_point, x, y, _ = setup_point(mus_label, body_labels[0],
                                            'origin')
        mus_const += [x, y]
        insert_point, x, y, _ = setup_point(mus_label, body_labels[-1],
                                            'insert')
        mus_const += [x, y]

        if pathway_type == 'Linear':
            pathway = me.LinearPathway(origin_point, insert_point)
        elif len(body_labels) > 2:
            middle_point, x, y, r = setup_point(mus_label, body_labels[1],
                                                'middle')
            if pathway_type == 'Obstacle':
                mus_const += [x, y]  # skip radius r
                pathway = me.ObstacleSetPathway(origin_point, middle_point,
                                                insert_point)
            elif pathway_type == 'Extensor':
                mus_const += [x, y, r]
                seg = get_segment_by_label(body_labels[-1])  # shin
                pathway = ExtensorPathway(
                    origin_point,
                    insert_point,
                    middle_point,
                    seg.inertial_frame.z,
                    origin_point.pos_from(middle_point),
                    insert_point.pos_from(middle_point),
                    seg.parent_reference_frame.y,
                    seg.parent_reference_frame.x,
                    seg.reference_frame.y,
                    seg.reference_frame.x,
                    r,  # radius
                    # a negative knee angle flexes the knee (so switch sign)
                    -seg.generalized_coordinate_symbol)

        act = bm.FirstOrderActivationDeGroote2016.with_defaults(mus_label)
        mus = bm.MusculotendonDeGroote2016.with_defaults(
            mus_label, pathway, act)
        muscles.append(mus)
        mus_loads += list(mus.to_loads())
        mus_excit.append(mus.e)
        mus_const += mus.constants[:]

    return mus_loads, mus_excit, mus_const, muscles


def derive_equations_of_motion(seat_force=False, gait_cycle_control=False,
                               include_muscles=False,
                               prevent_ground_penetration=False):
    """Returns the equations of motion for the planar walking model along with
    all of the constants, coordinates, speeds, joint torques, visualization
    frames, inertial reference frame, and origin point.

    Parameters
    ==========
    seat_force : boolean, optional, default=False
        If true, a contact force will be added to the hip joint to represent a
        surface higher than the ground to sit on.
    gait_cycle_control : boolean, optional, default=False
        If true, the specified forces and torques are replaced with a full
        state feeback controller summed with the forces and torques.
    include_muscles : boolean, optional, default=False
        If true, muscle actuators will be included in addition to the joint
        torque actuators.
    prevent_ground_penetration : boolean, optional
        If true, the ground force will be added to all joint centers as well as
        the feet to prevent the model from penetrating the ground at all.
        Otherwise, the force will only be applied to the feet bottoms.

    Returns
    =======
    symbolics: Symbolics
        Contains all symbolic model components, see :py:class:`Symbolics`.

    """

    print('Forming positions, velocities, accelerations and forces.')
    # reference frame label: Segment, segment name, distal joint name
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
            if prevent_ground_penetration:
                external_forces_torques.append((segment.joint,
                                                contact_force(segment.joint,
                                                              ground, origin)))

        # bodies
        bodies.append(segment.rigid_body)

        visualization_frames += segment.visualization_frames()

    if seat_force:
        # NOTE : The seat height is set to the length of the shank + the depth
        # of the foot.
        seat_level = origin.locatenew(
            'seat', (segments[2].length_symbol -
                     segments[3].foot_depth)*ground.y)
        external_forces_torques.append((segments[0].joint,
                                        contact_force(segments[0].joint,
                                                      ground, seat_level)))

    if prevent_ground_penetration:
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
    constants += list(sm.symbols('kc, cc, mu, vs', real=True, positive=True))

    if include_muscles:
        (mus_loads, mus_exc, mus_con, muscles) = generate_muscles(segments)
        external_forces_torques += mus_loads
        specified += mus_exc
        constants += mus_con

    # equations of motion
    print("Initializing Kane's Method.")
    kane = me.KanesMethod(ground, coordinates, speeds, kinematic_equations)
    print("Forming Kane's Equations.")
    fr, frstar = kane.kanes_equations(bodies, loads=external_forces_torques)

    if gait_cycle_control:
        repl, specified = generate_gait_cycle_torque_controller(coordinates,
                                                                speeds,
                                                                specified)

    sym_mod = Symbolics(
        kanes_method=kane,
        dyn_diff_eqs=fr + frstar,
        constants=constants,
        specifieds=specified,
        inertial_frame=ground,
        origin=origin,
        segments=segments,
        viz_frames=visualization_frames,
    )

    if include_muscles:
        sym_mod.muscles = muscles

    if gait_cycle_control:
        sym_mod.controller_repl = repl

    return sym_mod
