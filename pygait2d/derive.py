#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
class SymbolicModel():
    # states = [coordinates, speeds, activations]
    # specifieds = [forces, torques, excitations]
    # eoms = [kinematical, dynamical, muscle]

    kanes_method: me.KanesMethod
    dyn_diff_eqs: sm.Matrix
    constants: sm.Matrix
    specifieds: sm.Matrix
    inertial_frame: me.ReferenceFrame
    origin: me.Point
    segments: list
    viz_frames: list
    mus_diff_eqs: sm.Matrix = None
    activations: sm.Matrix = None
    excitations: sm.Matrix = None

    @property
    def equations_of_motion(self):
        eoms = self.kin_diff_eqs.col_join(self.dyn_diff_eqs)
        if self.mus_diff_eqs is not None:
            eoms = eoms.col_join(self.mus_diff_eqs)
        return eoms

    @property
    def states(self):
        states = self.coordinates.col_join(self.speeds)
        if self.activations is not None:
            states = states.col_join(self.activations)
        return states

    @property
    def coordinates(self):
        return self.kanes_method.q

    @property
    def speeds(self):
        return self.kanes_method.u

    @property
    def kin_diff_eqs(self):
        return sm.Matrix([k - v for k, v in
                          self.kanes_method.kindiffdict().items()])


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
        for seg in segments:
            if seg.reference_frame.name == label:
                return seg
        raise ValueError(f'No segment with label: {label}!')

    def setup_point(muscle_label, body_label, point_name):
        label = '_'.join([muscle_label, body_label, point_name])
        point = me.Point(label)
        # point will be fixed on this segment:
        seg = get_segment_by_label(body_label)
        print(seg)
        x, y = sm.symbols(label + '_x, ' + label + '_y', real=True)
        # the muscle points are defined using the body fixed unit vectors for
        # the body that the point is fixed in
        # TODO: Fix the y numerical values for points on thigh and shank
        # The numerical values for muscle attachments were measured from:
        # trunk: hip joint = segment.origin_joint/joint -> correct numbers
        # thigh: knee joint = segment.joint -> incorrect (should be hip joint)
        # shank: ankle joint = segment.joint -> incorrect (should be knee joint)
        # foot: ankle joint = segment.origin_joint/joint -> correct numbers
        point.set_pos(seg.origin_joint,
                      x*seg.reference_frame.x + y*seg.reference_frame.y)
        point.v2pt_theory(seg.origin_joint,
                          seg.inertial_frame,
                          seg.reference_frame)
        return point

    muscle_loads = []
    muscle_actvs = []
    muscle_states = []
    muscle_excit = []
    for muscle_label, pathway_data in muscle_descriptions.items():
        pathway_type = pathway_data[0]
        body_labels = pathway_data[1:]

        origin_point = setup_point(muscle_label, body_labels[0], 'origin')
        insert_point = setup_point(muscle_label, body_labels[-1], 'insert')

        if pathway_type == 'Linear':
            pathway = me.LinearPathway(origin_point, insert_point)
        elif len(body_labels) > 2:
            middle_point = setup_point(muscle_label, body_labels[1], 'middle')
            if pathway_type == 'Obstacle':
                pathway = me.ObstacleSetPathway(origin_point, middle_point,
                                                insert_point)
            elif pathway_type == 'Extensor':
                seg = get_segment_by_label(body_labels[-1])  # shin
                # TODO : Check the definition of the knee angle for this
                # pathway.
                pathway = ExtensorPathway(
                    origin_point,
                    insert_point,
                    middle_point,
                    seg.inertial_frame.z,
                    origin_point.pos_from(middle_point),
                    insert_point.pos_from(middle_point),
                    # TODO : move radius to data file with values
                    0.03,  # radius
                    # a negative knee angle flexes the knee (so switch sign)
                    -seg.generalized_coordinate_symbol)

        act = bm.FirstOrderActivationDeGroote2016.with_defaults(muscle_label)
        mus = bm.MusculotendonDeGroote2016.with_defaults(
            muscle_label, pathway, act)
        muscle_loads += list(mus.to_loads())
        muscle_states.append(mus.a)
        muscle_excit.append(mus.e)
        muscle_actvs.append(mus.a.diff() - mus.rhs()[0, 0])

    return muscle_loads, muscle_actvs, muscle_states, muscle_excit


def derive_equations_of_motion(seat_force=False, gait_cycle_control=False,
                               include_muscles=False):
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

    Returns
    =======
    mass_matrix : Matrix, shape(18, 18)
        Full mass matrix of the system to be multiplied by ``x' =
        [coordinates', speeds']``.
    forcing_vector : Matrix, shape(18, 1)
        Full forcing vector where: ``mass_matrix*x' = forcing vector``.
    kane : sympy.physics.mechanics.Kane
        A KanesMethod object in which the equations of motion have been
        derived. All symbolics are accessible from this object if needed.
    constants : list of Symbol
        The constants in the equations of motion.
    coordinates : list of Function(t)
        The generalized coordinates of the system.
    speeds : list of Function(t)
        The generalized speeds of the system.
    specified : list of Function(t), optional, default=None
        The specifed quantities of the system.
    visualization_frames : list of VizFrame
    ground : ReferenceFrame
        An inertial reference frame representing the Earth and a the direction
        of the uniform gravitational field.
    origin : Point
        A point fixed in the ground reference frame used for calculating
        translational velocities.
    segments : list of Segment
        All of the segment objects that make up the human.

    """

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

    states = coordinates + speeds
    if include_muscles:
        mus_loads, mus_actvs, mus_states, mus_exc = generate_muscles(segments)
        external_forces_torques += mus_loads
        states += mus_states
        specified += mus_exc

    # add contact model constants
    # TODO : these should be grabbed from the segments, not recreated.
    constants += list(sm.symbols('kc, cc, mu, vs', real=True, positive=True))

    # equations of motion
    print("Initializing Kane's Method.")
    kane = me.KanesMethod(ground, coordinates, speeds, kinematic_equations)
    print("Forming Kane's Equations.")
    fr, frstar = kane.kanes_equations(bodies, loads=external_forces_torques)
    mass_matrix = kane.mass_matrix_full
    forcing_vector = kane.forcing_full

    kin_diff_eqs = sm.Matrix([k - v for k, v in kane.kindiffdict().items()])
    equations_of_motion = kin_diff_eqs.col_join(fr + frstar)

    if include_muscles:
        equations_of_motion = equations_of_motion.col_join(
            sm.Matrix(mus_actvs))

    if gait_cycle_control:
        # joint_torques(phase) = mean_joint_torque + K*(joint_state_desired -
        # joint_state)
        # r = [Fax(t), Fay(t), Ta(t), Tb(t), Tc(t), Td(t), Te(t), Tf(t), Tg(t)]
        # x = [qax(t), qay(t), qa(t), qb(t), qc(t), qd(t), qe(t), qf(t), qg(t),
        #      uax(t), uay(t), ua(t), ub(t), uc(t), ud(t), ue(t), uf(t), ug(t)]
        # commanded states
        # xc = [qax_c(t), qay_c(t), qa_c(t), qb_c(t), qc_c(t), qd_c(t), qe_c(t), qf_c(t), qg_c(t)]
        #       uax_c(t), uay_c(t), ua_c(t), ub_c(t), uc_c(t), ud_c(t), ue_c(t), uf_c(t), ug_c(t)]
        # controlled joint torques
        # uc(t) = r(t) + K(t)*(xc(t) - x(t))
        # r(t) : force or torque
        # K(t) : time varying full state feedback gain matrix
        # xc(t) : commanded (desired) states
        # x(t) : states
        # K is, in general, 9 x 18
        # the first three rows and columns will be zero if hand of god is
        # absent, which effectively makes it a 6x18
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
                row.append(sm.Function('k_{}_{}'.format(ri.name, xi.name),
                                       real=True)(time_symbol))
            K.append(row)
        K = sm.Matrix(K)

        xc = []
        for xi in coordinates + speeds:
            xc.append(sm.Function('{}_c'.format(xi.name),
                                  real=True)(time_symbol))
        r = sm.Matrix(specified)
        xc = sm.Matrix(xc)
        x = sm.Matrix(coordinates + speeds)

        uc = r + K@(xc - x)

        repl = {k: v for k, v in zip(specified, uc)}

        forcing_vector = forcing_vector.xreplace(repl)

        specified += K[:]
        specified += xc[:]

    sym_mod = SymbolicModel(
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
        sym_mod.mus_diff_eqs = sm.Matrix(mus_actvs)
        sym_mod.activations = sm.Matrix(mus_states)
        sym_mod.excitations = sm.Matrix(mus_exc)

    return sym_mod
