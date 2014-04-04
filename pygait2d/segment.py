#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import symbols, exp
import sympy.physics.mechanics as me


class BodySegment(object):

    def __init__(self, label, description, parent_reference_frame,
                 origin_joint, joint_description, inertial_frame):
        """Initializes a body segment.

        Parameters
        ==========
        label : string
            A short label for the segment, like 'A'.
        description : string
            A short description of the segment, like 'Trunk'.
        parent_reference_frame : sympy.physics.vector.ReferenceFrame
            The parent reference frame for this segment.
        origin_joint : sympy.physics.vector.Point
            The joint which connects this segment to its parent.
        joint_description : string
            A short description of the new joint, e.g. 'knee'.
        inertial_frame : sympy.physics.mechanics.ReferenceFrame
            The inertial reference frame the segment is in. This is used to
            apply gravity to the segment (in the negative y direction of
            this frame).

        """

        self.label = label
        self.description = description
        self.parent_reference_frame = parent_reference_frame
        self.origin_joint = origin_joint
        self.joint_description = joint_description
        self.inertial_frame = inertial_frame

        self._create_symbols()
        self._kinematic_differential_equations()
        self._orient()
        self._set_angular_velocity()
        self._locate_joint()
        self._locate_mass_center()
        self._set_linear_velocities()
        self._inertia_dyadic()
        self._create_rigid_body()
        self._joint_torque()
        self._gravity()

    def _create_symbols(self):
        """Generates all of the SymPy symbols and functions of time
        associated with this segment."""

        subscript = self.label.lower()

        # constants
        self.g = symbols('g', real=True)
        self.mass_symbol = symbols('m{}'.format(subscript), real=True,
                                   positive=True)
        self.inertia_symbol = \
            symbols('i{}'.format(subscript), real=True, positive=True)
        self.length_symbol = \
            symbols('l{}'.format(subscript), real=True, positive=True)
        self.mass_center_x_symbol = symbols('x{}'.format(subscript), real=True)
        self.mass_center_y_symbol = symbols('y{}'.format(subscript))

        self.constants = [self.g, self.mass_symbol, self.inertia_symbol,
                          self.length_symbol, self.mass_center_x_symbol,
                          self.mass_center_y_symbol]

        # functions of time
        self.generalized_coordinate_symbol = \
            me.dynamicsymbols('q{}'.format(subscript))
        self.generalized_coordinate_derivative_symbol = \
            me.dynamicsymbols('q{}'.format(subscript), 1)
        self.generalized_speed_symbol = \
            me.dynamicsymbols('u{}'.format(subscript))
        self.joint_torque_symbol = me.dynamicsymbols('T{}'.format(subscript))

    def _kinematic_differential_equations(self):
        """Creates a list of the kinematic differential equations. We simply
        chose u = qdot."""
        self.kinematic_equations = \
            [self.generalized_coordinate_derivative_symbol -
             self.generalized_speed_symbol]

    def _orient(self):
        """Generates and orients the segment's reference frame relative to
        the parent reference frame by body fixed simple rotation about the
        generalized coordinate."""
        self.reference_frame = self.parent_reference_frame.orientnew(
            self.label, 'Axis', (self.generalized_coordinate_symbol,
                                 self.parent_reference_frame.z))

    def _set_angular_velocity(self):
        """Sets the angular velocity with the generalized speed."""
        self.reference_frame.set_ang_vel(self.parent_reference_frame,
                                         self.generalized_speed_symbol *
                                         self.parent_reference_frame.z)

    def _locate_joint(self):
        """Creates a point with respect to the origin joint for the next
        joint in the segment."""
        self.joint = self.origin_joint.locatenew(self.joint_description,
                                                 -self.length_symbol *
                                                 self.reference_frame.y)

    def _locate_mass_center(self):
        """Creates a point with respect the origin joint for the mass center
        of the segment."""
        self.mass_center = self.origin_joint.locatenew(
            '{} mass center'.format(self.description),
            self.mass_center_x_symbol * self.reference_frame.x +
            self.mass_center_y_symbol * self.reference_frame.y)

    def _set_linear_velocities(self):
        """Sets the linear velocities of the mass center and new joint."""
        self.mass_center.v2pt_theory(self.origin_joint, self.inertial_frame,
                                     self.parent_reference_frame)
        self.joint.v2pt_theory(self.origin_joint, self.inertial_frame,
                               self.parent_reference_frame)

    def _inertia_dyadic(self):
        """Creates an inertia dyadic for the segment."""
        self.inertia_dyadic = me.inertia(self.reference_frame, 0, 0,
                                         self.inertia_symbol)

    def _create_rigid_body(self):
        """Creates a rigid body for the segment."""
        self.rigid_body = me.RigidBody(self.description, self.mass_center,
                                       self.reference_frame,
                                       self.mass_symbol,
                                       (self.inertia_dyadic, self.mass_center))

    def _joint_torque(self):
        """Creates the joint torque vector acting on the segment."""
        self.torque = self.joint_torque_symbol * self.reference_frame.z

        # TODO : This is the torque vector which is applied to this segment,
        # but the negative of it should be applied to the parent segement.

    def _gravity(self):
        """Creates the gravitational force vector acting on the segment."""
        self.gravity = -self.mass_symbol * self.g * self.inertial_frame.y


class TrunkSegment(BodySegment):
    def __init__(self, *args):
        super(TrunkSegment, self).__init__(*args)
        self._trunk_extra_kinematic_equations()

    def _create_symbols(self):
        super(TrunkSegment, self)._create_symbols()
        # TODO : Format these with the subscript instead of a directly.
        self.qa = me.dynamicsymbols('qax, qay')
        self.ua = me.dynamicsymbols('uax, uay')

    def _trunk_extra_kinematic_equations(self):
        qaxd, qayd = me.dynamicsymbols('qax, qay', 1)
        self.kinematic_equations += [self.ua[0] - qaxd, self.ua[1] - qayd]

    def _locate_joint(self):
        """The trunk only has one joint, the hip, there is no other point."""
        self.joint = self.origin_joint.locatenew(self.joint_description,
                                                 self.qa[0] *
                                                 self.inertial_frame.x +
                                                 self.qa[1] *
                                                 self.inertial_frame.y)

    def _set_linear_velocities(self):
        """Sets the linear velocities of the mass center and new joint."""
        self.joint.set_vel(self.inertial_frame, self.ua[0] *
                           self.inertial_frame.x + self.ua[1] *
                           self.inertial_frame.y)
        self.mass_center.v2pt_theory(self.joint, self.inertial_frame,
                                     self.parent_reference_frame)


class FootSegment(BodySegment):
    def __init__(self, *args):
        super(FootSegment, self).__init__(*args)
        self._locate_foot_points()
        self._set_foot_linear_velocities()

    def _create_symbols(self):
        super(FootSegment, self)._create_symbols()
        self.heel_distance = symbols('hx{}'.format(self.label.lower()))
        self.toe_distance = symbols('tx{}'.format(self.label.lower()))
        self.foot_depth = symbols('fy{}'.format(self.label.lower()))

        self.constants += [self.heel_distance, self.toe_distance, self.foot_depth]

    def _locate_joint(self):
        """The foot has no joint."""
        pass

    def _locate_foot_points(self):

        self.heel = self.origin_joint.locatenew(
            '{} heel'.format(self.description), self.heel_distance *
            self.reference_frame.x + self.foot_depth *
            self.reference_frame.y)

        self.toe = self.origin_joint.locatenew(
            '{} toe'.format(self.description), self.toe_distance *
            self.reference_frame.x + self.foot_depth *
            self.reference_frame.y)

    def _set_linear_velocities(self):
        """There is no joint so pass this."""
        pass

    def _set_foot_linear_velocities(self):
        """Sets the linear velocities of the mass center and new joint."""
        self.mass_center.v2pt_theory(self.origin_joint, self.inertial_frame,
                                     self.parent_reference_frame)
        self.heel.v2pt_theory(self.origin_joint, self.inertial_frame,
                              self.parent_reference_frame)
        self.toe.v2pt_theory(self.origin_joint, self.inertial_frame,
                             self.parent_reference_frame)


def contact_force(point, ground, origin):
    """Returns a contact force vector acting on the given point made of
    friction along the contact surface and elastic force in the vertical
    direction.

    Parameters
    ==========
    point : sympy.physics.mechanics.Point
        The point which the contact force should be computed for.
    ground : sympy.physics.mechanics.ReferenceFrame
        A reference frame which represents the inerital ground in 2D space.
        The x axis defines the ground line and positive y is up.
    origin : sympy.physics.mechanics.Point
        An origin point located on the ground line.

    Returns
    =======
    force : sympy.physics.mechanics.Vector
        The contact force between the point and the ground.

    """
    # This is the "height" of the point above the ground, where a negative
    # value means that the point is below the ground.
    y_location = point.pos_from(origin).dot(ground.y)

    # The penetration into the ground is mathematically defined as:
    #
    #               { 0 if y_location > 0
    # deformation = {
    #               { abs(y_location) if y_location < 0
    #

    penetration = (abs(y_location) - y_location) / 2

    velocity = point.vel(ground)

    # The addition of "- y_location" here adds a small linear term to the
    # cubic stiffness and creates a light attractive force torwards the
    # ground. This is in place to ensure that gradients can be computed for
    # the optimization used in Ackermann and van den Bogert 2010.
    contact_stiffness, contact_damping = symbols('kc, cc')
    contact_friction_coefficient, friction_scaling_factor = symbols('mu, vs')

    vertical_force = (contact_stiffness * penetration ** 3 - y_location) * \
        (1 - contact_damping * velocity.dot(ground.y))

    friction = -contact_friction_coefficient * vertical_force * \
        ((2 / (1 + exp(-velocity.dot(ground.x) /
                       friction_scaling_factor))) - 1)

    return friction * ground.x + vertical_force * ground.y
