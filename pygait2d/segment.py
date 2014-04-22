#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import symbols, exp, acos, pi
import sympy.physics.mechanics as me
from pydy.viz import VisualizationFrame, Cylinder, Sphere


class BodySegment(object):

    viz_sphere_radius = 0.07
    viz_cylinder_radius = 0.035

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
            The point where this segment is connected to its parent.
        joint_description : string
            A short description of a new joint point, e.g. 'knee', for this
            segemt to attach to its child segment.
        inertial_frame : sympy.physics.mechanics.ReferenceFrame
            The global inertial reference frame of the system. This is used
            to apply gravity to the segment (in the negative y direction of
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
        self.g = symbols('g')
        self.mass_symbol = symbols('m{}'.format(subscript))
        self.inertia_symbol = \
            symbols('i{}'.format(subscript))
        self.length_symbol = \
            symbols('l{}'.format(subscript))
        self.mass_center_x_symbol = symbols('x{}'.format(subscript))
        self.mass_center_y_symbol = symbols('y{}'.format(subscript))

        self.constants = [self.g,
                          self.mass_symbol,
                          self.inertia_symbol,
                          self.length_symbol,
                          self.mass_center_x_symbol,
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
        """Creates a list of the kinematic differential equations. This is
        the simple definition:

        0 = \dot{q}_i - u_i

        """
        self.kinematic_equations = \
            [self.generalized_coordinate_derivative_symbol -
             self.generalized_speed_symbol]

    def _orient(self):
        """Generates and orients the segment's reference frame relative to
        the parent reference frame by body fixed simple rotation about the
        generalized coordinate."""
        self.reference_frame = \
            self.parent_reference_frame.orientnew(
                self.label, 'Axis', (self.generalized_coordinate_symbol,
                                     self.parent_reference_frame.z))

    def _set_angular_velocity(self):
        """Sets the angular velocity with the generalized speed."""
        self.reference_frame.set_ang_vel(self.parent_reference_frame,
                                         self.generalized_speed_symbol *
                                         self.parent_reference_frame.z)

    def _locate_joint(self):
        """Creates a point with respect to the origin joint for the next
        joint in the segment. This assumes that new joint is in the negative
        y direction with repect to the origin joint."""
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
                                     self.reference_frame)
        self.joint.v2pt_theory(self.origin_joint, self.inertial_frame,
                               self.reference_frame)

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
        # TODO : add in passive joint stiffness and damping

    def _gravity(self):
        """Creates the gravitational force vector acting on the segment."""
        self.gravity = -self.mass_symbol * self.g * self.inertial_frame.y

    def visualization_frames(self):
        """Returns visualization frames for the animation of the system.
        The requires numerical values of the cylinders and spheres."""

        viz_frames = []

        cylinder = Cylinder(color='red',
                            length=self.length_symbol,
                            radius=self.viz_cylinder_radius)

        center_point = self.origin_joint.locatenew('Cylinder Center',
                                                   -self.length_symbol / 2 *
                                                   self.reference_frame.y)

        viz_frames.append(VisualizationFrame('VizFrame',
                                             self.reference_frame,
                                             center_point, cylinder))

        viz_frames.append(VisualizationFrame('OriginJointFrame',
                                             self.reference_frame,
                                             self.origin_joint,
                                             Sphere(color='blue',
                                                    radius=self.viz_sphere_radius)))

        return viz_frames


class TrunkSegment(BodySegment):
    def __init__(self, *args):
        super(TrunkSegment, self).__init__(*args)
        self._trunk_extra_kinematic_equations()

    def _create_symbols(self):
        super(TrunkSegment, self)._create_symbols()
        # TODO : Format these with the subscript instead of a directly.
        self.qa = me.dynamicsymbols('qax, qay')
        self.ua = me.dynamicsymbols('uax, uay')
        self.constants.remove(self.length_symbol)
        del self.length_symbol

    def _trunk_extra_kinematic_equations(self):
        qaxd, qayd = me.dynamicsymbols('qax, qay', 1)
        self.kinematic_equations += [self.ua[0] - qaxd, self.ua[1] - qayd]

    def _locate_joint(self):
        """The trunk only has one joint, the hip, there is no other point."""
        # This locates the hip joint relative to the ground origin point.
        self.joint = self.origin_joint.locatenew(self.joint_description,
                                                 self.qa[0] *
                                                 self.inertial_frame.x +
                                                 self.qa[1] *
                                                 self.inertial_frame.y)

    def _locate_mass_center(self):
        """Creates a point with respect the hip joint for the mass center
        of the segment."""
        self.mass_center = self.joint.locatenew(
            '{} mass center'.format(self.description),
            self.mass_center_x_symbol * self.reference_frame.x +
            self.mass_center_y_symbol * self.reference_frame.y)

    def _set_linear_velocities(self):
        """Sets the linear velocities of the mass center and new joint."""
        # The joint is the hip. The origin joint is the ground's origin.
        self.joint.set_vel(self.inertial_frame, self.ua[0] *
                           self.inertial_frame.x + self.ua[1] *
                           self.inertial_frame.y)
        self.mass_center.v2pt_theory(self.joint, self.inertial_frame,
                                     self.reference_frame)

    def visualization_frames(self):
        """This should go from the hip to the mass center."""

        viz_frames = []

        hip_to_mc_vector = self.mass_center.pos_from(self.joint)

        cylinder = Cylinder(color='red', length=hip_to_mc_vector.magnitude(),
                            radius=self.viz_cylinder_radius)

        center_point = \
            self.joint.locatenew('Cylinder Center',
                                 hip_to_mc_vector.magnitude() / 2 *
                                 hip_to_mc_vector.normalize())

        viz_frames.append(VisualizationFrame('VizFrame',
                                             self.reference_frame,
                                             center_point, cylinder))

        viz_frames.append(VisualizationFrame('MassCenterFrame',
                                             self.reference_frame,
                                             self.mass_center,
                                             Sphere(color='blue',
                                                    radius=self.viz_sphere_radius)))

        return viz_frames


class FootSegment(BodySegment):

    viz_sphere_radius = 0.03
    viz_cylinder_radius = 0.01

    def __init__(self, *args):
        super(FootSegment, self).__init__(*args)
        self._locate_foot_points()
        self._set_foot_linear_velocities()

    def _create_symbols(self):
        super(FootSegment, self)._create_symbols()
        self.heel_distance = symbols('hx{}'.format(self.label.lower()))
        self.toe_distance = symbols('tx{}'.format(self.label.lower()))
        self.foot_depth = symbols('fy{}'.format(self.label.lower()))

        self.constants.remove(self.length_symbol)
        del self.length_symbol

        self.constants += [self.heel_distance, self.toe_distance,
                           self.foot_depth]

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
                                     self.reference_frame)
        self.heel.v2pt_theory(self.origin_joint, self.inertial_frame,
                              self.reference_frame)
        self.toe.v2pt_theory(self.origin_joint, self.inertial_frame,
                             self.reference_frame)

    def visualization_frames(self):
        """Returns a list of visualization frames needed to visualize the
        foot."""

        viz_frames = []

        heel_to_toe_length = self.toe.pos_from(self.heel).magnitude()
        bottom_cylinder = Cylinder(color='red',
                                   length=heel_to_toe_length,
                                   radius=self.viz_cylinder_radius)
        bottom_center_point = self.heel.locatenew('BottomCenter',
                                                  heel_to_toe_length / 2 *
                                                  self.reference_frame.x)
        # Creates a reference frame with the Y axis pointing from the heel
        # to the toe.
        bottom_rf = self.reference_frame.orientnew('Bottom', 'Axis',
                                                   (-pi / 2,
                                                    self.reference_frame.z))
        viz_frames.append(VisualizationFrame('BottomVizFrame',
                                             bottom_rf,
                                             bottom_center_point,
                                             bottom_cylinder))

        # top of foot
        ankle_to_toe_vector = self.toe.pos_from(self.origin_joint)
        top_cylinder = Cylinder(color='red',
                                length=ankle_to_toe_vector.magnitude(),
                                radius=self.viz_cylinder_radius)
        angle = -acos(ankle_to_toe_vector.normalize().dot(bottom_rf.y))
        top_foot_rf = bottom_rf.orientnew('Top', 'Axis', (angle,
                                                          bottom_rf.z))
        top_foot_center_point = \
            self.origin_joint.locatenew('TopCenter',
                                        ankle_to_toe_vector.magnitude() / 2
                                        * top_foot_rf.y)
        viz_frames.append(VisualizationFrame('TopVizFrame', top_foot_rf,
                                             top_foot_center_point,
                                             top_cylinder))

        # back of foot
        heel_to_ankle_vector = self.origin_joint.pos_from(self.heel)
        back_cylinder = Cylinder(color='red',
                                 length=heel_to_ankle_vector.magnitude(),
                                 radius=self.viz_cylinder_radius)
        angle = acos(heel_to_ankle_vector.normalize().dot(bottom_rf.y))
        back_foot_rf = bottom_rf.orientnew('Back', 'Axis', (angle,
                                                            bottom_rf.z))
        back_foot_center_point = \
            self.heel.locatenew('BackCenter',
                                heel_to_ankle_vector.magnitude() / 2 *
                                back_foot_rf.y)
        viz_frames.append(VisualizationFrame('BackVizFrame', back_foot_rf,
                                             back_foot_center_point,
                                             back_cylinder))

        # spheres for the ankle, toe, and heel
        viz_frames.append(VisualizationFrame('AnkleVizFrame',
                                             self.reference_frame,
                                             self.origin_joint,
                                             Sphere(color='blue',
                                                    radius=self.viz_sphere_radius)))
        viz_frames.append(VisualizationFrame('ToeVizFrame',
                                             self.reference_frame, self.toe,
                                             Sphere(color='blue',
                                                    radius=self.viz_sphere_radius)))
        viz_frames.append(VisualizationFrame('HeelVizFrame',
                                             self.reference_frame, self.heel,
                                             Sphere(color='blue',
                                                    radius=self.viz_sphere_radius)))

        return viz_frames


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
