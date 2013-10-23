#!/usr/bin/env python

from sympy import symbols as s
from sympy import exp
import sympy.physics.mechanics as me


class BodySegment(object):

    def __init__(self, label, description, parent_reference_frame,
                 origin_joint, joint_description, inertial_frame):
        """Initializes a body segment.

        Parameters
        ==========
        label : string
            A short label fro the segment, like 'A'.
        description : string
            A short description of the segement, like 'Trunk')
        parent_reference_frame : sympy.physics.mechanics.ReferenceFrame
            The parent reference frame for this segment.
        origin_joint : sympy.physics.mechanics.Point
            The joint which connects this segment to its parent.
        joint_description : string
            A short description of the new joint, e.g. 'knee'.
        inertial_frame : sympy.physics.mechanics.ReferenceFrame
            The inertial reference frame the segment is in.

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
        self._inertia_dyad()
        self._create_rigid_body()
        self._joint_torque()
        self._gravity()

    def _create_symbols(self):
        """Generates all of the SymPy symbols associated with this segment."""

        subscript = self.label.lower()
        self.g = s('g', real=True)
        self.mass_symbol = s('m{}'.format(subscript), real=True, positive=True)
        self.inertia_symbol = \
            s('i{}'.format(subscript), real=True, positive=True)
        self.length_symbol = \
            s('l{}'.format(subscript), real=True, positive=True)
        self.center_of_mass_x_symbol = s('x{}'.format(subscript), real=True)
        self.center_of_mass_y_symbol = s('y{}'.format(subscript))
        self.joint_torque_symbol = me.dynamicsymbols('M{}'.format(subscript))

        self.generalized_coordinate_symbol = \
            me.dynamicsymbols('q'.format(subscript))
        self.generalized_coordinate_derivative_symbol = \
            me.dynamicsymbols('q'.format(subscript), 1)
        self.generalized_speed_symbol = \
            me.dynamicsymbols('u{}'.format(subscript))

    def _kinematic_differential_equations(self):
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
        self.reference_frame.set_ang_vel(self.generalized_speed_symbol *
                                         self.parent_reference_frame.z)

    def _locate_joint(self):
        """Creates a point with respect to the origin joint for the next
        joint in the segment."""
        self.joint = self.origin_joint.locatenew(self.joint_description,
                                                 -self.length_symbol *
                                                 self.referenc_frame.y)

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

    def _inertia_dyad(self):
        """Creates an inertia dyad for the segment."""
        self.inertia = me.inertia(self.reference_frame, 0, 0,
                                  self.inertia_symbol)

    def _create_rigid_body(self):
        """Creates a rigid body for the segment."""
        self.rigid_body = me.RigidBody(self.description, self.mass_center,
                                       self.reference_frame,
                                       self.mass_symbol,
                                       (self.inertia, self.mass_center))

    def _joint_torque(self):
        """Creates the joint torque acting on the segment."""
        self.torque = self.joint_torque_symbol * self.reference_frame.z

    def _gravity(self):
        self.gravity = -self.mass_symbol * self.g * self.inertial.y


class TrunkSegment(BodySegment):
    def _create_symbols(self):
        super(TrunkSegment, self)._create_symbols()
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
        self.joint.set_vel(self.ua[0] * ground.x + self.ua[2] * ground.y)
        self.mass_center.v2pt_theory(self.joint, self.inertial_frame,
                                     self.parent_reference_frame)


class FootSegment(BodySegment):
    def __init__(self):
        super(FootSegment, self).__init__()
        self._locate_foot_points()
        self._set_foot_linear_velocities()

    def _create_symbols(self):
        super(FootSegment, self)._create_symbols()
        self.heel_distance = s('hx{}'.format(self.label.lower()))
        self.toe_distance = s('tx{}'.format(self.label.lower()))
        self.foot_depth = s('fy{}'.format(self.label.lower()))

    def _locate_joint(self):
        """The foot has no joint."""
        pass

    def _locate_foot_points(self):

        self.heel = self.origin_joint.locatenew(self.heel_distance *
                                                self.reference_frame.x +
                                                self.foot_depth *
                                                self.reference_frame.y)
        self.toe = self.origin_joint.locatenew(self.toe_distance *
                                               self.reference_frame.x +
                                               self.foot_depth *
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


def contact_force(point, ground):
    """Returns a contact force vector acting on the given point made of
    friction along the contact surface and elastic force in the vertical
    direction.

    Parameters
    ==========
    point : sympy.physics.mechanics.Point
    ground : sympy.physics.mehcanics.ReferenceFrame

    Returns
    =======
    force : sympy.physics.mechanics.Vector

    """
    # This is the "height" of the point above the ground, where a negative
    # value means that the point is below the ground.
    y_location = origin.pos_from(point).dot(ground.y)

    # The deformation, i.e., penetration into the ground is mathematically
    # defined as:
    #               { 0 if y_location > 0
    # deformation = {
    #               { abs(y_location) if y_location < 0
    # TODO: Ideally this would be translated to an if statement in the
    # generated numerical code but this is a hack to get around that.

    deformation = (abs(y_location) - y_location) / 2

    velocity = point.vel(ground)

    # The addition of "- y_location" here adds a small linear term to the
    # cubic stiffness and creates a light attractive force torwards the
    # ground. This is in place to ensure that gradients can be computed for
    # the optimization used in Ackermann and van den Bogert 2010.
    contact_stiffness, contact_damping = s('k, c')
    contact_friction_coefficient, friction_scaling_factor = s('mu, vs')

    vertical_force = (contact_stiffness * deformation ** 3 - y_location) * \
        (1 - contact_damping * velocity.dot(ground.y))

    friction = -contact_friction_coefficient * vertical_force * \
        ((2 / (1 + exp(-velocity.dot(ground.x) / friction_scaling_factor))) - 1)

    return friction * ground.x + vertical_force * ground.y

segment_descriptions = {'A': (TrunkSegment, 'Trunk', 'Hip'),
                        'B': (BodySegment, 'Right Thigh', 'Knee'),
                        'C': (BodySegment, 'Right Shank', 'Ankle'),
                        'D': (FootSegment, 'Right Foot', 'Heel'),
                        'E': (BodySegment, 'Left Thigh', 'Knee'),
                        'F': (BodySegment, 'Left Shank', 'Ankle'),
                        'G': (FootSegment, 'Left Foot', 'Heel')}


ground = me.ReferenceFrame('N')
origin = me.Point('O')
origin.set_vel(ground, 0)

segments = []
generalized_coordinates = []
generalized_speeds = []
kinematic_equations = []
generalized_forces = []
bodies = []
for label, (cla, desc, joint_desc) in segment_descriptions.items():
    if label == 'A':
        parent_reference_frame = ground
        origin_joint = origin
    elif label == 'E':
        parent_reference_frame = segments[0].reference_frame
        origin_joint = segments[0].joint
    else:
        parent_reference_frame = segments[-1].reference_frame
        origin_joint = segments[-1].joint

    segment = cla(label, desc, parent_reference_frame, origin_joint,
                  joint_desc, ground)
    segments.append(segment)

    # coordinates, speeds, kinematic differential equations
    generalized_coordinates.append(segment.generalized_coordinate_symbol)
    generalized_speeds.append(segment.generalized_speed_symbol)
    kinematic_equations += segment.kinematic_equations
    # gravity
    generalized_forces.append((segment.mass_center, segment.gravity))
    # joint torques
    generalized_forces.append((segment.reference_frame, segment.torque))
    generalized_forces.append((segment.parent_reference_frame,
                               -segment.torque))
    # contact
    generalized_forces.append((segment.joint, contact_force(segment.joint,
                                                            ground)))

    # bodies
    bodies.append(segment.rigid_body)


# Add contact force for trunk mass center.
generalized_forces.append((segments[0].mass_center,
                           contact_force(segments[0].mass_center, ground)))
# Add hand of god
trunk_force_x, trunk_force_y = s('Fax, Fay')
generalized_forces.append((segments[0].mass_center, trunk_force_x * ground.x
                           + trunk_force_y * ground.y))

# equations of motion
kane = me.KanesMethod(ground, generalized_coordinates, generalized_speeds,
                      kinematic_equations)
kane.kanes_equations(generalized_forces, bodies)

# outputs
#def ground_reaction_force():
    #vertical_force = heel_force + toe_force
    #moment = heel.pos_from(origin).cross(heel_force) + \
             #toe.pos_from(origin).cross(toe_force)
    #return vertical_force, moment
