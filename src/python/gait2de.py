#!/usr/bin/env python

from sympy import symbols as s
import sympy.physics.mechanics as me


class BodySegment(object):

    def __init__(self, label, description, origin_joint):

        self.label = label
        self.description = description

        self.generate_symbols(label.lower())

        self.reference_frame = me.ReferenceFrame(label)

        # origin point
        self.origin = me.Point(origin_joint)

        self.mass_center = me.Point('{} mass center'.format(description))

        self.mass_center.set_pos(self.origin,
                                 self.center_of_mass_x_symbol * self.reference_frame.x +
                                 self.center_of_mass_y_symbol * self.reference_frame.y)

    def _create_symbols(self, subscript):

        self.mass_symbol = s('m'.format(subscript), real=True, positive=True)
        self.inertia_symbol = s('i'.format(subscript), real=True, positive=True)
        self.length_symbol = s('l'.format(subscript), real=True, positive=True)
        self.center_of_mass_x_symbol = s('x'.format(subscript), real=True)
        self.center_of_mass_y_symbol = s('y'.format(subscript))

    @property
    def label(self, label):
        self.generate_symbols(label.lower())


trunk = BodySegment('A', 'Trunk', 'ma', 'ia', 'la', 'xa', 'ya')

body_descriptions = {'A': 'Trunk',
                     'B': 'Right Thigh',
                     'C': 'Right Shank',
                     'D': 'Right Foot',
                     'E': 'Left Thigh',
                     'F': 'Left Shank',
                     'G': 'Left Foot'}

q, u = me.dynamicsymbols(('q:9', 'u:9'), real=True)

frames = {}
parameters = {}
for body_label in sorted(body_descriptions.keys()):
    # mass, inertia, x/y center of mass location
    par = {}
    par['length'] = symbols('l{}'.format(body_label.lower()))
    par['mass'] = symbols('m{}'.format(body_label.lower()))
    par['inertia'] = symbols('i{}'.format(body_label.lower()))
    par['mass center'] = symbols('{}x, {}y'.format(body_label.lower()))
    parameters[body_label] = par
    frames[body_label] = me.ReferenceFrame(body_label)

#par__ContactY, par__ContactHeelX, par__ContactToeX
#par__ContactStiff, par__ContactDamp, par__ContactY0, par__ContactV0, par__ContactFric

ground = me.ReferenceFrame('N')

# orient the legs
for child, parent, angle in zip(['B', 'C', 'D', 'E', 'F', 'G'],
                                ['A', 'B', 'C', 'A', 'E', 'F'],
                                q[2:]):
    reference_frames[child].orient(parent, angle, 'axis',
                                   reference_frames[parent].z)

# locate joints
origin = me.Point('O')

hip = origin.locatenew('hip', q[0] * ground.x + q[1] * ground.y)

right_knee = hip.locatenew('right knee', -parameters['B']['length'] *
                           frames['B'].y)
right_ankle = right_knee.locatenew('right ankle', -parameters['C']['length']
                                   * frames['C'].y)

left_knee = hip.locatenew('left knee', -parameters['E']['length'] *
                          frames['E'].y)
left_ankle = left_knee.locatenew('left ankle', -parameters['F']['length'] *
                                 frames['F'].y)

# locate mass centers
x, y = parameters['A']['mass center']
trunk_mass_center = hip.locatenew('trunk mass center',
                                  x * frames['A'].x + y * frames['A'].y)

x, y = parameters['B']['mass center']
right_thigh_mass_center = right_knee.locatenew('right thigh mass center',
                                               x * frames['B'].x +
                                               y * frames['B'].y)

x, y = parameters['C']['mass center']
right_shank_mass_center = right_ankle.locatenew('right shank mass center',
                                                x * frames['C'].x +
                                                y * frames['C'].y)

x, y = parameters['D']['mass center']
right_foot_mass_center = right_ankle.locatenew('right foot mass center',
                                               x * frames['D'].x +
                                               y * frames['D'].y)

x, y = parameters['E']['mass center']
left_thigh_mass_center = left_knee.locatenew('left thigh mass center',
                                             x * frames['E'].x +
                                             y * frames['E'].y)

x, y = parameters['F']['mass center']
left_shank_mass_center = left_ankle.locatenew('left shank mass center',
                                              x * frames['F'].x +
                                              y * frames['F'].y)

x, y = parameters['G']['mass center']
left_foot_mass_center = left_ankle.locatenew('left foot mass center',
                                             x * frames['G'].x +
                                             y * frames['G'].y)

# angular velocities
for body_label, angular_rate in zip(['A', 'B', 'C', 'D', 'E'], u[2:]):
    frames[body_label].set_ang_vel(angular_rate * ground.z)
