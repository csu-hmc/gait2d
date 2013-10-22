#!/usr/bin/env python

from sympy import symbols as s
import sympy.physics.mechanics as me


class BodySegment(object):

    def __init__(self, label, description, mass, inertia, length,
                 center_of_mass):
        self.label = label
        self.description = description
        self.mass = mass
        self.inertia = inertia
        self.length = length
        self.center_of_mass_x = self.center_of_mass_x
        self.center_of_mass_y = self.center_of_mass_y

trunk = BodySegment('A', 'Trunk', s('ma'), s('ia'), s('la'), s('xa'), s('ya'))

body_descriptions = {'A': 'Trunk',
                     'B': 'Right Thigh',
                     'C': 'Right Shank',
                     'D': 'Right Foot',
                     'E': 'Left Thigh',
                     'F': 'Left Shank',
                     'G': 'Left Foot'}

q = me.dynamicsymbols('q:9')
u = me.dynamicsymbols('u:9')

reference_frames = {}
parameters = {}
for body_label in sorted(body_descriptions.keys()):
    # mass, inertia, x/y center of mass location
    par = {}
    par['length'] = symbols('l{}'.format(body_label.lower()))
    par['mass'] = symbols('m{}'.format(body_label.lower()))
    par['inertia'] = symbols('i{}'.format(body_label.lower()))
    par['mass center'] = symbols('{}x, {}y'.format(body_label.lower()))
    parameters[body_label] = par
    reference_frames[body_label] = me.ReferenceFrame(body_label)

#par__ContactY, par__ContactHeelX, par__ContactToeX
#par__ContactStiff, par__ContactDamp, par__ContactY0, par__ContactV0, par__ContactFric

ground = me.ReferenceFrame('N')

# orient the legs
for child, parent, angle in zip(['B', 'C', 'D', 'E', 'F', 'G'],
                                ['A', 'B', 'C', 'A', 'E', 'F'],
                                q[2:]):
    reference_frames[child].orient(parent, angle, 'axis',
                                   reference_frames[parent].z)

# locate
origin = me.Point('O')

hip = origin.locatenew('Pa', q[0] * ground.x + q[1] * ground.y)
right_knee = hip.locatenew('Pb', -parameters['B']['length'] * thigh.y)
right_ankle = right_knee.locatenew('Pc', -parameters['C']['length'] * shank.y )

leftt_knee = hip.locatenew('Pb', -parameters['E']['length'] * thigh.y)
left_ankle = left_knee.locatenew('Pc', -parameters['F']['length'] * shank.y )

x, y = parameters['B']['mass center']
right_thigh_mass_center = hip.locatenew('Pb', x * thigh.x + -y * thigh.y)
x, y = parameters['C']['mass center']
right_shank_mass_center = right_knee.locatenew('Pc', -parameters['C']['length'] * shank.y )

leftt_knee = hip.locatenew('Pb', -parameters['E']['length'] * thigh.y)
left_ankle = left_knee.locatenew('Pc', -parameters['F']['length'] * shank.y )


