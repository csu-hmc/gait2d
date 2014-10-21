#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This example simply simulates and visualizes the uncontrolled motion and
the model "falls down"."""

import numpy as np
from scipy.integrate import odeint
from pydy.codegen.code import generate_ode_function
from pydy.viz import Scene

from pygait2d import derive, simulate

(mass_matrix, forcing_vector, kane, constants, coordinates, speeds,
 specified, visualization_frames, ground, origin) = derive.derive_equations_of_motion()

rhs = generate_ode_function(mass_matrix, forcing_vector,
                            constants, coordinates, speeds,
                            specified=specified, generator='cython')

constant_values = simulate.load_constants(constants,
                                          'data/example_constants.yml')

args = {'constants': np.array(constant_values.values()),
        'specified': np.zeros(9)}

time_vector = np.linspace(0.0, 10.0, num=1000)
initial_conditions = np.zeros(18)
initial_conditions[1] = 1.0  # set hip above ground
initial_conditions[3] = np.deg2rad(5.0)  # right hip angle
initial_conditions[6] = -np.deg2rad(5.0)  # left hip angle
trajectories = odeint(rhs, initial_conditions, time_vector, args=(args,))

scene = Scene(ground, origin, *visualization_frames)

scene.generate_visualization_json(coordinates + speeds, constants,
                                  trajectories, args['constants'])

scene.display()
