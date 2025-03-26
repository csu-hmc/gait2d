#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This example simply simulates and visualizes the uncontrolled motion and
the model "falls down"."""

import os

import numpy as np
from scipy.integrate import odeint
from pydy.codegen.ode_function_generators import generate_ode_function
from pydy.viz import Scene

from pygait2d import derive, simulate

(mass_matrix, forcing_vector, kane, constants, coordinates, speeds, specified,
 visualization_frames, ground, origin, segments) = \
    derive.derive_equations_of_motion(ground_force_on_all_joints=False)

constant_values = simulate.load_constants(
    constants, os.path.join(os.path.dirname(__file__), '..',
                            'data/example_constants.yml'))

rhs = generate_ode_function(
    forcing_vector,
    coordinates,
    speeds,
    constants=list(constant_values.keys()),
    mass_matrix=mass_matrix,
    specifieds=specified,
    generator='cython',
    constants_arg_type='array',
    specifieds_arg_type='array',
)


args = (np.zeros(len(specified)), np.array(list(constant_values.values())))

time_vector = np.linspace(0.0, 10.0, num=1000)
initial_conditions = np.zeros(len(coordinates + speeds))
initial_conditions[1] = 1.0  # set hip above ground
initial_conditions[3] = np.deg2rad(5.0)  # right hip angle
initial_conditions[6] = -np.deg2rad(5.0)  # left hip angle
trajectories = odeint(rhs, initial_conditions, time_vector, args=args)

scene = Scene(ground, origin, *visualization_frames)
scene.states_symbols = coordinates + speeds
scene.constants = constant_values
scene.states_trajectories = trajectories
scene.times = time_vector
