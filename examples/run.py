#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This example simply simulates and visualizes the uncontrolled motion and
the model "falls down"."""

import os

import numpy as np
from scipy.integrate import odeint
from pydy.codegen.ode_function_generators import generate_ode_function
from pydy.viz import Scene
import matplotlib.pyplot as plt

from pygait2d import derive, simulate
from pygait2d.utils import generate_animation

symbolics = derive.derive_equations_of_motion()

constant_values = simulate.load_constants(
    symbolics.constants, os.path.join(os.path.dirname(__file__), '..',
                                      'data/example_constants.yml'))

rhs = generate_ode_function(
    symbolics.kanes_method.forcing_full,
    symbolics.coordinates,
    symbolics.speeds,
    constants=list(constant_values.keys()),
    mass_matrix=symbolics.kanes_method.mass_matrix_full,
    specifieds=symbolics.specifieds,
    generator='cython',
    constants_arg_type='array',
    specifieds_arg_type='array',
)

args = (np.zeros(len(symbolics.specifieds)),
        np.array(list(constant_values.values())))

time_vector = np.linspace(0.0, 2.0, num=1000)
initial_conditions = np.zeros(len(symbolics.states))
initial_conditions[1] = 1.0  # set hip above ground
initial_conditions[3] = np.deg2rad(5.0)  # right hip angle
initial_conditions[6] = -np.deg2rad(5.0)  # left hip angle
trajectories = odeint(rhs, initial_conditions, time_vector, args=args)

ani = generate_animation(symbolics,
                         time_vector,
                         trajectories,
                         np.zeros((len(time_vector),
                                   len(symbolics.specifieds))),
                         np.array(list(constant_values.values())))

plt.show()

scene = Scene(symbolics.inertial_frame, symbolics.origin,
              *symbolics.viz_frames)
scene.states_symbols = symbolics.states
scene.constants = constant_values
scene.states_trajectories = trajectories
scene.times = time_vector
