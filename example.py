#!/usr/bin/env python

"""This ensures that the PyDy model gives the same result as the Autolev
model."""

import yaml
import numpy as np
from numpy import testing
from pydy_code_gen.code import generate_ode_function

from algait2de.gait2de import evaluate_autolev_rhs as autolev_rhs
from pygait2d import derive, simulate

kane, constants, coordinates, speeds, specified = \
    derive.derive_equations_of_motion()

pydy_rhs = generate_ode_function(kane.mass_matrix_full, kane.forcing_full,
                                 constants, coordinates, speeds,
                                 specified=specified, generator='cython')

coordinate_values = np.random.random(9)
speed_values = np.random.random(9)
specified_values = np.random.random(9)

coordinate_values = np.arange(9.0)
speed_values = np.arange(9.0)
specified_values = np.arange(9.0)

constant_values = simulate.load_constants('data/example_constants.yml')
args = {'constants': np.array([constant_values[c] for c in constants]),
        'specified': specified_values}

pydy_xdot = pydy_rhs(np.hstack((coordinate_values, speed_values)), 0.0, args)

with open('data/example_constants.yml', 'r') as f:
    constants_dict = yaml.load(f)

constants_dict = simulate.map_values_to_autolev_symbols(constants_dict)

accelerations, grfs, animation_data = \
    autolev_rhs(coordinate_values, speed_values, specified_values,
                constants_dict)

testing.assert_allclose(pydy_xdot[9:], accelerations)
