#!/usr/bin/env python

"""This ensures that the PyDy model gives the same result as the Autolev
model."""

import os

import yaml
import numpy as np
from numpy import testing
from pydy.codegen.ode_function_generators import generate_ode_function
from algait2de.gait2de import evaluate_autolev_rhs as autolev_rhs

# local imports
from .. import derive, simulate

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')


def test_accelerations():
    (mass_matrix, forcing_vector, kane, constants, coordinates, speeds,
     specified, visualization_frames, ground, origin, segments) = \
        derive.derive_equations_of_motion()

    pydy_rhs = generate_ode_function(
        forcing_vector,
        coordinates,
        speeds,
        constants=constants,
        mass_matrix=mass_matrix,
        specifieds=specified,
        generator='cython',
        constants_arg_type='array',
        specifieds_arg_type='array',
    )

    coordinate_values = np.random.random(9)
    speed_values = np.random.random(9)
    specified_values = np.random.random(9)

    constant_map = simulate.load_constants(
        constants, os.path.join(ROOT, 'data/example_constants.yml'))
    args = (specified_values, np.array(list(constant_map.values())))

    x = np.hstack((coordinate_values, speed_values))
    pydy_xdot = pydy_rhs(x, 0.0, *args)

    with open(os.path.join(ROOT, 'data/example_constants.yml'), 'r') as f:
        constants_dict = yaml.load(f, Loader=yaml.SafeLoader)

    constants_dict = simulate.map_values_to_autolev_symbols(constants_dict)

    accelerations, grfs, animation_data = \
        autolev_rhs(coordinate_values, speed_values, specified_values,
                    constants_dict)

    testing.assert_allclose(pydy_xdot[9:], accelerations)
