#!/usr/bin/env python

"""This ensures that the PyDy model gives the same result as the Autolev
model."""

import os
import timeit

import numpy as np
import sympy as sm
import yaml
from algait2de.gait2de import evaluate_autolev_rhs as autolev_rhs
from pydy.codegen.ode_function_generators import generate_ode_function
from scipy.integrate import odeint

# local imports
from .. import derive, simulate, utils

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')


def test_accelerations():

    print('Deriving the symbolic dynamics model.')
    symbolics = derive.derive_equations_of_motion()
    print(symbolics)

    print('Generating the PyDy ODE numerical function.')
    pydy_rhs = generate_ode_function(
        symbolics.kanes_method.forcing,
        symbolics.coordinates,
        symbolics.speeds,
        constants=symbolics.constants,
        mass_matrix=symbolics.kanes_method.mass_matrix,
        coordinate_derivatives=sm.Matrix(symbolics.speeds),
        specifieds=symbolics.specifieds,
        generator='cython',
        # linear_sys_solve='sympy' is an undocumented option which calls
        # LUsolve on the cse'd Ax=b system before solving it. It takes a bit
        # longer to generate the function but results in fast evaluation.
        linear_sys_solver='sympy',
        constants_arg_type='array',
        specifieds_arg_type='array',
    )

    print('lambdifying the equations of motion')
    lam_M_F = sm.lambdify((symbolics.coordinates, symbolics.speeds,
                           symbolics.specifieds, symbolics.constants),
                          (symbolics.kanes_method.mass_matrix,
                           symbolics.kanes_method.forcing), cse=True)

    def lam_rhs(q, u, r, s):
        qdot = u
        M, F = lam_M_F(q, u, r, s)
        udot = np.linalg.solve(M, F.squeeze())
        return np.hstack((qdot, udot))

    coordinate_values = np.random.random(len(symbolics.coordinates))
    speed_values = np.random.random(len(symbolics.speeds))
    specified_values = np.random.random(len(symbolics.specifieds))

    constant_map = simulate.load_constants(
        symbolics.constants, os.path.join(ROOT, 'examples/example_constants.yml'))
    constant_values = np.array(list(constant_map.values()))

    args = (specified_values, constant_values)
    x = np.hstack((coordinate_values, speed_values))
    pydy_xdot = pydy_rhs(x, 0.0, *args)

    print('PyDy evaluation time:',
          timeit.timeit(lambda: pydy_rhs(x, 0.0, *args), number=1000))

    lam_xdot = lam_rhs(coordinate_values, speed_values, specified_values,
                       constant_values)

    print('Lambdify evaluation time:',
          timeit.timeit(lambda: lam_rhs(coordinate_values, speed_values,
                                        specified_values, constant_values),
                        number=1000))

    with open(os.path.join(ROOT, 'examples/example_constants.yml'), 'r') as f:
        constants_dict = yaml.load(f, Loader=yaml.SafeLoader)

    constants_dict = simulate.map_values_to_autolev_symbols(constants_dict)

    accelerations, grfs, animation_data = autolev_rhs(coordinate_values,
                                                      speed_values,
                                                      specified_values,
                                                      constants_dict)

    print('Autolev evaluation time:',
          timeit.timeit(lambda: autolev_rhs(coordinate_values, speed_values,
                                            specified_values, constants_dict),
                        number=1000))

    np.testing.assert_allclose(pydy_xdot[9:], accelerations)
    np.testing.assert_allclose(lam_xdot[9:], accelerations)


def test_with_control():
    symbolics = derive.derive_equations_of_motion(gait_cycle_control=True)

    pydy_rhs = generate_ode_function(
        symbolics.kanes_method.forcing_full,
        symbolics.coordinates,
        symbolics.speeds,
        constants=symbolics.constants,
        mass_matrix=symbolics.kanes_method.mass_matrix_full,
        specifieds=symbolics.specifieds,
        generator='cython',
        constants_arg_type='array',
        specifieds_arg_type='array',
    )

    coordinate_values = np.random.random(len(symbolics.coordinates))
    speed_values = np.random.random(len(symbolics.speeds))
    specified_values = np.random.random(len(symbolics.specifieds))

    constant_map = simulate.load_constants(
        symbolics.constants, os.path.join(ROOT, 'examples/example_constants.yml'))
    args = (specified_values, np.array(list(constant_map.values())))

    x = np.hstack((coordinate_values, speed_values))
    pydy_xdot = pydy_rhs(x, 0.0, *args)

    assert isinstance(pydy_xdot, np.ndarray)
    np.testing.assert_allclose(pydy_xdot[:9], speed_values)


def test_with_muscles(plot=False, animate=False):
    symbolics = derive.derive_equations_of_motion(include_muscles=True)
    print(symbolics)

    coordinate_values = np.random.random(len(symbolics.coordinates))
    muscle_values = np.random.random(len(symbolics.activations))
    speed_values = np.random.random(len(symbolics.speeds))
    specified_values = np.random.random(len(symbolics.specifieds))

    constant_map = simulate.load_constants(
        symbolics.constants, os.path.join(ROOT, 'examples/example_constants.yml'))
    constant_values = np.array(list(constant_map.values()))

    time_vector = np.linspace(0.0, 0.5, num=100)

    initial_conditions = np.zeros(len(symbolics.states))
    initial_conditions[1] = 1.0  # set hip above ground
    #initial_conditions[3] = np.deg2rad(120.0)  # right hip angle
    #initial_conditions[4] = -np.deg2rad(90.0)  # right knee angle
    #initial_conditions[5] = -np.deg2rad(25.0)  # right ankle angle
    #initial_conditions[6] = -np.deg2rad(40.0)  # left hip angle
    #initial_conditions[7] = -np.deg2rad(60.0)  # left knee angle
    #initial_conditions[8] = -np.deg2rad(25.0)  # left ankle angle

    if plot or animate:
        print('Generating the plot.')
        import matplotlib.pyplot as plt
        scene, fig, ax = utils.plot(symbolics, time_vector, initial_conditions,
                                    specified_values, constant_values)

    # TODO : Move the construction of M and F into Symbolics.
    num_simple = len(symbolics.coordinates) + len(symbolics.activations)
    num_dyn = len(symbolics.speeds)

    actd_zero_repl = dict(zip([a.diff() for a in symbolics.activations],
                              sm.zeros(len(symbolics.activations))))

    M = sm.BlockMatrix([[sm.eye(num_simple),
                         sm.zeros(num_simple, num_dyn)],
                        [sm.zeros(num_dyn, num_simple),
                         symbolics.kanes_method.mass_matrix]])
    F = sm.BlockMatrix([[sm.Matrix(symbolics.speeds)],
                        [-symbolics.mus_diff_eqs.xreplace(actd_zero_repl)],
                        [symbolics.kanes_method.forcing]])
    # NOTE: generate_ode_function fails if you leave these as block matrices
    # because iteration of a BlockMatrix returns nothing, unlike a normal SymPy
    # Matrix, so convert to normal matrices.
    M, F = sm.Matrix(M), sm.Matrix(F)

    print('Generating rhs function')
    pydy_rhs = generate_ode_function(
        F,
        symbolics.coordinates + symbolics.activations,
        symbolics.speeds,
        constants=symbolics.constants,
        mass_matrix=M,
        specifieds=symbolics.specifieds,
        generator='cython',
        # linear_sys_solve='sympy' is an undocumented option which calls
        # LUsolve on the cse'd Ax=b system before solving it. It takes a bit
        # longer to generate the function but results in fast evaluation.
        linear_sys_solver='sympy',
        constants_arg_type='array',
        specifieds_arg_type='array',
    )

    args = (specified_values, constant_values)

    x = np.hstack((coordinate_values, muscle_values, speed_values))
    pydy_xdot = pydy_rhs(x, 0.0, *args)

    assert isinstance(pydy_xdot, np.ndarray)
    np.testing.assert_allclose(pydy_xdot[:9], speed_values)

    args = (np.zeros(len(symbolics.specifieds)),
            np.array(list(constant_map.values())))

    trajectories = odeint(pydy_rhs, initial_conditions, time_vector, args=args)

    if plot and not animate:
        print('Showing the plot.')
        plt.show()
    elif animate:
        print('Generating the animation.')
        # NOTE : It is required to assign this to a variable if you want the
        # animation to run in the plot window.
        ani = utils.animate(scene, fig, time_vector, trajectories,
                            np.zeros((len(time_vector),
                                      len(symbolics.specifieds))),
                            np.array(list(constant_map.values())))
        print('Showing the animation.')
        plt.show()
