"""
Forward Simulation
==================

This example simply simulates and visualizes the uncontrolled motion and the
model falls down on the treadmill. It also compares the evaluation speed of
PyDy's and Autolev's models.
"""
import timeit

from algait2de.gait2de import evaluate_autolev_rhs
from pydy.codegen.ode_function_generators import generate_ode_function
from pygait2d import derive, simulate
from pygait2d import utils
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import sympy as sm
import yaml

# %%
# Derive the equations of motion, including a constant treadmill motion.
symbolics = derive.derive_equations_of_motion(treadmill=True)

# %%
# Load a parameter mapping from pygait2d symbol to numerical value, as well as
# a mappig of the symbol string to numerical value.
try:
    par_map = simulate.load_constants(symbolics.constants,
                                      'example_constants.yml')
    with open('example_constants.yml', 'r') as f:
        constants_dict = yaml.load(f, Loader=yaml.SafeLoader)
except FileNotFoundError:
    par_map = simulate.load_constants(symbolics.constants,
                                      'examples/example_constants.yml')
    with open('examples/example_constants.yml', 'r') as f:
        constants_dict = yaml.load(f, Loader=yaml.SafeLoader)

# %%
# Use PyDy to generate a function that can evaluate the right hand side of the
# ordinary differential equations of the multibody system. This uses PyDy's
# code geenration settings that result in the fastest numerical evaluation
# times at the cost of a slower code generation and compilation time.
rhs = generate_ode_function(
    symbolics.kanes_method.forcing,
    symbolics.coordinates,
    symbolics.speeds,
    constants=list(par_map.keys()),
    mass_matrix=symbolics.kanes_method.mass_matrix,
    coordinate_derivatives=sm.Matrix(symbolics.speeds),
    specifieds=symbolics.specifieds,
    generator='cython',
    constants_arg_type='array',
    specifieds_arg_type='array',
    linear_sys_solver='sympy',  # slowest code generation, fastest evaluation
)

# %%
# Prepare numerical arrays to be passed to the ODE functions.
specifieds_vals = np.zeros(len(symbolics.specifieds))
specifieds_vals[-1] = 1.0

args = (specifieds_vals, np.array(list(par_map.values())))

initial_conditions = np.zeros(len(symbolics.states))
initial_conditions[1] = 1.0  # set hip above ground
initial_conditions[3] = np.deg2rad(5.0)  # right hip angle
initial_conditions[6] = -np.deg2rad(5.0)  # left hip angle

# %%
# Time the average execuation of PyDy's ODE function evaluation.
print('PyDy evaluation time:',
      timeit.timeit(lambda: rhs(initial_conditions, 0.0, *args), number=1000))

# %%
# Time the average execuation of Autolev s ODE function evaluation.
autolev_constants_dict = simulate.map_values_to_autolev_symbols(constants_dict)
print('Autolev evaluation time:',
      timeit.timeit(lambda: evaluate_autolev_rhs(initial_conditions[:9],
                                                 initial_conditions[9:],
                                                 specifieds_vals,
                                                 autolev_constants_dict),
                    number=1000))

# %%
# Simulation the model for two seconds using the LSODA integrator (switches #
# between stiff and non-stiff modes).
time_vector = np.linspace(0.0, 2.0, num=201)
trajectories = odeint(rhs, initial_conditions, time_vector, args=args)

# %%
# Visualization the resulting forward simulation motion.
scene, fig, ax = utils.plot(symbolics, time_vector, initial_conditions,
                            args[0], args[1])
ax.set_xlim((-0.8, 0.8))
ax.set_ylim((-0.2, 1.4))
ani = utils.animate(scene, fig, time_vector, trajectories,
                    np.zeros((len(time_vector), len(symbolics.specifieds))),
                    args[1])

plt.show()
