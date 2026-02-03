"""
Forward Simulation
==================

This example simply simulates and visualizes the uncontrolled motion and the
model "falls down".
"""
from pygait2d import derive, simulate
from pygait2d import utils
import numpy as np
from pydy.codegen.ode_function_generators import generate_ode_function
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# %%
symbolics = derive.derive_equations_of_motion(treadmill=True)

# %%
constant_values = simulate.load_constants(symbolics.constants,
                                          'example_constants.yml')

# %%
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

# %%
specifieds_vals = np.zeros(len(symbolics.specifieds))
specifieds_vals[-1] = 1.0

args = (specifieds_vals, np.array(list(constant_values.values())))

time_vector = np.linspace(0.0, 2.0, num=1000)
initial_conditions = np.zeros(len(symbolics.states))
initial_conditions[1] = 1.0  # set hip above ground
initial_conditions[3] = np.deg2rad(5.0)  # right hip angle
initial_conditions[6] = -np.deg2rad(5.0)  # left hip angle
trajectories = odeint(rhs, initial_conditions, time_vector, args=args)

# %%
scene, fig, ax = utils.plot(symbolics, time_vector, initial_conditions,
                            args[0], args[1])
ax.set_xlim((-0.8, 0.8))
ax.set_ylim((-0.2, 1.4))
ani = utils.animate(scene, fig, time_vector, trajectories,
                    np.zeros((len(time_vector), len(symbolics.specifieds))),
                    args[1])

plt.show()
