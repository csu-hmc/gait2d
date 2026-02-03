"""
Joint Torque Gait Predictive Simulation
=======================================

This example replicates a similar predictive simulation solution as shown in
[Ackermann2010]_ using joint torques as inputs instead of muscle activations
[1]_.

The optimal control goal is to find the open-loop joint torque trajectories for
hip, knee, and ankle torques that generate a minimal mean-torque periodic
motion to ambulate at a specified average speed.

.. note::

   Requires `opty >= 1.5.0`::

      conda install -c conda-forge opty>=1.50

Import all necessary modules, functions, and classes:
"""
import numpy as np
import sympy as sm
from opty import Problem
import matplotlib.pyplot as plt
from pygait2d.derive import derive_equations_of_motion
from pygait2d.segment import time_symbol
from pygait2d.simulate import load_constants
from pygait2d.utils import plot, animate

# %%
# Derive the equations of motion using gait2d. Ground reaction forces are only
# needed on the feet points and no external forces and torques will act on the
# trunk.
symbolics = derive_equations_of_motion(
    prevent_ground_penetration=False,
    hand_of_god=False,
)

# %%
# The equations of motion have this many mathematical operations:
eom = symbolics.equations_of_motion
sm.count_ops(eom)

# %%
# :math:`t_f - t_0` needs to be available to compute the average speed in the
# instance constraint, so add an extra differential equation that is the time
# derivative of the difference in time.
#
# .. math::
#
#    \Delta_t(t) = \int_{t_0}^{t} d\tau
#
delt = sm.Function('delt', real=True)(time_symbol)
eom = eom.col_join(sm.Matrix([delt.diff(time_symbol) - 1]))

states = symbolics.states + [delt]
num_states = len(states)

# %%
# The generalized coordinates are the hip lateral position :math:`q_{ax}` and
# veritcal position :math:`q_{ay}`, the trunk angle with respect to vertical
# :math:`q_a` and the relative joint angles:
#
# - right: hip (b), knee (c), ankle (d)
# - left: hip (e), knee (f), ankle (g)
#
# Each joint has a joint torque acting between the adjacent bodies.
qax, qay, qa, qb, qc, qd, qe, qf, qg = symbolics.coordinates
uax, uay, ua, ub, uc, ud, ue, uf, ug = symbolics.speeds
Tb, Tc, Td, Te, Tf, Tg = symbolics.specifieds

# %%
# The constants are loaded from a file of realistic geometry, mass, inertia,
# and foot deformation properties of an adult human.
try:
    par_map = load_constants(symbolics.constants, 'example_constants.yml')
except FileNotFoundError:
    par_map = load_constants(symbolics.constants,
                             'examples/example_constants.yml')
par_map

# %%
# First solve a simple two-node standing solution to generate an initial guess
# for a walking solution. Set the average ambulation speed to zero and the
# number of discretization nodes for the half period to 2 and define the time
# step as a variable :math:`h`.
h = sm.symbols('h', real=True, positive=True)

static_num_nodes = 2
duration = (static_num_nodes - 1)*h

speed = sm.symbols('v', real=True)
par_map[speed] = 0.0

# %%
# Bound all the states to human realizable ranges.
#
# - The trunk should stay generally upright and be at a possible walking
#   height.
# - Only let the hip, knee, and ankle flex and extend to realistic limits.
# - Put a maximum on the peak torque values.
bounds = {
    h: (0.005, 0.1),
    delt: (0.0, 10.0),
    qax: (0.0, 10.0),
    qay: (0.5, 1.5),
    qa: np.deg2rad((-60.0, 60.0)),
    uax: (0.0, 10.0),
    uay: (-10.0, 10.0),
}
# hip
bounds.update({k: (-np.deg2rad(40.0), np.deg2rad(40.0))
               for k in [qb, qe]})
# knee
bounds.update({k: (-np.deg2rad(60.0), 0.0)
               for k in [qc, qf]})
# foot
bounds.update({k: (-np.deg2rad(30.0), np.deg2rad(30.0))
               for k in [qd, qg]})
# all rotational speeds
bounds.update({k: (-np.deg2rad(400.0), np.deg2rad(400.0))
               for k in [ua, ub, uc, ud, ue, uf, ug]})
# all joint torques
bounds.update({k: (-60.0, 60.0)
               for k in [Tb, Tc, Td, Te, Tf, Tg]})

# %%
# The average speed can be fixed by constraining the total distance traveled.
# To enforce a half period, set the right leg's angles at the initial time to
# be equal to the left leg's angles at the final time and vice versa. The same
# goes for the joint angular rates.
instance_constraints = (
    delt.func(0*h) - 0.0,
    qax.func(0*h) - 0.0,
    qax.func(duration) - speed*delt.func(duration),
    qay.func(0*h) - qay.func(duration),
    qa.func(0*h) - qa.func(duration),
    qb.func(0*h) - qe.func(duration),
    qc.func(0*h) - qf.func(duration),
    qd.func(0*h) - qg.func(duration),
    qe.func(0*h) - qb.func(duration),
    qf.func(0*h) - qc.func(duration),
    qg.func(0*h) - qd.func(duration),
    uax.func(0*h) - uax.func(duration),
    uay.func(0*h) - uay.func(duration),
    ua.func(0*h) - ua.func(duration),
    ub.func(0*h) - ue.func(duration),
    uc.func(0*h) - uf.func(duration),
    ud.func(0*h) - ug.func(duration),
    ue.func(0*h) - ub.func(duration),
    uf.func(0*h) - uc.func(duration),
    ug.func(0*h) - ud.func(duration),
)


# %%
# The objective is to minimize the mean of all joint torques.
def obj(prob, free):
    """Minimize the sum of the squares of the control torques."""
    interval = prob.extract_values(free, h)[0]
    torques = prob.extract_values(free, *symbolics.specifieds)
    return interval*np.sum(torques**2)


def obj_grad(prob, free):
    interval = prob.extract_values(free, h)[0]
    torques = prob.extract_values(free, *symbolics.specifieds)
    grad = np.zeros_like(free)
    prob.fill_free(grad, 2.0*interval*torques, *symbolics.specifieds)
    prob.fill_free(grad, np.sum(torques**2), h)
    return grad


# %%
# Create an optimization problem and solve it.
prob = Problem(
    obj,
    obj_grad,
    eom,
    states,
    static_num_nodes,
    h,
    known_parameter_map=par_map,
    instance_constraints=instance_constraints,
    bounds=bounds,
    time_symbol=time_symbol,
    parallel=True,
    tmp_dir='codegen',
)

# %%
# Find the optimal standing solution.
initial_guess = np.zeros(prob.num_free)
initial_guess[-1] = 0.01
solution, info = prob.solve(initial_guess)

# %%
# Repeat the standing solution for the number of nodes used in the walking
# solution and add some noise to the initial guess.
num_nodes = 50
xs, rs, _, h_val = prob.parse_free(solution)
x_rep = np.repeat(xs[:, 0:1], num_nodes, axis=1)
r_rep = np.repeat(rs[:, 0:1], num_nodes, axis=1)
initial_guess = np.hstack((x_rep.flatten(), r_rep.flatten(), 0.01))
initial_guess += np.random.normal(0.0, 0.01, size=initial_guess.shape)

# TODO : It would be nice if we could update the number of nodes and the
# instance constraints without creating a new problem.
duration = (num_nodes - 1)*h
instance_constraints = (
    delt.func(0*h) - 0.0,
    qax.func(0*h) - 0.0,
    qax.func(duration) - speed*delt.func(duration),
    qay.func(0*h) - qay.func(duration),
    qa.func(0*h) - qa.func(duration),
    qb.func(0*h) - qe.func(duration),
    qc.func(0*h) - qf.func(duration),
    qd.func(0*h) - qg.func(duration),
    qe.func(0*h) - qb.func(duration),
    qf.func(0*h) - qc.func(duration),
    qg.func(0*h) - qd.func(duration),
    uax.func(0*h) - uax.func(duration),
    uay.func(0*h) - uay.func(duration),
    ua.func(0*h) - ua.func(duration),
    ub.func(0*h) - ue.func(duration),
    uc.func(0*h) - uf.func(duration),
    ud.func(0*h) - ug.func(duration),
    ue.func(0*h) - ub.func(duration),
    uf.func(0*h) - uc.func(duration),
    ug.func(0*h) - ud.func(duration),
)

prob = Problem(
    obj,
    obj_grad,
    eom,
    states,
    num_nodes,
    h,
    known_parameter_map=par_map,
    instance_constraints=instance_constraints,
    bounds=bounds,
    time_symbol=time_symbol,
    parallel=True,
    tmp_dir='codegen',
)

# %%
# Solve the optimization problem for increasing walking speeds using the
# standing solution as the first initial guess and using the prior solution as
# the subsequent guesses.
solution = initial_guess
for new_speed in np.linspace(0.01, 1.3, num=8):
    par_map[speed] = new_speed
    print(f"Running optimization for walking speed: "
          f"{prob.collocator.known_parameter_map[speed]}")
    solution, info = prob.solve(solution)

# %%
# Extract the solution trajectories.
xs, rs, _, h_val = prob.parse_free(solution)
ps = np.array(list(par_map.values()))
times = prob.time_vector(solution=solution)

# TODO : plot() and animate() can only work with states and parameters defined
# in Symbolics, so we have to delete these extras added above. It would be
# better if plot() and animate() were not affected by this.
xs = xs[:-1, :]  # drop the extra time state
del par_map[speed]
ps = np.array(list(par_map.values()))

# %%
# Animate the motion.
scene, fig, ax = plot(symbolics, times, xs[:, 0], rs[:, 0], ps)
ax.set_xlim((-0.8, 0.8))
ax.set_ylim((-0.2, 1.4))
ani = animate(scene, fig, times, xs.T, rs.T, ps)

plt.show()
