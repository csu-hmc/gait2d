r"""
Joint Torque Predictive Simulation
==================================

.. note::

   Requires `opty >= 1.5.0`::

      conda install -c conda-forge opty>=1.50

Introduction
------------

This example replicates a similar predictive simulation solution as shown in
[Ackermann2010]_ using joint torques as inputs instead of muscle activations
[1]_.

The optimal control goal is to find the open-loop joint torque trajectories for
hip, knee, and ankle torques that generate a minimal mean-torque periodic
motion to ambulate at a specified average speed.

Import all necessary modules, functions, and classes:
"""
import numpy as np
import sympy as sm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from opty import Problem
from symmeplot.matplotlib import Scene3D

from pygait2d import derive, simulate
from pygait2d.segment import time_symbol, contact_force

# %%
# Derive the equations of motion using gait2d. Ground reaction forces are only
# # needed on the feet points and no external forces and torques will act on
# the trunk.
symbolics = derive.derive_equations_of_motion(
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
par_map = simulate.load_constants(symbolics.constants, 'example_constants.yml')
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
bounds.update({k: (-600.0, 600.0)
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
    torques = prob.extract_values(free, *symbolics.specifieds)
    interval = prob.extract_values(free, h)[0]
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
# Use symmeplot to make an animation of the motion.
xs, rs, _, h_val = prob.parse_free(solution)
times = prob.time_vector(solution=solution)
# TODO : Switch to pygait2d.utils.plot() and pygait2d.utils.animate()


def animate():

    origin = symbolics.origin
    ground = symbolics.inertial_frame
    trunk, rthigh, rshank, rfoot, lthigh, lshank, lfoot = symbolics.segments

    fig = plt.figure(figsize=(10.0, 4.0))

    ax3d = fig.add_subplot(1, 2, 1, projection='3d')
    ax2d = fig.add_subplot(1, 2, 2)

    hip_proj = origin.locatenew('m', qax*ground.x)
    scene = Scene3D(ground, hip_proj, ax=ax3d)

    # creates the stick person
    scene.add_line([
        rshank.joint,
        rfoot.toe,
        rfoot.heel,
        rshank.joint,
        rthigh.joint,
        trunk.joint,
        trunk.mass_center,
        trunk.joint,
        lthigh.joint,
        lshank.joint,
        lfoot.heel,
        lfoot.toe,
        lshank.joint,
    ], color="k")

    # creates a moving ground (many points to deal with matplotlib limitation)
    scene.add_line([origin.locatenew('gl', s*ground.x) for s in
                    np.linspace(-2.0, 2.0)], linestyle='--', color='tab:green',
                   axlim_clip=True)

    # adds CoM and unit vectors for each body segment
    for seg in symbolics.segments:
        scene.add_body(seg.rigid_body)

    # show ground reaction force vectors at the heels and toes, scaled to
    # visually reasonable length
    scene.add_vector(contact_force(rfoot.toe, ground, origin)/600.0,
                     rfoot.toe, color="tab:blue")
    scene.add_vector(contact_force(rfoot.heel, ground, origin)/600.0,
                     rfoot.heel, color="tab:blue")
    scene.add_vector(contact_force(lfoot.toe, ground, origin)/600.0,
                     lfoot.toe, color="tab:blue")
    scene.add_vector(contact_force(lfoot.heel, ground, origin)/600.0,
                     lfoot.heel, color="tab:blue")

    del par_map[speed]
    scene.lambdify_system(states + symbolics.specifieds + symbolics.constants)
    gait_cycle = np.vstack((
        xs,  # q, u shape(2n, N)
        rs,  # r, shape(q, N)
        np.repeat(np.atleast_2d(np.array(list(par_map.values()))).T,
                  len(times), axis=1),  # p, shape(r, N)
    ))
    scene.evaluate_system(*gait_cycle[:, 0])

    scene.axes.set_proj_type("ortho")
    scene.axes.view_init(90, -90, 0)
    scene.plot(prettify=False)

    ax3d.set_xlim((-0.8, 0.8))
    ax3d.set_ylim((-0.2, 1.4))
    ax3d.set_aspect('equal')
    for axis in (ax3d.xaxis, ax3d.yaxis, ax3d.zaxis):
        axis.set_ticklabels([])
        axis.set_ticks_position("none")

    eval_rforce = sm.lambdify(
        states + symbolics.specifieds + symbolics.constants,
        (contact_force(rfoot.toe, ground, origin) +
         contact_force(rfoot.heel, ground, origin)).to_matrix(ground),
        cse=True)

    eval_lforce = sm.lambdify(
        states + symbolics.specifieds + symbolics.constants,
        (contact_force(lfoot.toe, ground, origin) +
         contact_force(lfoot.heel, ground, origin)).to_matrix(ground),
        cse=True)

    rforces = np.array([eval_rforce(*gci).squeeze() for gci in gait_cycle.T])
    lforces = np.array([eval_lforce(*gci).squeeze() for gci in gait_cycle.T])

    ax2d.plot(times, rforces[:, :2], times, lforces[:, :2])
    ax2d.grid()
    ax2d.set_ylabel('Force [N]')
    ax2d.set_xlabel('Time [s]')
    ax2d.legend(['Horizontal GRF (r)', 'Vertical GRF (r)',
                 'Horizontal GRF (l)', 'Vertical GRF (l)'], loc='upper right')
    ax2d.set_title('Foot Ground Reaction Force Components')
    vline = ax2d.axvline(times[0], color='black')

    def update(i):
        scene.evaluate_system(*gait_cycle[:, i])
        scene.update()
        vline.set_xdata([times[i], times[i]])
        return scene.artists + (vline,)

    ani = FuncAnimation(
        fig,
        update,
        frames=range(len(times)),
        interval=h_val*1000,
    )

    return ani


animation = animate()


plt.show()
