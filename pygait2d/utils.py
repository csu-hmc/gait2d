#!/usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib.animation import FuncAnimation
from symmeplot.matplotlib import Scene3D
from sympy import srepr, Matrix, sympify
import matplotlib.pyplot as plt
import numpy as np
import sympy as sm
import sympy.physics.mechanics as me


def save_sympy_matrix(matrix, filename):
    """Writes a matrix to file in the SymPy representation (srepr)."""
    num_rows, num_cols = matrix.shape
    with open(filename, 'w') as f:
        f.write(str(num_rows) + "\n")
        f.write(str(num_cols) + "\n")
        for expr in matrix:
            f.write(srepr(expr) + "\n")


def load_sympy_matrix(filename):
    """Loads a matrix from file created with save_sympy_matrix."""
    exprs = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                num_rows = int(line.strip())
            elif i == 1:
                num_cols = int(line.strip())
            else:
                exprs.append(sympify(line.strip()))
    return Matrix(exprs).reshape(num_rows, num_cols)


class ExtensorPathway(me.PathwayBase):
    def __init__(self, origin, insertion, axis_point, axis, parent_axis,
                 child_axis, radius, coordinate):
        """A custom pathway that wraps a circular arc around a pin joint.  This
        is intended to be used for extensor muscles. For example, a triceps
        wrapping around the elbow joint to extend the upper arm at the elbow.

        Parameters
        ==========
        origin : Point
            Muscle origin point fixed on the parent body (A).
        insertion : Point
            Muscle insertion point fixed on the child body (B).
        axis_point : Point
            Pin joint location fixed in both the parent and child.
        axis : Vector
            Pin joint rotation axis.
        parent_axis : Vector
            Axis fixed in the parent frame (A) that is directed from the pin
            joint point to the muscle origin point.
        child_axis : Vector
            Axis fixed in the child frame (B) that is directed from the pin
            joint point to the muscle insertion point.
        radius : sympyfiable
            Radius of the arc that the muscle wraps around.
        coordinate : sympfiable function of time
            Joint angle, zero when parent and child frames align. Positive
            rotation about the pin joint axis, B with respect to A.

        Notes
        =====
        Only valid for coordinate >= 0.

        """
        super().__init__(origin, insertion)
        self.origin = origin
        self.insertion = insertion
        self.axis_point = axis_point
        self.axis = axis.normalize()
        self.parent_axis = parent_axis.normalize()
        self.child_axis = child_axis.normalize()
        self.radius = radius
        self.coordinate = coordinate
        self.origin_distance = axis_point.pos_from(origin).magnitude()
        self.insertion_distance = axis_point.pos_from(insertion).magnitude()
        self.origin_angle = sm.asin(self.radius/self.origin_distance)
        self.insertion_angle = sm.asin(self.radius/self.insertion_distance)

    @property
    def length(self):
        """Length of the pathway.
        Length of two fixed length line segments and a changing arc length
        of a circle.
        """
        angle = self.origin_angle + self.coordinate + self.insertion_angle
        arc_length = self.radius*angle
        origin_segment_length = self.origin_distance*sm.cos(self.origin_angle)
        insertion_segment_length = self.insertion_distance*sm.cos(
            self.insertion_angle)
        return origin_segment_length + arc_length + insertion_segment_length

    @property
    def extension_velocity(self):
        """Extension velocity of the pathway.
        Arc length of circle is the only thing that changes when the elbow
        flexes and extends.
        """
        return self.radius*self.coordinate.diff(me.dynamicsymbols._t)

    @property
    def plot_points(self):
        return (self.origin, self.parent_tangency_point, self.axis_point,
                self.child_tangency_point, self.insertion)

    def to_loads(self, force_magnitude):
        """Loads in the correct format to be supplied to `KanesMethod`.
        Forces applied to origin, insertion, and P from the muscle wrapped
        over circular arc of radius r.
        """
        # TODO: generate these points on init
        self.parent_tangency_point = me.Point('Aw')  # fixed in parent
        self.child_tangency_point = me.Point('Bw')  # fixed in child
        self.parent_tangency_point.set_pos(
            self.axis_point,
            self.radius*sm.cos(self.origin_angle)*self.parent_axis.cross(
                self.axis)
            + self.radius*sm.sin(self.origin_angle)*self.parent_axis,
        )
        self.child_tangency_point.set_pos(
            self.axis_point,
            -self.radius*sm.cos(self.insertion_angle)*self.child_axis.cross(
                self.axis)
            + self.radius*sm.sin(self.insertion_angle)*self.child_axis),
        parent_force_direction_vector = self.origin.pos_from(
            self.parent_tangency_point)
        child_force_direction_vector = self.insertion.pos_from(
            self.child_tangency_point)
        force_on_parent = (force_magnitude*
                           parent_force_direction_vector.normalize())
        force_on_child = (force_magnitude*
                          child_force_direction_vector.normalize())
        loads = [
            me.Force(self.origin, force_on_parent),
            me.Force(self.axis_point, -(force_on_parent + force_on_child)),
            me.Force(self.insertion, force_on_child),
        ]
        return loads


def generate_animation(symmod, times, xs, rs, ps):

    ground = symmod.inertial_frame
    origin = symmod.origin
    trunk, rthigh, rshank, rfoot, lthigh, lshank, lfoot = symmod.segments

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

    # hip_proj = origin.locatenew('m', qax*ground.x)
    # scene = Scene3D(ground, hip_proj, ax=ax3d)
    scene = Scene3D(ground, origin, ax=ax)

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

    if symmod.muscles is not None:
        for i, mus in enumerate(symmod.muscles):
            color = "C{}".format(i)
            try:
                scene.add_line(mus.pathway.plot_points, color=color,
                               marker='.')
                print(mus.pathway.plot_points, color)
            except AttributeError:
                scene.add_line(mus.pathway.attachments, color=color,
                               marker='.')
                print(mus.pathway.attachments, color)

    # creates a moving ground (many points to deal with matplotlib limitation)
    # ?? can we make the dashed line move to the left?
    scene.add_line([origin.locatenew('gl', s*ground.x) for s in
                    np.linspace(-2.0, 2.0)], linestyle='--',
                    color='tab:green', axlim_clip=True)

    # adds CoM and unit vectors for each body segment
    for seg in symmod.segments:
        scene.add_body(seg.rigid_body)

    # show ground reaction force vectors at the heels and toes, scaled to
    # visually reasonable length
    #scene.add_vector(contact_force(rfoot.toe, ground, origin, v)/600.0,
                        #rfoot.toe, color="tab:blue")
    #scene.add_vector(contact_force(rfoot.heel, ground, origin, v)/600.0,
                        #rfoot.heel, color="tab:blue")
    #scene.add_vector(contact_force(lfoot.toe, ground, origin, v)/600.0,
                        #lfoot.toe, color="tab:blue")
    #scene.add_vector(contact_force(lfoot.heel, ground, origin, v)/600.0,
                        #lfoot.heel, color="tab:blue")

    scene.lambdify_system(symmod.states[:] + symmod.specifieds[:] +
                          symmod.constants[:])
    gait_cycle = np.vstack((
        xs.T,  # q, u shape(2n, N)
        #np.zeros((3, len(times))),  # Fax, Fay, Ta (hand of god), shape(3, N)
        rs.T,  # r, shape(q, N)
        np.repeat(np.atleast_2d(ps).T, len(times), axis=1),  # p, shape(r, N)
    ))
    scene.evaluate_system(*gait_cycle[:, 0])

    scene.axes.set_proj_type("ortho")
    scene.axes.view_init(90, -90, 0)
    scene.plot()

    ax.set_xlim((-0.8, 0.8))
    ax.set_ylim((-0.2, 1.4))
    ax.set_aspect('equal')

    # TODO : return plot here so a single frame can be viewed
    # plt.show()

    def update(i):
        scene.evaluate_system(*gait_cycle[:, i])
        scene.update()
        return scene.artists

    ani = FuncAnimation(
        fig,
        update,
        frames=range(len(times)),
        interval=(times[1] - times[0])*1000,
    )

    #ani.save('muscles.gif', fps=int(1.0/(times[1] - times[0])))

    return ani
