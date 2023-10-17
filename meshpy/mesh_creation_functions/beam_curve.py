# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# MeshPy: A beam finite element input generator
#
# MIT License
#
# Copyright (c) 2018-2023
#     Ivo Steinbrecher
#     Institute for Mathematics and Computer-Based Simulation
#     Universitaet der Bundeswehr Muenchen
#     https://www.unibw.de/imcs-en
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# -----------------------------------------------------------------------------
"""
This file has functions to create a beam from a parametric curve.
"""

# Python packages.
import numpy as np

# Meshpy modules.
from ..conf import mpy
from ..rotation import Rotation, smallest_rotation
from .beam_generic import create_beam_mesh_function


def create_beam_mesh_curve(
    mesh,
    beam_object,
    material,
    function,
    interval,
    *,
    function_derivative=None,
    function_rotation=None,
    **kwargs
):
    """
    Generate a beam from a parametric curve. Integration along the beam is
    performed with scipy, and if the gradient is not explicitly provided, it is
    calculated with the numpy wrapper autograd.

    Args
    ----
    mesh: Mesh
        Mesh that the curve will be added to.
    beam_object: Beam
        Class of beam that will be used for this line.
    material: Material
        Material for this line.
    function: function
        3D-parametric curve that represents the beam axis. If only a 2D
        point is returned, the triad creation is simplified. If
        mathematical functions are used, they have to come from the wrapper
        autograd.numpy.
    interval: [start end]
        Start and end values for the parameter of the curve.
    function_derivative: function -> R3
        Explicitly provide the jacobian of the centerline position.
    function_rotation: function -> Rotation
        If this argument is given, the triads are computed with this
        function, on the same interval as the position function. Must
        return a Rotation object.
        If no function_rotation is given, the rotation of the first node
        is calculated automatically and all subsequent nodal rotations
        are calculated based on a smallest rotation mapping onto the curve
        tangent vector.

    **kwargs (for all of them look into create_beam_mesh_function)
    ----
    n_el: int
        Number of equally spaced beam elements along the line. Defaults to 1.
        Mutually exclusive with l_el.
    l_el: float
        Desired length of beam elements. This requires the option interval_length
        to be set. Mutually exclusive with n_el. Be aware, that this length
        might not be achieved, if the elements are warped after they are
        created.

    Return
    ----
    return_set: GeometryName
        Set with the 'start' and 'end' node of the curve. Also a 'line' set
        with all nodes of the curve.
    """

    print("start")

    # Packages for AD and numerical integration.
    from autograd import jacobian
    import autograd.numpy as npAD
    import scipy.integrate as integrate
    import scipy.optimize as optimize

    # Check size of position function
    if len(function(interval[0])) == 2:
        is_3d_curve = False
    elif len(function(interval[0])) == 3:
        is_3d_curve = True
    else:
        raise ValueError("Function must return either 2d or 3d curve!")

    # Check rotation function.
    if function_rotation is None:
        is_rot_funct = False
    else:
        is_rot_funct = True

    # Check that the position is an np.array
    if not isinstance(function(float(interval[0])), np.ndarray):
        raise TypeError(
            "Function must be of type np.ndarray, got {}!".format(
                type(function(float(interval[0])))
            )
        )

    # Get the derivative of the position function and the increment along
    # the curve.
    if function_derivative is None:
        rp = jacobian(function)
    else:
        rp = function_derivative

    # Check which one of the boundaries is larger.
    if interval[0] > interval[1]:
        # In this case rp needs to be negated.
        rp_positive = rp
        rp = lambda t: -(rp_positive(t))

    def ds(t):
        """Increment along the curve."""
        return npAD.linalg.norm(rp(t))

    def S(t, start_t=None, start_S=None):
        """
        Function that integrates the length until a parameter value.
        A speedup can be achieved by giving start_t and start_S, the
        parameter and Length at a known point.
        """
        if start_t is None and start_S is None:
            st = interval[0]
            sS = 0
        elif start_t is not None and start_S is not None:
            st = start_t
            sS = start_S
        else:
            raise ValueError("Input parameters are wrong!")
        return integrate.quad(ds, st, t)[0] + sS

    def get_t_along_curve(arc_length, t0, **kwargs):
        """
        Calculate the parameter t where the length along the curve is
        arc_length. t0 is the start point for the Newton iteration.
        """
        t_root = optimize.newton(
            lambda t: S(t, **kwargs) - arc_length, t0, fprime=ds, full_output=True
        )
        print(t_root[1])
        return t_root[0]

    class BeamFunctions:
        """
        This object manages the creation of functions which are used to create the
        beam nodes. By wrapping this in this object, we can store some data and speed
        up the numerical integration along the curve.
        """

        def __init__(self):
            """
            Initialize the object.
            """
            self.t_start_element = interval[0]
            self.last_triad = None

            if is_3d_curve:
                r_prime = rp(float(interval[0]))
                if abs(np.dot(r_prime, [0, 0, 1])) < abs(np.dot(r_prime, [0, 1, 0])):
                    t2_temp = [0, 0, 1]
                else:
                    t2_temp = [0, 1, 0]
                self.last_triad = Rotation.from_basis(r_prime, t2_temp)

        def __call__(self, length_a, length_b):
            """
            This object is called with the interval limits. This method returns a
            function that evaluates the position and rotation within this interval.
            """

            # Length of the beam element in physical space.
            L = length_b - length_a

            def get_beam_position_and_rotation_at_xi(xi):
                """
                Evaluate the beam position and rotation at xi. xi is the beam element
                parameter coordinate, i.e., xi = [-1, 1].
                """
                # Parameter for xi.
                t = get_t_along_curve(
                    length_a + 0.5 * (xi + 1) * L,
                    self.t_start_element,
                    start_t=self.t_start_element,
                    start_S=length_a,
                )

                # Position at xi.
                if is_3d_curve:
                    pos = function(t)
                else:
                    pos = np.zeros(3)
                    pos[:2] = function(t)

                # Rotation at xi.
                if is_rot_funct:
                    rot = function_rotation(t)
                else:
                    r_prime = rp(t)
                    if is_3d_curve:
                        # Create the next triad via the smallest rotation mapping based
                        # on the last triad.
                        rot = smallest_rotation(self.last_triad, r_prime)
                        self.last_triad = rot.copy()
                    else:
                        # The rotation simplifies in the 2d case.
                        rot = Rotation([0, 0, 1], np.arctan2(r_prime[1], r_prime[0]))

                if np.abs(xi - 1) < mpy.eps_pos:
                    # Set start values for the next element.
                    self.t_start_element = t

                # Return the needed values for beam creation.
                return (pos, rot)

            return get_beam_position_and_rotation_at_xi

    # Now create the beam.
    # Get the length of the whole segment.
    length = S(interval[1])

    # print(integrate.quad(ds, interval[0], interval[1]))

    # n = 10
    # full_integration_output = integrate.quad(
    #     ds,
    #     interval[0],
    #     interval[1],
    #     points=np.linspace(interval[0], interval[1], n),
    #     limit=n,
    #     full_output=True,
    # )
    # points = full_integration_output[2]["alist"]
    # print(points)
    # # points.append(interval[1])
    # interval_lengths = full_integration_output[2]["rlist"]
    # accumulated_lengths = interval_lengths
    # for i in range(1, len(accumulated_lengths)):
    #     accumulated_lengths[i] += accumulated_lengths[i - 1]

    # Create the beam in the mesh
    return create_beam_mesh_function(
        mesh,
        beam_object=beam_object,
        material=material,
        function_generator=BeamFunctions(),
        interval=[0.0, length],
        interval_length=length,
        **kwargs,
    )
