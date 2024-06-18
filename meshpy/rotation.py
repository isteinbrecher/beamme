# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# MeshPy: A beam finite element input generator
#
# MIT License
#
# Copyright (c) 2018-2024
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
This module defines a class that represents a rotation in 3D.
"""

# Python modules.
import numpy as np

# Meshpy modules.
from . import mpy


def skew_matrix(vector):
    """Return the skew matrix for the vector"""
    skew = np.zeros([3, 3])
    skew[0, 1] = -vector[2]
    skew[0, 2] = vector[1]
    skew[1, 0] = vector[2]
    skew[1, 2] = -vector[0]
    skew[2, 0] = -vector[1]
    skew[2, 1] = vector[0]
    return skew


class Rotation:
    """
    A class that represents a rotation of a coordinate system.
    Internally the rotations are stored as quaternions.
    """

    def __init__(self, *args):
        """
        Initialize the rotation object.

        Args
        ----
        *args:
            - Rotation()
                Create a identity rotation.
            - Rotation(axis, phi)
                Create a rotation around the vector axis with the angle phi.
            - Rotation([q0, q1, q2, q3])
                Create a rotation with the quaternion values q0...q3.
        """

        self.q = np.zeros(4)

        if len(args) == 0:
            # Identity element.
            self.q[0] = 1
        elif len(args) == 1 and len(args[0]) == 4:
            # Set directly from quaternion
            # To avoid error accumulation, normalize the quaternion here
            q = np.array(args[0])
            self.q[:] = q / np.linalg.norm(q)
        elif len(args) == 2:
            # Set from vector and rotation angle.
            vector = args[0]
            phi = args[1]
            norm = np.linalg.norm(vector)
            if np.abs(phi) < mpy.eps_quaternion:
                self.q[0] = 1
            elif norm < mpy.eps_quaternion:
                raise ValueError("The rotation axis can not be a zero vector!")
            else:
                self.q[0] = np.cos(0.5 * phi)
                self.q[1:] = np.sin(0.5 * phi) * np.array(vector) / norm
        else:
            raise ValueError(f"The given arguments {args} are invalid!")

    @classmethod
    def from_rotation_matrix(cls, R):
        """
        Create the object from a rotation matrix.
        The code is based on Spurriers algorithm:
            R. A. Spurrier (1978): “Comment on “Singularity-free extraction of a quaternion from a
            direction-cosine matrix”
        """

        q = np.zeros(4)
        trace = np.trace(R)
        values = [R[i, i] for i in range(3)]
        values.append(trace)
        arg_max = np.argmax(values)
        if arg_max == 3:
            q[0] = np.sqrt(trace + 1) * 0.5
            q[1] = (R[2, 1] - R[1, 2]) / (4 * q[0])
            q[2] = (R[0, 2] - R[2, 0]) / (4 * q[0])
            q[3] = (R[1, 0] - R[0, 1]) / (4 * q[0])
        else:
            i_index = arg_max
            j_index = (i_index + 1) % 3
            k_index = (i_index + 2) % 3
            q_i = np.sqrt(R[i_index, i_index] * 0.5 + (1 - trace) * 0.25)
            q[0] = (R[k_index, j_index] - R[j_index, k_index]) / (4 * q_i)
            q[i_index + 1] = q_i
            q[j_index + 1] = (R[j_index, i_index] + R[i_index, j_index]) / (4 * q_i)
            q[k_index + 1] = (R[k_index, i_index] + R[i_index, k_index]) / (4 * q_i)

        return cls(q)

    @classmethod
    def from_basis(cls, t1, t2):
        """
        Create the object from two basis vectors t1, t2.
        t2 will be orthogonalized on t1, and t3 will be calculated with the
        cross product.
        """

        t1_normal = t1 / np.linalg.norm(t1)
        t2_ortho = t2 - t1_normal * np.dot(t1_normal, t2)
        t2_normal = t2_ortho / np.linalg.norm(t2_ortho)
        t3_normal = np.cross(t1_normal, t2_normal)

        R = np.transpose([t1_normal, t2_normal, t3_normal])
        return cls.from_rotation_matrix(R)

    @classmethod
    def from_rotation_vector(cls, rotation_vector):
        """Create the object from a rotation vector."""

        phi = np.linalg.norm(rotation_vector)
        if np.abs(phi) < mpy.eps_quaternion:
            return cls([0, 0, 0], 0)
        else:
            return cls(rotation_vector, phi)

    def check(self):
        """Perform all checks for the rotation."""
        self.check_uniqueness()
        self.check_quaternion_constraint()

    def check_uniqueness(self):
        """
        We always want q0 to be positive -> the range for the rotational angle
        is [-pi, pi].
        """

        if self.q[0] < 0:
            self.q = -self.q

    def check_quaternion_constraint(self):
        """We want to check that q.q = 1."""

        if np.abs(1 - np.linalg.norm(self.q)) > mpy.eps_quaternion:
            raise ValueError(
                f"The rotation object is corrupted. q.q does not equal 1! q={self.q}"
            )

    def get_rotation_matrix(self):
        """
        Return the rotation matrix for this rotation.
        (Krenk (3.50))
        """
        q_skew = skew_matrix(self.q[1:])
        R = (
            (self.q[0] ** 2 - np.dot(self.q[1:], self.q[1:])) * np.eye(3)
            + 2 * self.q[0] * q_skew
            + 2 * ([self.q[1:]] * np.transpose([self.q[1:]]))
        )

        return R

    def get_quaternion(self):
        """
        Return the quaternion for this rotation, as tuple.
        """

        return np.array(self.q)

    def get_rotation_vector(self):
        """
        Return the rotation vector for this object.
        """

        self.check()

        norm = np.linalg.norm(self.q[1:])
        phi = 2 * np.arctan2(norm, self.q[0])

        # Check if effective rotation angle is 0.
        if np.abs(np.sin(phi / 2)) < mpy.eps_quaternion:
            return np.zeros(3)
        else:
            return phi * self.q[1:] / norm

    def inv(self):
        """
        Return the inverse of this rotation.
        """

        tmp_quaternion = self.q.copy()
        tmp_quaternion[0] *= -1.0
        return Rotation(tmp_quaternion)

    def __mul__(self, other):
        """
        Add this rotation to another, or apply it on a vector.
        """

        # Check if the other object is also a rotation.
        if isinstance(other, Rotation):
            # Get quaternions of the two objects.
            p = self.q
            q = other.q
            # Add the rotations.
            added_rotation = np.zeros_like(self.q)
            added_rotation[0] = p[0] * q[0] - np.dot(p[1:], q[1:])
            added_rotation[1:] = p[0] * q[1:] + q[0] * p[1:] + np.cross(p[1:], q[1:])
            return Rotation(added_rotation)
        elif isinstance(other, (list, np.ndarray)) and len(other) == 3:
            # Apply rotation to vector.
            return np.dot(self.get_rotation_matrix(), np.array(other))
        raise NotImplementedError("Error, not implemented, does not make sense anyway!")

    def __eq__(self, other):
        """
        Check if the other rotation is equal to this one
        """

        if isinstance(other, Rotation):
            return bool(
                (np.linalg.norm(self.q - other.q) < mpy.eps_quaternion)
                or (np.linalg.norm(self.q + other.q) < mpy.eps_quaternion)
            )
        else:
            return object.__eq__(self, other)

    def get_dat(self):
        """
        Return a string with the triad components for the .dat line
        """

        rotation_vector = self.get_rotation_vector()

        # The zeros are added to avoid negative zeros in the input file.
        return " ".join(
            [
                (
                    mpy.dat_precision.format(component + 0)
                    if np.abs(component) >= mpy.eps_quaternion
                    else "0"
                )
                for component in rotation_vector
            ]
        )

    def copy(self):
        """Return a deep copy of this object."""
        return Rotation(self.q)

    def __str__(self):
        """
        String representation of object.
        """

        self.check()
        return f"Rotation:\n    q0: {self.q[0]}\n    q: {self.q[1:]}"


def get_relative_rotation(rotation1, rotation2):
    """Return the rotation from rotation1 to rotation2."""
    return rotation2 * rotation1.inv()


def add_rotations(rotation_21, rotation_10):
    """
    Multiply a rotation onto another.

    Args
    ----
    rotation_10: np.ndarray
        Array with the dimensions n x 4 or 4 x 1.
        The first rotation that is applied.
    rotation_21: np.ndarray
        Array with the dimensions n x 4 or 4 x 1.
        The second rotation that is applied.

    Return
    ----
    rot_new: np.ndarray
        Array with the dimensions n x 4.
        This array contains the new quaternions.
    """

    # Transpose the arrays, to work with the following code.
    if isinstance(rotation_10, Rotation):
        rot1 = rotation_10.get_quaternion().transpose()
    else:
        rot1 = np.transpose(rotation_10)
    if isinstance(rotation_21, Rotation):
        rot2 = rotation_21.get_quaternion().transpose()
    else:
        rot2 = np.transpose(rotation_21)

    if rot1.size > rot2.size:
        rotnew = np.zeros_like(rot1)
    else:
        rotnew = np.zeros_like(rot2)

    # Multiply the two rotations (code is taken from /utility/rotation.nb).
    rotnew[0] = (
        rot1[0] * rot2[0] - rot1[1] * rot2[1] - rot1[2] * rot2[2] - rot1[3] * rot2[3]
    )
    rotnew[1] = (
        rot1[1] * rot2[0] + rot1[0] * rot2[1] + rot1[3] * rot2[2] - rot1[2] * rot2[3]
    )
    rotnew[2] = (
        rot1[2] * rot2[0] - rot1[3] * rot2[1] + rot1[0] * rot2[2] + rot1[1] * rot2[3]
    )
    rotnew[3] = (
        rot1[3] * rot2[0] + rot1[2] * rot2[1] - rot1[1] * rot2[2] + rot1[0] * rot2[3]
    )

    return rotnew.transpose()


def rotate_coordinates(coordinates, rotation, *, origin=None):
    """
    Rotate all given coordinates

    Args
    ----
    coordinates: np.array
        Array of 3D coordinates to be rotated
    rotation: Rotation, list(quaternions) (nx4)
        The rotation that will be applied to the coordinates. Can also be an
        array with a quaternion for each coordinate.
    origin: 3D vector
        If this is given, the mesh is rotated about this point. Default is
        (0,0,0)
    """

    if isinstance(rotation, Rotation):
        rotation = rotation.get_quaternion().transpose()

    # Check if origin has to be added
    if origin is None:
        origin = [0.0, 0.0, 0.0]

    # New position array
    coordinates_new = np.zeros_like(coordinates)

    # Evaluate the new positions using the numpy data structure
    # (code is taken from /utility/rotation.nb)
    rotation = rotation.transpose()

    q0_q0 = np.square(rotation[0])
    q0_q1_2 = 2.0 * rotation[0] * rotation[1]
    q0_q2_2 = 2.0 * rotation[0] * rotation[2]
    q0_q3_2 = 2.0 * rotation[0] * rotation[3]

    q1_q1 = np.square(rotation[1])
    q1_q2_2 = 2.0 * rotation[1] * rotation[2]
    q1_q3_2 = 2.0 * rotation[1] * rotation[3]

    q2_q2 = np.square(rotation[2])
    q2_q3_2 = 2.0 * rotation[2] * rotation[3]

    q3_q3 = np.square(rotation[3])

    coordinates_new[:, 0] = (
        (q0_q0 + q1_q1 - q2_q2 - q3_q3) * (coordinates[:, 0] - origin[0])
        + (q1_q2_2 - q0_q3_2) * (coordinates[:, 1] - origin[1])
        + (q0_q2_2 + q1_q3_2) * (coordinates[:, 2] - origin[2])
    )
    coordinates_new[:, 1] = (
        (q1_q2_2 + q0_q3_2) * (coordinates[:, 0] - origin[0])
        + (q0_q0 - q1_q1 + q2_q2 - q3_q3) * (coordinates[:, 1] - origin[1])
        + (-q0_q1_2 + q2_q3_2) * (coordinates[:, 2] - origin[2])
    )
    coordinates_new[:, 2] = (
        (-q0_q2_2 + q1_q3_2) * (coordinates[:, 0] - origin[0])
        + (q0_q1_2 + q2_q3_2) * (coordinates[:, 1] - origin[1])
        + (q0_q0 - q1_q1 - q2_q2 + q3_q3) * (coordinates[:, 2] - origin[2])
    )

    if origin is not None:
        coordinates_new += origin

    return coordinates_new


def smallest_rotation(q: Rotation, t):
    """
    Get the triad that results from the smallest rotation (rotation without twist) from
    the triad q such that the rotated first basis vector aligns with t. For more details
    see Christoph Meier's dissertation chapter 2.1.2.

    Args
    ----
    q: Rotation
        Starting triad.
    t: Vector in R3
        Direction of the first basis of the rotated triad.
    Return
    ----
    q_sr: Rotation
        The triad that results from a smallest rotation.
    """

    R_old = q.get_rotation_matrix()
    g1_old = R_old[:, 0]
    g2_old = R_old[:, 1]
    g3_old = R_old[:, 2]

    g1 = np.array(t) / np.linalg.norm(t)
    g2 = g2_old - np.dot(g2_old, g1) / (1 + np.dot(g1_old, g1)) * (g1 + g1_old)
    g3 = g3_old - np.dot(g3_old, g1) / (1 + np.dot(g1_old, g1)) * (g1 + g1_old)

    return Rotation.from_rotation_matrix(np.transpose([g1, g2, g3]))


class SpecialEuclideanGroup:
    """This class represents an object of the special Euclidean group (SE3)

    For more details have a look at:
    Sonneville, V., Cardona, A., and Brüls, O., 2014, “Geometrically Exact Beam Finite Element
    Formulated on the Special Euclidean Group SE(3),” Computer Methods in Applied Mechanics and
    Engineering, 268, pp. 451–474.
    """

    def __init__(self, position, rotation: Rotation) -> None:
        """Set the position and the rotation for this SE3 element"""
        self.position = position
        self.rotation = rotation

    @classmethod
    def from_matrix_representation(cls, H):
        """Return the object from a matrix representation"""
        if not H.shape == (4, 4):
            raise ValueError(f"Matrix shape is wrong, got {H.shape}")
        if not np.allclose(H[3, :], [0, 0, 0, 1]):
            raise ValueError(f"Last row in matrix is wrong {H[3, :]}")

        rotation = Rotation.from_rotation_matrix(H[:3, :3])
        position = H[:3, 3]
        return cls(position, rotation)

    def inv(self):
        """Return the inverse of this SE3 object."""
        return SpecialEuclideanGroup.from_matrix_representation(
            self.get_inv_matrix_representation()
        )

    def __mul__(self, other):
        """Multiply this SE3 item with another SE3 item"""

        # Check if the other object is also a SE3 item
        if isinstance(other, SpecialEuclideanGroup):
            H_self = self.get_matrix_representation()
            H_other = other.get_matrix_representation()
            return SpecialEuclideanGroup.from_matrix_representation(
                np.dot(H_self, H_other)
            )
        raise NotImplementedError("Error, not implemented, does not make sense anyway!")

    def get_matrix_representation(self):
        """Return the SE3 matrix representation"""
        H = np.zeros((4, 4))
        H[:3, :3] = self.rotation.get_rotation_matrix()
        H[:3, 3] = self.position
        H[3, 3] = 1.0
        return H

    def get_inv_matrix_representation(self):
        """Return the SE3 matrix representation"""
        H = np.zeros((4, 4))
        rotation_matrix_inv = self.rotation.get_rotation_matrix().transpose()
        H[:3, :3] = rotation_matrix_inv
        H[:3, 3] = -np.dot(rotation_matrix_inv, self.position)
        H[3, 3] = 1.0
        return H

    def Log_SE3(self):
        """Return the Log_SE3 of this SE3 element.
        See Section A.2 of the reference mentioned above."""

        omega = self.rotation.get_rotation_vector()

        Log_SE3 = H = np.zeros(6)
        Log_SE3[:3] = np.dot(
            self.get_transformation_matrix_inv(omega).transpose(), self.position
        )
        Log_SE3[3:] = omega

        return Log_SE3

    @classmethod
    def from_Exp_SE3(cls, h):
        """Return the exponential map Exp_SE3 for a given vector h.
        See Section A.2 of the reference mentioned above."""
        h_u = h[:3]
        h_omega = h[3:]
        T = cls.get_transformation_matrix(h_omega)
        return cls(np.dot(T.transpose(), h_u), Rotation.from_rotation_vector(h_omega))

    @staticmethod
    def get_transformation_matrix(rotation_vector):
        """Return the transformation matrix for a given rotation"""

        omega = rotation_vector
        omega_norm = np.linalg.norm(omega)
        omega_skew = skew_matrix(omega)

        # Note: We have to take the square here, since there there is a division by the
        # square of the angle in the following formula
        if omega_norm**2 > mpy.eps_quaternion:
            alpha = np.sin(omega_norm) / omega_norm
            beta = 2.0 * (1.0 - np.cos(omega_norm)) / omega_norm**2
            T_SO3 = (
                np.identity(3)
                - 0.5 * beta * omega_skew
                + (1.0 - alpha) / omega_norm**2 * (np.dot(omega_skew, omega_skew))
            )
        else:
            # Note: this is the constant part of the Taylor series expansion. If this function is
            # derived with automatic differentiation, we have to also add higher order terms.
            T_SO3 = np.identity(3)
        return T_SO3

    @staticmethod
    def get_transformation_matrix_inv(rotation_vector):
        """Return the inverse of the transformation matrix for a given rotation"""

        omega = rotation_vector
        omega_norm = np.linalg.norm(omega)
        omega_skew = skew_matrix(omega)

        # Note: We have to take the square here, since there there is a division by the
        # square of the angle in the following formula
        if omega_norm**2 > mpy.eps_quaternion:
            alpha = np.sin(omega_norm) / omega_norm
            beta = 2.0 * (1.0 - np.cos(omega_norm)) / omega_norm**2
            T_inv_SO3 = (
                np.identity(3)
                + 0.5 * omega_skew
                + 1.0
                / omega_norm**2
                * (1 - alpha / beta)
                * (np.dot(omega_skew, omega_skew))
            )
        else:
            # Note: this is the constant part of the Taylor series expansion. If this function is
            # derived with automatic differentiation, we have to also add higher order terms.
            T_inv_SO3 = np.identity(3)
        return T_inv_SO3


def SE3_interpolation_between_nodes(positions, rotations, parameter_coordinates):
    """Perform a spatial interpolation between elements of SE3.

    This is a mix between the approach presented in Sonneville (2014) and Crisfield, Jelenic (1999).
    We generate an averaged element in SE3 (arithmetic average of position and normed average of
    the quaternion components) and interpolate the relative logarithmic maps in SE3 for each node.

    Args
    ----
    positions: list(3D array)
        List of the nodal positions
    rotations: list(Rotation)
        List of the nodal rotations
    parameter_coordinates: list(float)
        List of the parameter coordinates for each node
    """

    from scipy.interpolate import lagrange

    # Get averaged position and rotation
    r_average = np.zeros(3)
    q_average = np.zeros(4)
    for position, rotation in zip(positions, rotations):
        r_average += position
        q_average += rotation.q
    r_average = r_average / len(position)
    q_average = q_average / np.linalg.norm(q_average)
    q_average = Rotation(q_average)
    H_average = SpecialEuclideanGroup(r_average, q_average)
    H_average_inv = H_average.inv()

    # Get the relative SE3 distance for each node
    d_nodes = np.zeros((len(positions), 6))
    for i_node, (position, rotation) in enumerate(zip(positions, rotations)):
        H_node = SpecialEuclideanGroup(position, rotation)
        H_rel = H_average_inv * H_node
        d_nodes[i_node] = H_rel.Log_SE3()

    # Get the interpolation function
    interpolation_functions = [
        lagrange(parameter_coordinates, d_nodes[:, i]) for i in range(6)
    ]

    def interpolate(xi):
        H = np.array([interpolation_functions[i](xi) for i in range(6)])
        H_interpolated = H_average * SpecialEuclideanGroup.from_Exp_SE3(H)
        return (H_interpolated.position, H_interpolated.rotation)

    return interpolate
