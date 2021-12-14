# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# MeshPy: A beam finite element input generator
#
# MIT License
#
# Copyright (c) 2021 Ivo Steinbrecher
#                    Institute for Mathematics and Computer-Based Simulation
#                    Universitaet der Bundeswehr Muenchen
#                    https://www.unibw.de/imcs-en
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
This module implements the class that represents one node in the Mesh.
"""

# Python modules.
import numpy as np

# Meshpy modules.
from .conf import mpy
from .base_mesh_item import BaseMeshItem


class Node(BaseMeshItem):
    """
    This object represents one node in the mesh. The node can have a rotation
    and can be rotated moved and so on.
    """

    def __init__(self, coordinates, rotation=None, is_middle_node=False,
                is_dat=False, **kwargs):
        super().__init__(data=None, is_dat=is_dat, **kwargs)

        # Coordinates and rotation of this node.
        self.coordinates = np.array(coordinates)
        if rotation is None:
            self.rotation = None
        else:
            self.rotation = rotation.copy()

        # If this node is at the end of a line or curve (by default only those
        # nodes are checked for overlapping nodes).
        self.is_end_node = False

        # If the node is in the middle of a beam element.
        self.is_middle_node = is_middle_node

        # Lists with the objects that this node is linked to.
        self.element_link = []
        self.node_sets_link = []
        self.element_partner_index = None
        self.coupling_link = None
        self.mesh = None

        # If this node is replaced, store a link to the remaining node.
        self.master_node = None

    @classmethod
    def from_dat(cls, input_line):
        """Create the Node object from a line in the input file."""

        # Split up the input line.
        line_split = input_line[0].split()

        # Convert the node coordinates into a Node object.
        return cls([float(line_split[i]) for i in range(3, 6)], is_dat=True,
            comments=input_line[1])

    def get_master_node(self):
        """
        Return the master node of this node. If the node has not been replaced,
        this object is returned.
        """

        if self.master_node is None:
            return self
        else:
            return self.master_node.get_master_node()

    def replace_with(self, master_node):
        """Replace this node with another node object."""

        # Replace the links to this node in the referenced objects.
        self.mesh.replace_node(self, master_node)
        for element in self.element_link:
            element.replace_node(self, master_node)
        for node_set in self.node_sets_link:
            node_set.replace_node(self, master_node)

        # Set link to master node.
        self.master_node = master_node.get_master_node()

    def unlink(self):
        """Reset the links to elements, node sets and global indices."""
        self.element_link = []
        self.node_sets_link = []
        self.coupling_link = None
        self.mesh = None
        self.n_global = None

    def rotate(self, rotation, origin=None, only_rotate_triads=False):
        """
        Rotate this node. By default the node is rotated around the origin
        (0,0,0), if the keyword argument origin is given, it is rotated around
        that point.
        If only_rotate_triads is True, then only the rotation is affected,
        the position of the node stays the same.
        """

        # If the node has a rotation, rotate it.
        if self.rotation is not None:
            self.rotation = rotation * self.rotation

            # Rotate the positions (around origin).
            if not only_rotate_triads:
                if origin is not None:
                    self.coordinates = self.coordinates - origin
                self.coordinates = rotation * self.coordinates
                if origin is not None:
                    self.coordinates = self.coordinates + origin

    def _get_dat(self):
        """
        Return the line that represents this node in the input file.
        """

        coordinate_string = ' '.join([
            mpy.dat_precision.format(component + 0)
            if np.abs(component) >= mpy.eps_pos
            else '0'
            for component in self.coordinates
            ])
        return 'NODE {} COORD {}'.format(self.n_global, coordinate_string)
