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
This module implements the class that represents one element in the Mesh.
"""

# Meshpy modules.
from .base_mesh_item import BaseMeshItem


class Element(BaseMeshItem):
    """A base class for an FEM element in the mesh."""

    def __init__(self, nodes=None, material=None, is_dat=False, **kwargs):
        super().__init__(data=None, is_dat=is_dat, **kwargs)

        # List of nodes that are connected to the element.
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

        # Material of this element.
        self.material = material

        # VTK cell data for this element.
        self.vtk_cell_data = {}

    @classmethod
    def from_dat(cls, input_line):
        """
        Check the string to decide which element to use. Nodes are linked with
        the 0 based index and will be connected to the Node objects after the
        whole input file is parsed.
        """

        # Import solid element classes for creation of the element.
        from .element_volume import (VolumeHEX8, SolidRigidSphere, VolumeHEX27,
            VolumeHEX20, VolumeTET10, VolumeTET4)

        # Split up input line and get pre node string.
        line_split = input_line[0].split()
        dat_pre_nodes = ' '.join(line_split[1:3])

        # Get a list of the element nodes.
        element_nodes = []
        for i, item in enumerate(line_split[3:]):
            if item.isdigit():
                element_nodes.append(int(item) - 1)
            else:
                break
        else:
            raise ValueError(('The input line:\n"{}"\ncould not be converted '
                + 'to a solid element!').format(input_line))

        # Get the post node string
        dat_post_nodes = ' '.join(line_split[3 + i:])

        # Depending on the number of nodes chose which solid element to return.
        n_nodes = len(element_nodes)
        if n_nodes == 8:
            return VolumeHEX8(nodes=element_nodes, dat_pre_nodes=dat_pre_nodes,
                dat_post_nodes=dat_post_nodes, comments=input_line[1])
        elif len(element_nodes) == 4:
            return VolumeTET4(nodes=element_nodes, dat_pre_nodes=dat_pre_nodes,
                dat_post_nodes=dat_post_nodes, comments=input_line[1])
        elif len(element_nodes) == 10:
            return VolumeTET10(nodes=element_nodes, dat_pre_nodes=dat_pre_nodes,
                dat_post_nodes=dat_post_nodes, comments=input_line[1])
        elif len(element_nodes) == 20:
            return VolumeHEX20(nodes=element_nodes, dat_pre_nodes=dat_pre_nodes,
                dat_post_nodes=dat_post_nodes, comments=input_line[1])
        elif len(element_nodes) == 27:
            return VolumeHEX27(nodes=element_nodes, dat_pre_nodes=dat_pre_nodes,
                dat_post_nodes=dat_post_nodes, comments=input_line[1])
        elif len(element_nodes) == 1:
            return SolidRigidSphere(nodes=element_nodes,
                dat_pre_nodes=dat_pre_nodes, dat_post_nodes=dat_post_nodes,
                comments=input_line[1])
        else:
            raise TypeError('Could not find a element type for '
                + '{}, with {} nodes'.format(dat_pre_nodes, n_nodes))

    def replace_node(self, old_node, new_node):
        """Replace old_node with new_node."""

        # Look for old_node and replace it. If it is not found, throw error.
        for i, node in enumerate(self.nodes):
            if node == old_node:
                self.nodes[i] = new_node
                break
        else:
            raise ValueError('The node that should be replaced is not in the '
                + 'current element')

    def get_vtk(self, vtk_writer_beam, vtk_writer_solid):
        """
        Add representation of this element to the vtk_writers for solid and
        beam.
        """
        raise NotImplementedError(
            'VTK output has to be implemented in the class!')
