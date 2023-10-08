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
This module defines classes and functions to create and edit a Baci input file.
"""

# Global configuration object.
from .conf import mpy

# Utility functions.
from .utility import clean_simulation_directory, flatten

# 3D rotations for nodes.
from .rotation import Rotation

# Mesh items.
from .base_mesh_item import BaseMeshItem
from .function import Function
from .material import (
    MaterialBeam,
    MaterialEulerBernoulli,
    MaterialKirchhoff,
    MaterialReissner,
    MaterialReissnerElastoplastic,
    MaterialStVenantKirchhoff,
)
from .element_beam import Beam3eb, Beam3k, Beam3rHerm2Line3, Beam3rLine2Line2
from .geometry_set import GeometrySet, GeometrySetNodes

# Boundary conditions and couplings for geometry in the mesh.
from .boundary_condition import BoundaryCondition
from .coupling import Coupling

# The mesh class itself and the input file classes.
from .mesh import Mesh
from .inputfile import InputFile, InputSection, InputSectionMultiKey

# Functions to set default header options.
from .header_functions import (
    get_comment,
    get_yes_no,
    set_beam_to_solid_meshtying,
    set_header_static,
    set_runtime_output,
)

# Sub modules to be included in the main module module
from . import abaqus
from . import mesh_creation_functions

# Define the items that will be exported by default.
__all__ = [
    # Option object.
    "mpy",
    # Basic stuff.
    "Rotation",
    "BaseMeshItem",
    "Function",
    "MaterialReissner",
    "MaterialReissnerElastoplastic",
    "MaterialKirchhoff",
    "MaterialBeam",
    "GeometrySet",
    "GeometrySetNodes",
    "BoundaryCondition",
    "Coupling",
    "MaterialEulerBernoulli",
    "MaterialStVenantKirchhoff",
    # Mesh items.
    "Beam3rHerm2Line3",
    "Beam3rLine2Line2",
    "Beam3k",
    "Mesh",
    "InputFile",
    "InputSection",
    "Beam3eb",
    # Header functions.
    "set_header_static",
    "set_runtime_output",
    "set_beam_to_solid_meshtying",
    # Sub modules.
    "abaqus",
    "mesh_creation_functions",
]
