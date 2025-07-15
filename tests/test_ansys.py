# The MIT License (MIT)
#
# Copyright (c) 2018-2025 BeamMe Authors
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
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""This script is used to test the creation of Ansys input files."""

import numpy as np
import pytest

from beamme.abaqus.beam import generate_abaqus_beam
from beamme.abaqus.input_file import AbaqusBeamNormalDefinition, AbaqusInputFile
from beamme.abaqus.material import AbaqusBeamMaterial
from beamme.core.geometry_set import GeometrySet
from beamme.core.mesh import Mesh
from beamme.core.rotation import Rotation
from beamme.mesh_creation_functions.beam_arc import create_beam_mesh_arc_segment_2d
from beamme.mesh_creation_functions.beam_line import create_beam_mesh_line

PYTEST_ABAQUS_NORMAL_DEFINITION_PARAMETRIZE = [
    ("normal_definition", "additional_identifier"),
    [
        (AbaqusBeamNormalDefinition.normal, "normal"),
        (AbaqusBeamNormalDefinition.normal_and_extra_node, "normal_and_extra_node"),
    ],
]
ELEMENT = ["element_type", ["B31H", "B32H"]]


@pytest.mark.parametrize(*PYTEST_ABAQUS_NORMAL_DEFINITION_PARAMETRIZE)
@pytest.mark.parametrize(*ELEMENT)
def test_ansys_frame(
    normal_definition,
    additional_identifier,
    element_type,
    assert_results_close,
    get_corresponding_reference_file_path,
):
    """Create a frame out of connected beams with different materials."""

    mesh = Mesh()
    mat_1 = AbaqusBeamMaterial("beam_material")
    mat_2 = AbaqusBeamMaterial("beam_material_2")
    beam_type_b23 = generate_abaqus_beam(element_type)
    beam_type_b33 = generate_abaqus_beam(element_type)

    create_beam_mesh_line(mesh, beam_type_b23, mat_1, [0, 0, 0], [0, 2, 0], n_el=2)
    set_2 = create_beam_mesh_line(
        mesh, beam_type_b23, mat_1, [0, 2, 0], [1, 2, 0], n_el=2
    )
    # mesh.rotate(Rotation([1, 0, 0], np.pi * 0.5))
    # create_beam_mesh_line(mesh, beam_type_b23, mat_2, [1, 0, 0], [1, 1, 0], n_el=2)
    # create_beam_mesh_line(mesh, beam_type_b33, mat_1, [1, 1, 0], [1, 1, 1], n_el=2)
    # create_beam_mesh_line(mesh, beam_type_b33, mat_2, [1, 1, 1], [0, 1, 1], n_el=2)
    mesh.couple_nodes()

    fix_set = GeometrySet(mesh.nodes[0], name="fix_node")
    load_set = set_2["line"]
    load_set.name = "load_node"
    mesh.add(fix_set, load_set)

    mesh.display_pyvista(beam_radius_for_display=0.075)

    input_file = AbaqusInputFile(mesh)
    input_file.write_input_file(
        f"elbow_{element_type}_{additional_identifier}.inp",
        normal_definition=normal_definition,
    )
    # assert_results_close(
    #     get_corresponding_reference_file_path(
    #         additional_identifier=additional_identifier, extension="inp"
    #     ),
    #     input_file.get_input_file_string(normal_definition),
    #     atol=1e-15,
    # )


@pytest.mark.parametrize(*PYTEST_ABAQUS_NORMAL_DEFINITION_PARAMETRIZE)
@pytest.mark.parametrize(*ELEMENT)
def test_ansys_curved_beam(
    normal_definition,
    additional_identifier,
    element_type,
    assert_results_close,
    get_corresponding_reference_file_path,
):
    """Create a frame out of connected beams with different materials."""

    mesh = Mesh()
    mat_1 = AbaqusBeamMaterial("beam_material")
    mat_2 = AbaqusBeamMaterial("beam_material_2")
    beam_type_b23 = generate_abaqus_beam(element_type)
    beam_type_b33 = generate_abaqus_beam(element_type)

    create_beam_mesh_arc_segment_2d(
        mesh, beam_type_b23, mat_1, [0, 0, 0], 1, 0, np.pi, n_el=6
    )

    for i_node, node in enumerate(mesh.nodes):
        node.rotation = node.rotation * Rotation(
            [1, 0, 0], 2.0 * np.pi * i_node / (len(mesh.nodes) - 1)
        )

    # set_2 = create_beam_mesh_line(
    #     mesh, beam_type_b23, mat_1, [0, 2, 0], [1, 2, 0], n_el=2
    # )
    # mesh.rotate(Rotation([1, 0, 0], np.pi * 0.5))
    # create_beam_mesh_line(mesh, beam_type_b23, mat_2, [1, 0, 0], [1, 1, 0], n_el=2)
    # create_beam_mesh_line(mesh, beam_type_b33, mat_1, [1, 1, 0], [1, 1, 1], n_el=2)
    # create_beam_mesh_line(mesh, beam_type_b33, mat_2, [1, 1, 1], [0, 1, 1], n_el=2)
    # mesh.couple_nodes()

    fix_set = GeometrySet(mesh.nodes[0], name="fix_node")
    # load_set = set_2["line"]
    # load_set.name = "load_node"
    # mesh.add(fix_set, load_set)

    mesh.display_pyvista(beam_radius_for_display=0.075)

    input_file = AbaqusInputFile(mesh)
    input_file.write_input_file(
        f"arc_{element_type}_{additional_identifier}.inp",
        normal_definition=normal_definition,
    )
    # assert_results_close(
    #     get_corresponding_reference_file_path(
    #         additional_identifier=additional_identifier, extension="inp"
    #     ),
    #     input_file.get_input_file_string(normal_definition),
    #     atol=1e-15,
    # )
