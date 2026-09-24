# The MIT License (MIT)
#
# Copyright (c) 2018-2026 BeamMe Authors
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
"""This script is used to test general functionality of the core mesh representation
class."""

from beamme.core.conf import bme
from beamme.core.element_beam import Beam3
from beamme.core.material import MaterialBeamBase
from beamme.core.mesh import Mesh
from beamme.mesh_creation_functions.beam_line import create_beam_mesh_line


def test_integration_core_mesh_representation_get_geometry_set_infos(
    assert_results_close,
):
    """Test the get_geometry_set_infos method."""
    mesh = Mesh()
    beam_set = create_beam_mesh_line(
        mesh, Beam3, MaterialBeamBase, [0, 0, 0], [1, 2, 3], n_el=3
    )
    beam_set["start"].name = "start"
    beam_set["line"].name = "line"
    mesh.add(beam_set["start"])
    mesh.add(beam_set["end"])
    mesh.add(beam_set["line"])

    mesh_representation = mesh.get_mesh_representation()[0]
    geometry_set_infos = mesh_representation.get_geometry_set_infos()

    assert len(geometry_set_infos) == 3

    for i_info, (
        name,
        geometry_type,
        i_global,
        point_flag_vector,
        cell_flag_vector,
    ) in enumerate(
        [
            ("start", bme.geo.point, 0, [1, 0, 0, 0, 0, 0, 0], None),
            (None, bme.geo.point, 1, [0, 0, 0, 0, 0, 0, 1], None),
            ("line", bme.geo.line, 2, None, [1, 1, 1]),
        ]
    ):
        geometry_set_info = geometry_set_infos[i_info]
        assert geometry_set_info.name == name
        assert geometry_set_info.i_global == i_global
        assert geometry_set_info.geometry_type == geometry_type
        assert_results_close(geometry_set_info.point_flag_vector, point_flag_vector)
        assert_results_close(geometry_set_info.cell_flag_vector, cell_flag_vector)
