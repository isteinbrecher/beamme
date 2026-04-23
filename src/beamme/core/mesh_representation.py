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
"""This module defines the mesh representation data structure."""

from dataclasses import dataclass as _dataclass
from typing import Any as _Any

import numpy as _np
import pyvista as _pv
from numpy.typing import NDArray as _NDArray

from beamme.core.conf import Geometry as _Geometry
from beamme.core.conf import bme as _bme

MESH_REPRESENTATION_MAPPINGS: dict[str, _Any] = {}
# fmt: off
MESH_REPRESENTATION_MAPPINGS[
    "element_type_and_n_nodes_to_connectivity_mapping_beamme_to_vtk"
] = {

    # Only list the non-standard mappings
    (_bme.element_type.solid, 20):
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15],
    (_bme.element_type.solid, 27):
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15,
         24, 22, 21, 23, 20, 25, 26],
}
# fmt: on
MESH_REPRESENTATION_MAPPINGS[
    "element_type_and_n_nodes_to_connectivity_mapping_vtk_to_beamme"
] = {
    # Only list the non-standard mappings
    (_bme.element_type.solid, 20): _np.argsort(
        MESH_REPRESENTATION_MAPPINGS[
            "element_type_and_n_nodes_to_connectivity_mapping_beamme_to_vtk"
        ][(_bme.element_type.solid, 20)]
    ),
    (_bme.element_type.solid, 27): _np.argsort(
        MESH_REPRESENTATION_MAPPINGS[
            "element_type_and_n_nodes_to_connectivity_mapping_beamme_to_vtk"
        ][(_bme.element_type.solid, 27)]
    ),
}

GEOMETRY_SET_INFO_PREFIX = "geometry_set"


@_dataclass
class GeometrySetInfo:
    """Data class that contains information of a geometry set."""

    geometry_type: _Geometry
    i_global: int
    name: str | None = None
    point_flag_vector: _NDArray | None = None
    cell_flag_vector: _NDArray | None = None

    def __str__(self):
        """Return a string representation of this geometry set info - the node and cell flags are not included."""
        name = self.name
        return (
            f"{GEOMETRY_SET_INFO_PREFIX}_{self.i_global}_{self.geometry_type.name}"
            + (f"_{name}" if name is not None else "")
        )


def string_to_geometry_set_info(name: str) -> GeometrySetInfo | None:
    """Extract the geometry set information from a given string."""
    if not name.startswith(GEOMETRY_SET_INFO_PREFIX):
        return None

    name_without_prefix = name[len(GEOMETRY_SET_INFO_PREFIX) + 1 :]
    split = name_without_prefix.split("_", 2)
    i_global = int(split[0])
    geometry_type = _Geometry[split[1]]
    return GeometrySetInfo(
        i_global=i_global,
        geometry_type=geometry_type,
        name=split[2] if len(split) == 3 else None,
    )


class MeshRepresentation:
    """Class representing a generic mesh."""

    def __init__(
        self,
        cell_connectivity: _NDArray[_np.integer] | None = None,
        cell_types: _NDArray[_np.integer] | None = None,
        points: _NDArray[_np.floating] | None = None,
        geometry_sets: list[GeometrySetInfo] | None = None,
        cell_data: dict[str, _NDArray | None] | None = None,
        point_data: dict[str, _NDArray | None] | None = None,
    ):
        if cell_connectivity is None or len(cell_connectivity) == 0:
            # In this case, we have an empty mesh representation.
            # Safety check that all provided data is None or empty
            if (
                (cell_types is not None and len(cell_types) > 0)
                or (points is not None and len(points) > 0)
                or (geometry_sets is not None and len(geometry_sets) > 0)
                or (cell_data is not None and len(cell_data) > 0)
                or (point_data is not None and len(point_data) > 0)
            ):
                raise ValueError(
                    "If cell_connectivity is None or empty, all other data must be None or empty as well."
                )
            self._grid = None
        else:
            # We use a unstructured vtk grid as main mesh data structure
            self._grid = _pv.UnstructuredGrid(cell_connectivity, cell_types, points)

            # Store node set information
            if geometry_sets is not None:
                for node_set in geometry_sets:
                    self._grid.point_data[str(node_set)] = node_set.point_flag_vector

            # Store the required data in the grid.
            for data, attribute_name in (
                (cell_data, "cell_data"),
                (point_data, "point_data"),
            ):
                if data is not None:
                    data_attribute = getattr(self._grid, attribute_name)
                    for key, value in data.items():
                        if value is not None:
                            data_attribute[key] = value

    def is_empty(self) -> bool:
        """Return True if this mesh representation is empty, i.e., it does not
        contain any grid."""
        return self._grid is None

    @property
    def grid(self) -> _pv.UnstructuredGrid:
        """Return the internal grid of this mesh representation."""
        if self.is_empty():
            raise ValueError("This mesh representation does not contain a grid.")
        return self._grid
