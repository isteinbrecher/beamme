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
from typing import Iterable as _Iterable

import numpy as _np
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
        def _convert_argument_numpy(
            argument: _NDArray | None,
            default_shape: tuple[int, ...],
            dtype: type,
        ) -> _NDArray:
            """Convert given array arguments so we can store them in this
            object."""
            if argument is None:
                return _np.empty(default_shape, dtype=dtype)
            else:
                return _np.asarray(argument, dtype=dtype)

        def _filter_none_entries(
            argument: dict[str, _NDArray | None] | None, expected_size: int
        ) -> dict[str, _NDArray]:
            """Check if a dictionary is given, and if so, filter None entries
            from it."""
            if argument is None:
                return {}
            else:
                data = {
                    key: _np.asarray(value)
                    for key, value in argument.items()
                    if value is not None
                }
                for key, value in data.items():
                    if len(value) != expected_size:
                        raise ValueError(
                            f"Data field {key} has size {len(value)}, but expected size is {expected_size}."
                        )
                return data

        self.cell_connectivity = _convert_argument_numpy(
            cell_connectivity, default_shape=(0,), dtype=int
        )
        self.cell_types = _convert_argument_numpy(
            cell_types, default_shape=(0,), dtype=int
        )
        self.points = _convert_argument_numpy(points, default_shape=(0, 3), dtype=float)

        self.cell_data = _filter_none_entries(cell_data, len(self.cell_types))
        self.point_data = _filter_none_entries(point_data, len(self.points))

        def _store_geometry_set_data(
            data_field: str,
            name: str,
            field_vector: _NDArray | None,
            expected_size: int,
        ) -> None:
            """Store point or cell data defining a geometry set."""
            if field_vector is not None:
                if len(field_vector) != expected_size:
                    raise ValueError(
                        f"Geometry set {name} has size {len(field_vector)},"
                        f" but expected size is {expected_size}."
                    )
                getattr(self, data_field)[str(name)] = field_vector

        if geometry_sets is not None:
            for geometry_set in geometry_sets:
                geometry_set_name = str(geometry_set)
                _store_geometry_set_data(
                    "point_data",
                    geometry_set_name,
                    geometry_set.point_flag_vector,
                    len(self.points),
                )
                _store_geometry_set_data(
                    "cell_data",
                    geometry_set_name,
                    geometry_set.cell_flag_vector,
                    len(self.cell_types),
                )

        # Get the offset array so we can iterate over the connectivity in a
        # performant manner.
        self.cell_connectivity_offsets = _np.zeros(len(self.cell_types) + 1, dtype=int)
        for i_cell in range(len(self.cell_types)):
            last_offset = self.cell_connectivity_offsets[i_cell]
            if last_offset >= len(self.cell_connectivity):
                raise ValueError(
                    "Invalid offset, is larger than the size of the connectivity."
                )
            cell_size = self.cell_connectivity[last_offset]
            if cell_size < 1:
                raise ValueError("Invalid offset, has to be at least 1.")
            self.cell_connectivity_offsets[i_cell + 1] = (
                cell_size + self.cell_connectivity_offsets[i_cell] + 1
            )
        if self.cell_connectivity_offsets[-1] != len(self.cell_connectivity):
            raise ValueError("Invalid cell connectivity offsets.")

    def connectivity_iterator(self) -> _Iterable:
        """Return an iterator over the cell connectivity.

        Returns:
            An iterator that returns the connectivity array for each cell.
        """
        for i in range(len(self.cell_types)):
            start = self.cell_connectivity_offsets[i] + 1
            end = self.cell_connectivity_offsets[i + 1]
            yield self.cell_connectivity[start:end]
