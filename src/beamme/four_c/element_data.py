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
"""Define a data class for 4C element data."""

from dataclasses import dataclass as _dataclass
from dataclasses import field as _field
from typing import Any as _Any

import numpy as _np
from numpy.typing import NDArray as _NDArray

from beamme.utils.data_structures import (
    compare_nested_dicts_or_lists as _compare_nested_dicts_or_lists,
)


@_dataclass(eq=False)
class FourCElementData:
    """Class that contains the data for a 4C element block."""

    four_c_type: str
    four_c_cell: str
    element_technology: dict[str, _Any] = _field(default_factory=dict)

    def get_legacy_dict(
        self, element_id, connectivity, element_material_id, additional_element_data
    ) -> dict:
        """Return the dictionary to write this element data to a legacy element
        definition in the input file."""

        return {
            "id": element_id + 1,
            "cell": {
                "type": self.four_c_cell,
                "connectivity": connectivity + 1,
            },
            "data": {
                "type": self.four_c_type,
                **(
                    {"MAT": element_material_id + 1}
                    if element_material_id != -1
                    else {}
                ),
                **self.element_technology,
                **(additional_element_data or {}),
            },
        }

    def __eq__(self, other) -> bool:
        """Check if two 4C element data objects are equal."""
        if not isinstance(other, FourCElementData):
            return False
        if self.four_c_cell != other.four_c_cell:
            return False
        if self.four_c_type != other.four_c_type:
            return False
        if not _compare_nested_dicts_or_lists(
            self.element_technology, other.element_technology
        ):
            return False
        return True


def four_c_element_data_from_legacy_dict(
    legacy_dict: dict,
) -> tuple[FourCElementData, int, _NDArray, int]:
    """Extract the 4C element data from a legacy element definition in the
    input file.

    Args:
        legacy_dict: The legacy element definition in the input file, will be modified in place.

    Returns:
        A tuple containing the 4C element data, the element ID, the connectivity, and the material ID.
    """

    element_id = legacy_dict["id"]
    connectivity = _np.array(legacy_dict["cell"]["connectivity"], dtype=int) - 1
    four_c_cell = legacy_dict["cell"]["type"]

    # Since we directly store `legacy_dict["data"]` as the element technology data,
    # the material ID and four_c_type have to be removed from the `legacy_dict["data"]`
    # dictionary, thus the `pop`.
    material_id = legacy_dict["data"].pop("MAT", -1)
    four_c_type = legacy_dict["data"].pop("type")

    data = FourCElementData(
        four_c_type=four_c_type,
        four_c_cell=four_c_cell,
        element_technology=legacy_dict["data"],
    )

    return data, element_id, connectivity, material_id
