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

from typing import Any as _Any


class FourCElementData:
    """Class that contains the data for a 4C element block."""

    def __init__(
        self,
        four_c_type: str,
        four_c_cell: str,
        element_technology: dict[str, _Any] | None = None,
    ):
        """Initialize the 4C element data."""
        self.four_c_type = four_c_type
        self.four_c_cell = four_c_cell
        if element_technology is None:
            element_technology = {}
        self.element_technology = element_technology

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
