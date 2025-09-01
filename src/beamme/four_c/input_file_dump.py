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
"""This file defines functions to dump mesh items for 4C."""

from beamme.four_c.four_c_types import BeamType as _BeamType


def dump_coupling(coupling):
    """Return the input file representation of the coupling condition."""

    if isinstance(coupling.data, dict):
        data = coupling.data
    else:
        # In this case we have to check which beams are connected to the node.
        # TODO: Coupling also makes sense for different beam types, this can
        # be implemented at some point.
        nodes = coupling.geometry_set.get_points()
        connected_elements = [
            element for node in nodes for element in node.element_link
        ]
        element_types = {type(element) for element in connected_elements}
        if len(element_types) > 1:
            raise TypeError(
                f"Expected a single connected type of beam elements, got {element_types}"
            )
        element_type = element_types.pop()
        if element_type.four_c_beam_type is _BeamType.kirchhoff:
            rotvec = {
                type(element).four_c_element_data["ROTVEC"]
                for element in connected_elements
            }
            if len(rotvec) > 1 or not rotvec.pop():
                raise TypeError(
                    "Couplings for Kirchhoff beams and rotvec==False not yet implemented."
                )

        data = element_type.get_coupling_dict(coupling.data)

    return {"E": coupling.geometry_set.i_global + 1, **data}
