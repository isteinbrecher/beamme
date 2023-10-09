# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import autograd.numpy as npAD
from meshpy import Rotation, InputFile, Beam3rHerm2Line3, MaterialReissner
from meshpy.mesh_creation_functions import create_beam_mesh_curve


from cProfile import Profile
from pstats import SortKey, Stats


# Create input file.
input_file = InputFile(maintainer="Ivo Steinbrecher")

# Add material and functions.
mat = MaterialReissner()

# Set parameters for the helix.
R = 2.0
tz = 4.0  # incline
n = 1  # number of turns
n_el = 5

# Create a helix with a parametric curve.
def helix(t):
    factor = 2
    t_trans = npAD.exp(factor * t / (2.0 * np.pi * n)) * t / npAD.exp(factor)
    return npAD.array(
        [
            R * npAD.cos(t_trans),
            R * npAD.sin(t_trans),
            t_trans * tz / (2 * np.pi),
        ]
    )


def meins():
    helix_set = create_beam_mesh_curve(
        input_file, Beam3rHerm2Line3, mat, helix, [0.0, 2.0 * np.pi * n], n_el=n_el
    )


def fib(n):
    return n if n < 2 else fib(n - 2) + fib(n - 1)


with Profile() as profile:
    print(f"{meins() = }")
    (Stats(profile).strip_dirs().sort_stats(SortKey.CALLS).print_stats())
