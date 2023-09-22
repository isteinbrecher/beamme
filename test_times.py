import numpy as np
from meshpy import mpy
from meshpy.geometric_search import find_close_points, FindClosePointAlgorithm
import time


# Import matplotlib in headless mode
# fmt: off
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# fmt: on


def get_times(n, algorithms):
    points = np.random.rand(n, 3)

    times = {}

    for algorithm in algorithms:
        if n < 500:
            n_iter = 500
        elif n < 1000:
            n_iter = 100
        elif n < 10000:
            n_iter = 10
        else:
            n_iter = 1

        start = time.time()
        for i in range(n_iter):
            find_close_points(points, algorithm=algorithm)
        search_time = (time.time() - start) / n_iter
        times[algorithm] = search_time

    return times


get_times(
    10,
    [
        FindClosePointAlgorithm.binning_cython,
        FindClosePointAlgorithm.boundary_volume_hierarchy_arborx,
        FindClosePointAlgorithm.brute_force_cython,
        FindClosePointAlgorithm.kd_tree_scipy,
    ],
)
mpy.geometric_search_max_nodes_brute_force = 1e10

times = {
    FindClosePointAlgorithm.binning_cython: [],
    FindClosePointAlgorithm.boundary_volume_hierarchy_arborx: [],
    FindClosePointAlgorithm.brute_force_cython: [],
    FindClosePointAlgorithm.kd_tree_scipy: [],
}

n_points = [10 * 2**i for i in range(14)] + [2**16, 2**18]
print(n_points)
for n in n_points:

    if n < 20000:
        algorithms = [
            FindClosePointAlgorithm.binning_cython,
            FindClosePointAlgorithm.boundary_volume_hierarchy_arborx,
            FindClosePointAlgorithm.brute_force_cython,
            FindClosePointAlgorithm.kd_tree_scipy,
        ]
    else:
        algorithms = [
            FindClosePointAlgorithm.binning_cython,
            FindClosePointAlgorithm.boundary_volume_hierarchy_arborx,
            FindClosePointAlgorithm.kd_tree_scipy,
        ]

    my_time = get_times(n, algorithms)
    for key in my_time.keys():
        times[key].append(my_time[key])


for key in times.keys():
    plt.loglog(n_points[: len(times[key])], times[key])
plt.savefig("times.png")
