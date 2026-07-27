# %% [markdown]
# # A moving VTK object
#
# This example creates a short sequence of VTK scenes in which a cube moves
# along the $x$-axis. The visualization can be stepped forwards and backwards.
#
# In Jupyter, both the static frames and the interactive 3D scenes are embedded
# in the notebook output. The same two views are exported when this notebook is
# included in the Sphinx documentation.

# %%
import pyvista as pv
import vtk

from beamme.utils.environment import is_testing
from beamme.utils.visualization import show_plotter_sequence

# %%
# A separate VTK data set and PyVista plotter are created for every time step.
# Keeping the camera and scene bounds fixed makes the cube's motion easy to
# compare between frames.
cube_positions = [-2.0, -1.0, 0.0, 1.0, 2.0]
plotters = []

for x_position in cube_positions:
    cube_source = vtk.vtkCubeSource()
    cube_source.SetCenter(x_position, 0.0, 0.0)
    cube_source.SetXLength(0.8)
    cube_source.SetYLength(0.8)
    cube_source.SetZLength(0.8)
    cube_source.Update()

    plotter = pv.Plotter(window_size=(900, 500))
    plotter.add_mesh(
        pv.wrap(cube_source.GetOutput()),
        color="royalblue",
        show_edges=True,
    )
    plotter.show_bounds(
        bounds=(-2.8, 2.8, -1.2, 1.2, -1.2, 1.2),
        xtitle="x",
        ytitle="y",
        ztitle="z",
    )
    plotter.camera_position = [
        (6.5, -8.0, 5.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
    ]
    plotters.append(plotter)

# %% [markdown]
# Use the buttons below the scene to move through the time steps. Switching to
# **Interactive 3D** allows each frame to be rotated, panned, and zoomed.

# %%
if not is_testing():
    show_plotter_sequence(
        plotters,
        labels=[f"x = {x_position:.1f}" for x_position in cube_positions],
    )
