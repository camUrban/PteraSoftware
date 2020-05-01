"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the solution to the aerodynamics problem.
"""

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt


def draw(aerodynamics_problem, points_type=None):
    """Draw the solution to the aerodynamics problem.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param aerodynamics_problem: AerodynamicsProblem
        Aerodynamics problem object whose solution will be plotted.
    :param points_type: str
        The name of the object whose points will be plotted.
    :return: None
    """

    # Initialize the plotter.
    plotter = pv.Plotter()

    # Make the airplane geometry.
    vertices = np.vstack((
        aerodynamics_problem.front_left_vertices,
        aerodynamics_problem.front_right_vertices,
        aerodynamics_problem.back_right_vertices,
        aerodynamics_problem.back_left_vertices
    ))
    faces = np.transpose(np.vstack((
        4 * np.ones(aerodynamics_problem.n_panels),
        np.arange(aerodynamics_problem.n_panels),
        np.arange(aerodynamics_problem.n_panels) + aerodynamics_problem.n_panels,
        np.arange(aerodynamics_problem.n_panels) + 2 * aerodynamics_problem.n_panels,
        np.arange(aerodynamics_problem.n_panels) + 3 * aerodynamics_problem.n_panels,
    )))
    faces = np.reshape(faces, (-1), order="C")
    wing_surfaces = pv.PolyData(vertices, faces)

    # If the delta coefficients of pressure have not been added, calculate the delta coefficients of pressure.
    if not hasattr(aerodynamics_problem, "delta_cp"):
        aerodynamics_problem.calculate_delta_cp()

    # Define the minimum and maximum delta coefficients of pressure to scale the panel colors.
    delta_cp_min = -1.5
    delta_cp_max = 1.5

    # Create the scalars, color map, mesh, and scalar bar.
    scalars = np.minimum(np.maximum(aerodynamics_problem.delta_cp, delta_cp_min), delta_cp_max)
    color_map = plt.cm.get_cmap("viridis")
    plotter.add_mesh(wing_surfaces, scalars=scalars, cmap=color_map, color='tan', show_edges=True,
                     smooth_shading=True)
    plotter.add_scalar_bar(title="Pressure Coefficient Differential", n_labels=5, shadow=True,
                           font_family="arial")

    # Check if the points_type has been set.
    if points_type is not None:
        # If so, add points to the plotter.
        points = getattr(aerodynamics_problem, points_type)
        plotter.add_points(points)

    # Set the background color and camera positions.
    plotter.set_background(color="black")
    plotter.show(cpos=(-1, -1, 1), full_screen=True)
