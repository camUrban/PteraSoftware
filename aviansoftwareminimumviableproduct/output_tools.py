from aviansoftwareminimumviableproduct import mesh_tools
import numpy as np
import pyvista as pv


# Adapted from:         vlm3.make_panels() in AeroSandbox
# Author:               Peter Sharpe
# Date of Retrieval:    03/28/2020
def draw(aerodynamics_problem, points_type=None):
    if aerodynamics_problem.verbose:
        print("Drawing...")

    # Initialize Plotter
    plotter = pv.Plotter()

    # Make airplane geometry
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
    faces = np.reshape(faces, (-1), order='C')
    wing_surfaces = pv.PolyData(vertices, faces)

    plotter.add_mesh(wing_surfaces, color='white', show_edges=True,
                     smooth_shading=True)

    # Points
    if points_type is not None:
        points = getattr(aerodynamics_problem, points_type)

        plotter.add_points(points)

    # Do the plotting
    plotter.set_background(color="black")
    plotter.show(cpos=(-1, -1, 1), full_screen=True)

    if aerodynamics_problem.verbose:
        print("Drawing finished!")


# ToDo: Properly comment and cite the code in the following method
def make_gif(aerodynamics_problem):
    plotter = pv.Plotter()

    for i in range(len(aerodynamics_problem.time_steps)):
        aerodynamics_problem.current_time = aerodynamics_problem.time_steps[i]
        aerodynamics_problem.current_sweep_angle = aerodynamics_problem.sweep_angles[i]
        if aerodynamics_problem.current_time == 0:
            mesh_tools.make_panels(aerodynamics_problem)

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
            faces = np.reshape(faces, (-1), order='C')
            wing_surfaces = pv.PolyData(vertices, faces)

            plotter.add_mesh(wing_surfaces, color='white', show_edges=True,
                             smooth_shading=False)

            print("Press \"q\" to close window and create an animation")

            plotter.show(auto_close=False)
            plotter.open_gif("flapping.gif")

            aerodynamics_problem.initial_front_left_vertices = aerodynamics_problem.front_left_vertices
            aerodynamics_problem.initial_front_right_vertices = aerodynamics_problem.front_right_vertices
            aerodynamics_problem.initial_back_left_vertices = aerodynamics_problem.back_left_vertices
            aerodynamics_problem.initial_back_right_vertices = aerodynamics_problem.back_right_vertices

            plotter.write_frame()
            plotter.clear()
        else:
            aerodynamics_problem.move_panels()

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
            faces = np.reshape(faces, (-1), order='C')
            wing_surfaces = pv.PolyData(vertices, faces)

            plotter.add_mesh(wing_surfaces, color='white', show_edges=True,
                             smooth_shading=False)

            plotter.write_frame()
            plotter.clear()

    # Close movie and delete object
    plotter.close()
