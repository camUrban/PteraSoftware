"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the geometry of an current_airplane object.
    draw_geometry: Draw the geometry of an current_airplane
"""

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt


# ToDo: Properly document this function.
def draw(airplane, show_delta_pressures):
    """Draw the geometry of an current_airplane object.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param airplane: Airplane
        This is the current_airplane object whose geometry is to be plotted.
    :param show_delta_pressures: bool
        Set this variable to true to show the change in pressure across the panels.
    :return: None
    """

    # Initialize the plotter.
    plotter = pv.Plotter()

    # Set the color map.
    color_map = plt.cm.get_cmap('plasma')

    # Initialize empty ndarrays to hold the things to plot.
    panel_vertices = np.empty((0, 3))
    panel_faces = np.empty(0)
    scalars = np.empty(0)

    # Increment through the current_airplane's wings.
    for wing in airplane.wings:
        # Increment through the wing's chordwise and spanwise positions.
        for chordwise_position in range(wing.num_chordwise_panels):
            for spanwise_position in range(wing.num_spanwise_panels):
                # Calculate the panel number, starting from zero.
                panel_num = (chordwise_position * wing.num_spanwise_panels) + spanwise_position

                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]

                # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
                # details.
                panel_vertices_to_add = np.vstack((
                    panel.front_left_vertex,
                    panel.front_right_vertex,
                    panel.back_right_vertex,
                    panel.back_left_vertex
                ))
                panel_face_to_add = np.array([4,
                                              (panel_num * 4),
                                              (panel_num * 4) + 1,
                                              (panel_num * 4) + 2,
                                              (panel_num * 4) + 3])

                # Stack this panel's vertices, faces, and scalars with the ndarray of all the vertices, faces, and
                # scalars.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))
                scalar_to_add = np.maximum(np.minimum(panel.delta_pressure, 1000), -1000)
                scalars = np.hstack((scalars, scalar_to_add))

        for wake_ring_vortex in np.ravel(wing.wake_ring_vortices):
            plotter.add_mesh(pv.Line(wake_ring_vortex.front_right_vertex, wake_ring_vortex.front_left_vertex),
                             show_edges=True, cmap=color_map, color='white')
            plotter.add_mesh(pv.Line(wake_ring_vortex.front_left_vertex, wake_ring_vortex.back_left_vertex),
                             show_edges=True, cmap=color_map, color='white')
            plotter.add_mesh(pv.Line(wake_ring_vortex.back_left_vertex, wake_ring_vortex.back_right_vertex),
                             show_edges=True, cmap=color_map, color='white')
            plotter.add_mesh(pv.Line(wake_ring_vortex.back_right_vertex, wake_ring_vortex.front_right_vertex),
                             show_edges=True, cmap=color_map, color='white')

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)
    plotter.add_mesh(panel_surface, show_edges=True, cmap=color_map, scalars=scalars, color='white',
                     smooth_shading=True)

    for wing in airplane.wings:
        streamline_points = np.reshape(wing.streamline_points, (-1, 1, 3))
        plotter.add_points(pv.PolyData(streamline_points))

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="black")
    plotter.show(cpos=(-1, -1, 1), full_screen=False)


# ToDo: Properly document this function.
def draw_geometry(airplane):
    """Draw the geometry of an current_airplane object.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param airplane: Airplane
        This is the current_airplane object whose geometry is to be plotted.
    :return: None
    """

    # Initialize the plotter.
    plotter = pv.Plotter()

    # Initialize empty ndarrays to hold the things to plot.
    panel_vertices = np.empty((0, 3))
    panel_faces = np.empty(0)

    # Increment through the current_airplane's wings.
    for wing in airplane.wings:
        # Increment through the wing's chordwise and spanwise positions.
        for chordwise_position in range(wing.num_chordwise_panels):
            for spanwise_position in range(wing.num_spanwise_panels):
                # Calculate the panel number, starting from zero.
                panel_num = (chordwise_position * wing.num_spanwise_panels) + spanwise_position

                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]

                # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
                # details.
                panel_vertices_to_add = np.vstack((
                    panel.front_left_vertex,
                    panel.front_right_vertex,
                    panel.back_right_vertex,
                    panel.back_left_vertex
                ))
                panel_face_to_add = np.array([4,
                                              (panel_num * 4),
                                              (panel_num * 4) + 1,
                                              (panel_num * 4) + 2,
                                              (panel_num * 4) + 3])

                # Stack this panel's vertices, faces, and scalars with the ndarray of all the vertices, faces, and
                # scalars.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)
    plotter.add_mesh(panel_surface, show_edges=True, color='white')

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="black")
    plotter.show(cpos=(-1, -1, 1), full_screen=False)


# ToDo: Properly document and cite this function.
def make_flapping_gif(movement):
    airplanes = movement.airplanes

    # Initialize the plotter.
    plotter = pv.Plotter()

    color_map = plt.cm.get_cmap('plasma', 256)

    # Initialize empty ndarrays to hold the things to plot.
    panel_vertices = np.empty((0, 3))
    panel_faces = np.empty(0)
    scalars = np.empty(0)

    # Increment through the current_airplane's wings.
    for wing in airplanes[0].wings:
        # Increment through the wing's chordwise and spanwise positions.
        for chordwise_position in range(wing.num_chordwise_panels):
            for spanwise_position in range(wing.num_spanwise_panels):
                # Calculate the panel number, starting from zero.
                panel_num = (chordwise_position * wing.num_spanwise_panels) + spanwise_position

                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]

                # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
                # details.
                panel_vertices_to_add = np.vstack((
                    panel.front_left_vertex,
                    panel.front_right_vertex,
                    panel.back_right_vertex,
                    panel.back_left_vertex
                ))
                panel_face_to_add = np.array([4,
                                              (panel_num * 4),
                                              (panel_num * 4) + 1,
                                              (panel_num * 4) + 2,
                                              (panel_num * 4) + 3])
                scalar_to_add = panel.center[2] / 2

                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))
                scalars = np.hstack((scalars, scalar_to_add))

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)

    plotter.add_mesh(panel_surface, show_edges=True, scalars=scalars, cmap=color_map, smooth_shading=True)
    plotter.update_scalar_bar_range(clim=[-2, 2])
    plotter.update_scalars(scalars)

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="black")

    print('Orient the view, then press "q" to close window and produce movie')

    # setup camera and close
    plotter.show(cpos=(-1, -1, 1), full_screen=False, auto_close=False)

    # Open a gif
    plotter.open_gif("flapping.gif")

    for airplane in airplanes:
        # Initialize empty ndarrays to hold the things to plot.
        panel_vertices = np.empty((0, 3))
        panel_faces = np.empty(0)
        scalars = np.empty(0)

        # Increment through the current_airplane's wings.
        for wing in airplane.wings:
            # Increment through the wing's chordwise and spanwise positions.
            for chordwise_position in range(wing.num_chordwise_panels):
                for spanwise_position in range(wing.num_spanwise_panels):
                    # Calculate the panel number, starting from zero.
                    panel_num = (chordwise_position * wing.num_spanwise_panels) + spanwise_position

                    # Pull the panel object out of the wing's list of panels.
                    panel = wing.panels[chordwise_position, spanwise_position]

                    # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
                    # details.
                    panel_vertices_to_add = np.vstack((
                        panel.front_left_vertex,
                        panel.front_right_vertex,
                        panel.back_right_vertex,
                        panel.back_left_vertex
                    ))
                    panel_face_to_add = np.array([4,
                                                  (panel_num * 4),
                                                  (panel_num * 4) + 1,
                                                  (panel_num * 4) + 2,
                                                  (panel_num * 4) + 3])
                    scalar_to_add = panel.center[2] / 2

                    panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                    panel_faces = np.hstack((panel_faces, panel_face_to_add))
                    scalars = np.hstack((scalars, scalar_to_add))

        # Initialize the panel surfaces and add the meshes to the plotter.
        panel_surface = pv.PolyData(panel_vertices, panel_faces)
        plotter.clear()

        plotter.add_mesh(panel_surface, show_edges=True, scalars=scalars, cmap=color_map, smooth_shading=True)
        plotter.update_scalar_bar_range(clim=[-2, 2])
        plotter.update_scalars(scalars)

        plotter.write_frame()

    # Close movie and delete object
    plotter.close()
