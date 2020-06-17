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
def draw(airplane, show_delta_pressures, show_streamlines):
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

                if show_delta_pressures:
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

    if show_delta_pressures:
        plotter.add_mesh(panel_surface, show_edges=True, cmap=color_map, scalars=scalars, color='white',
                         smooth_shading=True)
    else:
        plotter.add_mesh(panel_surface, show_edges=True, cmap=color_map, color='white', smooth_shading=True)

    if show_streamlines:
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
def make_flapping_gif(movement, show_delta_pressures):
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

                if show_delta_pressures:
                    scalar_to_add = np.maximum(np.minimum(panel.delta_pressure, 100), -100)
                else:
                    scalar_to_add = panel.center[2] / 2

                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))
                scalars = np.hstack((scalars, scalar_to_add))

        for wake_ring_vortex in np.flip(np.ravel(wing.wake_ring_vortices)):
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

    if show_delta_pressures:
        plotter.update_scalar_bar_range(clim=[-100, 100])
    else:
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

        plotter.clear()

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
                    if show_delta_pressures:
                        scalar_to_add = np.maximum(np.minimum(panel.delta_pressure, 100), -100)
                    else:
                        scalar_to_add = panel.center[2] / 2

                    panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                    panel_faces = np.hstack((panel_faces, panel_face_to_add))
                    scalars = np.hstack((scalars, scalar_to_add))

            for wake_ring_vortex in np.flip(np.ravel(wing.wake_ring_vortices)):
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
        if show_delta_pressures:
            plotter.update_scalar_bar_range(clim=[-100, 100])
        else:
            plotter.update_scalar_bar_range(clim=[-2, 2])
        plotter.update_scalars(scalars)

        plotter.write_frame()

    # Close movie and delete object
    plotter.close()


# ToDo: Properly document this function.
def plot_results_versus_time(movement):

    num_steps = movement.num_steps
    delta_time = movement.delta_time
    airplanes = movement.airplanes

    times = np.arange(0, num_steps * delta_time, delta_time)

    total_near_field_force_wind_axes = np.zeros((3, num_steps))
    total_near_field_force_coefficients_wind_axes = np.zeros((3, num_steps))
    total_near_field_moment_wind_axes = np.zeros((3, num_steps))
    total_near_field_moment_coefficients_wind_axes = np.zeros((3, num_steps))

    for step in range(num_steps):

        airplane = airplanes[step]

        total_near_field_force_wind_axes[:, step] = airplane.total_near_field_force_wind_axes
        total_near_field_force_coefficients_wind_axes[:, step] = airplane.total_near_field_force_coefficients_wind_axes
        total_near_field_moment_wind_axes[:, step] = airplane.total_near_field_moment_wind_axes
        total_near_field_moment_coefficients_wind_axes[:, step] = (airplane.total_near_field_moment_coefficients_wind_axes)

    force_figure, force_axes = plt.subplots()
    force_axes.plot(times, total_near_field_force_wind_axes[0], label='Induced Drag')
    force_axes.plot(times, total_near_field_force_wind_axes[1], label='Side Force')
    force_axes.plot(times, total_near_field_force_wind_axes[2], label='Lift')
    force_axes.set_xlabel('Time (s)')
    force_axes.set_ylabel('Force (N)')
    force_axes.set_title('Total Forces in Wind Axes versus Time')
    force_axes.legend()
    force_figure.show()

    force_coefficients_figure, force_coefficients_axes = plt.subplots()
    force_coefficients_axes.plot(
        times, total_near_field_force_coefficients_wind_axes[0], label='Coefficient of Induced Drag'
    )
    force_coefficients_axes.plot(
        times, total_near_field_force_coefficients_wind_axes[1], label='Coefficient of Side Force'
    )
    force_coefficients_axes.plot(
        times, total_near_field_force_coefficients_wind_axes[2], label='Coefficient of Lift'
    )
    force_coefficients_axes.set_xlabel('Time (s)')
    force_coefficients_axes.set_ylabel('Dimensionless')
    force_coefficients_axes.set_title('Total Force Coefficients in Wind Axes versus Time')
    force_coefficients_axes.legend()
    force_coefficients_figure.show()

    moment_figure, moment_axes = plt.subplots()
    moment_axes.plot(times, total_near_field_moment_wind_axes[0], label='Rolling Moment')
    moment_axes.plot(times, total_near_field_moment_wind_axes[1], label='Pitching Moment')
    moment_axes.plot(times, total_near_field_moment_wind_axes[2], label='Yawing Moment')
    moment_axes.set_xlabel('Time (s)')
    moment_axes.set_ylabel('Moment (Nm)')
    moment_axes.set_title('Total Moments in Wind Axes versus Time')
    moment_axes.legend()
    moment_figure.show()

    moment_coefficients_figure, moment_coefficients_axes = plt.subplots()
    moment_coefficients_axes.plot(
        times, total_near_field_moment_coefficients_wind_axes[0], label='Coefficient of Rolling Moment'
    )
    moment_coefficients_axes.plot(
        times, total_near_field_moment_coefficients_wind_axes[1], label='Coefficient of Pitching Moment'
    )
    moment_coefficients_axes.plot(
        times, total_near_field_moment_coefficients_wind_axes[2], label='Coefficient of Yawing Moment'
    )
    moment_coefficients_axes.set_xlabel('Time (s)')
    moment_coefficients_axes.set_ylabel('Dimensionless')
    moment_coefficients_axes.set_title('Total Moment Coefficients in Wind Axes versus Time')
    moment_coefficients_axes.legend()
    moment_coefficients_figure.show()
