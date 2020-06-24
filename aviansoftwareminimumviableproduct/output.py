"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the geometry of an airplane object.
    animate: Create an animation of a problem's movement.
"""

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


def draw(
    airplane,
    show_delta_pressures=False,
    show_streamlines=False,
    show_wake_vortices=False,
):
    """Draw the geometry of an airplane object.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param airplane: Airplane
        This is the airplane object whose geometry is to be plotted.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels. The default value is False.
    :param show_streamlines: bool, optional
        Set this variable to true to show the streamlines emanating from the back of the wings. The default value is
        False.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices. The default value is False.
    :return: None
    """

    # Initialize the plotter.
    plotter = pv.Plotter()

    # Set the color map.
    color_map = plt.cm.get_cmap("plasma")

    # Initialize empty ndarrays to hold the things to plot.
    panel_vertices = np.empty((0, 3), dtype=int)
    panel_faces = np.empty(0, dtype=int)
    scalars = np.empty(0, dtype=int)

    # Initialize a variable to keep track of how many panels have been added thus far.
    current_panel_num = 0

    # Increment through the current airplane's wings.
    for wing in airplane.wings:

        # Unravel the wing's panel matrix and iterate through it.
        panels = np.ravel(wing.panels)
        for panel in panels:

            # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
            # details.
            panel_vertices_to_add = np.vstack(
                (
                    panel.front_left_vertex,
                    panel.front_right_vertex,
                    panel.back_right_vertex,
                    panel.back_left_vertex,
                )
            )
            panel_face_to_add = np.array(
                [
                    4,
                    (current_panel_num * 4),
                    (current_panel_num * 4) + 1,
                    (current_panel_num * 4) + 2,
                    (current_panel_num * 4) + 3,
                ]
            )

            # Stack this panel's vertices, faces, and scalars with the ndarray of all the vertices, faces, and
            # scalars.
            panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
            panel_faces = np.hstack((panel_faces, panel_face_to_add))

            # If the user wants to plot the pressures, add the panel's delta pressure to the ndarray of scalars.
            if show_delta_pressures:
                scalar_to_add = panel.delta_pressure
                scalars = np.hstack((scalars, scalar_to_add))

            # Update the number of previous panels.
            current_panel_num += 1

        # Check if the user wants to show the wake vortices.
        if show_wake_vortices:

            # Iterate through the unraveled array of wake vortices for the given wing.
            for wake_ring_vortex in np.ravel(wing.wake_ring_vortices):
                # Add a line to make this wake ring vortex's front edge.
                plotter.add_mesh(
                    pv.Line(
                        wake_ring_vortex.front_right_vertex,
                        wake_ring_vortex.front_left_vertex,
                    ),
                    show_edges=True,
                    color="#E62128",
                    line_width=2,
                )

                # Add a line to make this wake ring vortex's left edge.
                plotter.add_mesh(
                    pv.Line(
                        wake_ring_vortex.front_left_vertex,
                        wake_ring_vortex.back_left_vertex,
                    ),
                    show_edges=True,
                    color="#E62128",
                    line_width=2,
                )

                # Add a line to make this wake ring vortex's back edge.
                plotter.add_mesh(
                    pv.Line(
                        wake_ring_vortex.back_left_vertex,
                        wake_ring_vortex.back_right_vertex,
                    ),
                    show_edges=True,
                    color="#E62128",
                    line_width=2,
                )

                # Add a line to make this wake ring vortex's right edge.
                plotter.add_mesh(
                    pv.Line(
                        wake_ring_vortex.back_right_vertex,
                        wake_ring_vortex.front_right_vertex,
                    ),
                    show_edges=True,
                    color="#E62128",
                    line_width=2,
                )

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the plotter with the pressure scalars.
    # Otherwise, add the panel surfaces without the pressure scalars.
    if show_delta_pressures:
        plotter.add_mesh(
            panel_surface,
            show_edges=True,
            cmap=color_map,
            scalars=scalars,
            smooth_shading=True,
        )
    else:
        plotter.add_mesh(
            panel_surface,
            show_edges=True,
            cmap=color_map,
            color="#86C552",
            smooth_shading=True,
        )

    # Check if the user wants to plot streamlines.
    if show_streamlines:

        # Iterate through the airplane's wings.
        for wing in airplane.wings:

            # Iterate through the wing's spanwise positions, and slice the wing's streamline point matrix into columns.
            for spanwise_position in range(wing.num_spanwise_panels):
                streamline_point_column = wing.streamline_points[
                    :, spanwise_position, :
                ]

                # Iterate through each streamline point column.
                for point_index in range(streamline_point_column.shape[0]):

                    # Skip the first point because it has not previous point with which to make a line.
                    if point_index != 0:
                        # Get the current, and the last point.
                        point = streamline_point_column[point_index, :]
                        last_point = streamline_point_column[point_index - 1, :]

                        # Add a line to make this segment of the streamline.
                        plotter.add_mesh(
                            pv.Line(last_point, point,),
                            show_edges=True,
                            color="#FEEEEE",
                            line_width=2,
                        )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="#000000")
    plotter.show(cpos=(-1, -1, 1), full_screen=False)


def animate(unsteady_solver, show_delta_pressures=False, show_wake_vortices=False):
    """Create an animation of a solver's geometry.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry is to be animated.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels. The default value is false.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices. The default value is false.
    :return: None
    """

    # Get this solver's problem's airplanes.

    airplanes = []

    for steady_problem in unsteady_solver.steady_problems:
        airplanes.append(steady_problem.airplane)

    # Initialize the plotter and get the plasma color map.
    plotter = pv.Plotter()
    color_map = plt.cm.get_cmap("plasma", 256)

    # Initialize empty ndarrays to hold the things to plot.
    panel_vertices = np.empty((0, 3), dtype=int)
    panel_faces = np.empty(0, dtype=int)
    scalars = np.empty(0, dtype=int)

    # Initialize a variable to keep track of how many panels have been added thus far.
    current_panel_num = 0

    # Increment through the current airplane's wings.
    for wing in airplanes[0].wings:

        # Unravel the wing's panel matrix and iterate through it.
        panels = np.ravel(wing.panels)
        for panel in panels:

            # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
            # details.
            panel_vertices_to_add = np.vstack(
                (
                    panel.front_left_vertex,
                    panel.front_right_vertex,
                    panel.back_right_vertex,
                    panel.back_left_vertex,
                )
            )
            panel_face_to_add = np.array(
                [
                    4,
                    (current_panel_num * 4),
                    (current_panel_num * 4) + 1,
                    (current_panel_num * 4) + 2,
                    (current_panel_num * 4) + 3,
                ]
            )

            # If the user wants to plot the pressures, add the panel's delta pressure to the ndarray of scalars.
            if show_delta_pressures:
                scalar_to_add = panel.delta_pressure
                scalars = np.hstack((scalars, scalar_to_add))

            # Stack this panel's vertices, faces, and scalars with the ndarray of all the vertices, faces, and
            # scalars.
            panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
            panel_faces = np.hstack((panel_faces, panel_face_to_add))

            # Update the number of previous panels.
            current_panel_num += 1

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the plotter with the pressure scalars.
    # Otherwise, add the panel surfaces without the pressure scalars.
    if show_delta_pressures:
        plotter.add_mesh(
            panel_surface,
            show_edges=True,
            cmap=color_map,
            scalars=scalars,
            smooth_shading=True,
        )
    else:
        plotter.add_mesh(
            panel_surface,
            show_edges=True,
            cmap=color_map,
            color="#86C552",
            smooth_shading=True,
        )

    # If the user wants to show the pressures, update the scalar bar range (so it is not automatically set each frame),
    # and update the scalars after doing so.
    if show_delta_pressures:
        plotter.update_scalar_bar_range(clim=[-1000, 1000])
        plotter.update_scalars(scalars)

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="black")

    # Print a message to the console on how to set up the window.
    print(
        'Orient the view, then press "q" to close the window and produce the animation.'
    )

    # Set up the camera and close the window.
    plotter.show(cpos=(-1, -1, 1), full_screen=False, auto_close=False)

    # Open a gif.
    plotter.open_gif("animation.gif")

    # Begin to iterate through all the other airplanes.
    for airplane in airplanes[1:]:

        # Clear the plotter.
        plotter.clear()

        # Initialize empty ndarrays to hold the things to plot.
        panel_vertices = np.empty((0, 3), dtype=int)
        panel_faces = np.empty(0, dtype=int)
        scalars = np.empty(0, dtype=int)

        # Initialize a variable to keep track of how many panels have been added thus far.
        current_panel_num = 0

        # Increment through the current airplane's wings.
        for wing in airplane.wings:

            # Unravel the wing's panel matrix and iterate through it.
            panels = np.ravel(wing.panels)
            for panel in panels:

                # Stack this panel's vertices, faces, and scalars. Look through the PolyData documentation for more
                # details.
                panel_vertices_to_add = np.vstack(
                    (
                        panel.front_left_vertex,
                        panel.front_right_vertex,
                        panel.back_right_vertex,
                        panel.back_left_vertex,
                    )
                )
                panel_face_to_add = np.array(
                    [
                        4,
                        (current_panel_num * 4),
                        (current_panel_num * 4) + 1,
                        (current_panel_num * 4) + 2,
                        (current_panel_num * 4) + 3,
                    ]
                )

                # Stack this panel's vertices, faces, and scalars with the ndarray of all the vertices, faces, and
                # scalars.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # If the user wants to plot the pressures, add the panel's delta pressure to the ndarray of scalars.
                if show_delta_pressures:
                    scalar_to_add = panel.delta_pressure
                    scalars = np.hstack((scalars, scalar_to_add))

                # Update the number of previous panels.
                current_panel_num += 1

            # Check if the user wants to show the wake vortices.
            if show_wake_vortices:

                # Iterate through the unraveled array of wake vortices for the given wing.
                for wake_ring_vortex in np.ravel(wing.wake_ring_vortices):
                    # Add a line to make this wake ring vortex's front edge.
                    plotter.add_mesh(
                        pv.Line(
                            wake_ring_vortex.front_right_vertex,
                            wake_ring_vortex.front_left_vertex,
                        ),
                        show_edges=True,
                        color="#E62128",
                        line_width=2,
                    )

                    # Add a line to make this wake ring vortex's left edge.
                    plotter.add_mesh(
                        pv.Line(
                            wake_ring_vortex.front_left_vertex,
                            wake_ring_vortex.back_left_vertex,
                        ),
                        show_edges=True,
                        color="#E62128",
                        line_width=2,
                    )

                    # Add a line to make this wake ring vortex's back edge.
                    plotter.add_mesh(
                        pv.Line(
                            wake_ring_vortex.back_left_vertex,
                            wake_ring_vortex.back_right_vertex,
                        ),
                        show_edges=True,
                        color="#E62128",
                        line_width=2,
                    )

                    # Add a line to make this wake ring vortex's right edge.
                    plotter.add_mesh(
                        pv.Line(
                            wake_ring_vortex.back_right_vertex,
                            wake_ring_vortex.front_right_vertex,
                        ),
                        show_edges=True,
                        color="#E62128",
                        line_width=2,
                    )

        # Initialize the panel surfaces and add the meshes to the plotter.
        panel_surface = pv.PolyData(panel_vertices, panel_faces)

        # Check if the user wants to plot pressures. If so, add the panel surfaces to the plotter with the pressure scalars.
        # Otherwise, add the panel surfaces without the pressure scalars.
        if show_delta_pressures:
            plotter.add_mesh(
                panel_surface,
                show_edges=True,
                cmap=color_map,
                scalars=scalars,
                smooth_shading=True,
            )
        else:
            plotter.add_mesh(
                panel_surface,
                show_edges=True,
                cmap=color_map,
                color="#86C552",
                smooth_shading=True,
            )

        # Write the current frame to the plotter.
        plotter.write_frame()

    # Close the animation and delete the plotter.
    plotter.close()


# ToDo: Properly document this function.
def plot_results_versus_time(movement, verbose=True):
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

        total_near_field_force_wind_axes[
            :, step
        ] = airplane.total_near_field_force_wind_axes
        total_near_field_force_coefficients_wind_axes[
            :, step
        ] = airplane.total_near_field_force_coefficients_wind_axes
        total_near_field_moment_wind_axes[
            :, step
        ] = airplane.total_near_field_moment_wind_axes
        total_near_field_moment_coefficients_wind_axes[
            :, step
        ] = airplane.total_near_field_moment_coefficients_wind_axes

    force_figure, force_axes = plt.subplots()
    force_axes.plot(times, total_near_field_force_wind_axes[0], label="Induced Drag")
    force_axes.plot(times, total_near_field_force_wind_axes[1], label="Side Force")
    force_axes.plot(times, total_near_field_force_wind_axes[2], label="Lift")
    force_axes.set_xlabel("Time (s)")
    force_axes.set_ylabel("Force (N)")
    force_axes.set_title("Total Forces in Wind Axes versus Time")
    force_axes.legend()
    if verbose:
        force_figure.show()

    force_coefficients_figure, force_coefficients_axes = plt.subplots()
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[0],
        label="Coefficient of Induced Drag",
    )
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[1],
        label="Coefficient of Side Force",
    )
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[2],
        label="Coefficient of Lift",
    )
    force_coefficients_axes.set_xlabel("Time (s)")
    force_coefficients_axes.set_ylabel("Dimensionless")
    force_coefficients_axes.set_title(
        "Total Force Coefficients in Wind Axes versus Time"
    )
    force_coefficients_axes.legend()
    if verbose:
        force_coefficients_figure.show()

    moment_figure, moment_axes = plt.subplots()
    moment_axes.plot(
        times, total_near_field_moment_wind_axes[0], label="Rolling Moment"
    )
    moment_axes.plot(
        times, total_near_field_moment_wind_axes[1], label="Pitching Moment"
    )
    moment_axes.plot(times, total_near_field_moment_wind_axes[2], label="Yawing Moment")
    moment_axes.set_xlabel("Time (s)")
    moment_axes.set_ylabel("Moment (Nm)")
    moment_axes.set_title("Total Moments in Wind Axes versus Time")
    moment_axes.legend()
    if verbose:
        moment_figure.show()

    moment_coefficients_figure, moment_coefficients_axes = plt.subplots()
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[0],
        label="Coefficient of Rolling Moment",
    )
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[1],
        label="Coefficient of Pitching Moment",
    )
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[2],
        label="Coefficient of Yawing Moment",
    )
    moment_coefficients_axes.set_xlabel("Time (s)")
    moment_coefficients_axes.set_ylabel("Dimensionless")
    moment_coefficients_axes.set_title(
        "Total Moment Coefficients in Wind Axes versus Time"
    )
    moment_coefficients_axes.legend()
    if verbose:
        moment_coefficients_figure.show()
