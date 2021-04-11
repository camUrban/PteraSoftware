"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the geometry of an airplane object.
    animate: Create an animation of a problem's movement.
    plot_results_versus_time: This method takes in an unsteady solver object,
    and plots the geometries' forces, moments,
                              and coefficients as a function of time.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

import pterasoftware as ps


def draw(
    solver,
    show_delta_pressures=False,
    show_streamlines=False,
    show_wake_vortices=False,
):
    """Draw the geometry of an airplane object.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry and attributes are to be plotted.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels.
        The default value is False.
    :param show_streamlines: bool, optional
        Set this variable to true to show the streamlines emanating from the back of
        the wings. The default value is
        False.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices.
        The default value is False.
    :return: None
    """

    # Initialize the plotter and get the plasma color map.
    plotter = pv.Plotter()
    color_map = plt.cm.get_cmap("plasma")

    # Get the solver's geometry.
    if isinstance(
        solver,
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        airplane = solver.steady_problems[-1].airplane
    else:
        airplane = solver.airplane

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

            # Stack this panel's vertices, faces, and scalars. Look through the
            # PolyData documentation for more
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

            # Stack this panel's vertices, faces, and scalars with the array of all
            # the vertices, faces, and
            # scalars.
            panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
            panel_faces = np.hstack((panel_faces, panel_face_to_add))

            # If the user wants to plot the pressures, add the panel's delta pressure
            # to the array of scalars.
            if show_delta_pressures:
                scalar_to_add = panel.delta_pressure
                scalars = np.hstack((scalars, scalar_to_add))

            # Update the number of previous panels.
            current_panel_num += 1

        # Check if the user wants to show the wake vortices.
        if show_wake_vortices:

            wake_ring_vortex_vertices = np.empty((0, 3), dtype=int)
            wake_ring_vortex_faces = np.empty(0, dtype=int)
            current_wake_ring_vortex_num = 0

            # Iterate through the unraveled array of wake vortices for the given wing.
            for wake_ring_vortex in np.ravel(wing.wake_ring_vortices):
                wake_ring_vortex_vertices_to_add = np.vstack(
                    (
                        wake_ring_vortex.front_left_vertex,
                        wake_ring_vortex.front_right_vertex,
                        wake_ring_vortex.back_right_vertex,
                        wake_ring_vortex.back_left_vertex,
                    )
                )
                wake_ring_vortex_face_to_add = np.array(
                    [
                        4,
                        (current_wake_ring_vortex_num * 4),
                        (current_wake_ring_vortex_num * 4) + 1,
                        (current_wake_ring_vortex_num * 4) + 2,
                        (current_wake_ring_vortex_num * 4) + 3,
                    ]
                )

                wake_ring_vortex_vertices = np.vstack(
                    (wake_ring_vortex_vertices, wake_ring_vortex_vertices_to_add)
                )
                wake_ring_vortex_faces = np.hstack(
                    (wake_ring_vortex_faces, wake_ring_vortex_face_to_add)
                )

                current_wake_ring_vortex_num += 1

            wake_ring_vortex_surface = pv.PolyData(
                wake_ring_vortex_vertices, wake_ring_vortex_faces
            )

            plotter.add_mesh(
                wake_ring_vortex_surface,
                show_edges=True,
                smooth_shading=True,
                color="white",
            )

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the
    # plotter with the pressure scalars.
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

        # Iterate through the spanwise positions in the solver's streamline point
        # matrix.
        for spanwise_position in range(solver.streamline_points.shape[1]):

            # Get the column of streamline points at this spanwise position.
            streamline_point_column = solver.streamline_points[:, spanwise_position, :]

            # Iterate through each streamline point column.
            for point_index in range(streamline_point_column.shape[0]):

                # Skip the first point because it has not previous point with which
                # to make a line.
                if point_index != 0:
                    # Get the current, and the last point.
                    point = streamline_point_column[point_index, :]
                    last_point = streamline_point_column[point_index - 1, :]

                    # Add a line to make this segment of the streamline.
                    plotter.add_mesh(
                        pv.Line(
                            last_point,
                            point,
                        ),
                        show_edges=True,
                        color="#EEEEEF",
                        line_width=2,
                    )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="#000000")
    plotter.show(cpos=(-1, -1, 1), full_screen=False)


def animate(
    unsteady_solver,
    show_delta_pressures=False,
    show_wake_vortices=False,
    keep_file=True,
):
    """Create an animation of a solver's geometry.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry is to be animated.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels.
        The default value is false.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices.
        The default value is false.
    :param keep_file: bool, optional
        Set this variable to false in order to not save the resulting GIF. The
        default value is true.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

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

            # Stack this panel's vertices, faces, and scalars. Look through the
            # PolyData documentation for more
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

            # If the user wants to plot the pressures, add the panel's delta pressure
            # to the array of scalars.
            if show_delta_pressures:
                scalar_to_add = panel.delta_pressure
                scalars = np.hstack((scalars, scalar_to_add))

            # Stack this panel's vertices, faces, and scalars with the array of all
            # the vertices, faces, and
            # scalars.
            panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
            panel_faces = np.hstack((panel_faces, panel_face_to_add))

            # Update the number of previous panels.
            current_panel_num += 1

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surface = pv.PolyData(panel_vertices, panel_faces)

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the
    # plotter with the pressure scalars.
    # Otherwise, add the panel surfaces without the pressure scalars.
    if show_delta_pressures and first_results_step == 0:
        plotter.add_mesh(
            panel_surface,
            show_edges=True,
            cmap=color_map,
            scalars=scalars,
            smooth_shading=True,
        )
        # As this is the first step with pressures update the scalar bar range (so it
        # is not automatically set each frame), and update the scalars after doing so.
        if show_delta_pressures:
            plotter.update_scalar_bar_range(clim=[-1000, 1000])
            plotter.update_scalars(scalars)
    else:
        plotter.add_mesh(
            panel_surface,
            show_edges=True,
            cmap=color_map,
            color="#86C552",
            smooth_shading=True,
        )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color="#000000")

    # Print a message to the console on how to set up the window.
    print(
        'Orient the view, then press "q" to close the window and produce the '
        "animation."
    )

    # Set up the camera and close the window.
    plotter.show(cpos=(-1, -1, 1), full_screen=False, auto_close=False)

    # Open a gif.
    plotter.open_gif("animation.gif")

    # Initialize a variable to keep track of which step we are on.
    current_step = 1

    # Begin to iterate through all the other airplanes.
    for airplane in airplanes[1:]:

        # Clear the plotter.
        plotter.clear()

        # Initialize empty ndarrays to hold the things to plot.
        panel_vertices = np.empty((0, 3), dtype=int)
        panel_faces = np.empty(0, dtype=int)
        scalars = np.empty(0, dtype=int)

        # Initialize a variable to keep track of how many panels have been added thus
        # far.
        current_panel_num = 0

        # Increment through the current airplane's wings.
        for wing in airplane.wings:

            # Unravel the wing's panel matrix and iterate through it.
            panels = np.ravel(wing.panels)
            for panel in panels:

                # Stack this panel's vertices, faces, and scalars. Look through the
                # PolyData documentation for more
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

                # Stack this panel's vertices, faces, and scalars with the array of
                # all the vertices, faces, and
                # scalars.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # If the user wants to plot the pressures, add the panel's delta
                # pressure to the array of scalars.
                if show_delta_pressures:
                    scalar_to_add = panel.delta_pressure
                    scalars = np.hstack((scalars, scalar_to_add))

                # Update the number of previous panels.
                current_panel_num += 1

            # Check if the user wants to show the wake vortices.
            if show_wake_vortices:

                wake_ring_vortex_vertices = np.empty((0, 3), dtype=int)
                wake_ring_vortex_faces = np.empty(0, dtype=int)
                current_wake_ring_vortex_num = 0

                # Iterate through the unraveled array of wake vortices for the given
                # wing.
                for wake_ring_vortex in np.ravel(wing.wake_ring_vortices):
                    wake_ring_vortex_vertices_to_add = np.vstack(
                        (
                            wake_ring_vortex.front_left_vertex,
                            wake_ring_vortex.front_right_vertex,
                            wake_ring_vortex.back_right_vertex,
                            wake_ring_vortex.back_left_vertex,
                        )
                    )
                    wake_ring_vortex_face_to_add = np.array(
                        [
                            4,
                            (current_wake_ring_vortex_num * 4),
                            (current_wake_ring_vortex_num * 4) + 1,
                            (current_wake_ring_vortex_num * 4) + 2,
                            (current_wake_ring_vortex_num * 4) + 3,
                        ]
                    )

                    wake_ring_vortex_vertices = np.vstack(
                        (wake_ring_vortex_vertices, wake_ring_vortex_vertices_to_add)
                    )
                    wake_ring_vortex_faces = np.hstack(
                        (wake_ring_vortex_faces, wake_ring_vortex_face_to_add)
                    )

                    current_wake_ring_vortex_num += 1

                wake_ring_vortex_surface = pv.PolyData(
                    wake_ring_vortex_vertices, wake_ring_vortex_faces
                )

                plotter.add_mesh(
                    wake_ring_vortex_surface,
                    show_edges=True,
                    smooth_shading=True,
                    color="white",
                )

        # Initialize the panel surfaces and add the meshes to the plotter.
        panel_surface = pv.PolyData(panel_vertices, panel_faces)

        # Check if the user wants to plot pressures and this step is equal to or
        # greater than the first step with calculated results. If so, add the panel
        # surfaces to the plotter with the pressure scalars. Otherwise, add the panel
        # surfaces without the pressure scalars.
        if show_delta_pressures and first_results_step <= current_step:
            plotter.add_mesh(
                panel_surface,
                show_edges=True,
                cmap=color_map,
                scalars=scalars,
                smooth_shading=True,
            )
            if first_results_step == current_step:

                # If this is the first step with pressures update the scalar bar
                # range (so it is not automatically set each frame), and update the
                # scalars after doing so.
                if show_delta_pressures:
                    plotter.update_scalar_bar_range(clim=[-1000, 1000])
                    plotter.update_scalars(scalars)
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
        plotter.clear()

        # Increment the step number tracker.
        current_step += 1

    # Close the animation and delete the plotter.
    plotter.close()

    # Delete the file if requested.
    if not keep_file:
        os.remove("animation.gif")


def plot_results_versus_time(unsteady_solver, testing=False):
    """This method takes in an unsteady solver object, and plots the geometries'
    forces, moments, and coefficients as a
    function of time.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose resulting forces, moments, and coefficients
        are to be plotted.
    :param testing: bool, Optional
        This boolean determines if the plots will be shown. If true, no plots will be
        shown. It is useful for testing,
        where the user wants to know that the plots were created without having to
        show them. It's default value is
        false.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

    # Get this solver's time step characteristics. Note that the first time step,
    # time step 0, occurs at 0 seconds.
    num_steps = unsteady_solver.num_steps
    delta_time = unsteady_solver.delta_time
    first_results_time_step_time = delta_time * first_results_step
    final_time_step_time = delta_time * (num_steps - 1)
    num_steps_with_results = num_steps - first_results_step

    # Create a 1D array with the time at each time step where results have been
    # calculated.
    times = np.linspace(
        first_results_time_step_time,
        final_time_step_time,
        num_steps_with_results,
        endpoint=True,
    )

    # Initialize matrices to hold the forces, moments, and coefficients at each time
    # step which has results.
    total_near_field_force_wind_axes = np.zeros((3, num_steps_with_results))
    total_near_field_force_coefficients_wind_axes = np.zeros(
        (3, num_steps_with_results)
    )
    total_near_field_moment_wind_axes = np.zeros((3, num_steps_with_results))
    total_near_field_moment_coefficients_wind_axes = np.zeros(
        (3, num_steps_with_results)
    )

    # Initialize a variable to track position in the results arrays.
    results_step = 0

    # Iterate through the time steps and add the results to their respective matrices.
    for step in range(first_results_step, num_steps):
        # Get the airplane from the problem at this step.
        airplane = unsteady_solver.steady_problems[step].airplane

        total_near_field_force_wind_axes[
            :, results_step
        ] = airplane.total_near_field_force_wind_axes
        total_near_field_force_coefficients_wind_axes[
            :, results_step
        ] = airplane.total_near_field_force_coefficients_wind_axes
        total_near_field_moment_wind_axes[
            :, results_step
        ] = airplane.total_near_field_moment_wind_axes
        total_near_field_moment_coefficients_wind_axes[
            :, results_step
        ] = airplane.total_near_field_moment_coefficients_wind_axes

        results_step += 1

    # Create and show the force plot.
    # Initialize the plot.
    force_figure, force_axes = plt.subplots()

    # Add each of the three components of the force.
    force_axes.plot(
        times,
        total_near_field_force_wind_axes[0],
        label="Induced Drag",
        color="#000000",
    )
    force_axes.plot(
        times, total_near_field_force_wind_axes[1], label="Side Force", color="#86C552"
    )
    force_axes.plot(
        times, total_near_field_force_wind_axes[2], label="Lift", color="#E62128"
    )

    # Name the axis labels and the title.
    force_axes.set_xlabel("Time (s)", color="#000000")
    force_axes.set_ylabel("Force (N)", color="#000000")
    force_axes.set_title("Total Forces in Wind Axes versus Time", color="#000000")

    # Set the plot's background color.
    force_figure.patch.set_facecolor("#EEEEEF")
    force_axes.set_facecolor("#EEEEEF")

    # Add a legend.
    force_axes.legend(facecolor="#EEEEEF")

    # Show the plot.
    if not testing:
        force_figure.show()

    # Create and show the force coefficient plot.
    # Initialize the plot.
    force_coefficients_figure, force_coefficients_axes = plt.subplots()

    # Add each of the three force coefficients.
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[0],
        label="Coefficient of Induced Drag",
        color="#000000",
    )
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[1],
        label="Coefficient of Side Force",
        color="#86C552",
    )
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[2],
        label="Coefficient of Lift",
        color="#E62128",
    )

    # Name the axis labels and the title.
    force_coefficients_axes.set_xlabel("Time (s)", color="#000000")
    force_coefficients_axes.set_ylabel("Dimensionless", color="#000000")
    force_coefficients_axes.set_title(
        "Total Force Coefficients in Wind Axes versus Time", color="#000000"
    )

    # Set the plot's background color.
    force_coefficients_figure.patch.set_facecolor("#EEEEEF")
    force_coefficients_axes.set_facecolor("#EEEEEF")

    # Add a legend.
    force_coefficients_axes.legend(facecolor="#EEEEEF")

    # Show the plot.
    if not testing:
        force_coefficients_figure.show()

    # Create and show the moment plot.
    # Initialize the plot.
    moment_figure, moment_axes = plt.subplots()

    # Add each of the three components of the moment.
    moment_axes.plot(
        times,
        total_near_field_moment_wind_axes[0],
        label="Rolling Moment",
        color="#000000",
    )
    moment_axes.plot(
        times,
        total_near_field_moment_wind_axes[1],
        label="Pitching Moment",
        color="#86C552",
    )
    moment_axes.plot(
        times,
        total_near_field_moment_wind_axes[2],
        label="Yawing Moment",
        color="#E62128",
    )

    # Name the axis labels and the title.
    moment_axes.set_xlabel("Time (s)", color="#000000")
    moment_axes.set_ylabel("Moment (Nm)", color="#000000")
    moment_axes.set_title("Total Moments in Wind Axes versus Time", color="#000000")

    # Set the plot's background color.
    moment_figure.patch.set_facecolor("#EEEEEF")
    moment_axes.set_facecolor("#EEEEEF")

    # Add a legend.
    moment_axes.legend(facecolor="#EEEEEF")

    # Show the plot.
    if not testing:
        moment_figure.show()

    # Create and show the moment coefficient plot.
    # Initialize the plot.
    moment_coefficients_figure, moment_coefficients_axes = plt.subplots()

    # Add each of the three moment coefficients.
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[0],
        label="Coefficient of Rolling Moment",
        color="#000000",
    )
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[1],
        label="Coefficient of Pitching Moment",
        color="#86C552",
    )
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[2],
        label="Coefficient of Yawing Moment",
        color="#E62128",
    )

    # Name the axis labels and the title.
    moment_coefficients_axes.set_xlabel("Time (s)", color="#000000")
    moment_coefficients_axes.set_ylabel("Dimensionless", color="#000000")
    moment_coefficients_axes.set_title(
        "Total Moment Coefficients in Wind Axes versus Time", color="#000000"
    )

    # Set the plot's background color.
    moment_coefficients_figure.patch.set_facecolor("#EEEEEF")
    moment_coefficients_axes.set_facecolor("#EEEEEF")

    # Add a legend.
    moment_coefficients_axes.legend(facecolor="#EEEEEF")

    # Show the plot.
    if not testing:
        moment_coefficients_figure.show()
