# ToDo: Update this module's documentation.
"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the geometry of an airplane object.

    animate: Create an animation of a problem's movement.

    plot_results_versus_time: This method takes in an unsteady solver object,
    and plots the geometries' forces, moments, and coefficients as a function of
    time.

    get_wake_ring_vortex_surfaces: This function returns the PolyData object for the
    surface of wake ring vortices at a given time step.

    get_scalars: This function gets the delta pressure values from a problem's
    airplane objects, and puts them into a 1D array to be used as scalars for display
    by other output methods.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from . import unsteady_ring_vortex_lattice_method

sequential_color_map = "speed"
diverging_color_map = "delta"
wake_vortex_color = "white"
panel_color = "chartreuse"
streamline_color = "orchid"
plotter_background_color = "black"
figure_background_color = "black"
figure_text_color = "white"

# For the figure lines, use the "Prism" qualitative color map from
# carto.com/carto-colors.
prism = [
    "#5F4690",
    "#1D6996",
    "#38A6A5",
    "#0F8554",
    "#73AF48",
    "#EDAD08",
    "#E17C05",
    "#CC503E",
    "#94346E",
    "#6F4070",
    "#994E95",
    "#666666",
]
[
    drag_color,
    side_color,
    lift_color,
    roll_color,
    pitch_color,
    yaw_color,
] = prism[3:9]

# Set the number of markers and the marker size for the results plots.
num_markers = 6
marker_size = 8

# Calculate the normalized spacing between the markers for the results plots.
marker_spacing = 1.0 / num_markers


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

    # Initialize the plotter.
    plotter = pv.Plotter()

    # Get the solver's geometry.
    if isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        airplanes = solver.steady_problems[-1].airplanes
        last_step = solver.num_steps - 1

        # If the user wants to show the wake ring vortices, then get their surfaces and
        # plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = get_wake_ring_vortex_surfaces(solver, last_step)
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=True,
                color=wake_vortex_color,
            )

    else:
        airplanes = solver.airplanes

    # Get the panel surfaces.
    panel_surfaces = get_panel_surfaces(airplanes)

    # Check if the user wants to plot pressures.
    if show_delta_pressures:

        # Get the scalars
        scalars = get_scalars(airplanes)

        # Choose the color map and set its limits based on if the min and max scalars
        # have the same sign (sequential color map) or if they have different signs
        # (diverging color map).
        if np.sign(np.min(scalars)) == np.sign(np.max(scalars)):
            color_map = sequential_color_map
            c_min = max(np.mean(scalars) - 2 * np.std(scalars), np.min(scalars))
            c_max = min(np.mean(scalars) + 2 * np.std(scalars), np.max(scalars))
        else:
            color_map = diverging_color_map
            c_min = -2 * np.std(scalars)
            c_max = 2 * np.std(scalars)

        # Add the panel surfaces to the plotter with the pressure scalars.
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            cmap=color_map,
            scalars=scalars,
            clim=([c_min, c_max]),
            smooth_shading=True,
            show_scalar_bar=False,
        )
        plotter.add_scalar_bar(
            title="Lifting Pressure (Pa)",
            bold=False,
            title_font_size=16,
            label_font_size=14,
            width=0.5,
            position_x=0.25,
            position_y=0.05,
            n_labels=2,
            italic=True,
            fmt="%.1f",
        )
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
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
                        color=streamline_color,
                        line_width=2,
                    )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color=plotter_background_color)
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

    # Initialize the plotter and get the color map.
    plotter = pv.Plotter()

    # Initialize values to hold the color map choice and its limits.
    c_min = 0
    c_max = 0
    color_map = None

    # Check if the user wants to show pressures.
    if show_delta_pressures:

        # Initialize an empty array to hold all of the problem's scalars.
        all_scalars = np.empty(0, dtype=int)

        # Now iterate through all the time steps to get all of the scalars. These
        # values will be used to configure the color map.
        for airplane in airplanes:
            scalars_to_add = get_scalars(airplane)
            all_scalars = np.hstack((all_scalars, scalars_to_add))

        # Choose the color map and set its limits based on if the min and max scalars
        # across all time steps have the same sign (sequential color map) or if they
        # have different signs (diverging color map).
        if np.sign(np.min(all_scalars)) == np.sign(np.max(all_scalars)):
            color_map = sequential_color_map
            c_min = max(
                np.mean(all_scalars) - 2 * np.std(all_scalars), np.min(all_scalars)
            )
            c_max = min(
                np.mean(all_scalars) + 2 * np.std(all_scalars), np.max(all_scalars)
            )
        else:
            color_map = diverging_color_map
            c_min = -2 * np.std(all_scalars)
            c_max = 2 * np.std(all_scalars)

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surfaces = get_panel_surfaces(airplanes[0])

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the
    # plotter with the pressure scalars. Otherwise, add the panel surfaces without
    # the pressure scalars.
    if show_delta_pressures and first_results_step == 0:
        scalars = get_scalars(airplanes[0])
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            cmap=color_map,
            clim=[c_min, c_max],
            scalars=scalars,
            smooth_shading=True,
            show_scalar_bar=False,
        )
        plotter.add_scalar_bar(
            title="Lifting Pressure (Pa)",
            bold=False,
            title_font_size=16,
            label_font_size=14,
            width=0.5,
            position_x=0.25,
            position_y=0.05,
            n_labels=2,
            italic=True,
            fmt="%.1f",
        )

        # Update the scalars is the user wants to show the pressures.
        plotter.update_scalars(scalars)
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
            smooth_shading=True,
        )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color=plotter_background_color)

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

        # Get the panel surfaces.
        panel_surfaces = get_panel_surfaces(airplane)

        # If the user wants to show the wake ring vortices, then get their surfaces and
        # plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = get_wake_ring_vortex_surfaces(
                unsteady_solver, current_step
            )
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=True,
                color=wake_vortex_color,
            )

        # Check if the user wants to plot pressures and this step is equal to or
        # greater than the first step with calculated results. If so, add the panel
        # surfaces to the plotter wit  h the pressure scalars. Otherwise, add the
        # panel surfaces without the pressure scalars.
        if show_delta_pressures and first_results_step <= current_step:
            scalars = get_scalars(airplane)
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                cmap=color_map,
                clim=[c_min, c_max],
                scalars=scalars,
                smooth_shading=True,
                show_scalar_bar=False,
            )
            plotter.add_scalar_bar(
                title="Lifting Pressure (Pa)",
                bold=False,
                title_font_size=16,
                label_font_size=14,
                width=0.5,
                position_x=0.25,
                position_y=0.05,
                n_labels=2,
                italic=True,
                fmt="%.1f",
            )

            if first_results_step == current_step:
                plotter.update_scalars(scalars)
        else:
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                color=panel_color,
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
    forces, moments, and coefficients as a function of time.

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

    # Remove the top and right spines.
    force_axes.spines.right.set_visible(False)
    force_axes.spines.top.set_visible(False)
    force_axes.spines.bottom.set_color(figure_text_color)
    force_axes.spines.left.set_color(figure_text_color)
    force_axes.xaxis.label.set_color(figure_text_color)
    force_axes.yaxis.label.set_color(figure_text_color)
    force_axes.tick_params(axis="x", colors=figure_text_color)
    force_axes.tick_params(axis="y", colors=figure_text_color)

    # Add each of the three components of the force.
    force_axes.plot(
        times,
        total_near_field_force_wind_axes[0],
        label="$\it{Induced\ Drag}$",
        color=drag_color,
        marker=".",
        markevery=(marker_spacing * 2 / 3, marker_spacing),
        markersize=marker_size,
    )
    force_axes.plot(
        times,
        total_near_field_force_wind_axes[1],
        label="$\it{Side\ Force}$",
        color=side_color,
        marker=".",
        markevery=(marker_spacing * 2 / 3, marker_spacing),
        markersize=marker_size,
    )
    force_axes.plot(
        times,
        total_near_field_force_wind_axes[2],
        label="$\it{Lift}$",
        color=lift_color,
        marker=".",
        markevery=(marker_spacing * 2 / 3, marker_spacing),
        markersize=marker_size,
    )

    # Name the axis labels and the title.
    force_axes.set_xlabel("$\it{Time\ (s)}$", color=figure_text_color)
    force_axes.set_ylabel("$\it{Force\ (N)}$", color=figure_text_color)
    force_axes.set_title("$\it{Forces\ vs.\ Time}$", color=figure_text_color)

    # Set the plot's background color.
    force_figure.patch.set_facecolor(figure_background_color)
    force_axes.set_facecolor(figure_background_color)

    # Add a legend.
    force_axes.legend(
        facecolor=figure_background_color,
        edgecolor=figure_background_color,
        labelcolor=figure_text_color,
    )

    # Show the plot.
    if not testing:
        force_figure.show()

    # Create and show the force coefficient plot.
    # Initialize the plot.
    force_coefficients_figure, force_coefficients_axes = plt.subplots()

    # Remove the top and right spines.
    force_coefficients_axes.spines.right.set_visible(False)
    force_coefficients_axes.spines.top.set_visible(False)
    force_coefficients_axes.spines.bottom.set_color(figure_text_color)
    force_coefficients_axes.spines.left.set_color(figure_text_color)
    force_coefficients_axes.xaxis.label.set_color(figure_text_color)
    force_coefficients_axes.yaxis.label.set_color(figure_text_color)
    force_coefficients_axes.tick_params(axis="x", colors=figure_text_color)
    force_coefficients_axes.tick_params(axis="y", colors=figure_text_color)

    # Add each of the three force coefficients.
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[0],
        label="$\it{Induced\ Drag}$",
        color=drag_color,
        marker=".",
        markevery=(marker_spacing * 0 / 3, marker_spacing),
        markersize=marker_size,
    )
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[1],
        label="$\it{Side\ Force}$",
        color=side_color,
        marker=".",
        markevery=(marker_spacing * 1 / 3, marker_spacing),
        markersize=marker_size,
    )
    force_coefficients_axes.plot(
        times,
        total_near_field_force_coefficients_wind_axes[2],
        label="$\it{Lift}$",
        color=lift_color,
        marker=".",
        markevery=(marker_spacing * 0 / 3, marker_spacing),
        markersize=marker_size,
    )

    # Name the axis labels and the title.
    force_coefficients_axes.set_xlabel("$\it{Time\ (s)}$", color=figure_text_color)
    force_coefficients_axes.set_ylabel(
        "$\it{Coefficient\ (Unitless)}$", color=figure_text_color
    )
    force_coefficients_axes.set_title(
        "$\it{Force\ Coefficients\ vs.\ Time}$",
        color=figure_text_color,
    )

    # Set the plot's background color.
    force_coefficients_figure.patch.set_facecolor(figure_background_color)
    force_coefficients_axes.set_facecolor(figure_background_color)

    # Add a legend.
    force_coefficients_axes.legend(
        facecolor=figure_background_color,
        edgecolor=figure_background_color,
        labelcolor=figure_text_color,
    )

    # Show the plot.
    if not testing:
        force_coefficients_figure.show()

    # Create and show the moment plot.
    # Initialize the plot.
    moment_figure, moment_axes = plt.subplots()

    # Remove the top and right spines.
    moment_axes.spines.right.set_visible(False)
    moment_axes.spines.top.set_visible(False)
    moment_axes.spines.bottom.set_color(figure_text_color)
    moment_axes.spines.left.set_color(figure_text_color)
    moment_axes.xaxis.label.set_color(figure_text_color)
    moment_axes.yaxis.label.set_color(figure_text_color)
    moment_axes.tick_params(axis="x", colors=figure_text_color)
    moment_axes.tick_params(axis="y", colors=figure_text_color)

    # Add each of the three components of the moment.
    moment_axes.plot(
        times,
        total_near_field_moment_wind_axes[0],
        label="$\it{Roll}$",
        color=roll_color,
        marker=".",
        markevery=(marker_spacing * 0 / 3, marker_spacing),
        markersize=marker_size,
    )
    moment_axes.plot(
        times,
        total_near_field_moment_wind_axes[1],
        label="$\it{Pitch}$",
        color=pitch_color,
        marker=".",
        markevery=(marker_spacing * 1 / 3, marker_spacing),
        markersize=marker_size,
    )
    moment_axes.plot(
        times,
        total_near_field_moment_wind_axes[2],
        label="$\it{Yaw}$",
        color=yaw_color,
        marker=".",
        markevery=(marker_spacing * 2 / 3, marker_spacing),
        markersize=marker_size,
    )

    # Name the axis labels and the title.
    moment_axes.set_xlabel("$\it{Time\ (s)}$", color=figure_text_color)
    moment_axes.set_ylabel("$\it{Moment\ (N\ m)}$", color=figure_text_color)
    moment_axes.set_title("$\it{Moments\ vs.\ Time}$", color=figure_text_color)

    # Set the plot's background color.
    moment_figure.patch.set_facecolor(figure_background_color)
    moment_axes.set_facecolor(figure_background_color)

    # Add a legend.
    moment_axes.legend(
        facecolor=figure_background_color,
        edgecolor=figure_background_color,
        labelcolor=figure_text_color,
    )

    # Show the plot.
    if not testing:
        moment_figure.show()

    # Create and show the moment coefficient plot.
    # Initialize the plot.
    moment_coefficients_figure, moment_coefficients_axes = plt.subplots()

    # Remove the top and right spines.
    moment_coefficients_axes.spines.right.set_visible(False)
    moment_coefficients_axes.spines.top.set_visible(False)
    moment_coefficients_axes.spines.bottom.set_color(figure_text_color)
    moment_coefficients_axes.spines.left.set_color(figure_text_color)
    moment_coefficients_axes.xaxis.label.set_color(figure_text_color)
    moment_coefficients_axes.yaxis.label.set_color(figure_text_color)
    moment_coefficients_axes.tick_params(axis="x", colors=figure_text_color)
    moment_coefficients_axes.tick_params(axis="y", colors=figure_text_color)

    # Add each of the three moment coefficients.
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[0],
        label="$\it{Roll}$",
        color=roll_color,
        marker=".",
        markevery=(marker_spacing * 0 / 3, marker_spacing),
        markersize=marker_size,
    )
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[1],
        label="$\it{Pitch}$",
        color=pitch_color,
        marker=".",
        markevery=(marker_spacing * 1 / 3, marker_spacing),
        markersize=marker_size,
    )
    moment_coefficients_axes.plot(
        times,
        total_near_field_moment_coefficients_wind_axes[2],
        label="$\it{Yaw}$",
        color=yaw_color,
        marker=".",
        markevery=(marker_spacing * 2 / 3, marker_spacing),
        markersize=marker_size,
    )

    # Name the axis labels and the title.
    moment_coefficients_axes.set_xlabel(
        "$\it{Time\ (s)}$",
        color=figure_text_color,
    )
    moment_coefficients_axes.set_ylabel(
        "$\it{Coefficient\ (Unitless)}$",
        color=figure_text_color,
    )
    moment_coefficients_axes.set_title(
        "$\it{Moment\ Coefficients\ vs.\ Time}$",
        color=figure_text_color,
    )

    # Set the plot's background color.
    moment_coefficients_figure.patch.set_facecolor(figure_background_color)
    moment_coefficients_axes.set_facecolor(figure_background_color)

    # Add a legend.
    moment_coefficients_axes.legend(
        facecolor=figure_background_color,
        edgecolor=figure_background_color,
        labelcolor=figure_text_color,
    )

    # Show the plot.
    if not testing:
        moment_coefficients_figure.show()


# ToDo: Document this method.
def print_steady_results(steady_solver):
    """

    :param steady_solver:
    :return:
    """

    for airplane_num, airplane in enumerate(steady_solver.airplanes):
        print("Airplane ", airplane_num + 1, ":", sep="")
        # Print out the total forces and moments.
        print("\tForces in Wind Axes:")
        print(
            "\t\tInduced Drag:\t\t\t",
            np.round(airplane.total_near_field_force_wind_axes[0], 3),
            " N",
            sep="",
        )
        print(
            "\t\tSide Force:\t\t\t\t",
            np.round(airplane.total_near_field_force_wind_axes[1], 3),
            " N",
            sep="",
        )
        print(
            "\t\tLift:\t\t\t\t\t",
            np.round(airplane.total_near_field_force_wind_axes[2], 3),
            " N",
            sep="",
        )
        print("\n\tMoments in Wind Axes:")
        print(
            "\t\tRolling Moment:\t\t\t",
            np.round(airplane.total_near_field_moment_wind_axes[0], 3),
            " Nm",
            sep="",
        )
        print(
            "\t\tPitching Moment:\t\t",
            np.round(airplane.total_near_field_moment_wind_axes[1], 3),
            " Nm",
            sep="",
        )
        print(
            "\t\tYawing Moment:\t\t\t",
            np.round(airplane.total_near_field_moment_wind_axes[2], 3),
            " Nm",
            sep="",
        )

        # Print out the coefficients.
        print("\n\tCoefficients in Wind Axes:")
        print(
            "\t\tCDi:\t\t\t\t\t",
            np.round(airplane.total_near_field_force_coefficients_wind_axes[0], 3),
            sep="",
        )
        print(
            "\t\tCY:\t\t\t\t\t\t",
            np.round(airplane.total_near_field_force_coefficients_wind_axes[1], 3),
            sep="",
        )
        print(
            "\t\tCL:\t\t\t\t\t\t",
            np.round(airplane.total_near_field_force_coefficients_wind_axes[2], 3),
            sep="",
        )
        print(
            "\t\tCl:\t\t\t\t\t\t",
            np.round(airplane.total_near_field_moment_coefficients_wind_axes[0], 3),
            sep="",
        )
        print(
            "\t\tCm:\t\t\t\t\t\t",
            np.round(airplane.total_near_field_moment_coefficients_wind_axes[1], 3),
            sep="",
        )
        print(
            "\t\tCn:\t\t\t\t\t\t",
            np.round(airplane.total_near_field_moment_coefficients_wind_axes[2], 3),
            sep="",
        )

        # If there's more airplane's whose results are going to be printed, print new
        # line to separate them.
        if (airplane_num + 1) < steady_solver.num_airplanes:
            print("")


# ToDo: Document this method.
def get_panel_surfaces(
    airplanes,
):
    """

    :param airplanes:
    :return:
    """

    # Initialize empty arrays to hold the panel vertices and faces.
    panel_vertices = np.empty((0, 3), dtype=int)
    panel_faces = np.empty(0, dtype=int)

    # Initialize a variable to keep track of how many panels have been added thus far.
    panel_num = 0

    # Increment through the airplanes' wings.
    for airplane in airplanes:
        for wing in airplane.wings:

            # Unravel the wing's panel matrix and iterate through it.
            panels = np.ravel(wing.panels)
            for panel in panels:

                # Stack this panel's vertices and faces. Look through the PolyData
                # documentation for more details.
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
                        (panel_num * 4),
                        (panel_num * 4) + 1,
                        (panel_num * 4) + 2,
                        (panel_num * 4) + 3,
                    ]
                )

                # Stack this panel's vertices and faces with the array of all the
                # vertices and faces.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # Update the number of previous panels.
                panel_num += 1

    # Return the panel surfaces.
    return pv.PolyData(panel_vertices, panel_faces)


def get_wake_ring_vortex_surfaces(solver, step):
    """This function returns the PolyData object for the surface of wake ring
    vortices at a given time step.

    :param solver: UnsteadyRingVortexLatticeMethodSolver
        This is the unsteady solver with the wake.
    :param step: int
        This is the step number at which to look at the wake.
    :return: PolyData
        This is the PolyData object of the wake surface that can be displayed by
        other output methods.
    """
    num_wake_ring_vortices = solver.num_wake_ring_vortices_list[step]
    wake_ring_vortex_front_right_vertices = (
        solver.wake_ring_vortex_front_right_vertices_list[step]
    )
    wake_ring_vortex_front_left_vertices = (
        solver.wake_ring_vortex_front_left_vertices_list[step]
    )
    wake_ring_vortex_back_left_vertices = (
        solver.wake_ring_vortex_back_left_vertices_list[step]
    )
    wake_ring_vortex_back_right_vertices = (
        solver.wake_ring_vortex_back_right_vertices_list[step]
    )

    # Initialize empty arrays to hold each wake ring vortex's vertices and its face.
    wake_ring_vortex_vertices = np.zeros((0, 3), dtype=int)
    wake_ring_vortex_faces = np.zeros(0, dtype=int)

    for wake_ring_vortex_num in range(num_wake_ring_vortices):

        this_front_right_vertex = wake_ring_vortex_front_right_vertices[
            wake_ring_vortex_num
        ]
        this_front_left_vertex = wake_ring_vortex_front_left_vertices[
            wake_ring_vortex_num
        ]
        this_back_left_vertex = wake_ring_vortex_back_left_vertices[
            wake_ring_vortex_num
        ]
        this_back_right_vertex = wake_ring_vortex_back_right_vertices[
            wake_ring_vortex_num
        ]

        wake_ring_vortex_vertices_to_add = np.vstack(
            (
                this_front_left_vertex,
                this_front_right_vertex,
                this_back_right_vertex,
                this_back_left_vertex,
            )
        )
        wake_ring_vortex_face_to_add = np.array(
            [
                4,
                (wake_ring_vortex_num * 4),
                (wake_ring_vortex_num * 4) + 1,
                (wake_ring_vortex_num * 4) + 2,
                (wake_ring_vortex_num * 4) + 3,
            ]
        )

        # Stack this wake ring vortex's vertices and faces.
        wake_ring_vortex_vertices = np.vstack(
            (wake_ring_vortex_vertices, wake_ring_vortex_vertices_to_add)
        )
        wake_ring_vortex_faces = np.hstack(
            (wake_ring_vortex_faces, wake_ring_vortex_face_to_add)
        )

        # Increment the wake ring vortex counter.
        wake_ring_vortex_num += 1

    # Return the vortex surfaces.
    return pv.PolyData(wake_ring_vortex_vertices, wake_ring_vortex_faces)


def get_scalars(
    airplanes,
):
    """This function gets the delta pressure values from a problem's airplane
    objects, and puts them into a 1D array to be used as scalars for display by other
    output methods.

    :param airplanes: list of Airplane objects
        This is the list of airplane objects with the scalars we are collecting.
    :return scalars: 1D array of ints
        This is the 1D array of integers for each panel's delta pressure values.
    """

    # Initialize an empty array to hold the scalars.
    scalars = np.empty(0, dtype=int)

    # Increment through the airplanes' wings.
    for airplane in airplanes:
        for wing in airplane.wings:

            # Unravel the wing's panel matrix and iterate through it.
            panels = np.ravel(wing.panels)
            for panel in panels:

                # Stack this panel's scalars.
                scalar_to_add = panel.delta_pressure
                scalars = np.hstack((scalars, scalar_to_add))

    # Return the resulting 1D array of scalars.
    return scalars
