"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the geometry of the airplanes in a solver object.

    animate: Create an animation of a problem's movement.

    animate_frames: This function creates an animation of the solver's geometries but
    it saves each time step as a separate PNG, instead of compiling them as frames in
    a GIF.

    plot_results_versus_time: This method takes in an unsteady solver object,
    and plots the geometries' forces, moments, force coefficients, and moment
    coefficients as a function of time.

    print_steady_results: This function prints the forces, moments,
    force coefficients, and moment coefficients calculated by a steady solver.

    print_unsteady_results: This function prints the averages of the forces, moments,
    force coefficients, and moment coefficients calculated by a unsteady solver.

    get_panel_surfaces: This function returns a PolyData representation of the wing
    panel surfaces associated with all the airplanes in a given list.

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

# Define the color and colormaps used by the visualization functions.
sequential_color_map = "speed"
diverging_color_map = "delta"
wake_vortex_color = "white"
panel_color = "chartreuse"
streamline_color = "orchid"
plotter_background_color = "black"
figure_background_color = "black"
text_color = "white"

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

# Set constants for the color maps, scalar bars, and text boxes.
color_map_num_sig = 3
bar_title_font_size = 20
bar_label_font_size = 14
bar_width = 0.5
bar_position_x = 0.25
bar_position_y = 0.05
bar_n_labels = 2
text_max_position = (0.85, 0.075)
text_min_position = (0.85, 0.05)
text_font_size = 7

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
    """Draw the geometry of the airplanes in a solver object.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry and attributes are to be plotted.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels.
        See the src.panel.update_pressure() method for more details on the pressure
        value. The default value is false.
    :param show_streamlines: bool, optional
        Set this variable to true to show the streamlines emanating from the back of
        the wings. The default value is False.
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
            c_min = max(
                np.mean(scalars) - color_map_num_sig * np.std(scalars), np.min(scalars)
            )
            c_max = min(
                np.mean(scalars) + color_map_num_sig * np.std(scalars), np.max(scalars)
            )
        else:
            color_map = diverging_color_map
            c_min = -color_map_num_sig * np.std(scalars)
            c_max = color_map_num_sig * np.std(scalars)

        # Add the panel surfaces to the plotter with the pressure scalars.
        scalar_bar_args = dict(
            title="Normal Pressure (Pa)",
            title_font_size=bar_title_font_size,
            label_font_size=bar_label_font_size,
            width=bar_width,
            position_x=bar_position_x,
            position_y=bar_position_y,
            n_labels=bar_n_labels,
            fmt="%.1f",
        )
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            cmap=color_map,
            clim=[c_min, c_max],
            scalars=scalars,
            smooth_shading=True,
            scalar_bar_args=scalar_bar_args,
        )
        plotter.add_text(
            text="Max: " + str(round(max(scalars), 1)) + " Pa",
            position=text_max_position,
            font_size=text_font_size,
            viewport=True,
        )
        plotter.add_text(
            text="Min: " + str(round(min(scalars), 1)) + " Pa",
            position=text_min_position,
            font_size=text_font_size,
            viewport=True,
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
    """Create an animation of a solver's geometries.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry is to be animated.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels.
        See the src.panel.update_pressure() method for more details on the pressure
        value. The default value is false.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices.
        The default value is false.
    :param keep_file: bool, optional
        Set this variable to false in order to not save the resulting GIF. The
        default value is true.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

    # Get this solver's problems' airplanes. This will become a list of lists,
    # with the first index being the time step and the second index identifying each
    # of the solver's airplanes at that time step.
    step_airplanes = []
    for steady_problem in unsteady_solver.steady_problems:
        step_airplanes.append(steady_problem.airplanes)

    # Initialize the plotter and get the color map. Also, turn off the lighting to avoid
    # making the scalar bar jittery.
    plotter = pv.Plotter(lighting="none")

    # Initialize values to hold the color map choice and its limits.
    c_min = 0
    c_max = 0
    color_map = None

    # Initialize an empty array to hold all of the problem's scalars.
    all_scalars = np.empty(0, dtype=int)

    # Check if the user wants to show pressures.
    if show_delta_pressures:

        # Now iterate through each time step and gather all of the scalars for its
        # list of airplanes. These values will be used to configure the color map.
        for airplanes in step_airplanes:
            scalars_to_add = get_scalars(airplanes)
            all_scalars = np.hstack((all_scalars, scalars_to_add))

        # Choose the color map and set its limits based on if the min and max scalars
        # across all time steps have the same sign (sequential color map) or if they
        # have different signs (diverging color map).
        if np.sign(np.min(all_scalars)) == np.sign(np.max(all_scalars)):
            color_map = sequential_color_map
            c_min = max(
                np.mean(all_scalars) - color_map_num_sig * np.std(all_scalars),
                np.min(all_scalars),
            )
            c_max = min(
                np.mean(all_scalars) + color_map_num_sig * np.std(all_scalars),
                np.max(all_scalars),
            )
        else:
            color_map = diverging_color_map
            c_min = -color_map_num_sig * np.std(all_scalars)
            c_max = color_map_num_sig * np.std(all_scalars)

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surfaces = get_panel_surfaces(step_airplanes[0])

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the
    # plotter with the pressure scalars. Otherwise, add the panel surfaces without
    # the pressure scalars.
    if show_delta_pressures and first_results_step == 0:
        scalars = get_scalars(step_airplanes[0])

        # Add the panel surfaces to the plotter with the pressure scalars.
        scalar_bar_args = dict(
            title="Normal Pressure (Pa)",
            title_font_size=bar_title_font_size,
            label_font_size=bar_label_font_size,
            width=bar_width,
            position_x=bar_position_x,
            position_y=bar_position_y,
            n_labels=bar_n_labels,
            fmt="%.1f",
        )
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            cmap=color_map,
            clim=[c_min, c_max],
            scalars=scalars,
            smooth_shading=True,
            scalar_bar_args=scalar_bar_args,
        )
        plotter.add_text(
            text="Max: " + str(round(max(all_scalars), 1)) + " Pa",
            position=text_max_position,
            font_size=text_font_size,
            viewport=True,
        )
        plotter.add_text(
            text="Min: " + str(round(min(all_scalars), 1)) + " Pa",
            position=text_min_position,
            font_size=text_font_size,
            viewport=True,
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

    # Begin to iterate through all the other steps' airplanes.
    for airplanes in step_airplanes[1:]:

        # Clear the plotter.
        plotter.clear()

        # Get the panel surfaces.
        panel_surfaces = get_panel_surfaces(airplanes)

        # If the user wants to show the wake ring vortices, then get their surfaces
        # and plot them.
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
        # surfaces to the plotter with the pressure scalars. Otherwise, add the panel
        # surfaces without the pressure scalars.
        if show_delta_pressures and first_results_step <= current_step:
            scalars = get_scalars(airplanes)
            # Add the panel surfaces to the plotter with the pressure scalars.
            scalar_bar_args = dict(
                title="Normal Pressure (Pa)",
                title_font_size=bar_title_font_size,
                label_font_size=bar_label_font_size,
                width=bar_width,
                position_x=bar_position_x,
                position_y=bar_position_y,
                n_labels=bar_n_labels,
                fmt="%.1f",
            )
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                cmap=color_map,
                clim=[c_min, c_max],
                scalars=scalars,
                smooth_shading=True,
                scalar_bar_args=scalar_bar_args,
            )
            plotter.add_text(
                text="Max: " + str(round(max(all_scalars), 1)) + " Pa",
                position=text_max_position,
                font_size=text_font_size,
                viewport=True,
            )
            plotter.add_text(
                text="Min: " + str(round(min(all_scalars), 1)) + " Pa",
                position=text_min_position,
                font_size=text_font_size,
                viewport=True,
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


def animate_frames(
    unsteady_solver,
    show_delta_pressures=False,
    show_wake_vortices=False,
):
    """This function creates an animation of the solver's geometries but it saves
    each time step as a separate PNG, instead of compiling them as frames in a GIF.

    Note: By uploading the images to an online tool (such as
    https://ezgif.com/webp-maker/), the frames can be compiled into an animated,
    transparent WebP file. This technique was used to create the first simulation
    image shown in README.md.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry is to be animated.
    :param show_delta_pressures: bool, optional
        Set this variable to true to show the change in pressure across the panels.
        See the src.panel.update_pressure() method for more details on the pressure
        value. The default value is false.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices.
        The default value is false.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

    # Get this solver's problems' airplanes. This will become a list of lists,
    # with the first index being the time step and the second index identifying each
    # of the solver's airplanes at that time step.
    step_airplanes = []
    for steady_problem in unsteady_solver.steady_problems:
        step_airplanes.append(steady_problem.airplanes)

    # Initialize the plotter and get the color map. Also, turn off the lighting to avoid
    # making the scalar bar jittery.
    plotter = pv.Plotter(lighting="none")
    pv.rcParams["transparent_background"] = True

    # Initialize values to hold the color map choice and its limits.
    c_min = 0
    c_max = 0
    color_map = None

    # Initialize an empty array to hold all of the problem's scalars.
    all_scalars = np.empty(0, dtype=int)

    # Check if the user wants to show pressures.
    if show_delta_pressures:

        # Now iterate through each time step and gather all of the scalars for its
        # list of airplanes. These values will be used to configure the color map.
        for airplanes in step_airplanes:
            scalars_to_add = get_scalars(airplanes)
            all_scalars = np.hstack((all_scalars, scalars_to_add))

        # Choose the color map and set its limits based on if the min and max scalars
        # across all time steps have the same sign (sequential color map) or if they
        # have different signs (diverging color map).
        if np.sign(np.min(all_scalars)) == np.sign(np.max(all_scalars)):
            color_map = sequential_color_map
            c_min = max(
                np.mean(all_scalars) - color_map_num_sig * np.std(all_scalars),
                np.min(all_scalars),
            )
            c_max = min(
                np.mean(all_scalars) + color_map_num_sig * np.std(all_scalars),
                np.max(all_scalars),
            )
        else:
            color_map = diverging_color_map
            c_min = -color_map_num_sig * np.std(all_scalars)
            c_max = color_map_num_sig * np.std(all_scalars)

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surfaces = get_panel_surfaces(step_airplanes[0])

    # Check if the user wants to plot pressures. If so, add the panel surfaces to the
    # plotter with the pressure scalars. Otherwise, add the panel surfaces without
    # the pressure scalars.
    if show_delta_pressures and first_results_step == 0:
        scalars = get_scalars(step_airplanes[0])

        # Add the panel surfaces to the plotter with the pressure scalars.
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            cmap=color_map,
            clim=[c_min, c_max],
            scalars=scalars,
            smooth_shading=True,
            show_scalar_bar=False,
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

    # Print a message to the console on how to set up the window.
    print(
        'Orient the view, then press "q" to close the window and produce the '
        "animation."
    )

    # Set up the camera and close the window.
    plotter.show(cpos=(-1, -1, 1), full_screen=False, auto_close=False)

    plotter.screenshot(filename="0", transparent_background=True)

    # Initialize a variable to keep track of which step we are on.
    current_step = 1

    # Begin to iterate through all the other steps' airplanes.
    for airplanes in step_airplanes[1:]:

        # Clear the plotter.
        plotter.clear()

        # Get the panel surfaces.
        panel_surfaces = get_panel_surfaces(airplanes)

        # If the user wants to show the wake ring vortices, then get their surfaces
        # and plot them.
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
        # surfaces to the plotter with the pressure scalars. Otherwise, add the panel
        # surfaces without the pressure scalars.
        if show_delta_pressures and first_results_step <= current_step:
            scalars = get_scalars(airplanes)

            # Add the panel surfaces to the plotter with the pressure scalars.
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                cmap=color_map,
                clim=[c_min, c_max],
                scalars=scalars,
                smooth_shading=True,
                show_scalar_bar=False,
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

        plotter.screenshot(filename=str(current_step), transparent_background=True)
        plotter.clear()

        # Increment the step number tracker.
        current_step += 1

    # Close the animation and delete the plotter.
    plotter.close()


def plot_results_versus_time(unsteady_solver, testing=False):
    """This method takes in an unsteady solver object, and plots the geometries'
    forces, moments, force coefficients, and moment coefficients as a function of time.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose resulting forces, moments, and coefficients
        are to be plotted.
    :param testing: bool, Optional
        This boolean determines if the plots will be shown. If true, no plots will be
        shown. It is useful for testing, where the user wants to know that the plots
        were created without having to show them. It's default value is false.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

    # Get this solver's time step characteristics. Note that the first time step (
    # time step 0), occurs at 0 seconds.
    num_steps = unsteady_solver.num_steps
    delta_time = unsteady_solver.delta_time
    num_airplanes = unsteady_solver.num_airplanes
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
    total_near_field_force_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )
    total_near_field_force_coefficients_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )
    total_near_field_moment_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )
    total_near_field_moment_coefficients_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )

    # Initialize a variable to track position in the results arrays.
    results_step = 0

    # Iterate through the time steps and add the results to their respective matrices.
    for step in range(first_results_step, num_steps):

        # Get the airplanes from the problem at this step.
        airplanes = unsteady_solver.steady_problems[step].airplanes

        # Iterate through this step's airplanes.
        for airplane_id, airplane in enumerate(airplanes):
            total_near_field_force_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_force_wind_axes
            total_near_field_force_coefficients_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_force_coefficients_wind_axes
            total_near_field_moment_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_moment_wind_axes
            total_near_field_moment_coefficients_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_moment_coefficients_wind_axes

        results_step += 1

    # Iterate through the airplane id's to plot each airplane's figures.
    for airplane_id in range(num_airplanes):

        # Initialize the four figures.
        force_figure, force_axes = plt.subplots()
        force_coefficients_figure, force_coefficients_axes = plt.subplots()
        moment_coefficients_figure, moment_coefficients_axes = plt.subplots()
        moment_figure, moment_axes = plt.subplots()

        # Remove all the plots' top and right spines.
        force_axes.spines.right.set_visible(False)
        force_axes.spines.top.set_visible(False)
        force_coefficients_axes.spines.right.set_visible(False)
        force_coefficients_axes.spines.top.set_visible(False)
        moment_axes.spines.right.set_visible(False)
        moment_axes.spines.top.set_visible(False)
        moment_coefficients_axes.spines.right.set_visible(False)
        moment_coefficients_axes.spines.top.set_visible(False)

        # Format all the plots' spine and label colors.
        force_axes.spines.bottom.set_color(text_color)
        force_axes.spines.left.set_color(text_color)
        force_axes.xaxis.label.set_color(text_color)
        force_axes.yaxis.label.set_color(text_color)
        force_coefficients_axes.spines.bottom.set_color(text_color)
        force_coefficients_axes.spines.left.set_color(text_color)
        force_coefficients_axes.xaxis.label.set_color(text_color)
        force_coefficients_axes.yaxis.label.set_color(text_color)
        moment_coefficients_axes.spines.bottom.set_color(text_color)
        moment_coefficients_axes.spines.left.set_color(text_color)
        moment_coefficients_axes.xaxis.label.set_color(text_color)
        moment_coefficients_axes.yaxis.label.set_color(text_color)
        moment_axes.spines.bottom.set_color(text_color)
        moment_axes.spines.left.set_color(text_color)
        moment_axes.xaxis.label.set_color(text_color)
        moment_axes.yaxis.label.set_color(text_color)

        # Format all the plots' tick colors.
        force_axes.tick_params(axis="x", colors=text_color)
        force_axes.tick_params(axis="y", colors=text_color)
        force_coefficients_axes.tick_params(axis="x", colors=text_color)
        force_coefficients_axes.tick_params(axis="y", colors=text_color)
        moment_coefficients_axes.tick_params(axis="x", colors=text_color)
        moment_coefficients_axes.tick_params(axis="y", colors=text_color)
        moment_axes.tick_params(axis="x", colors=text_color)
        moment_axes.tick_params(axis="y", colors=text_color)

        # Format all the plots' background colors.
        force_figure.patch.set_facecolor(figure_background_color)
        force_axes.set_facecolor(figure_background_color)
        force_coefficients_figure.patch.set_facecolor(figure_background_color)
        force_coefficients_axes.set_facecolor(figure_background_color)
        moment_figure.patch.set_facecolor(figure_background_color)
        moment_axes.set_facecolor(figure_background_color)
        moment_coefficients_figure.patch.set_facecolor(figure_background_color)
        moment_coefficients_axes.set_facecolor(figure_background_color)

        # Populate the plots.
        force_axes.plot(
            times,
            total_near_field_force_wind_axes[airplane_id, 0],
            label="$\it{Induced\ Drag}$",
            color=drag_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_axes.plot(
            times,
            total_near_field_force_wind_axes[airplane_id, 1],
            label="$\it{Side\ Force}$",
            color=side_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_axes.plot(
            times,
            total_near_field_force_wind_axes[airplane_id, 2],
            label="$\it{Lift}$",
            color=lift_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_coefficients_axes.plot(
            times,
            total_near_field_force_coefficients_wind_axes[airplane_id, 0],
            label="$\it{Induced\ Drag}$",
            color=drag_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_coefficients_axes.plot(
            times,
            total_near_field_force_coefficients_wind_axes[airplane_id, 1],
            label="$\it{Side\ Force}$",
            color=side_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_coefficients_axes.plot(
            times,
            total_near_field_force_coefficients_wind_axes[airplane_id, 2],
            label="$\it{Lift}$",
            color=lift_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_axes.plot(
            times,
            total_near_field_moment_wind_axes[airplane_id, 0],
            label="$\it{Roll}$",
            color=roll_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_axes.plot(
            times,
            total_near_field_moment_wind_axes[airplane_id, 1],
            label="$\it{Pitch}$",
            color=pitch_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_axes.plot(
            times,
            total_near_field_moment_wind_axes[airplane_id, 2],
            label="$\it{Yaw}$",
            color=yaw_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            total_near_field_moment_coefficients_wind_axes[airplane_id, 0],
            label="$\it{Roll}$",
            color=roll_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            total_near_field_moment_coefficients_wind_axes[airplane_id, 1],
            label="$\it{Pitch}$",
            color=pitch_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            total_near_field_moment_coefficients_wind_axes[airplane_id, 2],
            label="$\it{Yaw}$",
            color=yaw_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )

        # Find and format this airplane's name for use in the plot titles.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[airplane_id].name
        airplane_name_split = airplane_name.split(" ")
        title_prefix = "$\it{"
        for sub in airplane_name_split:
            title_prefix += sub + "\ "
        force_title = title_prefix + "Forces\ vs.\ Time}$"
        force_coefficient_title = title_prefix + "Force\ Coefficients\ vs.\ Time}$"
        moment_title = title_prefix + "Moments\ vs.\ Time}$"
        moment_coefficient_title = title_prefix + "Moment\ Coefficients\ vs.\ Time}$"

        # Name the plots' axis labels and titles.
        force_axes.set_xlabel("$\it{Time\ (s)}$", color=text_color)
        force_axes.set_ylabel("$\it{Force\ (N)}$", color=text_color)
        force_axes.set_title(force_title, color=text_color)
        force_coefficients_axes.set_xlabel("$\it{Time\ (s)}$", color=text_color)
        force_coefficients_axes.set_ylabel("$\it{Coefficient$", color=text_color)
        force_coefficients_axes.set_title(force_coefficient_title, color=text_color)
        moment_axes.set_xlabel("$\it{Time\ (s)}$", color=text_color)
        moment_axes.set_ylabel("$\it{Moment\ (N\ m)}$", color=text_color)
        moment_axes.set_title(moment_title, color=text_color)
        moment_coefficients_axes.set_xlabel("$\it{Time\ (s)}$", color=text_color)
        moment_coefficients_axes.set_ylabel("$\it{Coefficient$", color=text_color)
        moment_coefficients_axes.set_title(moment_coefficient_title, color=text_color)

        # Format the plots' legends.
        force_axes.legend(
            facecolor=figure_background_color,
            edgecolor=figure_background_color,
            labelcolor=text_color,
        )
        force_coefficients_axes.legend(
            facecolor=figure_background_color,
            edgecolor=figure_background_color,
            labelcolor=text_color,
        )
        moment_axes.legend(
            facecolor=figure_background_color,
            edgecolor=figure_background_color,
            labelcolor=text_color,
        )
        moment_coefficients_axes.legend(
            facecolor=figure_background_color,
            edgecolor=figure_background_color,
            labelcolor=text_color,
        )

        # Show this airplane's plots, if not testing.
        if not testing:
            force_figure.show()
            force_coefficients_figure.show()
            moment_figure.show()
            moment_coefficients_figure.show()


def print_steady_results(steady_solver):
    """This function prints the forces, moments, force coefficients, and moment
    coefficients calculated by a steady solver.

    :param steady_solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver
        This is the solver object with the results to be printed.
    :return: None
    """

    for airplane_num, airplane in enumerate(steady_solver.airplanes):
        print("Airplane ", airplane.name, ":", sep="")

        # Print out this airplane's forces.
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

        # Print out this airplane's moments.
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

        # Print out this airplane's force coefficients.
        print("\n\tForce Coefficients in Wind Axes:")
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

        # Print out this airplane's moment coefficients.
        print("\n\tMoment Coefficients in Wind Axes:")
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

        # If the results from more airplanes are going to be printed, print new line
        # to separate them.
        if (airplane_num + 1) < steady_solver.num_airplanes:
            print("")


def print_unsteady_results(unsteady_solver):
    """This function prints the averages of the forces, moments, force coefficients,
    and moment coefficients calculated by a unsteady solver.

    Note: This method averages the values for every time step that calculated
    results. Therefore, the averages are not necessarily the final-cycle averages.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver or
        This is the solver object with the results to be printed.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

    # Get this solver's time step characteristics. Note that the first time step (
    # time step 0), occurs at 0 seconds.
    num_steps = unsteady_solver.num_steps
    num_airplanes = unsteady_solver.num_airplanes
    num_steps_with_results = num_steps - first_results_step

    # Initialize matrices to hold the forces, moments, and coefficients at each of
    # the time steps that has results.
    total_near_field_force_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )
    total_near_field_force_coefficients_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )
    total_near_field_moment_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )
    total_near_field_moment_coefficients_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_with_results)
    )

    # Initialize a variable to track position in the results arrays.
    results_step = 0

    # Iterate through the time steps with results and add the results to their
    # respective matrices.
    for step in range(first_results_step, num_steps):

        # Get the airplanes from the problem at this step.
        airplanes = unsteady_solver.steady_problems[step].airplanes

        # Iterate through this step's airplanes.
        for airplane_id, airplane in enumerate(airplanes):
            total_near_field_force_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_force_wind_axes
            total_near_field_force_coefficients_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_force_coefficients_wind_axes
            total_near_field_moment_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_moment_wind_axes
            total_near_field_moment_coefficients_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_moment_coefficients_wind_axes

        results_step += 1

    # For each airplane object, calculate and print the average force, moment,
    # force coefficient, and moment coefficient values.
    for airplane_id, airplane in enumerate(
        unsteady_solver.steady_problems[0].airplanes
    ):
        # Calculate the average values.
        this_induced_drag = np.mean(total_near_field_force_wind_axes[airplane_id, 0, :])
        this_side_force = np.mean(total_near_field_force_wind_axes[airplane_id, 1, :])
        this_lift = np.mean(total_near_field_force_wind_axes[airplane_id, 2, :])
        this_rolling_moment = np.mean(
            total_near_field_moment_wind_axes[airplane_id, 0, :]
        )
        this_pitching_moment = np.mean(
            total_near_field_moment_wind_axes[airplane_id, 1, :]
        )
        this_yawing_moment = np.mean(
            total_near_field_moment_wind_axes[airplane_id, 2, :]
        )
        this_induced_drag_coefficient = np.mean(
            total_near_field_force_coefficients_wind_axes[airplane_id, 0, :]
        )
        this_side_force_coefficient = np.mean(
            total_near_field_force_coefficients_wind_axes[airplane_id, 1, :]
        )
        this_lift_coefficient = np.mean(
            total_near_field_force_coefficients_wind_axes[airplane_id, 2, :]
        )
        this_rolling_coefficient = np.mean(
            total_near_field_moment_coefficients_wind_axes[airplane_id, 0, :]
        )
        this_pitching_coefficient = np.mean(
            total_near_field_moment_coefficients_wind_axes[airplane_id, 1, :]
        )
        this_yawing_coefficient = np.mean(
            total_near_field_moment_coefficients_wind_axes[airplane_id, 2, :]
        )

        print(airplane.name, ":", sep="")

        # Print out this airplane's average forces.
        print("\tAverage Forces in Wind Axes:")
        print(
            "\t\tInduced Drag:\t\t\t",
            np.round(this_induced_drag, 3),
            " N",
            sep="",
        )
        print(
            "\t\tSide Force:\t\t\t\t",
            np.round(this_side_force, 3),
            " N",
            sep="",
        )
        print(
            "\t\tLift:\t\t\t\t\t",
            np.round(this_lift, 3),
            " N",
            sep="",
        )

        # Print out this airplane's average moments.
        print("\n\tAverage Moments in Wind Axes:")
        print(
            "\t\tRolling Moment:\t\t\t",
            np.round(this_rolling_moment, 3),
            " Nm",
            sep="",
        )
        print(
            "\t\tPitching Moment:\t\t",
            np.round(this_pitching_moment, 3),
            " Nm",
            sep="",
        )
        print(
            "\t\tYawing Moment:\t\t\t",
            np.round(this_yawing_moment, 3),
            " Nm",
            sep="",
        )

        # Print out this airplane's average force coefficients.
        print("\n\tAverage Force Coefficients in Wind Axes:")
        print(
            "\t\tCDi:\t\t\t\t\t",
            np.round(this_induced_drag_coefficient, 3),
            sep="",
        )
        print(
            "\t\tCY:\t\t\t\t\t\t",
            np.round(this_side_force_coefficient, 3),
            sep="",
        )
        print(
            "\t\tCL:\t\t\t\t\t\t",
            np.round(this_lift_coefficient, 3),
            sep="",
        )

        # Print out this airplane's average moment coefficients.
        print("\n\tAverage Moment Coefficients in Wind Axes:")
        print(
            "\t\tCl:\t\t\t\t\t\t",
            np.round(this_rolling_coefficient, 3),
            sep="",
        )
        print(
            "\t\tCm:\t\t\t\t\t\t",
            np.round(this_pitching_coefficient, 3),
            sep="",
        )
        print(
            "\t\tCn:\t\t\t\t\t\t",
            np.round(this_yawing_coefficient, 3),
            sep="",
        )

        # If the results from more airplanes are going to be printed, print new line
        # to separate them.
        if (airplane_id + 1) < num_airplanes:
            print("")


def get_panel_surfaces(
    airplanes,
):
    """This function returns a PolyData representation of the wing panel surfaces
    associated with all the airplanes in a given list.

    :param airplanes: list of Airplane objects
        This is a list of airplane objects whose wing panel surfaces will be returned.
    :return: pv.PolyData
        This is a PolyData representation of the airplanes' wing panel surfaces.
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

                # Stack this panel's scalars. The scalar is the pressure on the
                # panel. See the src.panel.update_pressure() method for more details
                # on the pressure value.
                scalars = np.hstack((scalars, panel.delta_pressure))

    # Return the resulting 1D array of scalars.
    return scalars
