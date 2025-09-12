"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    draw: Draw the geometry of the airplanes in a solver object.

    animate: Create an animation of a problem's movement.

    plot_results_versus_time: This method takes in an unsteady solver object,
    and plots the geometries' forces, moments, force coefficients, and moment
    coefficients as a function of time.

    print_steady_results: This function prints the forces, moments,
    force coefficients, and moment coefficients calculated by a steady solver.

    print_unsteady_results: This function prints the final cycle-averaged of the
    forces, moments, force coefficients, and moment coefficients calculated by an
    unsteady solver.

    get_panel_surfaces: This function returns a PolyData representation of the wing
    panel surfaces associated with all the airplanes in a given list.

    get_wake_ring_vortex_surfaces: This function returns the PolyData object for the
    surface of wake ring vortices at a given time step.

    get_scalars: This function gets the coefficient values from a problem's airplane
    objects, and puts them into a 1D array to be used as scalars for display by other
    output methods.

    plot_scalars: This function plots a scalar bar, the mesh panels with a particular
    set of scalars, and labels for the minimum and maximum scalar values."""

import math
import time

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import webp

from . import unsteady_ring_vortex_lattice_method

# Define the color and colormaps used by the visualization functions.
sequential_color_map = "speed"
diverging_color_map = "delta"
wake_vortex_color = "white"
panel_color = "chartreuse"
streamline_color = "orchid"
plotter_background_color = "black"
figure_background_color = "None"
text_color = "#818181"
quality = 75
window_size = [1024, 768]

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
bar_title_font_size = 30
bar_label_font_size = 21
bar_width = 0.5
bar_position_x = 0.25
bar_position_y = 0.05
bar_n_labels = 2
text_max_position = (0.85, 0.075)
text_min_position = (0.85, 0.050)
text_speed_position = (0.05, 0.075)
text_font_size = 11

# Set the number of markers and the marker size for the results plots.
num_markers = 6
marker_size = 8

# Calculate the normalized spacing between the markers for the results plots.
marker_spacing = 1.0 / num_markers


def draw(
    solver,
    scalar_type=None,
    show_streamlines=False,
    show_wake_vortices=False,
    save=False,
    testing=False,
):
    """Draw the geometry of the airplanes in a solver object.

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry and attributes are to be plotted.
    :param scalar_type: str or None, optional
        This variable determines which values will be used to color the wing panels,
        if any. The default value is None, which gives colors all wing panels
        uniformly. Acceptable alternatives are "induced drag", "side force",
        and "lift", which respectively use each panel's induced drag, side force,
        or lift coefficient.
    :param show_streamlines: bool, optional
        Set this variable to true to show the streamlines emanating from the back of
        the wings. The default value is False.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices.
        The default value is False.
    :param save: bool, optional
        Set this variable to True to save the image as a WebP. The default value is
        False.
    :param testing: bool, optional
        Set this variable to True to close the image after 1 second, which is useful
        for running test suites. The default value is False.
    :return: None
    """

    # Initialize the plotter and set it to use parallel projection (instead of
    # perspective).
    plotter = pv.Plotter(window_size=window_size, lighting=None)
    plotter.enable_parallel_projection()

    # Get the solver's geometry.
    if isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        draw_step = solver.num_steps - 1
        airplanes = solver.steady_problems[draw_step].airplanes

        # If the user wants to show the wake ring vortices, then get their surfaces and
        # plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = get_wake_ring_vortex_surfaces(solver, draw_step)
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=False,
                color=wake_vortex_color,
            )
    else:
        airplanes = solver.airplanes

    # Check if the user wants to show scalars on the wing panels.
    show_scalars = False
    if (
        scalar_type == "induced drag"
        or scalar_type == "side force"
        or scalar_type == "lift"
    ):
        show_scalars = True

    # Get the panel surfaces.
    panel_surfaces = get_panel_surfaces(airplanes)

    # Check if the user wants to plot any scalars.
    if show_scalars:

        # Get the scalars
        these_scalars = get_scalars(airplanes, scalar_type)
        min_scalar = round(min(these_scalars), 2)
        max_scalar = round(max(these_scalars), 2)

        # Choose the color map and set its limits based on if the min and max scalars
        # have the same sign (sequential color map) or if they have different signs
        # (diverging color map).
        if np.sign(np.min(these_scalars)) == np.sign(np.max(these_scalars)):
            color_map = sequential_color_map
            c_min = max(
                np.mean(these_scalars) - color_map_num_sig * np.std(these_scalars),
                np.min(these_scalars),
            )
            c_max = min(
                np.mean(these_scalars) + color_map_num_sig * np.std(these_scalars),
                np.max(these_scalars),
            )
        else:
            color_map = diverging_color_map
            c_min = -color_map_num_sig * np.std(these_scalars)
            c_max = color_map_num_sig * np.std(these_scalars)

        plot_scalars(
            plotter,
            these_scalars,
            scalar_type,
            min_scalar,
            max_scalar,
            color_map,
            c_min,
            c_max,
            panel_surfaces,
        )
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
            smooth_shading=False,
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
                        smooth_shading=False,
                    )

    # Set the plotter's background color.
    plotter.set_background(color=plotter_background_color)
    if not testing:
        # Show the plotter so the user can adjust the camera position and window.
        # When the user closes the window, the plotter object itself won't close so
        # that it can be saved as an image if the user wants.
        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )
    else:
        # Show the plotter for 1 second, then proceed automatically. This is useful
        # for testing.
        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # If the user wants to save the image, take a screenshot, convert it into an
    # image object, and save it as a WebP.
    if save:
        image = webp.Image.fromarray(
            plotter.screenshot(
                filename=None,
                transparent_background=True,
                return_img=True,
            )
        )
        webp.save_image(
            img=image, file_path="Draw.webp", lossless=False, quality=quality
        )

    # Close all the plotters.
    pv.close_all()


def animate(
    unsteady_solver,
    scalar_type=None,
    show_wake_vortices=False,
    save=False,
    testing=False,
):
    """Create an animation of a solver's geometries.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose geometry is to be animated.
    :param scalar_type: str or None, optional
        This variable determines which values will be used to color the wing panels,
        if any. The default value is None, which gives colors all wing panels
        uniformly. Acceptable alternatives are "induced drag", "side force",
        and "lift", which respectively use each panel's induced drag, side force,
        or lift coefficient.
    :param show_wake_vortices: bool, optional
        Set this variable to true to show the airplane object's wake ring vortices.
        The default value is false.
    :param save: bool, optional
        Set this variable to True in order to save the resulting WebP animation. The
        default value is False.
    :param testing: bool, optional
        Set this variable to True to close the image after 1 second, which is useful
        for running test suites. The default value is False.
    :return: None
    """

    first_results_step = unsteady_solver.first_results_step

    # Get this solver's problems' airplanes. This will become a list of lists,
    # with the first index being the time step and the second index identifying each
    # of the solver's airplanes at that time step.
    step_airplanes = []
    for steady_problem in unsteady_solver.steady_problems:
        step_airplanes.append(steady_problem.airplanes)

    # Scale down the true-speed frames per second to at most 50 fps. This is the
    # maximum speed some programs can render WebPs.
    requested_fps = 1 / unsteady_solver.delta_time
    speed = 1
    if requested_fps > 50:
        speed = 50 / requested_fps
    actual_fps = math.floor(requested_fps * speed)

    # Initialize the plotter and get the color map.
    plotter = pv.Plotter(window_size=window_size, lighting=None)

    # Initialize values to hold the color map choice and its limits.
    c_min = 0
    c_max = 0
    color_map = None

    if save:
        # Add text to the animation that display's its speed relative to the true
        # speed.
        plotter.add_text(
            text="Speed: " + str(round(100 * speed)) + "%",
            position=text_speed_position,
            font_size=text_font_size,
            viewport=True,
            color=text_color,
        )

    # Check if the user wants to show scalars on the wing panels.
    show_scalars = False
    if (
        scalar_type == "induced drag"
        or scalar_type == "side force"
        or scalar_type == "lift"
    ):
        show_scalars = True

    # Initialize variables to hold the problems' scalars and their attributes.
    all_scalars = np.empty(0, dtype=float)
    min_scalar = None
    max_scalar = None

    # Check if the user wants to show scalars on the wing panels.
    if show_scalars:

        # Now iterate through each time step and gather all the scalars for its list
        # of airplanes. These values will be used to configure the color map.
        for airplanes in step_airplanes:
            scalars_to_add = get_scalars(airplanes, scalar_type)
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

        min_scalar = round(min(all_scalars), 2)
        max_scalar = round(max(all_scalars), 2)

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surfaces = get_panel_surfaces(step_airplanes[0])

    # Check if the user wants to show any scalars. If so, add the panel surfaces to
    # the plotter with these scalars.
    if show_scalars and first_results_step == 0:
        these_scalars = get_scalars(step_airplanes[0], scalar_type)

        plot_scalars(
            plotter,
            these_scalars,
            scalar_type,
            min_scalar,
            max_scalar,
            color_map,
            c_min,
            c_max,
            panel_surfaces,
        )
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
            smooth_shading=False,
        )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color=plotter_background_color)

    if not testing:
        # If not testing, print a message to the console on how to set up the window.
        print(
            'Orient the view, then press "q" to close the window and produce the animation.'
        )

        # Show the plotter so the user can set up the camera. Then, they will close the
        # window, starting the animation, but the plotter object will stay open.
        plotter.show(
            title="Rendering speed not to scale.",
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )
    else:
        # If we are testing, show the plotter for 1 second, then start the animation.
        plotter.show(
            title="Rendering speed not to scale.",
            cpos=(-1, -1, 1),
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # Start a list which will hold a WebP image of each frame.
    images = [
        webp.Image.fromarray(
            plotter.screenshot(
                transparent_background=True,
                return_img=True,
            )
        )
    ]

    # Initialize a variable to keep track of which step we are on.
    current_step = 1

    # Begin to iterate through all the other steps' airplanes.
    for airplanes in step_airplanes[1:]:

        # Clear the plotter.
        plotter.clear()

        # Get the panel surfaces.
        panel_surfaces = get_panel_surfaces(airplanes)

        if save:
            plotter.add_text(
                text="Speed: " + str(round(100 * speed)) + "%",
                position=text_speed_position,
                font_size=text_font_size,
                viewport=True,
                color=text_color,
            )

        # If the user wants to show the wake ring vortices, then get their surfaces
        # and plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = get_wake_ring_vortex_surfaces(
                unsteady_solver, current_step
            )
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=False,
                color=wake_vortex_color,
            )

        # Check if the user wants to plot scalars and this step is equal to or
        # greater than the first step with calculated results. If so, add the panel
        # surfaces to the plotter with the scalars.
        if show_scalars and first_results_step <= current_step:

            these_scalars = get_scalars(airplanes, scalar_type)

            plot_scalars(
                plotter,
                these_scalars,
                scalar_type,
                min_scalar,
                max_scalar,
                color_map,
                c_min,
                c_max,
                panel_surfaces,
            )
        else:
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                color=panel_color,
                smooth_shading=False,
            )

        # Append a WebP image of this frame to the list of frame images if the user
        # wants to save an animation.
        if save:
            images.append(
                webp.Image.fromarray(
                    plotter.screenshot(
                        filename=None,
                        transparent_background=True,
                        return_img=True,
                    )
                )
            )

        # Increment the step number tracker.
        current_step += 1

    # If the user wants to save the file, save the list of images as an animated WebP.
    if save:
        # Convert the list of WebP images to a WebP animation.
        webp.save_images(
            images, "Animate.webp", fps=actual_fps, lossless=False, quality=quality
        )

    # Close all the plotters.
    pv.close_all()


def plot_results_versus_time(unsteady_solver, show=True, save=False):
    """This method takes in an unsteady solver object, and plots the geometries'
    forces, moments, force coefficients, and moment coefficients as a function of time.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver
        This is the solver object whose resulting forces, moments, and coefficients
        are to be plotted.
    :param show: bool, Optional
        This boolean determines if the plots will be shown. If False, no plots will be
        shown, which is useful for testing when the user wants to know that the plots
        were created without having to show them. Its default value is True.
    :param save: bool, Optional
        This boolean determines if the plots will be saved as WebP images. The
        default value is False.
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
    num_steps_to_average = num_steps - first_results_step

    # Create a 1D array with the time at each time step where results have been
    # calculated.
    times = np.linspace(
        first_results_time_step_time,
        final_time_step_time,
        num_steps_to_average,
        endpoint=True,
    )

    # Initialize matrices to hold the forces, moments, and coefficients at each time
    # step which has results.
    total_near_field_force_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_to_average)
    )
    total_near_field_force_coefficients_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_to_average)
    )
    total_near_field_moment_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_to_average)
    )
    total_near_field_moment_coefficients_wind_axes = np.zeros(
        (num_airplanes, 3, num_steps_to_average)
    )

    # Initialize a variable to track position in the results arrays.
    results_step = 0

    # Iterate through the time steps and add the results to their respective matrices.
    for step in range(first_results_step, num_steps):

        # Get the airplanes from the problem at this step.
        airplanes = unsteady_solver.steady_problems[step].airplanes

        # Iterate through this step's airplanes.
        for airplane_id, airplane in enumerate(airplanes):
            total_near_field_force_wind_axes[airplane_id, :, results_step] = (
                airplane.total_near_field_force_wind_axes
            )
            total_near_field_force_coefficients_wind_axes[
                airplane_id, :, results_step
            ] = airplane.total_near_field_force_coefficients_wind_axes
            total_near_field_moment_wind_axes[airplane_id, :, results_step] = (
                airplane.total_near_field_moment_wind_axes
            )
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
            label="Induced Drag",
            color=drag_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_axes.plot(
            times,
            total_near_field_force_wind_axes[airplane_id, 1],
            label="Side Force",
            color=side_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_axes.plot(
            times,
            total_near_field_force_wind_axes[airplane_id, 2],
            label="Lift",
            color=lift_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_coefficients_axes.plot(
            times,
            total_near_field_force_coefficients_wind_axes[airplane_id, 0],
            label="Induced Drag",
            color=drag_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_coefficients_axes.plot(
            times,
            total_near_field_force_coefficients_wind_axes[airplane_id, 1],
            label="Side Force",
            color=side_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        force_coefficients_axes.plot(
            times,
            total_near_field_force_coefficients_wind_axes[airplane_id, 2],
            label="Lift",
            color=lift_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_axes.plot(
            times,
            total_near_field_moment_wind_axes[airplane_id, 0],
            label="Roll",
            color=roll_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_axes.plot(
            times,
            total_near_field_moment_wind_axes[airplane_id, 1],
            label="Pitch",
            color=pitch_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_axes.plot(
            times,
            total_near_field_moment_wind_axes[airplane_id, 2],
            label="Yaw",
            color=yaw_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            total_near_field_moment_coefficients_wind_axes[airplane_id, 0],
            label="Roll",
            color=roll_color,
            marker=".",
            markevery=(marker_spacing * 0 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            total_near_field_moment_coefficients_wind_axes[airplane_id, 1],
            label="Pitch",
            color=pitch_color,
            marker=".",
            markevery=(marker_spacing * 1 / 3, marker_spacing),
            markersize=marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            total_near_field_moment_coefficients_wind_axes[airplane_id, 2],
            label="Yaw",
            color=yaw_color,
            marker=".",
            markevery=(marker_spacing * 2 / 3, marker_spacing),
            markersize=marker_size,
        )

        # Find and format this airplane's name for use in the plot titles.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[airplane_id].name
        force_title = airplane_name + " Forces vs. Time"
        force_coefficient_title = airplane_name + " Force Coefficients vs. Time"
        moment_title = airplane_name + " Moments vs. Time"
        moment_coefficient_title = airplane_name + " Moment Coefficients vs. Time"

        # Name the plots' axis labels and titles.
        force_axes.set_xlabel("Time (s)", color=text_color)
        force_axes.set_ylabel("Force (N)", color=text_color)
        force_axes.set_title(force_title, color=text_color)
        force_coefficients_axes.set_xlabel("Time (s)", color=text_color)
        force_coefficients_axes.set_ylabel("Coefficient", color=text_color)
        force_coefficients_axes.set_title(force_coefficient_title, color=text_color)
        moment_axes.set_xlabel("Time (s)", color=text_color)
        moment_axes.set_ylabel("Moment (N m)", color=text_color)
        moment_axes.set_title(moment_title, color=text_color)
        moment_coefficients_axes.set_xlabel("Time (s)", color=text_color)
        moment_coefficients_axes.set_ylabel("Coefficient", color=text_color)
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

        # Save the figures as PNGs if the user wants to do so.
        if save:
            force_figure.savefig(
                airplane_name + " Forces.png",
                dpi=300,
            )
            force_coefficients_figure.savefig(
                airplane_name + " Force Coefficients.png",
                dpi=300,
            )
            moment_figure.savefig(
                airplane_name + " Moments.png",
                dpi=300,
            )
            moment_coefficients_figure.savefig(
                airplane_name + " Moment Coefficients.png",
                dpi=300,
            )

        # If the user wants to show the plots, do so.
        if show:
            force_figure.show()
            force_coefficients_figure.show()
            moment_figure.show()
            moment_coefficients_figure.show()
        else:
            plt.close("all")


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
    """This function prints the final cycle-averaged of the forces, moments,
    force coefficients, and moment coefficients calculated by an unsteady solver.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver or
        This is the solver object with the results to be printed.
    :return: None
    """
    forces = unsteady_solver.unsteady_problem.final_mean_near_field_forces_wind_axes
    moments = unsteady_solver.unsteady_problem.final_mean_near_field_moments_wind_axes
    force_coefficients = (
        unsteady_solver.unsteady_problem.final_mean_near_field_force_coefficients_wind_axes
    )
    moment_coefficients = (
        unsteady_solver.unsteady_problem.final_mean_near_field_moment_coefficients_wind_axes
    )

    # For each airplane object, calculate and print the average force, moment,
    # force coefficient, and moment coefficient values.
    for airplane_id, airplane in enumerate(
        unsteady_solver.steady_problems[0].airplanes
    ):
        these_forces = forces[airplane_id]
        these_moments = moments[airplane_id]
        these_force_coefficients = force_coefficients[airplane_id]
        these_moment_coefficients = moment_coefficients[airplane_id]

        this_induced_drag = these_forces[0]
        this_side_force = these_forces[1]
        this_lift = these_forces[2]
        this_rolling_moment = these_moments[0]
        this_pitching_moment = these_moments[1]
        this_yawing_moment = these_moments[2]

        this_induced_drag_coefficient = these_force_coefficients[0]
        this_side_force_coefficient = these_force_coefficients[1]
        this_lift_coefficient = these_force_coefficients[2]
        this_rolling_moment_coefficient = these_moment_coefficients[0]
        this_pitching_moment_coefficient = these_moment_coefficients[1]
        this_yawing_moment_coefficient = these_moment_coefficients[2]

        print(airplane.name, ":", sep="")

        # Print out this airplane's average forces.
        print("\tFinal Forces in Wind Axes:")
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
        print("\n\tFinal Moments in Wind Axes:")
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
        print("\n\tFinal Force Coefficients in Wind Axes:")
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
        print("\n\tFinal Moment Coefficients in Wind Axes:")
        print(
            "\t\tCl:\t\t\t\t\t\t",
            np.round(this_rolling_moment_coefficient, 3),
            sep="",
        )
        print(
            "\t\tCm:\t\t\t\t\t\t",
            np.round(this_pitching_moment_coefficient, 3),
            sep="",
        )
        print(
            "\t\tCn:\t\t\t\t\t\t",
            np.round(this_yawing_moment_coefficient, 3),
            sep="",
        )

        # If the results from more airplanes are going to be printed, print new line
        # to separate them.
        if (airplane_id + 1) < unsteady_solver.num_airplanes:
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
    panel_vertices = np.empty((0, 3), dtype=float)
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
                        panel.Flpp_G_Cg,
                        panel.Frpp_G_Cg,
                        panel.Brpp_G_Cg,
                        panel.Blpp_G_Cg,
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
    wake_ring_vortex_vertices = np.zeros((0, 3), dtype=float)
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
    scalar_type,
):
    """This function gets the coefficient values from a problem's airplane objects,
    and puts them into a 1D array to be used as scalars for display by other output
    methods.

    :param airplanes: list of Airplane objects
        This is the list of airplane objects with the scalars we are collecting.
    :param scalar_type: str
        This variable determines which scalar values will be returned. Acceptable
        inputs are "induced drag", "side force", and "lift", which respectively
        return each panel's induced drag, side force, or lift coefficient.
    :return scalars: 1D array of floats
        This is the 1D array of floats for each panel's coefficient value.
    """

    # Initialize an empty array to hold the scalars.
    scalars = np.empty(0, dtype=float)

    # Increment through the airplanes' wings.
    for airplane in airplanes:
        for wing in airplane.wings:

            # Unravel the wing's panel matrix and iterate through it.
            panels = np.ravel(wing.panels)
            for panel in panels:

                # Stack this panel's scalars.
                if scalar_type == "induced drag":
                    scalars = np.hstack((scalars, panel.induced_drag_coefficient))
                if scalar_type == "side force":
                    scalars = np.hstack((scalars, panel.side_force_coefficient))
                if scalar_type == "lift":
                    scalars = np.hstack((scalars, panel.lift_coefficient))

    # Return the resulting 1D array of scalars.
    return scalars


def plot_scalars(
    plotter,
    these_scalars,
    scalar_type,
    min_scalar,
    max_scalar,
    color_map,
    c_min,
    c_max,
    panel_surfaces,
):
    """This function plots a scalar bar, the mesh panels with a particular set of
    scalars, and labels for the minimum and maximum scalar values.

    :param plotter: `pyvista.Plotter`
        The plotter object used for visualization.
    :param these_scalars: 1D array of floats
        This is the 1D array of floats for each panel's coefficient value.
    :param scalar_type: str
        This variable determines which scalar values will be returned. Acceptable
        inputs are "induced drag", "side force", and "lift", which respectively
        return each panel's induced drag, side force, or lift coefficient.
    :param min_scalar: float
        Minimum scalar value, which is displayed as text on the plotter.
    :param max_scalar: float
        Maximum scalar value, which is displayed as text on the plotter.
    :param color_map: str
        Name of the colormap to use for scalar visualization.
    :param c_min: float
        Lower bound for the colormap scaling.
    :param c_max: float
        Upper bound for the colormap scaling.
    :param panel_surfaces: `pyvista.PolyData`
        PolyData object representing the mesh panel surfaces.
    :return: None
    """
    scalar_bar_args = dict(
        title=scalar_type.title() + " Coefficient",
        title_font_size=bar_title_font_size,
        label_font_size=bar_label_font_size,
        width=bar_width,
        position_x=bar_position_x,
        position_y=bar_position_y,
        n_labels=bar_n_labels,
        fmt="%.2f",
        color=text_color,
    )
    plotter.add_mesh(
        panel_surfaces,
        show_edges=True,
        cmap=color_map,
        clim=[c_min, c_max],
        scalars=these_scalars,
        smooth_shading=False,
        scalar_bar_args=scalar_bar_args,
    )
    plotter.add_text(
        text="Max: " + str(max_scalar),
        position=text_max_position,
        font_size=text_font_size,
        viewport=True,
        color=text_color,
    )
    plotter.add_text(
        text="Min: " + str(min_scalar),
        position=text_min_position,
        font_size=text_font_size,
        viewport=True,
        color=text_color,
    )
