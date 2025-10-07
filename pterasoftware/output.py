# NOTE: I've started refactoring this module.
"""This module contains useful functions for visualizing solutions to problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

# TODO: Update the function descriptions.
This module contains the following functions:
    draw: Draw the geometry of the airplanes in a solver object.

    animate: Create an animation of a problem's movement.

    plot_results_versus_time: This method takes in an unsteady solver object,
    and plots the geometries' forces, moments, force coefficients, and moment
    coefficients as a function of time.

    print_results: This function prints the forces, and force coefficients (in
    wind axes) and the moments, and moment coefficients (in wind axes, relative to
    the CG) calculated by a steady solver.
"""

import math
import time

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import webp

from . import steady_horseshoe_vortex_lattice_method
from . import steady_ring_vortex_lattice_method
from . import unsteady_ring_vortex_lattice_method

# Define the color and colormaps used by the visualization functions.
_sequential_color_map = "speed"
_diverging_color_map = "delta"
_wake_vortex_color = "white"
_panel_color = "chartreuse"
_streamline_color = "orchid"
_plotter_background_color = "black"
_figure_background_color = "None"
_text_color = "#818181"
_quality = 75
_window_size = [1024, 768]

# For the figure lines, use the "Prism" qualitative color map from
# carto.com/carto-colors.
_prism = [
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
    _drag_color,
    _side_color,
    _lift_color,
    _roll_color,
    _pitch_color,
    _yaw_color,
] = _prism[3:9]

# Set constants for the color maps, scalar bars, and text boxes.
_color_map_num_sig = 3
_bar_title_font_size = 30
_bar_label_font_size = 21
_bar_width = 0.5
_bar_position_x = 0.25
_bar_position_y = 0.05
_bar_n_labels = 2
_text_max_position = (0.85, 0.075)
_text_min_position = (0.85, 0.050)
_text_speed_position = (0.05, 0.075)
_text_font_size = 11

# Set the number of markers and the marker size for the results plots.
_num_markers = 6
_marker_size = 8

# Calculate the normalized spacing between the markers for the results plots.
_marker_spacing = 1.0 / _num_markers


# NOTE: I haven't yet started refactoring this function.
def draw(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
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
    plotter = pv.Plotter(window_size=_window_size, lighting=None)
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
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces(
                solver, draw_step
            )
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=False,
                color=_wake_vortex_color,
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
    panel_surfaces = _get_panel_surfaces(airplanes)

    # Check if the user wants to plot any scalars.
    if show_scalars:

        # Get the scalars
        these_scalars = _get_scalars(airplanes, scalar_type)
        min_scalar = round(min(these_scalars), 2)
        max_scalar = round(max(these_scalars), 2)

        # Choose the color map and set its limits based on if the min and max scalars
        # have the same sign (sequential color map) or if they have different signs
        # (diverging color map).
        if np.sign(np.min(these_scalars)) == np.sign(np.max(these_scalars)):
            color_map = _sequential_color_map
            c_min = max(
                np.mean(these_scalars) - _color_map_num_sig * np.std(these_scalars),
                np.min(these_scalars),
            )
            c_max = min(
                np.mean(these_scalars) + _color_map_num_sig * np.std(these_scalars),
                np.max(these_scalars),
            )
        else:
            color_map = _diverging_color_map
            c_min = -_color_map_num_sig * np.std(these_scalars)
            c_max = _color_map_num_sig * np.std(these_scalars)

        _plot_scalars(
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
            color=_panel_color,
            smooth_shading=False,
        )

    # Check if the user wants to plot streamlines.
    if show_streamlines:
        # Iterate through the spanwise positions in the solver's streamline point
        # matrix.
        for spanwise_position in range(solver.stackStreamlinePoints_G_Cg.shape[1]):

            # Get the column of streamline points at this spanwise position.
            streamline_point_column = solver.stackStreamlinePoints_G_Cg[
                :, spanwise_position, :
            ]

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
                        color=_streamline_color,
                        line_width=2,
                        smooth_shading=False,
                    )

    # Set the plotter's background color.
    plotter.set_background(color=_plotter_background_color)
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
            img=image, file_path="Draw.webp", lossless=False, quality=_quality
        )

    # Close all the plotters.
    pv.close_all()


# NOTE: I haven't yet started refactoring this function.
def animate(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
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
    plotter = pv.Plotter(window_size=_window_size, lighting=None)

    # Initialize values to hold the color map choice and its limits.
    c_min = 0
    c_max = 0
    color_map = None

    if save:
        # Add text to the animation that display's its speed relative to the true
        # speed.
        plotter.add_text(
            text="Speed: " + str(round(100 * speed)) + "%",
            position=_text_speed_position,
            font_size=_text_font_size,
            viewport=True,
            color=_text_color,
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
            scalars_to_add = _get_scalars(airplanes, scalar_type)
            all_scalars = np.hstack((all_scalars, scalars_to_add))

        # Choose the color map and set its limits based on if the min and max scalars
        # across all time steps have the same sign (sequential color map) or if they
        # have different signs (diverging color map).
        if np.sign(np.min(all_scalars)) == np.sign(np.max(all_scalars)):
            color_map = _sequential_color_map
            c_min = max(
                np.mean(all_scalars) - _color_map_num_sig * np.std(all_scalars),
                np.min(all_scalars),
            )
            c_max = min(
                np.mean(all_scalars) + _color_map_num_sig * np.std(all_scalars),
                np.max(all_scalars),
            )
        else:
            color_map = _diverging_color_map
            c_min = -_color_map_num_sig * np.std(all_scalars)
            c_max = _color_map_num_sig * np.std(all_scalars)

        min_scalar = round(min(all_scalars), 2)
        max_scalar = round(max(all_scalars), 2)

    # Initialize the panel surfaces and add the meshes to the plotter.
    panel_surfaces = _get_panel_surfaces(step_airplanes[0])

    # Check if the user wants to show any scalars. If so, add the panel surfaces to
    # the plotter with these scalars.
    if show_scalars and first_results_step == 0:
        these_scalars = _get_scalars(step_airplanes[0], scalar_type)

        _plot_scalars(
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
            color=_panel_color,
            smooth_shading=False,
        )

    # Set the plotter background color and show the plotter.
    plotter.set_background(color=_plotter_background_color)

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
        panel_surfaces = _get_panel_surfaces(airplanes)

        if save:
            plotter.add_text(
                text="Speed: " + str(round(100 * speed)) + "%",
                position=_text_speed_position,
                font_size=_text_font_size,
                viewport=True,
                color=_text_color,
            )

        # If the user wants to show the wake ring vortices, then get their surfaces
        # and plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces(
                unsteady_solver, current_step
            )
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=False,
                color=_wake_vortex_color,
            )

        # Check if the user wants to plot scalars and this step is equal to or
        # greater than the first step with calculated results. If so, add the panel
        # surfaces to the plotter with the scalars.
        if show_scalars and first_results_step <= current_step:

            these_scalars = _get_scalars(airplanes, scalar_type)

            _plot_scalars(
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
                color=_panel_color,
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
            images, "Animate.webp", fps=actual_fps, lossless=False, quality=_quality
        )

    # Close all the plotters.
    pv.close_all()


# NOTE: I haven't yet started refactoring this function.
def plot_results_versus_time(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    show=True,
    save=False,
):
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
    forces_W = np.zeros((num_airplanes, 3, num_steps_to_average))
    forceCoefficients_W = np.zeros((num_airplanes, 3, num_steps_to_average))
    moments_W_Cg = np.zeros((num_airplanes, 3, num_steps_to_average))
    momentCoefficients_W_Cg = np.zeros((num_airplanes, 3, num_steps_to_average))

    # Initialize a variable to track position in the results arrays.
    results_step = 0

    # Iterate through the time steps and add the results to their respective matrices.
    for step in range(first_results_step, num_steps):

        # Get the airplanes from the problem at this step.
        airplanes = unsteady_solver.steady_problems[step].airplanes

        # Iterate through this step's airplanes.
        for airplane_id, airplane in enumerate(airplanes):
            forces_W[airplane_id, :, results_step] = airplane.forces_W
            forceCoefficients_W[airplane_id, :, results_step] = (
                airplane.forceCoefficients_W
            )
            moments_W_Cg[airplane_id, :, results_step] = airplane.moments_W_Cg
            momentCoefficients_W_Cg[airplane_id, :, results_step] = (
                airplane.momentCoefficients_W_Cg
            )

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
        force_axes.spines.bottom.set_color(_text_color)
        force_axes.spines.left.set_color(_text_color)
        force_axes.xaxis.label.set_color(_text_color)
        force_axes.yaxis.label.set_color(_text_color)
        force_coefficients_axes.spines.bottom.set_color(_text_color)
        force_coefficients_axes.spines.left.set_color(_text_color)
        force_coefficients_axes.xaxis.label.set_color(_text_color)
        force_coefficients_axes.yaxis.label.set_color(_text_color)
        moment_coefficients_axes.spines.bottom.set_color(_text_color)
        moment_coefficients_axes.spines.left.set_color(_text_color)
        moment_coefficients_axes.xaxis.label.set_color(_text_color)
        moment_coefficients_axes.yaxis.label.set_color(_text_color)
        moment_axes.spines.bottom.set_color(_text_color)
        moment_axes.spines.left.set_color(_text_color)
        moment_axes.xaxis.label.set_color(_text_color)
        moment_axes.yaxis.label.set_color(_text_color)

        # Format all the plots' tick colors.
        force_axes.tick_params(axis="x", colors=_text_color)
        force_axes.tick_params(axis="y", colors=_text_color)
        force_coefficients_axes.tick_params(axis="x", colors=_text_color)
        force_coefficients_axes.tick_params(axis="y", colors=_text_color)
        moment_coefficients_axes.tick_params(axis="x", colors=_text_color)
        moment_coefficients_axes.tick_params(axis="y", colors=_text_color)
        moment_axes.tick_params(axis="x", colors=_text_color)
        moment_axes.tick_params(axis="y", colors=_text_color)

        # Format all the plots' background colors.
        force_figure.patch.set_facecolor(_figure_background_color)
        force_axes.set_facecolor(_figure_background_color)
        force_coefficients_figure.patch.set_facecolor(_figure_background_color)
        force_coefficients_axes.set_facecolor(_figure_background_color)
        moment_figure.patch.set_facecolor(_figure_background_color)
        moment_axes.set_facecolor(_figure_background_color)
        moment_coefficients_figure.patch.set_facecolor(_figure_background_color)
        moment_coefficients_axes.set_facecolor(_figure_background_color)

        # Populate the plots.
        force_axes.plot(
            times,
            forces_W[airplane_id, 0],
            label="Induced Drag",
            color=_drag_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_axes.plot(
            times,
            forces_W[airplane_id, 1],
            label="Side Force",
            color=_side_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_axes.plot(
            times,
            forces_W[airplane_id, 2],
            label="Lift",
            color=_lift_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            forceCoefficients_W[airplane_id, 0],
            label="Induced Drag",
            color=_drag_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            forceCoefficients_W[airplane_id, 1],
            label="Side Force",
            color=_side_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            forceCoefficients_W[airplane_id, 2],
            label="Lift",
            color=_lift_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_Cg[airplane_id, 0],
            label="Roll",
            color=_roll_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_Cg[airplane_id, 1],
            label="Pitch",
            color=_pitch_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_Cg[airplane_id, 2],
            label="Yaw",
            color=_yaw_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_Cg[airplane_id, 0],
            label="Roll",
            color=_roll_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_Cg[airplane_id, 1],
            label="Pitch",
            color=_pitch_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_Cg[airplane_id, 2],
            label="Yaw",
            color=_yaw_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )

        # Find and format this airplane's name for use in the plot titles.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[airplane_id].name
        force_title = airplane_name + " Forces vs. Time"
        force_coefficient_title = airplane_name + " Force Coefficients vs. Time"
        moment_title = airplane_name + " Moments vs. Time"
        moment_coefficient_title = airplane_name + " Moment Coefficients vs. Time"

        # Name the plots' axis labels and titles.
        force_axes.set_xlabel("Time (s)", color=_text_color)
        force_axes.set_ylabel("Force (N)", color=_text_color)
        force_axes.set_title(force_title, color=_text_color)
        force_coefficients_axes.set_xlabel("Time (s)", color=_text_color)
        force_coefficients_axes.set_ylabel("Coefficient", color=_text_color)
        force_coefficients_axes.set_title(force_coefficient_title, color=_text_color)
        moment_axes.set_xlabel("Time (s)", color=_text_color)
        moment_axes.set_ylabel("Moment (N m)", color=_text_color)
        moment_axes.set_title(moment_title, color=_text_color)
        moment_coefficients_axes.set_xlabel("Time (s)", color=_text_color)
        moment_coefficients_axes.set_ylabel("Coefficient", color=_text_color)
        moment_coefficients_axes.set_title(moment_coefficient_title, color=_text_color)

        # Format the plots' legends.
        force_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color,
        )
        force_coefficients_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color,
        )
        moment_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color,
        )
        moment_coefficients_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color,
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


# TEST: Add unit tests for this method.
def print_results(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
):
    """This function prints the forces, and force coefficients (in wind axes) and the
    moments, and moment coefficients (in wind axes, relative to the CG) calculated by
    a solver.

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver

        This is the solver with the results to be printed.

    :return: None
    """
    if isinstance(
        solver,
        (
            steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
        ),
    ):
        these_airplanes = solver.airplanes
        solver_type = "steady"
    elif isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        these_airplanes = solver.current_airplanes
        if solver.unsteady_problem.movement.static:
            solver_type = "static geometry unsteady"
        else:
            solver_type = "variable geometry unsteady"
    else:
        raise TypeError(
            "solver must be a SteadyHorseshoeVortexLatticeMethodSolver, a SteadyRingVortexLatticeMethodSolver, or an UnsteadyRingVortexLatticeMethodSolver."
        )

    padding_spaces = 2

    col1 = [
        "FX_W",
        "FY_W",
        "FZ_W",
        "MX_W_Cg",
        "MY_W_Cg",
        "MZ_W_Cg",
        "cFX_W",
        "cFY_W",
        "cFZ_W",
        "cMX_W_Cg",
        "cMY_W_Cg",
        "cMZ_W_Cg",
    ]
    col1 = [label + ":" for label in col1]
    col1_space = max(len(elem) for elem in col1) + padding_spaces

    col3 = [
        "Drag",
        "Side Force",
        "Lift",
        "Rolling Moment",
        "Pitching Moment",
        "Yawing Moment",
        "CDi",
        "CY",
        "CL",
        "Cl",
        "Cm",
        "Cn",
    ]
    col3 = [label + ":" for label in col3]
    col3_space = max(len(elem) for elem in col3) + padding_spaces

    pad = " " * padding_spaces

    for airplane_num, airplane in enumerate(these_airplanes):
        title1 = None
        title2 = None
        title3 = None
        title4 = None
        these_forces_W = None
        these_moments_W_Cg = None
        these_forceCoefficients_W = None
        these_momentCoefficients_W_Cg = None

        match solver_type:
            case "steady":
                title1 = f"{pad}Forces (in wind axes):"
                title2 = f"{pad}Moments (in wind axes, relative to the CG):"
                title3 = f"{pad}Force Coefficients (in wind axes):"
                title4 = f"{pad}Moment Coefficients (in wind axes, relative to the CG):"
                these_forces_W = airplane.forces_W
                these_moments_W_Cg = airplane.moments_W_Cg
                these_forceCoefficients_W = airplane.forceCoefficients_W
                these_momentCoefficients_W_Cg = airplane.momentCoefficients_W_Cg
            case "static geometry unsteady":
                title1 = f"{pad}Final Forces (in wind axes):"
                title2 = f"{pad}Final Moments (in wind axes, relative to the CG):"
                title3 = f"{pad}Final Force Coefficients (in wind axes):"
                title4 = f"{pad}Final Moment Coefficients (in wind axes, relative to the CG):"
                these_forces_W = solver.unsteady_problem.finalForces_W[airplane_num]
                these_moments_W_Cg = solver.unsteady_problem.finalMoments_W_Cg[
                    airplane_num
                ]
                these_forceCoefficients_W = (
                    solver.unsteady_problem.finalForceCoefficients_W[airplane_num]
                )
                these_momentCoefficients_W_Cg = (
                    solver.unsteady_problem.finalMomentCoefficients_W_Cg[airplane_num]
                )
            case "variable geometry unsteady":
                title1 = f"{pad}Final Cycle-Averaged Forces (in wind axes):"
                title2 = f"{pad}Final Cycle-Averaged Moments (in wind axes, relative to the CG):"
                title3 = f"{pad}Final Cycle-Averaged Force Coefficients (in wind axes):"
                title4 = f"{pad}Final Cycle-Averaged Moment Coefficients (in wind axes, relative to the CG):"
                these_forces_W = solver.unsteady_problem.finalMeanForces_W[airplane_num]
                these_moments_W_Cg = solver.unsteady_problem.finalMeanMoments_W_Cg[
                    airplane_num
                ]
                these_forceCoefficients_W = (
                    solver.unsteady_problem.finalMeanForceCoefficients_W[airplane_num]
                )
                these_momentCoefficients_W_Cg = (
                    solver.unsteady_problem.finalMeanMomentCoefficients_W_Cg[
                        airplane_num
                    ]
                )
            case _:
                raise ValueError(f"Unknown solver type: {solver_type}")

        col2 = [
            these_forces_W[0],
            these_forces_W[1],
            these_forces_W[2],
            these_moments_W_Cg[0],
            these_moments_W_Cg[1],
            these_moments_W_Cg[2],
            these_forceCoefficients_W[0],
            these_forceCoefficients_W[1],
            these_forceCoefficients_W[2],
            these_momentCoefficients_W_Cg[0],
            these_momentCoefficients_W_Cg[1],
            these_momentCoefficients_W_Cg[2],
        ]
        col2 = [str(np.round(val, 3)) for val in col2]
        col2 = [
            val + " N" if i < 3 else val + " Nm" if i < 6 else val
            for i, val in enumerate(col2)
        ]
        col2_space = max(len(elem) for elem in col2) + 2 * padding_spaces

        col4 = [
            -these_forces_W[0],
            these_forces_W[1],
            -these_forces_W[2],
            these_moments_W_Cg[0],
            these_moments_W_Cg[1],
            these_moments_W_Cg[2],
            -these_forceCoefficients_W[0],
            these_forceCoefficients_W[1],
            -these_forceCoefficients_W[2],
            these_momentCoefficients_W_Cg[0],
            these_momentCoefficients_W_Cg[1],
            these_momentCoefficients_W_Cg[2],
        ]
        col4 = [str(np.round(val, 3)) for val in col4]
        col4 = [
            val + " N" if i < 3 else val + " Nm" if i < 6 else val
            for i, val in enumerate(col4)
        ]

        print(f'Airplane "{airplane.name}":')

        for i in range(len(col1)):
            if i == 0:
                print(title1)
            elif i == 3:
                print(title2)
            elif i == 6:
                print(title3)
            elif i == 9:
                print(title4)

            s = f"{2 * pad}{col1[i]:<{col1_space}}{col2[i]:<{col2_space}}{col3[i]:<{col3_space}}{col4[i]}"
            print(s)

        # If the results from more Airplanes are going to be printed, print two new
        # lines to separate them.
        if (airplane_num + 1) < solver.num_airplanes:
            print("\n")


# NOTE: I haven't yet started refactoring this function.
def _get_panel_surfaces(
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


# NOTE: I haven't yet started refactoring this function.
def _get_wake_ring_vortex_surfaces(solver, step):
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
    num_wake_ring_vortices = solver.list_num_wake_vortices[step]
    wake_ring_vortex_front_right_vertices = solver.listStackFrwrvp_G_Cg[step]
    wake_ring_vortex_front_left_vertices = solver.listStackFlwrvp_G_Cg[step]
    wake_ring_vortex_back_left_vertices = solver.listStackBlwrvp_G_Cg[step]
    wake_ring_vortex_back_right_vertices = solver.listStackBrwrvp_G_Cg[step]

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


# NOTE: I haven't yet started refactoring this function.
def _get_scalars(
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


# NOTE: I haven't yet started refactoring this function.
def _plot_scalars(
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
        title_font_size=_bar_title_font_size,
        label_font_size=_bar_label_font_size,
        width=_bar_width,
        position_x=_bar_position_x,
        position_y=_bar_position_y,
        n_labels=_bar_n_labels,
        fmt="%.2f",
        color=_text_color,
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
        position=_text_max_position,
        font_size=_text_font_size,
        viewport=True,
        color=_text_color,
    )
    plotter.add_text(
        text="Min: " + str(min_scalar),
        position=_text_min_position,
        font_size=_text_font_size,
        viewport=True,
        color=_text_color,
    )
