"""Contains useful functions for visualizing geometry and results.

This module contains the following classes:
    None

This module contains the following functions:
    draw: Draw a solver's Airplane(s).

    animate: Animate an UnsteadyRingVortexLatticeMethodSolver's Airplane(s).

    plot_results_versus_time: This function takes in an
    UnsteadyRingVortexLatticeMethodSolver, and plots the loads and load coefficients
    as a function of time.

    print_results: This function prints the load and load coefficients calculated by
    a solver.
"""

from __future__ import annotations

import math
import time
from typing import cast

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import webp

from . import geometry
from . import _parameter_validation
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
_quality = 75.0
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


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
def draw(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    scalar_type: str | None = None,
    show_streamlines: bool | np.bool_ = False,
    show_wake_vortices: bool | np.bool_ = False,
    save: bool | np.bool_ = False,
    testing: bool | np.bool_ = False,
) -> None:
    """Draw a solver's Airplane(s).

    Citation:
        Adapted from:         vlm3.draw in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver

        This is the solver whose Airplane(s) will be plotted.

    :param scalar_type: str or None, optional

        This variable determines how to color the Panels. Setting this to None colors
        the Panels uniformly. It can also be "induced drag", "side force", or "lift",
        which respectively use each Panel's induced drag, side force, and lift
        coefficient. The default value is None.

    :param show_streamlines: boolLike, optional

        Set this to True to show the streamlines emanating from the back of the
        Wings. It can be a bool or a numpy bool and will be converted
        internally to a bool. The default value is False.

    :param show_wake_vortices: boolLike, optional

        Set this to True to show any wake RingVortices. It can be a bool or a
        numpy bool and will be converted internally to a bool. The default
        value is False.

    :param save: boolLike, optional

        Set this to True to save the image as a WebP. It can be a bool or a NumPy
        bool and will be converted internally to a bool. The default value is
        False.

    :param testing: boolLike, optional

        Set this to True to close the image after 1 second, which is useful for
        running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default value is False.

    :return: None
    """
    if not isinstance(
        solver,
        (
            steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
            unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        ),
    ):
        raise TypeError(
            "solver must be a SteadyHorseshoeVortexLatticeMethodSolver, "
            "SteadyRingVortexLatticeMethodSolver, "
            "or UnsteadyRingVortexLatticeMethodSolver."
        )

    if scalar_type is not None:
        if not solver.ran:
            raise RuntimeError(
                "solver must have run before drawing with scalar_type not None."
            )

        scalar_type = _parameter_validation.str_return_str(scalar_type, "scalar_type")
        if scalar_type not in ("induced drag", "side force", "lift"):
            raise ValueError(
                'scalar_type must be None, "induced drag", "side force", or "lift".'
            )

    show_streamlines = _parameter_validation.boolLike_return_bool(
        show_streamlines, "show_streamlines"
    )
    if show_streamlines and not solver.ran:
        raise RuntimeError(
            "solver must have run before drawing with show_streamlines set to True."
        )
    if show_streamlines and len(solver.gridStreamlinePoints_GP1_CgP1) == 0:
        raise RuntimeError(
            "solver must have streamline points calculated before drawing with "
            "show_streamlines set to True."
        )

    show_wake_vortices = _parameter_validation.boolLike_return_bool(
        show_wake_vortices, "show_wake_vortices"
    )
    if show_wake_vortices and not isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        raise ValueError(
            "show_wake_vortices can only be True when drawing an "
            "UnsteadyRingVortexLatticeMethodSolver."
        )
    if show_wake_vortices and not solver.ran:
        raise RuntimeError(
            "solver must have run before drawing with show_wake_vortices set to True."
        )

    save = _parameter_validation.boolLike_return_bool(save, "save")
    testing = _parameter_validation.boolLike_return_bool(testing, "testing")

    # Create the Plotter and set it to use parallel projection (instead of perspective).
    plotter = pv.Plotter(window_size=_window_size, lighting=None)
    plotter.enable_parallel_projection()  # type: ignore[call-arg]

    # Get the solver's geometry.
    if isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        draw_step = solver.num_steps - 1

        airplanes = solver.steady_problems[draw_step].airplanes
        qInf__E = solver.steady_problems[draw_step].operating_point.qInf__E

        # If showing wake RingVortices, get their surfaces and plot them.
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
        qInf__E = solver.operating_point.qInf__E

    # Get the Panel surfaces.
    panel_surfaces = _get_panel_surfaces(airplanes)

    # Plot the Panels either with scalar coloring or with a uniform color.
    if scalar_type in ("induced drag", "side force", "lift"):
        these_scalars = _get_scalars(airplanes, scalar_type, qInf__E)
        min_scalar = round(min(these_scalars), 2)
        max_scalar = round(max(these_scalars), 2)

        # Choose the color map and set its limits based on if the min and max scalars
        # have the same sign (sequential color map) or if they have different signs
        # (diverging color map).
        if np.sign(np.min(these_scalars)) == np.sign(np.max(these_scalars)):
            color_map = _sequential_color_map
            c_min = max(
                float(np.mean(these_scalars))
                - _color_map_num_sig * float(np.std(these_scalars)),
                float(np.min(these_scalars)),
            )
            c_max = min(
                float(np.mean(these_scalars))
                + _color_map_num_sig * float(np.std(these_scalars)),
                float(np.max(these_scalars)),
            )
        else:
            color_map = _diverging_color_map
            c_min = -_color_map_num_sig * float(np.std(these_scalars))
            c_max = _color_map_num_sig * float(np.std(these_scalars))

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

    # If showing streamlines, plot them.
    if show_streamlines:
        # Iterate through the spanwise positions in the solver's streamline point
        # ndarray.
        for spanwise_position in range(solver.gridStreamlinePoints_GP1_CgP1.shape[1]):
            # Get the ndarray of streamline points at this spanwise position (in the
            # first Airplane's geometry axes, relative to the first Airplane's CG).
            stackStreamlinePoints_GP1_CgP1 = solver.gridStreamlinePoints_GP1_CgP1[
                :, spanwise_position, :
            ]

            # Iterate through the streamline points at this spanwise position.
            for point_index in range(stackStreamlinePoints_GP1_CgP1.shape[0]):

                # Skip the first point because it has no previous point with which
                # to make a line.
                if point_index != 0:
                    # Get the current and last point.
                    point = stackStreamlinePoints_GP1_CgP1[point_index, :]
                    last_point = stackStreamlinePoints_GP1_CgP1[point_index - 1, :]

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

    # Set the Plotter's background color.
    plotter.set_background(color=_plotter_background_color)  # type: ignore[call-arg]
    if not testing:
        # Show the Plotter so the user can adjust the camera position and window.
        # When the user closes the window, the Plotter still exists. Therefore,
        # it can later be saved as an image if desired.
        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )
    else:
        # Show the Plotter for 1 second, then proceed automatically. This is useful
        # for testing.
        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # If saving, take a screenshot, convert it to a ndarray, convert that to an Image,
    # and save it as a WebP.
    if save:
        image = webp.Image.fromarray(
            np.array(
                plotter.screenshot(
                    filename=None,
                    transparent_background=True,
                    return_img=True,
                )
            )
        )

        webp.save_image(
            img=image, file_path="Draw.webp", lossless=False, quality=_quality
        )

    # Close all Plotters.
    pv.close_all()


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
def animate(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    scalar_type: str | None = None,
    show_wake_vortices: bool | np.bool_ = False,
    save: bool | np.bool_ = False,
    testing: bool | np.bool_ = False,
) -> None:
    """Animate an UnsteadyRingVortexLatticeMethodSolver's Airplane(s).

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver

        This is the UnsteadyRingVortexLatticeMethodSolver whose Airplane(s) will be
        animated.

    :param scalar_type: str or None, optional

        This variable determines how to color the Panels. Setting this to None colors
        the Panels uniformly. It can also be "induced drag", "side force", or "lift",
        which respectively use each Panel's induced drag, side force, and lift
        coefficient. The default value is None.

    :param show_wake_vortices: boolLike, optional

        Set this to True to show any wake RingVortices. It can be a bool or a
        numpy bool and will be converted internally to a bool. The default
        value is False.

    :param save: boolLike, optional

        Set this to True to save the image as a WebP. It can be a bool or a NumPy
        bool and will be converted internally to a bool. The default value is
        False.

    :param testing: boolLike, optional

        Set this to True to close the image after 1 second, which is useful for
        running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default value is False.

    :return: None
    """
    if not isinstance(
        unsteady_solver,
        (unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,),
    ):
        raise TypeError(
            "unsteady_solver must be an UnsteadyRingVortexLatticeMethodSolver."
        )

    if scalar_type is not None:
        if not unsteady_solver.ran:
            raise RuntimeError(
                "unsteady_solver must have run before animating with scalar_type not "
                "None."
            )

        scalar_type = _parameter_validation.str_return_str(scalar_type, "scalar_type")
        if scalar_type not in ("induced drag", "side force", "lift"):
            raise ValueError(
                'scalar_type must be None, "induced drag", "side force", or "lift".'
            )

    show_wake_vortices = _parameter_validation.boolLike_return_bool(
        show_wake_vortices, "show_wake_vortices"
    )
    if show_wake_vortices and not unsteady_solver.ran:
        raise RuntimeError(
            "unsteady_solver must have run before animating with show_wake_vortices set"
            " to True."
        )

    save = _parameter_validation.boolLike_return_bool(save, "save")
    testing = _parameter_validation.boolLike_return_bool(testing, "testing")

    first_results_step = unsteady_solver.first_results_step

    # Get the solver's SteadyProblems' Airplanes. This will become a list of lists,
    # with the first index being the time step and the second index identifying each
    # Airplane at that time step.
    step_airplanes = []
    for steady_problem in unsteady_solver.steady_problems:
        step_airplanes.append(steady_problem.airplanes)

    # Scale down the true-speed frames per second to at most 50 fps. This is the
    # maximum speed at which some programs can render WebPs.
    requested_fps = 1.0 / unsteady_solver.delta_time
    speed = 1.0
    if requested_fps > 50.0:
        speed = 50.0 / requested_fps
    actual_fps = float(math.floor(requested_fps * speed))

    # Create the Plotter and set it to use parallel projection (instead of perspective).
    plotter = pv.Plotter(window_size=_window_size, lighting=None)
    plotter.enable_parallel_projection()  # type: ignore[call-arg]

    # Initialize values to hold the color map choice and its limits.
    c_min = 0.0
    c_max = 0.0
    color_map: str = ""

    # If saving the animation, add text that displays its speed.
    if save:
        plotter.add_text(
            text="Speed: " + str(round(100 * speed)) + "%",
            position=_text_speed_position,
            font_size=_text_font_size,
            viewport=True,
            color=_text_color,
        )

    # Initialize variables to hold the SteadyProblems' scalars and their attributes.
    all_scalars = np.empty(0, dtype=float)
    min_scalar = 0.0
    max_scalar = 0.0

    # If coloring the Panels based on scalars, gather all the scalars across all the
    # time steps and Airplanes. These will be used to set the color map limits.
    if scalar_type is not None:
        for step_id, airplanes in enumerate(step_airplanes):
            scalars_to_add = _get_scalars(
                airplanes,
                scalar_type,
                unsteady_solver.steady_problems[step_id].operating_point.qInf__E,
            )
            all_scalars = cast(np.ndarray, np.hstack((all_scalars, scalars_to_add)))

        # Choose the color map and set its limits based on if the min and max scalars
        # across all time steps have the same sign (sequential color map) or if they
        # have different signs (diverging color map).
        if np.sign(np.min(all_scalars)) == np.sign(np.max(all_scalars)):
            color_map = _sequential_color_map
            c_min = max(
                float(np.mean(all_scalars))
                - _color_map_num_sig * float(np.std(all_scalars)),
                float(np.min(all_scalars)),
            )
            c_max = min(
                float(np.mean(all_scalars))
                + _color_map_num_sig * float(np.std(all_scalars)),
                float(np.max(all_scalars)),
            )
        else:
            color_map = _diverging_color_map
            c_min = -_color_map_num_sig * float(np.std(all_scalars))
            c_max = _color_map_num_sig * float(np.std(all_scalars))

        min_scalar = round(min(all_scalars), 2)
        max_scalar = round(max(all_scalars), 2)

    # Get the Panel surfaces of the first time step's Airplane(s).
    panel_surfaces = _get_panel_surfaces(step_airplanes[0])

    # Plot the first time step's Airplanes' Panels either with scalar coloring or
    # with a uniform color.
    if scalar_type is not None and first_results_step == 0:
        these_scalars = _get_scalars(
            step_airplanes[0],
            scalar_type,
            unsteady_solver.steady_problems[0].operating_point.qInf__E,
        )

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

    # Set the Plotter's background color.
    plotter.set_background(color=_plotter_background_color)  # type: ignore[call-arg]

    # If not testing, show the Plotter with the first time step, and print a message
    # to the console on how to adjust the view and start the animation. If testing,
    # show the Plotter with the first time step for 1 second, and start the animation
    # with the current window view.
    if not testing:
        print(
            'Orient the view, then press "q" to close the window and produce the '
            "animation."
        )
        plotter.show(
            title="Rendering speed not to scale.",
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )
    else:
        plotter.show(
            title="Rendering speed not to scale.",
            cpos=(-1, -1, 1),
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # Start a list to hold a WebP Image of each frame. To start, take a screenshot,
    # convert it to a ndarray, and convert that to an Image.
    images = [
        webp.Image.fromarray(
            np.array(
                plotter.screenshot(
                    transparent_background=True,
                    return_img=True,
                )
            )
        )
    ]

    # Initialize a variable to keep track of the current time step.
    current_step = 1

    # Begin to iterate through the Airplane(s) from the subsequent time steps.
    for airplanes in step_airplanes[1:]:

        # Clear the Plotter.
        plotter.clear()

        # Get the Panel surfaces of this time step's Airplane(s).
        panel_surfaces = _get_panel_surfaces(airplanes)

        # If saving the animation, add text that displays its speed.
        if save:
            plotter.add_text(
                text="Speed: " + str(round(100 * speed)) + "%",
                position=_text_speed_position,
                font_size=_text_font_size,
                viewport=True,
                color=_text_color,
            )

        # If showing wake RingVortices, get their surfaces and plot them.
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

        # Plot the Panels either with a uniform color or, if the current time step
        # has results, with scalar coloring.
        if scalar_type is not None and first_results_step <= current_step:
            these_scalars = _get_scalars(
                airplanes,
                scalar_type,
                unsteady_solver.steady_problems[current_step].operating_point.qInf__E,
            )

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

        # If saving, append a WebP Image of this frame to the list of Images. To do
        # so, take a screenshot, convert it to a ndarray, and convert that to an Image.
        if save:
            images.append(
                webp.Image.fromarray(
                    np.array(
                        plotter.screenshot(
                            filename=None,
                            transparent_background=True,
                            return_img=True,
                        )
                    )
                )
            )

        # Increment the time step tracker.
        current_step += 1

    # If saving, save the list of Images as an animated WebP.
    if save:
        # Convert the list of WebP Images to an WebP animation.
        webp.save_images(
            images, "Animate.webp", fps=actual_fps, lossless=False, quality=_quality
        )

    # Close all the Plotters.
    pv.close_all()


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
def plot_results_versus_time(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    show: bool | np.bool_ = True,
    save: bool | np.bool_ = False,
) -> None:
    """This function takes in an UnsteadyRingVortexLatticeMethodSolver, and plots the
    loads and load coefficients as a function of time.

    :param unsteady_solver: UnsteadyRingVortexLatticeMethodSolver This is the
        UnsteadyRingVortexLatticeMethodSolver whose loads will be plotted.
    :param show: boolLike, optional Set this to True to show the plots. It can be a
        bool or a numpy bool and will be converted internally to a bool. The
        default is True.
    :param save: boolLike, Optional Set this to True to save the plots as PNGs. It can
        be a bool or a numpy bool and will be converted internally to a bool.
        The default is True.
    :return: None
    """
    if not isinstance(
        unsteady_solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        raise TypeError(
            "unsteady_solver must be an " "UnsteadyRingVortexLatticeMethodSolver."
        )
    show = _parameter_validation.boolLike_return_bool(show, "show")
    save = _parameter_validation.boolLike_return_bool(save, "save")

    if not unsteady_solver.ran:
        raise RuntimeError(
            "unsteady_solver must have run before plotting results versus time."
        )

    first_results_step = unsteady_solver.first_results_step

    # Get the time step characteristics. Note that the first time step (time step
    # 0), occurs at 0 seconds.
    num_steps = unsteady_solver.num_steps
    delta_time = unsteady_solver.delta_time
    num_airplanes = unsteady_solver.num_airplanes
    first_results_time_step_time = delta_time * first_results_step
    final_time_step_time = delta_time * (num_steps - 1)
    num_steps_to_average = num_steps - first_results_step

    # Create a 1D ndarray with the time at each time step where results have been
    # calculated.
    times = np.linspace(
        first_results_time_step_time,
        final_time_step_time,
        num_steps_to_average,
        endpoint=True,
    )

    # Initialize matrices to hold the loads and load coefficients at every time step
    # that has results.
    forces_W = np.zeros((num_airplanes, 3, num_steps_to_average), dtype=float)
    forceCoefficients_W = np.zeros(
        (num_airplanes, 3, num_steps_to_average), dtype=float
    )
    moments_W_CgP1 = np.zeros((num_airplanes, 3, num_steps_to_average), dtype=float)
    momentCoefficients_W_CgP1 = np.zeros(
        (num_airplanes, 3, num_steps_to_average), dtype=float
    )

    # Initialize a variable to track position in the results arrays.
    results_step = 0

    # Iterate through the time steps and add the results to their respective matrices.
    for step in range(first_results_step, num_steps):

        # Get the Airplanes from the SteadyProblem at this time step.
        airplanes = unsteady_solver.steady_problems[step].airplanes

        # Iterate through this time step's Airplanes.
        for airplane_id, airplane in enumerate(airplanes):
            forces_W[airplane_id, :, results_step] = airplane.forces_W
            forceCoefficients_W[airplane_id, :, results_step] = (
                airplane.forceCoefficients_W
            )
            moments_W_CgP1[airplane_id, :, results_step] = airplane.moments_W_CgP1
            momentCoefficients_W_CgP1[airplane_id, :, results_step] = (
                airplane.momentCoefficients_W_CgP1
            )

        results_step += 1

    # Iterate through the Airplane ID's to plot each Airplane's figures.
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
            -forces_W[airplane_id, 0],
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
            -forces_W[airplane_id, 2],
            label="Lift",
            color=_lift_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            -forceCoefficients_W[airplane_id, 0],
            label="Induced Drag Coefficient",
            color=_drag_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            forceCoefficients_W[airplane_id, 1],
            label="Side Force Coefficient",
            color=_side_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            -forceCoefficients_W[airplane_id, 2],
            label="Lift Coefficient",
            color=_lift_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_CgP1[airplane_id, 0],
            label="Roll",
            color=_roll_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_CgP1[airplane_id, 1],
            label="Pitch",
            color=_pitch_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_CgP1[airplane_id, 2],
            label="Yaw",
            color=_yaw_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_CgP1[airplane_id, 0],
            label="Roll Coefficient",
            color=_roll_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_CgP1[airplane_id, 1],
            label="Pitch Coefficient",
            color=_pitch_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_CgP1[airplane_id, 2],
            label="Yaw Coefficient",
            color=_yaw_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )

        # Find and format this Airplane's name for use in the plot titles.
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
        force_coefficients_axes.set_ylabel("Force Coefficient", color=_text_color)
        force_coefficients_axes.set_title(force_coefficient_title, color=_text_color)
        moment_axes.set_xlabel("Time (s)", color=_text_color)
        moment_axes.set_ylabel("Moment (N m)", color=_text_color)
        moment_axes.set_title(moment_title, color=_text_color)
        moment_coefficients_axes.set_xlabel("Time (s)", color=_text_color)
        moment_coefficients_axes.set_ylabel("Moment Coefficient", color=_text_color)
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

    # If the user wants to show the plots, do so. This is done outside the loop so
    # that plt.show() is only called once after all figures are created.
    if show:
        plt.show()
    else:
        plt.close("all")


# TEST: Consider adding unit tests for this function.
# TEST: Consider adding integration tests for this function.
def print_results(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
) -> None:
    """This function prints the load and load coefficients calculated by a solver.

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
            "solver must be a SteadyHorseshoeVortexLatticeMethodSolver, "
            "a SteadyRingVortexLatticeMethodSolver, or an "
            "UnsteadyRingVortexLatticeMethodSolver."
        )

    if not solver.ran:
        raise RuntimeError("solver must have run before printing results.")

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
        title1: str = ""
        title2: str = ""
        title3: str = ""
        title4: str = ""
        these_forces_W: np.ndarray = np.empty(0, dtype=float)
        these_moments_W_CgP1: np.ndarray = np.empty(0, dtype=float)
        these_forceCoefficients_W: np.ndarray = np.empty(0, dtype=float)
        these_momentCoefficients_W_CgP1: np.ndarray = np.empty(0, dtype=float)

        match solver_type:
            case "steady":
                title1 = f"{pad}Forces (in wind axes):"
                title2 = f"{pad}Moments (in wind axes, relative to the CG):"
                title3 = f"{pad}Force Coefficients (in wind axes):"
                title4 = f"{pad}Moment Coefficients (in wind axes, relative to the CG):"

                _forces_W = airplane.forces_W
                assert _forces_W is not None

                these_forces_W = _forces_W

                _moments_W_CgP1 = airplane.moments_W_CgP1
                assert _moments_W_CgP1 is not None

                these_moments_W_CgP1 = _moments_W_CgP1

                _forceCoefficients_W = airplane.forceCoefficients_W
                assert _forceCoefficients_W is not None

                these_forceCoefficients_W = _forceCoefficients_W

                _momentCoefficients_W_CgP1 = airplane.momentCoefficients_W_CgP1
                assert _momentCoefficients_W_CgP1 is not None

                these_momentCoefficients_W_CgP1 = _momentCoefficients_W_CgP1

            case "static geometry unsteady":
                assert isinstance(
                    solver,
                    unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
                )

                title1 = f"{pad}Final Forces (in wind axes):"
                title2 = f"{pad}Final Moments (in wind axes, relative to the CG):"
                title3 = f"{pad}Final Force Coefficients (in wind axes):"
                title4 = (
                    f"{pad}Final Moment Coefficients (in wind axes, relative to "
                    f"the CG):"
                )
                these_forces_W = solver.unsteady_problem.finalForces_W[airplane_num]
                these_moments_W_CgP1 = solver.unsteady_problem.finalMoments_W_CgP1[
                    airplane_num
                ]
                these_forceCoefficients_W = (
                    solver.unsteady_problem.finalForceCoefficients_W[airplane_num]
                )
                these_momentCoefficients_W_CgP1 = (
                    solver.unsteady_problem.finalMomentCoefficients_W_CgP1[airplane_num]
                )
            case "variable geometry unsteady":
                assert isinstance(
                    solver,
                    unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
                )

                title1 = f"{pad}Final Cycle-Averaged Forces (in wind axes):"
                title2 = (
                    f"{pad}Final Cycle-Averaged Moments (in wind axes, relative "
                    f"to the CG):"
                )
                title3 = f"{pad}Final Cycle-Averaged Force Coefficients (in wind axes):"
                title4 = (
                    f"{pad}Final Cycle-Averaged Moment Coefficients (in wind "
                    f"axes, relative to the CG):"
                )
                these_forces_W = solver.unsteady_problem.finalMeanForces_W[airplane_num]
                these_moments_W_CgP1 = solver.unsteady_problem.finalMeanMoments_W_CgP1[
                    airplane_num
                ]
                these_forceCoefficients_W = (
                    solver.unsteady_problem.finalMeanForceCoefficients_W[airplane_num]
                )
                these_momentCoefficients_W_CgP1 = (
                    solver.unsteady_problem.finalMeanMomentCoefficients_W_CgP1[
                        airplane_num
                    ]
                )
            case _:
                raise ValueError(f"Unknown solver type: {solver_type}")

        col2 = [
            these_forces_W[0],
            these_forces_W[1],
            these_forces_W[2],
            these_moments_W_CgP1[0],
            these_moments_W_CgP1[1],
            these_moments_W_CgP1[2],
            these_forceCoefficients_W[0],
            these_forceCoefficients_W[1],
            these_forceCoefficients_W[2],
            these_momentCoefficients_W_CgP1[0],
            these_momentCoefficients_W_CgP1[1],
            these_momentCoefficients_W_CgP1[2],
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
            these_moments_W_CgP1[0],
            these_moments_W_CgP1[1],
            these_moments_W_CgP1[2],
            -these_forceCoefficients_W[0],
            these_forceCoefficients_W[1],
            -these_forceCoefficients_W[2],
            these_momentCoefficients_W_CgP1[0],
            these_momentCoefficients_W_CgP1[1],
            these_momentCoefficients_W_CgP1[2],
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


# TEST: Consider adding unit tests for this function.
def _get_panel_surfaces(
    airplanes: list[geometry.airplane.Airplane],
) -> pv.PolyData:
    """This function returns a PolyData representation of the Wings' Panels' surfaces
    associated with all the Airplanes in a list of Airplanes.

    :param airplanes: list of Airplanes This is a list of Airplanes whose Wings' Panels'
        surfaces will be returned.
    :return: pv.PolyData This is a PolyData representation of the Airplanes' Wings'
        Panels' surfaces.
    """
    # Initialize empty ndarrays to hold the Panels' vertices and faces.
    panel_vertices = np.empty((0, 3), dtype=float)
    panel_faces = np.empty(0, dtype=int)

    # Initialize a variable to keep track of how many Panels have been added thus far.
    panel_num = 0

    # Increment through the Airplanes' Wing(s).
    for airplane in airplanes:
        for wing in airplane.wings:
            _panels = wing.panels
            assert _panels is not None

            # Unravel this Wing's ndarray of Panels iterate through it.
            panels = np.ravel(_panels)
            for panel in panels:
                # Arrange this Panel's vertices and faces into ndarrays in the
                # proper form to represent PolyData surfaces.
                panel_vertices_to_add = np.vstack(
                    (
                        panel.Flpp_GP1_CgP1,
                        panel.Frpp_GP1_CgP1,
                        panel.Brpp_GP1_CgP1,
                        panel.Blpp_GP1_CgP1,
                    )
                )
                panel_face_to_add = np.array(
                    [
                        4,
                        (panel_num * 4),
                        (panel_num * 4) + 1,
                        (panel_num * 4) + 2,
                        (panel_num * 4) + 3,
                    ],
                    dtype=int,
                )

                # Add this Panel's vertices and faces to the ndarray of all vertices
                # and faces.
                panel_vertices = cast(
                    np.ndarray, np.vstack((panel_vertices, panel_vertices_to_add))
                )
                panel_faces = cast(
                    np.ndarray, np.hstack((panel_faces, panel_face_to_add))
                )

                # Update the number of Panels.
                panel_num += 1

    # Return the Panels' surfaces.
    return pv.PolyData(panel_vertices, panel_faces)


# TEST: Consider adding unit tests for this function.
def _get_wake_ring_vortex_surfaces(
    solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    step: int,
) -> pv.PolyData:
    """This function returns the PolyData representation of surfaces an
    UnsteadyRingVortexLatticeMethodSolver's wake RingVortices at a given time step.

    :param solver: UnsteadyRingVortexLatticeMethodSolver This is the
        UnsteadyRingVortexLatticeMethodSolver with the wake RingVortices to process.
    :param step: int This is the time step number at which to process the wake
        RingVortices.
    :return: PolyData  This is the PolyData representation of the wake RingVortices.
    """
    num_wake_ring_vortices = solver.list_num_wake_vortices[step]
    stackFrwrvp_GP1_CgP1 = solver.listStackFrwrvp_GP1_CgP1[step]
    stackFlwrvp_GP1_CgP1 = solver.listStackFlwrvp_GP1_CgP1[step]
    stackBlwrvp_GP1_CgP1 = solver.listStackBlwrvp_GP1_CgP1[step]
    stackBrwrvp_GP1_CgP1 = solver.listStackBrwrvp_GP1_CgP1[step]

    # Initialize empty ndarrays to hold each wake RingVortex's vertices and face.
    wake_ring_vortex_vertices = np.zeros((0, 3), dtype=float)
    wake_ring_vortex_faces = np.zeros(0, dtype=int)

    for wake_ring_vortex_num in range(num_wake_ring_vortices):
        Frwrvp_GP1_CgP1 = stackFrwrvp_GP1_CgP1[wake_ring_vortex_num]
        Flwrvp_GP1_CgP1 = stackFlwrvp_GP1_CgP1[wake_ring_vortex_num]
        Blwrvp_GP1_CgP1 = stackBlwrvp_GP1_CgP1[wake_ring_vortex_num]
        Brwrvp_GP1_CgP1 = stackBrwrvp_GP1_CgP1[wake_ring_vortex_num]

        wake_ring_vortex_vertices_to_add = np.vstack(
            (
                Flwrvp_GP1_CgP1,
                Frwrvp_GP1_CgP1,
                Brwrvp_GP1_CgP1,
                Blwrvp_GP1_CgP1,
            )
        )
        wake_ring_vortex_face_to_add = np.array(
            [
                4,
                (wake_ring_vortex_num * 4),
                (wake_ring_vortex_num * 4) + 1,
                (wake_ring_vortex_num * 4) + 2,
                (wake_ring_vortex_num * 4) + 3,
            ],
            dtype=int,
        )

        # Stack this wake RingVortex's vertices and faces to the ndarrays of all wake
        # RingVortices' vertices and faces.
        wake_ring_vortex_vertices = cast(
            np.ndarray,
            np.vstack((wake_ring_vortex_vertices, wake_ring_vortex_vertices_to_add)),
        )
        wake_ring_vortex_faces = cast(
            np.ndarray,
            np.hstack((wake_ring_vortex_faces, wake_ring_vortex_face_to_add)),
        )

        # Increment the wake RingVortex counter.
        wake_ring_vortex_num += 1

    # Return the wake RingVortex surfaces.
    return pv.PolyData(wake_ring_vortex_vertices, wake_ring_vortex_faces)


# TEST: Consider adding unit tests for this function.
def _get_scalars(
    airplanes: list[geometry.airplane.Airplane],
    scalar_type: str,
    qInf__E: float,
) -> np.ndarray:
    """This function gets the coefficient values from a SteadyProblem's Airplanes'
    Wings' Panels.

    :param airplanes: list of Airplanes This is the list of Airplanes with the scalars
        to return.
    :param scalar_type: str This variable determines how which load coefficient to
        return as scalars. It can be "induced drag", "side force", or "lift", which
        respectively use each Panel's induced drag, side force, and lift coefficient.
    :param qInf__E: float This is the freestream dynamic pressure experienced by this
        SteadyProblem's Airplane(s) (observed in the Earth frame).
    :return scalars: (N,) ndarray of floats This is the (N,) ndarray of floats for the N
        Panels' load coefficients.
    """
    scalars = np.empty(0, dtype=float)

    # Iterate through the Airplanes' Wings.
    for airplane in airplanes:
        for wing in airplane.wings:
            _panels = wing.panels
            assert _panels is not None

            # Unravel this Wing's ndarray of Panels iterate through them.
            these_panels = np.ravel(_panels)
            for this_panel in these_panels:

                # Stack this Panel's scalars.
                if scalar_type == "induced drag":
                    this_induced_drag_coefficient = (
                        -this_panel.forces_W[0] / qInf__E / this_panel.area
                    )

                    scalars = cast(
                        np.ndarray, np.hstack((scalars, this_induced_drag_coefficient))
                    )

                if scalar_type == "side force":
                    this_side_force_coefficient = (
                        this_panel.forces_W[1] / qInf__E / this_panel.area
                    )

                    scalars = cast(
                        np.ndarray, np.hstack((scalars, this_side_force_coefficient))
                    )

                if scalar_type == "lift":
                    this_lift_coefficient = (
                        -this_panel.forces_W[2] / qInf__E / this_panel.area
                    )

                    scalars = cast(
                        np.ndarray, np.hstack((scalars, this_lift_coefficient))
                    )

    # Return the resulting ndarray of scalars.
    return scalars


# TEST: Consider adding unit tests for this function.
def _plot_scalars(
    plotter: pv.Plotter,
    these_scalars: np.ndarray,
    scalar_type: str,
    min_scalar: float,
    max_scalar: float,
    color_map: str,
    c_min: float,
    c_max: float,
    panel_surfaces: pv.PolyData,
) -> None:
    """This function plots a scalar bar, the surfaces of a set of Panels with a
    particular set of scalars, and labels for the minimum and maximum scalar values.

    :param plotter: pyvista.Plotter

        The Plotter used for visualization.

    :param these_scalars: (N,) ndarray of floats

        This is the ndarray of floats for each of the N Panels' coefficient values.

    :param scalar_type: str

        This variable determines how which load coefficient to return as scalars. It
        can be "induced drag", "side force", or "lift", which respectively use each
        Panel's induced drag, side force, and lift coefficient.

    :param min_scalar: float

        Minimum scalar value, which is displayed as text on the Plotter.

    :param max_scalar: float

        Maximum scalar value, which is displayed as text on the Plotter.

    :param color_map: str

        Name of the color map to use for scalar visualization.

    :param c_min: float

        Lower bound for the color map scaling.

    :param c_max: float

        Upper bound for the color map scaling.

    :param panel_surfaces: pyvista.PolyData

        PolyData representing the Panels' surfaces.

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
        cmap=color_map,  # type: ignore[arg-type]
        clim=[c_min, c_max],
        scalars=these_scalars,
        smooth_shading=False,
        scalar_bar_args=scalar_bar_args,  # type: ignore[arg-type]
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
