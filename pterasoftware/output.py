"""Contains functions for visualizing geometry and results.

**Contains the following classes:**

None

**Contains the following functions:**

draw: Draws a solver's Airplane(s).

animate: Animates the Airplane(s) of an UnsteadyRingVortexLatticeMethodSolver or one of
its subclasses (the aeroelastic or free flight solver).

plot_results_versus_time: Plots the loads and load coefficients of an
UnsteadyRingVortexLatticeMethodSolver or one of its subclasses (the aeroelastic or free
flight solver) as a function of time. For a free flight solver, it also plots the first
Airplane's six-degree-of-freedom state history.

log_results: Logs a solver's load and load coefficients, and, for a free flight solver,
the first Airplane's initial and final six-degree-of-freedom state.
"""

from __future__ import annotations

import math
import time

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import webp

from . import (
    _logging,
    _parameter_validation,
    _transformations,
    free_flight_unsteady_ring_vortex_lattice_method,
    geometry,
)
from . import operating_point as operating_point_mod
from . import (
    steady_horseshoe_vortex_lattice_method,
    steady_ring_vortex_lattice_method,
    unsteady_ring_vortex_lattice_method,
)
from .movements import free_flight_movement as free_flight_movement_mod

_logger = _logging.get_logger("output")

# Define the color and colormaps used by the visualization functions.
_sequential_color_map = "speed"
_diverging_color_map = "delta"
_wake_vortex_color = "white"
_panel_color = "chartreuse"
_streamline_color = "orchid"
_image_surface_opacity = 0.5
_image_surface_scale = 5.0
_image_reflection_mute_factor = 0.5
_image_surface_checker_size = 25
_image_surface_color_a = np.array([40, 40, 40], dtype=np.uint8)
_image_surface_color_b = np.array([80, 80, 80], dtype=np.uint8)
_plotter_background_color = "black"
_figure_background_color = "None"
_text_color = (129, 129, 129)
_text_color_normalized: tuple[float, float, float] = (
    _text_color[0] / 255,
    _text_color[1] / 255,
    _text_color[2] / 255,
)
_text_color_surface = (220, 220, 220)
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
    _alpha_color,
    _beta_color,
    _linear_x_color,
    _linear_y_color,
    _linear_z_color,
    _angular_x_color,
    _angular_y_color,
    _angular_z_color,
] = _prism[1:9]

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

# Define the camera's view-up direction for free flight visualizations. Earth axes have
# +z pointing down, so physical up is the -z direction. The free flight visualizations
# render geometry in Earth axes (so the body flies through the scene in its true pose)
# and use this view-up so that down appears downward on screen. This is a rendering
# setting, not an axis system.
_free_flight_view_up_E = np.array([0.0, 0.0, -1.0], dtype=float)

# Define the camera's view direction for free flight visualizations, given as the offset
# from the focal point to the camera position (in Earth axes). This views the scene
# obliquely from the South, West, and above (Earth -x, -y, and -z).
_free_flight_view_direction_E = np.array([1.0, -1.0, -1.0], dtype=float)
_free_flight_view_direction_E = _free_flight_view_direction_E / np.linalg.norm(
    _free_flight_view_direction_E
)


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
    """Draws a solver's Airplane(s).

    For a FreeFlightUnsteadyRingVortexLatticeMethodSolver, the final time step's
    geometry is rendered in Earth axes at its true position and orientation, so the
    body's flight attitude is visible, rather than being drawn body-fixed in geometry
    axes.

    **Citation:**

    Adapted from: vlm3.draw in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 03/28/2020

    :param solver: The solver whose Airplane(s) will be plotted. This can be a
        SteadyHorseshoeVortexLatticeMethodSolver, a SteadyRingVortexLatticeMethodSolver,
        or an UnsteadyRingVortexLatticeMethodSolver. The
        UnsteadyRingVortexLatticeMethodSolver's subclasses, the
        AeroelasticUnsteadyRingVortexLatticeMethodSolver and the
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, are also accepted.
    :param scalar_type: Determines how to color the Panels. Setting this to None colors
        the Panels uniformly. If the solver has been run, it can also be "induced drag",
        "side force", or "lift", which respectively use each Panel's induced drag, side
        force, and lift coefficient. The default is None.
    :param show_streamlines: Set this to True to show the streamlines emanating from the
        back of the Wings. If True, the solver's streamlines must have already been
        calculated. Can be a bool or a numpy bool and will be converted internally to a
        bool. The default is False.
    :param show_wake_vortices: Set this to True to show any wake ring vortices. If True,
        the solver must be an UnsteadyRingVortexLatticeMethodSolver and must have
        already been run. Can be a bool or a numpy bool and will be converted internally
        to a bool. The default is False.
    :param save: Set this to True to save the image as a WebP. It can be a bool or a
        numpy bool and will be converted internally to a bool. The default is False.
    :param testing: Set this to True to close the image after one second, which is
        useful for running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
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

    # For a free flight solver, geometry is rendered in its true Earth-frame pose so the
    # body flies through the scene. T_pas_GP1_CgP1_to_E_Eo holds the passive
    # transformation from the first Airplane's geometry axes (relative to its CG) to Earth
    # axes (relative to the Earth origin) for the drawn time step, and stays None for the
    # standard body-fixed rendering in geometry axes.
    is_free_flight = isinstance(
        solver,
        free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    )
    T_pas_GP1_CgP1_to_E_Eo: np.ndarray | None = None

    # Get the solver's geometry and OperatingPoint.
    if isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        draw_step = solver.num_steps - 1

        airplanes = solver.steady_problems[draw_step].airplanes
        draw_operating_point = solver.steady_problems[draw_step].operating_point
        qInf__E = draw_operating_point.qInf__E

        if is_free_flight:
            T_pas_GP1_CgP1_to_E_Eo = _get_T_pas_GP1_CgP1_to_E_Eo(draw_operating_point)

        # If showing wake ring vortices, get their surfaces and plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces(
                solver, draw_step
            )
            if T_pas_GP1_CgP1_to_E_Eo is not None:
                wake_ring_vortex_surfaces = _transform_mesh(
                    wake_ring_vortex_surfaces, T_pas_GP1_CgP1_to_E_Eo
                )
            plotter.add_mesh(
                wake_ring_vortex_surfaces,
                show_edges=True,
                smooth_shading=False,
                color=_wake_vortex_color,
            )
    else:
        airplanes = solver.airplanes
        draw_operating_point = solver.operating_point
        qInf__E = draw_operating_point.qInf__E

    # Get the Panel surfaces, mapping them into Earth axes for free flight.
    panel_surfaces = _get_panel_surfaces(airplanes)
    if T_pas_GP1_CgP1_to_E_Eo is not None:
        panel_surfaces = _transform_mesh(panel_surfaces, T_pas_GP1_CgP1_to_E_Eo)

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

        T_reflect = draw_operating_point.surfaceReflect_T_act_GP1_CgP1
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
            text_color=_text_color_surface if T_reflect is not None else _text_color,
        )
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=_panel_color,
            smooth_shading=False,
        )
        T_reflect = draw_operating_point.surfaceReflect_T_act_GP1_CgP1
    image_surface_mesh = None

    # For free flight, the active reflection is represented in geometry axes, but the
    # geometry has been mapped into Earth axes. Re-expressing the reflection in Earth axes
    # (a change of basis by the same passive transformation) lets the reflected-geometry
    # code below operate entirely in Earth axes.
    if T_pas_GP1_CgP1_to_E_Eo is not None and T_reflect is not None:
        T_reflect = (
            T_pas_GP1_CgP1_to_E_Eo
            @ T_reflect
            @ _transformations.invert_T_pas(T_pas_GP1_CgP1_to_E_Eo)
        )

    # If an image surface is defined, add reflected geometry. The image surface plane
    # is added later, after the geometry bounds are captured.
    if T_reflect is not None:
        mute = _image_reflection_mute_factor
        muted_edge_color = _mute_color("black", mute)

        # Add reflected Panel surfaces with muted coloring.
        reflected_panel_surfaces = _reflect_mesh(panel_surfaces, T_reflect)
        if scalar_type in ("induced drag", "side force", "lift"):
            plotter.add_mesh(
                reflected_panel_surfaces,
                show_edges=True,
                edge_color=muted_edge_color,
                cmap=_mute_colormap(color_map, mute),
                clim=[c_min, c_max],
                scalars=these_scalars,
                smooth_shading=False,
                show_scalar_bar=False,
            )
        else:
            plotter.add_mesh(
                reflected_panel_surfaces,
                show_edges=True,
                edge_color=muted_edge_color,
                color=_mute_color(_panel_color, mute),
                smooth_shading=False,
            )

        # Add reflected wake ring vortex surfaces if they are being shown.
        if show_wake_vortices:
            plotter.add_mesh(
                _reflect_mesh(wake_ring_vortex_surfaces, T_reflect),
                show_edges=True,
                edge_color=muted_edge_color,
                smooth_shading=False,
                color=_mute_color(_wake_vortex_color, mute),
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

                    # For free flight, map the segment into Earth axes.
                    if T_pas_GP1_CgP1_to_E_Eo is not None:
                        point = _transformations.apply_T_to_vectors(
                            T_pas_GP1_CgP1_to_E_Eo, point, is_position=True
                        )
                        last_point = _transformations.apply_T_to_vectors(
                            T_pas_GP1_CgP1_to_E_Eo, last_point, is_position=True
                        )

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

                    # If an image surface is defined, add the reflected streamline
                    # segment.
                    if T_reflect is not None:
                        reflected_point = _transformations.apply_T_to_vectors(
                            T_reflect,
                            point,
                            is_position=True,
                        )
                        reflected_last_point = _transformations.apply_T_to_vectors(
                            T_reflect,
                            last_point,
                            is_position=True,
                        )
                        plotter.add_mesh(
                            pv.Line(
                                reflected_last_point,
                                reflected_point,
                            ),
                            show_edges=True,
                            color=_mute_color(_streamline_color, mute),
                            line_width=2,
                            smooth_shading=False,
                        )

    # If an image surface is defined, save the geometry bounds (which now include
    # the reflected geometry but not the image surface plane), add the image surface
    # plane, then fit the camera to the saved bounds so the view is not dominated by
    # the much larger image surface plane. When an image surface is present, cpos is
    # not passed to show() because that would trigger an auto-fit to all actors
    # (including the image surface).
    if T_reflect is not None:
        geometry_bounds = plotter.bounds
        if T_pas_GP1_CgP1_to_E_Eo is not None:
            # The image surface helper builds the plane from geometry-axis quantities, so
            # it needs geometry-axis bounds. Build the plane there, then map it into Earth
            # axes to match the rendered geometry.
            geometry_axis_bounds = _get_panel_surfaces(airplanes).bounds
            image_surface_result = _get_image_surface_mesh_and_texture(
                draw_operating_point, geometry_axis_bounds
            )
            assert image_surface_result is not None
            image_surface_mesh, image_surface_texture = image_surface_result
            image_surface_mesh = _transform_mesh(
                image_surface_mesh, T_pas_GP1_CgP1_to_E_Eo
            )
        else:
            image_surface_result = _get_image_surface_mesh_and_texture(
                draw_operating_point, geometry_bounds
            )
            assert image_surface_result is not None
            image_surface_mesh, image_surface_texture = image_surface_result
        plotter.add_mesh(
            image_surface_mesh,
            texture=image_surface_texture,
            opacity=_image_surface_opacity,
            smooth_shading=True,
        )

        # For the standard body-fixed rendering, fit the camera to the geometry bounds
        # so the much larger image surface plane does not dominate the view. Free flight
        # uses its own Earth-axes camera, computed below.
        if T_pas_GP1_CgP1_to_E_Eo is None:
            plotter.camera.position = (-1, -1, 1)
            plotter.camera.focal_point = (0, 0, 0)
            plotter.camera.up = (0, 0, 1)
            plotter.reset_camera(bounds=geometry_bounds)  # type: ignore[call-arg]

    # Choose the camera position. Free flight frames the body in Earth axes with physical
    # up as Earth -z. The standard rendering views geometry axes from (-1, -1, 1), unless
    # an image surface is present, in which case the camera was already fitted above and
    # cpos is left None so show() does not auto-fit to the large image surface plane.
    draw_cpos: tuple | None
    if is_free_flight:
        # Aim the camera along the Earth-axes view direction with physical up, then fit
        # the parallel projection to the rendered geometry. The camera is set explicitly
        # here (rather than via cpos) so that show() preserves the fit, since passing a
        # full position to cpos would not set the parallel projection's scale.
        airplane_bounds = np.array(panel_surfaces.bounds, dtype=float)
        center_E_Eo = 0.5 * (airplane_bounds[1::2] + airplane_bounds[::2])
        airplane_diagonal = float(
            np.linalg.norm(airplane_bounds[1::2] - airplane_bounds[::2])
        )
        plotter.camera.focal_point = tuple(center_E_Eo)
        plotter.camera.position = tuple(
            center_E_Eo + 3.0 * airplane_diagonal * _free_flight_view_direction_E
        )
        plotter.camera.up = _free_flight_view_up_E
        plotter.reset_camera()  # type: ignore[call-arg]
        draw_cpos = None
    elif image_surface_mesh is None:
        draw_cpos = (-1, -1, 1)
    else:
        draw_cpos = None

    # Set the Plotter's background color.
    plotter.set_background(color=_plotter_background_color)  # type: ignore[call-arg]
    if not testing:
        # Show the Plotter so the user can adjust the camera position and window.
        # When the user closes the window, the Plotter still exists. Therefore,
        # it can later be saved as an image if desired.
        plotter.show(
            title="Orient the view, then press any key to continue.",
            cpos=draw_cpos,
            full_screen=False,
            auto_close=False,
        )
    else:
        # Show the Plotter for 1 second, then proceed automatically. This is useful
        # for testing.
        plotter.show(
            cpos=draw_cpos,
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
            img=image, file_path="draw.webp", lossless=False, quality=_quality
        )

    # Close all Plotters.
    pv.close_all()


def animate(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    scalar_type: str | None = None,
    show_wake_vortices: bool | np.bool_ = False,
    save: bool | np.bool_ = False,
    testing: bool | np.bool_ = False,
) -> None:
    """Animates the Airplane(s) of an UnsteadyRingVortexLatticeMethodSolver or one of
    its subclasses (the aeroelastic or free flight solver).

    For a FreeFlightUnsteadyRingVortexLatticeMethodSolver, each time step's geometry is
    rendered in Earth axes at its true position and orientation, so the body flies
    through the scene along its trajectory rather than staying fixed while the wake
    streams behind it. The camera frames the whole trajectory.

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver whose Airplane(s)
        will be animated. Its subclasses, the
        AeroelasticUnsteadyRingVortexLatticeMethodSolver and the
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, are also accepted.
    :param scalar_type: Determines how to color the Panels. Setting this to None colors
        the Panels uniformly. If the solver has been run, it can also be "induced drag",
        "side force", or "lift", which respectively use each Panel's induced drag, side
        force, and lift coefficient. The default is None.
    :param show_wake_vortices: Set this to True to show any wake ring vortices. If True,
        the solver must have already been run. Can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :param save: Set this to True to save the animation as an animated WebP. It can be a
        bool or a numpy bool and will be converted internally to a bool. The default is
        False.
    :param testing: Set this to True to start the animation after one second, which is
        useful for running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :return: None
    """
    if not isinstance(
        unsteady_solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
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

    # For a free flight solver, each time step's geometry is rendered in its true Earth-
    # frame pose so the body flies through the scene. step_transforms holds, per time
    # step, the passive transformation from the first Airplane's geometry axes (relative
    # to its CG) to Earth axes (relative to the Earth origin). It stays empty for the
    # standard body-fixed rendering in geometry axes.
    is_free_flight = isinstance(
        unsteady_solver,
        free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    )
    step_transforms: list[np.ndarray] = []
    if is_free_flight:
        step_transforms = [
            _get_T_pas_GP1_CgP1_to_E_Eo(steady_problem.operating_point)
            for steady_problem in unsteady_solver.steady_problems
        ]

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
            all_scalars = np.hstack((all_scalars, scalars_to_add))

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

    # Pre-compute the image surface mesh and reflection matrix from the last time
    # step's geometry so that the plane is large enough to encompass the fully
    # developed wake. The mesh, texture, and reflection matrix are static and reused
    # for every frame. The last step's geometry bounds (including reflected geometry
    # but not the image surface plane) are also saved so the camera can be fitted to
    # the geometry rather than the larger image surface.
    last_step = len(step_airplanes) - 1
    last_step_operating_point = unsteady_solver.steady_problems[
        last_step
    ].operating_point
    T_reflect = last_step_operating_point.surfaceReflect_T_act_GP1_CgP1
    animate_text_color = _text_color_surface if T_reflect is not None else _text_color
    image_surface_geometry_bounds = None
    if is_free_flight:
        if T_reflect is not None:
            # The image surface is fixed in the world, so build its plane once from the
            # last step's geometry-axis quantities, then map it into Earth axes. Re-
            # expressing the reflection in Earth axes (a change of basis by the same
            # passive transformation) lets the reflected geometry be built from the
            # Earth-axes panels each frame.
            T_pas_last = step_transforms[last_step]
            geometry_axis_bounds = _get_panel_surfaces(step_airplanes[last_step]).bounds
            image_surface_result = _get_image_surface_mesh_and_texture(
                last_step_operating_point, geometry_axis_bounds
            )
            assert image_surface_result is not None
            image_surface_mesh, image_surface_texture = image_surface_result
            image_surface_mesh = _transform_mesh(image_surface_mesh, T_pas_last)
            T_reflect = (
                T_pas_last @ T_reflect @ _transformations.invert_T_pas(T_pas_last)
            )
        else:
            image_surface_mesh = None
            image_surface_texture = None
    elif T_reflect is not None:
        last_step_panel_surfaces = _get_panel_surfaces(step_airplanes[last_step])
        reflected_last_step_panel_surfaces = _reflect_mesh(
            last_step_panel_surfaces, T_reflect
        )
        if show_wake_vortices:
            last_step_wake_surfaces = _get_wake_ring_vortex_surfaces(
                unsteady_solver, last_step
            )
            reflected_last_step_wake_surfaces = _reflect_mesh(
                last_step_wake_surfaces, T_reflect
            )
            combined = (
                last_step_panel_surfaces.merge(last_step_wake_surfaces)
                .merge(reflected_last_step_panel_surfaces)
                .merge(reflected_last_step_wake_surfaces)
            )
            image_surface_geometry_bounds = combined.bounds
        else:
            combined = last_step_panel_surfaces.merge(
                reflected_last_step_panel_surfaces
            )
            image_surface_geometry_bounds = combined.bounds
        image_surface_result = _get_image_surface_mesh_and_texture(
            last_step_operating_point, image_surface_geometry_bounds
        )
        assert image_surface_result is not None
        image_surface_mesh, image_surface_texture = image_surface_result
    else:
        image_surface_mesh = None
        image_surface_texture = None

    # For free flight, compute a fixed camera that frames the whole trajectory. The body
    # moves through the scene, so the camera is centered on the trajectory's midpoint and
    # the parallel scale is fit to the projected extent of the trajectory's geometry so
    # the whole glide stays in view. The clipping range is sized later, after the user has
    # oriented the view, so it tracks the user's chosen camera rather than the default.
    free_flight_cpos: list | None = None
    free_flight_parallel_scale = 0.0
    free_flight_clip_meshes: list = []
    if is_free_flight:
        initialPosition_E_Eo = unsteady_solver.steady_problems[
            0
        ].operating_point.CgP1_E_Eo
        finalPosition_E_Eo = last_step_operating_point.CgP1_E_Eo
        trajectoryMidpoint_E_Eo = 0.5 * (initialPosition_E_Eo + finalPosition_E_Eo)
        trajectory_extent = float(
            np.linalg.norm(finalPosition_E_Eo - initialPosition_E_Eo)
        )

        # Map the first and last frames' Panel surfaces into Earth axes. These two frames
        # bound the trajectory, so their combined extent frames the whole glide.
        first_step_panel_surfaces = _transform_mesh(
            _get_panel_surfaces(step_airplanes[0]), step_transforms[0]
        )
        last_step_panel_surfaces = _transform_mesh(
            _get_panel_surfaces(step_airplanes[last_step]),
            step_transforms[last_step],
        )
        airplane_bounds = np.array(first_step_panel_surfaces.bounds, dtype=float)
        airplane_diagonal = float(
            np.linalg.norm(airplane_bounds[1::2] - airplane_bounds[::2])
        )

        # Aim the camera along the Earth-axes view direction with physical up, centered on
        # the trajectory's midpoint and far enough back to clear the geometry at both ends.
        padding = max(2.0 * airplane_diagonal, 0.5 * trajectory_extent)
        camera_distance = trajectory_extent + padding
        cameraPosition_E_Eo = (
            trajectoryMidpoint_E_Eo + camera_distance * _free_flight_view_direction_E
        )
        free_flight_cpos = [
            tuple(cameraPosition_E_Eo),
            tuple(trajectoryMidpoint_E_Eo),
            _free_flight_view_up_E,
        ]

        # Collect the geometry that frames the trajectory: the body at both ends, plus the
        # last frame's wake (the largest) if it is shown. The last frame's surfaces are
        # reused below to size the clipping range after the user orients the view.
        framing_meshes = [first_step_panel_surfaces, last_step_panel_surfaces]
        free_flight_clip_meshes = [last_step_panel_surfaces]
        if show_wake_vortices:
            last_step_wake_surfaces = _transform_mesh(
                _get_wake_ring_vortex_surfaces(unsteady_solver, last_step),
                step_transforms[last_step],
            )
            if last_step_wake_surfaces.n_points > 0:
                framing_meshes.append(last_step_wake_surfaces)
                free_flight_clip_meshes.append(last_step_wake_surfaces)

        # Fit the parallel scale (half the viewport height in world units, since the
        # projection is parallel) to the projected extent of that geometry about the focal
        # point. This frames the glide snugly; the user can rescale interactively before
        # the animation is captured.
        free_flight_parallel_scale = _free_flight_fit_parallel_scale(
            framing_meshes,
            trajectoryMidpoint_E_Eo,
            _free_flight_view_direction_E,
            _free_flight_view_up_E,
        )

    # If saving the animation, add text that displays its speed.
    if save:
        plotter.add_text(
            text="Speed: " + str(round(100 * speed)) + "%",
            position=_text_speed_position,
            font_size=_text_font_size,
            viewport=True,
            color=animate_text_color,
        )

    # Get the Panel surfaces of the first time step's Airplane(s), mapping them into
    # Earth axes for free flight.
    panel_surfaces = _get_panel_surfaces(step_airplanes[0])
    if is_free_flight:
        panel_surfaces = _transform_mesh(panel_surfaces, step_transforms[0])

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
            text_color=_text_color_surface if T_reflect is not None else _text_color,
        )
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=_panel_color,
            smooth_shading=False,
        )

    # If an image surface is defined, add reflected geometry, plot the pre-computed
    # plane, set the camera direction, and fit the camera to the last time step's
    # geometry bounds so the view is not dominated by the much larger image surface
    # plane. When an image surface is present, cpos is not passed to show() because
    # that would trigger an auto-fit to all actors (including the image surface).
    if T_reflect is not None:
        assert image_surface_mesh is not None
        mute = _image_reflection_mute_factor
        muted_edge_color = _mute_color("black", mute)
        muted_panel_color = _mute_color(_panel_color, mute)
        muted_wake_color = _mute_color(_wake_vortex_color, mute)
        if color_map:
            muted_color_map = _mute_colormap(color_map, mute)
        else:
            muted_color_map = None

        # Add reflected Panel surfaces with muted coloring.
        reflected_panel_surfaces = _reflect_mesh(panel_surfaces, T_reflect)
        if scalar_type is not None and first_results_step == 0:
            plotter.add_mesh(
                reflected_panel_surfaces,
                show_edges=True,
                edge_color=muted_edge_color,
                cmap=muted_color_map,
                clim=[c_min, c_max],
                scalars=these_scalars,
                smooth_shading=False,
                show_scalar_bar=False,
            )
        else:
            plotter.add_mesh(
                reflected_panel_surfaces,
                show_edges=True,
                edge_color=muted_edge_color,
                color=muted_panel_color,
                smooth_shading=False,
            )

        # Add the image surface plane.
        plotter.add_mesh(
            image_surface_mesh,
            texture=image_surface_texture,
            opacity=_image_surface_opacity,
            smooth_shading=True,
        )

        # For the standard body-fixed rendering, fit the camera to the geometry bounds
        # so the much larger image surface plane does not dominate the view. Free flight
        # uses its own trajectory-framing camera, applied below.
        if not is_free_flight:
            plotter.camera.position = (-1, -1, 1)
            plotter.camera.focal_point = (0, 0, 0)
            plotter.camera.up = (0, 0, 1)
            plotter.reset_camera(bounds=image_surface_geometry_bounds)  # type: ignore[call-arg]

    # Choose the camera position. Free flight sets its trajectory-framing camera
    # explicitly (rather than through show()'s cpos) so it survives the mesh additions
    # above and becomes the default for the held first frame; cpos is then left None so
    # show() does not reset it, which lets the user's interactive reorientation and
    # rescaling carry through to the animation. The standard rendering views geometry axes
    # from (-1, -1, 1), unless an image surface is present, in which case the camera was
    # already fitted above and cpos is left None so show() does not auto-fit to the large
    # image surface plane.
    animate_cpos: tuple | list | None
    if is_free_flight:
        assert free_flight_cpos is not None
        plotter.camera_position = free_flight_cpos
        plotter.camera.parallel_scale = free_flight_parallel_scale
        animate_cpos = None
    elif image_surface_mesh is None:
        animate_cpos = (-1, -1, 1)
    else:
        animate_cpos = None

    # Set the Plotter's background color.
    plotter.set_background(color=_plotter_background_color)  # type: ignore[call-arg]

    # If not testing, show the Plotter with the first time step so the user can
    # orient the view. When the user presses any key, set the title back to the
    # animation title and proceed. If testing, show the Plotter with the first time
    # step for 1 second, and start the animation with the current window view.
    if not testing:
        plotter.show(
            title="Orient the view, then press any key to produce the animation.",
            cpos=animate_cpos,
            full_screen=False,
            auto_close=False,
        )
        assert plotter.ren_win is not None
        plotter.ren_win.SetWindowName("Rendering speed not to scale.")
    else:
        plotter.show(
            title="Rendering speed not to scale.",
            cpos=animate_cpos,
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # The user may have reoriented or rescaled the view during the held first frame.
    # Preserve that camera and only size the clipping range so every frame stays visible:
    # temporarily add the last frame's geometry (the first frame's is already present), fit
    # the clipping range to both, then remove the temporary actors. The body moves through
    # the scene, so a clipping range fit to the first frame alone would clip later frames.
    if is_free_flight:
        temporary_actors = [
            plotter.add_mesh(clip_mesh) for clip_mesh in free_flight_clip_meshes
        ]
        plotter.reset_camera_clipping_range()
        free_flight_clipping_range = plotter.camera.clipping_range
        for temporary_actor in temporary_actors:
            plotter.remove_actor(temporary_actor)
        plotter.camera.clipping_range = free_flight_clipping_range

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

        # Get the Panel surfaces of this time step's Airplane(s), mapping them into Earth
        # axes for free flight.
        panel_surfaces = _get_panel_surfaces(airplanes)
        if is_free_flight:
            panel_surfaces = _transform_mesh(
                panel_surfaces, step_transforms[current_step]
            )

        # If saving the animation, add text that displays its speed.
        if save:
            plotter.add_text(
                text="Speed: " + str(round(100 * speed)) + "%",
                position=_text_speed_position,
                font_size=_text_font_size,
                viewport=True,
                color=animate_text_color,
            )

        # If showing wake ring vortices, get their surfaces and plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces(
                unsteady_solver, current_step
            )
            if is_free_flight:
                wake_ring_vortex_surfaces = _transform_mesh(
                    wake_ring_vortex_surfaces, step_transforms[current_step]
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
                text_color=(
                    _text_color_surface if T_reflect is not None else _text_color
                ),
            )
        else:
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                color=_panel_color,
                smooth_shading=False,
            )

        # If an image surface is defined, add reflected geometry and the pre-computed
        # image surface plane.
        if T_reflect is not None:
            assert image_surface_mesh is not None

            # Add reflected Panel surfaces with muted coloring.
            reflected_panel_surfaces = _reflect_mesh(panel_surfaces, T_reflect)
            if scalar_type is not None and first_results_step <= current_step:
                plotter.add_mesh(
                    reflected_panel_surfaces,
                    show_edges=True,
                    edge_color=muted_edge_color,
                    cmap=muted_color_map,
                    clim=[c_min, c_max],
                    scalars=these_scalars,
                    smooth_shading=False,
                    show_scalar_bar=False,
                )
            else:
                plotter.add_mesh(
                    reflected_panel_surfaces,
                    show_edges=True,
                    edge_color=muted_edge_color,
                    color=muted_panel_color,
                    smooth_shading=False,
                )

            # Add reflected wake ring vortex surfaces if they are being shown.
            if show_wake_vortices:
                plotter.add_mesh(
                    _reflect_mesh(wake_ring_vortex_surfaces, T_reflect),
                    show_edges=True,
                    edge_color=muted_edge_color,
                    smooth_shading=False,
                    color=muted_wake_color,
                )

            # Add the image surface plane.
            plotter.add_mesh(
                image_surface_mesh,
                texture=image_surface_texture,
                opacity=_image_surface_opacity,
                smooth_shading=True,
            )

        # If an image surface is present, force VTK to recalculate the scalar bar
        # layout. Adding the image surface mesh with opacity causes VTK's
        # UnconstrainedFontSize layout to misposition the left label (PyVista
        # issue #7516).
        if T_reflect is not None:
            for scalar_bar_actor in plotter.scalar_bars.values():
                scalar_bar_actor.Modified()
            plotter.render()

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
            images, "animate.webp", fps=actual_fps, lossless=False, quality=_quality
        )

    # Close all the Plotters.
    pv.close_all()


def plot_results_versus_time(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    show: bool | np.bool_ = True,
    save: bool | np.bool_ = False,
) -> None:
    """Plots the loads and load coefficients of an UnsteadyRingVortexLatticeMethodSolver
    or one of its subclasses (the aeroelastic or free flight solver) as a function of
    time.

    For a FreeFlightUnsteadyRingVortexLatticeMethodSolver, this also plots the first
    Airplane's six-degree-of-freedom state history: its position, velocity, orientation,
    angular velocity, and aerodynamic angles versus time. These describe the first
    Airplane, the rigid body the dynamics integrate, so they are plotted once for the
    whole simulation rather than per Airplane.

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver whose loads and
        load coefficients will be plotted. Its subclasses, the
        AeroelasticUnsteadyRingVortexLatticeMethodSolver and the
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, are also accepted. For a
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, the first Airplane's state
        history is plotted as well.
    :param show: Set this to True to show the plots. It can be a bool or a numpy bool
        and will be converted internally to a bool. The default is True.
    :param save: Set this to True to save the plots as PNGs. It can be a bool or a numpy
        bool and will be converted internally to a bool. The default is False.
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
        force_axes.spines.bottom.set_color(_text_color_normalized)
        force_axes.spines.left.set_color(_text_color_normalized)
        force_axes.xaxis.label.set_color(_text_color_normalized)
        force_axes.yaxis.label.set_color(_text_color_normalized)
        force_coefficients_axes.spines.bottom.set_color(_text_color_normalized)
        force_coefficients_axes.spines.left.set_color(_text_color_normalized)
        force_coefficients_axes.xaxis.label.set_color(_text_color_normalized)
        force_coefficients_axes.yaxis.label.set_color(_text_color_normalized)
        moment_coefficients_axes.spines.bottom.set_color(_text_color_normalized)
        moment_coefficients_axes.spines.left.set_color(_text_color_normalized)
        moment_coefficients_axes.xaxis.label.set_color(_text_color_normalized)
        moment_coefficients_axes.yaxis.label.set_color(_text_color_normalized)
        moment_axes.spines.bottom.set_color(_text_color_normalized)
        moment_axes.spines.left.set_color(_text_color_normalized)
        moment_axes.xaxis.label.set_color(_text_color_normalized)
        moment_axes.yaxis.label.set_color(_text_color_normalized)

        # Format all the plots' tick colors.
        force_axes.tick_params(axis="x", colors=_text_color_normalized)
        force_axes.tick_params(axis="y", colors=_text_color_normalized)
        force_coefficients_axes.tick_params(axis="x", colors=_text_color_normalized)
        force_coefficients_axes.tick_params(axis="y", colors=_text_color_normalized)
        moment_coefficients_axes.tick_params(axis="x", colors=_text_color_normalized)
        moment_coefficients_axes.tick_params(axis="y", colors=_text_color_normalized)
        moment_axes.tick_params(axis="x", colors=_text_color_normalized)
        moment_axes.tick_params(axis="y", colors=_text_color_normalized)

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
            color=_linear_x_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_axes.plot(
            times,
            forces_W[airplane_id, 1],
            label="Side Force",
            color=_linear_y_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_axes.plot(
            times,
            -forces_W[airplane_id, 2],
            label="Lift",
            color=_linear_z_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            -forceCoefficients_W[airplane_id, 0],
            label="Induced Drag Coefficient",
            color=_linear_x_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            forceCoefficients_W[airplane_id, 1],
            label="Side Force Coefficient",
            color=_linear_y_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        force_coefficients_axes.plot(
            times,
            -forceCoefficients_W[airplane_id, 2],
            label="Lift Coefficient",
            color=_linear_z_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_CgP1[airplane_id, 0],
            label="Rolling Moment",
            color=_angular_x_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_CgP1[airplane_id, 1],
            label="Pitching Moment",
            color=_angular_y_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_axes.plot(
            times,
            moments_W_CgP1[airplane_id, 2],
            label="Yawing Moment",
            color=_angular_z_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_CgP1[airplane_id, 0],
            label="Rolling Moment Coefficient",
            color=_angular_x_color,
            marker=".",
            markevery=(_marker_spacing * 0 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_CgP1[airplane_id, 1],
            label="Pitching Moment Coefficient",
            color=_angular_y_color,
            marker=".",
            markevery=(_marker_spacing * 1 / 3, _marker_spacing),
            markersize=_marker_size,
        )
        moment_coefficients_axes.plot(
            times,
            momentCoefficients_W_CgP1[airplane_id, 2],
            label="Yawing Moment Coefficient",
            color=_angular_z_color,
            marker=".",
            markevery=(_marker_spacing * 2 / 3, _marker_spacing),
            markersize=_marker_size,
        )

        # Find and format this Airplane's name for use in the plot titles.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[airplane_id].name
        airplane_name_snake = airplane_name.lower().replace(" ", "_")
        force_title = airplane_name + " Forces vs. Time"
        force_coefficient_title = airplane_name + " Force Coefficients vs. Time"
        moment_title = airplane_name + " Moments vs. Time"
        moment_coefficient_title = airplane_name + " Moment Coefficients vs. Time"

        # Name the plots' axis labels and titles.
        force_axes.set_xlabel("Time (s)", color=_text_color_normalized)
        force_axes.set_ylabel("Force (N)", color=_text_color_normalized)
        force_axes.set_title(force_title, color=_text_color_normalized)
        force_coefficients_axes.set_xlabel("Time (s)", color=_text_color_normalized)
        force_coefficients_axes.set_ylabel(
            "Force Coefficient", color=_text_color_normalized
        )
        force_coefficients_axes.set_title(
            force_coefficient_title, color=_text_color_normalized
        )
        moment_axes.set_xlabel("Time (s)", color=_text_color_normalized)
        moment_axes.set_ylabel("Moment (N m)", color=_text_color_normalized)
        moment_axes.set_title(moment_title, color=_text_color_normalized)
        moment_coefficients_axes.set_xlabel("Time (s)", color=_text_color_normalized)
        moment_coefficients_axes.set_ylabel(
            "Moment Coefficient", color=_text_color_normalized
        )
        moment_coefficients_axes.set_title(
            moment_coefficient_title, color=_text_color_normalized
        )

        # Format the plots' legends.
        force_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color_normalized,
        )
        force_coefficients_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color_normalized,
        )
        moment_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color_normalized,
        )
        moment_coefficients_axes.legend(
            facecolor=_figure_background_color,
            edgecolor=_figure_background_color,
            labelcolor=_text_color_normalized,
        )

        # Save the figures as PNGs if the user wants to do so.
        if save:
            force_figure.savefig(
                airplane_name_snake + "_forces.png",
                dpi=300,
            )
            force_coefficients_figure.savefig(
                airplane_name_snake + "_force_coefficients.png",
                dpi=300,
            )
            moment_figure.savefig(
                airplane_name_snake + "_moments.png",
                dpi=300,
            )
            moment_coefficients_figure.savefig(
                airplane_name_snake + "_moment_coefficients.png",
                dpi=300,
            )

    # For a free flight solver, also plot the first Airplane's six-degree-of-freedom
    # state history. This is plotted once for the whole simulation, since the state
    # describes the first Airplane, the single rigid body the dynamics integrate.
    if isinstance(
        unsteady_solver,
        free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    ):
        # Narrow the movement to a FreeFlightMovement so its operating point history is
        # typed. The solver type guarantees this, so the assert documents the invariant.
        movement = unsteady_solver.unsteady_problem.movement
        assert isinstance(movement, free_flight_movement_mod.FreeFlightMovement)
        operating_points = movement.operating_point_movement.operating_points

        # The state history covers every time step, starting at time step 0, so it uses
        # its own time array rather than the results-averaging window above.
        num_state_steps = len(operating_points)
        state_times = np.linspace(
            0.0,
            delta_time * (num_state_steps - 1),
            num_state_steps,
            endpoint=True,
        )

        # Initialize matrices to hold the state quantities at every time step.
        positions_E_Eo = np.zeros((3, num_state_steps), dtype=float)
        velocities_E__E = np.zeros((3, num_state_steps), dtype=float)
        anglesDeg_E_to_BP1_izyx = np.zeros((3, num_state_steps), dtype=float)
        omegasDeg_BP1__E = np.zeros((3, num_state_steps), dtype=float)
        alphas = np.zeros(num_state_steps, dtype=float)
        betas = np.zeros(num_state_steps, dtype=float)

        # Iterate through the time steps and extract each step's state.
        for step, this_operating_point in enumerate(operating_points):
            positions_E_Eo[:, step] = this_operating_point.CgP1_E_Eo
            velocities_E__E[:, step] = _velocity_E__E_from_operating_point(
                this_operating_point
            )
            anglesDeg_E_to_BP1_izyx[:, step] = this_operating_point.angles_E_to_BP1_izyx
            omegasDeg_BP1__E[:, step] = this_operating_point.omegas_BP1__E
            alphas[step] = this_operating_point.alpha
            betas[step] = this_operating_point.beta

        # The state describes the first Airplane (the rigid body MuJoCo integrates), so
        # the plot titles and file names use the first Airplane's name.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[0].name
        airplane_name_snake = airplane_name.lower().replace(" ", "_")

        _plot_state_history(
            state_times,
            [positions_E_Eo[0], positions_E_Eo[1], positions_E_Eo[2]],
            ["X Component", "Y Component", "Z Component"],
            [_linear_x_color, _linear_y_color, _linear_z_color],
            airplane_name
            + " Position (in Earth Axes, Relative to the Earth Origin) vs. Time",
            "Position (m)",
            save,
            airplane_name_snake + "_position.png",
        )
        _plot_state_history(
            state_times,
            [velocities_E__E[0], velocities_E__E[1], velocities_E__E[2]],
            ["X Component", "Y Component", "Z Component"],
            [_linear_x_color, _linear_y_color, _linear_z_color],
            airplane_name
            + " Velocity (in Earth Axes, Observed from the Earth Frame) vs. Time",
            "Velocity (m/s)",
            save,
            airplane_name_snake + "_velocity.png",
        )
        _plot_state_history(
            state_times,
            [
                anglesDeg_E_to_BP1_izyx[0],
                anglesDeg_E_to_BP1_izyx[1],
                anglesDeg_E_to_BP1_izyx[2],
            ],
            ["Roll Angle", "Pitch Angle", "Yaw Angle"],
            [_angular_x_color, _angular_y_color, _angular_z_color],
            airplane_name
            + " Orientation (from Earth Axes to the First Airplane's Body Axes Using "
            + "an Intrinsic zy'x\" Sequence) vs. Time",
            "Orientation (deg)",
            save,
            airplane_name_snake + "_orientation.png",
        )
        _plot_state_history(
            state_times,
            [omegasDeg_BP1__E[0], omegasDeg_BP1__E[1], omegasDeg_BP1__E[2]],
            ["Roll Rate", "Pitch Rate", "Yaw Rate"],
            [_angular_x_color, _angular_y_color, _angular_z_color],
            airplane_name
            + " Angular Velocity (of the First Airplane's Body Axes, Observed from the "
            + "Earth Frame) vs. Time",
            "Angular Velocity (deg/s)",
            save,
            airplane_name_snake + "_angular_velocity.png",
        )
        _plot_state_history(
            state_times,
            [alphas, betas],
            ["Angle of Attack", "Sideslip Angle"],
            [_alpha_color, _beta_color],
            airplane_name + " Aerodynamic Angles vs. Time",
            "Angle (deg)",
            save,
            airplane_name_snake + "_aerodynamic_angles.png",
        )

    # If the user wants to show the plots, do so. This is done outside the loop so
    # that plt.show() is only called once after all figures are created.
    if show:
        plt.show()
    else:
        plt.close("all")


# TEST: Consider adding integration tests for this function.
def log_results(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
) -> None:
    """Logs a solver's load and load coefficients.

    The logging level must be set to INFO or lower in order to see results. See
    set_up_logging for details on configuring the logging level.

    For a FreeFlightUnsteadyRingVortexLatticeMethodSolver, this also logs the first
    Airplane's initial and final six-degree-of-freedom state: its position, velocity,
    orientation, angular velocity, and aerodynamic angles.

    :param solver: The solver whose load and load coefficients will be logged. This can
        be a SteadyHorseshoeVortexLatticeMethodSolver, a
        SteadyRingVortexLatticeMethodSolver, or an
        UnsteadyRingVortexLatticeMethodSolver. The
        UnsteadyRingVortexLatticeMethodSolver's subclasses, the
        AeroelasticUnsteadyRingVortexLatticeMethodSolver and the
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, are also accepted. For a
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, the first Airplane's initial
        and final state is logged as well.
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
            "solver must be a SteadyHorseshoeVortexLatticeMethodSolver, a "
            "SteadyRingVortexLatticeMethodSolver, or an "
            "UnsteadyRingVortexLatticeMethodSolver."
        )

    if not solver.ran:
        raise RuntimeError("solver must have run before logging results.")

    padding_spaces = 2

    col1 = [
        "fX_W",
        "fY_W",
        "fZ_W",
        "mX_W_CgP1",
        "mY_W_CgP1",
        "mZ_W_CgP1",
        "cFX_W",
        "cFY_W",
        "cFZ_W",
        "cMX_W_CgP1",
        "cMY_W_CgP1",
        "cMZ_W_CgP1",
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
                title1 = f"{pad}Forces (in Wind Axes):"
                title2 = (
                    f"{pad}Moments (in Wind Axes, Relative to the First Airplane's CG):"
                )
                title3 = f"{pad}Force Coefficients (in Wind Axes):"
                title4 = (
                    f"{pad}Moment Coefficients (in Wind Axes, Relative to the First "
                    f"Airplane's CG):"
                )

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

                title1 = f"{pad}Final Forces (in Wind Axes):"
                title2 = (
                    f"{pad}Final Moments (in Wind Axes, Relative to the First "
                    f"Airplane's CG):"
                )
                title3 = f"{pad}Final Force Coefficients (in Wind Axes):"
                title4 = (
                    f"{pad}Final Moment Coefficients (in Wind Axes, Relative to "
                    f"the First Airplane's CG):"
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

                title1 = f"{pad}Final Cycle-Averaged Forces (in Wind Axes):"
                title2 = (
                    f"{pad}Final Cycle-Averaged Moments (in Wind Axes, Relative to the "
                    f"First Airplane's CG):"
                )
                title3 = f"{pad}Final Cycle-Averaged Force Coefficients (in Wind Axes):"
                title4 = (
                    f"{pad}Final Cycle-Averaged Moment Coefficients (in Wind Axes, "
                    f"Relative to the First Airplane's CG):"
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

        _logger.info(f'Airplane "{airplane.name}":')

        # Display the Reynolds number for steady solvers.
        if solver_type == "steady":
            assert isinstance(
                solver,
                (
                    steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
                    steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
                ),
            )
            re = solver.reynolds_numbers[airplane_num]
            _logger.info(f"{pad}Reynolds Number: {re:.2e}")

        for i in range(len(col1)):
            if i == 0:
                _logger.info(title1)
            elif i == 3:
                _logger.info(title2)
            elif i == 6:
                _logger.info(title3)
            elif i == 9:
                _logger.info(title4)

            s = f"{2 * pad}{col1[i]:<{col1_space}}{col2[i]:<{col2_space}}{col3[i]:<{col3_space}}{col4[i]}"
            _logger.info(s)

        # If the results from more Airplanes are going to be logged, log a blank
        # line to separate them.
        if (airplane_num + 1) < solver.num_airplanes:
            _logger.info("")

    # For a free flight solver, also log the first Airplane's initial and final
    # six-degree-of-freedom state. This is logged once, since the state describes the
    # first Airplane, the single rigid body the dynamics integrate.
    if isinstance(
        solver,
        free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    ):
        # Narrow the movement to a FreeFlightMovement so its operating point history is
        # typed. The solver type guarantees this, so the assert documents the invariant.
        movement = solver.unsteady_problem.movement
        assert isinstance(movement, free_flight_movement_mod.FreeFlightMovement)
        operating_points = movement.operating_point_movement.operating_points

        final_time = solver.delta_time * (len(operating_points) - 1)

        _logger.info("")
        _logger.info("The First Airplane's Free Flight State History:")

        # Each vector state quantity is broken into one row per component, mirroring
        # the per-component force and moment rows above, and each component row is
        # labeled with its variable-convention name. The four group headers are logged
        # before their first component (at flat indices 0, 3, 6, and 9), and the
        # component labels are padded to a common width so the values align.
        state_component_labels = [
            "cgP1X_E_Eo",
            "cgP1Y_E_Eo",
            "cgP1Z_E_Eo",
            "angleX_E_to_BP1_izyx",
            "angleY_E_to_BP1_izyx",
            "angleZ_E_to_BP1_izyx",
            "vCgP1X_E__E",
            "vCgP1Y_E__E",
            "vCgP1Z_E__E",
            "omegaX_BP1__E",
            "omegaY_BP1__E",
            "omegaZ_BP1__E",
        ]
        state_component_labels = [label + ":" for label in state_component_labels]
        state_component_units = ["m"] * 3 + ["deg"] * 3 + ["m/s"] * 3 + ["deg/s"] * 3
        state_component_space = (
            max(len(label) for label in state_component_labels) + 2 * padding_spaces
        )

        state_group_header_position = (
            f"{2 * pad}Position (of the First Airplane's CG, in Earth Axes, Relative "
            f"to the Earth Origin):"
        )
        state_group_header_orientation = (
            f"{2 * pad}Orientation (from Earth Axes to the First Airplane's Body Axes "
            f"Using an Intrinsic zy'x\" Sequence):"
        )
        state_group_header_velocity = (
            f"{2 * pad}Velocity (of the First Airplane's CG, in Earth Axes, Observed "
            f"from the Earth Frame):"
        )
        state_group_header_angular_velocity = (
            f"{2 * pad}Angular Velocity (of the First Airplane's Body Axes, Observed "
            f"from the Earth Frame):"
        )

        # Log the initial state (at time step 0) and the final state.
        for state_label, state_time, this_operating_point in [
            ("Initial State", 0.0, operating_points[0]),
            ("Final State", final_time, operating_points[-1]),
        ]:
            CgP1_E_Eo = this_operating_point.CgP1_E_Eo
            angles_E_to_BP1_izyx = this_operating_point.angles_E_to_BP1_izyx
            velocity_E__E = _velocity_E__E_from_operating_point(this_operating_point)
            omegas_BP1__E = this_operating_point.omegas_BP1__E

            state_component_values = [
                CgP1_E_Eo[0],
                CgP1_E_Eo[1],
                CgP1_E_Eo[2],
                angles_E_to_BP1_izyx[0],
                angles_E_to_BP1_izyx[1],
                angles_E_to_BP1_izyx[2],
                velocity_E__E[0],
                velocity_E__E[1],
                velocity_E__E[2],
                omegas_BP1__E[0],
                omegas_BP1__E[1],
                omegas_BP1__E[2],
            ]
            state_component_values = [
                f"{np.round(value, 3)} {unit}"
                for value, unit in zip(state_component_values, state_component_units)
            ]

            _logger.info(f"{pad}{state_label} (at t = {state_time:.3f} s):")

            for i in range(len(state_component_labels)):
                if i == 0:
                    _logger.info(state_group_header_position)
                elif i == 3:
                    _logger.info(state_group_header_orientation)
                elif i == 6:
                    _logger.info(state_group_header_velocity)
                elif i == 9:
                    _logger.info(state_group_header_angular_velocity)

                _logger.info(
                    f"{3 * pad}{state_component_labels[i]:<{state_component_space}}"
                    f"{state_component_values[i]}"
                )

            # The aerodynamic angles are scalars, so each is logged as a single row
            # labeled with its variable-convention name. The labels are padded to a
            # common width so the two values align with each other.
            alpha_label = "Angle of Attack (alpha):"
            beta_label = "Sideslip Angle (beta):"
            aerodynamic_angle_space = (
                max(len(alpha_label), len(beta_label)) + padding_spaces
            )
            _logger.info(
                f"{2 * pad}{alpha_label:<{aerodynamic_angle_space}}"
                f"{np.round(this_operating_point.alpha, 3)} deg"
            )
            _logger.info(
                f"{2 * pad}{beta_label:<{aerodynamic_angle_space}}"
                f"{np.round(this_operating_point.beta, 3)} deg"
            )


def _velocity_E__E_from_operating_point(
    this_operating_point: operating_point_mod.OperatingPoint,
) -> np.ndarray:
    """Returns the first Airplane's CG velocity (in Earth axes, observed from the Earth
    frame) for a free flight OperatingPoint.

    The CG velocity is the negative of the freestream velocity, since the freestream (a
    still airmass) is entirely due to the first Airplane's motion. The OperatingPoint
    stores the freestream velocity in the first Airplane's geometry axes, so it is
    rotated into Earth axes here.

    :param this_operating_point: The OperatingPoint whose CG velocity will be returned.
    :return: A (3,) ndarray of floats representing the first Airplane's CG velocity (in
        Earth axes, observed from the Earth frame) in meters per second.
    """
    vInf_E__E = _transformations.apply_T_to_vectors(
        this_operating_point.T_pas_GP1_CgP1_to_E_CgP1,
        this_operating_point.vInf_GP1__E,
        is_position=False,
    )
    return -vInf_E__E


def _plot_state_history(
    times: np.ndarray,
    series: list[np.ndarray],
    labels: list[str],
    colors: list[str],
    title: str,
    y_label: str,
    save: bool,
    save_name: str,
) -> None:
    """Plots one free flight state-history figure: a set of series sharing a y-axis,
    plotted against time and styled to match the load plots of plot_results_versus_time.

    :param times: A (num_steps,) ndarray of floats representing the time, in seconds, at
        each time step.
    :param series: A list of (num_steps,) ndarrays of floats, one per line to plot.
    :param labels: A list of the legend labels, one per series.
    :param colors: A list of the line colors, one per series.
    :param title: The figure's title.
    :param y_label: The figure's y-axis label.
    :param save: Set this to True to save the figure as a PNG.
    :param save_name: The file name to save the figure under if save is True.
    :return: None
    """
    figure, axes = plt.subplots()

    # Remove the plot's top and right spines.
    axes.spines.right.set_visible(False)
    axes.spines.top.set_visible(False)

    # Format the plot's spine and label colors.
    axes.spines.bottom.set_color(_text_color_normalized)
    axes.spines.left.set_color(_text_color_normalized)
    axes.xaxis.label.set_color(_text_color_normalized)
    axes.yaxis.label.set_color(_text_color_normalized)

    # Format the plot's tick colors.
    axes.tick_params(axis="x", colors=_text_color_normalized)
    axes.tick_params(axis="y", colors=_text_color_normalized)

    # Format the plot's background colors.
    figure.patch.set_facecolor(_figure_background_color)
    axes.set_facecolor(_figure_background_color)

    # Populate the plot, staggering the markers across the series.
    num_series = len(series)
    for series_id, (this_series, label, color) in enumerate(
        zip(series, labels, colors)
    ):
        axes.plot(
            times,
            this_series,
            label=label,
            color=color,
            marker=".",
            markevery=(_marker_spacing * series_id / num_series, _marker_spacing),
            markersize=_marker_size,
        )

    # Name the plot's axis labels and title.
    axes.set_xlabel("Time (s)", color=_text_color_normalized)
    axes.set_ylabel(y_label, color=_text_color_normalized)
    axes.set_title(title, color=_text_color_normalized)

    # Format the plot's legend.
    axes.legend(
        facecolor=_figure_background_color,
        edgecolor=_figure_background_color,
        labelcolor=_text_color_normalized,
    )

    # Save the figure as a PNG if the user wants to do so.
    if save:
        figure.savefig(save_name, dpi=300)


def _get_panel_surfaces(
    airplanes: tuple[geometry.airplane.Airplane, ...],
) -> pv.PolyData:
    """Returns a PolyData representation of the Wings' Panels' surfaces associated with
    all the Airplanes in a tuple of Airplanes.

    :param airplanes: The tuple of Airplanes whose Wings' Panels' surfaces will be
        returned.
    :return: A PolyData representation of the Airplanes' Wings' Panels' surfaces.
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
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # Update the number of Panels.
                panel_num += 1

    # Return the Panels' surfaces.
    return pv.PolyData(panel_vertices, panel_faces)


def _get_image_surface_mesh_and_texture(
    this_operating_point: operating_point_mod.OperatingPoint,
    geometry_bounds: tuple[float, float, float, float, float, float],
) -> tuple[pv.PolyData, pv.Texture] | None:
    """Returns a PolyData plane mesh and checkerboard Texture representing the image
    surface, or None if no image surface is defined.

    The plane is centered at the projection of the geometry bounding box center onto the
    image surface, and sized proportionally to the bounding box diagonal so that it
    appears large relative to the geometry.

    :param this_operating_point: The OperatingPoint that may define an image surface.
    :param geometry_bounds: The (xmin, xmax, ymin, ymax, zmin, zmax) bounding box of the
        geometry used to determine the plane's center and size.
    :return: A tuple of (PolyData plane mesh, checkerboard Texture) representing the
        image surface, or None if no image surface is defined.
    """
    surface_normal = this_operating_point.surfaceNormal_GP1
    surface_point = this_operating_point.surfacePoint_GP1_CgP1

    if surface_normal is None or surface_point is None:
        return None

    # Compute the bounding box center and diagonal length.
    bounds = np.array(geometry_bounds, dtype=float)
    bbox_center = np.array(
        [
            0.5 * (bounds[0] + bounds[1]),
            0.5 * (bounds[2] + bounds[3]),
            0.5 * (bounds[4] + bounds[5]),
        ],
        dtype=float,
    )
    bbox_diagonal = float(
        np.linalg.norm(
            np.array(
                [
                    bounds[1] - bounds[0],
                    bounds[3] - bounds[2],
                    bounds[5] - bounds[4],
                ],
                dtype=float,
            )
        )
    )

    # Project the bounding box center onto the image surface to get the plane's center.
    offset = np.dot(bbox_center - surface_point, surface_normal)
    plane_center = bbox_center - offset * surface_normal

    # Size the plane proportionally to the bounding box diagonal.
    plane_size = _image_surface_scale * bbox_diagonal

    mesh = pv.Plane(
        center=plane_center,
        direction=surface_normal,
        i_size=plane_size,
        j_size=plane_size,
    )

    # Build a checkerboard texture image. Each cell is one pixel, so a 25 x 25
    # checkerboard is a 25 x 25 x 3 RGB image.
    n = _image_surface_checker_size
    row = np.arange(n, dtype=int)
    col = np.arange(n, dtype=int)
    rr, cc = np.meshgrid(row, col, indexing="ij")
    is_dark = (rr + cc) % 2 == 0
    image = np.where(
        is_dark[:, :, np.newaxis], _image_surface_color_a, _image_surface_color_b
    )
    texture = pv.numpy_to_texture(image)

    return mesh, texture


def _get_T_pas_GP1_CgP1_to_E_Eo(
    this_operating_point: operating_point_mod.OperatingPoint,
) -> np.ndarray:
    """Returns the passive transformation from the first Airplane's geometry axes,
    relative to the first Airplane's CG, to Earth axes, relative to the Earth origin.

    Free flight visualizations render each time step's geometry in its true Earth-frame
    position and orientation, so the body translates and rotates through the scene. This
    transformation chains the per-step geometry-to-Earth rotation with the translation
    from the first Airplane's CG to the Earth origin.

    :param this_operating_point: The OperatingPoint whose Earth-frame position and
        orientation define the transformation.
    :return: A (4,4) ndarray of floats representing the passive transformation from the
        first Airplane's geometry axes, relative to the first Airplane's CG, to Earth
        axes, relative to the Earth origin.
    """
    # Translate from Earth axes relative to the first Airplane's CG to Earth axes
    # relative to the Earth origin. For a passive translation, the parameter is the
    # position of the final reference point (the Earth origin) relative to the initial
    # one (the first Airplane's CG), which is the negative of the CG's position relative
    # to the Earth origin.
    T_pas_E_CgP1_to_E_Eo = _transformations.generate_trans_T(
        translations=-this_operating_point.CgP1_E_Eo,
        passive=True,
    )

    # Chain geometry-to-Earth (relative to the CG) with the CG-to-origin translation.
    return _transformations.compose_T_pas(
        this_operating_point.T_pas_GP1_CgP1_to_E_CgP1,
        T_pas_E_CgP1_to_E_Eo,
    )


def _transform_mesh(
    mesh: pv.PolyData,
    T_pas: np.ndarray,
) -> pv.PolyData:
    """Returns a copy of a PolyData mesh with its points mapped through a passive
    transformation.

    :param mesh: The PolyData mesh to transform.
    :param T_pas: A (4,4) ndarray of floats representing the passive transformation to
        apply to the mesh's points.
    :return: A new PolyData mesh with all points mapped through the transformation.
    """
    transformed = mesh.copy()
    transformed.points = _transformations.apply_T_to_vectors(
        T_pas,
        mesh.points,
        is_position=True,
    )
    return transformed


def _reflect_mesh(
    mesh: pv.PolyData,
    T_reflect: np.ndarray,
) -> pv.PolyData:
    """Returns a copy of a PolyData mesh with its points reflected across the image
    surface.

    :param mesh: The PolyData mesh to reflect.
    :param T_reflect: A (4,4) ndarray of floats representing the active reflection
        transformation matrix (in the first Airplane's geometry axes, relative to the
        first Airplane's CG).
    :return: A new PolyData mesh with all points reflected across the image surface.
    """
    return _transform_mesh(mesh, T_reflect)


def _free_flight_fit_parallel_scale(
    meshes: list[pv.PolyData],
    focalPoint_E_Eo: np.ndarray,
    viewDirection_E: np.ndarray,
    viewUp_E: np.ndarray,
    margin: float = 1.15,
) -> float:
    """Returns a parallel-projection scale that frames a set of meshes about a focal
    point.

    The scale is half the viewport's height in world units (the convention for a
    parallel projection). It is sized to the largest projection of the meshes' bounding-
    box corners onto the camera's screen-right and screen-up axes, measured from the
    focal point, so everything stays in view regardless of the viewport's aspect ratio.
    A margin leaves a little space around the geometry.

    :param meshes: The PolyData meshes to frame.
    :param focalPoint_E_Eo: A (3,) ndarray of floats locating the camera's focal point
        (in Earth axes, relative to the Earth origin). Extents are measured from this
        point, since it projects to the center of the viewport.
    :param viewDirection_E: A (3,) ndarray of floats giving the offset from the focal
        point to the camera position (in Earth axes). The camera looks back along it.
    :param viewUp_E: A (3,) ndarray of floats giving the camera's up direction (in Earth
        axes).
    :param margin: A factor (at least 1.0) by which to pad the fitted scale. The default
        is 1.15.
    :return: The parallel-projection scale.
    """
    # Build the camera's screen right and up axes in Earth axes. The camera looks from its
    # position back toward the focal point, i.e. along the negative view direction. The up
    # axis is the supplied up made orthogonal to that look direction.
    lookDirection_E = -viewDirection_E
    lookDirection_E = lookDirection_E / np.linalg.norm(lookDirection_E)
    upDirection_E = viewUp_E - np.dot(viewUp_E, lookDirection_E) * lookDirection_E
    upDirection_E = upDirection_E / np.linalg.norm(upDirection_E)
    rightDirection_E = np.cross(lookDirection_E, upDirection_E)
    rightDirection_E = rightDirection_E / np.linalg.norm(rightDirection_E)

    # Collect every mesh's eight bounding-box corners, measured from the focal point.
    corners_E_Eo: list[list[float]] = []
    for mesh in meshes:
        bounds = np.array(mesh.bounds, dtype=float)
        corners_E_Eo.extend(
            [float(x), float(y), float(z)]
            for x in bounds[0:2]
            for y in bounds[2:4]
            for z in bounds[4:6]
        )
    corners_E = np.array(corners_E_Eo) - focalPoint_E_Eo

    # Each matrix-vector product is a batch of dot products, one per corner: the (N, 3)
    # array of corners times a (3,) screen axis gives an (N,) array whose entries are each
    # corner's signed distance from the focal point along that axis (the axes are unit
    # vectors, so the dot product is the scalar projection). Take the largest magnitude
    # along either axis so neither screen dimension is clipped.
    half_extent_right = float(np.abs(corners_E @ rightDirection_E).max())
    half_extent_up = float(np.abs(corners_E @ upDirection_E).max())
    return margin * max(half_extent_right, half_extent_up)


def _mute_color(
    color: str | tuple[float, ...],
    factor: float,
) -> tuple[float, float, float]:
    """Returns a muted version of a color by linearly interpolating it toward middle
    gray.

    :param color: Any color that PyVista can parse (name, hex string, RGB tuple, etc.).
    :param factor: The muting factor in [0, 1]. 0 means no change, 1 means fully gray.
    :return: A (R, G, B) tuple of floats in [0, 1].
    """
    rgb = np.array(pv.Color(color).float_rgb)
    gray = np.full(3, 0.5)
    muted = rgb + factor * (gray - rgb)
    return float(muted[0]), float(muted[1]), float(muted[2])


def _mute_colormap(
    cmap_name: str,
    factor: float,
) -> matplotlib.colors.ListedColormap:
    """Returns a muted version of a named colormap by linearly interpolating each color
    toward middle gray.

    :param cmap_name: The name of a Matplotlib or cmocean colormap.
    :param factor: The muting factor in [0, 1]. 0 means no change, 1 means fully gray.
    :return: A ListedColormap with muted colors.
    """
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        cmap = plt.get_cmap("cmo." + cmap_name)
    colors = cmap(np.linspace(0, 1, 256))
    gray = 0.5
    colors[:, :3] = colors[:, :3] + factor * (gray - colors[:, :3])
    return matplotlib.colors.ListedColormap(colors)


def _get_wake_ring_vortex_surfaces(
    solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    step: int,
) -> pv.PolyData:
    """Returns the PolyData representation of the surfaces of an
    UnsteadyRingVortexLatticeMethodSolver's wake ring vortices at a given time step.

    :param solver: The UnsteadyRingVortexLatticeMethodSolver with the wake ring vortices
        to process.
    :param step: The time step at which to process the wake ring vortices.
    :return: The PolyData representation of the wake ring vortices.
    """
    num_wake_ring_vortices = solver.list_num_wake_vortices[step]
    stackFrwrvp_GP1_CgP1 = solver.listStackFrwrvp_GP1_CgP1[step]
    stackFlwrvp_GP1_CgP1 = solver.listStackFlwrvp_GP1_CgP1[step]
    stackBlwrvp_GP1_CgP1 = solver.listStackBlwrvp_GP1_CgP1[step]
    stackBrwrvp_GP1_CgP1 = solver.listStackBrwrvp_GP1_CgP1[step]

    # Initialize empty ndarrays to hold each wake ring vortex's vertices and face.
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

        # Stack this wake ring vortex's vertices and faces to the ndarrays of all wake
        # ring vortices' vertices and faces.
        wake_ring_vortex_vertices = np.vstack(
            (wake_ring_vortex_vertices, wake_ring_vortex_vertices_to_add)
        )
        wake_ring_vortex_faces = np.hstack(
            (wake_ring_vortex_faces, wake_ring_vortex_face_to_add)
        )

    # Return the wake ring vortex surfaces.
    return pv.PolyData(wake_ring_vortex_vertices, wake_ring_vortex_faces)


def _get_scalars(
    airplanes: tuple[geometry.airplane.Airplane, ...],
    scalar_type: str,
    qInf__E: float,
) -> np.ndarray:
    """Returns the load coefficient values from a SteadyProblem's Airplanes' Wings'
    Panels.

    :param airplanes: The tuple of Airplanes with the scalars to return.
    :param scalar_type: Determines which load coefficient to return as scalars. Can be
        "induced drag", "side force", or "lift", which respectively use each Panel's
        induced drag, side force, and lift coefficient.
    :param qInf__E: The current freestream dynamic pressure experienced by this
        SteadyProblem's Airplane(s) (observed in the Earth frame). The units are in
        Pascals.
    :return: A (N,) ndarray of floats representing the N Panels' load coefficients.
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

                    scalars = np.hstack((scalars, this_induced_drag_coefficient))

                if scalar_type == "side force":
                    this_side_force_coefficient = (
                        this_panel.forces_W[1] / qInf__E / this_panel.area
                    )

                    scalars = np.hstack((scalars, this_side_force_coefficient))

                if scalar_type == "lift":
                    this_lift_coefficient = (
                        -this_panel.forces_W[2] / qInf__E / this_panel.area
                    )

                    scalars = np.hstack((scalars, this_lift_coefficient))

    # Return the resulting ndarray of scalars.
    return scalars


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
    text_color: tuple[int, int, int] = _text_color,
) -> None:
    """Plots a scalar bar, the surfaces of a set of Panels with particular scalars, and
    labels for the minimum and maximum scalar values.

    :param plotter: The Plotter used for visualization.
    :param these_scalars: A (N,) ndarray of floats representing the N Panels' load
        coefficients.
    :param scalar_type: Which load coefficient is represented by the scalars. Can be
        "induced drag", "side force", or "lift".
    :param min_scalar: Minimum scalar value, which is displayed as text on the Plotter.
    :param max_scalar: Maximum scalar value, which is displayed as text on the Plotter.
    :param color_map: Name of the color map to use for scalar visualization. Check the
        pyvista.add_mesh documentation for the list of acceptable values.
    :param c_min: Lower bound for the color map scaling.
    :param c_max: Upper bound for the color map scaling.
    :param panel_surfaces: PolyData representing the Panels' surfaces.
    :param text_color: The color used for the scalar bar and label text. The default is
        _text_color.
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
        color=text_color,
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
        color=text_color,
    )
    plotter.add_text(
        text="Min: " + str(min_scalar),
        position=_text_min_position,
        font_size=_text_font_size,
        viewport=True,
        color=text_color,
    )
