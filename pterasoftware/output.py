"""Contains functions for visualizing geometry and results.

**Contains the following classes:**

None

**Contains the following functions:**

draw: Draws a solver's Airplane(s).

animate: Animates an UnsteadyRingVortexLatticeMethodSolver's Airplane(s).

animate_free_flight: Animates a CoupledUnsteadyRingVortexLatticeMethodSolver's Airplane.

plot_results_versus_time: Plots an UnsteadyRingVortexLatticeMethodSolver's loads and
load coefficients as a function of time.

print_results: Prints a solver's load and load coefficients.
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
    _parameter_validation,
    _transformations,
    coupled_unsteady_ring_vortex_lattice_method,
    geometry,
)
from . import operating_point as operating_point_mod
from . import (
    steady_horseshoe_vortex_lattice_method,
    steady_ring_vortex_lattice_method,
    unsteady_ring_vortex_lattice_method,
)

# Define the color and colormaps used by the visualization functions.
_sequential_color_map = "speed"
_diverging_color_map = "delta"
_wake_vortex_color = "white"
_panel_color = "chartreuse"
_streamline_color = "orchid"
_image_surface_opacity = 0.75
_image_surface_scale = 5.0
_image_reflection_mute_factor = 0.5
_image_surface_checker_size = 25
_image_surface_color_a = np.array([40, 40, 40], dtype=np.uint8)
_image_surface_color_b = np.array([80, 80, 80], dtype=np.uint8)
_plotter_background_color = "black"
_body_mesh_color = "cornflowerblue"
_ground_plane_checker_light = [200, 200, 200]
_ground_plane_checker_dark = [120, 120, 120]
_ground_plane_checker_block_size = 32
_ground_plane_checker_resolution = 256
_ground_plane_opacity = 0.5
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
        | coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver
    ),
    scalar_type: str | None = None,
    show_streamlines: bool | np.bool_ = False,
    show_wake_vortices: bool | np.bool_ = False,
    save: bool | np.bool_ = False,
    testing: bool | np.bool_ = False,
    body_mesh_path: str | None = None,
    body_mesh_scale: float | int = 1.0,
    ground_plane_z_E: float | int | None = None,
) -> None:
    """Draws a solver's Airplane(s).

    **Citation:**

    Adapted from: vlm3.draw in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 03/28/2020

    :param solver: The solver whose Airplane(s) will be plotted.
    :param scalar_type: Determines how to color the Panels. Setting this to None colors
        the Panels uniformly. If the solver has been run, it can also be "induced drag",
        "side force", or "lift", which respectively use each Panel's induced drag, side
        force, and lift coefficient. The default is None.
    :param show_streamlines: Set this to True to show the streamlines emanating from the
        back of the Wings. If True, the solver's streamlines must have already been
        calculated. Can be a bool or a numpy bool and will be converted internally to a
        bool. The default is False.
    :param show_wake_vortices: Set this to True to show any wake RingVortices. If True,
        the solver must be an UnsteadyRingVortexLatticeMethodSolver and must have
        already been run. Can be a bool or a numpy bool and will be converted internally
        to a bool. The default is False.
    :param save: Set this to True to save the image as a WebP. It can be a bool or a
        numpy bool and will be converted internally to a bool. The default is False.
    :param testing: Set this to True to close the image after one second, which is
        useful for running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :param body_mesh_path: The path to an STL file containing the body mesh to display.
        The mesh is assumed to be in the first Airplane's body axes, relative to the
        first Airplane's CG. Only used when the solver is a
        CoupledUnsteadyRingVortexLatticeMethodSolver. Setting this to None disables body
        mesh rendering. The default is None.
    :param body_mesh_scale: The scale factor to apply to the body mesh before rendering.
        Only used when body_mesh_path is not None. The default is 1.0.
    :param ground_plane_z_E: The z coordinate of the ground plane (in Earth axes, in
        meters). Only used when the solver is a
        CoupledUnsteadyRingVortexLatticeMethodSolver. Setting this to None disables
        ground plane rendering. The default is None.
    :return: None
    """
    if not isinstance(
        solver,
        (
            steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
            unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
            coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
        ),
    ):
        raise TypeError(
            "solver must be a SteadyHorseshoeVortexLatticeMethodSolver, "
            "a SteadyRingVortexLatticeMethodSolver, "
            "an UnsteadyRingVortexLatticeMethodSolver, "
            "or a CoupledUnsteadyRingVortexLatticeMethodSolver."
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
    if show_streamlines and isinstance(
        solver,
        coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    ):
        raise ValueError(
            "show_streamlines can't be True when drawing a "
            "CoupledUnsteadyRingVortexLatticeMethodSolver."
        )
    if show_streamlines:
        assert isinstance(
            solver,
            (
                steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
                steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
                unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
            ),
        )
        if len(solver.gridStreamlinePoints_GP1_CgP1) == 0:
            raise RuntimeError(
                "solver must have streamline points calculated before drawing with "
                "show_streamlines set to True."
            )

    show_wake_vortices = _parameter_validation.boolLike_return_bool(
        show_wake_vortices, "show_wake_vortices"
    )
    if show_wake_vortices and not isinstance(
        solver,
        (
            unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
            coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
        ),
    ):
        raise ValueError(
            "show_wake_vortices can only be True when drawing an "
            "UnsteadyRingVortexLatticeMethodSolver or a "
            "CoupledUnsteadyRingVortexLatticeMethodSolver."
        )
    if show_wake_vortices and not solver.ran:
        raise RuntimeError(
            "solver must have run before drawing with show_wake_vortices set to True."
        )

    save = _parameter_validation.boolLike_return_bool(save, "save")
    testing = _parameter_validation.boolLike_return_bool(testing, "testing")

    if body_mesh_path is not None:
        if not isinstance(
            solver,
            coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
        ):
            raise ValueError(
                "body_mesh_path can only be used when the solver is a "
                "CoupledUnsteadyRingVortexLatticeMethodSolver."
            )
        body_mesh_path = _parameter_validation.str_return_str(
            body_mesh_path, "body_mesh_path"
        )
    body_mesh_scale = _parameter_validation.number_in_range_return_float(
        body_mesh_scale, "body_mesh_scale"
    )
    if ground_plane_z_E is not None:
        if not isinstance(
            solver,
            coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
        ):
            raise ValueError(
                "ground_plane_z_E can only be used when the solver is a "
                "CoupledUnsteadyRingVortexLatticeMethodSolver."
            )
        ground_plane_z_E = _parameter_validation.number_in_range_return_float(
            ground_plane_z_E, "ground_plane_z_E"
        )

    # If a body mesh path was provided, load the STL and pre-scale its points.
    body_mesh_BP1_CgP1: pv.PolyData | None = None
    if body_mesh_path is not None:
        body_mesh_BP1_CgP1 = pv.read(body_mesh_path)
        body_mesh_BP1_CgP1.points = body_mesh_BP1_CgP1.points * body_mesh_scale

    # Create the Plotter and set it to use parallel projection (instead of perspective).
    plotter = pv.Plotter(window_size=_window_size, lighting=None)
    plotter.enable_parallel_projection()  # type: ignore[call-arg]

    # Get the solver's geometry and OperatingPoint.
    if isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        draw_step = solver.num_steps - 1

        airplanes = solver.steady_problems[draw_step].airplanes
        draw_operating_point = solver.steady_problems[draw_step].operating_point
        qInf__E = draw_operating_point.qInf__E

        # If showing wake RingVortices, get their surfaces and plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces(
                solver, draw_step
            )
            if wake_ring_vortex_surfaces.n_points > 0:
                plotter.add_mesh(
                    wake_ring_vortex_surfaces,
                    show_edges=True,
                    smooth_shading=False,
                    color=_wake_vortex_color,
                )
    elif isinstance(
        solver,
        coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    ):
        draw_step = solver.num_steps - 1

        airplanes = (solver.current_airplane,)
        qInf__E = solver.current_coupled_operating_point.qInf__E

        # If showing wake RingVortices, get their surfaces and plot them.
        if show_wake_vortices:
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces_free_flight(
                solver, draw_step
            )
            if wake_ring_vortex_surfaces.n_points > 0:
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

    # Get the Panel surfaces.
    if isinstance(
        solver,
        coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    ):
        panel_surfaces = _get_panel_surfaces_free_flight(
            airplanes[0],
            position_E_E=solver.stackPosition_E_E[-1, :],
            R_pas_E_to_BP1=solver.stackR_pas_E_to_BP1[-1, :],
        )
    else:
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

    # Add the body mesh, ground plane, and camera setup for a coupled solver.
    draw_cpos: tuple | None = None
    if isinstance(
        solver,
        coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    ):
        # Compute the airplane's bounding diagonal to use as a scale reference for
        # the camera and the ground plane texture tiling.
        airplane_bounds = np.array(panel_surfaces.bounds)
        airplane_diagonal = float(
            np.linalg.norm(airplane_bounds[1::2] - airplane_bounds[::2])
        )

        # Add the body mesh at the final time step.
        if body_mesh_BP1_CgP1 is not None:
            body_surface = _get_body_mesh_surface_free_flight(
                body_mesh_BP1_CgP1,
                solver.stackPosition_E_E[-1, :],
                solver.stackR_pas_E_to_BP1[-1, :],
            )
            plotter.add_mesh(
                body_surface,
                color=_body_mesh_color,
                smooth_shading=True,
            )

        # Add the ground plane.
        if ground_plane_z_E is not None:
            initialPosition_E_E = solver.stackPosition_E_E[0]
            finalPosition_E_E = solver.stackPosition_E_E[-1]
            trajectoryMidpoint_E_E = (initialPosition_E_E + finalPosition_E_E) / 2.0
            trajectory_extent = float(
                np.linalg.norm(finalPosition_E_E - initialPosition_E_E)
            )
            ground_plane_surface = _get_ground_plane_surface(
                ground_plane_z_E,
                trajectory_extent,
                trajectoryMidpoint_E_E,
                airplane_diagonal,
            )
            ground_plane_texture = _get_ground_plane_texture()
            plotter.add_mesh(
                ground_plane_surface,
                texture=ground_plane_texture,
                opacity=_ground_plane_opacity,
                smooth_shading=False,
            )

            # Set up lighting and shadows.
            _setup_free_flight_lighting(plotter)

        # Compute a camera position based on the airplane geometry (not the ground
        # plane) so that the view is zoomed to the airplane's scale.
        T_pas_E_to_V = _transformations.generate_rot_T(
            angles=(0.0, 180.0, 0.0),
            passive=True,
            intrinsic=True,
            order="xyz",
        )
        airplaneCenter_V = np.array(
            [
                (airplane_bounds[0] + airplane_bounds[1]) / 2.0,
                (airplane_bounds[2] + airplane_bounds[3]) / 2.0,
                (airplane_bounds[4] + airplane_bounds[5]) / 2.0,
            ]
        )
        padding = airplane_diagonal * 2.0
        camera_distance = airplane_diagonal + padding
        cameraDirection_V = np.array([-1.0, -1.0, 1.0])
        cameraDirection_V = cameraDirection_V / np.linalg.norm(cameraDirection_V)
        cameraPosition_V = airplaneCenter_V + camera_distance * cameraDirection_V
        draw_cpos = (
            tuple(cameraPosition_V),
            tuple(airplaneCenter_V),
            (0.0, 0.0, 1.0),
        )

    T_reflect = draw_operating_point.surfaceReflect_T_act_GP1_CgP1
    image_surface_mesh = None

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

        # Add reflected wake RingVortex surfaces if they are being shown.
        if show_wake_vortices:
            plotter.add_mesh(
                _reflect_mesh(wake_ring_vortex_surfaces, T_reflect),
                show_edges=True,
                edge_color=muted_edge_color,
                smooth_shading=False,
                color=_mute_color(_wake_vortex_color, mute),
            )
    # If showing streamlines, plot them.
    if show_streamlines and not isinstance(
        solver,
        coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    ):
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

                    # If an image surface is defined, add the reflected streamline
                    # segment.
                    if T_reflect is not None:
                        reflected_point = _transformations.apply_T_to_vectors(
                            T_reflect,
                            point,
                            has_point=True,
                        )
                        reflected_last_point = _transformations.apply_T_to_vectors(
                            T_reflect,
                            last_point,
                            has_point=True,
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
        plotter.camera.position = (-1, -1, 1)
        plotter.camera.focal_point = (0, 0, 0)
        plotter.camera.up = (0, 0, 1)
        plotter.reset_camera(bounds=geometry_bounds)  # type: ignore[call-arg]

    # Set the Plotter's background color.
    plotter.set_background(color=_plotter_background_color)  # type: ignore[call-arg]

    # Use the computed camera position for coupled solvers (which accounts for the
    # airplane's actual size). When an image surface is defined, use None (let PyVista
    # auto-fit to include the reflected geometry). Otherwise, use the default direction.
    if draw_cpos is not None:
        cpos_arg = draw_cpos
    elif image_surface_mesh is not None:
        cpos_arg = None
    else:
        cpos_arg = (-1, -1, 1)

    if not testing:
        # Show the Plotter so the user can adjust the camera position and window.
        # When the user closes the window, the Plotter still exists. Therefore,
        # it can later be saved as an image if desired.
        plotter.show(
            cpos=cpos_arg,
            full_screen=False,
            auto_close=False,
        )
    else:
        # Show the Plotter for 1 second, then proceed automatically. This is useful
        # for testing.
        plotter.show(
            cpos=cpos_arg,
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
    """Animates an UnsteadyRingVortexLatticeMethodSolver's Airplane(s).

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver whose Airplane(s)
        will be animated.
    :param scalar_type: Determines how to color the Panels. Setting this to None colors
        the Panels uniformly. If the solver has been run, it can also be "induced drag",
        "side force", or "lift", which respectively use each Panel's induced drag, side
        force, and lift coefficient. The default is None.
    :param show_wake_vortices: Set this to True to show any wake RingVortices. If True,
        the solver must have already been run. Can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :param save: Set this to True to save the image as a WebP. It can be a bool or a
        numpy bool and will be converted internally to a bool. The default is False.
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

    # Scale down the true-speed frames per second to at most 50 fps. This is the
    # maximum speed at which some programs can render WebPs.
    requested_fps = 1.0 / unsteady_solver.delta_time
    speed = 1.0
    if requested_fps > 50.0:
        speed = 50.0 / requested_fps
    actual_fps = float(math.floor(requested_fps * speed))
    # actual_fps = 30

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
    if T_reflect is not None:
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
        image_surface_geometry_bounds = None

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
        plotter.camera.position = (-1, -1, 1)
        plotter.camera.focal_point = (0, 0, 0)
        plotter.camera.up = (0, 0, 1)
        plotter.reset_camera(bounds=image_surface_geometry_bounds)  # type: ignore[call-arg]

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
            cpos=(-1, -1, 1) if image_surface_mesh is None else None,
            full_screen=False,
            auto_close=False,
        )
    else:
        plotter.show(
            title="Rendering speed not to scale.",
            cpos=(-1, -1, 1) if image_surface_mesh is None else None,
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
            if wake_ring_vortex_surfaces.n_points > 0:
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

            # Add reflected wake RingVortex surfaces if they are being shown.
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
def animate_free_flight(
    coupled_solver: coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    scalar_type: str | None = None,
    show_wake_vortices: bool = False,
    save: bool = False,
    testing: bool = False,
    body_mesh_path: str | None = None,
    body_mesh_scale: float | int = 1.0,
    ground_plane_z_E: float | int | None = None,
) -> None:
    """Animates a CoupledUnsteadyRingVortexLatticeMethodSolver's Airplane.

    :param coupled_solver: The CoupledUnsteadyRingVortexLatticeMethodSolver whose
        Airplane will be animated.
    :param scalar_type: Determines how to color the Panels. Setting this to None colors
        the Panels uniformly. If the solver has been run, it can also be "induced drag",
        "side force", or "lift", which respectively use each Panel's induced drag, side
        force, and lift coefficient. The default is None.
    :param show_wake_vortices: Set this to True to show any wake RingVortices. If True,
        the solver must have already been run. Can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :param save: Set this to True to save the image as a WebP. It can be a bool or a
        numpy bool and will be converted internally to a bool. The default is False.
    :param testing: Set this to True to start the animation after one second, which is
        useful for running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :param body_mesh_path: The path to an STL file containing the body mesh to display.
        The mesh is assumed to be in the first Airplane's body axes, relative to the
        first Airplane's CG. Setting this to None disables body mesh rendering. The
        default is None.
    :param body_mesh_scale: The scale factor to apply to the body mesh before rendering.
        The default is 1.0.
    :param ground_plane_z_E: The z coordinate of the ground plane (in Earth axes, in
        meters). Setting this to None disables ground plane rendering. The default is
        None.
    :return: None
    """
    if not isinstance(
        coupled_solver,
        coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    ):
        raise TypeError(
            "coupled_solver must be a CoupledUnsteadyRingVortexLatticeMethodSolver."
        )

    if scalar_type is not None:
        scalar_type = _parameter_validation.str_return_str(scalar_type, "scalar_type")
        if scalar_type not in ("induced drag", "side force", "lift"):
            raise ValueError(
                'scalar_type must be None, "induced drag", "side force", or "lift".'
            )

    show_wake_vortices = _parameter_validation.boolLike_return_bool(
        show_wake_vortices, "show_wake_vortices"
    )
    save = _parameter_validation.boolLike_return_bool(save, "save")
    testing = _parameter_validation.boolLike_return_bool(testing, "testing")

    if body_mesh_path is not None:
        body_mesh_path = _parameter_validation.str_return_str(
            body_mesh_path, "body_mesh_path"
        )
    body_mesh_scale = _parameter_validation.number_in_range_return_float(
        body_mesh_scale, "body_mesh_scale"
    )
    if ground_plane_z_E is not None:
        ground_plane_z_E = _parameter_validation.number_in_range_return_float(
            ground_plane_z_E, "ground_plane_z_E"
        )

    # If a body mesh path was provided, load the STL and pre-scale its points.
    body_mesh_BP1_CgP1: pv.PolyData | None = None
    if body_mesh_path is not None:
        body_mesh_BP1_CgP1 = pv.read(body_mesh_path)
        body_mesh_BP1_CgP1.points = body_mesh_BP1_CgP1.points * body_mesh_scale

    # Get the Airplane at each time step from the solver's CoupledSteadyProblems.
    # These are the Airplanes with forces computed on their Panels, as opposed to the
    # geometry only Airplanes from the CoupledMovement.
    airplanes = tuple(csp.airplane for csp in coupled_solver.coupled_steady_problems)

    # Scale down the true-speed frames per second to at most 50 fps. This is the
    # maximum speed at which some programs can render WebPs.
    requested_fps = 1.0 / coupled_solver.delta_time
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

    # Initialize variables to hold the CoupledSteadyProblems' scalars and their
    # attributes.
    all_scalars = np.empty(0, dtype=float)
    min_scalar = 0.0
    max_scalar = 0.0

    # If coloring the Panels based on scalars, gather all the scalars across all the
    # time steps and Airplanes. These will be used to set the color map limits.
    if scalar_type is not None:
        for step_id, airplane in enumerate(airplanes):
            scalars_to_add = _get_scalars(
                (airplane,),
                scalar_type,
                coupled_solver.coupled_steady_problems[
                    step_id
                ].coupled_operating_point.qInf__E,
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

    # Compute a camera position that can see the entire trajectory. Find the midpoint
    # and extent of the trajectory, then position the camera far enough back to see it
    # all.
    initialPosition_E_E = coupled_solver.stackPosition_E_E[0]
    finalPosition_E_E = coupled_solver.stackPosition_E_E[-1]
    trajectoryMidpoint_E_E = (initialPosition_E_E + finalPosition_E_E) / 2.0
    trajectory_extent = float(np.linalg.norm(finalPosition_E_E - initialPosition_E_E))

    # Compute the airplane's bounding diagonal at the first time step to use as a
    # scale reference for the camera. This ensures the padding adapts to the size of
    # the aircraft rather than using a hardcoded minimum.
    first_panel_surfaces_for_scale = _get_panel_surfaces_free_flight(
        airplanes[0],
        coupled_solver.stackPosition_E_E[0],
        coupled_solver.stackR_pas_E_to_BP1[0],
    )
    airplane_bounds = np.array(first_panel_surfaces_for_scale.bounds)
    airplane_diagonal = float(
        np.linalg.norm(airplane_bounds[1::2] - airplane_bounds[::2])
    )

    # If a ground plane z coordinate was provided, compute the ground plane surface
    # once. It is static and will be re-added to the Plotter each frame.
    ground_plane_surface: pv.PolyData | None = None
    ground_plane_texture: pv.Texture | None = None
    if ground_plane_z_E is not None:
        ground_plane_surface = _get_ground_plane_surface(
            ground_plane_z_E,
            trajectory_extent,
            trajectoryMidpoint_E_E,
            airplane_diagonal,
        )
        ground_plane_texture = _get_ground_plane_texture()

    # Add some padding to ensure we can see the whole airplane and wake at each end.
    # Use the airplane's bounding diagonal as the minimum scale reference.
    padding = float(max(airplane_diagonal * 2.0, trajectory_extent * 0.5))
    camera_distance = trajectory_extent + padding

    # Position the camera along the direction (-1, -1, 1) (in PyVista axes) from the
    # trajectory midpoint.
    cameraDirection_V = np.array([-1.0, -1.0, 1.0])
    cameraDirection_V = cameraDirection_V / np.linalg.norm(cameraDirection_V)

    # Transform from Earth axes to PyVista axes. Earth axes have +Z pointing down, but
    # PyVista expects +Z pointing up. A 180 degree rotation about the Y axis flips both
    # X and Z, giving a right handed system with Z up.
    T_pas_E_to_V = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )
    focalPoint_V_E = _transformations.apply_T_to_vectors(
        T_pas_E_to_V,
        trajectoryMidpoint_E_E,
        has_point=True,
    )

    cameraPosition_V_E = focalPoint_V_E + camera_distance * cameraDirection_V
    viewUp_V = (0.0, 0.0, 1.0)
    cpos = [tuple(cameraPosition_V_E), tuple(focalPoint_V_E), viewUp_V]

    # For parallel projection, set the parallel_scale to control the viewport size.
    parallel_scale = camera_distance * 0.6

    # To compute the correct camera clipping range, we need to show PyVista the full
    # extent of the geometry (first and last frames). Add meshes at both positions,
    # set the camera, compute the clipping range, then clear and restart with just
    # the first frame.
    first_panel_surfaces = _get_panel_surfaces_free_flight(
        airplanes[0],
        coupled_solver.stackPosition_E_E[0],
        coupled_solver.stackR_pas_E_to_BP1[0],
    )
    last_step = len(airplanes) - 1
    last_panel_surfaces = _get_panel_surfaces_free_flight(
        airplanes[last_step],
        coupled_solver.stackPosition_E_E[last_step],
        coupled_solver.stackR_pas_E_to_BP1[last_step],
    )

    # Add both first and last frame meshes to compute clipping range.
    plotter.add_mesh(first_panel_surfaces, show_edges=True, color=_panel_color)
    plotter.add_mesh(last_panel_surfaces, show_edges=True, color=_panel_color)

    # Add body mesh at first and last positions to extend the geometry bounds.
    if body_mesh_BP1_CgP1 is not None:
        first_body_surface = _get_body_mesh_surface_free_flight(
            body_mesh_BP1_CgP1,
            coupled_solver.stackPosition_E_E[0],
            coupled_solver.stackR_pas_E_to_BP1[0],
        )
        last_body_surface = _get_body_mesh_surface_free_flight(
            body_mesh_BP1_CgP1,
            coupled_solver.stackPosition_E_E[last_step],
            coupled_solver.stackR_pas_E_to_BP1[last_step],
        )
        plotter.add_mesh(first_body_surface, color=_body_mesh_color)
        plotter.add_mesh(last_body_surface, color=_body_mesh_color)

    # Add wake at last frame if showing wake (this extends the geometry bounds).
    if show_wake_vortices:
        wake_surfaces = _get_wake_ring_vortex_surfaces_free_flight(
            coupled_solver, last_step
        )
        if wake_surfaces.n_points > 0:
            plotter.add_mesh(wake_surfaces, show_edges=True, color=_wake_vortex_color)

    # Set camera and compute clipping range.
    plotter.camera.position = cpos[0]
    plotter.camera.focal_point = cpos[1]
    plotter.camera.up = cpos[2]
    plotter.camera.parallel_scale = parallel_scale
    plotter.reset_camera_clipping_range()
    stored_clipping_range = plotter.camera.clipping_range

    # Clear the plotter and set up the first frame properly.
    plotter.clear()

    # Get the Panel surfaces of the first time step's Airplane in PyVista axes,
    # relative to the Earth origin.
    panel_surfaces = _get_panel_surfaces_free_flight(
        airplanes[0],
        coupled_solver.stackPosition_E_E[0],
        coupled_solver.stackR_pas_E_to_BP1[0],
    )

    # Plot the first time step's Airplanes' Panels either with scalar coloring or
    # with a uniform color.
    if scalar_type is not None:
        these_scalars = _get_scalars(
            (airplanes[0],),
            scalar_type,
            coupled_solver.coupled_steady_problems[0].coupled_operating_point.qInf__E,
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

    # Add the body mesh for the first frame.
    if body_mesh_BP1_CgP1 is not None:
        body_surface = _get_body_mesh_surface_free_flight(
            body_mesh_BP1_CgP1,
            coupled_solver.stackPosition_E_E[0],
            coupled_solver.stackR_pas_E_to_BP1[0],
        )
        plotter.add_mesh(
            body_surface,
            color=_body_mesh_color,
            smooth_shading=True,
        )

    # Add the ground plane for the first frame.
    if ground_plane_surface is not None:
        plotter.add_mesh(
            ground_plane_surface,
            texture=ground_plane_texture,
            opacity=_ground_plane_opacity,
            smooth_shading=False,
        )

    # Set up lighting and shadows when a ground plane is present.
    if ground_plane_surface is not None:
        _setup_free_flight_lighting(plotter)

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
            title="Free Flight Animation - Rendering speed not to scale.",
            cpos=cpos,
            full_screen=False,
            auto_close=False,
        )
    else:
        plotter.show(
            title="Free Flight Animation - Rendering speed not to scale.",
            cpos=cpos,
            full_screen=False,
            interactive=False,
            auto_close=False,
        )

    # Apply the pre-computed clipping range to ensure all geometry throughout the
    # animation is visible. The clipping range defines the near and far planes -
    # objects outside this range won't be rendered.
    plotter.camera.position = cpos[0]
    plotter.camera.focal_point = cpos[1]
    plotter.camera.up = cpos[2]
    plotter.camera.parallel_scale = parallel_scale
    plotter.camera.clipping_range = stored_clipping_range
    time.sleep(1)

    # Start a list to hold a WebP Image of each frame.
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

    # Begin to iterate through the Airplanes from the subsequent time steps.
    for airplane in airplanes[1:]:

        # Clear the Plotter.
        plotter.clear()

        # Get the Panel surfaces of this time step's Airplane in Earth axes,
        # relative to the Earth origin.
        panel_surfaces = _get_panel_surfaces_free_flight(
            airplane,
            coupled_solver.stackPosition_E_E[current_step],
            coupled_solver.stackR_pas_E_to_BP1[current_step],
        )

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
            wake_ring_vortex_surfaces = _get_wake_ring_vortex_surfaces_free_flight(
                coupled_solver, current_step
            )
            if wake_ring_vortex_surfaces.n_points > 0:
                plotter.add_mesh(
                    wake_ring_vortex_surfaces,
                    show_edges=True,
                    smooth_shading=False,
                    color=_wake_vortex_color,
                )

        # Plot the Panels either with a uniform color or with scalar coloring.
        if scalar_type is not None:
            these_scalars = _get_scalars(
                (airplane,),
                scalar_type,
                coupled_solver.coupled_steady_problems[
                    current_step
                ].coupled_operating_point.qInf__E,
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

        # Add the body mesh for this frame.
        if body_mesh_BP1_CgP1 is not None:
            body_surface = _get_body_mesh_surface_free_flight(
                body_mesh_BP1_CgP1,
                coupled_solver.stackPosition_E_E[current_step],
                coupled_solver.stackR_pas_E_to_BP1[current_step],
            )
            plotter.add_mesh(
                body_surface,
                color=_body_mesh_color,
                smooth_shading=True,
            )

        # Add the ground plane for this frame.
        if ground_plane_surface is not None:
            plotter.add_mesh(
                ground_plane_surface,
                texture=ground_plane_texture,
                opacity=_ground_plane_opacity,
                smooth_shading=False,
            )

        # Re-enable lighting and shadows after clearing the Plotter.
        if ground_plane_surface is not None:
            _setup_free_flight_lighting(plotter)

        # If saving, append a WebP Image of this frame to the list of Images.
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
            images,
            "AnimateFreeFlight.webp",
            fps=actual_fps,
            lossless=False,
            quality=_quality,
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
    """Plots an UnsteadyRingVortexLatticeMethodSolver's loads and load coefficients as a
    function of time.

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver whose loads and
        load coefficients will be plotted.
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
    """Prints a solver's load and load coefficients.

    :param solver: The solver whose load and load coefficients will be printed.
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
        col2 = [f"{val:.4e}" for val in col2]
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
        col4 = [f"{val:.4e}" for val in col4]
        col4 = [
            val + " N" if i < 3 else val + " Nm" if i < 6 else val
            for i, val in enumerate(col4)
        ]

        print(f'Airplane "{airplane.name}":')

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
            print(f"{pad}Reynolds Number: {re:.2e}")

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
def _get_panel_surfaces_free_flight(
    airplane: geometry.airplane.Airplane,
    position_E_E: np.ndarray,
    R_pas_E_to_BP1: np.ndarray,
) -> pv.PolyData:
    """Returns a PolyData representation of the Wings' Panels' surfaces, in PyVista
    axes, relative to the Earth origin, for free flight visualization.

    :param airplane: The Airplane whose Wings' Panels' surfaces will be returned.
    :param position_E_E: A (3,) ndarray of floats representing the current position of
        the Airplane's CG (in Earth axes, relative to the Earth origin). The units are
        in meters.
    :param R_pas_E_to_BP1: A (3,3) ndarray of floats representing the current
        orientation of the Airplane, expressed as a passive rotation matrix from Earth
        axes to the Airplane's body axes.
    :return: The PolyData representation of the Airplane's Wings' Panels' surfaces (in
        PyVista axes, relative to the Earth origin).
    """
    # Initialize empty ndarrays to hold the Panels' vertices and faces.
    stackVertices_GP1_CgP1 = np.empty((0, 3), dtype=float)
    faces = np.empty(0, dtype=int)

    # Initialize a variable to keep track of how many Panels have been added thus far.
    panel_num = 0

    # Increment through the Airplane's Wing(s).
    for wing in airplane.wings:
        _panels = wing.panels
        assert _panels is not None

        # Unravel this Wing's ndarray of Panels iterate through it.
        panels = np.ravel(_panels)
        for panel in panels:
            # Arrange this Panel's vertices and faces into ndarrays in the
            # proper form to represent PolyData surfaces.
            vertices_GP1_CgP1 = np.vstack(
                (
                    panel.Flpp_GP1_CgP1,
                    panel.Frpp_GP1_CgP1,
                    panel.Brpp_GP1_CgP1,
                    panel.Blpp_GP1_CgP1,
                )
            )
            face_to_add = np.array(
                [
                    4,
                    (panel_num * 4),
                    (panel_num * 4) + 1,
                    (panel_num * 4) + 2,
                    (panel_num * 4) + 3,
                ],
                dtype=int,
            )

            # Add this Panel's vertices and faces to the ndarray of all vertices and
            # faces.
            stackVertices_GP1_CgP1 = np.vstack(
                (stackVertices_GP1_CgP1, vertices_GP1_CgP1)
            )
            faces = np.hstack((faces, face_to_add))

            # Update the number of Panels.
            panel_num += 1

    # Geometry axes point opposite to body axes (rotated 180 degrees about y-axis).
    # Create the passive rotation transformation from the first Airplane's geometry
    # axes, relative to the first Airplane's CG, to the first Airplane's body axes,
    # relative to the first Airplane's CG.
    T_pas_GP1_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )

    # Create the passive rotation transformation from Earth axes, relative to the
    # first Airplane's CG, to the first Airplane's body axes, relative to the first
    # Airplane's CG. R_pas_E_to_BP1 is a (3, 3) matrix, so we need to embed it in a (
    # 4, 4) homogeneous transformation matrix.
    T_pas_E_CgP1_to_BP1_CgP1 = np.eye(4, dtype=float)
    T_pas_E_CgP1_to_BP1_CgP1[:3, :3] = R_pas_E_to_BP1

    # Invert to get the transformation from the first Airplane's body axes, relative
    # to the first Airplane's CG, to Earth axes, relative to the first Airplane's CG.
    T_pas_BP1_CgP1_to_E_CgP1 = _transformations.invert_T_pas(T_pas_E_CgP1_to_BP1_CgP1)

    # Create the passive translation transformation from the Earth axes, relative to
    # the first Airplane's CG, to Earth axes, relative to the Earth origin.
    #
    # For generate_trans_T with passive=True, the translations parameter should be the
    # position of the final reference point (the Earth origin) relative to the initial
    # reference point (the first Airplane's CG). Since position_E_E is the position of
    # the first Airplane's CG relative to Earth origin, the position of the Earth origin
    # relative to the first Airplane's CG is -position_E_E.
    T_pas_E_CgP1_to_E_E = _transformations.generate_trans_T(
        translations=-position_E_E,
        passive=True,
    )

    # Compose the transformations: GP1_CgP1 -> BP1_CgP1 -> E_CgP1 -> E_E.
    # Note: For passive transformations, compose_T_pas takes them in path order.
    T_pas_GP1_CgP1_to_E_E = _transformations.compose_T_pas(
        T_pas_GP1_CgP1_to_BP1_CgP1,
        T_pas_BP1_CgP1_to_E_CgP1,
        T_pas_E_CgP1_to_E_E,
    )

    # Apply transformation to all vertices.
    stackVertices_E_E = _transformations.apply_T_to_vectors(
        T_pas_GP1_CgP1_to_E_E,
        stackVertices_GP1_CgP1,
        has_point=True,
    )

    # REFACTOR: Add section to AXES_POINTS_AND_FRAMES.md about PyVista axes.
    # Transform from Earth axes to PyVista axes. Earth axes have +Z pointing down, but
    # PyVista expects +Z pointing up. A 180 degree rotation about the Y axis flips both
    # X and Z, giving a right handed system with Z up.
    T_pas_E_to_V = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )
    stackVertices_V_E = _transformations.apply_T_to_vectors(
        T_pas_E_to_V,
        stackVertices_E_E,
        has_point=True,
    )

    # Return the Panels' surfaces.
    return pv.PolyData(stackVertices_V_E, faces)


# TEST: Consider adding unit tests for this function.
def _get_body_mesh_surface_free_flight(
    body_mesh_BP1_CgP1: pv.PolyData,
    position_E_E: np.ndarray,
    R_pas_E_to_BP1: np.ndarray,
) -> pv.PolyData:
    """Returns a PolyData representation of the body STL mesh, in PyVista axes, relative
    to the Earth origin, for free flight visualization.

    The body mesh is assumed to already be in the first Airplane's body axes, relative
    to the first Airplane's CG. The transformation chain is: BP1_CgP1 -> E_CgP1 -> E_E
    -> V_E.

    :param body_mesh_BP1_CgP1: A PolyData representation of the body mesh (in the first
        Airplane's body axes, relative to the first Airplane's CG).
    :param position_E_E: A (3,) ndarray of floats representing the current position of
        the Airplane's CG (in Earth axes, relative to the Earth origin). The units are
        in meters.
    :param R_pas_E_to_BP1: A (3,3) ndarray of floats representing the current
        orientation of the Airplane, expressed as a passive rotation matrix from Earth
        axes to the Airplane's body axes.
    :return: A PolyData representation of the body mesh (in PyVista axes, relative to
        the Earth origin).
    """
    # Get the body mesh's vertices in body axes, relative to the first Airplane's CG.
    stackVertices_BP1_CgP1 = np.array(body_mesh_BP1_CgP1.points, dtype=float)

    # Create the passive rotation transformation from Earth axes, relative to the
    # first Airplane's CG, to the first Airplane's body axes, relative to the first
    # Airplane's CG. R_pas_E_to_BP1 is a (3, 3) matrix, so we need to embed it in a
    # (4, 4) homogeneous transformation matrix.
    T_pas_E_CgP1_to_BP1_CgP1 = np.eye(4, dtype=float)
    T_pas_E_CgP1_to_BP1_CgP1[:3, :3] = R_pas_E_to_BP1

    # Invert to get the transformation from the first Airplane's body axes, relative
    # to the first Airplane's CG, to Earth axes, relative to the first Airplane's CG.
    T_pas_BP1_CgP1_to_E_CgP1 = _transformations.invert_T_pas(T_pas_E_CgP1_to_BP1_CgP1)

    # Create the passive translation transformation from Earth axes, relative to
    # the first Airplane's CG, to Earth axes, relative to the Earth origin.
    T_pas_E_CgP1_to_E_E = _transformations.generate_trans_T(
        translations=-position_E_E,
        passive=True,
    )

    # Compose the transformations: BP1_CgP1 -> E_CgP1 -> E_E.
    T_pas_BP1_CgP1_to_E_E = _transformations.compose_T_pas(
        T_pas_BP1_CgP1_to_E_CgP1,
        T_pas_E_CgP1_to_E_E,
    )

    # Apply transformation to all vertices.
    stackVertices_E_E = _transformations.apply_T_to_vectors(
        T_pas_BP1_CgP1_to_E_E,
        stackVertices_BP1_CgP1,
        has_point=True,
    )

    # Transform from Earth axes to PyVista axes. Earth axes have +Z pointing down, but
    # PyVista expects +Z pointing up. A 180 degree rotation about the Y axis flips both
    # X and Z, giving a right handed system with Z up.
    T_pas_E_to_V = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )
    stackVertices_V_E = _transformations.apply_T_to_vectors(
        T_pas_E_to_V,
        stackVertices_E_E,
        has_point=True,
    )

    # Return a new PolyData with the transformed points and the original faces.
    return pv.PolyData(stackVertices_V_E, body_mesh_BP1_CgP1.faces)


# TEST: Consider adding unit tests for this function.
def _get_ground_plane_surface(
    ground_plane_z_E: float,
    trajectory_extent: float,
    trajectoryMidpoint_E_E: np.ndarray,
    airplane_diagonal: float,
) -> pv.PolyData:
    """Returns a PolyData representation of a ground plane, in PyVista axes, relative to
    the Earth origin.

    The plane is centered at the trajectory's x/y midpoint in Earth axes, at z =
    ground_plane_z_E. It is sized to always extend beyond the trajectory. The returned
    plane has its texture coordinates scaled so that a checkerboard texture tiles
    proportionally to the airplane's size.

    :param ground_plane_z_E: The z coordinate of the ground plane (in Earth axes). The
        units are in meters.
    :param trajectory_extent: The distance between the initial and final positions of
        the trajectory. The units are in meters.
    :param trajectoryMidpoint_E_E: A (3,) ndarray of floats representing the midpoint of
        the trajectory (in Earth axes, relative to the Earth origin). The units are in
        meters.
    :param airplane_diagonal: The bounding diagonal of the airplane's panel surfaces (in
        PyVista axes). Used to scale the checkerboard texture tiles. The units are in
        meters.
    :return: A PolyData representation of the ground plane (in PyVista axes, relative to
        the Earth origin).
    """
    # Compute the plane size to extend beyond the trajectory. Use the airplane
    # diagonal as the minimum scale reference.
    plane_size = max(airplane_diagonal * 8.0, trajectory_extent * 2.0)

    # Build the plane center in Earth axes, using the trajectory's x/y midpoint and
    # the specified ground z coordinate.
    planeCenter_E_E = np.array(
        [trajectoryMidpoint_E_E[0], trajectoryMidpoint_E_E[1], ground_plane_z_E],
        dtype=float,
    )

    # Transform the center from Earth axes to PyVista axes. Earth axes have +Z
    # pointing down, but PyVista expects +Z pointing up. A 180 degree rotation about
    # the Y axis flips both X and Z, giving a right handed system with Z up.
    T_pas_E_to_V = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )
    planeCenter_V_E = _transformations.apply_T_to_vectors(
        T_pas_E_to_V,
        planeCenter_E_E,
        has_point=True,
    )

    # Create the ground plane. The ground plane is perpendicular to Earth's z axis.
    # Under the 180 degree Y rotation (E to V), Earth's z axis (0,0,1) maps to
    # PyVista's (0,0,-1). The ground plane normal should point "up" in PyVista, so
    # the direction is (0, 0, 1).
    plane = pv.Plane(
        center=planeCenter_V_E,
        direction=(0.0, 0.0, 1.0),
        i_size=plane_size,
        j_size=plane_size,
    )

    # Scale the texture coordinates so the checkerboard tiles proportionally to the
    # physical size. Each tile covers one unit of plane_size, so the number of tiles
    # equals the plane size divided by the airplane diagonal.
    num_tiles = (
        max(4, int(plane_size / airplane_diagonal)) if airplane_diagonal > 0 else 4
    )
    plane.active_texture_coordinates = plane.active_texture_coordinates * num_tiles

    return plane


def _setup_free_flight_lighting(plotter: pv.Plotter) -> None:
    """Adds a directional light and enables shadow mapping on the given Plotter.

    The light direction is chosen to illuminate the scene from above and to the side,
    producing visible shadows on the ground plane.

    :param plotter: The Plotter to configure.
    :return: None
    """
    light = pv.Light(
        position=(1.0, -1.0, -1.0),
        focal_point=(0.0, 0.0, 0.0),
        color="white",
        intensity=1.0,
    )
    plotter.add_light(light)

    # Add a dim ambient fill light so shadowed regions are not completely black.
    fill = pv.Light(
        position=(0.0, 0.0, -1.0),
        focal_point=(0.0, 0.0, 0.0),
        color="white",
        intensity=0.3,
    )
    plotter.add_light(fill)

    plotter.enable_shadows()


def _get_ground_plane_texture() -> pv.Texture:
    """Returns a checkerboard Texture for the ground plane.

    The texture is a repeating checkerboard pattern. The number of visible tiles depends
    on the texture coordinate scaling applied to the ground plane surface by
    ``_get_ground_plane_surface``.

    :return: A PyVista Texture containing a checkerboard pattern with wrap mode set to
        repeat.
    """
    n = _ground_plane_checker_resolution
    block = _ground_plane_checker_block_size
    checker = np.zeros((n, n, 3), dtype=np.uint8)
    for i in range(n):
        for j in range(n):
            if (i // block + j // block) % 2 == 0:
                checker[i, j] = _ground_plane_checker_light
            else:
                checker[i, j] = _ground_plane_checker_dark
    tex = pv.numpy_to_texture(checker)
    tex.SetRepeat(True)
    return tex


# TEST: Consider adding unit tests for this function.
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
    reflected = mesh.copy()
    reflected.points = _transformations.apply_T_to_vectors(
        T_reflect,
        mesh.points,
        has_point=True,
    )
    return reflected


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


# TEST: Consider adding unit tests for this function.
def _get_wake_ring_vortex_surfaces_free_flight(
    coupled_solver: coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
    step: int,
) -> pv.PolyData:
    """Returns the PolyData representation of surfaces a
    CoupledUnsteadyRingVortexLatticeMethodSolver's wake RingVortices at a given time
    step, in PyVista axes, relative to the Earth origin.

    :param coupled_solver: The CoupledUnsteadyRingVortexLatticeMethodSolver with the
        wake RingVortices to process.
    :param step: The time step number at which to process the wake RingVortices.
    :return: The PolyData representation of the wake RingVortices, in PyVista axes,
        relative to the Earth origin.
    """
    position_E_E = coupled_solver.stackPosition_E_E[step]
    R_pas_E_to_BP1 = coupled_solver.stackR_pas_E_to_BP1[step]

    num_wake_ring_vortices = coupled_solver.list_num_wake_vortices[step]
    stackFrwrvp_GP1_CgP1 = coupled_solver.listStackFrwrvp_GP1_CgP1[step]
    stackFlwrvp_GP1_CgP1 = coupled_solver.listStackFlwrvp_GP1_CgP1[step]
    stackBlwrvp_GP1_CgP1 = coupled_solver.listStackBlwrvp_GP1_CgP1[step]
    stackBrwrvp_GP1_CgP1 = coupled_solver.listStackBrwrvp_GP1_CgP1[step]

    # Initialize empty ndarrays to hold each wake RingVortex's vertices and face.
    stackVertices_GP1_CgP1 = np.zeros((0, 3), dtype=float)
    faces = np.zeros(0, dtype=int)

    for wake_ring_vortex_num in range(num_wake_ring_vortices):
        Frwrvp_GP1_CgP1 = stackFrwrvp_GP1_CgP1[wake_ring_vortex_num]
        Flwrvp_GP1_CgP1 = stackFlwrvp_GP1_CgP1[wake_ring_vortex_num]
        Blwrvp_GP1_CgP1 = stackBlwrvp_GP1_CgP1[wake_ring_vortex_num]
        Brwrvp_GP1_CgP1 = stackBrwrvp_GP1_CgP1[wake_ring_vortex_num]

        newVertices_GP1_CgP1 = np.vstack(
            (
                Flwrvp_GP1_CgP1,
                Frwrvp_GP1_CgP1,
                Brwrvp_GP1_CgP1,
                Blwrvp_GP1_CgP1,
            )
        )
        new_faces = np.array(
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
        stackVertices_GP1_CgP1 = np.vstack(
            (stackVertices_GP1_CgP1, newVertices_GP1_CgP1)
        )
        faces = np.hstack((faces, new_faces))

    # Geometry axes point opposite to body axes (rotated 180 degrees about y-axis).
    # Create the passive rotation transformation from the first Airplane's geometry
    # axes, relative to the first Airplane's CG, to the first Airplane's body axes,
    # relative to the first Airplane's CG.
    T_pas_GP1_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )

    # Create the passive rotation transformation from Earth axes, relative to the
    # first Airplane's CG, to the first Airplane's body axes, relative to the first
    # Airplane's CG. R_pas_E_to_BP1 is a (3, 3) matrix, so we need to embed it in a (
    # 4, 4) homogeneous transformation matrix.
    T_pas_E_CgP1_to_BP1_CgP1 = np.eye(4, dtype=float)
    T_pas_E_CgP1_to_BP1_CgP1[:3, :3] = R_pas_E_to_BP1

    # Invert to get the transformation from the first Airplane's body axes, relative
    # to the first Airplane's CG, to Earth axes, relative to the first Airplane's CG.
    T_pas_BP1_CgP1_to_E_CgP1 = _transformations.invert_T_pas(T_pas_E_CgP1_to_BP1_CgP1)

    # Create the passive translation transformation from the Earth axes, relative to
    # the first Airplane's CG, to Earth axes, relative to the Earth origin.
    #
    # For generate_trans_T with passive=True, the translations parameter should be the
    # position of the final reference point (the Earth origin) relative to the initial
    # reference point (the first Airplane's CG). Since position_E_E is the position of
    # the first Airplane's CG relative to the Earth origin, the position of the Earth
    # origin relative to the first Airplane's CG is -position_E_E.
    T_pas_E_CgP1_to_E_E = _transformations.generate_trans_T(
        translations=-position_E_E,
        passive=True,
    )

    # Compose the transformations: GP1_CgP1 -> BP1_CgP1 -> E_CgP1 -> E_E.
    # Note: For passive transformations, compose_T_pas takes them in path order.
    T_pas_GP1_CgP1_to_E_E = _transformations.compose_T_pas(
        T_pas_GP1_CgP1_to_BP1_CgP1,
        T_pas_BP1_CgP1_to_E_CgP1,
        T_pas_E_CgP1_to_E_E,
    )

    # Apply transformation to all vertices.
    stackVertices_E_E = _transformations.apply_T_to_vectors(
        T_pas_GP1_CgP1_to_E_E,
        stackVertices_GP1_CgP1,
        has_point=True,
    )

    # Transform from Earth axes to PyVista's viewing convention. Earth axes have +Z
    # pointing down, but PyVista expects +Z pointing up. A 180-degree rotation about
    # the Y-axis flips both X and Z, giving a right-handed system with Z up.
    T_pas_E_to_V = _transformations.generate_rot_T(
        angles=(0.0, 180.0, 0.0),
        passive=True,
        intrinsic=True,
        order="xyz",
    )
    stackVertices_V_E = _transformations.apply_T_to_vectors(
        T_pas_E_to_V,
        stackVertices_E_E,
        has_point=True,
    )

    # Return the wake RingVortex surfaces.
    return pv.PolyData(stackVertices_V_E, faces)


# TEST: Consider adding unit tests for this function.
def _get_wake_ring_vortex_surfaces(
    solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    step: int,
) -> pv.PolyData:
    """Returns the PolyData representation of the surfaces of an
    UnsteadyRingVortexLatticeMethodSolver's wake RingVortices at a given time step.

    :param solver: The UnsteadyRingVortexLatticeMethodSolver with the wake RingVortices
        to process.
    :param step: The time step at which to process the wake RingVortices.
    :return: The PolyData representation of the wake RingVortices.
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
        wake_ring_vortex_vertices = np.vstack(
            (wake_ring_vortex_vertices, wake_ring_vortex_vertices_to_add)
        )
        wake_ring_vortex_faces = np.hstack(
            (wake_ring_vortex_faces, wake_ring_vortex_face_to_add)
        )

    # Return the wake RingVortex surfaces.
    return pv.PolyData(wake_ring_vortex_vertices, wake_ring_vortex_faces)


# TEST: Consider adding unit tests for this function.
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
