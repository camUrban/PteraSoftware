"""Contains functions for visualizing geometry and results."""

from __future__ import annotations

import math
import os.path
import time
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import matplotlib.colors
import matplotlib.legend_handler
import matplotlib.pyplot as plt
import matplotlib.typing
import numpy as np
import pyvista as pv
import webp
from tqdm import tqdm

from . import (
    _free_flight_unsteady_ring_vortex_lattice_method,
    _logging,
    _mujoco_model,
    _output_plotting,
    _output_rendering,
    _parameter_validation,
    _steady_horseshoe_vortex_lattice_method,
    _steady_ring_vortex_lattice_method,
    _transformations,
    _unsteady_ring_vortex_lattice_method,
)
from ._movements import free_flight_movement as free_flight_movement_mod

logger = _logging.get_logger("output")

# Define the Plotter's appearance. The streamline line width is in pixels and is tuned
# for _output_rendering.REFERENCE_WINDOW_SIZE, so it is scaled wherever it is used, as
# the font sizes in _output_rendering are.
STREAMLINE_COLOR = "orchid"
STREAMLINE_LINE_WIDTH = 2.0
IMAGE_SURFACE_OPACITY = 0.5
ANIMATE_PREVIEW_FIRST_OPACITY = 0.10
ANIMATE_PREVIEW_LAST_OPACITY = 0.35

# Define the number of samples used for multisample anti-aliasing. PyVista defaults to
# 8, whose resolve is not reproducible on every driver: rendering one scene twice can
# differ by a few intensity levels along an anti-aliased edge, which makes a saved WebP
# vary between runs. Four samples renders indistinguishably and varied less than 8 in
# testing, so the visualizations pin it rather than take the default. It is not fully
# reproducible either, since a pixel along an edge can still differ by one intensity
# level between runs.
MULTI_SAMPLES = 4

# Define the colors of the series in the results plots, which are named colors from
# Matplotlib's xkcd color survey table.
ALPHA_COLOR = "xkcd:peacock blue"
BETA_COLOR = "xkcd:sea"
LINEAR_X_COLOR = "xkcd:dark sea green"
LINEAR_Y_COLOR = "xkcd:nasty green"
LINEAR_Z_COLOR = "xkcd:squash"
ANGULAR_X_COLOR = "xkcd:pumpkin"
ANGULAR_Y_COLOR = "xkcd:dark coral"
ANGULAR_Z_COLOR = "xkcd:red violet"

# Define the text that the results outputs share. Every figure's legend labels,
# subtitle, and y axis label are named once here because the other two outputs restate
# them: a CSV header is the transformed form of the same three pieces, taking its
# quantity from the legend label, its axes, point, and frame from the subtitle, and its
# unit from the y axis label, while a logged group header pairs a quantity with the same
# subtitle. Naming them once is what keeps the three describing a quantity the same way.
FORCE_LABELS = ["Induced Drag", "Crosswind Force", "Lift"]
FORCE_COEFFICIENT_LABELS = [
    "Induced Drag Coefficient",
    "Crosswind Force Coefficient",
    "Lift Coefficient",
]
MOMENT_LABELS = ["Rolling Moment", "Pitching Moment", "Yawing Moment"]
MOMENT_COEFFICIENT_LABELS = [
    "Rolling Moment Coefficient",
    "Pitching Moment Coefficient",
    "Yawing Moment Coefficient",
]

# The position and velocity figures label their series by component alone, since their
# titles name the quantity. A CSV column has no title, so it takes the quantity from the
# y axis label instead.
COMPONENT_LABELS = ["X Component", "Y Component", "Z Component"]
ORIENTATION_LABELS = ["Roll Angle", "Pitch Angle", "Yaw Angle"]
ANGULAR_VELOCITY_LABELS = ["Roll Rate", "Pitch Rate", "Yaw Rate"]
AERODYNAMIC_ANGLE_LABELS = ["Angle of Attack", "Sideslip Angle"]

WIND_AXES_SUBTITLE = "(in Wind Axes)"
WIND_AXES_CG_SUBTITLE = "(in Wind Axes, Relative to the CG)"

# The logged results report the loads in each Airplane's own geometry axes as well,
# which no figure plots.
GEOMETRY_AXES_SUBTITLE = "(in Geometry Axes)"
GEOMETRY_AXES_CG_SUBTITLE = "(in Geometry Axes, Relative to the CG)"
POSITION_SUBTITLE = (
    "(of the First Airplane's CG, in Earth Axes, Relative to the Earth Origin)"
)
VELOCITY_SUBTITLE = (
    "(of the First Airplane's CG, in Earth Axes, Observed from the Earth Frame)"
)
ORIENTATION_SUBTITLE = (
    "(of the First Airplane's Body Axes Relative to Earth Axes Using an Intrinsic "
    "zy'x\" Sequence)"
)
ANGULAR_VELOCITY_SUBTITLE = (
    "(in the First Airplane's Body Axes, Observed from the Earth Frame)"
)
AERODYNAMIC_ANGLE_SUBTITLE = ""

FORCE_Y_LABEL = "Force (N)"
FORCE_COEFFICIENT_Y_LABEL = "Force Coefficient"
MOMENT_Y_LABEL = "Moment (N m)"
MOMENT_COEFFICIENT_Y_LABEL = "Moment Coefficient"
POSITION_Y_LABEL = "Position (m)"
VELOCITY_Y_LABEL = "Velocity (m/s)"
ORIENTATION_Y_LABEL = "Orientation (deg)"
ANGULAR_VELOCITY_Y_LABEL = "Angular Velocity (deg/s)"
AERODYNAMIC_ANGLE_Y_LABEL = "Angle (deg)"

# Define the camera's view-up direction for free flight visualizations. Earth axes have
# +z pointing down, so physical up is the -z direction. The free flight visualizations
# render geometry in Earth axes (so the body flies through the scene in its true pose)
# and use this view-up so that down appears downward on screen. This is a rendering
# setting, not an axis system.
freeFlightViewUp_E = np.array([0.0, 0.0, -1.0], dtype=float)

# Define the camera's offset direction for free flight visualizations, which is the unit
# vector from the focal point toward the camera position (in Earth axes). The camera
# looks along the negative of this vector. This views the scene obliquely from the
# North, West, and above (Earth +x, -y, and -z). The normalization is computed in plain
# Python because the documentation build imports this module with NumPy replaced by a
# stand-in that does not support arithmetic.
freeFlightCameraOffsetDirection_E = np.array(
    [component / math.sqrt(3.0) for component in (1.0, -1.0, -1.0)], dtype=float
)


def set_preview_opacity(actors: list[pv.Actor], opacity: float) -> None:
    """Sets the opacity of temporary animation-preview actors.

    :param actors: The actors whose opacity should be changed.
    :param opacity: The opacity to apply to the actors.
    :return: None.
    """
    for actor in actors:
        if isinstance(actor, pv.Actor):
            actor.prop.opacity = opacity


def create_plotter(
    window_width: int, window_height: int, off_screen: bool
) -> pv.Plotter:
    """Creates a Plotter with a parallel projection, multisample anti-aliasing, and the
    background color.

    The background color is set before anything realizes the window, so that the window
    appears in its final color rather than flashing white first.

    :param window_width: The width, in pixels, of the render window.
    :param window_height: The height, in pixels, of the render window.
    :param off_screen: Set this to True to render off screen, without showing a window.
    :return: The new Plotter.
    """
    plotter = pv.Plotter(
        window_size=[window_width, window_height],
        off_screen=off_screen,
        lighting=None,
    )
    plotter.enable_parallel_projection()  # type: ignore[call-arg]
    plotter.enable_anti_aliasing("msaa", multi_samples=MULTI_SAMPLES)
    plotter.set_background(  # type: ignore[call-arg]
        color=_output_rendering.PLOTTER_BACKGROUND_COLOR
    )
    return plotter


def create_preview_plotter(
    window_width: int, window_height: int
) -> tuple[pv.Plotter, float]:
    """Creates the Plotter whose window shows a drawing or an animation's preview,
    scaled down when necessary to fit on screen with the requested aspect ratio.

    A window manager will not grant an on-screen render window larger than the display
    less any docks or bars and less the window's own title bar, and VTK silently shrinks
    one that asks for more. A window that was shrunk is replaced with one scaled down
    uniformly to fit within the size that was granted, so it keeps the requested aspect
    ratio and frames the scene the same way a saved file of the requested size does. The
    sizes compared are all ones the render window reports about itself, which accounts
    for the title bar and the other window decorations without measuring them. Only a
    granted size smaller than the request counts as a shrink, since a display that
    scales its pixels reports a larger one. Off screen, the window is never shrunk.

    :param window_width: The requested width, in pixels, of the render window.
    :param window_height: The requested height, in pixels, of the render window.
    :return: A tuple of the new Plotter and the factor by which to scale the font sizes
        and line widths in its window, as returned by _output_rendering.get_window_scale
        for the window's size.
    """
    plotter = create_plotter(window_width, window_height, pv.OFF_SCREEN)
    if pv.OFF_SCREEN:
        return plotter, _output_rendering.get_window_scale(window_width, window_height)

    # Rendering the empty scene realizes the window, which is what makes the granted
    # size readable.
    render_window = plotter.ren_win
    assert render_window is not None
    render_window.Render()
    render_window.SetWindowName("Assembling the scene. Please wait.")
    render_window.Render()
    granted_width, granted_height = render_window.GetSize()
    fit = min(granted_width / window_width, granted_height / window_height)
    if fit >= 1.0:
        return plotter, _output_rendering.get_window_scale(window_width, window_height)

    # The window was shrunk, so replace it with one scaled down uniformly to fit within
    # the granted size. Its outer frame is then no larger than the one the window
    # manager just granted.
    plotter.close()
    fitted_width = max(1, int(window_width * fit))
    fitted_height = max(1, int(window_height * fit))
    plotter = create_plotter(fitted_width, fitted_height, False)
    render_window = plotter.ren_win
    assert render_window is not None
    render_window.Render()
    render_window.SetWindowName("Assembling the scene. Please wait.")
    render_window.Render()
    granted_width, granted_height = render_window.GetSize()
    if granted_width < fitted_width or granted_height < fitted_height:
        plotter.close()
        raise ValueError(
            f"window_size {window_width} by {window_height} does not fit on screen, "
            f"and the window manager granted only {granted_width} by "
            f"{granted_height} pixels to a window scaled down to {fitted_width} by "
            f"{fitted_height} pixels to fit. Request a smaller window, or render off "
            f"screen by setting pyvista.OFF_SCREEN to True."
        )
    return plotter, _output_rendering.get_window_scale(fitted_width, fitted_height)


def add_scene(
    plotter: pv.Plotter,
    panel_surfaces: pv.PolyData,
    wake_ring_vortex_surfaces: pv.PolyData | None,
    coloring: _output_rendering.ScalarColoring | None,
    reflect_T_act: np.ndarray | None,
    window_scale: float,
    streamline_surfaces: pv.PolyData | None,
    image_surface_mesh: pv.PolyData | None,
    image_surface_texture: pv.Texture | None,
    worldbody_geoms: list[_mujoco_model.RenderGeom],
    body_geoms: list[_mujoco_model.RenderGeom],
    T_pas_BP1_CgP1_to_E_Eo: np.ndarray | None,
) -> None:
    """Adds the actors of a drawing's scene, or of one animation frame's scene, to a
    Plotter.

    The actors are added in a fixed order: the wake and the Panels, the streamlines, the
    image surface plane, and then the MuJoCo geometry. The meshes arrive already mapped
    into whichever axes they are being rendered in. This function does not touch the
    camera, and it does not add an animation's text overlays.

    :param plotter: The Plotter to add the scene to.
    :param panel_surfaces: The PolyData representation of the Panel surfaces.
    :param wake_ring_vortex_surfaces: The PolyData representation of the wake ring
        vortex surfaces, or None to omit the wake.
    :param coloring: The scalar coloring to apply to the Panels, or None to color them
        uniformly.
    :param reflect_T_act: A (4,4) ndarray of floats representing the active
        transformation (in whichever axes the meshes are rendered in) that reflects
        geometry across the image surface, or None when no image surface is defined, in
        which case no reflected geometry is added.
    :param window_scale: The factor by which to scale the line widths and font sizes, as
        returned by _output_rendering.get_window_scale.
    :param streamline_surfaces: The PolyData representation of the streamlines, or None
        to omit them.
    :param image_surface_mesh: The image surface plane's mesh, or None to omit it. It
        must not be None when reflect_T_act is not None.
    :param image_surface_texture: The image surface plane's Texture, or None when
        image_surface_mesh is None.
    :param worldbody_geoms: The RenderGeoms attached to the worldbody.
    :param body_geoms: The RenderGeoms attached to the body.
    :param T_pas_BP1_CgP1_to_E_Eo: A (4,4) ndarray of floats representing the passive
        transformation from the first Airplane's body axes, relative to the first
        Airplane's CG, to Earth axes, relative to the Earth origin, at the scene's time
        step, or None to omit the MuJoCo geometry.
    :return: None
    """
    # Add the wake, the Panels, and, if an image surface is defined, their reflections.
    _output_rendering.add_frame_geometry(
        plotter,
        panel_surfaces,
        wake_ring_vortex_surfaces,
        coloring,
        reflect_T_act,
        window_scale,
    )

    # If showing streamlines, plot them.
    if streamline_surfaces is not None:
        plotter.add_mesh(
            streamline_surfaces,
            show_edges=True,
            color=STREAMLINE_COLOR,
            line_width=STREAMLINE_LINE_WIDTH * window_scale,
            smooth_shading=False,
            lighting=False,
            render=False,
        )

        # If an image surface is defined, add the reflected streamlines, muted toward
        # gray so that they read as a reflection rather than as more streamlines.
        if reflect_T_act is not None:
            plotter.add_mesh(
                _output_rendering.transform_mesh(streamline_surfaces, reflect_T_act),
                show_edges=True,
                color=_output_rendering.mute_color(
                    STREAMLINE_COLOR, _output_rendering.IMAGE_REFLECTION_MUTE_FACTOR
                ),
                line_width=STREAMLINE_LINE_WIDTH * window_scale,
                smooth_shading=False,
                lighting=False,
                render=False,
            )

    # If an image surface is defined, add its plane.
    if image_surface_mesh is not None:
        plotter.add_mesh(
            image_surface_mesh,
            texture=image_surface_texture,
            opacity=IMAGE_SURFACE_OPACITY,
            smooth_shading=True,
            lighting=False,
            render=False,
        )

    # If showing MuJoCo geometry, add it posed at the scene's time step.
    if T_pas_BP1_CgP1_to_E_Eo is not None:
        _output_rendering.add_mujoco_geometry(
            plotter,
            worldbody_geoms,
            body_geoms,
            T_pas_BP1_CgP1_to_E_Eo,
            reflect_T_act,
        )


def get_draw_geometry_meshes(
    panel_surfaces: pv.PolyData,
    wake_ring_vortex_surfaces: pv.PolyData | None,
    streamline_surfaces: pv.PolyData | None,
    reflect_T_act: np.ndarray | None,
) -> list[pv.PolyData]:
    """Returns the meshes of draw's geometry, which are the wake, the Panels, and the
    streamlines, along with their reflections when an image surface is defined.

    Empty meshes are left out, since PyVista cannot add them to a Plotter, so the
    returned meshes' bounds match the bounds of the actors draw adds for its geometry.

    :param panel_surfaces: The PolyData representation of the Panel surfaces.
    :param wake_ring_vortex_surfaces: The PolyData representation of the wake ring
        vortex surfaces, or None when the wake is omitted.
    :param streamline_surfaces: The PolyData representation of the streamlines, or None
        when they are omitted.
    :param reflect_T_act: A (4,4) ndarray of floats representing the active
        transformation (in whichever axes the meshes are rendered in) that reflects
        geometry across the image surface, or None when no image surface is defined, in
        which case no reflected meshes are returned.
    :return: A list of the non empty geometry meshes.
    """
    geometry_meshes = [panel_surfaces]
    if wake_ring_vortex_surfaces is not None:
        geometry_meshes.append(wake_ring_vortex_surfaces)
    if streamline_surfaces is not None:
        geometry_meshes.append(streamline_surfaces)
    if reflect_T_act is not None:
        geometry_meshes += [
            _output_rendering.transform_mesh(geometry_mesh, reflect_T_act)
            for geometry_mesh in geometry_meshes
        ]
    return [
        geometry_mesh for geometry_mesh in geometry_meshes if geometry_mesh.n_points > 0
    ]


def draw(
    solver: (
        _steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | _steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    scalar_type: str | None = None,
    show_streamlines: bool | np.bool = False,
    show_wake_vortices: bool | np.bool = False,
    show_mujoco_geometry: bool | np.bool = False,
    window_size: Sequence[int] = (1024, 768),
    save: bool | np.bool = False,
    path: str | Path = "draw.webp",
    quality: int | float = 75.0,
    testing: bool | np.bool = False,
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
        "crosswind force", or "lift", which respectively use each Panel's induced drag,
        crosswind force, and lift coefficient. The default is None.
    :param show_streamlines: Set this to True to show the streamlines emanating from the
        back of the Wings. If True, the solver's streamlines must have already been
        calculated. Can be a bool or a numpy bool and will be converted internally to a
        bool. The default is False.
    :param show_wake_vortices: Set this to True to show any wake ring vortices. If True,
        the solver must be an UnsteadyRingVortexLatticeMethodSolver and must have
        already been run. Can be a bool or a numpy bool and will be converted internally
        to a bool. The default is False.
    :param show_mujoco_geometry: Set this to True to show the MuJoCo geometry that
        extra_xml and mujoco_assets inject into the solver's model, with body geoms
        posed at the drawn time step and worldbody geoms static in Earth axes. If True,
        the solver must be a FreeFlightUnsteadyRingVortexLatticeMethodSolver. Can be a
        bool or a numpy bool and will be converted internally to a bool. The default is
        False.
    :param window_size: The width and height, in pixels, of the render window. This also
        sets the resolution of the saved WebP. It must be a sequence of two positive
        ints. The text and line widths scale with it, so a larger or smaller window is
        legible rather than being drawn with the same pixel counts as the default. When
        rendering on screen, a window_size that does not fit in the area the window
        manager grants is shown scaled down to fit, keeping its aspect ratio, with its
        text and line widths scaled down to match, while a saved WebP keeps the
        requested size. The default is (1024, 768).
    :param save: Set this to True to save the image as a WebP. The image is rendered
        separately, off screen, at window_size and from the view the window ends with,
        so resizing the window does not change the saved image's size or styling.
        However, a window resized to a different aspect ratio shows more or less of the
        scene horizontally than the saved image does. It can be a bool or a numpy bool
        and will be converted internally to a bool. The default is False.
    :param path: The file path to save the image to. It can be a str or a Path, must end
        with ".webp", and its directory must already exist. This has no effect unless
        save is True. The default is "draw.webp".
    :param quality: The quality of the saved WebP, where 0.0 is the smallest file with
        the most compression artifacts and 100.0 is the largest file with the fewest. It
        can be an int or a float and will be converted internally to a float. This has
        no effect unless save is True. The default is 75.0.
    :param testing: Set this to True to close the image after one second, which is
        useful for running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :return: None
    """
    if not isinstance(
        solver,
        (
            _steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            _steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
            _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
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
        if scalar_type not in _output_rendering.VALID_SCALAR_TYPES:
            valid_types = ", ".join(
                f'"{t}"' for t in _output_rendering.VALID_SCALAR_TYPES
            )
            raise ValueError(
                f'scalar_type must be None, {valid_types}, got "{scalar_type}".'
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
        _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        raise ValueError(
            "show_wake_vortices can only be True when drawing an "
            "UnsteadyRingVortexLatticeMethodSolver."
        )
    if show_wake_vortices and not solver.ran:
        raise RuntimeError(
            "solver must have run before drawing with show_wake_vortices set to True."
        )

    show_mujoco_geometry = _parameter_validation.boolLike_return_bool(
        show_mujoco_geometry, "show_mujoco_geometry"
    )
    if show_mujoco_geometry and not isinstance(
        solver,
        _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    ):
        raise ValueError(
            "show_mujoco_geometry can only be True when drawing a "
            "FreeFlightUnsteadyRingVortexLatticeMethodSolver."
        )

    if not isinstance(window_size, Sequence) or len(window_size) != 2:
        raise ValueError("window_size must be a sequence of two ints.")
    window_width = _parameter_validation.int_in_range_return_int(
        window_size[0], "window_size[0]", 0, False
    )
    window_height = _parameter_validation.int_in_range_return_int(
        window_size[1], "window_size[1]", 0, False
    )

    save = _parameter_validation.boolLike_return_bool(save, "save")
    path = _parameter_validation.pathLike_return_path(path, "path", (".webp",))
    quality = _parameter_validation.number_in_range_return_float(
        quality, "quality", 0.0, True, 100.0, True
    )
    testing = _parameter_validation.boolLike_return_bool(testing, "testing")

    # Create the Plotter whose window shows the drawing, along with the scale of its
    # text and line widths. Those differ from the saved WebP's when the requested window
    # does not fit on screen and is scaled down.
    plotter, preview_window_scale = create_preview_plotter(window_width, window_height)
    window_scale = _output_rendering.get_window_scale(window_width, window_height)

    # For a free flight solver, geometry is rendered in its true Earth-frame pose so the
    # body flies through the scene. T_pas_GP1_CgP1_to_E_Eo holds the passive
    # transformation from the first Airplane's geometry axes (relative to its CG) to
    # Earth axes (relative to the Earth origin) for the drawn time step, and stays None
    # for the standard body-fixed rendering in geometry axes.
    is_free_flight = isinstance(
        solver,
        _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    )
    T_pas_GP1_CgP1_to_E_Eo: np.ndarray | None = None

    # Get the solver's geometry and OperatingPoint, along with the wake ring vortex
    # surfaces when they are being shown. The wake stays None otherwise, which omits it
    # from the scene.
    wake_ring_vortex_surfaces: pv.PolyData | None = None
    if isinstance(
        solver,
        _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        draw_step = solver.num_steps - 1

        airplanes = solver.steady_problems[draw_step].airplanes
        draw_operating_point = solver.steady_problems[draw_step].operating_point
        qInf__E = draw_operating_point.qInf__E

        if is_free_flight:
            T_pas_GP1_CgP1_to_E_Eo = _output_rendering.get_free_flight_transformation(
                draw_operating_point
            )

        if show_wake_vortices:
            wake_ring_vortex_surfaces = _output_rendering.get_wake_ring_vortex_surfaces(
                solver, draw_step
            )
            if T_pas_GP1_CgP1_to_E_Eo is not None:
                wake_ring_vortex_surfaces = _output_rendering.transform_mesh(
                    wake_ring_vortex_surfaces, T_pas_GP1_CgP1_to_E_Eo
                )
    else:
        airplanes = solver.airplanes
        draw_operating_point = solver.operating_point
        qInf__E = draw_operating_point.qInf__E

    # Get the Panel surfaces, mapping them into Earth axes for free flight.
    panel_surfaces = _output_rendering.get_panel_surfaces(airplanes)
    if T_pas_GP1_CgP1_to_E_Eo is not None:
        panel_surfaces = _output_rendering.transform_mesh(
            panel_surfaces, T_pas_GP1_CgP1_to_E_Eo
        )

    reflect_T_act = draw_operating_point.surfaceReflect_T_act_GP1_CgP1
    image_surface_mesh: pv.PolyData | None = None
    image_surface_texture: pv.Texture | None = None

    # For free flight, the active reflection is represented in geometry axes, but the
    # geometry has been mapped into Earth axes. Re-expressing the reflection in Earth
    # axes (a change of basis by the same passive transformation) lets the
    # reflected-geometry code below operate entirely in Earth axes.
    if T_pas_GP1_CgP1_to_E_Eo is not None and reflect_T_act is not None:
        reflect_T_act = (
            T_pas_GP1_CgP1_to_E_Eo
            @ reflect_T_act
            @ _transformations.invert_T_pas(T_pas_GP1_CgP1_to_E_Eo)
        )

    # Choose the scalar coloring for the Panels, leaving it None to color them
    # uniformly.
    coloring: _output_rendering.ScalarColoring | None = None
    if scalar_type in ("induced drag", "crosswind force", "lift"):
        these_scalars = _output_rendering.get_scalars(airplanes, scalar_type, qInf__E)
        color_map, c_min, c_max = _output_rendering.choose_color_map(these_scalars)
        coloring = _output_rendering.ScalarColoring(
            scalars=these_scalars,
            scalar_type=scalar_type,
            min_scalar=float(min(these_scalars)),
            max_scalar=float(max(these_scalars)),
            color_map=color_map,
            muted_color_map=_output_rendering.mute_colormap(
                color_map, _output_rendering.IMAGE_REFLECTION_MUTE_FACTOR
            ),
            c_min=c_min,
            c_max=c_max,
        )

    # If showing MuJoCo geometry, gather the geoms that extra_xml injects, along with
    # the drawn time step's body axes to Earth axes transformation. The worldbody geoms
    # are left out of the camera's framing fit below, while the body geoms join it.
    worldbody_geoms: list[_mujoco_model.RenderGeom] = []
    body_geoms: list[_mujoco_model.RenderGeom] = []
    T_pas_BP1_CgP1_to_E_Eo: np.ndarray | None = None
    if show_mujoco_geometry:
        assert isinstance(
            solver,
            _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
        )
        worldbody_geoms, body_geoms = _output_rendering.get_mujoco_render_geometry(
            solver
        )
        if worldbody_geoms or body_geoms:
            T_pas_BP1_CgP1_to_E_Eo = (
                _output_rendering.get_free_flight_body_transformation(
                    draw_operating_point
                )
            )

    # If showing streamlines, get their surfaces, mapping them into Earth axes for free
    # flight. They stay None otherwise, which omits them from the scene.
    streamline_surfaces: pv.PolyData | None = None
    if show_streamlines:
        streamline_surfaces = _output_rendering.get_streamline_surfaces(
            solver.gridStreamlinePoints_GP1_CgP1
        )
        if T_pas_GP1_CgP1_to_E_Eo is not None:
            streamline_surfaces = _output_rendering.transform_mesh(
                streamline_surfaces, T_pas_GP1_CgP1_to_E_Eo
            )

    # Find the bounds of the geometry. They are found from the meshes rather than from a
    # Plotter's actors, so the image surface plane they size is known before any actor
    # is added. The scene meshes start as the geometry meshes and gain the image surface
    # plane when one is defined, which gives the free flight camera its framing bounds.
    geometry_meshes = get_draw_geometry_meshes(
        panel_surfaces, wake_ring_vortex_surfaces, streamline_surfaces, reflect_T_act
    )
    geometry_bounds = pv.MultiBlock(geometry_meshes).bounds
    scene_meshes = list(geometry_meshes)

    # If an image surface is defined, build its plane.
    if reflect_T_act is not None:
        if T_pas_GP1_CgP1_to_E_Eo is not None:
            # The image surface helper builds the plane from geometry-axis quantities,
            # so it needs geometry-axis bounds. Build the plane there, then map it into
            # Earth axes to match the rendered geometry.
            geometry_axis_bounds = _output_rendering.get_panel_surfaces(
                airplanes
            ).bounds
            image_surface_result = _output_rendering.get_image_surface_mesh_and_texture(
                draw_operating_point, geometry_axis_bounds
            )
            assert image_surface_result is not None
            image_surface_mesh, image_surface_texture = image_surface_result
            image_surface_mesh = _output_rendering.transform_mesh(
                image_surface_mesh, T_pas_GP1_CgP1_to_E_Eo
            )
        else:
            image_surface_result = _output_rendering.get_image_surface_mesh_and_texture(
                draw_operating_point, geometry_bounds
            )
            assert image_surface_result is not None
            image_surface_mesh, image_surface_texture = image_surface_result
        scene_meshes.append(image_surface_mesh)

    # Add the scene's actors.
    add_scene(
        plotter,
        panel_surfaces,
        wake_ring_vortex_surfaces,
        coloring,
        reflect_T_act,
        preview_window_scale,
        streamline_surfaces,
        image_surface_mesh,
        image_surface_texture,
        worldbody_geoms,
        body_geoms,
        T_pas_BP1_CgP1_to_E_Eo,
    )

    # For the standard body-fixed rendering with an image surface, fit the camera to the
    # geometry bounds so the view is not dominated by the much larger image surface
    # plane. When an image surface is present, cpos is not passed to show() because that
    # would trigger an auto-fit to all actors (including the image surface). Free flight
    # uses its own Earth-axes camera, computed below.
    if reflect_T_act is not None and T_pas_GP1_CgP1_to_E_Eo is None:
        plotter.camera.position = (-1, -1, 1)
        plotter.camera.focal_point = (0, 0, 0)
        plotter.camera.up = (0, 0, 1)
        plotter.reset_camera(bounds=geometry_bounds)  # type: ignore[call-arg]

    # Choose the camera position. Free flight frames the body in Earth axes with
    # physical up as Earth -z. The standard rendering views geometry axes from (-1, -1,
    # 1), unless an image surface is present, in which case the camera was already
    # fitted above and cpos is left None so show() does not auto-fit to the large image
    # surface plane.
    draw_cpos: tuple | None
    if is_free_flight:
        # Place the camera along the Earth-axes camera offset direction from the
        # geometry's center, with physical up, then fit the parallel projection to the
        # rendered geometry. The camera is set explicitly here (rather than via cpos) so
        # that show() preserves the fit, since passing a full position to cpos would not
        # set the parallel projection's scale.
        airplane_bounds = np.array(panel_surfaces.bounds, dtype=float)
        center_E_Eo = 0.5 * (airplane_bounds[1::2] + airplane_bounds[::2])
        airplane_diagonal = float(
            np.linalg.norm(airplane_bounds[1::2] - airplane_bounds[::2])
        )
        plotter.camera.focal_point = tuple(center_E_Eo)
        plotter.camera.position = tuple(
            center_E_Eo + 3.0 * airplane_diagonal * freeFlightCameraOffsetDirection_E
        )
        plotter.camera.up = freeFlightViewUp_E

        if T_pas_BP1_CgP1_to_E_Eo is not None:
            # Fit the camera to explicit framing bounds: the geometry and the image
            # surface plane, plus the body geoms posed at the drawn time step, plus
            # their reflections when an image surface is defined. The worldbody geoms
            # stay out of the fit, so a worldbody geom that is much larger than the
            # body, like a ground plane, cannot dominate it. The fit re-centers the
            # focal point on those bounds and keeps the view direction and up set above.
            posed_body_geom_meshes = [
                _output_rendering.transform_mesh(
                    render_geom.mesh, T_pas_BP1_CgP1_to_E_Eo
                )
                for render_geom in body_geoms
            ]
            if reflect_T_act is not None:
                posed_body_geom_meshes += [
                    _output_rendering.transform_mesh(
                        posed_body_geom_mesh, reflect_T_act
                    )
                    for posed_body_geom_mesh in posed_body_geom_meshes
                ]
            all_bounds = np.array(
                [pv.MultiBlock(scene_meshes).bounds]
                + [
                    posed_body_geom_mesh.bounds
                    for posed_body_geom_mesh in posed_body_geom_meshes
                ],
                dtype=float,
            )
            framing_bounds = np.empty(6, dtype=float)
            framing_bounds[0::2] = all_bounds[:, 0::2].min(axis=0)
            framing_bounds[1::2] = all_bounds[:, 1::2].max(axis=0)
            plotter.reset_camera(bounds=tuple(framing_bounds))  # type: ignore[call-arg]

            # Re-fit only the clipping range to all actors, which keeps the fitted
            # framing while ensuring the worldbody geoms are not cut off.
            plotter.reset_camera_clipping_range()
        else:
            plotter.reset_camera()  # type: ignore[call-arg]
        draw_cpos = None
    elif image_surface_mesh is None:
        draw_cpos = (-1, -1, 1)
    else:
        draw_cpos = None

    # The scene is translucent when it holds the image surface plane or a geom whose
    # rgba has an alpha below one. Settle the scalar bar layout before such a drawing is
    # displayed. This is the first of the two passes it takes, and show below is the
    # second, so the labels are in place by the time the user sees anything.
    scene_is_translucent = reflect_T_act is not None or any(
        float(render_geom.rgba[3]) < 1.0 for render_geom in worldbody_geoms + body_geoms
    )
    if scene_is_translucent:
        _output_rendering.settle_scalar_bar_layout(plotter)

    if not testing:
        # Show the Plotter so the user can adjust the camera position and window. When
        # the user closes the window, the Plotter still exists. Therefore, it can later
        # be saved as an image if desired.
        plotter.show(
            title="Orient the view, then press any key to continue.",
            cpos=draw_cpos,
            full_screen=False,
            auto_close=False,
        )
    else:
        # Show the Plotter for 1 second, then proceed automatically. This is useful for
        # testing.
        plotter.show(
            title="Testing. Please wait.",
            cpos=draw_cpos,
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # If saving, render the scene again in a separate off-screen Plotter, from the
    # camera the window ended with, and save that as a WebP. The user can resize the
    # window in ways that cannot be prevented on every platform, and the window's text
    # and line widths stay scaled to the size it opened at, so capturing the window
    # could save an image of a different size with mismatched styling. The off-screen
    # Plotter is always the requested size. Its scalar bar layout is settled whether or
    # not its scene is translucent, since an opaque scene's layout is already settled
    # and the extra pass leaves it unchanged.
    if save:
        save_plotter = create_plotter(window_width, window_height, True)
        add_scene(
            save_plotter,
            panel_surfaces,
            wake_ring_vortex_surfaces,
            coloring,
            reflect_T_act,
            window_scale,
            streamline_surfaces,
            image_surface_mesh,
            image_surface_texture,
            worldbody_geoms,
            body_geoms,
            T_pas_BP1_CgP1_to_E_Eo,
        )
        save_plotter.camera.DeepCopy(plotter.camera)
        save_plotter.renderer.camera_set = True
        _output_rendering.settle_scalar_bar_layout(save_plotter)
        image = _output_rendering.screenshot_image(save_plotter)
        save_plotter.close()

        # webp annotates file_path as a str, so the Path is converted at the boundary.
        webp.save_image(
            img=image,
            file_path=str(path),
            lossless=False,
            quality=quality,
            method=_output_rendering.WEBP_METHOD,
        )

    # Close the Plotter.
    plotter.close()


def animate(
    unsteady_solver: _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    scalar_type: str | None = None,
    show_wake_vortices: bool | np.bool = False,
    show_mujoco_geometry: bool | np.bool = False,
    window_size: Sequence[int] = (1024, 768),
    save: bool | np.bool = False,
    path: str | Path = "animate.webp",
    quality: int | float = 75.0,
    speed: int | float | None = None,
    testing: bool | np.bool = False,
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
        "crosswind force", or "lift", which respectively use each Panel's induced drag,
        crosswind force, and lift coefficient. The default is None.
    :param show_wake_vortices: Set this to True to show any wake ring vortices. If True,
        the solver must have already been run. Can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :param show_mujoco_geometry: Set this to True to show the MuJoCo geometry that
        extra_xml and mujoco_assets inject into the solver's model, with body geoms re-
        posed every time step and worldbody geoms static in Earth axes. If True, the
        unsteady_solver must be a FreeFlightUnsteadyRingVortexLatticeMethodSolver. Can
        be a bool or a numpy bool and will be converted internally to a bool. The
        default is False.
    :param window_size: The width and height, in pixels, of the render window. This also
        sets the resolution of the saved WebP. It must be a sequence of two positive
        ints. The text and line widths scale with it, so a larger or smaller window is
        legible rather than being drawn with the same pixel counts as the default. When
        rendering on screen, a window_size that does not fit in the area the window
        manager grants is shown scaled down to fit, keeping its aspect ratio, with its
        text and line widths scaled down to match, while a saved WebP keeps the
        requested size. The default is (1024, 768).
    :param save: Set this to True to save the animation as an animated WebP. Once the
        view is oriented, the window closes and the frames are rendered separately, off
        screen, at window_size and from the view the window ended with, while a progress
        bar reports the rendering's progress. Resizing the window therefore does not
        change the saved frames' size or styling. However, a window resized to a
        different aspect ratio shows more or less of the scene horizontally than the
        saved frames do. If False, the animation plays in the window instead. It can be
        a bool or a numpy bool and will be converted internally to a bool. The default
        is False.
    :param path: The file path to save the animation to. It can be a str or a Path, must
        end with ".webp", and its directory must already exist. This has no effect
        unless save is True. The default is "animate.webp".
    :param quality: The quality of the saved WebP, where 0.0 is the smallest file with
        the most compression artifacts and 100.0 is the largest file with the fewest. It
        can be an int or a float and will be converted internally to a float. This has
        no effect unless save is True. The default is 75.0.
    :param speed: The playback speed of the saved animation, as a multiple of real time,
        where 1.0 plays the simulation at true speed and 0.5 plays it at half speed.
        Setting this to None plays at true speed, slowing the animation down only when
        true speed would need more than 50.0 frames per second, which is the fastest
        some programs will render a WebP. Any speed that would need more than 50.0
        frames per second is instead reached by saving only every Nth frame, which
        trades the animation's temporal resolution for its speed. It can be an int or a
        float, must be positive, and will be converted internally to a float. A speed so
        slow that the animation would save fewer than one frame per second of playback,
        or so fast that it would save fewer than two frames in total, is rejected. This
        has no effect unless save is True. The default is None.
    :param testing: Set this to True to start the animation after one second, which is
        useful for running test suites. It can be a bool or a numpy bool and will be
        converted internally to a bool. The default is False.
    :return: None
    """
    if not isinstance(
        unsteady_solver,
        _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
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
        if scalar_type not in _output_rendering.VALID_SCALAR_TYPES:
            valid_types = ", ".join(
                f'"{t}"' for t in _output_rendering.VALID_SCALAR_TYPES
            )
            raise ValueError(
                f'scalar_type must be None, {valid_types}, got "{scalar_type}".'
            )

    show_wake_vortices = _parameter_validation.boolLike_return_bool(
        show_wake_vortices, "show_wake_vortices"
    )
    if show_wake_vortices and not unsteady_solver.ran:
        raise RuntimeError(
            "unsteady_solver must have run before animating with show_wake_vortices set"
            " to True."
        )

    show_mujoco_geometry = _parameter_validation.boolLike_return_bool(
        show_mujoco_geometry, "show_mujoco_geometry"
    )
    if show_mujoco_geometry and not isinstance(
        unsteady_solver,
        _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    ):
        raise ValueError(
            "show_mujoco_geometry can only be True when animating a "
            "FreeFlightUnsteadyRingVortexLatticeMethodSolver."
        )

    if not isinstance(window_size, Sequence) or len(window_size) != 2:
        raise ValueError("window_size must be a sequence of two ints.")
    window_width = _parameter_validation.int_in_range_return_int(
        window_size[0], "window_size[0]", 0, False
    )
    window_height = _parameter_validation.int_in_range_return_int(
        window_size[1], "window_size[1]", 0, False
    )

    save = _parameter_validation.boolLike_return_bool(save, "save")
    path = _parameter_validation.pathLike_return_path(path, "path", (".webp",))
    quality = _parameter_validation.number_in_range_return_float(
        quality, "quality", 0.0, True, 100.0, True
    )
    if speed is not None:
        speed = _parameter_validation.number_in_range_return_float(
            speed, "speed", 0.0, False
        )
    testing = _parameter_validation.boolLike_return_bool(testing, "testing")

    first_results_step = unsteady_solver.first_results_step

    # Get the solver's SteadyProblems' Airplanes. This will become a list of lists, with
    # the first index being the time step and the second index identifying each Airplane
    # at that time step.
    step_airplanes = []
    for steady_problem in unsteady_solver.steady_problems:
        step_airplanes.append(steady_problem.airplanes)

    # For a free flight solver, each time step's geometry is rendered in its true
    # Earth-frame pose so the body flies through the scene. step_transforms holds, per
    # time step, the passive transformation from the first Airplane's geometry axes
    # (relative to its CG) to Earth axes (relative to the Earth origin). It stays empty
    # for the standard body-fixed rendering in geometry axes.
    is_free_flight = isinstance(
        unsteady_solver,
        _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    )
    step_transforms: list[np.ndarray] = []
    if is_free_flight:
        step_transforms = [
            _output_rendering.get_free_flight_transformation(
                steady_problem.operating_point
            )
            for steady_problem in unsteady_solver.steady_problems
        ]

    # If showing MuJoCo geometry, gather the geoms that extra_xml injects, along with
    # each time step's body axes to Earth axes transformation. The worldbody geoms stay
    # fixed in Earth axes while the body geoms are re-posed every frame.
    worldbody_geoms: list[_mujoco_model.RenderGeom] = []
    body_geoms: list[_mujoco_model.RenderGeom] = []
    step_body_transforms: list[np.ndarray] = []
    if show_mujoco_geometry:
        assert isinstance(
            unsteady_solver,
            _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
        )
        worldbody_geoms, body_geoms = _output_rendering.get_mujoco_render_geometry(
            unsteady_solver
        )
        step_body_transforms = [
            _output_rendering.get_free_flight_body_transformation(
                steady_problem.operating_point
            )
            for steady_problem in unsteady_solver.steady_problems
        ]

    # Resolve the playback speed into the frame stride, frame rate, and text overlays
    # that describe how the saved animation steps through the time steps.
    playback = _output_rendering.resolve_playback(unsteady_solver, speed, save)

    # Create the Plotter whose window shows the preview, along with the scale of its
    # text and line widths. Those differ from the saved WebP's when the requested window
    # does not fit on screen and is scaled down.
    plotter, preview_window_scale = create_preview_plotter(window_width, window_height)
    window_scale = _output_rendering.get_window_scale(window_width, window_height)

    # Initialize values to hold the color map choice and its limits.
    c_min = 0.0
    c_max = 0.0
    color_map: matplotlib.colors.Colormap | None = None
    muted_color_map: matplotlib.colors.Colormap | None = None

    # Initialize variables to hold the SteadyProblems' scalars and their attributes.
    all_scalars = np.empty(0, dtype=float)
    min_scalar = 0.0
    max_scalar = 0.0

    # If coloring the Panels based on scalars, gather all the scalars across all the
    # time steps and Airplanes. These will be used to set the color map limits.
    if scalar_type is not None:
        for step_id, airplanes in enumerate(step_airplanes):
            scalars_to_add = _output_rendering.get_scalars(
                airplanes,
                scalar_type,
                unsteady_solver.steady_problems[step_id].operating_point.qInf__E,
            )
            all_scalars = np.hstack((all_scalars, scalars_to_add))

        # Choose the color map from the scalars across all the time steps, so that one
        # map and one pair of limits apply to every frame. Mute it once here as well,
        # since building the muted copy walks a 256 entry table and every frame that
        # carries reflected geometry needs it.
        color_map, c_min, c_max = _output_rendering.choose_color_map(all_scalars)
        muted_color_map = _output_rendering.mute_colormap(
            color_map, _output_rendering.IMAGE_REFLECTION_MUTE_FACTOR
        )

        min_scalar = float(min(all_scalars))
        max_scalar = float(max(all_scalars))

    # Pre-compute the image surface plane and the reflection that maps geometry onto its
    # far side. Both are static across the animation, so they are built once here from
    # the last time step's geometry and reused for every frame.
    last_step = len(step_airplanes) - 1
    last_step_operating_point = unsteady_solver.steady_problems[
        last_step
    ].operating_point
    (
        image_surface_mesh,
        image_surface_texture,
        reflect_T_act,
        image_surface_geometry_bounds,
    ) = _output_rendering.get_animation_image_surface(
        unsteady_solver,
        step_airplanes,
        step_transforms,
        is_free_flight,
        show_wake_vortices,
    )
    animate_text_color = (
        _output_rendering.TEXT_COLOR_SURFACE
        if reflect_T_act is not None
        else _output_rendering.TEXT_COLOR
    )

    # For free flight, compute a fixed camera that frames the whole trajectory. The body
    # moves through the scene, so the camera is centered on the trajectory's midpoint
    # and the parallel scale is fit to the projected extent of the trajectory's geometry
    # so the whole glide stays in view. The clipping range is sized later, after the
    # user has oriented the view, so it tracks the user's chosen camera rather than the
    # default.
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

        # Map the first and last frames' Panel surfaces into Earth axes. These two
        # frames bound the trajectory, so their combined extent frames the whole glide.
        first_step_panel_surfaces = _output_rendering.transform_mesh(
            _output_rendering.get_panel_surfaces(step_airplanes[0]), step_transforms[0]
        )
        last_step_panel_surfaces = _output_rendering.transform_mesh(
            _output_rendering.get_panel_surfaces(step_airplanes[last_step]),
            step_transforms[last_step],
        )
        airplane_bounds = np.array(first_step_panel_surfaces.bounds, dtype=float)
        airplane_diagonal = float(
            np.linalg.norm(airplane_bounds[1::2] - airplane_bounds[::2])
        )

        # Place the camera along the Earth-axes camera offset direction from the
        # trajectory's midpoint, with physical up, far enough back to clear the geometry
        # at both ends.
        padding = max(2.0 * airplane_diagonal, 0.5 * trajectory_extent)
        camera_distance = trajectory_extent + padding
        cameraPosition_E_Eo = (
            trajectoryMidpoint_E_Eo
            + camera_distance * freeFlightCameraOffsetDirection_E
        )
        free_flight_cpos = [
            tuple(cameraPosition_E_Eo),
            tuple(trajectoryMidpoint_E_Eo),
            freeFlightViewUp_E,
        ]

        # Collect the geometry that frames the trajectory: the body at both ends, plus
        # the last frame's wake (the largest) if it is shown. The last frame's surfaces
        # are reused below to size the clipping range after the user orients the view.
        framing_meshes = [first_step_panel_surfaces, last_step_panel_surfaces]
        free_flight_clip_meshes = [last_step_panel_surfaces]
        if show_wake_vortices:
            last_step_wake_surfaces = _output_rendering.transform_mesh(
                _output_rendering.get_wake_ring_vortex_surfaces(
                    unsteady_solver, last_step
                ),
                step_transforms[last_step],
            )
            if last_step_wake_surfaces.n_points > 0:
                framing_meshes.append(last_step_wake_surfaces)
                free_flight_clip_meshes.append(last_step_wake_surfaces)

        # The body geoms fly with the body, so they join the framing fit posed at both
        # ends of the trajectory, keeping a geom that extends past the Panel surfaces
        # from being cropped. Only the last time step's posed copies join the clipping
        # meshes, since the first time step's geom actors are already present when the
        # clipping range is sized. The worldbody geoms join neither: their actors are
        # also present when the clipping range is sized, which keeps them visible
        # without letting a worldbody geom that is much larger than the body, like a
        # ground plane, dominate the framing fit.
        first_step_body_geom_meshes = [
            _output_rendering.transform_mesh(render_geom.mesh, step_body_transforms[0])
            for render_geom in body_geoms
        ]
        last_step_body_geom_meshes = [
            _output_rendering.transform_mesh(
                render_geom.mesh, step_body_transforms[last_step]
            )
            for render_geom in body_geoms
        ]
        framing_meshes += first_step_body_geom_meshes + last_step_body_geom_meshes
        free_flight_clip_meshes += last_step_body_geom_meshes

        # Fit the parallel scale (half the viewport height in world units, since the
        # projection is parallel) to the projected extent of that geometry about the
        # focal point. This frames the glide snugly. The user can rescale interactively
        # before the animation is captured.
        free_flight_parallel_scale = (
            _output_rendering.get_free_flight_fit_parallel_scale(
                framing_meshes,
                trajectoryMidpoint_E_Eo,
                freeFlightCameraOffsetDirection_E,
                freeFlightViewUp_E,
            )
        )

    # If saving the animation, add the text overlays that describe its playback.
    if save:
        _output_rendering.add_playback_overlays(
            plotter, playback, preview_window_scale, animate_text_color
        )

    # Show the first and last time steps together during framing. This gives the user
    # both the initial body pose and the final wake/trajectory extent without
    # multiplying actors across every time step of a long simulation.
    first_panel_surfaces = _output_rendering.get_panel_surfaces(step_airplanes[0])
    if is_free_flight:
        first_panel_surfaces = _output_rendering.transform_mesh(
            first_panel_surfaces, step_transforms[0]
        )
    first_frame_coloring: _output_rendering.ScalarColoring | None = None
    if scalar_type is not None and first_results_step == 0:
        assert color_map is not None
        assert muted_color_map is not None
        first_frame_coloring = _output_rendering.ScalarColoring(
            _output_rendering.get_scalars(
                step_airplanes[0],
                scalar_type,
                unsteady_solver.steady_problems[0].operating_point.qInf__E,
            ),
            scalar_type,
            min_scalar,
            max_scalar,
            color_map,
            muted_color_map,
            c_min,
            c_max,
        )

    first_preview_actors = _output_rendering.add_frame_geometry(
        plotter,
        first_panel_surfaces,
        None,
        first_frame_coloring,
        reflect_T_act,
        preview_window_scale,
    )
    if last_step != 0:
        set_preview_opacity(first_preview_actors, ANIMATE_PREVIEW_FIRST_OPACITY)

    if show_mujoco_geometry:
        first_mujoco_actors = _output_rendering.add_mujoco_geometry(
            plotter, worldbody_geoms, body_geoms, step_body_transforms[0], reflect_T_act
        )
        if last_step != 0:
            set_preview_opacity(first_mujoco_actors, ANIMATE_PREVIEW_FIRST_OPACITY)

    if last_step != 0:
        last_panel_surfaces = _output_rendering.get_panel_surfaces(
            step_airplanes[last_step]
        )
        if is_free_flight:
            last_panel_surfaces = _output_rendering.transform_mesh(
                last_panel_surfaces, step_transforms[last_step]
            )

        last_wake_surfaces = None
        if show_wake_vortices:
            last_wake_surfaces = _output_rendering.get_wake_ring_vortex_surfaces(
                unsteady_solver, last_step
            )
            if is_free_flight:
                last_wake_surfaces = _output_rendering.transform_mesh(
                    last_wake_surfaces, step_transforms[last_step]
                )

        last_coloring: _output_rendering.ScalarColoring | None = None
        if scalar_type is not None:
            assert color_map is not None
            assert muted_color_map is not None
            last_coloring = _output_rendering.ScalarColoring(
                _output_rendering.get_scalars(
                    step_airplanes[last_step],
                    scalar_type,
                    unsteady_solver.steady_problems[last_step].operating_point.qInf__E,
                ),
                scalar_type,
                min_scalar,
                max_scalar,
                color_map,
                muted_color_map,
                c_min,
                c_max,
            )
        last_preview_actors = _output_rendering.add_frame_geometry(
            plotter,
            last_panel_surfaces,
            last_wake_surfaces,
            last_coloring,
            reflect_T_act,
            preview_window_scale,
        )
        set_preview_opacity(last_preview_actors, ANIMATE_PREVIEW_LAST_OPACITY)

        if show_mujoco_geometry:
            last_mujoco_actors = _output_rendering.add_mujoco_geometry(
                plotter,
                [],
                body_geoms,
                step_body_transforms[last_step],
                reflect_T_act,
            )
            set_preview_opacity(last_mujoco_actors, ANIMATE_PREVIEW_LAST_OPACITY)

    # If an image surface is defined, plot the pre-computed plane, set the camera
    # direction, and fit the camera to the last time step's geometry bounds so the view
    # is not dominated by the much larger image surface plane. When an image surface is
    # present, cpos is not passed to show() because that would trigger an auto-fit to
    # all actors (including the image surface).
    if reflect_T_act is not None:
        assert image_surface_mesh is not None

        # Add the image surface plane.
        plotter.add_mesh(
            image_surface_mesh,
            texture=image_surface_texture,
            opacity=IMAGE_SURFACE_OPACITY,
            smooth_shading=True,
            lighting=False,
            render=False,
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
    # rescaling carry through to the animation. The standard rendering views geometry
    # axes from (-1, -1, 1), unless an image surface is present, in which case the
    # camera was already fitted above and cpos is left None so show() does not auto-fit
    # to the large image surface plane.
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

    # The frames are translucent when the scene holds the image surface plane or a geom
    # whose rgba has an alpha below one, and the held first frame is translucent
    # whenever it has ghosts as well. Give a translucent held frame the scalar bar
    # layout pass that translucent frames in the loop get, so the view held for the user
    # matches the animation that follows it. This is the first of the two passes the
    # layout takes to settle, and show below is the second, so the labels are in place
    # by the time the user sees anything.
    scene_is_translucent = reflect_T_act is not None or any(
        float(render_geom.rgba[3]) < 1.0 for render_geom in worldbody_geoms + body_geoms
    )
    if scene_is_translucent or last_step != 0:
        _output_rendering.settle_scalar_bar_layout(plotter)

    # If not testing, show the Plotter with the first time step so the user can orient
    # the view. When the user presses any key, set the title back to the animation title
    # and proceed. If testing, show the Plotter with the first time step for 1 second,
    # and start the animation with the current window view.
    if not testing:
        title = "Orient the view, then press any key."
        if last_step != 0:
            title += " The ghosts preview the animation and are not saved."
        plotter.show(
            title=title,
            cpos=animate_cpos,
            full_screen=False,
            auto_close=False,
        )

    else:
        plotter.show(
            title="Testing. Please wait.",
            cpos=animate_cpos,
            full_screen=False,
            interactive=False,
            auto_close=False,
        )
        time.sleep(1)

    # If saving, render the frames in a separate off-screen Plotter, from the camera the
    # window ended with, and close the window. The user can resize the window in ways
    # that cannot be prevented on every platform, and the window's text and line widths
    # stay scaled to the size it opened at, so frames captured from it could change size
    # partway through an animation and have mismatched styling. The off-screen Plotter
    # is always the requested size, and it starts with the playback text overlays. If
    # not saving, the frames play in the window, which is cleared of the preview first,
    # and keep its scale. Saving or not, the first frame then starts from a Plotter
    # without any of the preview, so it is built in full, with the image surface plane.
    if save:
        frame_plotter = create_plotter(window_width, window_height, True)
        frame_window_scale = window_scale
        frame_plotter.camera.DeepCopy(plotter.camera)
        frame_plotter.renderer.camera_set = True
        _output_rendering.add_playback_overlays(
            frame_plotter, playback, window_scale, animate_text_color
        )
        plotter.close()
    else:
        frame_plotter = plotter
        frame_window_scale = preview_window_scale
        assert plotter.ren_win is not None
        plotter.ren_win.SetWindowName(
            "Rendering animation. Please leave the window open until rendering "
            "finishes."
        )
        plotter.render()
        plotter.clear()

    # Build the first frame as the actual animation frame, and settle its scalar bar
    # layout. The layout is settled whether or not the frame is translucent, since an
    # opaque frame's layout is already settled and the extra pass leaves it unchanged.
    # The first time step has not shed a wake yet.
    first_frame_panel_surfaces = _output_rendering.get_panel_surfaces(step_airplanes[0])
    if is_free_flight:
        first_frame_panel_surfaces = _output_rendering.transform_mesh(
            first_frame_panel_surfaces, step_transforms[0]
        )
    add_scene(
        frame_plotter,
        first_frame_panel_surfaces,
        None,
        first_frame_coloring,
        reflect_T_act,
        frame_window_scale,
        None,
        image_surface_mesh,
        image_surface_texture,
        worldbody_geoms,
        body_geoms,
        step_body_transforms[0] if show_mujoco_geometry else None,
    )
    _output_rendering.settle_scalar_bar_layout(frame_plotter)

    # The user may have reoriented or rescaled the view during the held first frame.
    # Preserve that camera and only size the clipping range so every frame stays
    # visible: temporarily add the last frame's geometry (the first frame's is already
    # present), fit the clipping range to both, then remove the temporary actors. The
    # body moves through the scene, so a clipping range fit to the first frame alone
    # would clip later frames.
    if is_free_flight:
        temporary_actors = [
            frame_plotter.add_mesh(clip_mesh, lighting=False, render=False)
            for clip_mesh in free_flight_clip_meshes
        ]
        frame_plotter.reset_camera_clipping_range()
        free_flight_clipping_range = frame_plotter.camera.clipping_range
        for temporary_actor in temporary_actors:
            frame_plotter.remove_actor(temporary_actor, render=False)
        frame_plotter.camera.clipping_range = free_flight_clipping_range

    frame_plotter.render()

    # Show the rendering's progress, one time step at a time. When saving, the window
    # has closed, so this is the only sign that the animation is still being rendered.
    # When not saving, the window itself shows the progress, so the bar is disabled.
    progress_bar = tqdm(
        total=len(step_airplanes),
        unit="",
        ncols=100,
        desc="Rendering",
        disable=not save,
        bar_format="{desc}:{percentage:3.0f}% |{bar}| Elapsed: {elapsed}, "
        "Remaining: {remaining}",
    )

    # If saving, start the writer that encodes each frame into the WebP as it is
    # captured, so the frames never accumulate in memory, and hand it this first frame.
    animation_writer = None
    if save:
        # Ask the driver not to wait for the display's vertical blank before each buffer
        # swap. With the wait in place, every captured frame costs at least one refresh
        # period however little work it holds, so capture speed follows the display's
        # refresh rate. The frames are read back from the buffer, so the tearing the
        # wait prevents never reaches the file. This is a request the platform may
        # ignore, in which case capture simply keeps the refresh rate's pace. It is not
        # made when not saving, since the window is then what the user watches, and the
        # wait is the only pacing the playback has.
        assert frame_plotter.ren_win is not None
        frame_plotter.ren_win.SetSwapControl(0)

        animation_writer = _output_rendering.AnimationWriter(
            path, playback.frame_rate, quality
        )
        animation_writer.add_frame(_output_rendering.screenshot_image(frame_plotter))
    progress_bar.update()

    # Initialize a variable to keep track of the current time step.
    current_step = 1

    # Begin to iterate through the Airplane(s) from the subsequent time steps.
    for airplanes in step_airplanes[1:]:

        # Clear the Plotter.
        frame_plotter.clear()

        # Get the Panel surfaces of this time step's Airplane(s), mapping them into
        # Earth axes for free flight.
        panel_surfaces = _output_rendering.get_panel_surfaces(airplanes)
        if is_free_flight:
            panel_surfaces = _output_rendering.transform_mesh(
                panel_surfaces, step_transforms[current_step]
            )

        # If saving the animation, add the text overlays that describe its playback.
        if save:
            _output_rendering.add_playback_overlays(
                frame_plotter, playback, window_scale, animate_text_color
            )

        # If showing wake ring vortices, get their surfaces, mapping them into Earth
        # axes for free flight. They stay None otherwise, which omits them from the
        # scene.
        wake_ring_vortex_surfaces = None
        if show_wake_vortices:
            wake_ring_vortex_surfaces = _output_rendering.get_wake_ring_vortex_surfaces(
                unsteady_solver, current_step
            )
            if is_free_flight:
                wake_ring_vortex_surfaces = _output_rendering.transform_mesh(
                    wake_ring_vortex_surfaces, step_transforms[current_step]
                )

        # Choose this time step's scalar coloring, leaving it None to color the Panels
        # uniformly. Time steps before the first one with results have no scalars.
        coloring = None
        if scalar_type is not None and first_results_step <= current_step:
            assert color_map is not None
            assert muted_color_map is not None
            coloring = _output_rendering.ScalarColoring(
                _output_rendering.get_scalars(
                    airplanes,
                    scalar_type,
                    unsteady_solver.steady_problems[
                        current_step
                    ].operating_point.qInf__E,
                ),
                scalar_type,
                min_scalar,
                max_scalar,
                color_map,
                muted_color_map,
                c_min,
                c_max,
            )

        # Add this time step's scene, which includes the pre-computed image surface
        # plane when one is defined.
        add_scene(
            frame_plotter,
            panel_surfaces,
            wake_ring_vortex_surfaces,
            coloring,
            reflect_T_act,
            frame_window_scale,
            None,
            image_surface_mesh,
            image_surface_texture,
            worldbody_geoms,
            body_geoms,
            step_body_transforms[current_step] if show_mujoco_geometry else None,
        )

        # If the frame is translucent, settle the scalar bar layout before it is
        # displayed, leaving the second of its two passes to the render below.
        if scene_is_translucent:
            _output_rendering.settle_scalar_bar_layout(frame_plotter)

        # Render the assembled frame. Every add above is made with render=False, so
        # without this the frame would never reach the screen. Rendering once the frame
        # is whole is invisible, unlike the renders the adds used to trigger.
        frame_plotter.render()

        # If saving, hand a WebP Image of this frame to the writer. Only the time steps
        # that are multiples of the stride are saved, so a speed the maximum frame rate
        # cannot carry drops the ones in between.
        if animation_writer is not None and current_step % playback.keep_every == 0:
            animation_writer.add_frame(
                _output_rendering.screenshot_image(frame_plotter)
            )

        # Increment the time step tracker, and advance the progress bar.
        current_step += 1
        progress_bar.update()

    progress_bar.close()

    # If saving, finish the animation and write it to its file.
    if animation_writer is not None:
        animation_writer.close()

    # Close the Plotter the frames were rendered in. When saving, the window's Plotter
    # was already closed before the frames were rendered.
    frame_plotter.close()


def plot_results_versus_time(
    unsteady_solver: _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    show: bool | np.bool = True,
    figure_size_in: Sequence[int | float] = (6.4, 4.8),
    save: bool | np.bool = False,
    save_csv: bool | np.bool = False,
    directory: str | Path = ".",
    prefix: str = "",
    resolution_dpi: int | float = 300.0,
    file_format: str = "png",
    show_titles: bool | np.bool = True,
    font_size: int | float | None = None,
    text_color: str | Sequence[float] | None = None,
    line_width: int | float | None = None,
) -> None:
    """Plots the loads and load coefficients of an UnsteadyRingVortexLatticeMethodSolver
    or one of its subclasses (the aeroelastic or free flight solver) as a function of
    time.

    For a FreeFlightUnsteadyRingVortexLatticeMethodSolver, this also plots the first
    Airplane's six-degree-of-freedom state history: its position, velocity, orientation,
    angular velocity, and aerodynamic angles versus time. These describe the first
    Airplane, the rigid body the dynamics integrate, so they are plotted once for the
    whole simulation rather than per Airplane.

    Each file is named after the Airplane it describes, so an Airplane whose name holds
    a path separator is rejected when save or save_csv is True. Such a name would
    compose a destination outside the directory that was asked for.

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver whose loads and
        load coefficients will be plotted. Its subclasses, the
        AeroelasticUnsteadyRingVortexLatticeMethodSolver and the
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, are also accepted. For a
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, the first Airplane's state
        history is plotted as well.
    :param show: Set this to True to show the plots. It can be a bool or a numpy bool
        and will be converted internally to a bool. The default is True.
    :param figure_size_in: The width and height, in inches, of each figure. Multiplying
        this by resolution_dpi gives the resolution of each saved PNG. It must be a
        sequence of two positive numbers. The default is (6.4, 4.8), which is
        Matplotlib's own default.
    :param save: Set this to True to save the plots in the format named by file_format.
        It can be a bool or a numpy bool and will be converted internally to a bool. The
        default is False.
    :param save_csv: Set this to True to save the plotted data as CSVs, which is
        independent of save, so the data can be exported without rendering any images.
        One file holds the loads and load coefficients of each Airplane, and, for a
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, a second holds the first
        Airplane's state history. The two are separate files because the load histories
        begin at the solver's first results step while the state history begins at time
        step 0, so they have different numbers of rows. It can be a bool or a numpy bool
        and will be converted internally to a bool. The default is False.
    :param directory: The directory to save the plots and CSVs in. It can be a str or a
        Path and must already exist. This has no effect unless save or save_csv is True.
        The default is ".", the current working directory.
    :param prefix: A prefix to prepend to each file's name, which distinguishes one
        run's output from another's. It must be a str and must be a file name component
        rather than a path. With the default, an empty string, each file is named after
        its Airplane and the quantity it holds, as in "example_airplane_forces.png".
        With a non-empty prefix, each name becomes
        "<prefix>_<airplane>_<quantity>.<ext>". The Airplane's name is included either
        way, since dropping it would collide across the Airplanes of a formation
        simulation. This has no effect unless save or save_csv is True.
    :param resolution_dpi: The dots per inch at which to save each PNG. It can be an int
        or a float and will be converted internally to a float. This has no effect
        unless save is True and file_format is "png". The default is 300.0.
    :param file_format: The file format to save the plots in, which is also each plot
        file's extension. It must be "png", "svg", or "pdf". A PNG is a raster image
        with a transparent background. An SVG and a PDF are vector images, with their
        text kept as selectable text in an embedded font. Use a PDF for documents, such
        as LaTeX papers, and an SVG for the web and for slides. This has no effect
        unless save is True. The default is "png".
    :param show_titles: Set this to False to omit each figure's title and subtitle, for
        instance when a caption in a paper will carry them instead. It can be a bool or
        a numpy bool and will be converted internally to a bool. The default is True.
    :param font_size: The size, in points, of every piece of text in each figure. Text
        is sized against figure_size_in, so for a paper, set figure_size_in to the width
        the figure will print at and font_size to the size the publisher asks for. It
        can be an int or a float and will be converted internally to a float. It must be
        positive. Pass None to size each piece of text by Matplotlib's defaults, with
        the subtitle smaller than the rest. The default is None.
    :param text_color: The color of the text, the axis spines, and the ticks in each
        figure. It can be any color Matplotlib accepts, such as a name like "black", a
        hex string like "#333333", or a sequence of three or four numbers from 0.0 to
        1.0. Pass None to use the gray the rendered visualizations' text uses. The
        default is None.
    :param line_width: The middle line width, in points, in each figure. A figure's
        lines are drawn from thickest to thinnest so that all remain visible where they
        overlap, with widths spread evenly from 1.4 times this width down to 0.6 times
        it, and the legend draws every line at this width. Like font_size, it is
        measured against figure_size_in, so a figure printed at a small width usually
        calls for thinner lines. It can be an int or a float and will be converted
        internally to a float. It must be positive. Pass None to use a width of 2.5,
        which spreads the lines from 3.5 to 1.5. The default is None.
    :return: None
    """
    if not isinstance(
        unsteady_solver,
        _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        raise TypeError(
            "unsteady_solver must be an " "UnsteadyRingVortexLatticeMethodSolver."
        )
    show = _parameter_validation.boolLike_return_bool(show, "show")

    if not isinstance(figure_size_in, Sequence) or len(figure_size_in) != 2:
        raise ValueError("figure_size_in must be a sequence of two numbers.")
    figure_width_in = _parameter_validation.number_in_range_return_float(
        figure_size_in[0], "figure_size_in[0]", 0.0, False
    )
    figure_height_in = _parameter_validation.number_in_range_return_float(
        figure_size_in[1], "figure_size_in[1]", 0.0, False
    )

    save = _parameter_validation.boolLike_return_bool(save, "save")
    save_csv = _parameter_validation.boolLike_return_bool(save_csv, "save_csv")

    # A missing directory is an error rather than something to create, matching what
    # pathLike_return_path does for draw's and animate's single file destinations.
    if not isinstance(directory, (str, Path)):
        raise TypeError("directory must be a str or a Path.")
    directory = Path(directory)
    if directory.exists() and not directory.is_dir():
        raise ValueError(f'directory must be a directory, got file "{directory}".')
    if not directory.is_dir():
        raise ValueError(
            f'directory "{directory}" does not exist. Create it first, or choose a '
            f"destination that already exists."
        )

    # A prefix carrying a separator would compose a path into a subdirectory that this
    # function never checked for, so the missing directory would surface as a
    # FileNotFoundError from Matplotlib instead of as the error above. The test goes
    # through os.path rather than Path, since Path normalizes a "." away and so would
    # reject it while accepting "..", even though both compose a valid file name.
    prefix = _parameter_validation.str_return_str(prefix, "prefix")
    if prefix != os.path.basename(prefix):
        raise ValueError(
            f'prefix must be a file name component rather than a path, got "{prefix}".'
        )

    resolution_dpi = _parameter_validation.number_in_range_return_float(
        resolution_dpi, "resolution_dpi", 0.0, False
    )

    file_format = _parameter_validation.str_return_str(file_format, "file_format")
    if file_format not in _output_plotting.VALID_FILE_FORMATS:
        valid_formats = ", ".join(f'"{f}"' for f in _output_plotting.VALID_FILE_FORMATS)
        raise ValueError(
            f'file_format must be one of {valid_formats}, got "{file_format}".'
        )

    show_titles = _parameter_validation.boolLike_return_bool(show_titles, "show_titles")
    validated_font_size = None
    if font_size is not None:
        validated_font_size = _parameter_validation.number_in_range_return_float(
            font_size, "font_size", 0.0, False
        )
    validated_line_width = None
    if line_width is not None:
        validated_line_width = _parameter_validation.number_in_range_return_float(
            line_width, "line_width", 0.0, False
        )
    text_rgba = None
    if text_color is not None:
        if not matplotlib.colors.is_color_like(text_color):
            raise ValueError(
                f'text_color must be a color Matplotlib accepts, got "{text_color}".'
            )
        # Matplotlib's type hints accept only tuples of exact lengths, while a color can
        # arrive as any sequence, which is_color_like has just vetted.
        text_rgba = matplotlib.colors.to_rgba(
            cast(matplotlib.typing.ColorType, text_color)
        )

    if not unsteady_solver.ran:
        raise RuntimeError(
            "unsteady_solver must have run before plotting results versus time."
        )

    # The Airplanes' names reach the composed file paths, so they are held to the same
    # rule as prefix: a separator in one would compose a destination outside the
    # directory the caller named. Every Airplane is checked before any figure is drawn,
    # so a bad name never leaves a partly written set behind.
    if save or save_csv:
        for airplane in unsteady_solver.steady_problems[0].airplanes:
            airplane_name_snake = airplane.name.lower().replace(" ", "_")
            if airplane_name_snake != os.path.basename(airplane_name_snake):
                raise ValueError(
                    f'An Airplane\'s name, "{airplane.name}", cannot be used in a file '
                    f"name, since it contains a path separator."
                )

    first_results_step = unsteady_solver.first_results_step

    # Get the time step characteristics. Note that the first time step (time step 0),
    # occurs at 0 seconds.
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

    # Initialize matrices to hold the named loads and named load coefficients at every
    # time step that has results.
    namedForces_W = np.zeros((num_airplanes, 3, num_steps_to_average), dtype=float)
    namedForceCoefficients_W = np.zeros(
        (num_airplanes, 3, num_steps_to_average), dtype=float
    )
    namedMoments_W_Cg = np.zeros((num_airplanes, 3, num_steps_to_average), dtype=float)
    namedMomentCoefficients_W_Cg = np.zeros(
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
            namedForces_W[airplane_id, 0, results_step] = airplane.inducedDrag_W
            namedForces_W[airplane_id, 1, results_step] = airplane.crosswindForce_W
            namedForces_W[airplane_id, 2, results_step] = airplane.lift_W
            namedForceCoefficients_W[airplane_id, 0, results_step] = (
                airplane.inducedDragCoefficient_W
            )
            namedForceCoefficients_W[airplane_id, 1, results_step] = (
                airplane.crosswindForceCoefficient_W
            )
            namedForceCoefficients_W[airplane_id, 2, results_step] = (
                airplane.liftCoefficient_W
            )
            namedMoments_W_Cg[airplane_id, 0, results_step] = (
                airplane.rollingMoment_W_Cg
            )
            namedMoments_W_Cg[airplane_id, 1, results_step] = (
                airplane.pitchingMoment_W_Cg
            )
            namedMoments_W_Cg[airplane_id, 2, results_step] = airplane.yawingMoment_W_Cg
            namedMomentCoefficients_W_Cg[airplane_id, 0, results_step] = (
                airplane.rollingMomentCoefficient_W_Cg
            )
            namedMomentCoefficients_W_Cg[airplane_id, 1, results_step] = (
                airplane.pitchingMomentCoefficient_W_Cg
            )
            namedMomentCoefficients_W_Cg[airplane_id, 2, results_step] = (
                airplane.yawingMomentCoefficient_W_Cg
            )

        results_step += 1

    # Iterate through the Airplane ID's to plot each Airplane's figures.
    for airplane_id in range(num_airplanes):

        # Find and format this Airplane's name for use in the plot titles and file
        # names.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[airplane_id].name
        airplane_name_snake = airplane_name.lower().replace(" ", "_")

        # Compose the stem that each of this Airplane's file names begins with.
        file_stem = airplane_name_snake
        if prefix:
            file_stem = prefix + "_" + airplane_name_snake

        # Plot this Airplane's four load figures.
        _output_plotting.plot_time_history(
            times,
            [
                namedForces_W[airplane_id, 0],
                namedForces_W[airplane_id, 1],
                namedForces_W[airplane_id, 2],
            ],
            FORCE_LABELS,
            [LINEAR_X_COLOR, LINEAR_Y_COLOR, LINEAR_Z_COLOR],
            airplane_name + " Forces",
            WIND_AXES_SUBTITLE,
            FORCE_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_forces." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            times,
            [
                namedForceCoefficients_W[airplane_id, 0],
                namedForceCoefficients_W[airplane_id, 1],
                namedForceCoefficients_W[airplane_id, 2],
            ],
            FORCE_COEFFICIENT_LABELS,
            [LINEAR_X_COLOR, LINEAR_Y_COLOR, LINEAR_Z_COLOR],
            airplane_name + " Force Coefficients",
            WIND_AXES_SUBTITLE,
            FORCE_COEFFICIENT_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_force_coefficients." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            times,
            [
                namedMoments_W_Cg[airplane_id, 0],
                namedMoments_W_Cg[airplane_id, 1],
                namedMoments_W_Cg[airplane_id, 2],
            ],
            MOMENT_LABELS,
            [ANGULAR_X_COLOR, ANGULAR_Y_COLOR, ANGULAR_Z_COLOR],
            airplane_name + " Moments",
            WIND_AXES_CG_SUBTITLE,
            MOMENT_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_moments." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            times,
            [
                namedMomentCoefficients_W_Cg[airplane_id, 0],
                namedMomentCoefficients_W_Cg[airplane_id, 1],
                namedMomentCoefficients_W_Cg[airplane_id, 2],
            ],
            MOMENT_COEFFICIENT_LABELS,
            [ANGULAR_X_COLOR, ANGULAR_Y_COLOR, ANGULAR_Z_COLOR],
            airplane_name + " Moment Coefficients",
            WIND_AXES_CG_SUBTITLE,
            MOMENT_COEFFICIENT_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_moment_coefficients." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )

        # Write this Airplane's twelve plotted load series to one CSV. A reader
        # comparing induced drag against pitching moment should read one table rather
        # than join four, and a table has none of a plot's limited visual capacity, so
        # the four figures collapse to a single file here. The columns repeat the
        # figures' sign conventions and the headers are derived from the same four
        # figures' text, so a row reads the way the figures do.
        if save_csv:
            _output_plotting.write_time_history_csv(
                times,
                _output_plotting.csv_headers(
                    FORCE_LABELS, WIND_AXES_SUBTITLE, FORCE_Y_LABEL
                )
                + _output_plotting.csv_headers(
                    FORCE_COEFFICIENT_LABELS,
                    WIND_AXES_SUBTITLE,
                    FORCE_COEFFICIENT_Y_LABEL,
                )
                + _output_plotting.csv_headers(
                    MOMENT_LABELS, WIND_AXES_CG_SUBTITLE, MOMENT_Y_LABEL
                )
                + _output_plotting.csv_headers(
                    MOMENT_COEFFICIENT_LABELS,
                    WIND_AXES_CG_SUBTITLE,
                    MOMENT_COEFFICIENT_Y_LABEL,
                ),
                [
                    namedForces_W[airplane_id, 0],
                    namedForces_W[airplane_id, 1],
                    namedForces_W[airplane_id, 2],
                    namedForceCoefficients_W[airplane_id, 0],
                    namedForceCoefficients_W[airplane_id, 1],
                    namedForceCoefficients_W[airplane_id, 2],
                    namedMoments_W_Cg[airplane_id, 0],
                    namedMoments_W_Cg[airplane_id, 1],
                    namedMoments_W_Cg[airplane_id, 2],
                    namedMomentCoefficients_W_Cg[airplane_id, 0],
                    namedMomentCoefficients_W_Cg[airplane_id, 1],
                    namedMomentCoefficients_W_Cg[airplane_id, 2],
                ],
                directory / (file_stem + "_loads.csv"),
            )

    # For a free flight solver, also plot the first Airplane's six-degree-of-freedom
    # state history. This is plotted once for the whole simulation, since the state
    # describes the first Airplane, the single rigid body the dynamics integrate.
    if isinstance(
        unsteady_solver,
        _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
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
        angles_E_to_BP1_izyx = np.zeros((3, num_state_steps), dtype=float)
        omegas_BP1__E = np.zeros((3, num_state_steps), dtype=float)
        alphas = np.zeros(num_state_steps, dtype=float)
        betas = np.zeros(num_state_steps, dtype=float)

        # Iterate through the time steps and extract each step's state.
        for step, this_operating_point in enumerate(operating_points):
            positions_E_Eo[:, step] = this_operating_point.CgP1_E_Eo
            velocities_E__E[:, step] = _output_plotting.get_operating_point_velocity(
                this_operating_point
            )
            angles_E_to_BP1_izyx[:, step] = this_operating_point.angles_E_to_BP1_izyx
            omegas_BP1__E[:, step] = this_operating_point.omegas_BP1__E
            alphas[step] = this_operating_point.alpha
            betas[step] = this_operating_point.beta

        # The state describes the first Airplane (the rigid body MuJoCo integrates), so
        # the plot titles and file names use the first Airplane's name.
        airplane_name = unsteady_solver.steady_problems[0].airplanes[0].name
        airplane_name_snake = airplane_name.lower().replace(" ", "_")

        # Compose the stem that each of the state figures' file names begins with.
        file_stem = airplane_name_snake
        if prefix:
            file_stem = prefix + "_" + airplane_name_snake

        _output_plotting.plot_time_history(
            state_times,
            [positions_E_Eo[0], positions_E_Eo[1], positions_E_Eo[2]],
            COMPONENT_LABELS,
            [LINEAR_X_COLOR, LINEAR_Y_COLOR, LINEAR_Z_COLOR],
            airplane_name + " Position",
            POSITION_SUBTITLE,
            POSITION_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_position." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            state_times,
            [velocities_E__E[0], velocities_E__E[1], velocities_E__E[2]],
            COMPONENT_LABELS,
            [LINEAR_X_COLOR, LINEAR_Y_COLOR, LINEAR_Z_COLOR],
            airplane_name + " Velocity",
            VELOCITY_SUBTITLE,
            VELOCITY_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_velocity." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            state_times,
            [
                angles_E_to_BP1_izyx[0],
                angles_E_to_BP1_izyx[1],
                angles_E_to_BP1_izyx[2],
            ],
            ORIENTATION_LABELS,
            [ANGULAR_X_COLOR, ANGULAR_Y_COLOR, ANGULAR_Z_COLOR],
            airplane_name + " Orientation",
            ORIENTATION_SUBTITLE,
            ORIENTATION_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_orientation." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            state_times,
            [omegas_BP1__E[0], omegas_BP1__E[1], omegas_BP1__E[2]],
            ANGULAR_VELOCITY_LABELS,
            [ANGULAR_X_COLOR, ANGULAR_Y_COLOR, ANGULAR_Z_COLOR],
            airplane_name + " Angular Velocity",
            ANGULAR_VELOCITY_SUBTITLE,
            ANGULAR_VELOCITY_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_angular_velocity." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )
        _output_plotting.plot_time_history(
            state_times,
            [alphas, betas],
            AERODYNAMIC_ANGLE_LABELS,
            [ALPHA_COLOR, BETA_COLOR],
            airplane_name + " Aerodynamic Angles",
            AERODYNAMIC_ANGLE_SUBTITLE,
            AERODYNAMIC_ANGLE_Y_LABEL,
            (figure_width_in, figure_height_in),
            save,
            directory / (file_stem + "_aerodynamic_angles." + file_format),
            resolution_dpi,
            show_titles,
            validated_font_size,
            text_rgba,
            validated_line_width,
        )

        # Write the state history to its own CSV rather than into the loads file. The
        # loads begin at the solver's first results step while the state begins at time
        # step 0, so the two have different numbers of rows and cannot share a table.
        # The headers are derived from the same five figures' text, so a column and the
        # figure that plots it describe the quantity the same way.
        if save_csv:
            _output_plotting.write_time_history_csv(
                state_times,
                _output_plotting.csv_headers(
                    COMPONENT_LABELS, POSITION_SUBTITLE, POSITION_Y_LABEL
                )
                + _output_plotting.csv_headers(
                    COMPONENT_LABELS, VELOCITY_SUBTITLE, VELOCITY_Y_LABEL
                )
                + _output_plotting.csv_headers(
                    ORIENTATION_LABELS, ORIENTATION_SUBTITLE, ORIENTATION_Y_LABEL
                )
                + _output_plotting.csv_headers(
                    ANGULAR_VELOCITY_LABELS,
                    ANGULAR_VELOCITY_SUBTITLE,
                    ANGULAR_VELOCITY_Y_LABEL,
                )
                + _output_plotting.csv_headers(
                    AERODYNAMIC_ANGLE_LABELS,
                    AERODYNAMIC_ANGLE_SUBTITLE,
                    AERODYNAMIC_ANGLE_Y_LABEL,
                ),
                [
                    positions_E_Eo[0],
                    positions_E_Eo[1],
                    positions_E_Eo[2],
                    velocities_E__E[0],
                    velocities_E__E[1],
                    velocities_E__E[2],
                    angles_E_to_BP1_izyx[0],
                    angles_E_to_BP1_izyx[1],
                    angles_E_to_BP1_izyx[2],
                    omegas_BP1__E[0],
                    omegas_BP1__E[1],
                    omegas_BP1__E[2],
                    alphas,
                    betas,
                ],
                directory / (file_stem + "_state.csv"),
            )

    # If the user wants to show the plots, do so. This is done outside the loop so that
    # plt.show() is only called once after all figures are created.
    if show:
        plt.show()
    else:
        plt.close("all")


def log_results(
    solver: (
        _steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | _steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
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
            _steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            _steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
        ),
    ):
        these_airplanes = solver.airplanes
        solver_type = "steady"
    elif isinstance(
        solver,
        _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
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

    # Labels for the loads in each Airplane's own geometry axes, followed by the loads
    # in wind axes.
    col1 = [
        "fX_G",
        "fY_G",
        "fZ_G",
        "mX_G_Cg",
        "mY_G_Cg",
        "mZ_G_Cg",
        "cFX_G",
        "cFY_G",
        "cFZ_G",
        "cMX_G_Cg",
        "cMY_G_Cg",
        "cMZ_G_Cg",
        "fX_W",
        "fY_W",
        "fZ_W",
        "mX_W_Cg",
        "mY_W_Cg",
        "mZ_W_Cg",
        "cFX_W",
        "cFY_W",
        "cFZ_W",
        "cMX_W_Cg",
        "cMY_W_Cg",
        "cMZ_W_Cg",
    ]
    col1 = [label + ":" for label in col1]
    col1_space = max(len(elem) for elem in col1) + padding_spaces

    # Named load labels for the wind axes rows only, because names like induced drag and
    # lift are wind axes concepts with no geometry axes counterparts. All four load
    # groups are named the way the figures that plot them are.
    col3 = (
        FORCE_LABELS
        + MOMENT_LABELS
        + FORCE_COEFFICIENT_LABELS
        + MOMENT_COEFFICIENT_LABELS
    )
    col3 = [label + ":" for label in col3]
    col3_space = max(len(elem) for elem in col3) + padding_spaces

    for airplane_num, airplane in enumerate(these_airplanes):
        title_prefix: str = ""
        theseForces_G: np.ndarray = np.empty(0, dtype=float)
        theseMoments_G_Cg: np.ndarray = np.empty(0, dtype=float)
        theseForceCoefficients_G: np.ndarray = np.empty(0, dtype=float)
        theseMomentCoefficients_G_Cg: np.ndarray = np.empty(0, dtype=float)
        theseForces_W: np.ndarray = np.empty(0, dtype=float)
        theseMoments_W_Cg: np.ndarray = np.empty(0, dtype=float)
        theseForceCoefficients_W: np.ndarray = np.empty(0, dtype=float)
        theseMomentCoefficients_W_Cg: np.ndarray = np.empty(0, dtype=float)
        theseNamedForces_W: list[float | None] = []
        theseNamedMoments_W_Cg: list[float | None] = []
        theseNamedForceCoefficients_W: list[float | None] = []
        theseNamedMomentCoefficients_W_Cg: list[float | None] = []

        match solver_type:
            case "steady":
                title_prefix = ""

                _forces_G = airplane.forces_G
                assert _forces_G is not None

                theseForces_G = _forces_G

                _moments_G_Cg = airplane.moments_G_Cg
                assert _moments_G_Cg is not None

                theseMoments_G_Cg = _moments_G_Cg

                _forceCoefficients_G = airplane.forceCoefficients_G
                assert _forceCoefficients_G is not None

                theseForceCoefficients_G = _forceCoefficients_G

                _momentCoefficients_G_Cg = airplane.momentCoefficients_G_Cg
                assert _momentCoefficients_G_Cg is not None

                theseMomentCoefficients_G_Cg = _momentCoefficients_G_Cg

                _forces_W = airplane.forces_W
                assert _forces_W is not None

                theseForces_W = _forces_W

                _moments_W_Cg = airplane.moments_W_Cg
                assert _moments_W_Cg is not None

                theseMoments_W_Cg = _moments_W_Cg

                _forceCoefficients_W = airplane.forceCoefficients_W
                assert _forceCoefficients_W is not None

                theseForceCoefficients_W = _forceCoefficients_W

                _momentCoefficients_W_Cg = airplane.momentCoefficients_W_Cg
                assert _momentCoefficients_W_Cg is not None

                theseMomentCoefficients_W_Cg = _momentCoefficients_W_Cg

                theseNamedForces_W = [
                    airplane.inducedDrag_W,
                    airplane.crosswindForce_W,
                    airplane.lift_W,
                ]
                theseNamedMoments_W_Cg = [
                    airplane.rollingMoment_W_Cg,
                    airplane.pitchingMoment_W_Cg,
                    airplane.yawingMoment_W_Cg,
                ]
                theseNamedForceCoefficients_W = [
                    airplane.inducedDragCoefficient_W,
                    airplane.crosswindForceCoefficient_W,
                    airplane.liftCoefficient_W,
                ]
                theseNamedMomentCoefficients_W_Cg = [
                    airplane.rollingMomentCoefficient_W_Cg,
                    airplane.pitchingMomentCoefficient_W_Cg,
                    airplane.yawingMomentCoefficient_W_Cg,
                ]

            case "static geometry unsteady":
                assert isinstance(
                    solver,
                    _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
                )

                unsteady_problem = solver.unsteady_problem

                title_prefix = "Final "
                theseForces_G = unsteady_problem.finalForces_G[airplane_num]
                theseMoments_G_Cg = unsteady_problem.finalMoments_G_Cg[airplane_num]
                theseForceCoefficients_G = unsteady_problem.finalForceCoefficients_G[
                    airplane_num
                ]
                theseMomentCoefficients_G_Cg = (
                    unsteady_problem.finalMomentCoefficients_G_Cg[airplane_num]
                )
                theseForces_W = unsteady_problem.finalForces_W[airplane_num]
                theseMoments_W_Cg = unsteady_problem.finalMoments_W_Cg[airplane_num]
                theseForceCoefficients_W = unsteady_problem.finalForceCoefficients_W[
                    airplane_num
                ]
                theseMomentCoefficients_W_Cg = (
                    unsteady_problem.finalMomentCoefficients_W_Cg[airplane_num]
                )
                theseNamedForces_W = [
                    unsteady_problem.finalInducedDrags_W[airplane_num],
                    unsteady_problem.finalCrosswindForces_W[airplane_num],
                    unsteady_problem.finalLifts_W[airplane_num],
                ]
                theseNamedMoments_W_Cg = [
                    unsteady_problem.finalRollingMoments_W_Cg[airplane_num],
                    unsteady_problem.finalPitchingMoments_W_Cg[airplane_num],
                    unsteady_problem.finalYawingMoments_W_Cg[airplane_num],
                ]
                theseNamedForceCoefficients_W = [
                    unsteady_problem.finalInducedDragCoefficients_W[airplane_num],
                    unsteady_problem.finalCrosswindForceCoefficients_W[airplane_num],
                    unsteady_problem.finalLiftCoefficients_W[airplane_num],
                ]
                theseNamedMomentCoefficients_W_Cg = [
                    unsteady_problem.finalRollingMomentCoefficients_W_Cg[airplane_num],
                    unsteady_problem.finalPitchingMomentCoefficients_W_Cg[airplane_num],
                    unsteady_problem.finalYawingMomentCoefficients_W_Cg[airplane_num],
                ]
            case "variable geometry unsteady":
                assert isinstance(
                    solver,
                    _unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
                )

                unsteady_problem = solver.unsteady_problem

                title_prefix = "Final Cycle-Averaged "
                theseForces_G = unsteady_problem.finalMeanForces_G[airplane_num]
                theseMoments_G_Cg = unsteady_problem.finalMeanMoments_G_Cg[airplane_num]
                theseForceCoefficients_G = (
                    unsteady_problem.finalMeanForceCoefficients_G[airplane_num]
                )
                theseMomentCoefficients_G_Cg = (
                    unsteady_problem.finalMeanMomentCoefficients_G_Cg[airplane_num]
                )
                theseForces_W = unsteady_problem.finalMeanForces_W[airplane_num]
                theseMoments_W_Cg = unsteady_problem.finalMeanMoments_W_Cg[airplane_num]
                theseForceCoefficients_W = (
                    unsteady_problem.finalMeanForceCoefficients_W[airplane_num]
                )
                theseMomentCoefficients_W_Cg = (
                    unsteady_problem.finalMeanMomentCoefficients_W_Cg[airplane_num]
                )
                theseNamedForces_W = [
                    unsteady_problem.finalMeanInducedDrags_W[airplane_num],
                    unsteady_problem.finalMeanCrosswindForces_W[airplane_num],
                    unsteady_problem.finalMeanLifts_W[airplane_num],
                ]
                theseNamedMoments_W_Cg = [
                    unsteady_problem.finalMeanRollingMoments_W_Cg[airplane_num],
                    unsteady_problem.finalMeanPitchingMoments_W_Cg[airplane_num],
                    unsteady_problem.finalMeanYawingMoments_W_Cg[airplane_num],
                ]
                theseNamedForceCoefficients_W = [
                    unsteady_problem.finalMeanInducedDragCoefficients_W[airplane_num],
                    unsteady_problem.finalMeanCrosswindForceCoefficients_W[
                        airplane_num
                    ],
                    unsteady_problem.finalMeanLiftCoefficients_W[airplane_num],
                ]
                theseNamedMomentCoefficients_W_Cg = [
                    unsteady_problem.finalMeanRollingMomentCoefficients_W_Cg[
                        airplane_num
                    ],
                    unsteady_problem.finalMeanPitchingMomentCoefficients_W_Cg[
                        airplane_num
                    ],
                    unsteady_problem.finalMeanYawingMomentCoefficients_W_Cg[
                        airplane_num
                    ],
                ]
            case _:
                raise ValueError(f"Unknown solver type: {solver_type}")

        # One title per three-row group, in the order the groups are logged: the loads
        # in this Airplane's own geometry axes, then the loads in wind axes.
        titles = [
            _logging.indent(1) + title_prefix + f"Forces {GEOMETRY_AXES_SUBTITLE}:",
            _logging.indent(1) + title_prefix + f"Moments {GEOMETRY_AXES_CG_SUBTITLE}:",
            _logging.indent(1)
            + title_prefix
            + f"Force Coefficients {GEOMETRY_AXES_SUBTITLE}:",
            _logging.indent(1)
            + title_prefix
            + f"Moment Coefficients {GEOMETRY_AXES_CG_SUBTITLE}:",
            _logging.indent(1) + title_prefix + f"Forces {WIND_AXES_SUBTITLE}:",
            _logging.indent(1) + title_prefix + f"Moments {WIND_AXES_CG_SUBTITLE}:",
            _logging.indent(1)
            + title_prefix
            + f"Force Coefficients {WIND_AXES_SUBTITLE}:",
            _logging.indent(1)
            + title_prefix
            + f"Moment Coefficients {WIND_AXES_CG_SUBTITLE}:",
        ]

        col2 = [
            theseForces_G[0],
            theseForces_G[1],
            theseForces_G[2],
            theseMoments_G_Cg[0],
            theseMoments_G_Cg[1],
            theseMoments_G_Cg[2],
            theseForceCoefficients_G[0],
            theseForceCoefficients_G[1],
            theseForceCoefficients_G[2],
            theseMomentCoefficients_G_Cg[0],
            theseMomentCoefficients_G_Cg[1],
            theseMomentCoefficients_G_Cg[2],
            theseForces_W[0],
            theseForces_W[1],
            theseForces_W[2],
            theseMoments_W_Cg[0],
            theseMoments_W_Cg[1],
            theseMoments_W_Cg[2],
            theseForceCoefficients_W[0],
            theseForceCoefficients_W[1],
            theseForceCoefficients_W[2],
            theseMomentCoefficients_W_Cg[0],
            theseMomentCoefficients_W_Cg[1],
            theseMomentCoefficients_W_Cg[2],
        ]
        col2 = [f"{val:#10.3G}" for val in col2]
        col2 = [
            val + " N" if i % 12 < 3 else val + " Nm" if i % 12 < 6 else val
            for i, val in enumerate(col2)
        ]
        col2_space = max(len(elem) for elem in col2) + 2 * padding_spaces

        # Assemble the named load values for the wind axes rows. Each branch above reads
        # these from named load properties, so the mapping from raw components to named
        # loads lives with the loads classes rather than here.
        col4_values = (
            theseNamedForces_W
            + theseNamedMoments_W_Cg
            + theseNamedForceCoefficients_W
            + theseNamedMomentCoefficients_W_Cg
        )
        col4 = [f"{val:#10.3G}" for val in col4_values]
        col4 = [
            val + " N" if i < 3 else val + " Nm" if i < 6 else val
            for i, val in enumerate(col4)
        ]

        logger.info(_logging.indent() + f'Airplane "{airplane.name}":')

        # Display the Reynolds number for steady solvers.
        if solver_type == "steady":
            assert isinstance(
                solver,
                (
                    _steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
                    _steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
                ),
            )
            re = solver.reynolds_numbers[airplane_num]
            logger.info(_logging.indent(1) + f"Reynolds Number: {re:#.3G}")

        for i in range(len(col1)):
            if i % 3 == 0:
                logger.info(titles[i // 3])

            if i < 12:
                # The geometry axes rows have no named load columns.
                s = _logging.indent(2) + f"{col1[i]:<{col1_space}}{col2[i]}"
            else:
                j = i - 12
                s = (
                    _logging.indent(2)
                    + f"{col1[i]:<{col1_space}}{col2[i]:<{col2_space}}{col3[j]:<{col3_space}}{col4[j]}"
                )
            logger.info(s)

    # For a free flight solver, also log the first Airplane's initial and final
    # six-degree-of-freedom state. This is logged once, since the state describes the
    # first Airplane, the single rigid body the dynamics integrate.
    if isinstance(
        solver,
        _free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    ):
        # Narrow the movement to a FreeFlightMovement so its operating point history is
        # typed. The solver type guarantees this, so the assert documents the invariant.
        movement = solver.unsteady_problem.movement
        assert isinstance(movement, free_flight_movement_mod.FreeFlightMovement)
        operating_points = movement.operating_point_movement.operating_points

        final_time = solver.delta_time * (len(operating_points) - 1)

        logger.info(
            _logging.indent() + "The First Airplane's Free Flight State History:"
        )

        # Each vector state quantity is broken into one row per component, mirroring the
        # per-component force and moment rows above, and each component row is labeled
        # with its variable-convention name. The four group headers are logged before
        # their first component (at flat indices 0, 3, 6, and 9), and the component
        # labels are padded to a common width so the values align.
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
            _logging.indent(2) + f"Position {POSITION_SUBTITLE}:"
        )
        state_group_header_orientation = (
            _logging.indent(2) + f"Orientation {ORIENTATION_SUBTITLE}:"
        )
        state_group_header_velocity = (
            _logging.indent(2) + f"Velocity {VELOCITY_SUBTITLE}:"
        )
        state_group_header_angular_velocity = (
            _logging.indent(2) + f"Angular Velocity {ANGULAR_VELOCITY_SUBTITLE}:"
        )

        # Log the initial state (at time step 0) and the final state.
        for state_label, state_time, this_operating_point in [
            ("Initial State", 0.0, operating_points[0]),
            ("Final State", final_time, operating_points[-1]),
        ]:
            CgP1_E_Eo = this_operating_point.CgP1_E_Eo
            angles_E_to_BP1_izyx = this_operating_point.angles_E_to_BP1_izyx
            velocity_E__E = _output_plotting.get_operating_point_velocity(
                this_operating_point
            )
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
                f"{value:#10.3G} {unit}"
                for value, unit in zip(state_component_values, state_component_units)
            ]

            logger.info(
                _logging.indent(1) + f"{state_label} (at t = {state_time:#.3G} s):"
            )

            for i in range(len(state_component_labels)):
                if i == 0:
                    logger.info(state_group_header_position)
                elif i == 3:
                    logger.info(state_group_header_orientation)
                elif i == 6:
                    logger.info(state_group_header_velocity)
                elif i == 9:
                    logger.info(state_group_header_angular_velocity)

                logger.info(
                    _logging.indent(3)
                    + f"{state_component_labels[i]:<{state_component_space}}"
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
            logger.info(
                _logging.indent(2) + f"{alpha_label:<{aerodynamic_angle_space}}"
                f"{this_operating_point.alpha:#10.3G} deg"
            )
            logger.info(
                _logging.indent(2) + f"{beta_label:<{aerodynamic_angle_space}}"
                f"{this_operating_point.beta:#10.3G} deg"
            )
