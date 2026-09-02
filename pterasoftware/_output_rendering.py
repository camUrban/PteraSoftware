"""Contains the styling, geometry building, and scene assembly used to render the
visualizations with PyVista."""

from __future__ import annotations

import math
import queue
import threading
from pathlib import Path
from typing import NamedTuple, cast

import matplotlib.colors
import numpy as np
import pyvista as pv
import webp

from . import (
    _colormaps,
    _logging,
    _mujoco_model,
    _private_access,
    _transformations,
    free_flight_unsteady_ring_vortex_lattice_method,
    geometry,
)
from . import operating_point as operating_point_mod
from . import (
    problems,
    unsteady_ring_vortex_lattice_method,
)

_logger = _logging.get_logger("output")

# Define the colors, sizes, and positions used when rendering the geometry. The color
# maps and color palettes live in the _colormaps module. The edge line widths are in
# pixels, so they are tuned for REFERENCE_WINDOW_SIZE and scaled by get_window_scale
# wherever they are used, as the font sizes below are.
_WAKE_VORTEX_COLOR = "white"
_WAKE_VORTEX_EDGE_LINE_WIDTH = 1.0
_PANEL_COLOR = "chartreuse"
_PANEL_EDGE_LINE_WIDTH = 1.0
_IMAGE_SURFACE_SCALE = 5.0
IMAGE_REFLECTION_MUTE_FACTOR = 0.5
_IMAGE_SURFACE_CHECKER_SIZE = 25
_IMAGE_SURFACE_COLOR_A = np.array([40, 40, 40], dtype=np.uint8)
_IMAGE_SURFACE_COLOR_B = np.array([80, 80, 80], dtype=np.uint8)
TEXT_COLOR = (129, 129, 129)
TEXT_COLOR_SURFACE = (220, 220, 220)
PLOTTER_BACKGROUND_COLOR = "black"

# Define the valid scalar types for coloring Panels.
VALID_SCALAR_TYPES = ("induced drag", "side force", "lift")

# Define the shading coefficients for the MuJoCo geom actors, which are the only lit
# actors in any scene. Every other actor is added with lighting disabled, so the
# headlight that add_mujoco_geometry brings into a scene shades the geoms alone, and a
# scene without MuJoCo geoms never carries a light at all. An actor with lighting
# disabled renders identically whether or not the scene holds a light, so the geoms'
# shading leaves every other actor's appearance untouched. The ambient floor keeps the
# faces the headlight misses from going black.
_MUJOCO_GEOMETRY_AMBIENT = 0.3
_MUJOCO_GEOMETRY_DIFFUSE = 0.7

# Define the window length used to measure the largest window a window manager will
# grant. X11 window dimensions are 16 bit, so this is the largest a window can ask to be
# and is certain to be shrunk to the maximum on any display.
_OVERSIZED_WINDOW_LENGTH = 32767

# Set constants for the color maps, scalar bars, and text boxes. The positions are
# fractions of the render window, so they track its size on their own. The font sizes
# are in pixels and do not, so they are tuned for REFERENCE_WINDOW_SIZE and scaled by
# get_window_scale wherever they are used.
_COLOR_MAP_NUM_SIG = 3
_BAR_TITLE_FONT_SIZE = 30
_BAR_LABEL_FONT_SIZE = 21
_BAR_WIDTH = 0.5
_BAR_POSITION_X = 0.25
_BAR_POSITION_Y = 0.05
_BAR_N_LABELS = 2
# The maximum and minimum scalar labels are right justified, so these positions anchor
# their right edges against the right margin and the text grows leftward. No label
# length can overflow the window that way.
_TEXT_MAX_POSITION = (0.99, 0.080)
_TEXT_MIN_POSITION = (0.99, 0.045)
TEXT_FONT_SIZE = 10

# An animation's playback overlays sit against the left margin, a step apart, mirroring
# the maximum and minimum scalar text on the other side of the window. The scalar bar
# occupies the lower of the two rows as well, and its widget's box hides anything drawn
# behind it, so both overlays are held close to the margin to leave the lower one as
# much room as _BAR_POSITION_X allows.
_TEXT_SPEED_POSITION = (0.01, 0.080)
_TEXT_DROPPED_FRAMES_POSITION = (0.01, 0.045)

# Define the render window size that every font size and line width in the
# visualizations is tuned against.
REFERENCE_WINDOW_SIZE = (1024, 768)

# Define the highest frame rate, in frames per second, that an animation is saved at.
# Some programs will not render a WebP any faster than this.
_MAX_FRAME_RATE = 50.0

# Define the libwebp compression method that every saved WebP is encoded with, from 0
# (the fastest) to 6 (the slowest). The method sets how hard the encoder works to shrink
# the file at a given quality, so it trades encode time for file size and leaves the
# image quality alone. The fastest method encodes about twice as fast as the default of
# 4 for a file about a fifth larger, which is the better side of that trade for
# animations, whose save time is otherwise dominated by the encode.
WEBP_METHOD = 0

# Define the number of captured frames an AnimationWriter holds while they wait to be
# encoded. When frames are captured faster than they are encoded, the capture loop
# blocks once this many are waiting, which bounds the raw frames in memory no matter how
# many time steps an animation has. Each slot costs one raw frame of memory, and the
# slower of the two stages sets the pace whatever the depth, so a deeper queue buys
# nothing. Two slots let one frame wait while another is encoded.
_ANIMATION_WRITER_QUEUE_DEPTH = 2

# Define the fewest frames per cycle of the fastest prescribed motion that a saved
# animation can show before it stops resolving that motion. The Nyquist floor is 2.0
# frames per cycle, below which the apparent motion can reverse direction or stand
# still, but an animation is already jerky and misleading well above that floor. One
# threshold here covers both regimes.
_MIN_FRAMES_PER_PERIOD = 10.0


def get_window_scale(window_width: int, window_height: int) -> float:
    """Returns the factor by which to scale the font sizes and line widths that are
    tuned for REFERENCE_WINDOW_SIZE.

    VTK sizes text and lines in pixels, so without scaling they keep the same pixel
    count at every window size. That crowds a small window, where the scalar bar's
    labels overlap the bar itself, and it leaves them nearly unreadable in a large one.

    The smaller of the two ratios is used rather than the height ratio alone. The scalar
    labels grow leftward from the right margin toward the centered scalar bar, so a wide
    but short window runs out of horizontal room before it runs out of vertical room.

    :param window_width: The render window's width in pixels.
    :param window_height: The render window's height in pixels.
    :return: The scale factor, which is exactly 1.0 at REFERENCE_WINDOW_SIZE.
    """
    return min(
        window_width / REFERENCE_WINDOW_SIZE[0],
        window_height / REFERENCE_WINDOW_SIZE[1],
    )


def get_largest_window_size() -> tuple[int, int]:
    """Returns the largest on-screen render window a window manager will grant.

    A window cannot be larger than the area the window manager grants, which is the
    display less any docks or bars and less the window's own title bar, and VTK silently
    shrinks one that asks for more. GetScreenSize reports neither reduction, and the
    granted size only becomes readable once a window has been realized, so this measures
    it rather than computing it.

    The measurement is taken on a throwaway render window asking to be far larger than
    any display, which is shrunk in both dimensions and so reports the maximum in both.
    draw and animate detect a shrink on their own render window and call this only to
    name the ceiling in the resulting error, so the throwaway window is created on the
    error path alone.

    :return: A tuple of the largest grantable width and height in pixels. Off-screen
        rendering has no such ceiling, so this reports the requested size there.
    """
    probe = pv.Plotter(
        window_size=[_OVERSIZED_WINDOW_LENGTH, _OVERSIZED_WINDOW_LENGTH], lighting=None
    )
    try:
        # Match the background the visualizations render on so the probe, which is
        # briefly visible, does not flash a bright window against a dark one.
        probe.set_background(color=PLOTTER_BACKGROUND_COLOR)  # type: ignore[call-arg]
        probe_render_window = probe.ren_win
        assert probe_render_window is not None
        probe_render_window.Render()
        granted_width, granted_height = probe_render_window.GetSize()
    finally:
        probe.close()
    return granted_width, granted_height


class Playback(NamedTuple):
    """How a saved animation steps through a simulation's time steps and how fast it
    plays them back.

    keep_every is a whole number of time steps, so an animation saves the time steps at
    its multiples and drops the rest. frame_rate is deliberately not a whole number of
    frames per second: a WebP stores each frame's duration as a whole number of
    milliseconds, but the encoder reaches a frame rate by giving each frame a cumulative
    timestamp rather than a fixed duration, so a fractional rate is played back to
    within a small fraction of a percent instead of being quantized.
    """

    keep_every: int
    frame_rate: float
    overlay_texts: list[tuple[str, tuple[float, float]]]


def resolve_playback(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    speed: float | None,
    save: bool,
) -> Playback:
    """Resolves a requested playback speed into a frame stride, a frame rate, and the
    text overlays that describe them.

    A speed that the maximum frame rate cannot carry on its own is reached by saving
    only every keep_every-th frame, which trades the animation's temporal resolution for
    its speed. The stride is 1 whenever no frames need to be dropped, and it is also
    what holds the frame rate at or below the maximum.

    Dropping frames undersamples the prescribed motion, so this warns when the frames
    that survive no longer resolve the fastest motion in the simulation. The shortest
    non zero period is what governs that, since aliasing is set by the highest frequency
    present rather than by the composite motion's repeat interval.

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver being animated.
    :param speed: The requested playback speed as a multiple of real time, or None to
        play at true speed and scale that speed down only when true speed would need
        more than the maximum frame rate. It must be positive.
    :param save: Whether the animation is being saved. Frames are dropped from the saved
        file alone, so the aliasing warning is skipped when nothing is saved.
    :return: The resolved Playback.
    :raises ValueError: If the speed is so slow that the animation would save fewer than
        one frame per second of playback, or so fast that it would save fewer than two
        frames in total.
    """
    delta_time = unsteady_solver.delta_time
    num_steps = unsteady_solver.num_steps
    if speed is None:
        speed = min(1.0, _MAX_FRAME_RATE * delta_time)

    # The clamp only guards against the division landing an ulp above the maximum, since
    # the stride already holds the quotient at or below it.
    keep_every = math.ceil(speed / (_MAX_FRAME_RATE * delta_time))
    frame_rate = min(_MAX_FRAME_RATE, speed / (keep_every * delta_time))
    achieved_speed = frame_rate * keep_every * delta_time

    # Frames are saved at every time step that is a multiple of the stride, so time step
    # 0 is always saved and the last time step may not be.
    num_saved_frames = (num_steps - 1) // keep_every + 1

    # Guard the two speeds that cannot produce a usable animation. The default reaches
    # neither, except on a simulation whose delta_time exceeds one second, where playing
    # at true speed genuinely does need less than one frame per second.
    if frame_rate < 1.0:
        raise ValueError(
            f"speed is too slow to animate a simulation whose delta_time is "
            f"{delta_time} seconds, because it saves fewer than one frame per second "
            f"of playback. The slowest speed this simulation can be animated at is "
            f"{delta_time}."
        )
    if num_saved_frames < 2:
        largest_speed = (num_steps - 1) * _MAX_FRAME_RATE * delta_time
        raise ValueError(
            f"speed is too fast to animate a simulation of {num_steps} time steps, "
            f"because it saves fewer than two frames. The fastest speed this "
            f"simulation can be animated at is {largest_speed}."
        )

    # The speed overlay reports the speed the animation achieves, which is the requested
    # speed wherever the clamp above leaves the frame rate alone. The second overlay
    # appears only when frames are actually dropped.
    overlay_texts = [(f"Speed: {100 * achieved_speed:#.4G}%", _TEXT_SPEED_POSITION)]
    if keep_every > 1:
        overlay_texts.append(
            (f"Frames: Every 1 of {keep_every}", _TEXT_DROPPED_FRAMES_POSITION)
        )

    # A static geometry has nothing to alias, a stride of 1 drops nothing, and an
    # animation that is not saved keeps every frame, so none of the three is tested.
    movement = unsteady_solver.unsteady_problem.movement
    if save and keep_every > 1 and not movement.static:
        frames_per_period = movement.min_period / (keep_every * delta_time)
        if frames_per_period < _MIN_FRAMES_PER_PERIOD:
            _logger.warning(
                _logging.indent() + f"The animation saves {frames_per_period:.1f} "
                f"frames per cycle of its fastest motion, whose period is "
                f"{movement.min_period} seconds. Below {_MIN_FRAMES_PER_PERIOD} frames "
                f"per cycle the motion looks jerky, and below the Nyquist floor of 2.0 "
                f"it can appear to reverse direction or stand still. Animate at a "
                f"slower speed to save more frames."
            )

    return Playback(keep_every, frame_rate, overlay_texts)


def add_playback_overlays(
    plotter: pv.Plotter,
    playback: Playback,
    window_scale: float,
    text_color: tuple[int, int, int],
) -> None:
    """Adds a saved animation's playback text overlays to a Plotter.

    :param plotter: The Plotter to add the overlays to.
    :param playback: The Playback whose overlays will be added.
    :param window_scale: The render window's size relative to REFERENCE_WINDOW_SIZE,
        which the font size is scaled by.
    :param text_color: The RGB color to draw the overlays in.
    :return: None
    """
    for overlay_text, overlay_position in playback.overlay_texts:
        plotter.add_text(
            text=overlay_text,
            position=overlay_position,
            font_size=round(TEXT_FONT_SIZE * window_scale),
            viewport=True,
            color=text_color,
            render=False,
        )


def get_panel_surfaces(
    airplanes: tuple[geometry.airplane.Airplane, ...],
) -> pv.PolyData:
    """Returns a PolyData representation of the Wings' Panels' surfaces associated with
    all the Airplanes in a tuple of Airplanes.

    :param airplanes: The tuple of Airplanes whose Wings' Panels' surfaces will be
        returned.
    :return: A PolyData representation of the Airplanes' Wings' Panels' surfaces.
    """
    # Gather every Wing's Panels in order, so the vertex and face ndarrays can be sized
    # once rather than grown by stacking a Panel at a time, which is called once per
    # time step of an animation.
    wing_panels = []
    for airplane in airplanes:
        for wing in airplane.wings:
            _panels = wing.panels
            assert _panels is not None
            wing_panels.append(np.ravel(_panels))
    num_panels = sum(panels.size for panels in wing_panels)

    # Fill each Panel's four vertices, wound front left, front right, back right, and
    # back left so the cell traces the Panel's outline.
    panel_vertices = np.empty((num_panels * 4, 3), dtype=float)
    panel_num = 0
    for panels in wing_panels:
        for panel in panels:
            base_vertex = panel_num * 4
            panel_vertices[base_vertex] = panel.Flpp_GP1_CgP1
            panel_vertices[base_vertex + 1] = panel.Frpp_GP1_CgP1
            panel_vertices[base_vertex + 2] = panel.Brpp_GP1_CgP1
            panel_vertices[base_vertex + 3] = panel.Blpp_GP1_CgP1
            panel_num += 1

    # Return the Panels' surfaces.
    return pv.PolyData(panel_vertices, _get_quadrilateral_faces(num_panels))


def _get_quadrilateral_faces(num_quadrilaterals: int) -> np.ndarray:
    """Returns the faces ndarray of a PolyData made of consecutive quadrilaterals.

    The faces are in PolyData's padded form, where each face is its vertex count
    followed by its vertex indices. Face i covers vertices 4i through 4i + 3, so this is
    the faces ndarray for any mesh whose vertices are stored four to a cell in cell
    order.

    :param num_quadrilaterals: The number of quadrilateral faces.
    :return: A (num_quadrilaterals * 5,) ndarray of ints representing the faces.
    """
    base_vertices = np.arange(num_quadrilaterals, dtype=int) * 4
    faces = np.empty((num_quadrilaterals, 5), dtype=int)
    faces[:, 0] = 4
    faces[:, 1] = base_vertices
    faces[:, 2] = base_vertices + 1
    faces[:, 3] = base_vertices + 2
    faces[:, 4] = base_vertices + 3
    return np.ravel(faces)


def get_streamline_surfaces(
    gridStreamlinePoints_GP1_CgP1: np.ndarray,
) -> pv.PolyData:
    """Returns a PolyData representation of a solver's streamlines.

    Each streamline becomes one polyline through its own points, so the whole set is a
    single mesh. Adding one mesh per segment instead, which is what this replaced, gives
    VTK an actor per segment and costs seconds once there are a few hundred of them.

    :param gridStreamlinePoints_GP1_CgP1: A (num_points,num_streamlines,3) ndarray of
        floats representing the points along each streamline (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).
    :return: A PolyData representation of the streamlines.
    """
    num_points, num_streamlines = gridStreamlinePoints_GP1_CgP1.shape[:2]

    # Order the points streamline by streamline so that each one is a contiguous run.
    streamline_points = np.swapaxes(gridStreamlinePoints_GP1_CgP1, 0, 1).reshape(-1, 3)

    # Describe each streamline as one polyline over its run, in the flat cell format VTK
    # expects: a point count followed by that many point indices.
    first_indices = np.arange(num_streamlines) * num_points
    point_indices = first_indices[:, None] + np.arange(num_points)[None, :]
    lines = np.hstack(
        (np.full((num_streamlines, 1), num_points, dtype=int), point_indices)
    ).ravel()

    return pv.PolyData(streamline_points, lines=lines)


def get_image_surface_mesh_and_texture(
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
    plane_size = _IMAGE_SURFACE_SCALE * bbox_diagonal

    mesh = pv.Plane(
        center=plane_center,
        direction=surface_normal,
        i_size=plane_size,
        j_size=plane_size,
    )

    # Build a checkerboard texture image. Each cell is one pixel, so a 25 x 25
    # checkerboard is a 25 x 25 x 3 RGB image.
    n = _IMAGE_SURFACE_CHECKER_SIZE
    row = np.arange(n, dtype=int)
    col = np.arange(n, dtype=int)
    rr, cc = np.meshgrid(row, col, indexing="ij")
    is_dark = (rr + cc) % 2 == 0
    image = np.where(
        is_dark[:, :, np.newaxis], _IMAGE_SURFACE_COLOR_A, _IMAGE_SURFACE_COLOR_B
    )
    texture = pv.numpy_to_texture(image)

    return mesh, texture


def get_animation_image_surface(
    unsteady_solver: unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    step_airplanes: list[tuple[geometry.airplane.Airplane, ...]],
    step_transforms: list[np.ndarray],
    is_free_flight: bool,
    show_wake_vortices: bool,
) -> tuple[
    pv.PolyData | None,
    pv.Texture | None,
    np.ndarray | None,
    tuple[float, float, float, float, float, float] | None,
]:
    """Builds the image surface plane that an animation reuses for every frame, and
    returns it alongside the reflection that maps geometry onto its far side.

    The plane is built from the last time step's geometry so that it is large enough to
    encompass the fully developed wake. It is static across the animation, as is the
    reflection, so both are computed once here rather than per frame.

    For the standard body-fixed rendering, the plane is sized to the last step's
    geometry together with its reflected copy, and the wake is included in that extent
    when it is being shown. For free flight, the plane is fixed in the world while the
    geometry is rendered in Earth axes, so it is built from the last step's geometry
    axis bounds and then mapped into Earth axes. The reflection is re-expressed in Earth
    axes by the same passive transformation, which lets each frame's reflected geometry
    be built directly from the panels already mapped into Earth axes.

    :param unsteady_solver: The UnsteadyRingVortexLatticeMethodSolver being animated.
    :param step_airplanes: A list holding each time step's tuple of Airplanes.
    :param step_transforms: A list holding, per time step, the (4,4) ndarray of floats
        representing the passive transformation from the first Airplane's geometry axes
        (relative to its CG) to Earth axes (relative to the Earth origin). It is empty
        when is_free_flight is False.
    :param is_free_flight: Set this to True when the solver is a
        FreeFlightUnsteadyRingVortexLatticeMethodSolver, whose geometry is rendered in
        Earth axes.
    :param show_wake_vortices: Set this to True when the wake ring vortices are being
        shown, which includes them in the extent the plane is sized to.
    :return: A tuple holding the image surface plane's mesh, its Texture, the (4,4)
        ndarray of floats representing the active transformation that reflects geometry
        across the image surface, and the (xmin, xmax, ymin, ymax, zmin, zmax) bounding
        box to fit the camera to. The bounding box spans the geometry and its reflection
        but not the much larger plane, so fitting to it keeps the plane from dominating
        the view. All four are None when no image surface is defined, and the bounding
        box is also None for free flight, which frames its own camera to the trajectory.
    """
    last_step = len(step_airplanes) - 1
    last_step_operating_point = unsteady_solver.steady_problems[
        last_step
    ].operating_point
    T_reflect = last_step_operating_point.surfaceReflect_T_act_GP1_CgP1

    if T_reflect is None:
        return None, None, None, None

    last_step_panel_surfaces = get_panel_surfaces(step_airplanes[last_step])

    if is_free_flight:
        # Size the plane from the geometry axis bounds, since the image surface helper
        # builds it from geometry axis quantities, then map it into Earth axes to match
        # the rendered geometry.
        image_surface_result = get_image_surface_mesh_and_texture(
            last_step_operating_point, last_step_panel_surfaces.bounds
        )
        assert image_surface_result is not None
        image_surface_mesh, image_surface_texture = image_surface_result

        T_pas_last = step_transforms[last_step]
        image_surface_mesh = transform_mesh(image_surface_mesh, T_pas_last)
        T_reflect = T_pas_last @ T_reflect @ _transformations.invert_T_pas(T_pas_last)

        # Free flight frames its camera to the whole trajectory, so it needs no geometry
        # bounding box from here.
        return image_surface_mesh, image_surface_texture, T_reflect, None

    # Size the plane to the last step's geometry together with its reflected copy, so
    # that it spans both sides of the surface.
    reflected_last_step_panel_surfaces = _reflect_mesh(
        last_step_panel_surfaces, T_reflect
    )
    if show_wake_vortices:
        last_step_wake_surfaces = get_wake_ring_vortex_surfaces(
            unsteady_solver, last_step
        )
        combined = (
            last_step_panel_surfaces.merge(last_step_wake_surfaces)
            .merge(reflected_last_step_panel_surfaces)
            .merge(_reflect_mesh(last_step_wake_surfaces, T_reflect))
        )
    else:
        combined = last_step_panel_surfaces.merge(reflected_last_step_panel_surfaces)

    image_surface_geometry_bounds = combined.bounds
    image_surface_result = get_image_surface_mesh_and_texture(
        last_step_operating_point, image_surface_geometry_bounds
    )
    assert image_surface_result is not None
    image_surface_mesh, image_surface_texture = image_surface_result

    return (
        image_surface_mesh,
        image_surface_texture,
        T_reflect,
        image_surface_geometry_bounds,
    )


def get_free_flight_transformation(
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


def get_free_flight_body_transformation(
    this_operating_point: operating_point_mod.OperatingPoint,
) -> np.ndarray:
    """Returns the passive transformation from the first Airplane's body axes, relative
    to the first Airplane's CG, to Earth axes, relative to the Earth origin.

    Body-attached MuJoCo geoms are rigid in the first Airplane's body axes, so posing
    them for a time step chains the constant body axes to geometry axes rotation with
    the same geometry axes to Earth axes transformation that the Panel meshes are mapped
    through.

    :param this_operating_point: The OperatingPoint whose Earth-frame position and
        orientation define the transformation.
    :return: A (4,4) ndarray of floats representing the passive transformation from the
        first Airplane's body axes, relative to the first Airplane's CG, to Earth axes,
        relative to the Earth origin.
    """
    return _transformations.compose_T_pas(
        this_operating_point.T_pas_BP1_CgP1_to_GP1_CgP1,
        get_free_flight_transformation(this_operating_point),
    )


def get_mujoco_render_geometry(
    free_flight_solver: free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
) -> tuple[list[_mujoco_model.RenderGeom], list[_mujoco_model.RenderGeom]]:
    """Returns a free flight solver's MuJoCo render geometry, split into worldbody geoms
    and body geoms.

    The geometry comes from the compiled MuJoCo model held by the solver's
    FreeFlightUnsteadyProblem, reached through the registration pattern in
    _private_access.

    :param free_flight_solver: The FreeFlightUnsteadyRingVortexLatticeMethodSolver whose
        MuJoCo geometry will be returned.
    :return: A tuple of two lists of RenderGeoms. The first holds the worldbody geoms,
        whose meshes are in Earth axes, relative to the Earth origin. The second holds
        the body geoms, whose meshes are in the first Airplane's body axes, relative to
        the first Airplane's CG.
    """
    # The solver stores its problem widened to the core type, and its initialization
    # validates that the problem is a FreeFlightUnsteadyProblem, so the cast narrows
    # without a runtime check.
    free_flight_unsteady_problem = cast(
        problems.FreeFlightUnsteadyProblem, free_flight_solver.unsteady_problem
    )
    extracted_render_geoms = _private_access.get_mujoco_model(
        free_flight_unsteady_problem
    ).get_render_geometry()

    # Compute each geom's shading normals here, once, splitting the sharp edges so the
    # shading stays crisp across creases. The split duplicates the points along edges
    # sharper than PyVista's feature angle, letting each face shade with its own normal,
    # where plain smooth shading would average one normal across the crease and smear
    # it. The geoms are rigid, so transform_mesh carries the normals along with the
    # points, and add_mujoco_geometry hands PyVista a mesh that already has them rather
    # than having it split the edges again for every frame of an animation, which costs
    # more than the rest of the frame put together for a detailed mesh. The split leaves
    # behind a point id array that nothing reads, so it is dropped.
    render_geoms = []
    for render_geom in extracted_render_geoms:
        mesh = render_geom.mesh.compute_normals(
            cell_normals=False,
            split_vertices=True,
            feature_angle=pv.global_theme.sharp_edges_feature_angle,
        )
        if "pyvistaOriginalPointIds" in mesh.point_data:
            del mesh.point_data["pyvistaOriginalPointIds"]
        render_geoms.append(
            _mujoco_model.RenderGeom(
                mesh=mesh,
                rgba=render_geom.rgba,
                body_attached=render_geom.body_attached,
            )
        )

    worldbody_geoms = [
        render_geom for render_geom in render_geoms if not render_geom.body_attached
    ]
    body_geoms = [
        render_geom for render_geom in render_geoms if render_geom.body_attached
    ]
    return worldbody_geoms, body_geoms


def transform_mesh(
    mesh: pv.PolyData,
    T_pas: np.ndarray,
) -> pv.PolyData:
    """Returns a copy of a PolyData mesh with its points mapped through a passive
    transformation.

    :param mesh: The PolyData mesh to transform.
    :param T_pas: A (4,4) ndarray of floats representing the passive transformation to
        apply to the mesh's points.
    :return: A new PolyData mesh with all points mapped through the transformation,
        along with its active point normals, if it has any, mapped as directions.
    """
    transformed = mesh.copy()
    transformed.points = _transformations.apply_T_to_vectors(
        T_pas,
        mesh.points,
        is_position=True,
    )

    # Carry any shading normals along as directions, so they stay perpendicular to the
    # faces they shade. Under a reflection, this maps an outward normal to the reflected
    # face's outward normal, which is what the shading needs.
    normals_name = mesh.point_data.active_normals_name
    if normals_name is not None:
        transformed.point_data[normals_name] = _transformations.apply_T_to_vectors(
            T_pas,
            mesh.point_data[normals_name],
            is_position=False,
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
    return transform_mesh(mesh, T_reflect)


def get_free_flight_fit_parallel_scale(
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
    # Build the camera's screen right and up axes in Earth axes. The camera looks from
    # its position back toward the focal point, i.e. along the negative view direction.
    # The up axis is the supplied up made orthogonal to that look direction.
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
    # array of corners times a (3,) screen axis gives an (N,) array whose entries are
    # each corner's signed distance from the focal point along that axis (the axes are
    # unit vectors, so the dot product is the scalar projection). Take the largest
    # magnitude along either axis so neither screen dimension is clipped.
    half_extent_right = float(np.abs(corners_E @ rightDirection_E).max())
    half_extent_up = float(np.abs(corners_E @ upDirection_E).max())
    return margin * max(half_extent_right, half_extent_up)


def mute_color(
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


def mute_colormap(
    cmap: matplotlib.colors.Colormap,
    factor: float,
) -> matplotlib.colors.ListedColormap:
    """Returns a muted version of a colormap by linearly interpolating each color toward
    middle gray.

    :param cmap: The colormap to mute.
    :param factor: The muting factor in [0, 1]. 0 means no change, 1 means fully gray.
    :return: A ListedColormap with muted colors.
    """
    colors = cmap(np.linspace(0, 1, 256))
    gray = 0.5
    colors[:, :3] = colors[:, :3] + factor * (gray - colors[:, :3])
    return matplotlib.colors.ListedColormap(colors)


def get_wake_ring_vortex_surfaces(
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

    # Interleave the four corner stacks so each wake ring vortex's vertices are wound
    # front left, front right, back right, and back left. Every corner stack is a
    # (num_wake_ring_vortices, 3) ndarray, so four strided assignments build the whole
    # mesh at once rather than a wake ring vortex at a time, which matters because this
    # is called once per time step of an animation.
    wake_ring_vortex_vertices = np.empty((num_wake_ring_vortices * 4, 3), dtype=float)
    wake_ring_vortex_vertices[0::4] = stackFlwrvp_GP1_CgP1[:num_wake_ring_vortices]
    wake_ring_vortex_vertices[1::4] = stackFrwrvp_GP1_CgP1[:num_wake_ring_vortices]
    wake_ring_vortex_vertices[2::4] = stackBrwrvp_GP1_CgP1[:num_wake_ring_vortices]
    wake_ring_vortex_vertices[3::4] = stackBlwrvp_GP1_CgP1[:num_wake_ring_vortices]

    # Return the wake ring vortex surfaces.
    return pv.PolyData(
        wake_ring_vortex_vertices, _get_quadrilateral_faces(num_wake_ring_vortices)
    )


def get_scalars(
    airplanes: tuple[geometry.airplane.Airplane, ...],
    scalar_type: str,
    qInf__E: float,
) -> np.ndarray:
    """Returns the load coefficient values from a SteadyProblem's Airplanes' Wings'
    Panels.

    The scalar type is taken as given rather than checked, because the public output
    functions reject an unrecognized one before any geometry is walked, and this is
    called once per time step of an animation. A type outside the three named ones
    therefore contributes no scalars rather than raising.

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

    # Map the scalar type string to the corresponding Panel named force attribute.
    panel_force_attributes = {
        "induced drag": "inducedDrag_W",
        "side force": "sideForce_W",
        "lift": "lift_W",
    }

    attribute_name = panel_force_attributes.get(scalar_type)
    if attribute_name is None:
        return scalars

    # Iterate through the Airplanes' Wings.
    for airplane in airplanes:
        for wing in airplane.wings:
            _panels = wing.panels
            assert _panels is not None

            # Unravel this Wing's ndarray of Panels and iterate through them.
            these_panels = np.ravel(_panels)
            for this_panel in these_panels:
                force = getattr(this_panel, attribute_name)
                assert force is not None
                coefficient = force / qInf__E / this_panel.area
                scalars = np.hstack((scalars, coefficient))

    # Return the resulting ndarray of scalars.
    return scalars


def choose_color_map(
    scalars: np.ndarray,
) -> tuple[matplotlib.colors.Colormap, float, float]:
    """Chooses the color map to color a set of scalars with and finds the limits to
    apply to it.

    Scalars that keep one sign get a sequential color map, since their values run in a
    single direction. Scalars that change sign get a diverging color map with limits
    placed symmetrically about zero, so the map's midpoint marks where the scalar
    changes sign. In both cases the limits are held within a fixed number of standard
    deviations of the mean so that a few outlying Panels cannot flatten the contrast
    across the rest of the geometry.

    :param scalars: A (N,) ndarray of floats representing the scalar value at each of
        the N Panels being colored.
    :return: A tuple holding the color map, the lower limit to apply to it, and the
        upper limit to apply to it.
    """
    if np.sign(np.min(scalars)) == np.sign(np.max(scalars)):
        color_map: matplotlib.colors.Colormap = _colormaps.sequential_color_map
        c_min = max(
            float(np.mean(scalars)) - _COLOR_MAP_NUM_SIG * float(np.std(scalars)),
            float(np.min(scalars)),
        )
        c_max = min(
            float(np.mean(scalars)) + _COLOR_MAP_NUM_SIG * float(np.std(scalars)),
            float(np.max(scalars)),
        )
    else:
        color_map = _colormaps.diverging_color_map
        c_min = -_COLOR_MAP_NUM_SIG * float(np.std(scalars))
        c_max = _COLOR_MAP_NUM_SIG * float(np.std(scalars))

    return color_map, c_min, c_max


# The webp package re-exports Pillow's Image module under the name webp.Image, so the
# Image class itself is webp.Image.Image. Naming the class through webp keeps this
# module off Pillow, which Ptera Software only depends on transitively through webp.
def screenshot_image(plotter: pv.Plotter) -> webp.Image.Image:
    """Takes a screenshot of a Plotter and returns it as an Image.

    :param plotter: The Plotter to capture.
    :return: The captured frame as an Image with a transparent background.
    """
    return webp.Image.fromarray(
        np.array(
            plotter.screenshot(
                filename=None,
                transparent_background=True,
                return_img=True,
            )
        )
    )


class AnimationWriter:
    """Encodes an animation's frames into a WebP file as they are captured.

    **Contains the following methods:**

    add_frame: Hands a captured frame to the writer.

    close: Finishes the animation and writes it to its file.

    **Notes:**

    The frames are handed to a background thread that feeds them to libwebp's animation
    encoder one at a time, so the raw frames never accumulate. The queue between the
    capture loop and the thread is bounded, so a capture loop that outruns the encode
    blocks rather than piling frames up. The encode releases the global interpreter
    lock, so it overlaps with the rendering on the main thread rather than following it.

    The encoder lives on the background thread alone, from its creation through the
    assembly of the file, so no two threads ever touch it. The thread is a daemon, so a
    capture loop that raises and abandons the writer leaves nothing that could hold the
    interpreter open at exit.

    The file this writes is identical to the one the webp package's save_images produces
    from the same frames, frame rate, quality, and compression method. It gives each
    frame the same cumulative timestamp and uses the same encoder options.
    """

    def __init__(self, path: Path, frame_rate: float, quality: float) -> None:
        """The initialization method.

        :param path: The path of the WebP file to write. Its directory must already
            exist.
        :param frame_rate: The frame rate to play the animation back at, in frames per
            second. It must be positive.
        :param quality: The quality of the saved WebP, where 0.0 is the smallest file
            with the most compression artifacts and 100.0 is the largest file with the
            fewest.
        :return: None
        """
        self._path = path
        self._frame_rate = frame_rate
        self._quality = quality

        # The queue carries the captured frames to the encoding thread, and None marks
        # the end of the animation.
        self._queue: queue.Queue[webp.Image.Image | None] = queue.Queue(
            maxsize=_ANIMATION_WRITER_QUEUE_DEPTH
        )

        # An exception the encoding thread raises is held here until close re-raises it
        # on the calling thread.
        self._error: Exception | None = None

        self._thread = threading.Thread(
            target=self._encode, name="animation-writer", daemon=True
        )
        self._thread.start()

    def add_frame(self, image: webp.Image.Image) -> None:
        """Hands a captured frame to the writer.

        This blocks while the queue of frames waiting to be encoded is full. Every frame
        must have the same size as the first.

        :param image: The frame as an Image with a transparent background.
        :return: None
        """
        self._queue.put(image)

    def close(self) -> None:
        """Finishes the animation and writes it to its file.

        This blocks until every frame has been encoded and the file has been written.
        Any exception the encoding raised is re-raised here.

        :return: None
        """
        self._queue.put(None)
        self._thread.join()
        if self._error is not None:
            raise self._error

    def _encode(self) -> None:
        """Encodes the frames from the queue into the WebP file, on the background
        thread.

        An exception is held for close to re-raise. If it was raised before the queue's
        end marker was read, the queue is then drained through that marker so that
        add_frame never blocks on a queue that nothing is emptying.

        :return: None
        """
        ended = False
        try:
            config = webp.WebPConfig.new(
                lossless=False, quality=self._quality, method=WEBP_METHOD
            )
            encoder = None
            frame_size: tuple[int, int] | None = None
            num_frames = 0
            while True:
                image = self._queue.get()
                if image is None:
                    ended = True
                    break

                # The encoder is sized to the first frame, since the frames are the size
                # of the render window, which can be larger than the size that was asked
                # for on a display that scales its pixels.
                picture = webp.WebPPicture.from_pil(image)
                picture_size = (int(picture.ptr.width), int(picture.ptr.height))
                if encoder is None:
                    encoder = webp.WebPAnimEncoder.new(
                        picture_size[0],
                        picture_size[1],
                        webp.WebPAnimEncoderOptions.new(),
                    )
                    frame_size = picture_size
                elif picture_size != frame_size:
                    assert frame_size is not None
                    raise ValueError(
                        f"Every frame must have the same size as the first, which "
                        f"was {frame_size[0]} by {frame_size[1]} pixels, but frame "
                        f"{num_frames} was {picture_size[0]} by {picture_size[1]}."
                    )

                # Each frame starts at a cumulative timestamp, in whole milliseconds,
                # rather than lasting a fixed duration, which plays a fractional frame
                # rate back to within a small fraction of a percent.
                encoder.encode_frame(
                    picture, round((num_frames * 1000) / self._frame_rate), config
                )
                num_frames += 1

            if encoder is None:
                raise ValueError("An animation must have at least one frame.")

            animation = encoder.assemble(round((num_frames * 1000) / self._frame_rate))
            with open(self._path, "wb") as animation_file:
                animation_file.write(animation.buffer())
        except Exception as error:
            self._error = error
            while not ended:
                ended = self._queue.get() is None


def settle_scalar_bar_layout(plotter: pv.Plotter) -> None:
    """Runs the render pass that settles a Plotter's scalar bar layout, without
    displaying its result.

    Any actor with an opacity below one in the scene, such as the image surface plane, a
    preview ghost, or a translucent MuJoCo geom, leaves VTK's UnconstrainedFontSize
    layout with the scalar bar's labels off their bar edges on the first render after
    the scene is assembled (PyVista issue #7516). The layout only settles once the bar
    has been rendered again, so this marks the bars modified and renders. The result of
    that render is the mispositioned layout, so the buffers are held back from swapping,
    which keeps it off the screen and leaves the settled layout to the caller's next
    render. The render window is driven directly rather than through the Plotter, since
    the Plotter's render does nothing before its first show.

    :param plotter: The Plotter whose scalar bar layout to settle.
    :return: None
    """
    scalar_bar_actors = list(plotter.scalar_bars.values())
    if not scalar_bar_actors:
        return

    for scalar_bar_actor in scalar_bar_actors:
        scalar_bar_actor.Modified()

    render_window = plotter.ren_win
    assert render_window is not None
    render_window.SwapBuffersOff()
    render_window.Render()
    render_window.SwapBuffersOn()


def _plot_scalars(
    plotter: pv.Plotter,
    these_scalars: np.ndarray,
    scalar_type: str,
    min_scalar: float,
    max_scalar: float,
    color_map: matplotlib.colors.Colormap,
    c_min: float,
    c_max: float,
    panel_surfaces: pv.PolyData,
    window_scale: float,
    text_color: tuple[int, int, int] = TEXT_COLOR,
) -> list[pv.Actor]:
    """Plots a scalar bar, the surfaces of a set of Panels with particular scalars, and
    labels for the minimum and maximum scalar values.

    :param plotter: The Plotter used for visualization.
    :param these_scalars: A (N,) ndarray of floats representing the N Panels' load
        coefficients.
    :param scalar_type: Which load coefficient is represented by the scalars. Can be
        "induced drag", "side force", or "lift".
    :param min_scalar: Minimum scalar value, which is displayed as text on the Plotter.
    :param max_scalar: Maximum scalar value, which is displayed as text on the Plotter.
    :param color_map: The color map to use for scalar visualization.
    :param c_min: Lower bound for the color map scaling.
    :param c_max: Upper bound for the color map scaling.
    :param panel_surfaces: PolyData representing the Panels' surfaces.
    :param window_scale: The factor by which to scale the font sizes, as returned by
        get_window_scale.
    :param text_color: The color used for the scalar bar and label text. The default is
        TEXT_COLOR.
    :return: A list of the actors added to the plotter.
    """
    scalar_bar_args = dict(
        title=scalar_type.title() + " Coefficient",
        title_font_size=round(_BAR_TITLE_FONT_SIZE * window_scale),
        label_font_size=round(_BAR_LABEL_FONT_SIZE * window_scale),
        width=_BAR_WIDTH,
        position_x=_BAR_POSITION_X,
        position_y=_BAR_POSITION_Y,
        n_labels=_BAR_N_LABELS,
        fmt="%#.3G",
        color=text_color,
        # Suppress the render that adding the scalar bar would otherwise trigger. The
        # add_mesh call below forwards these arguments to add_scalar_bar, which renders
        # whenever the Plotter has already been shown unless render is passed
        # explicitly, and it ignores the render argument given to add_mesh. Without
        # this, each of an animation's frames is rendered half-assembled and then again
        # once whole, which reads on screen as the scalar bar and the labels flashing.
        render=False,
    )
    panel_actor = plotter.add_mesh(
        panel_surfaces,
        show_edges=True,
        line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
        cmap=color_map,
        clim=[c_min, c_max],
        scalars=these_scalars,
        smooth_shading=False,
        scalar_bar_args=scalar_bar_args,  # type: ignore[arg-type]
        lighting=False,
        render=False,
    )

    max_label = plotter.add_text(
        text=f"Max: {max_scalar:#.3G}",
        position=_TEXT_MAX_POSITION,
        font_size=round(TEXT_FONT_SIZE * window_scale),
        viewport=True,
        color=text_color,
        render=False,
    )
    min_label = plotter.add_text(
        text=f"Min: {min_scalar:#.3G}",
        position=_TEXT_MIN_POSITION,
        font_size=round(TEXT_FONT_SIZE * window_scale),
        viewport=True,
        color=text_color,
        render=False,
    )
    for label in (max_label, min_label):
        label.prop.justification_horizontal = "right"

    return [panel_actor, max_label, min_label]


class ScalarColoring(NamedTuple):
    """Everything needed to color one frame's Panels by a scalar, and to color their
    reflected copies to match.

    Only the scalars themselves change from frame to frame during an animation. The
    remaining fields are constant across the whole animation, because the color map and
    its limits are chosen once from the scalars of every time step so that a given color
    means the same value in every frame. muted_color_map is carried here, rather than
    derived where it is used, because building it walks a 256 entry table and it would
    otherwise be rebuilt on every frame.
    """

    scalars: np.ndarray
    scalar_type: str
    min_scalar: float
    max_scalar: float
    color_map: matplotlib.colors.Colormap
    muted_color_map: matplotlib.colors.Colormap
    c_min: float
    c_max: float


def add_frame_geometry(
    plotter: pv.Plotter,
    panel_surfaces: pv.PolyData,
    wake_surfaces: pv.PolyData | None,
    coloring: ScalarColoring | None,
    T_reflect: np.ndarray | None,
    window_scale: float,
) -> list[pv.Actor]:
    """Adds one frame's geometry to a Plotter, which is the wake, the Panels, and, when
    an image surface is defined, a muted reflected copy of both.

    A single drawing and every frame of an animation are the same scene, so both are
    assembled here. The meshes arrive already mapped into whichever axes they are being
    rendered in, so this function does not need to know whether it is drawing a free
    flight body in Earth axes or a body-fixed one in geometry axes.

    The image surface plane itself is not added here. A drawing sizes its plane from the
    bounds of the geometry added above, so it cannot be added until afterward, and it
    must in any case be added last so it renders over the geometry it is blended with.

    :param plotter: The Plotter to add the geometry to.
    :param panel_surfaces: The PolyData representation of the frame's Panel surfaces.
    :param wake_surfaces: The PolyData representation of the frame's wake ring vortex
        surfaces, or None to omit the wake.
    :param coloring: The scalar coloring to apply to the Panels, or None to color them
        uniformly.
    :param T_reflect: A (4,4) ndarray of floats representing the active transformation
        that reflects geometry across the image surface, or None when no image surface
        is defined, in which case no reflected geometry is added.
    :param window_scale: The factor by which to scale the scalar bar and label font
        sizes, as returned by get_window_scale.
    :return: A list of the actors added to the plotter.
    """
    actors: list[pv.Actor] = []
    # Add the wake ring vortex surfaces if they are being shown.
    if wake_surfaces is not None:
        actors.append(
            plotter.add_mesh(
                wake_surfaces,
                show_edges=True,
                line_width=_WAKE_VORTEX_EDGE_LINE_WIDTH * window_scale,
                smooth_shading=False,
                color=_WAKE_VORTEX_COLOR,
                lighting=False,
                render=False,
            )
        )

    # Plot the Panels either with scalar coloring or with a uniform color.
    if coloring is not None:
        actors.extend(
            _plot_scalars(
                plotter,
                coloring.scalars,
                coloring.scalar_type,
                coloring.min_scalar,
                coloring.max_scalar,
                coloring.color_map,
                coloring.c_min,
                coloring.c_max,
                panel_surfaces,
                window_scale,
                text_color=TEXT_COLOR_SURFACE if T_reflect is not None else TEXT_COLOR,
            )
        )
    else:
        actors.append(
            plotter.add_mesh(
                panel_surfaces,
                show_edges=True,
                line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
                color=_PANEL_COLOR,
                smooth_shading=False,
                lighting=False,
                render=False,
            )
        )

    if T_reflect is None:
        return actors

    # An image surface is defined, so add the reflected geometry. It is muted toward
    # gray so that it reads as a reflection rather than as more geometry.
    mute = IMAGE_REFLECTION_MUTE_FACTOR
    muted_edge_color = mute_color("black", mute)

    # Add reflected Panel surfaces with muted coloring.
    reflected_panel_surfaces = _reflect_mesh(panel_surfaces, T_reflect)
    if coloring is not None:
        actors.append(
            plotter.add_mesh(
                reflected_panel_surfaces,
                show_edges=True,
                line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
                edge_color=muted_edge_color,
                cmap=coloring.muted_color_map,
                clim=[coloring.c_min, coloring.c_max],
                scalars=coloring.scalars,
                smooth_shading=False,
                show_scalar_bar=False,
                lighting=False,
                render=False,
            )
        )

    else:
        actors.append(
            plotter.add_mesh(
                reflected_panel_surfaces,
                show_edges=True,
                line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
                edge_color=muted_edge_color,
                color=mute_color(_PANEL_COLOR, mute),
                smooth_shading=False,
                lighting=False,
                render=False,
            )
        )

    # Add reflected wake ring vortex surfaces if they are being shown.
    if wake_surfaces is not None:
        actors.append(
            plotter.add_mesh(
                _reflect_mesh(wake_surfaces, T_reflect),
                show_edges=True,
                line_width=_WAKE_VORTEX_EDGE_LINE_WIDTH * window_scale,
                edge_color=muted_edge_color,
                smooth_shading=False,
                color=mute_color(_WAKE_VORTEX_COLOR, mute),
                lighting=False,
                render=False,
            )
        )
    return actors


def add_mujoco_geometry(
    plotter: pv.Plotter,
    worldbody_geoms: list[_mujoco_model.RenderGeom],
    body_geoms: list[_mujoco_model.RenderGeom],
    T_pas_BP1_CgP1_to_E_Eo: np.ndarray,
    T_reflect: np.ndarray | None,
) -> list[pv.Actor]:
    """Adds one frame's MuJoCo geometry to a Plotter, alongside muted reflected copies
    when an image surface is defined.

    Worldbody geoms arrive with their meshes already in Earth axes, relative to the
    Earth origin, where MuJoCo holds them static, so they are added as they are. Body
    geoms arrive rigid in the first Airplane's body axes, relative to the first
    Airplane's CG, so each frame maps them through that time step's transformation to
    fly them along with the Panel meshes. The reflected copies are muted toward gray,
    matching the Panel treatment, so they read as reflections rather than as more
    geometry.

    The geom actors are the only lit actors in any scene. When there is geometry to add,
    a headlight comes with it, and every other actor disables its lighting, which
    renders identically with or without a light present, so the shading reaches the
    geoms alone and every scene without MuJoCo geometry stays free of lights entirely.

    :param plotter: The Plotter to add the geometry to.
    :param worldbody_geoms: The RenderGeoms attached to the worldbody.
    :param body_geoms: The RenderGeoms attached to the body.
    :param T_pas_BP1_CgP1_to_E_Eo: A (4,4) ndarray of floats representing the passive
        transformation from the first Airplane's body axes, relative to the first
        Airplane's CG, to Earth axes, relative to the Earth origin, at the frame's time
        step.
    :param T_reflect: A (4,4) ndarray of floats representing the active transformation
        (in Earth axes) that reflects geometry across the image surface, or None when no
        image surface is defined, in which case no reflected copies are added.
    :return: A list of the actors added to the plotter.
    """
    actors: list[pv.Actor] = []
    posed_meshes = [
        (render_geom.mesh, render_geom.rgba) for render_geom in worldbody_geoms
    ] + [
        (transform_mesh(render_geom.mesh, T_pas_BP1_CgP1_to_E_Eo), render_geom.rgba)
        for render_geom in body_geoms
    ]
    if not posed_meshes:
        return actors

    # Add the headlight that shades the geom actors, which follows the camera so the
    # geoms stay legible from any orientation the user chooses. Every non-geom actor is
    # added with lighting disabled, which renders identically with or without a light in
    # the scene, so the headlight changes nothing else. An animation clears the
    # Plotter's lights along with its actors between frames, so the light is added here,
    # once per assembled frame.
    if not plotter.renderer.lights:
        plotter.add_light(pv.Light(light_type="headlight"))

    for posed_mesh, rgba in posed_meshes:
        color = (float(rgba[0]), float(rgba[1]), float(rgba[2]))
        opacity = float(rgba[3])

        # The mesh carries its shading normals, computed with the sharp edges split when
        # get_mujoco_render_geometry extracted it, so PyVista passes them through here
        # rather than recomputing them for every frame.
        actor = plotter.add_mesh(
            posed_mesh,
            color=color,
            opacity=opacity,
            smooth_shading=True,
            ambient=_MUJOCO_GEOMETRY_AMBIENT,
            diffuse=_MUJOCO_GEOMETRY_DIFFUSE,
            render=False,
        )
        actors.append(actor)
        if T_reflect is not None:
            # A reflection inverts the faces' winding. Flip the reflected copy's faces
            # so they wind outward again, which the shading needs even though the
            # reflected normals it carries already point outward.
            actor = plotter.add_mesh(
                _reflect_mesh(posed_mesh, T_reflect).flip_faces(),
                color=mute_color(color, IMAGE_REFLECTION_MUTE_FACTOR),
                opacity=opacity,
                smooth_shading=True,
                ambient=_MUJOCO_GEOMETRY_AMBIENT,
                diffuse=_MUJOCO_GEOMETRY_DIFFUSE,
                render=False,
            )
            actors.append(actor)
    return actors
