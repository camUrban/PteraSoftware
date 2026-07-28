"""Contains the styling, geometry building, and scene assembly used to render the
visualizations with PyVista."""

from __future__ import annotations

from typing import NamedTuple

import matplotlib.colors
import numpy as np
import pyvista as pv
import webp

from . import _colormaps, _transformations, geometry
from . import operating_point as operating_point_mod
from . import unsteady_ring_vortex_lattice_method

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
_TEXT_MAX_POSITION = (0.85, 0.075)
_TEXT_MIN_POSITION = (0.85, 0.050)
TEXT_FONT_SIZE = 11

# Define the render window size that every font size and line width in the
# visualizations is tuned against.
REFERENCE_WINDOW_SIZE = (1024, 768)


def get_window_scale(window_width: int, window_height: int) -> float:
    """Returns the factor by which to scale the font sizes and line widths that are
    tuned for REFERENCE_WINDOW_SIZE.

    VTK sizes text and lines in pixels, so without scaling they keep the same pixel
    count at every window size. That crowds a small window, where the scalar bar's
    labels overlap the bar itself, and it leaves them nearly unreadable in a large one.

    The smaller of the two ratios is used rather than the height ratio alone. The scalar
    labels are anchored near the right edge and grow rightward, so a wide but short
    window runs out of horizontal room before it runs out of vertical room.

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


def get_panel_surfaces(
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
                # Arrange this Panel's vertices and faces into ndarrays in the proper
                # form to represent PolyData surfaces.
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

                # Add this Panel's vertices and faces to the ndarray of all vertices and
                # faces.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # Update the number of Panels.
                panel_num += 1

    # Return the Panels' surfaces.
    return pv.PolyData(panel_vertices, panel_faces)


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


def transform_mesh(
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


def get_scalars(
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
    :param color_map: The color map to use for scalar visualization.
    :param c_min: Lower bound for the color map scaling.
    :param c_max: Upper bound for the color map scaling.
    :param panel_surfaces: PolyData representing the Panels' surfaces.
    :param window_scale: The factor by which to scale the font sizes, as returned by
        get_window_scale.
    :param text_color: The color used for the scalar bar and label text. The default is
        TEXT_COLOR.
    :return: None
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
    )
    plotter.add_mesh(
        panel_surfaces,
        show_edges=True,
        line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
        cmap=color_map,
        clim=[c_min, c_max],
        scalars=these_scalars,
        smooth_shading=False,
        scalar_bar_args=scalar_bar_args,  # type: ignore[arg-type]
    )

    plotter.add_text(
        text=f"Max: {max_scalar:#.3G}",
        position=_TEXT_MAX_POSITION,
        font_size=round(TEXT_FONT_SIZE * window_scale),
        viewport=True,
        color=text_color,
    )
    plotter.add_text(
        text=f"Min: {min_scalar:#.3G}",
        position=_TEXT_MIN_POSITION,
        font_size=round(TEXT_FONT_SIZE * window_scale),
        viewport=True,
        color=text_color,
    )


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
) -> None:
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
    :return: None
    """
    # Add the wake ring vortex surfaces if they are being shown.
    if wake_surfaces is not None:
        plotter.add_mesh(
            wake_surfaces,
            show_edges=True,
            line_width=_WAKE_VORTEX_EDGE_LINE_WIDTH * window_scale,
            smooth_shading=False,
            color=_WAKE_VORTEX_COLOR,
        )

    # Plot the Panels either with scalar coloring or with a uniform color.
    if coloring is not None:
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
    else:
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
            color=_PANEL_COLOR,
            smooth_shading=False,
        )

    if T_reflect is None:
        return

    # An image surface is defined, so add the reflected geometry. It is muted toward
    # gray so that it reads as a reflection rather than as more geometry.
    mute = IMAGE_REFLECTION_MUTE_FACTOR
    muted_edge_color = mute_color("black", mute)

    # Add reflected Panel surfaces with muted coloring.
    reflected_panel_surfaces = _reflect_mesh(panel_surfaces, T_reflect)
    if coloring is not None:
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
        )
    else:
        plotter.add_mesh(
            reflected_panel_surfaces,
            show_edges=True,
            line_width=_PANEL_EDGE_LINE_WIDTH * window_scale,
            edge_color=muted_edge_color,
            color=mute_color(_PANEL_COLOR, mute),
            smooth_shading=False,
        )

    # Add reflected wake ring vortex surfaces if they are being shown.
    if wake_surfaces is not None:
        plotter.add_mesh(
            _reflect_mesh(wake_surfaces, T_reflect),
            show_edges=True,
            line_width=_WAKE_VORTEX_EDGE_LINE_WIDTH * window_scale,
            edge_color=muted_edge_color,
            smooth_shading=False,
            color=mute_color(_WAKE_VORTEX_COLOR, mute),
        )
