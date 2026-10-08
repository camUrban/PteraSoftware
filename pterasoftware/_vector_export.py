"""Contains the classes and functions that export the visualizations' scenes as svg and
pdf files.

A scene is recorded while it is assembled, in diagram axes, relative to the diagram
origin, and then exported through a parallel camera. The export projects every point to
display coordinates, which hold each point's x and y position in pixels, measured from
the window's bottom left corner, and its depth in pixels along the camera's view
direction, so that a larger depth is farther from the camera. Display coordinates are
the diagram axes rotated onto the camera's basis directions, scaled, and shifted, so
they keep every plane flat and every intersection where it was.

Visibility is resolved exactly in display coordinates. The opaque fills are split and
ordered by a binary space partitioning (BSP) tree, which paints them from back to front.
The strokes are clipped against every nearer triangle, which removes their hidden
stretches, and are painted after the fills. Matplotlib then writes the result, which
makes the svg and pdf writers the same ones the results plots use.
"""

from __future__ import annotations

import io
import re
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

import matplotlib
import matplotlib.axes
import matplotlib.collections
import matplotlib.colors
import matplotlib.figure
import matplotlib.font_manager
import matplotlib.patches
import matplotlib.path
import matplotlib.typing
import numpy as np
import scipy.spatial

from . import _fonts

# Define the tolerance, in pixels, within which a vertex counts as lying on a splitting
# plane, and by which a stroke wins a depth tie against a triangle, so that an edge
# shared by two faces is never hidden by the face that doesn't own it.
_PLANE_TOLERANCE = 1.0e-3

# Define the area, in square pixels, below which a fragment split off by the BSP is
# dropped rather than painted. Such a fragment is a sliver that no pixel can show.
_SLIVER_AREA = 1.0e-6

# Define the tolerance below which a polygon's normal or a triangle's footprint counts
# as degenerate.
_DEGENERATE_TOLERANCE = 1.0e-12

# Define the number of candidate splitting polygons the BSP auditions at each node, and
# the seed of the generator that samples them. The candidate that cuts the fewest other
# polygons is chosen, which splits far fewer polygons on a sheet-like mesh than always
# taking the first polygon. Seeding the generator keeps every export of the same scene
# identical.
_SPLITTER_CANDIDATES = 8
_SPLITTER_SEED = 7

# Define the side, in pixels, of the cells of the grid on screen that the triangles are
# registered in when clipping strokes, so that each stroke is only tested against the
# triangles near it.
_GRID_CELL_SIZE = 32.0

# Define the shortest fraction of a stroke that is still drawn once its hidden stretches
# are removed.
_MIN_VISIBLE_FRACTION = 1.0e-9

# Define the resolution the figures are laid out at. Each pixel of the window becomes
# one unit of the figure's axes, and is this many points wide.
_FIGURE_DPI = 100.0
POINTS_PER_PIXEL = 72.0 / _FIGURE_DPI

# Define the width, in pixels, of the stroke each fill is outlined with in its own
# color. Abutting fills are anti-aliased separately by the program that displays the
# file, which leaves faint seams between them unless each covers its edges this way.
_FILL_SEAM_LINE_WIDTH = 0.4

# Define the number of decimal places the coordinates of an svg's paths are rounded to.
# The coordinates are in points, so this rounds them to the nearest hundredth of a
# point, far finer than any display or print can show.
_SVG_DECIMAL_PLACES = 2


class VectorCamera(NamedTuple):
    """The parallel camera a scene is exported through, along with the size of the
    window it renders into.

    :param position_D_Do: A (3,) ndarray of floats holding the camera's position (in
        diagram axes, relative to the diagram origin). The units are in meters.
    :param focalPoint_D_Do: A (3,) ndarray of floats holding the camera's focal point
        (in diagram axes, relative to the diagram origin), which is shown at the center
        of the window. The units are in meters.
    :param viewUp_D: A (3,) ndarray of floats holding the camera's view up direction (in
        diagram axes). It needn't be perpendicular to the view direction or of unit
        length, since only its part perpendicular to the view direction is used.
    :param parallel_scale: Half of the window's height, in meters. It must be positive.
    :param window_width: The window's width, in pixels. It must be positive.
    :param window_height: The window's height, in pixels. It must be positive.
    """

    position_D_Do: np.ndarray
    focalPoint_D_Do: np.ndarray
    viewUp_D: np.ndarray
    parallel_scale: float
    window_width: float
    window_height: float

    def fit_to_page(self, page_width_in: float, page_height_in: float) -> VectorCamera:
        """Returns a camera that shows this camera's view scaled to fit a page.

        The returned camera's window is the page, in pixels at the resolution the
        figures are laid out at. The view is scaled by the largest factor that keeps all
        of it on the page, and stays centered on the focal point, so on a page whose
        aspect ratio differs from the window's, more of the scene shows along one
        direction.

        :param page_width_in: The page's width, in inches. It must be positive.
        :param page_height_in: The page's height, in inches. It must be positive.
        :return: The VectorCamera whose window is the page.
        """
        page_width = page_width_in * _FIGURE_DPI
        page_height = page_height_in * _FIGURE_DPI
        scale = min(page_width / self.window_width, page_height / self.window_height)
        return self._replace(
            parallel_scale=self.parallel_scale
            * page_height
            / (scale * self.window_height),
            window_width=page_width,
            window_height=page_height,
        )

    def to_display(self, stackPoints_D_Do: np.ndarray) -> np.ndarray:
        """Projects points to display coordinates.

        :param stackPoints_D_Do: A (N,3) ndarray of floats holding the points' positions
            (in diagram axes, relative to the diagram origin). The units are in meters.
        :return: A (N,3) ndarray of floats holding each point's x and y position, in
            pixels from the window's bottom left corner, and its depth, in pixels along
            the view direction from the focal point.
        """
        viewDirection_D = self.focalPoint_D_Do - self.position_D_Do
        viewDirection_D = viewDirection_D / np.linalg.norm(viewDirection_D)
        rightDirection_D = np.cross(viewDirection_D, self.viewUp_D)
        rightDirection_D = rightDirection_D / np.linalg.norm(rightDirection_D)
        upDirection_D = np.cross(rightDirection_D, viewDirection_D)

        pixels_per_meter = self.window_height / (2.0 * self.parallel_scale)
        stackOffsets_D = stackPoints_D_Do - self.focalPoint_D_Do
        return np.column_stack(
            [
                0.5 * self.window_width
                + pixels_per_meter * (stackOffsets_D @ rightDirection_D),
                0.5 * self.window_height
                + pixels_per_meter * (stackOffsets_D @ upDirection_D),
                pixels_per_meter * (stackOffsets_D @ viewDirection_D),
            ]
        )


class VectorText(NamedTuple):
    """A run of text placed in the window, drawn over the whole scene.

    :param x: The x position of the text's anchor, in pixels from the window's left
        edge.
    :param y: The y position of the text's anchor, in pixels from the window's bottom
        edge.
    :param text: The text. Any part between a pair of dollar signs is written as math
        with Matplotlib's mathtext.
    :param font_properties: The FontProperties that select the text's font. They must
        set both the font's file, which selects the font, and its family name, since an
        svg names its font by family. Their size is ignored in favor of font_size.
    :param font_size: The text's size, in pixels.
    :param color: The text's color, as any color Matplotlib accepts.
    :param horizontal_alignment: Which part of the text sits at the anchor horizontally.
        It must be "left", "center", or "right".
    :param vertical_alignment: Which part of the text sits at the anchor vertically. It
        must be "bottom", "baseline", "center", "center_baseline", or "top".
    :param math_font_family: The font set that the text's math is written in, as any
        font set Matplotlib's mathtext accepts, such as "stix". None leaves it at
        Matplotlib's default. The default is None.
    :param background_box: The box filled behind the text, as a tuple of its left,
        right, bottom, and top edges, in pixels from the window's bottom left corner, or
        None to fill nothing behind the text. The default is None.
    :param background_color: The color of the box filled behind the text, as any color
        Matplotlib accepts. It has no effect if background_box is None. The default is
        "white".
    """

    x: float
    y: float
    text: str
    font_properties: matplotlib.font_manager.FontProperties
    font_size: float
    color: matplotlib.typing.ColorType
    horizontal_alignment: str
    vertical_alignment: str
    math_font_family: str | None = None
    background_box: tuple[float, float, float, float] | None = None
    background_color: matplotlib.typing.ColorType = "white"


class _BspNode:
    """A node of a BSP tree, holding its splitting plane, the polygons that lie in that
    plane, and the subtrees on either side of it."""

    __slots__ = ("plane_normal_display", "plane_offset", "coplanar", "front", "back")

    def __init__(self) -> None:
        """The initialization method.

        :return: None
        """
        self.plane_normal_display = np.array([0.0, 0.0, 1.0], dtype=float)
        self.plane_offset = 0.0
        self.coplanar: list[tuple[np.ndarray, int]] = []
        self.front: _BspNode | None = None
        self.back: _BspNode | None = None


def _get_newell_normal(vertices_display: np.ndarray) -> np.ndarray:
    """Returns a planar polygon's normal, scaled by twice its area, by Newell's method.

    :param vertices_display: A (K,3) ndarray of floats holding the polygon's vertices,
        in order around it (in display coordinates).
    :return: A (3,) ndarray of floats holding the polygon's normal (in display
        coordinates), whose length is twice the polygon's area.
    """
    rolled_display = np.roll(vertices_display, -1, axis=0)
    return np.array(
        [
            np.sum(
                (vertices_display[:, 1] - rolled_display[:, 1])
                * (vertices_display[:, 2] + rolled_display[:, 2])
            ),
            np.sum(
                (vertices_display[:, 2] - rolled_display[:, 2])
                * (vertices_display[:, 0] + rolled_display[:, 0])
            ),
            np.sum(
                (vertices_display[:, 0] - rolled_display[:, 0])
                * (vertices_display[:, 1] + rolled_display[:, 1])
            ),
        ],
        dtype=float,
    )


def _get_plane(vertices_display: np.ndarray) -> tuple[np.ndarray, float]:
    """Returns the plane a polygon lies in.

    A degenerate polygon has no plane of its own, so it is given the plane of constant
    depth through its first vertex.

    :param vertices_display: A (K,3) ndarray of floats holding the polygon's vertices,
        in order around it (in display coordinates).
    :return: A tuple of the plane's unit normal, as a (3,) ndarray of floats (in display
        coordinates), and its offset, which is the dot product of that normal with any
        point in the plane, in pixels.
    """
    normal_display = _get_newell_normal(vertices_display)
    length = float(np.linalg.norm(normal_display))
    if length < _DEGENERATE_TOLERANCE:
        return np.array([0.0, 0.0, 1.0], dtype=float), float(vertices_display[0, 2])
    normal_display = normal_display / length
    return normal_display, float(normal_display @ vertices_display[0])


def _split_spanning_polygon(
    vertices_display: np.ndarray, distances: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Splits a convex polygon that spans a plane into its parts on either side.

    **Citation:**

    Adapted from: CSG.Plane.prototype.splitPolygon in csg.js

    Author: Evan Wallace

    :param vertices_display: A (K,3) ndarray of floats holding the polygon's vertices,
        in order around it (in display coordinates).
    :param distances: A (K,) ndarray of floats holding each vertex's signed distance
        from the plane, in pixels, which is positive on the plane's front side.
    :return: A tuple of two ndarrays of floats, holding the vertices of the part in
        front of the plane and of the part behind it, each in order around it (in
        display coordinates).
    """
    sides = np.zeros(distances.shape[0], dtype=int)
    sides[distances > _PLANE_TOLERANCE] = 1
    sides[distances < -_PLANE_TOLERANCE] = -1

    front_vertices: list[np.ndarray] = []
    back_vertices: list[np.ndarray] = []
    num_vertices = vertices_display.shape[0]
    for vertex_id in range(num_vertices):
        next_vertex_id = (vertex_id + 1) % num_vertices
        side = sides[vertex_id]
        next_side = sides[next_vertex_id]
        if side >= 0:
            front_vertices.append(vertices_display[vertex_id])
        if side <= 0:
            back_vertices.append(vertices_display[vertex_id])

        # Where the edge to the next vertex crosses the plane, both parts gain the
        # crossing point.
        if side * next_side < 0:
            fraction = distances[vertex_id] / (
                distances[vertex_id] - distances[next_vertex_id]
            )
            crossing_display = vertices_display[vertex_id] + fraction * (
                vertices_display[next_vertex_id] - vertices_display[vertex_id]
            )
            front_vertices.append(crossing_display)
            back_vertices.append(crossing_display)
    return np.array(front_vertices, dtype=float), np.array(back_vertices, dtype=float)


def _stack_polygon_vertices(
    polygons: Sequence[tuple[np.ndarray, int]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stacks the vertices of a sequence of polygons, so that they can all be measured
    against a plane at once.

    :param polygons: The polygons, each a tuple of a (K,3) ndarray of floats holding its
        vertices (in display coordinates) and the int ID of the fill it was cut from.
    :return: A tuple of three ndarrays. The first is a (V,3) ndarray of floats holding
        every polygon's vertices in turn (in display coordinates). The second and third
        are (P,) ndarrays of ints holding each polygon's number of vertices and the row
        of its first vertex.
    """
    counts = np.array([vertices.shape[0] for vertices, _ in polygons], dtype=int)
    starts = np.concatenate([np.zeros(1, dtype=int), np.cumsum(counts)[:-1]])
    return np.concatenate([vertices for vertices, _ in polygons]), counts, starts


def _choose_splitter(
    polygons: list[tuple[np.ndarray, int]], rng: np.random.Generator
) -> tuple[np.ndarray, float]:
    """Chooses a node's splitting polygon, moves it to the front of the node's polygons,
    and returns its plane.

    :param polygons: The node's polygons, each a tuple of a (K,3) ndarray of floats
        holding its vertices (in display coordinates) and the int ID of the fill it was
        cut from. It is reordered in place so that the splitter comes first.
    :param rng: The Generator that samples the candidate splitters.
    :return: A tuple of the splitting plane's unit normal, as a (3,) ndarray of floats
        (in display coordinates), and its offset, in pixels.
    """
    if len(polygons) <= 2:
        return _get_plane(polygons[0][0])

    candidate_ids = rng.choice(
        len(polygons), size=min(_SPLITTER_CANDIDATES, len(polygons)), replace=False
    )
    planes = [_get_plane(polygons[candidate_id][0]) for candidate_id in candidate_ids]
    stackNormals_display = np.array([normal for normal, _ in planes], dtype=float)
    offsets = np.array([offset for _, offset in planes], dtype=float)

    # Count the polygons each candidate's plane cuts, which are those with vertices on
    # both of its sides.
    vertices_display, _, starts = _stack_polygon_vertices(polygons)
    distances = vertices_display @ stackNormals_display.T - offsets
    max_distances = np.maximum.reduceat(distances, starts, axis=0)
    min_distances = np.minimum.reduceat(distances, starts, axis=0)
    num_cut = np.count_nonzero(
        (max_distances > _PLANE_TOLERANCE) & (min_distances < -_PLANE_TOLERANCE),
        axis=0,
    )

    best_id = int(np.argmin(num_cut))
    splitter_id = int(candidate_ids[best_id])
    polygons[0], polygons[splitter_id] = polygons[splitter_id], polygons[0]
    return planes[best_id]


def _split_node(
    node: _BspNode, polygons: list[tuple[np.ndarray, int]], rng: np.random.Generator
) -> tuple[list[tuple[np.ndarray, int]], list[tuple[np.ndarray, int]]]:
    """Chooses a node's splitting plane and sorts its polygons against it.

    The polygons in the plane are kept by the node. Those that span it are split, and
    any part that is a sliver is dropped.

    :param node: The node, whose plane and coplanar polygons are filled in.
    :param polygons: The node's polygons, each a tuple of a (K,3) ndarray of floats
        holding its vertices (in display coordinates) and the int ID of the fill it was
        cut from.
    :param rng: The Generator that samples the candidate splitters.
    :return: A tuple of two lists of polygons, which are those in front of the plane and
        those behind it.
    """
    node.plane_normal_display, node.plane_offset = _choose_splitter(polygons, rng)
    node.coplanar.append(polygons[0])

    vertices_display, counts, starts = _stack_polygon_vertices(polygons)
    distances = vertices_display @ node.plane_normal_display - node.plane_offset
    has_front = np.maximum.reduceat(distances, starts) > _PLANE_TOLERANCE
    has_back = np.minimum.reduceat(distances, starts) < -_PLANE_TOLERANCE

    front: list[tuple[np.ndarray, int]] = []
    back: list[tuple[np.ndarray, int]] = []
    for polygon_id in range(1, len(polygons)):
        polygon = polygons[polygon_id]
        if has_front[polygon_id] and has_back[polygon_id]:
            start = starts[polygon_id]
            front_display, back_display = _split_spanning_polygon(
                polygon[0], distances[start : start + counts[polygon_id]]
            )
            for part_display, destination in (
                (front_display, front),
                (back_display, back),
            ):
                if 0.5 * np.linalg.norm(_get_newell_normal(part_display)) >= (
                    _SLIVER_AREA
                ):
                    destination.append((part_display, polygon[1]))
        elif has_front[polygon_id]:
            front.append(polygon)
        elif has_back[polygon_id]:
            back.append(polygon)
        else:
            node.coplanar.append(polygon)
    return front, back


def get_paint_order(
    polygons: Sequence[tuple[np.ndarray, int]],
) -> list[tuple[np.ndarray, int]]:
    """Splits convex polygons where they cross, and orders the pieces to be painted from
    back to front.

    The polygons are partitioned by a BSP tree and then read back in the order a viewer
    at infinite negative depth sees them, so that every piece is painted after every
    piece behind it. The tree is built and read with explicit stacks rather than
    recursion, so a deep tree cannot exceed Python's recursion limit.

    **Citation:**

    Adapted from: CSG.Node.prototype.build in csg.js

    Author: Evan Wallace

    :param polygons: The convex polygons, each a tuple of a (K,3) ndarray of floats
        holding its vertices, in order around it (in display coordinates), and an int ID
        that its pieces carry.
    :return: The pieces, in the order to paint them, each a tuple of a (K,3) ndarray of
        floats holding its vertices, in order around it (in display coordinates), and
        the ID of the polygon it was cut from.
    """
    if not polygons:
        return []

    rng = np.random.default_rng(_SPLITTER_SEED)
    root = _BspNode()
    build_stack: list[tuple[_BspNode, list[tuple[np.ndarray, int]]]] = [
        (root, list(polygons))
    ]
    while build_stack:
        node, node_polygons = build_stack.pop()
        front, back = _split_node(node, node_polygons, rng)
        if front:
            node.front = _BspNode()
            build_stack.append((node.front, front))
        if back:
            node.back = _BspNode()
            build_stack.append((node.back, back))

    # The viewer is on a plane's front side exactly when the plane's normal points back
    # toward it, which is toward negative depth. Each node paints the subtree on the far
    # side of its plane, then its own polygons, then the subtree on the near side. The
    # stack is last in, first out, so they are pushed in the reverse order.
    order: list[tuple[np.ndarray, int]] = []
    traverse_stack: list[_BspNode | list[tuple[np.ndarray, int]] | None] = [root]
    while traverse_stack:
        item = traverse_stack.pop()
        if item is None:
            continue
        if isinstance(item, list):
            order.extend(item)
            continue
        if item.plane_normal_display[2] < 0.0:
            near, far = item.front, item.back
        else:
            near, far = item.back, item.front
        traverse_stack.extend([near, item.coplanar, far])
    return order


class _Occluders(NamedTuple):
    """The quantities that clipping strokes against a set of triangles needs, computed
    once for all of them.

    :param mins_display: A (T,3) ndarray of floats holding the smallest x, y, and depth
        of each triangle's vertices (in display coordinates).
    :param maxs_display: A (T,3) ndarray of floats holding the largest x, y, and depth
        of each triangle's vertices (in display coordinates).
    :param wound_display: A (T,3,3) ndarray of floats holding each triangle's vertices
        (in display coordinates), ordered counterclockwise on screen.
    :param normals_display: A (T,3) ndarray of floats holding each triangle's normal (in
        display coordinates), which needn't be of unit length.
    :param offsets: A (T,) ndarray of floats holding the dot product of each triangle's
        normal with its vertices.
    :param usable: A (T,) ndarray of bools that is False for each triangle that can't
        hide anything, which is one whose footprint on screen has no area.
    """

    mins_display: np.ndarray
    maxs_display: np.ndarray
    wound_display: np.ndarray
    normals_display: np.ndarray
    offsets: np.ndarray
    usable: np.ndarray


def _prepare_occluders(triangles_display: np.ndarray) -> _Occluders:
    """Computes the quantities that clipping strokes against a set of triangles needs.

    :param triangles_display: A (T,3,3) ndarray of floats holding each triangle's
        vertices (in display coordinates).
    :return: The _Occluders holding those quantities.
    """
    firstEdges_display = triangles_display[:, 1] - triangles_display[:, 0]
    secondEdges_display = triangles_display[:, 2] - triangles_display[:, 0]
    doubled_areas = (
        firstEdges_display[:, 0] * secondEdges_display[:, 1]
        - firstEdges_display[:, 1] * secondEdges_display[:, 0]
    )
    normals_display = np.cross(firstEdges_display, secondEdges_display)
    return _Occluders(
        mins_display=triangles_display.min(axis=1),
        maxs_display=triangles_display.max(axis=1),
        wound_display=np.where(
            (doubled_areas > 0.0)[:, np.newaxis, np.newaxis],
            triangles_display,
            triangles_display[:, ::-1, :],
        ),
        normals_display=normals_display,
        offsets=np.einsum("ij,ij->i", normals_display, triangles_display[:, 0]),
        usable=np.abs(doubled_areas) >= _DEGENERATE_TOLERANCE,
    )


def _clip_to_footprint(
    start_display: np.ndarray, direction_display: np.ndarray, wound_display: np.ndarray
) -> tuple[float, float] | None:
    """Finds the part of a stroke that lies within a triangle's footprint on screen.

    The stroke is clipped against the half plane inside each of the triangle's edges in
    turn.

    :param start_display: A (3,) ndarray of floats holding the stroke's start (in
        display coordinates).
    :param direction_display: A (3,) ndarray of floats holding the vector from the
        stroke's start to its end (in display coordinates).
    :param wound_display: A (3,3) ndarray of floats holding the triangle's vertices (in
        display coordinates), ordered counterclockwise on screen.
    :return: A tuple of the fractions along the stroke where the part within the
        footprint starts and ends, or None if the stroke misses the footprint.
    """
    low = 0.0
    high = 1.0
    for vertex_id in range(3):
        corner = wound_display[vertex_id, :2]
        along = wound_display[(vertex_id + 1) % 3, :2] - corner

        # The stroke is inside this edge's half plane where the cross product of the
        # edge with the vector from its corner to the stroke is positive, and that cross
        # product is linear in the fraction along the stroke.
        constant = along[0] * (start_display[1] - corner[1]) - along[1] * (
            start_display[0] - corner[0]
        )
        slope = along[0] * direction_display[1] - along[1] * direction_display[0]
        if abs(slope) < _DEGENERATE_TOLERANCE:
            if constant < 0.0:
                return None
            continue
        crossing = -constant / slope
        if slope > 0.0:
            low = max(low, crossing)
        else:
            high = min(high, crossing)
        if low >= high:
            return None
    return low, high


def _get_occluded_range(
    start_display: np.ndarray,
    direction_display: np.ndarray,
    normal_display: np.ndarray,
    offset: float,
    low: float,
    high: float,
) -> tuple[float, float] | None:
    """Finds the part of a stroke, within a triangle's footprint, that the triangle's
    plane is nearer than.

    Both the stroke's depth and the plane's depth beneath it are linear in the fraction
    along the stroke, so their difference crosses zero at most once. The stroke wins a
    tie within _PLANE_TOLERANCE.

    :param start_display: A (3,) ndarray of floats holding the stroke's start (in
        display coordinates).
    :param direction_display: A (3,) ndarray of floats holding the vector from the
        stroke's start to its end (in display coordinates).
    :param normal_display: A (3,) ndarray of floats holding the triangle's normal (in
        display coordinates), whose depth component mustn't be zero.
    :param offset: The dot product of the triangle's normal with its vertices.
    :param low: The fraction along the stroke where its part within the footprint
        starts.
    :param high: The fraction along the stroke where its part within the footprint ends.
    :return: A tuple of the fractions along the stroke where the hidden part starts and
        ends, or None if the triangle hides none of it.
    """
    gap_constant = (
        start_display[2]
        - (
            offset
            - normal_display[0] * start_display[0]
            - normal_display[1] * start_display[1]
        )
        / normal_display[2]
        - _PLANE_TOLERANCE
    )
    gap_slope = (
        direction_display[2]
        + (
            normal_display[0] * direction_display[0]
            + normal_display[1] * direction_display[1]
        )
        / normal_display[2]
    )
    if abs(gap_slope) < _DEGENERATE_TOLERANCE:
        if gap_constant <= 0.0:
            return None
        hidden_low, hidden_high = low, high
    else:
        crossing = -gap_constant / gap_slope
        if gap_slope > 0.0:
            hidden_low, hidden_high = max(low, crossing), high
        else:
            hidden_low, hidden_high = low, min(high, crossing)
    if hidden_low < hidden_high:
        return hidden_low, hidden_high
    return None


def _subtract_interval(
    intervals: list[tuple[float, float]], low: float, high: float
) -> list[tuple[float, float]]:
    """Removes an interval from a list of disjoint intervals.

    :param intervals: The disjoint intervals, each a tuple of its start and end.
    :param low: The start of the interval to remove.
    :param high: The end of the interval to remove.
    :return: The remaining disjoint intervals, each a tuple of its start and end.
    """
    remaining: list[tuple[float, float]] = []
    for start, end in intervals:
        if high <= start or low >= end:
            remaining.append((start, end))
            continue
        if low > start:
            remaining.append((start, low))
        if high < end:
            remaining.append((high, end))
    return remaining


def get_visible_intervals(
    strokes_display: np.ndarray,
    stroke_owners: Sequence[set[int]],
    triangles_display: np.ndarray,
) -> list[list[tuple[float, float]]]:
    """Finds the parts of each stroke that no triangle hides.

    A triangle hides the part of a stroke that lies within its footprint on screen and
    behind its plane. A stroke is never hidden by the triangles that own it, which are
    the faces it outlines, and it wins a depth tie against any other triangle. A
    triangle seen exactly edge on has no footprint, so it hides nothing.

    Visibility is judged along each stroke's centerline, so a stroke whose centerline
    runs along the boundary of what hides it keeps or loses its whole width at once.

    :param strokes_display: A (N,2,3) ndarray of floats holding each stroke's start and
        end (in display coordinates).
    :param stroke_owners: A sequence of N sets of ints, holding the indices of the
        triangles that own each stroke.
    :param triangles_display: A (T,3,3) ndarray of floats holding each triangle's
        vertices (in display coordinates).
    :return: A list of N lists, each holding a stroke's visible parts as tuples of the
        fractions along it where they start and end.
    """
    occluders = _prepare_occluders(triangles_display)

    # Register each triangle that can hide anything in the cells of a square grid on
    # screen that its bounding box overlaps, so each stroke only needs to consider the
    # triangles registered in the cells its own bounding box overlaps.
    grid: dict[tuple[int, int], list[int]] = {}
    cell_mins = np.floor(occluders.mins_display[:, :2] / _GRID_CELL_SIZE).astype(int)
    cell_maxs = np.floor(occluders.maxs_display[:, :2] / _GRID_CELL_SIZE).astype(int)
    for triangle_id in np.flatnonzero(occluders.usable):
        for cell_x in range(cell_mins[triangle_id, 0], cell_maxs[triangle_id, 0] + 1):
            for cell_y in range(
                cell_mins[triangle_id, 1], cell_maxs[triangle_id, 1] + 1
            ):
                grid.setdefault((cell_x, cell_y), []).append(int(triangle_id))

    visible_intervals: list[list[tuple[float, float]]] = []
    for stroke_display, owners in zip(strokes_display, stroke_owners):
        start_display = stroke_display[0]
        direction_display = stroke_display[1] - stroke_display[0]
        stroke_mins_display = stroke_display.min(axis=0)
        stroke_maxs_display = stroke_display.max(axis=0)

        registered_ids: set[int] = set()
        stroke_cell_mins = np.floor(stroke_mins_display[:2] / _GRID_CELL_SIZE)
        stroke_cell_maxs = np.floor(stroke_maxs_display[:2] / _GRID_CELL_SIZE)
        for cell_x in range(int(stroke_cell_mins[0]), int(stroke_cell_maxs[0]) + 1):
            for cell_y in range(int(stroke_cell_mins[1]), int(stroke_cell_maxs[1]) + 1):
                registered_ids.update(grid.get((cell_x, cell_y), ()))
        registered = np.array(sorted(registered_ids), dtype=int)

        # Only a triangle whose bounding box overlaps the stroke's on screen, and whose
        # nearest vertex is nearer than the stroke's farthest point, can hide any of it.
        candidate_ids = registered[
            (occluders.mins_display[registered, 0] <= stroke_maxs_display[0])
            & (occluders.maxs_display[registered, 0] >= stroke_mins_display[0])
            & (occluders.mins_display[registered, 1] <= stroke_maxs_display[1])
            & (occluders.maxs_display[registered, 1] >= stroke_mins_display[1])
            & (
                occluders.mins_display[registered, 2]
                < stroke_maxs_display[2] - _PLANE_TOLERANCE
            )
        ]

        intervals = [(0.0, 1.0)]
        for triangle_id in candidate_ids:
            if not intervals:
                break
            if int(triangle_id) in owners:
                continue
            footprint = _clip_to_footprint(
                start_display, direction_display, occluders.wound_display[triangle_id]
            )
            if footprint is None:
                continue
            hidden = _get_occluded_range(
                start_display,
                direction_display,
                occluders.normals_display[triangle_id],
                float(occluders.offsets[triangle_id]),
                *footprint,
            )
            if hidden is not None:
                intervals = _subtract_interval(intervals, *hidden)
        visible_intervals.append(intervals)
    return visible_intervals


def _get_ribbon_triangles(
    segments_display: np.ndarray, half_widths: np.ndarray
) -> np.ndarray:
    """Returns the ribbons that straight segments of lines cover on screen, each split
    into two triangles.

    Each ribbon runs along its segment and extends half its width to either side of it,
    perpendicular to it on screen, at the depths of the segment beneath it, so it lies
    in the plane through the segment that faces the camera. A segment that is a point on
    screen gives a ribbon with no area, which hides nothing.

    :param segments_display: A (S,2,3) ndarray of floats holding each segment's start
        and end (in display coordinates).
    :param half_widths: A (S,) ndarray of floats holding half of each ribbon's width, in
        pixels.
    :return: A (2S,3,3) ndarray of floats holding the triangles' vertices (in display
        coordinates). Segment i's ribbon is split into triangles 2i and 2i + 1.
    """
    alongs = segments_display[:, 1, :2] - segments_display[:, 0, :2]
    lengths = np.linalg.norm(alongs, axis=1)
    scales = np.divide(
        half_widths,
        lengths,
        out=np.zeros(lengths.shape[0], dtype=float),
        where=lengths > 0.0,
    )
    acrosses_display = np.zeros((segments_display.shape[0], 3), dtype=float)
    acrosses_display[:, 0] = -alongs[:, 1] * scales
    acrosses_display[:, 1] = alongs[:, 0] * scales
    first_left = segments_display[:, 0] + acrosses_display
    second_left = segments_display[:, 1] + acrosses_display
    second_right = segments_display[:, 1] - acrosses_display
    first_right = segments_display[:, 0] - acrosses_display
    return np.stack(
        [
            np.stack([first_left, second_left, second_right], axis=1),
            np.stack([first_left, second_right, first_right], axis=1),
        ],
        axis=1,
    ).reshape(-1, 3, 3)


def _round_svg_path_data(svg: str) -> str:
    """Returns an svg with the coordinates in its paths' data rounded to
    _SVG_DECIMAL_PLACES decimal places, with their trailing zeros removed.

    Matplotlib writes every coordinate to six decimal places, which makes up most of a
    large diagram's file without changing what is shown. Only the d attributes of path
    elements are rounded, so the text, its embedded fonts, and the styles are left
    alone.

    :param svg: The svg's contents.
    :return: The svg's contents, with its paths' coordinates rounded.
    """

    def round_number(number_match: re.Match[str]) -> str:
        # Rounding a small negative number to zero gives negative zero, which is written
        # without its sign.
        value = round(float(number_match.group()), _SVG_DECIMAL_PLACES) + 0.0
        return f"{value:.{_SVG_DECIMAL_PLACES}f}".rstrip("0").rstrip(".")

    return re.sub(
        r'<path d="[^"]*"',
        lambda path_match: re.sub(r"-?\d+\.\d+", round_number, path_match.group()),
        svg,
    )


class VectorLayer:
    """A layer of a scene, holding its opaque fills, the strokes that outline and run
    among them, its convex occluders, its dots, and its translucent polygons.

    Visibility is resolved within a layer, and the layers are painted in the order they
    were added, each over the ones before it. When a layer is drawn, each convex
    occluder becomes opaque fills and the strokes that outline it. The layer is then
    painted in four passes: its fills from back to front, then its strokes with their
    hidden parts removed, then its dots, and then its translucent polygons. The fills
    hide each other and the strokes behind them. The strokes also hide the strokes
    behind them, across their own widths, so where strokes cross, the nearer one shows.
    Where two strokes are equally deep, neither hides the other, and the one added later
    is painted over the other, which matches VTK, whose depth test lets a later fragment
    replace an equally deep one. The dots and the translucent polygons hide nothing and
    are hidden by nothing.
    """

    def __init__(self) -> None:
        """The initialization method.

        :return: None
        """
        self._listFillTriangles_D_Do: list[np.ndarray] = []
        self._fill_colors: list[np.ndarray] = []

        # Each stroke is one straight segment. A stroke that outlines fills records the
        # indices of the fill triangles it outlines, and the strokes are keyed by their
        # ends, rounded, so an edge shared by neighboring fills is outlined once.
        self._listStrokes_D_Do: list[np.ndarray] = []
        self._stroke_owners: list[set[int]] = []
        self._stroke_colors: list[np.ndarray] = []
        self._stroke_widths: list[float] = []
        self._edge_ids: dict[tuple[tuple[float, ...], ...], int] = {}

        # The polylines' strokes are keyed by their ends, rounded, along with their
        # color and width, so a stroke repeated in the same color and width, such as an
        # edge shared by neighboring closed polylines, is drawn once.
        self._polyline_stroke_keys: set[
            tuple[tuple[tuple[float, ...], ...], tuple[float, ...], float]
        ] = set()

        # Each convex occluder is a tuple of its points (in diagram axes, relative to
        # the diagram origin), its faces, its fill color, its outline color, its outline
        # width, and the indices of the faces whose boundaries are outlined while they
        # face the camera.
        self._occluders: list[
            tuple[
                np.ndarray, list[np.ndarray], np.ndarray, np.ndarray, float, list[int]
            ]
        ] = []

        self._listDots_D_Do: list[np.ndarray] = []
        self._dot_colors: list[np.ndarray] = []
        self._dot_diameters: list[float] = []

        self._listTranslucentPolygons_D_Do: list[np.ndarray] = []
        self._translucent_colors: list[np.ndarray] = []

    def add_quadrilaterals(
        self,
        gridQuadrilaterals_D_Do: np.ndarray,
        fill_colors: np.ndarray,
        edge_color: matplotlib.typing.ColorType,
        edge_width: float,
    ) -> None:
        """Adds opaque quadrilateral fills, outlined along their edges.

        Each quadrilateral is split into two triangles along its first diagonal. An edge
        shared by neighboring quadrilaterals in this layer is outlined once, in the
        color and width it was first added with.

        :param gridQuadrilaterals_D_Do: A (M,4,3) ndarray of floats holding each
            quadrilateral's corners, in order around it (in diagram axes, relative to
            the diagram origin). The units are in meters.
        :param fill_colors: A (M,4) ndarray of floats holding each quadrilateral's RGBA
            fill color, with each channel in the range [0.0, 1.0].
        :param edge_color: The color of the outlines, as any color Matplotlib accepts.
        :param edge_width: The width of the outlines, in pixels.
        :return: None
        """
        edge_rgba = np.array(matplotlib.colors.to_rgba(edge_color), dtype=float)
        for quadrilateral_D_Do, fill_color in zip(gridQuadrilaterals_D_Do, fill_colors):
            first_triangle_id = len(self._listFillTriangles_D_Do)
            self._listFillTriangles_D_Do += [
                quadrilateral_D_Do[[0, 1, 2]],
                quadrilateral_D_Do[[0, 2, 3]],
            ]
            self._fill_colors += [np.array(fill_color, dtype=float)] * 2
            owners = {first_triangle_id, first_triangle_id + 1}
            for corner_id in range(4):
                stroke_D_Do = quadrilateral_D_Do[[corner_id, (corner_id + 1) % 4]]
                key = tuple(sorted(tuple(np.round(end, 9)) for end in stroke_D_Do))
                edge_id = self._edge_ids.get(key)
                if edge_id is not None:
                    self._stroke_owners[edge_id].update(owners)
                    continue
                self._edge_ids[key] = len(self._listStrokes_D_Do)
                self._listStrokes_D_Do.append(stroke_D_Do)
                self._stroke_owners.append(set(owners))
                self._stroke_colors.append(edge_rgba)
                self._stroke_widths.append(edge_width)

    def add_polylines(
        self,
        listPolylines_D_Do: Sequence[np.ndarray],
        color: matplotlib.typing.ColorType,
        width: float,
        closed: bool = False,
    ) -> None:
        """Adds polylines, which the layer's fills and convex occluders can hide.

        A segment that repeats one already added by this method in the same color and
        width, such as an edge shared by neighboring closed polylines, is added once.

        :param listPolylines_D_Do: A sequence of (K,3) ndarrays of floats, each holding
            a polyline's points in order along it (in diagram axes, relative to the
            diagram origin). The units are in meters.
        :param color: The color of the polylines, as any color Matplotlib accepts.
        :param width: The width of the polylines, in pixels.
        :param closed: Determines whether each polyline is closed, joining its last
            point back to its first. The default is False.
        :return: None
        """
        rgba = np.array(matplotlib.colors.to_rgba(color), dtype=float)
        for polyline_D_Do in listPolylines_D_Do:
            if closed:
                polyline_D_Do = np.vstack([polyline_D_Do, polyline_D_Do[:1]])
            for point_id in range(polyline_D_Do.shape[0] - 1):
                stroke_D_Do = polyline_D_Do[point_id : point_id + 2]
                key = (
                    tuple(sorted(tuple(np.round(end, 9)) for end in stroke_D_Do)),
                    tuple(rgba),
                    width,
                )
                if key in self._polyline_stroke_keys:
                    continue
                self._polyline_stroke_keys.add(key)
                self._listStrokes_D_Do.append(stroke_D_Do)
                self._stroke_owners.append(set())
                self._stroke_colors.append(rgba)
                self._stroke_widths.append(width)

    def add_convex_occluder(
        self,
        stackPoints_D_Do: np.ndarray,
        faces: Sequence[np.ndarray],
        fill_color: matplotlib.typing.ColorType,
        outline_color: matplotlib.typing.ColorType,
        outline_width: float,
        outlined_face_ids: Sequence[int] = (),
    ) -> None:
        """Adds an opaque convex polyhedron, which is outlined along its silhouette and
        hides the fills and strokes behind it.

        Its silhouette is the boundary of its footprint on screen, which, since it is
        convex, is the convex hull of its points' positions on screen. The boundaries of
        the faces in outlined_face_ids are outlined too, while those faces face the
        camera. Its faces are ordered against the layer's fills, and its outlines
        against the layer's strokes, so it is hidden by whatever is nearer.

        :param stackPoints_D_Do: A (P,3) ndarray of floats holding the polyhedron's
            points (in diagram axes, relative to the diagram origin). The units are in
            meters.
        :param faces: A sequence of (K,) ndarrays of ints, each holding the indices of a
            face's points, in order around it. The faces must be convex.
        :param fill_color: The polyhedron's fill color, as any color Matplotlib accepts.
        :param outline_color: The color of its outlines, as any color Matplotlib
            accepts.
        :param outline_width: The width of its outlines, in pixels.
        :param outlined_face_ids: The indices of the faces whose boundaries are outlined
            while they face the camera. The default is an empty sequence.
        :return: None
        """
        self._occluders.append(
            (
                np.array(stackPoints_D_Do, dtype=float),
                [np.array(face, dtype=int) for face in faces],
                np.array(matplotlib.colors.to_rgba(fill_color), dtype=float),
                np.array(matplotlib.colors.to_rgba(outline_color), dtype=float),
                outline_width,
                list(outlined_face_ids),
            )
        )

    def add_dots(
        self,
        stackDots_D_Do: np.ndarray,
        color: matplotlib.typing.ColorType,
        diameter: float,
    ) -> None:
        """Adds round dots of a fixed size on screen.

        :param stackDots_D_Do: A (N,3) ndarray of floats holding each dot's center (in
            diagram axes, relative to the diagram origin). The units are in meters.
        :param color: The color of the dots, as any color Matplotlib accepts.
        :param diameter: The diameter of the dots, in pixels.
        :return: None
        """
        rgba = np.array(matplotlib.colors.to_rgba(color), dtype=float)
        for dot_D_Do in stackDots_D_Do:
            self._listDots_D_Do.append(np.array(dot_D_Do, dtype=float))
            self._dot_colors.append(rgba)
            self._dot_diameters.append(diameter)

    def add_translucent_polygons(
        self, listPolygons_D_Do: Sequence[np.ndarray], colors: np.ndarray
    ) -> None:
        """Adds translucent polygons, which are painted over the rest of the layer.

        :param listPolygons_D_Do: A sequence of N (K,3) ndarrays of floats, each holding
            a polygon's corners in order around it (in diagram axes, relative to the
            diagram origin). The units are in meters.
        :param colors: A (N,4) ndarray of floats holding each polygon's RGBA color, with
            each channel in the range [0.0, 1.0].
        :return: None
        """
        for polygon_D_Do, color in zip(listPolygons_D_Do, colors):
            self._listTranslucentPolygons_D_Do.append(
                np.array(polygon_D_Do, dtype=float)
            )
            self._translucent_colors.append(np.array(color, dtype=float))

    def draw(
        self,
        axes: matplotlib.axes.Axes,
        camera: VectorCamera,
        zorder: float,
        line_width_scale: float = 1.0,
    ) -> float:
        """Draws the layer onto a figure's axes, which are laid out one unit per pixel
        of the window.

        :param axes: The Axes to draw the layer onto.
        :param camera: The VectorCamera to project the layer through.
        :param zorder: The zorder of the layer's first pass. Each pass takes the next
            whole zorder.
        :param line_width_scale: The factor that the widths of the strokes, including
            the convex occluders' outlines, are scaled by, both where they are drawn and
            where they hide the strokes behind them. It must be positive. The default is
            1.0.
        :return: The zorder for whatever is drawn after the layer.
        """
        (
            listFillTriangles_D_Do,
            fill_colors,
            listStrokes_D_Do,
            base_stroke_owners,
            base_stroke_colors,
            base_stroke_widths,
        ) = self._get_fills_and_strokes(camera)

        fill_triangles_display = np.zeros((0, 3, 3), dtype=float)
        if listFillTriangles_D_Do:
            fill_triangles_display = camera.to_display(
                np.concatenate(listFillTriangles_D_Do)
            ).reshape(-1, 3, 3)
            order = get_paint_order(
                [
                    (triangle_display, triangle_id)
                    for triangle_id, triangle_display in enumerate(
                        fill_triangles_display
                    )
                ]
            )
            # Each run of consecutive pieces of the same color is written as one
            # compound path, rather than one path per piece, which keeps the file small.
            # Merging a run doesn't change what is painted, since its pieces are all
            # opaque and the same color. The files are filled by the nonzero rule, so
            # every piece is wound counterclockwise, which keeps overlapping pieces from
            # cancelling each other out.
            fill_paths: list[matplotlib.path.Path] = []
            fill_path_colors: list[np.ndarray] = []
            run_pieces: list[matplotlib.path.Path] = []
            for piece_id, (piece_display, fill_id) in enumerate(order):
                corners_display = piece_display[:, :2]
                signed_double_area = np.dot(
                    corners_display[:, 0], np.roll(corners_display[:, 1], -1)
                ) - np.dot(np.roll(corners_display[:, 0], -1), corners_display[:, 1])
                if signed_double_area < 0.0:
                    corners_display = corners_display[::-1]
                run_pieces.append(
                    matplotlib.path.Path(
                        np.vstack([corners_display, corners_display[:1]]), closed=True
                    )
                )
                color = fill_colors[fill_id]
                if piece_id == len(order) - 1 or not np.array_equal(
                    color, fill_colors[order[piece_id + 1][1]]
                ):
                    fill_paths.append(
                        matplotlib.path.Path.make_compound_path(*run_pieces)
                    )
                    fill_path_colors.append(color)
                    run_pieces = []
            axes.add_collection(
                matplotlib.collections.PathCollection(
                    fill_paths,
                    facecolors=fill_path_colors,
                    edgecolors=fill_path_colors,
                    linewidths=_FILL_SEAM_LINE_WIDTH * POINTS_PER_PIXEL,
                    zorder=zorder,
                )
            )
        zorder += 1.0

        if listStrokes_D_Do:
            strokes_display = camera.to_display(
                np.concatenate(listStrokes_D_Do)
            ).reshape(-1, 2, 3)
            stroke_widths = line_width_scale * np.array(base_stroke_widths, dtype=float)
            num_strokes = strokes_display.shape[0]

            # Every stroke hides the strokes behind it with a ribbon of its own width,
            # in the plane through it that faces the camera. The ribbons follow the
            # fills' triangles, two per stroke. A stroke is never hidden by its own
            # ribbon, or by the ribbon of a stroke that shares an end with it, such as
            # its neighbor along a polyline, which would otherwise clip it where they
            # meet.
            first_stroke_ribbon_id = fill_triangles_display.shape[0]
            stroke_ids_by_end: dict[tuple[float, ...], list[int]] = {}
            for stroke_id, stroke_display in enumerate(strokes_display):
                for end_display in stroke_display:
                    stroke_ids_by_end.setdefault(
                        tuple(np.round(end_display, 6)), []
                    ).append(stroke_id)
            stroke_owners: list[set[int]] = []
            for stroke_id, stroke_display in enumerate(strokes_display):
                connected_ids = {stroke_id}
                for end_display in stroke_display:
                    connected_ids.update(
                        stroke_ids_by_end[tuple(np.round(end_display, 6))]
                    )
                stroke_owners.append(
                    base_stroke_owners[stroke_id]
                    | {
                        first_stroke_ribbon_id + 2 * connected_id + half
                        for connected_id in connected_ids
                        for half in (0, 1)
                    }
                )

            # A clipped stroke's round end reaches half its width past where it is cut,
            # so the ribbons that clip it are widened by half its width, which keeps
            # that end from showing over the stroke in front. The strokes are clipped
            # one width at a time, against ribbons widened for that width.
            visible_intervals: list[list[tuple[float, float]]] = [[]] * num_strokes
            for width in np.unique(stroke_widths):
                width_stroke_ids = np.flatnonzero(stroke_widths == width)
                width_intervals = get_visible_intervals(
                    strokes_display[width_stroke_ids],
                    [stroke_owners[stroke_id] for stroke_id in width_stroke_ids],
                    np.concatenate(
                        [
                            fill_triangles_display,
                            _get_ribbon_triangles(
                                strokes_display, 0.5 * (stroke_widths + width)
                            ),
                        ]
                    ),
                )
                for width_stroke_id, intervals in zip(
                    width_stroke_ids, width_intervals
                ):
                    visible_intervals[int(width_stroke_id)] = intervals

            # Where two strokes are equally deep, neither hides the other, and the one
            # added later is painted over the other, as VTK draws them.
            segments: list[np.ndarray] = []
            segment_colors: list[np.ndarray] = []
            segment_widths: list[float] = []
            for stroke_id, intervals in enumerate(visible_intervals):
                start_display = strokes_display[stroke_id, 0]
                direction_display = (
                    strokes_display[stroke_id, 1] - strokes_display[stroke_id, 0]
                )
                for low, high in intervals:
                    if high - low <= _MIN_VISIBLE_FRACTION:
                        continue
                    segments.append(
                        np.array(
                            [
                                start_display[:2] + low * direction_display[:2],
                                start_display[:2] + high * direction_display[:2],
                            ],
                            dtype=float,
                        )
                    )
                    segment_colors.append(base_stroke_colors[stroke_id])
                    segment_widths.append(
                        float(stroke_widths[stroke_id]) * POINTS_PER_PIXEL
                    )

            # Each run of consecutive opaque segments of the same color and width is
            # written as one compound path, rather than one path per segment, which
            # keeps the file small. Translucent segments are each written alone, since a
            # compound path paints where its segments overlap only once. Matplotlib
            # simplifies a long path of straight pieces by dropping some of its points,
            # which is turned off so that every segment is kept.
            stroke_paths: list[matplotlib.path.Path] = []
            stroke_path_colors: list[np.ndarray] = []
            stroke_path_widths: list[float] = []
            run_segments: list[np.ndarray] = []
            for segment_id, segment_display in enumerate(segments):
                run_segments.append(segment_display)
                color = segment_colors[segment_id]
                width = segment_widths[segment_id]
                if (
                    segment_id == len(segments) - 1
                    or color[3] < 1.0
                    or not np.array_equal(color, segment_colors[segment_id + 1])
                    or width != segment_widths[segment_id + 1]
                ):
                    stroke_path = matplotlib.path.Path(
                        np.concatenate(run_segments),
                        np.tile(
                            [matplotlib.path.Path.MOVETO, matplotlib.path.Path.LINETO],
                            len(run_segments),
                        ),
                    )
                    stroke_path.should_simplify = False
                    stroke_paths.append(stroke_path)
                    stroke_path_colors.append(color)
                    stroke_path_widths.append(width)
                    run_segments = []
            axes.add_collection(
                matplotlib.collections.PathCollection(
                    stroke_paths,
                    facecolors="none",
                    edgecolors=stroke_path_colors,
                    linewidths=stroke_path_widths,
                    capstyle="round",
                    joinstyle="round",
                    zorder=zorder,
                )
            )
        zorder += 1.0

        if self._listDots_D_Do:
            dots_display = camera.to_display(np.array(self._listDots_D_Do, dtype=float))
            axes.add_collection(
                matplotlib.collections.PatchCollection(
                    [
                        matplotlib.patches.Circle(
                            (float(dot_display[0]), float(dot_display[1])),
                            0.5 * diameter,
                        )
                        for dot_display, diameter in zip(
                            dots_display, self._dot_diameters
                        )
                    ],
                    facecolors=self._dot_colors,
                    edgecolors="none",
                    zorder=zorder,
                )
            )
        zorder += 1.0

        if self._listTranslucentPolygons_D_Do:
            axes.add_collection(
                matplotlib.collections.PolyCollection(
                    [
                        camera.to_display(polygon_D_Do)[:, :2]
                        for polygon_D_Do in self._listTranslucentPolygons_D_Do
                    ],
                    facecolors=self._translucent_colors,
                    edgecolors="none",
                    linewidths=0.0,
                    zorder=zorder,
                )
            )
        return zorder + 1.0

    def _get_fills_and_strokes(self, camera: VectorCamera) -> tuple[
        list[np.ndarray],
        list[np.ndarray],
        list[np.ndarray],
        list[set[int]],
        list[np.ndarray],
        list[float],
    ]:
        """Returns the layer's fills and strokes, with its convex occluders converted
        into more of each as seen through a camera.

        Each occluder's faces, split into triangles, become opaque fills, which the BSP
        then orders exactly against the layer's other fills and against each other,
        however they interpenetrate. Its outline becomes strokes, owned by its own
        triangles, so they are hidden by whatever is nearer, but never by the occluder
        they outline. The outline runs along its silhouette, which is the boundary of
        its footprint on screen and, since it is convex, the convex hull of its points'
        positions on screen, and along the boundaries of its outlined faces that face
        the camera. The layer itself isn't changed.

        :param camera: The VectorCamera that decides which of the occluders' faces face
            the camera, and where their silhouettes are.
        :return: A tuple of six lists. The first holds each fill triangle as a (3,3)
            ndarray of floats (in diagram axes, relative to the diagram origin), and the
            second holds each fill triangle's RGBA color as a (4,) ndarray of floats.
            The third holds each stroke's ends as a (2,3) ndarray of floats (in diagram
            axes, relative to the diagram origin), the fourth holds the set of indices
            of the fill triangles each stroke outlines, the fifth holds each stroke's
            RGBA color as a (4,) ndarray of floats, and the sixth holds each stroke's
            width, in pixels. The units of the positions are in meters.
        """
        listFillTriangles_D_Do = list(self._listFillTriangles_D_Do)
        fill_colors = list(self._fill_colors)
        listStrokes_D_Do = list(self._listStrokes_D_Do)
        stroke_owners = [set(owners) for owners in self._stroke_owners]
        stroke_colors = list(self._stroke_colors)
        stroke_widths = list(self._stroke_widths)

        for (
            stackPoints_D_Do,
            faces,
            fill_rgba,
            outline_rgba,
            outline_width,
            outlined_face_ids,
        ) in self._occluders:
            first_triangle_id = len(listFillTriangles_D_Do)
            for face in faces:
                for corner_id in range(1, face.shape[0] - 1):
                    listFillTriangles_D_Do.append(
                        stackPoints_D_Do[
                            [face[0], face[corner_id], face[corner_id + 1]]
                        ]
                    )
                    fill_colors.append(fill_rgba)
            owners = set(range(first_triangle_id, len(listFillTriangles_D_Do)))

            points_display = camera.to_display(stackPoints_D_Do)
            try:
                hull = scipy.spatial.ConvexHull(points_display[:, :2])
            except scipy.spatial.QhullError:
                # A polyhedron whose footprint on screen has no area shows no outline.
                continue
            outlines = [hull.vertices]

            # A face of a convex polyhedron faces the camera where its outward normal
            # points back toward negative depth. Each face's normal is turned outward,
            # away from the polyhedron's center, so the faces may wind either way.
            center_display = points_display.mean(axis=0)
            for face_id in outlined_face_ids:
                face_display = points_display[faces[face_id]]
                normal_display = _get_newell_normal(face_display)
                if normal_display @ (face_display.mean(axis=0) - center_display) < 0.0:
                    normal_display = -normal_display
                if normal_display[2] < 0.0:
                    outlines.append(faces[face_id])

            # Each closed outline becomes a stroke between each pair of its consecutive
            # points.
            for outline in outlines:
                for point_id in range(len(outline)):
                    listStrokes_D_Do.append(
                        stackPoints_D_Do[
                            [outline[point_id], outline[(point_id + 1) % len(outline)]]
                        ]
                    )
                    stroke_owners.append(set(owners))
                    stroke_colors.append(outline_rgba)
                    stroke_widths.append(outline_width)

        return (
            listFillTriangles_D_Do,
            fill_colors,
            listStrokes_D_Do,
            stroke_owners,
            stroke_colors,
            stroke_widths,
        )


class VectorScene:
    """A scene recorded for export as an svg or pdf file.

    A scene holds its layers, which are painted in the order they were added, followed
    by its screen polygons and then its texts, which are placed in the window rather
    than in diagram axes.
    """

    def __init__(self) -> None:
        """The initialization method.

        :return: None
        """
        self._layers: list[VectorLayer] = []
        self._listScreenPolygons: list[np.ndarray] = []
        self._screen_polygon_colors: list[np.ndarray] = []
        self._texts: list[VectorText] = []

    def add_layer(self) -> VectorLayer:
        """Adds a layer, which is painted over every layer added before it.

        :return: The new VectorLayer.
        """
        layer = VectorLayer()
        self._layers.append(layer)
        return layer

    def add_screen_polygons(
        self, listPolygons: Sequence[np.ndarray], colors: np.ndarray
    ) -> None:
        """Adds opaque polygons placed in the window, which are painted over every
        layer.

        Each polygon is outlined in its own color, so abutting polygons show no seams
        between them.

        :param listPolygons: A sequence of N (K,2) ndarrays of floats, each holding a
            polygon's corners in order around it, in pixels from the window's bottom
            left corner.
        :param colors: A (N,4) ndarray of floats holding each polygon's RGBA color, with
            each channel in the range [0.0, 1.0].
        :return: None
        """
        for polygon, color in zip(listPolygons, colors):
            self._listScreenPolygons.append(np.array(polygon, dtype=float))
            self._screen_polygon_colors.append(np.array(color, dtype=float))

    def add_text(self, text: VectorText) -> None:
        """Adds a run of text, which is drawn over everything else.

        :param text: The VectorText to add.
        :return: None
        """
        self._texts.append(text)

    def save(
        self,
        path: Path,
        camera: VectorCamera,
        background_color: matplotlib.typing.ColorType | None,
        line_width_scale: float = 1.0,
        selectable_text: bool = True,
    ) -> None:
        """Saves the scene as an svg or pdf file, as seen through a camera.

        The format follows the path's suffix. A pdf embeds its fonts as TrueType, so its
        text stays selectable and searchable. By default, an svg writes its text as text
        with the fonts it uses embedded in it, as the results plots do, which also keeps
        it selectable and searchable. Some programs that import an svg ignore its
        embedded fonts, though, and draw its text in other fonts, which misplaces the
        parts of any math. An svg can instead draw each character of its text as a
        filled outline, which looks the same in every program, but can't be selected or
        searched.

        :param path: The path of the file to write. It must end with ".svg" or ".pdf",
            and its directory must already exist.
        :param camera: The VectorCamera to project the scene through, which also sets
            the size of the page.
        :param background_color: The color of the page, as any color Matplotlib accepts,
            or None to leave it transparent.
        :param line_width_scale: The factor that the widths of every layer's strokes and
            convex occluders' outlines are scaled by. It must be positive. The default
            is 1.0.
        :param selectable_text: Determines whether an svg writes its text as text, with
            its fonts embedded, rather than as filled outlines. It has no effect on a
            pdf. The default is True.
        :return: None
        """
        with matplotlib.rc_context(
            {
                "pdf.fonttype": 42,
                "svg.fonttype": "none" if selectable_text else "path",
            }
        ):
            figure = matplotlib.figure.Figure(
                figsize=(
                    camera.window_width / _FIGURE_DPI,
                    camera.window_height / _FIGURE_DPI,
                ),
                dpi=_FIGURE_DPI,
            )
            axes = figure.add_axes((0.0, 0.0, 1.0, 1.0))
            axes.set_axis_off()
            axes.set_xlim(0.0, camera.window_width)
            axes.set_ylim(0.0, camera.window_height)

            zorder = 1.0
            for layer in self._layers:
                zorder = layer.draw(axes, camera, zorder, line_width_scale)

            if self._listScreenPolygons:
                axes.add_collection(
                    matplotlib.collections.PolyCollection(
                        self._listScreenPolygons,
                        facecolors=self._screen_polygon_colors,
                        edgecolors=self._screen_polygon_colors,
                        linewidths=_FILL_SEAM_LINE_WIDTH * POINTS_PER_PIXEL,
                        zorder=zorder,
                    )
                )
            zorder += 1.0

            # Each text takes its own zorder, after its background box's, so a later
            # text's box covers an earlier text it overlaps, as it does on screen.
            for text in self._texts:
                if text.background_box is not None:
                    left, right, bottom, top = text.background_box
                    axes.add_patch(
                        matplotlib.patches.Rectangle(
                            (left, bottom),
                            right - left,
                            top - bottom,
                            facecolor=text.background_color,
                            edgecolor="none",
                            linewidth=0.0,
                            zorder=zorder,
                        )
                    )
                font_properties = text.font_properties.copy()
                font_properties.set_size(text.font_size * POINTS_PER_PIXEL)
                matplotlib_text = axes.text(
                    text.x,
                    text.y,
                    text.text,
                    fontproperties=font_properties,
                    color=text.color,
                    horizontalalignment=text.horizontal_alignment,
                    verticalalignment=text.vertical_alignment,
                    clip_on=False,
                    zorder=zorder + 0.5,
                )
                if text.math_font_family is not None:
                    matplotlib_text.set_math_fontfamily(text.math_font_family)
                zorder += 1.0

            # An svg is written to a buffer first, so the fonts can be embedded in it
            # and its paths' coordinates rounded before it reaches the file. An svg
            # whose text is drawn as outlines has no text to embed fonts for.
            facecolor = "none" if background_color is None else background_color
            if path.suffix.lower() == ".svg":
                svg_buffer = io.BytesIO()
                figure.savefig(svg_buffer, format="svg", facecolor=facecolor)
                svg = svg_buffer.getvalue().decode("utf-8")
                if selectable_text:
                    svg = _fonts.embed_fonts_in_svg(svg)
                svg = _round_svg_path_data(svg)
                path.write_bytes(svg.encode("utf-8"))
            else:
                figure.savefig(path, facecolor=facecolor)
