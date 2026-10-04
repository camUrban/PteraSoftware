"""This module contains classes to test the vector export classes and functions."""

import tempfile
import unittest
from pathlib import Path

import matplotlib.collections
import matplotlib.colors
import matplotlib.figure
import matplotlib.font_manager
import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _fonts, _vector_export
from tests.unit.fixtures import vector_export_fixtures


def _get_depth_at(piece_display: np.ndarray, x: float, y: float) -> float:
    """Returns the depth of a planar piece's plane beneath a point on screen.

    :param piece_display: A (K,3) ndarray of floats holding the piece's vertices (in
        display coordinates).
    :param x: The point's x position, in pixels.
    :param y: The point's y position, in pixels.
    :return: The depth of the piece's plane beneath the point, in pixels.
    """
    normal_display = np.cross(
        piece_display[1] - piece_display[0], piece_display[2] - piece_display[0]
    )
    return float(
        piece_display[0, 2]
        - (
            normal_display[0] * (x - piece_display[0, 0])
            + normal_display[1] * (y - piece_display[0, 1])
        )
        / normal_display[2]
    )


def _covers(piece_display: np.ndarray, x: float, y: float) -> bool:
    """Returns whether a convex piece's footprint on screen strictly contains a point.

    :param piece_display: A (K,3) ndarray of floats holding the piece's vertices (in
        display coordinates), in order around it.
    :param x: The point's x position, in pixels.
    :param y: The point's y position, in pixels.
    :return: True if the point is strictly inside the footprint, and False otherwise.
    """
    corners = piece_display[:, :2]
    edges = np.roll(corners, -1, axis=0) - corners
    to_point = np.array([x, y], dtype=float) - corners
    crosses = edges[:, 0] * to_point[:, 1] - edges[:, 1] * to_point[:, 0]
    return bool(np.all(crosses > 1.0e-9) or np.all(crosses < -1.0e-9))


class TestVectorCameraToDisplay(unittest.TestCase):
    """This class contains methods for testing
    _vector_export.VectorCamera.to_display."""

    def setUp(self) -> None:
        """Set up the camera.

        :return: None
        """
        self.camera = vector_export_fixtures.make_top_down_camera_fixture()

    def test_shows_the_focal_point_at_the_window_center(self) -> None:
        """Test that the focal point projects to the window's center at zero depth."""
        npt.assert_allclose(
            self.camera.to_display(np.zeros((1, 3), dtype=float)),
            [[100.0, 50.0, 0.0]],
        )

    def test_scales_by_the_window_height_over_the_parallel_scale(self) -> None:
        """Test that a meter spans the window's height over twice the parallel scale."""
        npt.assert_allclose(
            self.camera.to_display(np.array([[1.0, 2.0, 0.0]], dtype=float)),
            [[110.0, 70.0, 0.0]],
        )

    def test_depth_grows_away_from_the_camera(self) -> None:
        """Test that a point farther from the camera has a larger depth."""
        npt.assert_allclose(
            self.camera.to_display(np.array([[0.0, 0.0, -3.0]], dtype=float)),
            [[100.0, 50.0, 30.0]],
        )

    def test_uses_only_the_view_up_perpendicular_to_the_view(self) -> None:
        """Test that a view up direction tilted toward the view direction, and not of
        unit length, projects the same as the perpendicular unit one."""
        tilted_camera = self.camera._replace(
            viewUp_D=np.array([0.0, 3.0, 4.0], dtype=float)
        )
        stackPoints_D_Do = np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, -1.0]], dtype=float)
        npt.assert_allclose(
            tilted_camera.to_display(stackPoints_D_Do),
            self.camera.to_display(stackPoints_D_Do),
        )


class TestVectorCameraFitToPage(unittest.TestCase):
    """This class contains methods for testing
    _vector_export.VectorCamera.fit_to_page."""

    def setUp(self) -> None:
        """Set up the camera.

        :return: None
        """
        self.camera = vector_export_fixtures.make_top_down_camera_fixture()

    def test_scales_the_view_to_a_page_of_the_same_aspect_ratio(self) -> None:
        """Test that a page twice the window's size, at 100 pixels per inch, shows the
        same view at twice the pixels per meter, centered on the focal point."""
        page_camera = self.camera.fit_to_page(4.0, 2.0)
        self.assertEqual(
            (page_camera.window_width, page_camera.window_height), (400.0, 200.0)
        )
        npt.assert_allclose(
            page_camera.to_display(
                np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0]], dtype=float)
            ),
            [[200.0, 100.0, 0.0], [220.0, 140.0, 0.0]],
        )

    def test_fits_a_wider_view_to_a_taller_page(self) -> None:
        """Test that a view wider than its page is scaled so that its width spans the
        page, and is centered on the page vertically.

        The window spans 20.0 m by 10.0 m, and the page is 200 by 200 pixels, so the
        view keeps its 10.0 pixels per meter, and its right edge, 10.0 m from the focal
        point, lands on the page's right edge.
        """
        page_camera = self.camera.fit_to_page(2.0, 2.0)
        npt.assert_allclose(
            page_camera.to_display(
                np.array([[10.0, 0.0, 0.0], [0.0, 5.0, 0.0]], dtype=float)
            ),
            [[200.0, 100.0, 0.0], [100.0, 150.0, 0.0]],
        )


class TestGetPaintOrder(unittest.TestCase):
    """This class contains methods for testing _vector_export.get_paint_order."""

    def assert_painted_back_to_front(
        self, order: list[tuple[np.ndarray, int]], num_samples: int
    ) -> None:
        """Asserts that, at sampled points on screen, the pieces covering each point are
        painted in order of decreasing depth.

        The samples are drawn uniformly within each piece in turn, so every piece is
        sampled.

        :param order: The pieces in paint order, each a tuple of a (K,3) ndarray of
            floats holding its vertices (in display coordinates) and its int ID.
        :param num_samples: The number of samples to draw within each piece.
        :return: None
        """
        rng = np.random.default_rng(1)
        for piece_display, _ in order:
            for _ in range(num_samples):
                weights = rng.dirichlet(np.ones(piece_display.shape[0]))
                x, y = weights @ piece_display[:, :2]
                depths = [
                    _get_depth_at(other_display, x, y)
                    for other_display, _ in order
                    if _covers(other_display, x, y)
                ]
                for earlier, later in zip(depths[:-1], depths[1:]):
                    self.assertLessEqual(later, earlier + 1.0e-6)

    def test_returns_nothing_for_no_polygons(self) -> None:
        """Test that no polygons give no pieces."""
        self.assertEqual(_vector_export.get_paint_order([]), [])

    def test_splits_crossing_polygons(self) -> None:
        """Test that two triangles that cross each other are split, and that each piece
        keeps the ID of the triangle it was cut from."""
        order = _vector_export.get_paint_order(
            vector_export_fixtures.make_crossing_triangles_fixture()
        )
        self.assertGreater(len(order), 2)
        self.assertEqual({piece_id for _, piece_id in order}, {0, 1})

    def test_paints_crossing_polygons_back_to_front(self) -> None:
        """Test that the pieces of two crossing triangles are painted from back to front
        everywhere they overlap."""
        order = _vector_export.get_paint_order(
            vector_export_fixtures.make_crossing_triangles_fixture()
        )
        self.assert_painted_back_to_front(order, num_samples=20)

    def test_paints_random_polygons_back_to_front(self) -> None:
        """Test that the pieces of many overlapping and crossing triangles are painted
        from back to front everywhere they overlap."""
        order = _vector_export.get_paint_order(
            vector_export_fixtures.make_random_triangles_fixture()
        )
        self.assert_painted_back_to_front(order, num_samples=5)

    def test_keeps_the_area_of_every_polygon(self) -> None:
        """Test that each triangle's pieces cover the same area on screen as the
        triangle did."""
        triangles = vector_export_fixtures.make_random_triangles_fixture()
        order = _vector_export.get_paint_order(triangles)
        for triangle_display, triangle_id in triangles:
            piece_area = sum(
                0.5 * np.linalg.norm(_vector_export._get_newell_normal(piece_display))
                for piece_display, piece_id in order
                if piece_id == triangle_id
            )
            self.assertAlmostEqual(
                piece_area,
                0.5
                * float(
                    np.linalg.norm(_vector_export._get_newell_normal(triangle_display))
                ),
                delta=1.0e-3,
            )

    def test_is_deterministic(self) -> None:
        """Test that the same polygons are always split and ordered the same way."""
        triangles = vector_export_fixtures.make_random_triangles_fixture()
        first_order = _vector_export.get_paint_order(triangles)
        second_order = _vector_export.get_paint_order(triangles)
        self.assertEqual(len(first_order), len(second_order))
        for (first_display, first_id), (second_display, second_id) in zip(
            first_order, second_order
        ):
            self.assertEqual(first_id, second_id)
            npt.assert_array_equal(first_display, second_display)


class TestGetVisibleIntervals(unittest.TestCase):
    """This class contains methods for testing _vector_export.get_visible_intervals."""

    def setUp(self) -> None:
        """Set up the triangle that hides the strokes.

        :return: None
        """
        self.triangles_display = vector_export_fixtures.make_flat_triangle_fixture()

    def get_intervals(
        self, stroke_display: list[list[float]], owners: set[int] | None = None
    ) -> list[tuple[float, float]]:
        """Returns the visible parts of one stroke behind or in front of the triangle.

        :param stroke_display: The stroke's start and end (in display coordinates).
        :param owners: The indices of the triangles that own the stroke. None means it
            has no owners.
        :return: The stroke's visible parts, as tuples of the fractions along it where
            they start and end.
        """
        return _vector_export.get_visible_intervals(
            np.array([stroke_display], dtype=float),
            [owners if owners is not None else set()],
            self.triangles_display,
        )[0]

    def test_keeps_a_stroke_in_front(self) -> None:
        """Test that a stroke nearer than the triangle is wholly visible."""
        intervals = self.get_intervals([[10.0, 10.0, 20.0], [30.0, 10.0, 20.0]])
        self.assertEqual(intervals, [(0.0, 1.0)])

    def test_hides_a_stroke_behind(self) -> None:
        """Test that a stroke behind the triangle, within its footprint, is wholly
        hidden."""
        intervals = self.get_intervals([[10.0, 10.0, 80.0], [30.0, 10.0, 80.0]])
        self.assertEqual(intervals, [])

    def test_hides_the_part_of_a_stroke_within_the_footprint(self) -> None:
        """Test that a stroke behind the triangle stays visible where it leaves the
        triangle's footprint."""
        intervals = self.get_intervals([[-50.0, 10.0, 80.0], [50.0, 10.0, 80.0]])
        self.assertEqual(len(intervals), 1)
        npt.assert_allclose(intervals[0], (0.0, 0.5))

    def test_hides_the_part_of_a_stroke_behind_the_plane(self) -> None:
        """Test that a stroke that passes through the triangle is hidden from where it
        passes behind it.

        The stroke wins depth ties within a small tolerance, so it stays visible a
        sliver past the exact crossing, which is halfway along it.
        """
        intervals = self.get_intervals([[10.0, 10.0, 20.0], [30.0, 10.0, 80.0]])
        self.assertEqual(len(intervals), 1)
        npt.assert_allclose(intervals[0], (0.0, 0.5), atol=1.0e-4)

    def test_never_hides_a_stroke_behind_its_owner(self) -> None:
        """Test that a stroke is never hidden by a triangle that owns it."""
        intervals = self.get_intervals(
            [[10.0, 10.0, 80.0], [30.0, 10.0, 80.0]], owners={0}
        )
        self.assertEqual(intervals, [(0.0, 1.0)])

    def test_a_stroke_wins_a_depth_tie(self) -> None:
        """Test that a stroke lying in the triangle's plane stays visible."""
        intervals = self.get_intervals([[10.0, 10.0, 50.0], [30.0, 10.0, 50.0]])
        self.assertEqual(intervals, [(0.0, 1.0)])

    def test_an_edge_on_triangle_hides_nothing(self) -> None:
        """Test that a triangle seen exactly edge on hides nothing."""
        self.triangles_display = np.array(
            [[[0.0, 10.0, 0.0], [100.0, 10.0, 0.0], [0.0, 10.0, 100.0]]], dtype=float
        )
        intervals = self.get_intervals([[10.0, 10.0, 80.0], [30.0, 10.0, 80.0]])
        self.assertEqual(intervals, [(0.0, 1.0)])


class TestVectorLayer(unittest.TestCase):
    """This class contains methods for testing _vector_export.VectorLayer."""

    def setUp(self) -> None:
        """Set up the layer, the camera, and the Axes to draw onto.

        :return: None
        """
        self.layer = _vector_export.VectorLayer()
        self.camera = vector_export_fixtures.make_top_down_camera_fixture()
        figure = matplotlib.figure.Figure()
        self.axes = figure.add_axes((0.0, 0.0, 1.0, 1.0))

    def get_drawn_collection(
        self, collection_type: type[matplotlib.collections.Collection]
    ) -> matplotlib.collections.Collection:
        """Draws the layer and returns the one collection it drew of a type.

        :param collection_type: The type of the collection to return.
        :return: The collection.
        """
        self.layer.draw(self.axes, self.camera, 1.0)
        collections = [
            collection
            for collection in self.axes.collections
            if type(collection) is collection_type
        ]
        self.assertEqual(len(collections), 1)
        return collections[0]

    def test_outlines_a_shared_edge_once(self) -> None:
        """Test that neighboring quadrilaterals outline the edge they share once, and
        that it is owned by both of their triangles."""
        self.layer.add_quadrilaterals(
            np.array(
                [
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [1.0, 1.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ],
                    [
                        [1.0, 0.0, 0.0],
                        [2.0, 0.0, 0.0],
                        [2.0, 1.0, 0.0],
                        [1.0, 1.0, 0.0],
                    ],
                ],
                dtype=float,
            ),
            np.ones((2, 4), dtype=float),
            "black",
            1.0,
        )
        self.assertEqual(len(self.layer._listStrokes_D_Do), 7)
        self.assertIn({0, 1, 2, 3}, self.layer._stroke_owners)

    def test_closes_a_closed_polyline(self) -> None:
        """Test that a closed polyline gains a stroke from its last point back to its
        first, and an open one doesn't."""
        square_D_Do = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            dtype=float,
        )
        self.layer.add_polylines([square_D_Do], "black", 1.0)
        self.assertEqual(len(self.layer._listStrokes_D_Do), 3)
        self.layer.add_polylines([square_D_Do], "black", 1.0, closed=True)
        self.assertEqual(len(self.layer._listStrokes_D_Do), 7)
        npt.assert_array_equal(self.layer._listStrokes_D_Do[-1], square_D_Do[[3, 0]])

    def test_a_convex_occluder_hides_a_polyline_behind_it(self) -> None:
        """Test that a polyline passing behind a convex occluder is drawn as the two
        parts on either side of it.

        The cube's silhouette is at x = 90.0 and x = 110.0 pixels, and its 1.0 pixel
        wide outline hides the polyline for another pixel beyond it on each side, which
        is half the outline's width plus half the polyline's, so the polyline's round
        ends stop short of the outline rather than covering it.
        """
        self.layer.add_convex_occluder(
            vector_export_fixtures.make_cube_points_fixture(),
            vector_export_fixtures.make_cube_faces_fixture(),
            "white",
            "black",
            1.0,
        )
        self.layer.add_polylines(
            [np.array([[-3.0, 0.0, -5.0], [3.0, 0.0, -5.0]], dtype=float)],
            "black",
            1.0,
        )
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in self.get_drawn_collection(
                    matplotlib.collections.PathCollection
                ).get_paths()
            ]
        )
        self.assertEqual(len(segments), 2)
        npt.assert_allclose(segments[0], [[70.0, 50.0], [89.0, 50.0]])
        npt.assert_allclose(segments[1], [[111.0, 50.0], [130.0, 50.0]])

    def test_scales_the_line_widths_where_drawn_and_where_they_hide(self) -> None:
        """Test that a line width scale widens the strokes and the convex occluders'
        outlines both as drawn and where they hide the strokes behind them.

        With a scale of 2.0, the cube's outline and the polyline are each 2.0 pixels
        wide, so the outline hides the polyline for 2.0 pixels beyond the cube's
        silhouette at x = 90.0 and x = 110.0 pixels, which is half the outline's width
        plus half the polyline's.
        """
        self.layer.add_convex_occluder(
            vector_export_fixtures.make_cube_points_fixture(),
            vector_export_fixtures.make_cube_faces_fixture(),
            "white",
            "black",
            1.0,
        )
        self.layer.add_polylines(
            [np.array([[-3.0, 0.0, -5.0], [3.0, 0.0, -5.0]], dtype=float)],
            "black",
            1.0,
        )
        self.layer.draw(self.axes, self.camera, 1.0, line_width_scale=2.0)
        stroke_collection = next(
            collection
            for collection in self.axes.collections
            if type(collection) is matplotlib.collections.PathCollection
        )
        polygon_collection = next(
            collection
            for collection in self.axes.collections
            if type(collection) is matplotlib.collections.PolyCollection
        )
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in stroke_collection.get_paths()
            ]
        )
        self.assertEqual(len(segments), 2)
        npt.assert_allclose(segments[0], [[70.0, 50.0], [88.0, 50.0]])
        npt.assert_allclose(segments[1], [[112.0, 50.0], [130.0, 50.0]])
        npt.assert_allclose(
            np.array(stroke_collection.get_linewidth(), dtype=float),
            np.full(1, 2.0 * _vector_export.POINTS_PER_PIXEL, dtype=float),
        )
        npt.assert_allclose(
            np.array(polygon_collection.get_linewidth(), dtype=float),
            np.full(1, 2.0 * _vector_export.POINTS_PER_PIXEL, dtype=float),
        )

    def test_a_convex_occluder_leaves_a_polyline_in_front_of_it(self) -> None:
        """Test that a polyline passing in front of a convex occluder is drawn whole."""
        self.layer.add_convex_occluder(
            vector_export_fixtures.make_cube_points_fixture(),
            vector_export_fixtures.make_cube_faces_fixture(),
            "white",
            "black",
            1.0,
        )
        self.layer.add_polylines(
            [np.array([[-3.0, 0.0, 5.0], [3.0, 0.0, 5.0]], dtype=float)],
            "black",
            1.0,
        )
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in self.get_drawn_collection(
                    matplotlib.collections.PathCollection
                ).get_paths()
            ]
        )
        self.assertEqual(len(segments), 1)
        npt.assert_allclose(segments[0], [[70.0, 50.0], [130.0, 50.0]])

    def test_outlines_an_outlined_face_only_while_it_faces_the_camera(self) -> None:
        """Test that a convex occluder is outlined along its silhouette, plus the
        boundary of each outlined face that faces the camera, but not of one that faces
        away."""
        for outlined_face_id, num_paths in ((1, 2), (0, 1)):
            with self.subTest(outlined_face_id=outlined_face_id):
                self.setUp()
                self.layer.add_convex_occluder(
                    vector_export_fixtures.make_cube_points_fixture(),
                    vector_export_fixtures.make_cube_faces_fixture(),
                    "white",
                    "black",
                    1.0,
                    outlined_face_ids=[outlined_face_id],
                )
                paths = self.get_drawn_collection(
                    matplotlib.collections.PolyCollection
                ).get_paths()
                self.assertEqual(len(paths), num_paths)

    def test_returns_the_zorder_after_its_passes(self) -> None:
        """Test that drawing a layer returns the zorder after its five passes."""
        self.assertEqual(self.layer.draw(self.axes, self.camera, 3.0), 8.0)

    def test_paints_the_later_of_two_coincident_strokes_on_top(self) -> None:
        """Test that where two strokes lie on top of each other, neither hides the
        other, and the one added later is painted last, so it shows, as it does in
        VTK."""
        stroke_D_Do = np.array([[-3.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=float)
        self.layer.add_polylines([stroke_D_Do], "red", 1.0)
        self.layer.add_polylines([stroke_D_Do], "blue", 1.0)
        stroke_collection = self.get_drawn_collection(
            matplotlib.collections.PathCollection
        )
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in stroke_collection.get_paths()
            ]
        )
        self.assertEqual(len(segments), 2)
        npt.assert_array_equal(
            stroke_collection.get_edgecolor()[-1], matplotlib.colors.to_rgba("blue")
        )

    def test_a_nearer_stroke_hides_a_crossing_stroke_behind_it(self) -> None:
        """Test that where two strokes cross on screen, the nearer one hides the farther
        one across its width.

        The farther stroke runs along the x axis and the nearer one along the y axis,
        1.0 meter closer to the camera, crossing it at x = 100.0 pixels. Both are 1.0
        pixel wide, so the farther stroke is hidden for 1.0 pixel to either side, which
        is half the nearer stroke's width plus half its own.
        """
        self.layer.add_polylines(
            [np.array([[-3.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=float)], "red", 1.0
        )
        self.layer.add_polylines(
            [np.array([[0.0, -2.0, 1.0], [0.0, 2.0, 1.0]], dtype=float)], "blue", 1.0
        )
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in self.get_drawn_collection(
                    matplotlib.collections.PathCollection
                ).get_paths()
            ]
        )
        self.assertEqual(len(segments), 3)
        npt.assert_allclose(segments[0], [[70.0, 50.0], [99.0, 50.0]])
        npt.assert_allclose(segments[1], [[101.0, 50.0], [130.0, 50.0]])
        npt.assert_allclose(segments[2], [[100.0, 30.0], [100.0, 70.0]])

    def test_equally_deep_crossing_strokes_hide_nothing(self) -> None:
        """Test that two strokes crossing on screen at the same depth are both drawn
        whole."""
        self.layer.add_polylines(
            [np.array([[-3.0, 0.0, 0.0], [3.0, 0.0, 0.0]], dtype=float)], "red", 1.0
        )
        self.layer.add_polylines(
            [np.array([[0.0, -2.0, 0.0], [0.0, 2.0, 0.0]], dtype=float)], "blue", 1.0
        )
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in self.get_drawn_collection(
                    matplotlib.collections.PathCollection
                ).get_paths()
            ]
        )
        self.assertEqual(len(segments), 2)
        npt.assert_allclose(segments[0], [[70.0, 50.0], [130.0, 50.0]])
        npt.assert_allclose(segments[1], [[100.0, 30.0], [100.0, 70.0]])

    def test_strokes_that_share_an_end_never_hide_each_other(self) -> None:
        """Test that a polyline's segments don't hide each other where they meet.

        The polyline turns sharply back on itself at its middle point, so its second
        segment starts within the first segment's width on screen and runs away from the
        camera, behind it. Since the two share an end, neither hides the other, and both
        are drawn whole.
        """
        polyline_D_Do = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.05, -1.0]], dtype=float
        )
        self.layer.add_polylines([polyline_D_Do], "black", 1.0)
        segments = np.concatenate(
            [
                np.reshape(path.vertices, (-1, 2, 2))
                for path in self.get_drawn_collection(
                    matplotlib.collections.PathCollection
                ).get_paths()
            ]
        )
        self.assertEqual(len(segments), 2)
        npt.assert_allclose(segments[0], [[100.0, 50.0], [110.0, 50.0]])
        npt.assert_allclose(segments[1], [[110.0, 50.0], [105.0, 50.5]])


class TestVectorSceneSave(unittest.TestCase):
    """This class contains methods for testing _vector_export.VectorScene.save."""

    def setUp(self) -> None:
        """Set up a scene with a filled quadrilateral and a run of text, and a temporary
        directory to save it into.

        :return: None
        """
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temporary_path = Path(self.temporary_directory.name)
        self.camera = vector_export_fixtures.make_top_down_camera_fixture()

        self.scene = _vector_export.VectorScene()
        self.scene.add_layer().add_quadrilaterals(
            np.array(
                [
                    [
                        [-1.0, -1.0, 0.0],
                        [1.0, -1.0, 0.0],
                        [1.0, 1.0, 0.0],
                        [-1.0, 1.0, 0.0],
                    ]
                ],
                dtype=float,
            ),
            np.array([[0.5, 1.0, 0.0, 1.0]], dtype=float),
            "black",
            1.0,
        )
        self.scene.add_text(
            _vector_export.VectorText(
                x=100.0,
                y=10.0,
                text="Lift Coefficient",
                font_properties=matplotlib.font_manager.FontProperties(
                    family=_fonts.FONT_FAMILY, fname=_fonts.FONT_PATH
                ),
                font_size=15.0,
                color="black",
                horizontal_alignment="center",
                vertical_alignment="baseline",
            )
        )

    def tearDown(self) -> None:
        """Remove the temporary directory and any files the test wrote into it.

        :return: None
        """
        self.temporary_directory.cleanup()

    def test_saves_an_svg(self) -> None:
        """Test that a path ending with ".svg" saves an svg file sized to the window."""
        path = self.temporary_path / "scene.svg"
        self.scene.save(path, self.camera, None)
        contents = path.read_text()
        self.assertIn("<svg", contents)
        self.assertIn('width="144pt" height="72pt"', contents)

    def test_saves_an_svg_whose_text_carries_its_font(self) -> None:
        """Test that a saved svg writes its text as text and embeds the vendored font
        that the text names."""
        path = self.temporary_path / "scene.svg"
        self.scene.save(path, self.camera, None)
        contents = path.read_text(encoding="utf-8")
        self.assertEqual(contents.count("@font-face"), 1)
        self.assertIn(f"font-family: '{_fonts.FONT_FAMILY}'", contents)
        self.assertIn(">Lift Coefficient</text>", contents)

    def test_saves_a_pdf_with_truetype_fonts(self) -> None:
        """Test that a path ending with ".pdf" saves a pdf file whose fonts are embedded
        as TrueType."""
        path = self.temporary_path / "scene.pdf"
        self.scene.save(path, self.camera, "white")
        contents = path.read_bytes()
        self.assertTrue(contents.startswith(b"%PDF"))
        self.assertIn(b"/FontFile2", contents)
        self.assertNotIn(b"/Type3", contents)

    def test_fills_a_texts_background_box_behind_it(self) -> None:
        """Test that a text's background box is filled in its color, before the text is
        drawn over it."""
        self.scene.add_text(
            _vector_export.VectorText(
                x=50.0,
                y=50.0,
                text="Boxed",
                font_properties=matplotlib.font_manager.FontProperties(
                    family=_fonts.FONT_FAMILY, fname=_fonts.FONT_PATH
                ),
                font_size=15.0,
                color="black",
                horizontal_alignment="center",
                vertical_alignment="center",
                background_box=(30.0, 70.0, 40.0, 60.0),
                background_color="red",
            )
        )
        path = self.temporary_path / "scene.svg"
        self.scene.save(path, self.camera, None)
        contents = path.read_text(encoding="utf-8")
        self.assertIn("fill: #ff0000", contents)
        self.assertLess(contents.index("fill: #ff0000"), contents.index(">Boxed<"))

    def test_writes_math_in_the_texts_math_font_family(self) -> None:
        """Test that a text's math is written in its math font family's fonts."""
        self.scene.add_text(
            _vector_export.VectorText(
                x=50.0,
                y=50.0,
                text=r"$\hat{x}$",
                font_properties=matplotlib.font_manager.FontProperties(
                    family=_fonts.MONO_FONT_FAMILY, fname=_fonts.MONO_FONT_PATH
                ),
                font_size=15.0,
                color="black",
                horizontal_alignment="center",
                vertical_alignment="center",
                math_font_family="stix",
            )
        )
        path = self.temporary_path / "scene.svg"
        self.scene.save(path, self.camera, None)
        self.assertIn("font-family: 'STIXGeneral'", path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
