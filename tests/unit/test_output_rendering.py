"""This module contains classes to test the output rendering functions.

The functions that drive a Plotter or a render window are covered by the integration
tests instead, as are the wake ring vortex surfaces, which are built from a history that
only a solved simulation carries. The classes here cover the computation and the
geometry building that feed them, which are settled before any rendering begins.
"""

import unittest

import matplotlib.colors
import numpy as np
import numpy.testing as npt
import pyvista as pv

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _colormaps, _output_rendering, _transformations
from tests.unit.fixtures import operating_point_fixtures, output_rendering_fixtures


class TestGetWindowScale(unittest.TestCase):
    """This class contains methods for testing _output_rendering.get_window_scale."""

    def test_is_one_at_the_reference_window_size(self) -> None:
        """Test that the reference window size scales the font sizes by exactly 1.0.

        Every font size and line width in the visualizations is tuned against that
        window, so it is the one size that must leave them alone.
        """
        self.assertEqual(
            _output_rendering.get_window_scale(
                _output_rendering.REFERENCE_WINDOW_SIZE[0],
                _output_rendering.REFERENCE_WINDOW_SIZE[1],
            ),
            1.0,
        )

    def test_grows_with_a_larger_window(self) -> None:
        """Test that a window twice the reference size doubles the scale."""
        self.assertEqual(_output_rendering.get_window_scale(2048, 1536), 2.0)

    def test_takes_the_height_ratio_for_a_wide_window(self) -> None:
        """Test that a wide but short window is scaled by its height ratio."""
        self.assertEqual(_output_rendering.get_window_scale(4096, 384), 0.5)

    def test_takes_the_width_ratio_for_a_tall_window(self) -> None:
        """Test that a tall but narrow window is scaled by its width ratio.

        The scalar labels are anchored near the right edge and grow rightward, so the
        horizontal room is what runs out first in such a window.
        """
        self.assertEqual(_output_rendering.get_window_scale(512, 3072), 0.5)


class TestResolvePlayback(unittest.TestCase):
    """This class contains methods for testing _output_rendering.resolve_playback."""

    playback_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )
    long_step_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up the shared test fixtures."""
        cls.playback_solver = output_rendering_fixtures.make_playback_solver_fixture()
        cls.long_step_solver = (
            output_rendering_fixtures.make_long_step_playback_solver_fixture()
        )

    def test_default_speed_keeps_every_frame(self) -> None:
        """Test that the default speed is the fastest one that drops no frames.

        A time step of 0.01 seconds carries only 0.5 seconds of simulation per second of
        playback at the maximum frame rate, so the default plays at half speed rather
        than dropping frames to reach true speed.
        """
        playback = _output_rendering.resolve_playback(self.playback_solver, None, True)
        self.assertEqual(playback.keep_every, 1)
        self.assertEqual(playback.frame_rate, 50.0)

    def test_default_speed_is_true_speed_when_the_frame_rate_allows(self) -> None:
        """Test that the default speed is true speed when the frame rate can carry it.

        A time step of 0.05 seconds needs only 20 frames per second of playback to run
        at true speed, which is within the maximum, so the default asks for no less.
        """
        playback = _output_rendering.resolve_playback(self.long_step_solver, None, True)
        self.assertEqual(playback.keep_every, 1)
        self.assertEqual(playback.frame_rate, 20.0)
        self.assertEqual(playback.overlay_texts[0][0], "Speed: 100.0%")

    def test_a_slower_speed_lowers_the_frame_rate(self) -> None:
        """Test that a speed within the frame rate's reach drops no frames."""
        playback = _output_rendering.resolve_playback(self.playback_solver, 0.25, True)
        self.assertEqual(playback.keep_every, 1)
        self.assertEqual(playback.frame_rate, 25.0)

    def test_a_faster_speed_drops_frames(self) -> None:
        """Test that a speed beyond the frame rate's reach is met by dropping frames.

        True speed needs 100 frames per second of playback here, which is twice the
        maximum, so every other time step is dropped and the rest play at the maximum.
        """
        playback = _output_rendering.resolve_playback(self.playback_solver, 1.0, True)
        self.assertEqual(playback.keep_every, 2)
        self.assertEqual(playback.frame_rate, 50.0)

    def test_the_speed_overlay_reports_the_achieved_speed(self) -> None:
        """Test that the speed overlay reports the speed the animation reaches."""
        playback = _output_rendering.resolve_playback(self.playback_solver, 0.25, True)
        self.assertEqual(
            playback.overlay_texts[0],
            ("Speed: 25.00%", _output_rendering._TEXT_SPEED_POSITION),
        )

    def test_only_the_speed_overlay_appears_when_no_frames_are_dropped(self) -> None:
        """Test that an animation that keeps every frame carries one overlay."""
        playback = _output_rendering.resolve_playback(self.playback_solver, 0.25, True)
        self.assertEqual(len(playback.overlay_texts), 1)

    def test_a_second_overlay_reports_the_dropped_frames(self) -> None:
        """Test that an animation that drops frames says so in a second overlay."""
        playback = _output_rendering.resolve_playback(self.playback_solver, 1.0, True)
        self.assertEqual(
            playback.overlay_texts[1],
            ("Frames: Every 1 of 2", _output_rendering._TEXT_DROPPED_FRAMES_POSITION),
        )

    def test_rejects_a_speed_that_saves_less_than_one_frame_per_second(self) -> None:
        """Test that a speed too slow to fill a second of playback is rejected."""
        with self.assertRaises(ValueError) as context:
            _output_rendering.resolve_playback(self.playback_solver, 0.005, True)
        self.assertIn("too slow", str(context.exception))

    def test_rejects_a_speed_that_saves_fewer_than_two_frames(self) -> None:
        """Test that a speed too fast to save a second frame is rejected.

        The error names the fastest speed the simulation can be animated at, which is
        the speed that saves its first and last time steps and nothing between them.
        """
        with self.assertRaises(ValueError) as context:
            _output_rendering.resolve_playback(self.playback_solver, 6.0, True)
        self.assertIn("too fast", str(context.exception))
        self.assertIn("5.0", str(context.exception))


class TestResolvePlaybackAliasingWarning(unittest.TestCase):
    """This class contains methods for testing resolve_playback's aliasing warning."""

    fast_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )
    static_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up the shared test fixtures."""
        cls.fast_solver = (
            output_rendering_fixtures.make_fast_motion_playback_solver_fixture()
        )
        cls.static_solver = (
            output_rendering_fixtures.make_static_playback_solver_fixture()
        )

    def test_warns_when_the_saved_frames_stop_resolving_the_motion(self) -> None:
        """Test that undersampling the fastest prescribed motion is warned about.

        Dropping every other frame leaves 5 frames per cycle of a 0.1 second period,
        which is below the floor the warning is held to.
        """
        with self.assertLogs("pterasoftware.output", level="WARNING") as context:
            _output_rendering.resolve_playback(self.fast_solver, 1.0, True)
        self.assertIn("frames per cycle", context.output[0])

    def test_does_not_warn_when_no_frames_are_dropped(self) -> None:
        """Test that an animation that keeps every frame is not warned about."""
        with self.assertNoLogs("pterasoftware.output", level="WARNING"):
            _output_rendering.resolve_playback(self.fast_solver, 0.25, True)

    def test_does_not_warn_when_the_animation_is_not_saved(self) -> None:
        """Test that an animation that is only shown is not warned about.

        Frames are dropped from the saved file alone, so an animation that saves nothing
        undersamples nothing.
        """
        with self.assertNoLogs("pterasoftware.output", level="WARNING"):
            _output_rendering.resolve_playback(self.fast_solver, 1.0, False)

    def test_does_not_warn_for_a_static_geometry(self) -> None:
        """Test that a simulation whose geometry never moves is not warned about."""
        with self.assertNoLogs("pterasoftware.output", level="WARNING"):
            _output_rendering.resolve_playback(self.static_solver, 1.0, True)


class TestGetFreeFlightTransformation(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_free_flight_transformation."""

    def test_maps_the_cg_onto_its_earth_position(self) -> None:
        """Test that the first Airplane's CG lands at its position in Earth axes.

        The CG is the origin of the axes the transformation maps out of, so wherever the
        first Airplane has flown to, the CG maps onto that position.
        """
        operating_point = (
            operating_point_fixtures.make_with_cg_position_operating_point_fixture()
        )
        T_pas = _output_rendering.get_free_flight_transformation(operating_point)
        Cg_E_Eo = _transformations.apply_T_to_vectors(
            T_pas, np.zeros(3, dtype=float), is_position=True
        )
        npt.assert_allclose(Cg_E_Eo, operating_point.CgP1_E_Eo, atol=1e-12)

    def test_reduces_to_the_rotation_at_the_earth_origin(self) -> None:
        """Test that a first Airplane at the Earth origin needs only the rotation.

        With nothing to translate, the transformation is the geometry axes to Earth axes
        rotation it is built on.
        """
        operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        npt.assert_allclose(
            _output_rendering.get_free_flight_transformation(operating_point),
            operating_point.T_pas_GP1_CgP1_to_E_CgP1,
            atol=1e-12,
        )

    def test_leaves_directions_untranslated(self) -> None:
        """Test that the translation reaches positions alone.

        A direction has no reference point, so mapping one through this transformation
        must give what the geometry axes to Earth axes rotation alone gives.
        """
        operating_point = (
            operating_point_fixtures.make_with_cg_position_operating_point_fixture()
        )
        direction_GP1 = np.array([1.0, 2.0, 3.0], dtype=float)
        npt.assert_allclose(
            _transformations.apply_T_to_vectors(
                _output_rendering.get_free_flight_transformation(operating_point),
                direction_GP1,
                is_position=False,
            ),
            _transformations.apply_T_to_vectors(
                operating_point.T_pas_GP1_CgP1_to_E_CgP1,
                direction_GP1,
                is_position=False,
            ),
            atol=1e-12,
        )


class TestGetPanelSurfaces(unittest.TestCase):
    """This class contains methods for testing _output_rendering.get_panel_surfaces."""

    def test_builds_one_quadrilateral_per_panel(self) -> None:
        """Test that every Panel becomes one four sided cell."""
        airplanes = output_rendering_fixtures.make_placed_airplanes_fixture()
        assert airplanes[0].wings[0].panels is not None
        panels = np.ravel(airplanes[0].wings[0].panels)
        surfaces = _output_rendering.get_panel_surfaces(airplanes)
        self.assertEqual(surfaces.n_cells, len(panels))
        self.assertEqual(surfaces.n_points, 4 * len(panels))

    def test_walks_a_panels_corners_in_order(self) -> None:
        """Test that a cell's points are its Panel's corners, front left first.

        The corners are wound front left, front right, back right, back left, so that
        the cell traces the Panel's outline rather than crossing it.
        """
        airplanes = output_rendering_fixtures.make_placed_airplanes_fixture()
        assert airplanes[0].wings[0].panels is not None
        panel = np.ravel(airplanes[0].wings[0].panels)[0]
        surfaces = _output_rendering.get_panel_surfaces(airplanes)
        npt.assert_allclose(
            surfaces.points[:4],
            [
                panel.Flpp_GP1_CgP1,
                panel.Frpp_GP1_CgP1,
                panel.Brpp_GP1_CgP1,
                panel.Blpp_GP1_CgP1,
            ],
        )

    def test_gathers_every_airplanes_panels(self) -> None:
        """Test that a formation's Airplanes all reach the same mesh."""
        airplanes = output_rendering_fixtures.make_formation_airplanes_fixture()
        first_surfaces = _output_rendering.get_panel_surfaces((airplanes[0],))
        follower_surfaces = _output_rendering.get_panel_surfaces((airplanes[1],))
        both_surfaces = _output_rendering.get_panel_surfaces(airplanes)
        self.assertEqual(
            both_surfaces.n_cells, first_surfaces.n_cells + follower_surfaces.n_cells
        )


class TestGetStreamlineSurfaces(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_streamline_surfaces."""

    def test_builds_one_polyline_per_streamline(self) -> None:
        """Test that each streamline becomes one cell holding all of its points.

        The whole set is one mesh rather than one mesh per segment, which is what keeps
        VTK from being given an actor per segment.
        """
        grid = output_rendering_fixtures.make_streamline_points_fixture()
        num_points, num_streamlines = grid.shape[:2]
        surfaces = _output_rendering.get_streamline_surfaces(grid)
        self.assertEqual(surfaces.n_cells, num_streamlines)
        self.assertEqual(surfaces.n_points, num_points * num_streamlines)

    def test_runs_each_polyline_along_its_own_streamline(self) -> None:
        """Test that a polyline's points are its streamline's, in order."""
        grid = output_rendering_fixtures.make_streamline_points_fixture()
        num_points, num_streamlines = grid.shape[:2]
        surfaces = _output_rendering.get_streamline_surfaces(grid)
        for streamline_num in range(num_streamlines):
            first_point = streamline_num * num_points
            npt.assert_allclose(
                surfaces.points[first_point : first_point + num_points],
                grid[:, streamline_num, :],
            )

    def test_describes_each_polyline_as_a_count_then_its_indices(self) -> None:
        """Test that the cell array is in the flat format VTK expects."""
        grid = output_rendering_fixtures.make_streamline_points_fixture()
        surfaces = _output_rendering.get_streamline_surfaces(grid)
        npt.assert_array_equal(
            surfaces.lines, [4, 0, 1, 2, 3, 4, 4, 5, 6, 7, 4, 8, 9, 10, 11]
        )


class TestGetImageSurfaceMeshAndTexture(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_image_surface_mesh_and_texture."""

    def test_returns_nothing_without_an_image_surface(self) -> None:
        """Test that an OperatingPoint with no image surface builds no plane."""
        operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        self.assertIsNone(
            _output_rendering.get_image_surface_mesh_and_texture(
                operating_point,
                output_rendering_fixtures.make_geometry_bounds_fixture(),
            )
        )

    def test_centers_the_plane_on_the_projected_bounding_box_center(self) -> None:
        """Test that the plane is centered where the geometry sits over the surface."""
        operating_point = (
            operating_point_fixtures.make_with_ground_surface_operating_point_fixture()
        )
        result = _output_rendering.get_image_surface_mesh_and_texture(
            operating_point, output_rendering_fixtures.make_geometry_bounds_fixture()
        )
        self.assertIsNotNone(result)
        assert result is not None
        mesh, _ = result
        npt.assert_allclose(mesh.center, [0.0, 0.0, -10.0], atol=1e-9)

    def test_lays_the_plane_on_the_image_surface(self) -> None:
        """Test that every one of the plane's points lies on the image surface."""
        operating_point = (
            operating_point_fixtures.make_with_ground_surface_operating_point_fixture()
        )
        result = _output_rendering.get_image_surface_mesh_and_texture(
            operating_point, output_rendering_fixtures.make_geometry_bounds_fixture()
        )
        self.assertIsNotNone(result)
        assert result is not None
        mesh, _ = result
        surface_normal = operating_point.surfaceNormal_GP1
        surface_point = operating_point.surfacePoint_GP1_CgP1
        assert surface_normal is not None
        assert surface_point is not None
        offsets = (np.array(mesh.points) - surface_point) @ surface_normal
        npt.assert_allclose(offsets, np.zeros(mesh.n_points), atol=1e-9)

    def test_sizes_the_plane_to_the_bounding_box_diagonal(self) -> None:
        """Test that the plane is a fixed multiple of the geometry's diagonal.

        Sizing it this way is what keeps it looking large next to the geometry however
        large the geometry is.
        """
        operating_point = (
            operating_point_fixtures.make_with_ground_surface_operating_point_fixture()
        )
        result = _output_rendering.get_image_surface_mesh_and_texture(
            operating_point, output_rendering_fixtures.make_geometry_bounds_fixture()
        )
        self.assertIsNotNone(result)
        assert result is not None
        mesh, _ = result
        expected_size = float(_output_rendering._IMAGE_SURFACE_SCALE * np.sqrt(56.0))
        self.assertAlmostEqual(
            mesh.bounds.x_max - mesh.bounds.x_min, expected_size, places=4
        )

    def test_builds_a_two_color_checkerboard_texture(self) -> None:
        """Test that the texture alternates between the two image surface colors."""
        operating_point = (
            operating_point_fixtures.make_with_ground_surface_operating_point_fixture()
        )
        result = _output_rendering.get_image_surface_mesh_and_texture(
            operating_point, output_rendering_fixtures.make_geometry_bounds_fixture()
        )
        self.assertIsNotNone(result)
        assert result is not None
        _, texture = result
        image = texture.to_array()
        checker_size = _output_rendering._IMAGE_SURFACE_CHECKER_SIZE
        self.assertEqual(image.shape, (checker_size, checker_size, 3))
        npt.assert_array_equal(image[0, 0], _output_rendering._IMAGE_SURFACE_COLOR_A)
        npt.assert_array_equal(image[0, 1], _output_rendering._IMAGE_SURFACE_COLOR_B)
        npt.assert_array_equal(image[1, 0], _output_rendering._IMAGE_SURFACE_COLOR_B)
        npt.assert_array_equal(image[1, 1], _output_rendering._IMAGE_SURFACE_COLOR_A)


class TestTransformMesh(unittest.TestCase):
    """This class contains methods for testing _output_rendering.transform_mesh."""

    def test_maps_the_points_through_the_transformation(self) -> None:
        """Test that the mesh's points are mapped as positions."""
        mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        T_pas = _transformations.generate_trans_T(
            translations=np.array([1.0, 2.0, 3.0], dtype=float), passive=True
        )
        transformed = _output_rendering.transform_mesh(mesh, T_pas)
        npt.assert_allclose(
            transformed.points,
            np.array(mesh.points) - np.array([1.0, 2.0, 3.0]),
            atol=1e-12,
        )

    def test_leaves_the_original_mesh_alone(self) -> None:
        """Test that the mesh handed in is not the one that comes back changed.

        An animation transforms the same source geometry once per frame, so a
        transformation that wrote back into it would compound across the frames.
        """
        mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        original_points = np.array(mesh.points)
        T_pas = _transformations.generate_trans_T(
            translations=np.array([1.0, 2.0, 3.0], dtype=float), passive=True
        )
        _output_rendering.transform_mesh(mesh, T_pas)
        npt.assert_array_equal(mesh.points, original_points)

    def test_keeps_the_faces(self) -> None:
        """Test that only the points move, leaving the mesh's topology intact."""
        mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        transformed = _output_rendering.transform_mesh(mesh, np.eye(4, dtype=float))
        npt.assert_array_equal(transformed.faces, mesh.faces)


class TestGetFreeFlightFitParallelScale(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_free_flight_fit_parallel_scale."""

    def test_the_margin_pads_the_scale(self) -> None:
        """Test that the scale is the geometry's half extent, padded by the margin."""
        mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        scale = _output_rendering.get_free_flight_fit_parallel_scale(
            [mesh],
            np.zeros(3, dtype=float),
            np.array([0.0, 0.0, -10.0], dtype=float),
            np.array([0.0, -1.0, 0.0], dtype=float),
            margin=2.0,
        )
        self.assertAlmostEqual(scale, 2.0, places=12)

    def test_ignores_the_extent_along_the_view_direction(self) -> None:
        """Test that depth does not enlarge a parallel projection.

        A parallel projection's scale is half its viewport's height, so only the extents
        across the camera's screen axes can widen it.
        """
        cube_scale = _output_rendering.get_free_flight_fit_parallel_scale(
            [output_rendering_fixtures.make_cube_mesh_fixture()],
            np.zeros(3, dtype=float),
            np.array([0.0, 0.0, -10.0], dtype=float),
            np.array([0.0, -1.0, 0.0], dtype=float),
        )
        deep_box_scale = _output_rendering.get_free_flight_fit_parallel_scale(
            [output_rendering_fixtures.make_deep_box_mesh_fixture()],
            np.zeros(3, dtype=float),
            np.array([0.0, 0.0, -10.0], dtype=float),
            np.array([0.0, -1.0, 0.0], dtype=float),
        )
        self.assertAlmostEqual(deep_box_scale, cube_scale, places=12)

    def test_measures_the_extents_from_the_focal_point(self) -> None:
        """Test that geometry off to one side of the focal point stays in view.

        The focal point projects to the center of the viewport, so a mesh that sits away
        from it needs a scale large enough to reach back across that gap.
        """
        mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        scale = _output_rendering.get_free_flight_fit_parallel_scale(
            [mesh],
            np.array([5.0, 0.0, 0.0], dtype=float),
            np.array([0.0, 0.0, -10.0], dtype=float),
            np.array([0.0, -1.0, 0.0], dtype=float),
            margin=1.0,
        )
        self.assertAlmostEqual(scale, 6.0, places=12)

    def test_frames_every_mesh_it_is_given(self) -> None:
        """Test that the scale covers the whole set of meshes rather than one."""
        near_mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        far_mesh = output_rendering_fixtures.make_cube_mesh_fixture()
        far_mesh.points = np.array(far_mesh.points) + np.array([10.0, 0.0, 0.0])
        scale = _output_rendering.get_free_flight_fit_parallel_scale(
            [near_mesh, far_mesh],
            np.zeros(3, dtype=float),
            np.array([0.0, 0.0, -10.0], dtype=float),
            np.array([0.0, -1.0, 0.0], dtype=float),
            margin=1.0,
        )
        self.assertAlmostEqual(scale, 11.0, places=12)


class TestMuteColor(unittest.TestCase):
    """This class contains methods for testing _output_rendering.mute_color."""

    def test_a_factor_of_zero_leaves_the_color_alone(self) -> None:
        """Test that muting a color by nothing returns the color."""
        npt.assert_allclose(_output_rendering.mute_color("red", 0.0), (1.0, 0.0, 0.0))

    def test_a_factor_of_one_returns_middle_gray(self) -> None:
        """Test that muting a color fully returns middle gray."""
        npt.assert_allclose(_output_rendering.mute_color("red", 1.0), (0.5, 0.5, 0.5))

    def test_accepts_an_rgb_tuple(self) -> None:
        """Test that a color given as components is muted halfway to middle gray.

        The muting is a linear interpolation, and a color given as components is parsed
        the same way a name is.
        """
        npt.assert_allclose(
            _output_rendering.mute_color((1.0, 0.0, 0.0), 0.5), (0.75, 0.25, 0.25)
        )

    def test_returns_python_floats(self) -> None:
        """Test that the muted color is returned as three Python floats."""
        muted = _output_rendering.mute_color("red", 0.5)
        self.assertEqual(len(muted), 3)
        for component in muted:
            self.assertIsInstance(component, float)


class TestMuteColormap(unittest.TestCase):
    """This class contains methods for testing _output_rendering.mute_colormap."""

    def test_a_factor_of_zero_leaves_the_colors_alone(self) -> None:
        """Test that muting a color map by nothing returns its colors."""
        muted = _output_rendering.mute_colormap(_colormaps.sequential_color_map, 0.0)
        npt.assert_allclose(muted(0.0), _colormaps.sequential_color_map(0.0))
        npt.assert_allclose(muted(1.0), _colormaps.sequential_color_map(1.0))

    def test_a_factor_of_one_returns_middle_gray(self) -> None:
        """Test that muting a color map fully leaves every color middle gray."""
        muted = _output_rendering.mute_colormap(_colormaps.sequential_color_map, 1.0)
        npt.assert_allclose(muted(0.0)[:3], (0.5, 0.5, 0.5))
        npt.assert_allclose(muted(0.5)[:3], (0.5, 0.5, 0.5))
        npt.assert_allclose(muted(1.0)[:3], (0.5, 0.5, 0.5))

    def test_returns_a_listed_color_map_of_256_colors(self) -> None:
        """Test that the muted color map is a ListedColormap sampled at 256 colors."""
        muted = _output_rendering.mute_colormap(_colormaps.sequential_color_map, 0.5)
        self.assertIsInstance(muted, matplotlib.colors.ListedColormap)
        self.assertEqual(muted.N, 256)

    def test_leaves_the_alpha_channel_alone(self) -> None:
        """Test that muting reaches the colors rather than their opacity."""
        muted = _output_rendering.mute_colormap(_colormaps.sequential_color_map, 1.0)
        self.assertEqual(muted(0.0)[3], 1.0)


class TestGetScalars(unittest.TestCase):
    """This class contains methods for testing _output_rendering.get_scalars."""

    def test_induced_drag_negates_the_wind_axes_x_force(self) -> None:
        """Test that induced drag reads as positive against the freestream.

        Wind axes x points opposite the drag a Panel produces, so the coefficient
        carries the negated force.
        """
        airplanes = output_rendering_fixtures.make_loaded_airplanes_fixture()
        assert airplanes[0].wings[0].panels is not None
        panels = np.ravel(airplanes[0].wings[0].panels)
        scalars = _output_rendering.get_scalars(airplanes, "induced drag", 2.0)
        expected = [-panel.forces_W[0] / 2.0 / panel.area for panel in panels]
        npt.assert_allclose(scalars, expected)

    def test_side_force_keeps_the_wind_axes_y_force(self) -> None:
        """Test that side force reads in the same direction as wind axes y."""
        airplanes = output_rendering_fixtures.make_loaded_airplanes_fixture()
        assert airplanes[0].wings[0].panels is not None
        panels = np.ravel(airplanes[0].wings[0].panels)
        scalars = _output_rendering.get_scalars(airplanes, "side force", 2.0)
        expected = [panel.forces_W[1] / 2.0 / panel.area for panel in panels]
        npt.assert_allclose(scalars, expected)

    def test_lift_negates_the_wind_axes_z_force(self) -> None:
        """Test that lift reads as positive upward.

        Wind axes z points down, so the coefficient carries the negated force.
        """
        airplanes = output_rendering_fixtures.make_loaded_airplanes_fixture()
        assert airplanes[0].wings[0].panels is not None
        panels = np.ravel(airplanes[0].wings[0].panels)
        scalars = _output_rendering.get_scalars(airplanes, "lift", 2.0)
        expected = [-panel.forces_W[2] / 2.0 / panel.area for panel in panels]
        npt.assert_allclose(scalars, expected)

    def test_an_unrecognized_scalar_type_contributes_nothing(self) -> None:
        """Test that a scalar type outside the three named ones yields no scalars.

        The public output functions reject such a type before reaching here, so this is
        what an internal caller that skipped that check would see.
        """
        airplanes = output_rendering_fixtures.make_loaded_airplanes_fixture()
        self.assertEqual(
            _output_rendering.get_scalars(airplanes, "not a load", 2.0).shape, (0,)
        )


class TestChooseColorMap(unittest.TestCase):
    """This class contains methods for testing _output_rendering.choose_color_map."""

    def test_single_signed_scalars_get_the_sequential_color_map(self) -> None:
        """Test that scalars that keep one sign run in a single direction."""
        scalars = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
        color_map, _, _ = _output_rendering.choose_color_map(scalars)
        self.assertIs(color_map, _colormaps.sequential_color_map)

    def test_sign_changing_scalars_get_the_diverging_color_map(self) -> None:
        """Test that scalars that change sign are colored about their midpoint."""
        scalars = np.array([-2.0, -1.0, 1.0, 2.0], dtype=float)
        color_map, _, _ = _output_rendering.choose_color_map(scalars)
        self.assertIs(color_map, _colormaps.diverging_color_map)

    def test_the_diverging_limits_sit_symmetrically_about_zero(self) -> None:
        """Test that the diverging map's midpoint marks where the scalar changes
        sign."""
        scalars = np.array([-2.0, -1.0, 1.0, 2.0], dtype=float)
        _, c_min, c_max = _output_rendering.choose_color_map(scalars)
        self.assertAlmostEqual(c_min, -c_max, places=12)
        self.assertAlmostEqual(
            c_max,
            _output_rendering._COLOR_MAP_NUM_SIG * float(np.std(scalars)),
            places=12,
        )

    def test_a_tight_distribution_keeps_its_own_limits(self) -> None:
        """Test that limits never reach beyond the scalars they color.

        A distribution whose whole range sits within the sigma bound is colored across
        exactly that range.
        """
        scalars = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
        _, c_min, c_max = _output_rendering.choose_color_map(scalars)
        self.assertEqual(c_min, 1.0)
        self.assertEqual(c_max, 4.0)

    def test_an_outlier_cannot_flatten_the_contrast(self) -> None:
        """Test that a far outlier is held outside the limits.

        Coloring across the whole range would leave every other Panel in the bottom of
        the color map, so the upper limit is held within a fixed number of standard
        deviations of the mean.
        """
        scalars = output_rendering_fixtures.make_outlier_scalars_fixture()
        _, c_min, c_max = _output_rendering.choose_color_map(scalars)
        self.assertLess(c_max, float(np.max(scalars)))
        self.assertAlmostEqual(
            c_max,
            float(np.mean(scalars))
            + _output_rendering._COLOR_MAP_NUM_SIG * float(np.std(scalars)),
            places=12,
        )
        self.assertEqual(c_min, float(np.min(scalars)))


class TestGetAnimationImageSurface(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_animation_image_surface."""

    image_surface_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )
    plain_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up the shared test fixtures."""
        cls.image_surface_solver = (
            output_rendering_fixtures.make_image_surface_solver_fixture()
        )
        cls.plain_solver = output_rendering_fixtures.make_playback_solver_fixture()

    def test_builds_nothing_without_an_image_surface(self) -> None:
        """Test that a simulation with no image surface builds no plane."""
        step_airplanes = [
            steady_problem.airplanes
            for steady_problem in self.plain_solver.steady_problems
        ]
        self.assertEqual(
            _output_rendering.get_animation_image_surface(
                self.plain_solver, step_airplanes, [], False, False
            ),
            (None, None, None, None),
        )

    def test_builds_the_plane_from_the_last_time_step(self) -> None:
        """Test that a simulation with an image surface builds all four pieces.

        The plane is built from the last time step, whose wake is the most developed, so
        that one plane is large enough for every frame.
        """
        step_airplanes = [
            steady_problem.airplanes
            for steady_problem in self.image_surface_solver.steady_problems
        ]
        mesh, texture, T_reflect, bounds = (
            _output_rendering.get_animation_image_surface(
                self.image_surface_solver, step_airplanes, [], False, False
            )
        )
        self.assertIsInstance(mesh, pv.PolyData)
        self.assertIsInstance(texture, pv.Texture)
        self.assertIsNotNone(T_reflect)
        self.assertIsNotNone(bounds)

    def test_reflects_across_the_last_time_steps_image_surface(self) -> None:
        """Test that the reflection is the last time step's."""
        step_airplanes = [
            steady_problem.airplanes
            for steady_problem in self.image_surface_solver.steady_problems
        ]
        _, _, T_reflect, _ = _output_rendering.get_animation_image_surface(
            self.image_surface_solver, step_airplanes, [], False, False
        )
        last_operating_point = self.image_surface_solver.steady_problems[
            -1
        ].operating_point
        expected_T_reflect = last_operating_point.surfaceReflect_T_act_GP1_CgP1
        assert T_reflect is not None
        assert expected_T_reflect is not None
        npt.assert_allclose(T_reflect, expected_T_reflect)

    def test_the_bounding_box_spans_the_geometry_and_its_reflection(self) -> None:
        """Test that the box a camera fits to holds both copies of the geometry.

        The box leaves out the much larger plane, since fitting to that would leave the
        geometry too small to see.
        """
        step_airplanes = [
            steady_problem.airplanes
            for steady_problem in self.image_surface_solver.steady_problems
        ]
        _, _, T_reflect, bounds = _output_rendering.get_animation_image_surface(
            self.image_surface_solver, step_airplanes, [], False, False
        )
        panel_surfaces = _output_rendering.get_panel_surfaces(step_airplanes[-1])
        assert T_reflect is not None
        assert bounds is not None
        reflected_surfaces = _output_rendering.transform_mesh(panel_surfaces, T_reflect)
        self.assertAlmostEqual(
            bounds[4], min(panel_surfaces.bounds[4], reflected_surfaces.bounds[4])
        )
        self.assertAlmostEqual(
            bounds[5], max(panel_surfaces.bounds[5], reflected_surfaces.bounds[5])
        )

    def test_free_flight_maps_the_plane_into_earth_axes(self) -> None:
        """Test that free flight moves the plane by the last time step's transformation.

        Free flight renders its geometry in Earth axes and frames its camera to the
        whole trajectory, so it takes a plane in those axes and no bounding box.
        """
        step_airplanes = [
            steady_problem.airplanes
            for steady_problem in self.image_surface_solver.steady_problems
        ]
        T_pas = _transformations.generate_trans_T(
            translations=np.array([1.0, 2.0, 3.0], dtype=float), passive=True
        )
        step_transforms = [T_pas for _ in step_airplanes]
        mesh, _, _, bounds = _output_rendering.get_animation_image_surface(
            self.image_surface_solver, step_airplanes, step_transforms, True, False
        )
        body_fixed_mesh, _, _, _ = _output_rendering.get_animation_image_surface(
            self.image_surface_solver, step_airplanes, [], False, False
        )
        self.assertIsNone(bounds)
        assert mesh is not None
        assert body_fixed_mesh is not None

        # The two planes are sized differently, since free flight sizes its plane to the
        # geometry alone while a body-fixed drawing sizes it to the geometry and its
        # reflection together. PyVista stores a plane's points as 32 bit floats, so two
        # planes this far apart in size round their centers a little differently.
        npt.assert_allclose(
            mesh.center,
            np.array(body_fixed_mesh.center) - np.array([1.0, 2.0, 3.0]),
            atol=1e-5,
        )
