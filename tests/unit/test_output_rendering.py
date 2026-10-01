"""This module contains classes to test the output rendering functions.

The functions that drive a Plotter or a render window are covered by the integration
tests instead, as are the wake ring vortex surfaces, which are built from a history that
only a solved simulation carries. The classes here cover the computation and the
geometry building that feed them, which are settled before any rendering begins.
"""

import tempfile
import unittest
from pathlib import Path

import matplotlib
import matplotlib.colors
import numpy as np
import numpy.testing as npt
import pyvista as pv

# Load PyVista's plotting package, which registers VTK's Matplotlib backend for math
# text, as it is when a diagram creates its Plotter.
import pyvista.plotting  # noqa: F401
import webp
from vtkmodules.vtkRenderingFreeType import vtkMathTextUtilities

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _colormaps, _output_rendering, _transformations
from tests.unit.fixtures import (
    geometry_fixtures,
    operating_point_fixtures,
    output_rendering_fixtures,
    solver_fixtures,
)


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


class TestAnimationWriter(unittest.TestCase):
    """This class contains methods for testing _output_rendering.AnimationWriter."""

    def setUp(self) -> None:
        """Create a temporary directory to hold this test's animation."""
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.animation_path = Path(self.temporary_directory.name) / "animation.webp"

    def tearDown(self) -> None:
        """Remove the temporary directory and the animation it holds."""
        self.temporary_directory.cleanup()

    def test_writes_every_frame_at_the_frame_rate(self) -> None:
        """Test that the file holds every frame, each ending at its cumulative
        timestamp.

        At 25 frames per second, each frame lasts 40 milliseconds, so the five frames
        end at 40, 80, 120, 160, and 200 milliseconds.
        """
        frames = output_rendering_fixtures.make_animation_frames_fixture(5)

        writer = _output_rendering.AnimationWriter(self.animation_path, 25.0, 75.0)
        for frame in frames:
            writer.add_frame(frame)
        writer.close()

        with open(self.animation_path, "rb") as animation_file:
            animation_data = webp.WebPData.from_buffer(animation_file.read())
        decoder = webp.WebPAnimDecoder.new(animation_data)
        self.assertEqual(decoder.anim_info.frame_count, 5)
        self.assertEqual(
            [timestamp for _, timestamp in decoder.frames()], [40, 80, 120, 160, 200]
        )

    def test_matches_the_webp_packages_own_encoder(self) -> None:
        """Test that the file is byte for byte what the webp package's save_images
        writes from the same frames.

        save_images is what the animations were saved with before the writer existed, so
        this is what keeps the saved animations unchanged.
        """
        frames = output_rendering_fixtures.make_animation_frames_fixture(5)

        writer = _output_rendering.AnimationWriter(self.animation_path, 30.0, 50.0)
        for frame in frames:
            writer.add_frame(frame)
        writer.close()

        expected_path = self.animation_path.with_name("expected.webp")
        webp.save_images(
            frames,
            str(expected_path),
            fps=30.0,
            lossless=False,
            quality=50.0,
            method=_output_rendering.WEBP_METHOD,
        )

        self.assertEqual(self.animation_path.read_bytes(), expected_path.read_bytes())

    def test_rejects_a_frame_of_a_different_size(self) -> None:
        """Test that a frame whose size differs from the first's is rejected when the
        writer is closed, and that the frames after it do not block the caller.

        More frames follow the bad one than the writer's queue can hold, so the test
        would hang rather than fail if the writer stopped emptying the queue.
        """
        frames = output_rendering_fixtures.make_animation_frames_fixture(2)
        wrong_size_frame = output_rendering_fixtures.make_animation_frames_fixture(
            1, width=32, height=24
        )[0]
        trailing_frames = output_rendering_fixtures.make_animation_frames_fixture(
            _output_rendering._ANIMATION_WRITER_QUEUE_DEPTH + 2
        )

        writer = _output_rendering.AnimationWriter(self.animation_path, 25.0, 75.0)
        for frame in frames:
            writer.add_frame(frame)
        writer.add_frame(wrong_size_frame)
        for frame in trailing_frames:
            writer.add_frame(frame)

        with self.assertRaises(ValueError) as context:
            writer.close()
        self.assertIn("16 by 12", str(context.exception))
        self.assertIn("32 by 24", str(context.exception))
        self.assertFalse(self.animation_path.exists())

    def test_rejects_an_animation_without_frames(self) -> None:
        """Test that closing a writer that was given no frames is rejected."""
        writer = _output_rendering.AnimationWriter(self.animation_path, 25.0, 75.0)

        with self.assertRaises(ValueError):
            writer.close()
        self.assertFalse(self.animation_path.exists())


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


class TestGetFreeFlightBodyTransformation(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_free_flight_body_transformation."""

    def test_maps_the_cg_onto_its_earth_position(self) -> None:
        """Test that the first Airplane's CG lands at its position in Earth axes.

        The CG is the origin of the body axes the transformation maps out of, so
        wherever the first Airplane has flown to, the CG maps onto that position.
        """
        operating_point = (
            operating_point_fixtures.make_with_cg_position_operating_point_fixture()
        )
        T_pas = _output_rendering.get_free_flight_body_transformation(operating_point)
        Cg_E_Eo = _transformations.apply_T_to_vectors(
            T_pas, np.zeros(3, dtype=float), is_position=True
        )
        npt.assert_allclose(Cg_E_Eo, operating_point.CgP1_E_Eo, atol=1e-12)

    def test_aligns_body_axes_with_earth_axes_at_zero_body_angles(self) -> None:
        """Test that the transformation is the identity for a first Airplane at the
        Earth origin with zero body angles.

        At zero body angles, the geometry axes to Earth axes rotation is the same 180
        degree rotation about the y axis as the constant body axes to geometry axes
        rotation, so chaining the two cancels them and the body axes coincide with Earth
        axes.
        """
        operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        npt.assert_allclose(
            _output_rendering.get_free_flight_body_transformation(operating_point),
            np.eye(4, dtype=float),
            atol=1e-12,
        )

    def test_leaves_directions_untranslated(self) -> None:
        """Test that the translation reaches positions alone.

        A direction has no reference point, so mapping one through this transformation
        must give what chaining the body axes to geometry axes and geometry axes to
        Earth axes rotations alone gives.
        """
        operating_point = (
            operating_point_fixtures.make_with_cg_position_operating_point_fixture()
        )
        direction_BP1 = np.array([1.0, 2.0, 3.0], dtype=float)
        direction_GP1 = _transformations.apply_T_to_vectors(
            operating_point.T_pas_BP1_CgP1_to_GP1_CgP1,
            direction_BP1,
            is_position=False,
        )
        npt.assert_allclose(
            _transformations.apply_T_to_vectors(
                _output_rendering.get_free_flight_body_transformation(operating_point),
                direction_BP1,
                is_position=False,
            ),
            _transformations.apply_T_to_vectors(
                operating_point.T_pas_GP1_CgP1_to_E_CgP1,
                direction_GP1,
                is_position=False,
            ),
            atol=1e-12,
        )


class TestGetMuJoCoRenderGeometry(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_mujoco_render_geometry."""

    def test_splits_worldbody_and_body_geoms(self) -> None:
        """Test that worldbody geoms and body geoms land in their own lists."""
        solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture(
            extra_xml={
                "worldbody": '<geom name="ground" type="plane" size="5 4 0.1"/>',
                "body": '<geom name="box_geom" type="box" size="0.1 0.2 0.3"/>',
            }
        )
        worldbody_geoms, body_geoms = _output_rendering.get_mujoco_render_geometry(
            solver
        )

        self.assertEqual(len(worldbody_geoms), 1)
        self.assertFalse(worldbody_geoms[0].body_attached)
        self.assertEqual(len(body_geoms), 1)
        self.assertTrue(body_geoms[0].body_attached)

    def test_meshes_carry_their_shading_normals(self) -> None:
        """Test that every geom's mesh comes back with active point normals.

        The normals are computed once here so that adding a geom to a frame does not
        recompute them, which is the cost that dominates a frame for a detailed mesh.
        """
        solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture(
            extra_xml={
                "worldbody": '<geom name="ground" type="plane" size="5 4 0.1"/>',
                "body": '<geom name="box_geom" type="box" size="0.1 0.2 0.3"/>',
            }
        )
        worldbody_geoms, body_geoms = _output_rendering.get_mujoco_render_geometry(
            solver
        )
        for render_geom in worldbody_geoms + body_geoms:
            self.assertIsNotNone(render_geom.mesh.point_data.active_normals)

    def test_returns_empty_lists_without_extra_geometry(self) -> None:
        """Test that a solver whose model carries no extra geoms returns two empty
        lists."""
        solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()
        self.assertEqual(_output_rendering.get_mujoco_render_geometry(solver), ([], []))


class TestMathTextAvailability(unittest.TestCase):
    """This class contains methods for testing that VTK can render the diagrams' math
    labels."""

    def test_vtk_can_render_math_text(self) -> None:
        """Test that VTK's Matplotlib backend for math text is available.

        PyVista's plotting package registers the backend when it loads, which happens
        before any diagram creates its Plotter. The backend is only present in VTK
        builds that include VTK's Matplotlib module, and without it the math labels
        would not render as math.
        """
        math_text_utilities = vtkMathTextUtilities.GetInstance()
        self.assertIsNotNone(math_text_utilities)
        assert math_text_utilities is not None
        self.assertTrue(math_text_utilities.IsAvailable())

    def test_vtk_renders_times_math_text_in_stix(self) -> None:
        """Test that VTK renders math text in the Times font family with Matplotlib's
        STIX font set.

        VTK sets Matplotlib's mathtext font set from a Label's font family each time it
        renders math text, and the math labels rely on it choosing STIX for the Times
        family. The font set is first set to something else, so the test only passes if
        rendering changes it. The rc_context restores Matplotlib's settings afterward.
        """
        with matplotlib.rc_context({"mathtext.fontset": "cm"}):
            plotter = pv.Plotter(off_screen=True)
            label = pv.Label(
                text=r"$\hat{\mathbfit{x}}^{\mathrm{G}}$", position=(0.0, 0.0, 0.0)
            )
            label.prop.font_family = "times"
            plotter.add_actor(label)
            plotter.screenshot(return_img=True)
            plotter.close()
            self.assertEqual(matplotlib.rcParams["mathtext.fontset"], "stix")


class TestGetMathAxesLabel(unittest.TestCase):
    """This class contains methods for testing _output_rendering.get_math_axes_label."""

    def test_writes_a_unit_vector_whose_superscript_lists_the_axes_id(self) -> None:
        """Test that each axes ID the diagrams use becomes a bold italic unit vector for
        its basis direction, with the ID's abbreviations and numbers in its
        superscript."""
        cases = [
            ("E", "X", r"\hat{\mathbfit{x}}^{\mathrm{E}}"),
            ("G", "Y", r"\hat{\mathbfit{y}}^{\mathrm{G}}"),
            ("GP1", "Z", r"\hat{\mathbfit{z}}^{\mathrm{G}, \mathrm{P}1}"),
            ("Wn", "X", r"\hat{\mathbfit{x}}^{\mathrm{Wn}}"),
            ("Wn1P2", "X", r"\hat{\mathbfit{x}}^{\mathrm{Wn}1, \mathrm{P}2}"),
            ("Wcsp", "Y", r"\hat{\mathbfit{y}}^{\mathrm{Wcsp}}"),
            ("Wcs1Wn2", "Z", r"\hat{\mathbfit{z}}^{\mathrm{Wcs}1, \mathrm{Wn}2}"),
            ("A", "X", r"\hat{\mathbfit{x}}^{\mathrm{A}}"),
            (
                "AWcs1Wn2P3",
                "Y",
                r"\hat{\mathbfit{y}}^{\mathrm{A}, \mathrm{Wcs}1, \mathrm{Wn}2, "
                r"\mathrm{P}3}",
            ),
        ]
        for axes_id, component_letter, expected_label in cases:
            with self.subTest(axes_id=axes_id, component_letter=component_letter):
                self.assertEqual(
                    _output_rendering.get_math_axes_label(axes_id, component_letter),
                    expected_label,
                )

    def test_rejects_an_invalid_axes_id(self) -> None:
        """Test that an axes ID with an unknown abbreviation, a row and column, or
        characters that aren't abbreviations or numbers raises a ValueError."""
        for axes_id in ["Q", "GQ1", "Wnr1c2", "G-1", "", "g"]:
            with self.subTest(axes_id=axes_id):
                with self.assertRaises(ValueError):
                    _output_rendering.get_math_axes_label(axes_id, "X")


class TestGetMathPointLabel(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_math_point_label."""

    def test_writes_the_name_with_a_subscript_listing_its_owners(self) -> None:
        """Test that each point ID the diagrams use becomes the point's name, followed
        by its own number, a subscript listing its owners, and a Panel point's row and
        column."""
        cases = [
            ("Eo", r"\mathrm{EO}"),
            ("Cg", r"\mathrm{CG}"),
            ("CgP1", r"\mathrm{CG}_{\mathrm{P}1}"),
            ("Ler", r"\mathrm{LER}"),
            ("Ler1", r"\mathrm{LER}1"),
            ("Ler1P2", r"\mathrm{LER}1_{\mathrm{P}2}"),
            ("Lp", r"\mathrm{LP}"),
            ("Lpp", r"\mathrm{LPP}"),
            ("Lp1Wn2", r"\mathrm{LP}1_{\mathrm{Wn}2}"),
            ("Lp1Wn2P3", r"\mathrm{LP}1_{\mathrm{Wn}2, \mathrm{P}3}"),
            ("Cppr3c2", r"\mathrm{CPP}(3, 2)"),
            ("Cppr3c2Wn1", r"\mathrm{CPP}_{\mathrm{Wn}1}(3, 2)"),
            ("Cppr13c21Wn1P2", r"\mathrm{CPP}_{\mathrm{Wn}1, \mathrm{P}2}(13, 21)"),
        ]
        for point_id, expected_label in cases:
            with self.subTest(point_id=point_id):
                self.assertEqual(
                    _output_rendering.get_math_point_label(point_id), expected_label
                )

    def test_rejects_an_invalid_point_id(self) -> None:
        """Test that a point ID with an unknown name, an unknown owner, an owner with a
        row and column, or characters that aren't abbreviations or numbers raises a
        ValueError."""
        for point_id in ["Q1", "CgQ1", "Cppr1c2Wnr1c2", "Cg-1", "", "cg"]:
            with self.subTest(point_id=point_id):
                with self.assertRaises(ValueError):
                    _output_rendering.get_math_point_label(point_id)


class TestGetWingCrossSectionAirfoilLines(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_wing_cross_section_airfoil_lines."""

    def test_lays_the_outline_in_the_xz_plane_scaled_by_chord(self) -> None:
        """Test that the outline's airfoil axes x and y components become its wing cross
        section axes x and z components, scaled by the chord."""
        wing_cross_section = geometry_fixtures.make_basic_wing_cross_section_fixture()
        airfoilOutline_A_Lp = wing_cross_section.airfoil.outline_A_Lp
        airfoilOutline_Wcs_Lp, _ = (
            _output_rendering.get_wing_cross_section_airfoil_lines(wing_cross_section)
        )
        npt.assert_allclose(
            airfoilOutline_Wcs_Lp[:, 0],
            wing_cross_section.chord * airfoilOutline_A_Lp[:, 0],
        )
        npt.assert_array_equal(airfoilOutline_Wcs_Lp[:, 1], 0.0)
        npt.assert_allclose(
            airfoilOutline_Wcs_Lp[:, 2],
            wing_cross_section.chord * airfoilOutline_A_Lp[:, 1],
        )

    def test_lays_the_mcl_in_the_xz_plane_scaled_by_chord(self) -> None:
        """Test that the mean camber line's airfoil axes x and y components become its
        wing cross section axes x and z components, scaled by the chord."""
        wing_cross_section = geometry_fixtures.make_basic_wing_cross_section_fixture()
        airfoilMcl_A_Lp = wing_cross_section.airfoil.mcl_A_Lp
        assert airfoilMcl_A_Lp is not None
        _, airfoilMcl_Wcs_Lp = _output_rendering.get_wing_cross_section_airfoil_lines(
            wing_cross_section
        )
        npt.assert_allclose(
            airfoilMcl_Wcs_Lp[:, 0], wing_cross_section.chord * airfoilMcl_A_Lp[:, 0]
        )
        npt.assert_array_equal(airfoilMcl_Wcs_Lp[:, 1], 0.0)
        npt.assert_allclose(
            airfoilMcl_Wcs_Lp[:, 2], wing_cross_section.chord * airfoilMcl_A_Lp[:, 1]
        )


class TestGetCollocationPoints(unittest.TestCase):
    """This class contains methods for testing
    _output_rendering.get_collocation_points."""

    def test_returns_every_panel_row_by_row(self) -> None:
        """Test that every Panel's collocation point is returned, row by row."""
        wing = output_rendering_fixtures.make_placed_airplanes_fixture()[0].wings[0]
        assert wing.panels is not None
        num_rows, num_columns = wing.panels.shape
        ids, _, _, _ = _output_rendering.get_collocation_points(
            wing, "", np.eye(4, dtype=float)
        )
        self.assertEqual(
            ids,
            [
                f"Cppr{row_num}c{column_num}"
                for row_num in range(1, num_rows + 1)
                for column_num in range(1, num_columns + 1)
            ],
        )

    def test_numbers_panels_from_one_with_the_suffix(self) -> None:
        """Test that each ID numbers its Panel's row and column from one and ends with
        the suffix."""
        wing = output_rendering_fixtures.make_placed_airplanes_fixture()[0].wings[0]
        assert wing.panels is not None
        num_columns = wing.panels.shape[1]
        ids, _, _, _ = _output_rendering.get_collocation_points(
            wing, "Wn1", np.eye(4, dtype=float)
        )
        self.assertEqual(ids[1], "Cppr1c2Wn1")
        self.assertEqual(ids[num_columns], "Cppr2c1Wn1")

    def test_maps_each_collocation_point_into_diagram_axes(self) -> None:
        """Test that each position is its Panel's collocation point, mapped through the
        transformation."""
        wing = output_rendering_fixtures.make_placed_airplanes_fixture()[0].wings[0]
        assert wing.panels is not None
        T_pas_G_Cg_to_D_Do = _transformations.generate_trans_T(
            np.array([1.0, -2.0, 0.5]), passive=True
        )
        num_columns = wing.panels.shape[1]
        _, listCollocationPoints_D_Do, _, _ = _output_rendering.get_collocation_points(
            wing, "", T_pas_G_Cg_to_D_Do
        )
        npt.assert_allclose(
            listCollocationPoints_D_Do[num_columns + 2],
            _transformations.apply_T_to_vectors(
                T_pas_G_Cg_to_D_Do, wing.panels[1, 2].Cpp_G_Cg, is_position=True
            ),
        )

    def test_lays_each_cross_along_its_panels_diagonals(self) -> None:
        """Test that each cross's arms are unit vectors along its Panel's diagonals."""
        wing = output_rendering_fixtures.make_placed_airplanes_fixture()[0].wings[0]
        assert wing.panels is not None
        panel = wing.panels[0, 0]
        _, _, _, listCrossDirections_D = _output_rendering.get_collocation_points(
            wing, "", np.eye(4, dtype=float)
        )
        firstDiagonal_G = panel.Brpp_G_Cg - panel.Flpp_G_Cg
        secondDiagonal_G = panel.Blpp_G_Cg - panel.Frpp_G_Cg
        npt.assert_allclose(
            listCrossDirections_D[0],
            [
                firstDiagonal_G / np.linalg.norm(firstDiagonal_G),
                secondDiagonal_G / np.linalg.norm(secondDiagonal_G),
            ],
        )

    def test_returns_nothing_for_an_unmeshed_wing(self) -> None:
        """Test that a Wing without Panels gives empty lists."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        self.assertIsNone(wing.panels)
        self.assertEqual(
            _output_rendering.get_collocation_points(wing, "", np.eye(4, dtype=float)),
            ([], [], [], []),
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

    def test_maps_the_normals_as_directions(self) -> None:
        """Test that a mesh's active point normals are rotated but not translated.

        A body mesh carries its shading normals across the frames of an animation, so
        they have to turn with the faces they shade while ignoring where those faces
        move to.
        """
        mesh = output_rendering_fixtures.make_cube_mesh_fixture().compute_normals(
            cell_normals=False
        )
        T_pas = _transformations.compose_T_pas(
            _transformations.generate_rot_T(
                angles=np.array([0.0, 0.0, np.pi / 2], dtype=float),
                passive=True,
                intrinsic=True,
                order="xyz",
            ),
            _transformations.generate_trans_T(
                translations=np.array([1.0, 2.0, 3.0], dtype=float), passive=True
            ),
        )
        transformed = _output_rendering.transform_mesh(mesh, T_pas)
        npt.assert_allclose(
            transformed.point_data["Normals"],
            np.array(mesh.point_data["Normals"]) @ T_pas[:3, :3].T,
            atol=1e-12,
        )


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
        muted = _output_rendering.mute_colormap(_colormaps.SEQUENTIAL_COLOR_MAP, 0.0)
        npt.assert_allclose(muted(0.0), _colormaps.SEQUENTIAL_COLOR_MAP(0.0))
        npt.assert_allclose(muted(1.0), _colormaps.SEQUENTIAL_COLOR_MAP(1.0))

    def test_a_factor_of_one_returns_middle_gray(self) -> None:
        """Test that muting a color map fully leaves every color middle gray."""
        muted = _output_rendering.mute_colormap(_colormaps.SEQUENTIAL_COLOR_MAP, 1.0)
        npt.assert_allclose(muted(0.0)[:3], (0.5, 0.5, 0.5))
        npt.assert_allclose(muted(0.5)[:3], (0.5, 0.5, 0.5))
        npt.assert_allclose(muted(1.0)[:3], (0.5, 0.5, 0.5))

    def test_returns_a_listed_color_map_of_256_colors(self) -> None:
        """Test that the muted color map is a ListedColormap sampled at 256 colors."""
        muted = _output_rendering.mute_colormap(_colormaps.SEQUENTIAL_COLOR_MAP, 0.5)
        self.assertIsInstance(muted, matplotlib.colors.ListedColormap)
        self.assertEqual(muted.N, 256)

    def test_leaves_the_alpha_channel_alone(self) -> None:
        """Test that muting reaches the colors rather than their opacity."""
        muted = _output_rendering.mute_colormap(_colormaps.SEQUENTIAL_COLOR_MAP, 1.0)
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

    def test_crosswind_force_negates_the_wind_axes_y_force(self) -> None:
        """Test that crosswind force reads as positive toward the Airplane's left.

        Wind axes y points to the Airplane's right, so the coefficient carries the
        negated force.
        """
        airplanes = output_rendering_fixtures.make_loaded_airplanes_fixture()
        assert airplanes[0].wings[0].panels is not None
        panels = np.ravel(airplanes[0].wings[0].panels)
        scalars = _output_rendering.get_scalars(airplanes, "crosswind force", 2.0)
        expected = [-panel.forces_W[1] / 2.0 / panel.area for panel in panels]
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
        self.assertIs(color_map, _colormaps.SEQUENTIAL_COLOR_MAP)

    def test_sign_changing_scalars_get_the_diverging_color_map(self) -> None:
        """Test that scalars that change sign are colored about their midpoint."""
        scalars = np.array([-2.0, -1.0, 1.0, 2.0], dtype=float)
        color_map, _, _ = _output_rendering.choose_color_map(scalars)
        self.assertIs(color_map, _colormaps.DIVERGING_COLOR_MAP)

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
