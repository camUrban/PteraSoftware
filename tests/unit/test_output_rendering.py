"""This module contains classes to test the output rendering functions.

The functions that show a Plotter or drive a render window are covered by the
integration tests instead, as are the wake ring vortex surfaces, which are built from a
history that only a solved simulation carries. The classes here cover the computation
and the geometry building that feed them, which are settled before any rendering begins.
They also cover add_vortices, which only adds meshes to a Plotter, so its meshes can be
checked on an off screen Plotter that never renders. Likewise, they cover
add_axes_and_points, whose meshes and labels are checked the same way, and whose label
placement and dragging are checked by rendering that off screen Plotter and sending its
interactor style mouse events.
"""

import math
import tempfile
import unittest
from pathlib import Path
from typing import cast

import matplotlib
import matplotlib.colors
import numpy as np
import numpy.testing as npt
import pyvista as pv

# Load PyVista's plotting package, which registers VTK's Matplotlib backend for math
# text, as it is when a diagram creates its Plotter.
import pyvista.plotting  # noqa: F401
import webp
from vtkmodules.vtkCommonCore import reference as vtk_reference
from vtkmodules.vtkRenderingFreeType import vtkMathTextUtilities

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _colormaps, _fonts, _output_rendering, _transformations
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


def _add_vortices(
    plotter: pv.Plotter,
    ring_vortices: tuple[np.ndarray, ...] | None = None,
    wake_ring_vortices: tuple[np.ndarray, ...] | None = None,
    horseshoe_vortices: tuple[np.ndarray, ...] | None = None,
    horseshoe_vortices_are_wake: bool = False,
    largest_chord: float = 1.0,
    simplify: bool = False,
) -> None:
    """Adds vortices to a Plotter with _output_rendering.add_vortices, filling in empty
    stacks for every kind of vortex that isn't given.

    :param plotter: The Plotter to add the vortices to.
    :param ring_vortices: A tuple of five (N,3) ndarrays of floats holding the ring
        vortices' front right, front left, back left, and back right points (in diagram
        axes, relative to the diagram origin), and their Panels' unit normals (in
        diagram axes), or None to pass no ring vortices. The default is None.
    :param wake_ring_vortices: A tuple of five (K,3) ndarrays of floats holding the wake
        ring vortices' points and unit normals in the same order, or None to pass no
        wake ring vortices. The default is None.
    :param horseshoe_vortices: A tuple of five (M,3) ndarrays of floats holding the
        horseshoe vortices' points and unit normals in the same order, or None to pass
        no horseshoe vortices. The default is None.
    :param horseshoe_vortices_are_wake: Determines whether the horseshoe vortices are
        drawn as wake vortices. The default is False.
    :param largest_chord: The largest chord, which scales the trailing legs' overhang
        and their dashes. The units are in meters. The default is 1.0.
    :param simplify: Determines whether to simplify the vortices. The default is False.
    :return: None
    """
    noVortexPoints_D_Do = output_rendering_fixtures.make_no_vortex_points_fixture()
    no_vortices = (noVortexPoints_D_Do,) * 5
    (
        stackFrrvp_D_Do,
        stackFlrvp_D_Do,
        stackBlrvp_D_Do,
        stackBrrvp_D_Do,
        stackRingUnitNormals_D,
    ) = (
        ring_vortices or no_vortices
    )
    (
        stackFrwrvp_D_Do,
        stackFlwrvp_D_Do,
        stackBlwrvp_D_Do,
        stackBrwrvp_D_Do,
        stackWakeRingUnitNormals_D,
    ) = (
        wake_ring_vortices or no_vortices
    )
    (
        stackFrhvp_D_Do,
        stackFlhvp_D_Do,
        stackBlhvp_D_Do,
        stackBrhvp_D_Do,
        stackHorseshoeUnitNormals_D,
    ) = (
        horseshoe_vortices or no_vortices
    )
    _output_rendering.add_vortices(
        plotter,
        stackFrrvp_D_Do=stackFrrvp_D_Do,
        stackFlrvp_D_Do=stackFlrvp_D_Do,
        stackBlrvp_D_Do=stackBlrvp_D_Do,
        stackBrrvp_D_Do=stackBrrvp_D_Do,
        stackRingUnitNormals_D=stackRingUnitNormals_D,
        stackFrwrvp_D_Do=stackFrwrvp_D_Do,
        stackFlwrvp_D_Do=stackFlwrvp_D_Do,
        stackBlwrvp_D_Do=stackBlwrvp_D_Do,
        stackBrwrvp_D_Do=stackBrwrvp_D_Do,
        stackWakeRingUnitNormals_D=stackWakeRingUnitNormals_D,
        stackFrhvp_D_Do=stackFrhvp_D_Do,
        stackFlhvp_D_Do=stackFlhvp_D_Do,
        stackBlhvp_D_Do=stackBlhvp_D_Do,
        stackBrhvp_D_Do=stackBrhvp_D_Do,
        stackHorseshoeUnitNormals_D=stackHorseshoeUnitNormals_D,
        horseshoe_vortices_are_wake=horseshoe_vortices_are_wake,
        largest_chord=largest_chord,
        simplify=simplify,
    )


def _get_line_actors(plotter: pv.Plotter) -> list[tuple[pv.PolyData, pv.Actor]]:
    """Returns the meshes made of lines in a Plotter, along with their Actors, in the
    order they were added.

    An arrow tip's filled faces and its outline are left out. The faces' mesh has no
    lines, and the outline's Actor draws a filter's output, so its mapper's mesh is
    empty before the Plotter renders. Labels are left out too, since they aren't Actors
    and draw no mesh.

    :param plotter: The Plotter whose meshes to return.
    :return: A list of tuples, each holding a mesh made of lines and the Actor that
        draws it.
    """
    line_actors: list[tuple[pv.PolyData, pv.Actor]] = []
    for actor in plotter.actors.values():
        if not isinstance(actor, pv.Actor):
            continue
        mesh = cast(pv.DataSetMapper, actor.mapper).dataset
        if isinstance(mesh, pv.PolyData) and mesh.n_lines > 0:
            line_actors.append((mesh, actor))
    return line_actors


class TestAddVortices(unittest.TestCase):
    """This class contains methods for testing _output_rendering.add_vortices."""

    def setUp(self) -> None:
        """Create an off screen Plotter to add this test's vortices to."""
        self.plotter = pv.Plotter(off_screen=True)

    def tearDown(self) -> None:
        """Close the Plotter."""
        self.plotter.close()

    def test_adds_nothing_without_vortices(self) -> None:
        """Test that a solver that places no vortices adds nothing to the Plotter."""
        _add_vortices(self.plotter)
        self.assertEqual(len(self.plotter.actors), 0)

    def test_closes_each_exact_ring_vortex(self) -> None:
        """Test that an exact ring vortex is one closed polyline through its corners.

        The polyline runs front right, front left, back left, back right, and front
        right again, so it repeats its first corner to close itself.
        """
        ring_vortices = output_rendering_fixtures.make_square_ring_vortex_fixture()
        stackFrrvp_D_Do, stackFlrvp_D_Do, stackBlrvp_D_Do, stackBrrvp_D_Do, _ = (
            ring_vortices
        )
        _add_vortices(self.plotter, ring_vortices=ring_vortices)
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 1)
        mesh, _ = line_actors[0]
        npt.assert_array_equal(
            mesh.points,
            np.vstack(
                [
                    stackFrrvp_D_Do,
                    stackFlrvp_D_Do,
                    stackBlrvp_D_Do,
                    stackBrrvp_D_Do,
                    stackFrrvp_D_Do,
                ]
            ),
        )
        npt.assert_array_equal(mesh.lines, [5, 0, 1, 2, 3, 4])

    def test_exact_vortices_get_no_vorticity_arrows(self) -> None:
        """Test that an exact ring vortex adds its polyline and nothing else."""
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
        )
        self.assertEqual(len(self.plotter.actors), 1)

    def test_draws_one_mesh_per_color(self) -> None:
        """Test that the ring vortices sharing a color share one mesh.

        The bound ring vortices come first, in one mesh of one polyline each drawn in
        the vortex color at the vortex line width, and the wake ring vortices follow in
        a mesh of their own drawn in the wake vortex color.
        """
        _add_vortices(
            self.plotter,
            ring_vortices=(
                output_rendering_fixtures.make_neighboring_ring_vortices_fixture()
            ),
            wake_ring_vortices=(
                output_rendering_fixtures.make_square_ring_vortex_fixture()
            ),
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 2)
        bound_mesh, bound_actor = line_actors[0]
        wake_mesh, wake_actor = line_actors[1]
        self.assertEqual(bound_mesh.n_cells, 2)
        self.assertEqual(bound_actor.prop.color, _output_rendering._VORTEX_COLOR)
        self.assertEqual(
            bound_actor.prop.line_width, _output_rendering._VORTEX_LINE_WIDTH
        )
        self.assertEqual(wake_mesh.n_cells, 1)
        self.assertEqual(wake_actor.prop.color, _output_rendering._VORTEX_WAKE_COLOR)

    def test_ends_trailing_legs_a_fixed_overhang_past_the_bounding_box(self) -> None:
        """Test that an exact horseshoe vortex's solid polyline runs from its right
        trailing leg's solid end, along its finite leg, to its left trailing leg's solid
        end.

        The vortices' bounding box is flat in x, so each trailing leg ends 1.0 meter
        downstream of the finite leg, and its solid part ends 0.5 meters downstream,
        where its dashed length begins. The trailing legs' far away back points don't
        affect where they end.
        """
        horseshoe_vortices = (
            output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
        )
        _add_vortices(self.plotter, horseshoe_vortices=horseshoe_vortices)
        solid_mesh, _ = _get_line_actors(self.plotter)[0]
        npt.assert_allclose(
            solid_mesh.points,
            [
                [0.5, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
            ],
            atol=1e-12,
        )
        npt.assert_array_equal(solid_mesh.lines, [4, 0, 1, 2, 3])

    def test_dashes_the_end_of_each_trailing_leg(self) -> None:
        """Test that each trailing leg's dashed length is drawn as evenly spaced dashes.

        The dashed length runs from 0.5 to 1.0 meters downstream, and each dash and gap
        is 0.05 meters long, so each trailing leg gets five dashes, starting every 0.1
        meters from 0.5 meters. The right trailing leg's dashes come first.
        """
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
            ),
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 2)
        dash_mesh, _ = line_actors[1]
        dash_starts = 0.5 + 0.1 * np.arange(5, dtype=float)
        dash_distances = np.ravel(np.column_stack([dash_starts, dash_starts + 0.05]))
        expectedRightDashVertices_D_Do = np.column_stack(
            [dash_distances, np.ones(10, dtype=float), np.zeros(10, dtype=float)]
        )
        expectedLeftDashVertices_D_Do = np.column_stack(
            [dash_distances, np.zeros(10, dtype=float), np.zeros(10, dtype=float)]
        )
        self.assertEqual(dash_mesh.n_cells, 10)
        npt.assert_allclose(
            dash_mesh.points,
            np.vstack([expectedRightDashVertices_D_Do, expectedLeftDashVertices_D_Do]),
            atol=1e-12,
        )

    def test_lines_up_dashes_of_overlapping_trailing_legs(self) -> None:
        """Test that two trailing legs that lie on top of each other get matching
        dashes, even though they start at different points.

        The dashes are spaced by distance along the trailing legs' shared direction,
        rather than from each trailing leg's start, so the first horseshoe vortex's
        right trailing leg and the second one's left trailing leg draw the same dashes.
        Each of the four trailing legs gets five dashes, ordered by horseshoe vortex and
        then right before left.
        """
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_staggered_horseshoe_vortices_fixture()
            ),
        )
        dash_mesh, _ = _get_line_actors(self.plotter)[1]
        self.assertEqual(dash_mesh.n_cells, 20)
        npt.assert_allclose(dash_mesh.points[0:10], dash_mesh.points[30:40], atol=1e-12)

    def test_scales_the_overhang_with_the_largest_chord(self) -> None:
        """Test that doubling the largest chord doubles the trailing legs' overhang and
        their dashed length.

        With a largest chord of 2.0 meters, each trailing leg ends 2.0 meters downstream
        of the finite leg, and its solid part ends 1.0 meter downstream.
        """
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
            ),
            largest_chord=2.0,
        )
        solid_mesh, _ = _get_line_actors(self.plotter)[0]
        npt.assert_allclose(solid_mesh.points[0], [1.0, 1.0, 0.0], atol=1e-12)
        npt.assert_allclose(solid_mesh.points[-1], [1.0, 0.0, 0.0], atol=1e-12)

    def test_runs_trailing_legs_toward_their_back_points(self) -> None:
        """Test that each trailing leg runs along the unit vector from its front point
        toward its back point, whatever direction that is."""
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_slanted_horseshoe_vortex_fixture()
            ),
        )
        solid_mesh, _ = _get_line_actors(self.plotter)[0]
        trailingDirection_D = np.array([1.0, 0.0, 1.0], dtype=float) / math.sqrt(2.0)
        npt.assert_allclose(
            solid_mesh.points[0],
            np.array([0.0, 1.0, 0.0], dtype=float) + 0.5 * trailingDirection_D,
            atol=1e-12,
        )
        npt.assert_allclose(
            solid_mesh.points[-1], 0.5 * trailingDirection_D, atol=1e-12
        )

    def test_draws_bound_horseshoe_vortices_in_the_vortex_color(self) -> None:
        """Test that bound horseshoe vortices' solid polylines and dashes are drawn in
        the vortex color."""
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
            ),
            horseshoe_vortices_are_wake=False,
        )
        for _, actor in _get_line_actors(self.plotter):
            self.assertEqual(actor.prop.color, _output_rendering._VORTEX_COLOR)

    def test_draws_wake_horseshoe_vortices_in_the_wake_vortex_color(self) -> None:
        """Test that wake horseshoe vortices' solid polylines and dashes are drawn in
        the wake vortex color, while the bound ring vortices drawn with them keep the
        vortex color.

        This is the combination a steady ring vortex lattice method solver draws, with
        its bound ring vortices first and its wake horseshoe vortices' solid polylines
        and dashes after them.
        """
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
            horseshoe_vortices=(
                output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
            ),
            horseshoe_vortices_are_wake=True,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 3)
        self.assertEqual(line_actors[0][1].prop.color, _output_rendering._VORTEX_COLOR)
        self.assertEqual(
            line_actors[1][1].prop.color, _output_rendering._VORTEX_WAKE_COLOR
        )
        self.assertEqual(
            line_actors[2][1].prop.color, _output_rendering._VORTEX_WAKE_COLOR
        )

    def test_simplified_vortices_use_the_simplified_colors(self) -> None:
        """Test that simplified bound and wake ring vortices are drawn in their
        simplified colors."""
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
            wake_ring_vortices=(
                output_rendering_fixtures.make_square_ring_vortex_fixture()
            ),
            simplify=True,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(
            line_actors[0][1].prop.color, _output_rendering._VORTEX_SIMPLIFIED_COLOR
        )
        self.assertEqual(
            line_actors[1][1].prop.color,
            _output_rendering._VORTEX_WAKE_SIMPLIFIED_COLOR,
        )

    def test_shrinks_and_rounds_a_simplified_ring_vortex(self) -> None:
        """Test that a simplified ring vortex is shrunk toward its center and has its
        corners rounded.

        The shrunk ring vortex's corners sit 0.45 meters from its center along the x and
        y axes, so its legs span 0.05 to 0.95 meters. Each of its four corners is
        replaced by a curve of eight points, and the closed polyline repeats its first
        point at its end. The rounded polyline never reaches the shrunk corners.
        """
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
            simplify=True,
        )
        solid_mesh, _ = _get_line_actors(self.plotter)[0]
        num_corner_points = _output_rendering._VORTEX_SIMPLIFIED_CORNER_NUM_POINTS
        self.assertEqual(solid_mesh.n_points, 4 * num_corner_points + 1)
        npt.assert_allclose(solid_mesh.points[0], solid_mesh.points[-1], atol=1e-12)
        npt.assert_allclose(
            solid_mesh.bounds, (0.05, 0.95, 0.05, 0.95, 0.0, 0.0), atol=1e-6
        )
        for shrunkCorner_D_Do in [
            [0.05, 0.95, 0.0],
            [0.05, 0.05, 0.0],
            [0.95, 0.05, 0.0],
            [0.95, 0.95, 0.0],
        ]:
            self.assertGreater(
                float(
                    np.min(
                        np.linalg.norm(solid_mesh.points - shrunkCorner_D_Do, axis=1)
                    )
                ),
                1e-3,
            )

    def test_gives_each_leg_of_a_simplified_ring_vortex_an_arrow(self) -> None:
        """Test that a simplified ring vortex gets one vorticity arrow per leg.

        The arrows' shafts share one mesh of four half circular polylines, drawn at the
        axes arrows' line width, and their tips share one outline and one set of filled
        faces.
        """
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
            simplify=True,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 2)
        arc_mesh, arc_actor = line_actors[1]
        self.assertEqual(arc_mesh.n_cells, 4)
        self.assertEqual(
            arc_mesh.n_points, 4 * _output_rendering._VORTEX_VORTICITY_ARROW_NUM_POINTS
        )
        self.assertEqual(arc_actor.prop.line_width, _output_rendering._AXES_LINE_WIDTH)
        self.assertEqual(
            arc_actor.prop.color, _output_rendering._VORTEX_SIMPLIFIED_COLOR
        )
        self.assertEqual(len(self.plotter.actors), 2 + 2)

    def test_centers_each_vorticity_arrow_on_its_legs_midpoint(self) -> None:
        """Test that a vorticity arrow is a half circle around its leg's midpoint, in
        the plane perpendicular to the leg.

        The first leg runs from the shrunk front right corner to the shrunk front left
        corner, so its midpoint is (0.05, 0.5, 0.0). Every leg is 0.9 meters long, so
        each arrow's radius is 0.15 times that.
        """
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
            simplify=True,
        )
        arc_mesh, _ = _get_line_actors(self.plotter)[1]
        num_arc_points = _output_rendering._VORTEX_VORTICITY_ARROW_NUM_POINTS
        firstArcPoints_D_Do = arc_mesh.points[:num_arc_points]
        legMidpoint_D_Do = np.array([0.05, 0.5, 0.0], dtype=float)
        arrow_radius = _output_rendering._VORTEX_VORTICITY_ARROW_RADIUS * 0.9
        npt.assert_allclose(
            np.linalg.norm(firstArcPoints_D_Do - legMidpoint_D_Do, axis=1),
            arrow_radius,
            atol=1e-6,
        )
        npt.assert_allclose(firstArcPoints_D_Do[:, 1], 0.5, atol=1e-6)

    def test_sweeps_each_vorticity_arrow_across_its_legs_inward_side(self) -> None:
        """Test that a vorticity arrow starts on its Panel's upper side, sweeps across
        the inside of its vortex, and stops where its tip begins.

        The first leg runs along the negative y direction, so the vorticity for a
        negative vortex strength points along the positive y direction. Turning about
        that direction carries the upper side, along the positive z direction, toward
        the inside of the ring vortex, along the positive x direction. The tip takes up
        a fifth of the half circle's length, so the shaft stops at 0.8 * pi radians.
        """
        _add_vortices(
            self.plotter,
            ring_vortices=output_rendering_fixtures.make_square_ring_vortex_fixture(),
            simplify=True,
        )
        arc_mesh, _ = _get_line_actors(self.plotter)[1]
        num_arc_points = _output_rendering._VORTEX_VORTICITY_ARROW_NUM_POINTS
        firstArcPoints_D_Do = arc_mesh.points[:num_arc_points]
        legMidpoint_D_Do = np.array([0.05, 0.5, 0.0], dtype=float)
        arrow_radius = _output_rendering._VORTEX_VORTICITY_ARROW_RADIUS * 0.9
        shaft_end_angle = math.pi - _output_rendering._AXES_TIP_LENGTH * math.pi
        npt.assert_allclose(
            firstArcPoints_D_Do[0],
            legMidpoint_D_Do + arrow_radius * np.array([0.0, 0.0, 1.0]),
            atol=1e-6,
        )
        npt.assert_allclose(
            firstArcPoints_D_Do[-1],
            legMidpoint_D_Do
            + arrow_radius
            * np.array(
                [math.sin(shaft_end_angle), 0.0, math.cos(shaft_end_angle)],
                dtype=float,
            ),
            atol=1e-6,
        )
        self.assertTrue(np.all(firstArcPoints_D_Do[:, 0] >= 0.05 - 1e-6))

    def test_shrinks_the_finite_leg_of_a_simplified_horseshoe_vortex(self) -> None:
        """Test that a simplified horseshoe vortex's finite leg is shrunk toward its
        midpoint, with its trailing legs moving along with the finite leg's ends.

        The finite leg's ends move to 0.05 and 0.95 meters along the y axis. The
        trailing legs still end where the unshrunk vortices' bounding box puts them, so
        their solid parts end 0.5 meters downstream. The open polyline keeps its two end
        points, and each of its two corners is replaced by a curve of eight points.
        """
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
            ),
            simplify=True,
        )
        solid_mesh, _ = _get_line_actors(self.plotter)[0]
        num_corner_points = _output_rendering._VORTEX_SIMPLIFIED_CORNER_NUM_POINTS
        self.assertEqual(solid_mesh.n_points, 2 * num_corner_points + 2)
        npt.assert_allclose(solid_mesh.points[0], [0.5, 0.95, 0.0], atol=1e-6)
        npt.assert_allclose(solid_mesh.points[-1], [0.5, 0.05, 0.0], atol=1e-6)

    def test_gives_each_solid_leg_of_a_simplified_horseshoe_vortex_an_arrow(
        self,
    ) -> None:
        """Test that a simplified horseshoe vortex gets vorticity arrows on its finite
        leg and its two trailing legs' solid parts.

        The finite leg's arrow comes first. The inside of a horseshoe vortex lies
        downstream of its finite leg, so that arrow starts on the Panel's upper side and
        sweeps downstream, along the positive x direction. Its radius is 0.15 times the
        shrunk finite leg's length of 0.9 meters.
        """
        _add_vortices(
            self.plotter,
            horseshoe_vortices=(
                output_rendering_fixtures.make_straight_horseshoe_vortex_fixture()
            ),
            simplify=True,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 3)
        arc_mesh, _ = line_actors[2]
        self.assertEqual(arc_mesh.n_cells, 3)
        num_arc_points = _output_rendering._VORTEX_VORTICITY_ARROW_NUM_POINTS
        finiteLegArcPoints_D_Do = arc_mesh.points[:num_arc_points]
        arrow_radius = _output_rendering._VORTEX_VORTICITY_ARROW_RADIUS * 0.9
        npt.assert_allclose(
            finiteLegArcPoints_D_Do[0], [0.0, 0.5, arrow_radius], atol=1e-6
        )
        self.assertTrue(np.all(finiteLegArcPoints_D_Do[:, 0] >= -1e-6))


def _get_tip_fill_actors(plotter: pv.Plotter) -> list[tuple[pv.PolyData, pv.Actor]]:
    """Returns the meshes of the arrow tips' filled faces in a Plotter, along with their
    Actors, in the order they were added.

    These are the only meshes made of faces alone. The tips' outlines are left out,
    since their Actors draw a filter's output, so their mappers' meshes are empty before
    the Plotter renders.

    :param plotter: The Plotter whose meshes to return.
    :return: A list of tuples, each holding a mesh of arrow tips' filled faces and the
        Actor that draws it.
    """
    tip_fill_actors: list[tuple[pv.PolyData, pv.Actor]] = []
    for actor in plotter.actors.values():
        if not isinstance(actor, pv.Actor):
            continue
        mesh = cast(pv.DataSetMapper, actor.mapper).dataset
        if (
            isinstance(mesh, pv.PolyData)
            and mesh.n_cells > 0
            and mesh.n_lines == 0
            and mesh.n_verts == 0
        ):
            tip_fill_actors.append((mesh, actor))
    return tip_fill_actors


def _get_dot_actors(plotter: pv.Plotter) -> list[tuple[pv.PolyData, pv.Actor]]:
    """Returns the meshes of points drawn as dots in a Plotter, along with their Actors,
    in the order they were added.

    :param plotter: The Plotter whose meshes to return.
    :return: A list of tuples, each holding a mesh made of vertices and the Actor that
        draws it.
    """
    dot_actors: list[tuple[pv.PolyData, pv.Actor]] = []
    for actor in plotter.actors.values():
        if not isinstance(actor, pv.Actor):
            continue
        mesh = cast(pv.DataSetMapper, actor.mapper).dataset
        if isinstance(mesh, pv.PolyData) and mesh.n_verts > 0:
            dot_actors.append((mesh, actor))
    return dot_actors


def _get_labels(plotter: pv.Plotter) -> dict[str, pv.Label]:
    """Returns the Labels in a Plotter, keyed by their names.

    add_axes_and_points names each arrow's Label "arrow label " followed by its plain
    text, and each point's Label "point label " followed by its plain text, even when
    the Label displays math.

    :param plotter: The Plotter whose Labels to return.
    :return: A dict mapping each Label's name to the Label.
    """
    return {
        name: actor
        for name, actor in plotter.actors.items()
        if isinstance(actor, pv.Label)
    }


def _get_display_point(plotter: pv.Plotter, point_D_Do: np.ndarray) -> np.ndarray:
    """Returns the display coordinates of a point in a Plotter's renderer.

    :param plotter: The Plotter whose renderer maps the point.
    :param point_D_Do: A (3,) ndarray of floats representing the point's position (in
        diagram axes, relative to the diagram origin). The units are in meters.
    :return: A (3,) ndarray of floats holding the point's x and y display coordinates,
        in pixels, and its depth.
    """
    renderer = plotter.renderer
    renderer.SetWorldPoint(*point_D_Do, 1.0)
    renderer.WorldToDisplay()
    return np.array(renderer.GetDisplayPoint(), dtype=float)


def _drag_mouse(
    plotter: pv.Plotter, start_display: tuple[int, int], end_display: tuple[int, int]
) -> None:
    """Drags the mouse with its left button held down across a Plotter's render window.

    The events are sent to the interactor style, which is what add_axes_and_points
    observes, so the drag reaches its observers the way a real one does.

    :param plotter: The Plotter whose render window the mouse is dragged across.
    :param start_display: The x and y display coordinates, in pixels, at which the left
        button is pressed.
    :param end_display: The x and y display coordinates, in pixels, to which the mouse
        is moved before the left button is released.
    :return: None
    """
    assert plotter.iren is not None
    interactor = plotter.iren.interactor
    interactor_style = interactor.GetInteractorStyle()
    interactor.SetEventPosition(*start_display)
    interactor_style.InvokeEvent("LeftButtonPressEvent")
    interactor.SetEventPosition(*end_display)
    interactor_style.InvokeEvent("MouseMoveEvent")
    interactor_style.InvokeEvent("LeftButtonReleaseEvent")


class TestAddAxesAndPoints(unittest.TestCase):
    """This class contains methods for testing _output_rendering.add_axes_and_points."""

    def setUp(self) -> None:
        """Create an off screen Plotter to add this test's axes and points to."""
        self.plotter = pv.Plotter(off_screen=True)

    def tearDown(self) -> None:
        """Close the Plotter."""
        self.plotter.close()

    def test_sets_the_diagram_background_color(self) -> None:
        """Test that the Plotter's background is set to the diagram background color."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        self.assertEqual(
            self.plotter.background_color,
            pv.Color(_output_rendering._DIAGRAM_BACKGROUND_COLOR),
        )

    def test_draws_one_colored_shaft_per_basis_direction(self) -> None:
        """Test that an axes set's x, y, and z shafts are each drawn as one line, in
        red, green, and blue, from its point to where its tip begins.

        Each shaft starts at the axes set's point and runs along its basis direction, as
        the transformation gives them. The tip takes up a fifth of each arrow's length,
        so each shaft is 0.8 times the arrows' length.
        """
        T_pas_G_Cg_to_D_Do = (
            output_rendering_fixtures.make_offset_axes_transformation_fixture()
        )
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[T_pas_G_Cg_to_D_Do],
            axes_scale=axes_scale,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 3)
        Cg_D_Do = T_pas_G_Cg_to_D_Do[:3, 3]
        shaft_length = (1.0 - _output_rendering._AXES_TIP_LENGTH) * axes_scale
        for component_id, (mesh, actor) in enumerate(line_actors):
            with self.subTest(component_id=component_id):
                npt.assert_allclose(
                    mesh.points,
                    [
                        Cg_D_Do,
                        Cg_D_Do + shaft_length * T_pas_G_Cg_to_D_Do[:3, component_id],
                    ],
                    atol=1e-6,
                )
                self.assertEqual(
                    actor.prop.color,
                    pv.Color(_output_rendering._AXES_COLORS[component_id]),
                )
                self.assertEqual(
                    actor.prop.line_width, _output_rendering._AXES_LINE_WIDTH
                )

    def test_draws_only_the_x_and_y_arrows_of_two_dimensional_axes(self) -> None:
        """Test that a two dimensional axes set gets only its x and y shafts, tips, and
        labels."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["A"],
            point_ids=["Lp"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
            two_dimensional_axes_ids=["A"],
        )
        self.assertEqual(len(_get_line_actors(self.plotter)), 2)
        self.assertEqual(len(_get_tip_fill_actors(self.plotter)), 2)
        self.assertEqual(
            set(_get_labels(self.plotter)),
            {"arrow label AX", "arrow label AY", "point label Lp"},
        )

    def test_merges_arrows_one_basis_direction_at_a_time(self) -> None:
        """Test that two axes sets sharing a point and turned 90 degrees about their
        shared x basis direction merge only the arrows that coincide.

        The second axes set's x arrow coincides with the first one's x arrow, and its y
        arrow coincides with the first one's z arrow, so only its z arrow is drawn on
        its own. A merged arrow keeps the color of the first arrow merged into it, so
        the shared z and y arrow is blue.
        """
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G", "Wn"],
            point_ids=["Cg", "Ler"],
            transformations=[
                np.eye(4, dtype=float),
                output_rendering_fixtures.make_x_turned_axes_transformation_fixture(),
            ],
            axes_scale=1.0,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 4)
        self.assertEqual(
            [actor.prop.color for _, actor in line_actors],
            [pv.Color("red"), pv.Color("green"), pv.Color("blue"), pv.Color("blue")],
        )
        self.assertEqual(
            set(_get_labels(self.plotter)),
            {
                "arrow label GX/WnX",
                "arrow label GY",
                "arrow label GZ/WnY",
                "arrow label WnZ",
                "point label Cg/Ler",
            },
        )

    def test_does_not_repeat_the_ids_of_a_repeated_axes_set(self) -> None:
        """Test that an axes set passed twice is labeled with its IDs only once."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G", "G"],
            point_ids=["Cg", "Cg"],
            transformations=[np.eye(4, dtype=float), np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        self.assertEqual(
            set(_get_labels(self.plotter)),
            {"arrow label GX", "arrow label GY", "arrow label GZ", "point label Cg"},
        )

    def test_anchors_each_arrow_label_beyond_its_tip(self) -> None:
        """Test that each arrow's label is anchored along its arrow, 0.15 times the
        arrows' length beyond its tip."""
        T_pas_G_Cg_to_D_Do = (
            output_rendering_fixtures.make_offset_axes_transformation_fixture()
        )
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[T_pas_G_Cg_to_D_Do],
            axes_scale=axes_scale,
        )
        labels = _get_labels(self.plotter)
        label_distance = (1.0 + _output_rendering._AXES_ARROW_LABEL_OFFSET) * axes_scale
        for component_id, component_letter in enumerate(["X", "Y", "Z"]):
            with self.subTest(component_letter=component_letter):
                npt.assert_allclose(
                    labels[f"arrow label G{component_letter}"].position,
                    T_pas_G_Cg_to_D_Do[:3, 3]
                    + label_distance * T_pas_G_Cg_to_D_Do[:3, component_id],
                    atol=1e-12,
                )

    def test_offsets_each_point_label_away_from_its_arrows(self) -> None:
        """Test that a point's label is anchored along the negative sum of its axes
        set's basis directions, 0.3 times the arrows' length from the point, in the
        octant none of the arrows enter."""
        T_pas_G_Cg_to_D_Do = (
            output_rendering_fixtures.make_offset_axes_transformation_fixture()
        )
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[T_pas_G_Cg_to_D_Do],
            axes_scale=axes_scale,
        )
        labelDirection_D = (
            T_pas_G_Cg_to_D_Do[:3, :3]
            @ np.array([-1.0, -1.0, -1.0], dtype=float)
            / math.sqrt(3.0)
        )
        npt.assert_allclose(
            _get_labels(self.plotter)["point label Cg"].position,
            T_pas_G_Cg_to_D_Do[:3, 3]
            + _output_rendering._AXES_POINT_LABEL_OFFSET
            * axes_scale
            * labelDirection_D,
            atol=1e-12,
        )

    def test_marks_the_points_with_axes_with_dots(self) -> None:
        """Test that every point with axes is marked with a black dot, all of them in
        one mesh, drawn as spheres."""
        T_pas_Wn_Ler_to_D_Do = (
            output_rendering_fixtures.make_offset_axes_transformation_fixture()
        )
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G", "Wn"],
            point_ids=["Cg", "Ler"],
            transformations=[np.eye(4, dtype=float), T_pas_Wn_Ler_to_D_Do],
            axes_scale=1.0,
        )
        dot_actors = _get_dot_actors(self.plotter)
        self.assertEqual(len(dot_actors), 1)
        dot_mesh, dot_actor = dot_actors[0]
        npt.assert_allclose(
            dot_mesh.points,
            [np.zeros(3, dtype=float), T_pas_Wn_Ler_to_D_Do[:3, 3]],
            atol=1e-6,
        )
        self.assertEqual(dot_actor.prop.color, pv.Color("black"))
        self.assertEqual(dot_actor.prop.point_size, _output_rendering._AXES_POINT_SIZE)
        self.assertTrue(dot_actor.prop.render_points_as_spheres)

    def test_marks_each_extra_point_with_a_cross(self) -> None:
        """Test that an extra point is marked with a black cross, without arrows or a
        dot, whose two arms are centered on it and lie along its cross directions.

        Each arm is 0.1 times the arrows' length, so with arrows 2.0 meters long its
        ends sit 0.1 meters to either side of the point. The first arm's two ends come
        first, each starting from its negative side.
        """
        Lp_D_Do = np.array([1.0, 2.0, 3.0], dtype=float)
        crossDirections_D = output_rendering_fixtures.make_cross_directions_fixture()
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=[],
            point_ids=[],
            transformations=[],
            axes_scale=axes_scale,
            extra_point_ids=["Lp"],
            listExtraPoints_D_Do=[Lp_D_Do],
            listExtraPointBasisDirections_D=[np.eye(3, dtype=float)],
            listExtraPointCrossDirections_D=[crossDirections_D],
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 1)
        cross_mesh, cross_actor = line_actors[0]
        cross_half_length = 0.5 * _output_rendering._AXES_CROSS_SIZE * axes_scale
        npt.assert_allclose(
            cross_mesh.points,
            [
                Lp_D_Do - cross_half_length * crossDirections_D[0],
                Lp_D_Do + cross_half_length * crossDirections_D[0],
                Lp_D_Do - cross_half_length * crossDirections_D[1],
                Lp_D_Do + cross_half_length * crossDirections_D[1],
            ],
            atol=1e-6,
        )
        npt.assert_array_equal(cross_mesh.lines, [2, 0, 1, 2, 2, 3])
        self.assertEqual(cross_actor.prop.color, pv.Color("black"))
        self.assertEqual(
            cross_actor.prop.line_width, _output_rendering._AXES_LINE_WIDTH
        )
        self.assertEqual(len(_get_dot_actors(self.plotter)), 0)
        self.assertEqual(len(_get_tip_fill_actors(self.plotter)), 0)

    def test_offsets_each_extra_point_label_in_its_own_axes(self) -> None:
        """Test that an extra point's label is anchored along the unit vector (-1.0,
        1.0, 1.0) / sqrt(3.0) in the axes its basis directions give, 0.15 times the
        arrows' length from the point."""
        Lp_D_Do = np.array([1.0, 2.0, 3.0], dtype=float)
        extraPointBasisDirections_D = (
            output_rendering_fixtures.make_offset_axes_transformation_fixture()[:3, :3]
        )
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=[],
            point_ids=[],
            transformations=[],
            axes_scale=axes_scale,
            extra_point_ids=["Lp"],
            listExtraPoints_D_Do=[Lp_D_Do],
            listExtraPointBasisDirections_D=[extraPointBasisDirections_D],
            listExtraPointCrossDirections_D=[
                output_rendering_fixtures.make_cross_directions_fixture()
            ],
        )
        labelDirection_D = (
            extraPointBasisDirections_D
            @ np.array([-1.0, 1.0, 1.0], dtype=float)
            / math.sqrt(3.0)
        )
        npt.assert_allclose(
            _get_labels(self.plotter)["point label Lp"].position,
            Lp_D_Do
            + _output_rendering._AXES_CROSS_LABEL_OFFSET
            * axes_scale
            * labelDirection_D,
            atol=1e-12,
        )

    def test_merges_an_extra_point_into_a_point_with_axes(self) -> None:
        """Test that an extra point at an axes set's point is marked with that point's
        dot instead of a cross, and shares a label with it, anchored where that point's
        label is."""
        T_pas_G_Cg_to_D_Do = (
            output_rendering_fixtures.make_offset_axes_transformation_fixture()
        )
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[T_pas_G_Cg_to_D_Do],
            axes_scale=1.0,
            extra_point_ids=["Lp"],
            listExtraPoints_D_Do=[T_pas_G_Cg_to_D_Do[:3, 3]],
            listExtraPointBasisDirections_D=[np.eye(3, dtype=float)],
            listExtraPointCrossDirections_D=[
                output_rendering_fixtures.make_cross_directions_fixture()
            ],
        )
        self.assertEqual(len(_get_line_actors(self.plotter)), 3)
        self.assertEqual(len(_get_dot_actors(self.plotter)), 1)
        labels = _get_labels(self.plotter)
        self.assertNotIn("point label Cg", labels)
        self.assertNotIn("point label Lp", labels)
        labelDirection_D = (
            T_pas_G_Cg_to_D_Do[:3, :3]
            @ np.array([-1.0, -1.0, -1.0], dtype=float)
            / math.sqrt(3.0)
        )
        npt.assert_allclose(
            labels["point label Cg/Lp"].position,
            T_pas_G_Cg_to_D_Do[:3, 3]
            + _output_rendering._AXES_POINT_LABEL_OFFSET * labelDirection_D,
            atol=1e-12,
        )

    def test_orients_a_shared_cross_by_the_first_extra_point(self) -> None:
        """Test that two extra points at the same position share one cross, along the
        first one's cross directions, and one label."""
        Lp_D_Do = np.array([1.0, 2.0, 3.0], dtype=float)
        crossDirections_D = output_rendering_fixtures.make_cross_directions_fixture()
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=[],
            point_ids=[],
            transformations=[],
            axes_scale=axes_scale,
            extra_point_ids=["Lp1", "Lp2"],
            listExtraPoints_D_Do=[Lp_D_Do, Lp_D_Do.copy()],
            listExtraPointBasisDirections_D=[
                np.eye(3, dtype=float),
                np.eye(3, dtype=float),
            ],
            listExtraPointCrossDirections_D=[
                crossDirections_D,
                np.eye(3, dtype=float)[:2],
            ],
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 1)
        cross_mesh, _ = line_actors[0]
        cross_half_length = 0.5 * _output_rendering._AXES_CROSS_SIZE * axes_scale
        npt.assert_allclose(
            cross_mesh.points,
            [
                Lp_D_Do - cross_half_length * crossDirections_D[0],
                Lp_D_Do + cross_half_length * crossDirections_D[0],
                Lp_D_Do - cross_half_length * crossDirections_D[1],
                Lp_D_Do + cross_half_length * crossDirections_D[1],
            ],
            atol=1e-6,
        )
        self.assertEqual(set(_get_labels(self.plotter)), {"point label Lp1/Lp2"})

    def test_leaves_the_extra_points_unlabeled_when_asked(self) -> None:
        """Test that with label_extra_points set to False, an extra point on its own
        gets a cross but no label, and an extra point at an axes set's point leaves that
        point's label alone."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
            extra_point_ids=["Lp1", "Lp2"],
            listExtraPoints_D_Do=[
                np.zeros(3, dtype=float),
                np.array([1.0, 2.0, 3.0], dtype=float),
            ],
            listExtraPointBasisDirections_D=[
                np.eye(3, dtype=float),
                np.eye(3, dtype=float),
            ],
            listExtraPointCrossDirections_D=[
                output_rendering_fixtures.make_cross_directions_fixture(),
                output_rendering_fixtures.make_cross_directions_fixture(),
            ],
            label_extra_points=False,
        )
        line_actors = _get_line_actors(self.plotter)
        self.assertEqual(len(line_actors), 4)
        cross_mesh, _ = line_actors[3]
        self.assertEqual(cross_mesh.n_points, 4)
        self.assertEqual(
            set(_get_labels(self.plotter)),
            {"arrow label GX", "arrow label GY", "arrow label GZ", "point label Cg"},
        )

    def test_styles_the_plain_labels(self) -> None:
        """Test that plain labels are black, on an opaque background matching the
        Plotter's, and set in the monospaced font at the plain label font size."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        labels = _get_labels(self.plotter)
        self.assertEqual(len(labels), 4)
        for name, label in labels.items():
            with self.subTest(name=name):
                self.assertEqual(label.input, name.split(" ")[-1])
                self.assertEqual(label.size, _output_rendering._AXES_LABEL_FONT_SIZE)
                self.assertEqual(label.prop.color, pv.Color("black"))
                self.assertEqual(
                    label.prop.background_color, self.plotter.background_color
                )
                self.assertEqual(label.prop.background_opacity, 1.0)
                self.assertEqual(label.prop.GetFontFile(), str(_fonts.MONO_FONT_PATH))

    def test_writes_the_math_labels_as_math(self) -> None:
        """Test that math labels join their merged labels' math inside one pair of
        dollar signs, and are set in the Times font family at the math label font
        size."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G", "Wn"],
            point_ids=["Cg", "Ler"],
            transformations=[np.eye(4, dtype=float), np.eye(4, dtype=float)],
            axes_scale=1.0,
            math_labels=True,
        )
        labels = _get_labels(self.plotter)
        self.assertEqual(
            labels["arrow label GX/WnX"].input,
            "$"
            + _output_rendering.get_math_axes_label("G", "X")
            + r" \,/\, "
            + _output_rendering.get_math_axes_label("Wn", "X")
            + "$",
        )
        self.assertEqual(
            labels["point label Cg/Ler"].input,
            r"$\mathrm{CG} \,/\, \mathrm{LER}$",
        )
        for name, label in labels.items():
            with self.subTest(name=name):
                self.assertEqual(
                    label.size, _output_rendering._AXES_MATH_LABEL_FONT_SIZE
                )
                self.assertEqual(label.prop.font_family, "times")

    def test_rejects_an_id_that_cannot_be_written_as_math(self) -> None:
        """Test that math labels raise a ValueError for an axes ID or a point ID that
        isn't valid."""
        cases = [(["Q"], ["Cg"]), (["G"], ["Q1"])]
        for axes_ids, point_ids in cases:
            with self.subTest(axes_ids=axes_ids, point_ids=point_ids):
                plotter = pv.Plotter(off_screen=True)
                self.addCleanup(plotter.close)
                with self.assertRaises(ValueError):
                    _output_rendering.add_axes_and_points(
                        plotter,
                        axes_ids=axes_ids,
                        point_ids=point_ids,
                        transformations=[np.eye(4, dtype=float)],
                        axes_scale=1.0,
                        math_labels=True,
                    )

    def test_builds_each_tip_as_pyvista_builds_its_cone(self) -> None:
        """Test that each tip has the same points, in the same order, as a cone PyVista
        builds for it alone, and as many faces covering the same area.

        Each tip is 0.2 times the arrows' length long, with a base radius of 0.1 times
        their length, and sits at the end of its arrow. VTK turns a cone toward a
        direction with a negative x component differently than toward one with a
        positive x component. The oblique axes set's basis directions all have negative
        x components and the aligned axes set's don't, so each mesh mixes the two ways.
        The areas only match if each tip's faces use its own points.
        """
        transformations = [
            np.eye(4, dtype=float),
            output_rendering_fixtures.make_oblique_axes_transformation_fixture(),
        ]
        axes_scale = 2.0
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G", "Wn"],
            point_ids=["Cg", "Ler"],
            transformations=transformations,
            axes_scale=axes_scale,
        )
        tip_fill_actors = _get_tip_fill_actors(self.plotter)
        self.assertEqual(len(tip_fill_actors), 3)
        tip_length = _output_rendering._AXES_TIP_LENGTH * axes_scale
        for component_id, (tip_mesh, _) in enumerate(tip_fill_actors):
            with self.subTest(component_id=component_id):
                cones = [
                    pv.Cone(
                        center=transformation[:3, 3]
                        + (axes_scale - 0.5 * tip_length)
                        * transformation[:3, component_id],
                        direction=transformation[:3, component_id],
                        height=tip_length,
                        radius=_output_rendering._AXES_TIP_RADIUS * axes_scale,
                        resolution=_output_rendering._AXES_TIP_RESOLUTION,
                    )
                    for transformation in transformations
                ]
                npt.assert_allclose(
                    tip_mesh.points,
                    np.vstack([cone.points for cone in cones]),
                    atol=1e-6,
                )
                self.assertEqual(tip_mesh.n_cells, sum(cone.n_cells for cone in cones))
                self.assertAlmostEqual(
                    tip_mesh.area, sum(cone.area for cone in cones), places=5
                )

    def test_fills_the_tips_with_the_background_pushed_back(self) -> None:
        """Test that the tips' filled faces are unlit, match the background, and are
        pushed away from the camera by the tip fill polygon offset."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        for tip_mesh, tip_fill_actor in _get_tip_fill_actors(self.plotter):
            with self.subTest(n_points=tip_mesh.n_points):
                self.assertEqual(
                    tip_fill_actor.prop.color, self.plotter.background_color
                )
                self.assertFalse(tip_fill_actor.prop.lighting)
                tip_fill_mapper = tip_fill_actor.mapper
                offset_factor = vtk_reference(0.0)
                offset_units = vtk_reference(0.0)
                tip_fill_mapper.GetRelativeCoincidentTopologyPolygonOffsetParameters(
                    offset_factor, offset_units
                )
                self.assertEqual(
                    offset_factor, _output_rendering._AXES_TIP_FILL_OFFSET_FACTOR
                )
                self.assertEqual(
                    offset_units, _output_rendering._AXES_TIP_FILL_OFFSET_UNITS
                )

    def test_turns_on_polygon_offset_only_while_rendering(self) -> None:
        """Test that VTK's coincident topology resolution mode is polygon offset while
        the Plotter renders, and is restored to its previous mode afterward.

        The mode is shared by every mapper in the process, so it is first set to a mode
        other than its default, which shows that the previous mode is restored rather
        than the default, and it is reset when the test ends. The recording observer is
        added with a lower priority than the default, so it runs after the observer that
        turns polygon offset on.
        """
        resolve_mapper = pv.DataSetMapper()
        self.addCleanup(setattr, resolve_mapper, "resolve", resolve_mapper.resolve)
        resolve_mapper.resolve = "shift_zbuffer"
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        resolve_modes_while_rendering: list[str] = []

        def record_resolve_mode(caller: object, event: str) -> None:
            """Records VTK's coincident topology resolution mode.

            :param caller: The object that invoked the event, which is unused.
            :param event: The name of the event, which is unused.
            :return: None
            """
            resolve_modes_while_rendering.append(pv.DataSetMapper().resolve)

        self.plotter.renderer.AddObserver("StartEvent", record_resolve_mode, -1.0)
        self.plotter.screenshot(return_img=True)
        self.assertEqual(set(resolve_modes_while_rendering), {"polygon_offset"})
        self.assertEqual(resolve_mapper.resolve, "shift_zbuffer")

    def test_justifies_the_labels_away_from_what_they_label(self) -> None:
        """Test that rendering justifies each label so its text extends away from what
        it labels on screen.

        The camera looks down the negative z direction with parallel projection, so the
        x arrow's label extends right, the y arrow's label extends up, and the z arrow's
        label, whose anchor lies straight in front of the point, is centered. The
        point's label lies along (-1.0, -1.0, -1.0), so it extends left and down.
        """
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        self.plotter.camera.parallel_projection = True
        self.plotter.camera_position = "xy"
        self.plotter.screenshot(return_img=True)
        labels = _get_labels(self.plotter)
        expected_justifications = {
            "arrow label GX": ("left", "center"),
            "arrow label GY": ("center", "bottom"),
            "arrow label GZ": ("center", "center"),
            "point label Cg": ("right", "top"),
        }
        for name, expected_justification in expected_justifications.items():
            with self.subTest(name=name):
                self.assertEqual(
                    (
                        labels[name].prop.justification_horizontal,
                        labels[name].prop.justification_vertical,
                    ),
                    expected_justification,
                )

    def test_justifies_the_labels_of_reversed_arrows(self) -> None:
        """Test that the labels of arrows pointing left and down on screen extend left
        and down, rather than back over their arrows.

        The axes set is turned 180 degrees about the z axis, so its point's label lies
        along (1.0, 1.0, -1.0) and extends right and up.
        """
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[
                output_rendering_fixtures.make_reversed_axes_transformation_fixture()
            ],
            axes_scale=1.0,
        )
        self.plotter.camera.parallel_projection = True
        self.plotter.camera_position = "xy"
        self.plotter.screenshot(return_img=True)
        labels = _get_labels(self.plotter)
        expected_justifications = {
            "arrow label GX": ("right", "center"),
            "arrow label GY": ("center", "top"),
            "point label Cg": ("left", "bottom"),
        }
        for name, expected_justification in expected_justifications.items():
            with self.subTest(name=name):
                self.assertEqual(
                    (
                        labels[name].prop.justification_horizontal,
                        labels[name].prop.justification_vertical,
                    ),
                    expected_justification,
                )

    def test_dragging_a_label_moves_it_with_the_mouse(self) -> None:
        """Test that a label dragged with the left mouse button follows the mouse on
        screen, keeping its depth, while the camera stays still."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        self.plotter.camera_position = "xy"
        self.plotter.screenshot(return_img=True)
        label = _get_labels(self.plotter)["arrow label GX"]
        camera_position = self.plotter.camera.position
        anchor_display = _get_display_point(
            self.plotter, np.array(label.position, dtype=float)
        )

        # Press just inside the label, which extends right from its anchor, and drag it
        # 50 pixels to the right.
        start_display = (round(anchor_display[0]) + 2, round(anchor_display[1]) + 2)
        _drag_mouse(
            self.plotter, start_display, (start_display[0] + 50, start_display[1])
        )
        npt.assert_allclose(
            _get_display_point(self.plotter, np.array(label.position, dtype=float)),
            anchor_display + np.array([50.0, 0.0, 0.0], dtype=float),
            atol=1e-6,
        )
        self.assertEqual(self.plotter.camera.position, camera_position)

    def test_a_dragged_label_returns_when_the_camera_moves(self) -> None:
        """Test that a dragged label keeps its position through renders until the camera
        moves, and then returns to its anchor."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        self.plotter.camera_position = "xy"
        self.plotter.screenshot(return_img=True)
        label = _get_labels(self.plotter)["arrow label GX"]
        anchor_display = _get_display_point(
            self.plotter, np.array(label.position, dtype=float)
        )
        start_display = (round(anchor_display[0]) + 2, round(anchor_display[1]) + 2)
        _drag_mouse(
            self.plotter, start_display, (start_display[0] + 50, start_display[1])
        )
        draggedLabel_D_Do = np.array(label.position, dtype=float)

        self.plotter.screenshot(return_img=True)
        npt.assert_array_equal(label.position, draggedLabel_D_Do)

        self.plotter.camera_position = "xz"
        self.plotter.screenshot(return_img=True)
        npt.assert_allclose(
            label.position,
            [1.0 + _output_rendering._AXES_ARROW_LABEL_OFFSET, 0.0, 0.0],
            atol=1e-12,
        )

    def test_dragging_away_from_the_labels_rotates_the_camera(self) -> None:
        """Test that a left button drag that doesn't start on a label rotates the
        camera, as it would without the labels, and moves none of them."""
        _output_rendering.add_axes_and_points(
            self.plotter,
            axes_ids=["G"],
            point_ids=["Cg"],
            transformations=[np.eye(4, dtype=float)],
            axes_scale=1.0,
        )
        self.plotter.camera_position = "xy"
        self.plotter.screenshot(return_img=True)
        labels = _get_labels(self.plotter)
        listAnchors_D_Do = [
            np.array(label.position, dtype=float) for label in labels.values()
        ]
        camera_position = self.plotter.camera.position

        # The window's bottom left corner is far from every label.
        _drag_mouse(self.plotter, (1, 1), (40, 40))
        self.assertFalse(np.allclose(self.plotter.camera.position, camera_position))
        for label, anchor_D_Do in zip(labels.values(), listAnchors_D_Do):
            npt.assert_array_equal(label.position, anchor_D_Do)
