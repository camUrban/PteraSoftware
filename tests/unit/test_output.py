"""Unit tests for camera parameter validation and camera fitting."""

import unittest
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt

from pterasoftware import _output


class TestValidateCameraInputs(unittest.TestCase):
    """Tests validation of the shared camera arguments."""

    def test_valid_camera_inputs_are_normalized_to_float_arrays(self) -> None:
        position, focal_point, up, scale, return_state = _output.validate_camera_inputs(
            [3, 2, 1], [0, 0, 0], [0, 0, 1], 4, True
        )
        npt.assert_array_equal(position, [3.0, 2.0, 1.0])
        npt.assert_array_equal(focal_point, [0.0, 0.0, 0.0])
        npt.assert_array_equal(up, [0.0, 0.0, 1.0])
        self.assertEqual(scale, 4.0)
        self.assertIs(return_state, True)

    def test_invalid_vector_shape_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "camera_position"):
            _output.validate_camera_inputs([1, 2], None, None, None, False)

    def test_non_positive_parallel_scale_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "camera_parallel_scale"):
            _output.validate_camera_inputs(None, None, None, 0, False)

    def test_non_boolean_return_flag_is_rejected(self) -> None:
        with self.assertRaisesRegex(TypeError, "return_camera_state"):
            _output.validate_camera_inputs(None, None, None, None, 1)


class TestApplyCameraParameters(unittest.TestCase):
    """Tests geometric camera validation and scale fitting without rendering."""

    def setUp(self) -> None:
        self.plotter = SimpleNamespace(
            camera=SimpleNamespace(
                position=(3.0, -3.0, 3.0),
                focal_point=(0.0, 0.0, 0.0),
                up=(0.0, 0.0, 1.0),
                parallel_scale=1.0,
            ),
            renderer=SimpleNamespace(camera_set=False),
            reset_camera_clipping_range=lambda: None,
        )
        self.bounds = (-1.0, 1.0, -2.0, 2.0, -3.0, 3.0)

    def _apply(
        self,
        position: np.ndarray | None = None,
        focal_point: np.ndarray | None = None,
        up: np.ndarray | None = None,
        scale: float | None = None,
    ) -> None:
        _output.apply_camera_parameters(
            self.plotter,
            position,
            focal_point,
            up,
            scale,
            self.bounds,
            aspect_ratio=4 / 3,
        )

    def test_automatic_scale_is_positive_and_camera_values_are_applied(self) -> None:
        position = np.array([4.0, -3.0, 2.0])
        focal_point = np.array([0.0, 0.0, 0.0])
        up = np.array([0.0, 0.0, 1.0])
        self._apply(position, focal_point, up)
        npt.assert_allclose(self.plotter.camera.position, position)
        npt.assert_allclose(self.plotter.camera.focal_point, focal_point)
        npt.assert_allclose(self.plotter.camera.up, up)
        self.assertAlmostEqual(
            self.plotter.camera.parallel_scale,
            4.057443139169,
            places=10,
        )

    def test_explicit_parallel_scale_is_preserved(self) -> None:
        self._apply(scale=7.25)
        self.assertEqual(self.plotter.camera.parallel_scale, 7.25)

    def test_equal_position_and_focal_point_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not be equal"):
            self._apply(
                position=np.array([1.0, 2.0, 3.0]),
                focal_point=np.array([1.0, 2.0, 3.0]),
            )

    def test_zero_up_vector_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-zero length"):
            self._apply(up=np.array([0.0, 0.0, 0.0]))

    def test_small_nonzero_up_vector_is_accepted(self) -> None:
        # The magnitude of an up vector does not change its direction.
        self._apply(up=np.array([0.0, 0.0, 1e-10]))
        npt.assert_allclose(self.plotter.camera.up, [0.0, 0.0, 1e-10])

    def test_small_nonzero_view_direction_is_accepted(self) -> None:
        # A small but nonzero position-to-focal distance is still a direction.
        self._apply(
            position=np.array([0.0, 0.0, 1e-10]),
            focal_point=np.array([0.0, 0.0, 0.0]),
            up=np.array([0.0, 1.0, 0.0]),
        )
        npt.assert_allclose(self.plotter.camera.position, [0.0, 0.0, 1e-10])

    def test_partial_override_keeps_unspecified_camera_vectors(self) -> None:
        original_focal_point = self.plotter.camera.focal_point
        original_up = self.plotter.camera.up
        self._apply(position=np.array([5.0, 4.0, 3.0]))
        npt.assert_allclose(self.plotter.camera.position, [5.0, 4.0, 3.0])
        npt.assert_allclose(self.plotter.camera.focal_point, original_focal_point)
        npt.assert_allclose(self.plotter.camera.up, original_up)
        self.assertGreater(self.plotter.camera.parallel_scale, 0.0)

    def test_non_finite_camera_position_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "finite view direction"):
            self._apply(position=np.array([np.inf, 0.0, 1.0]))

    def test_non_finite_up_vector_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "finite, non-zero length"):
            self._apply(up=np.array([0.0, np.nan, 1.0]))

    def test_antiparallel_up_vector_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "parallel or antiparallel"):
            self._apply(
                position=np.array([0.0, 0.0, 1.0]),
                focal_point=np.array([0.0, 0.0, 0.0]),
                up=np.array([0.0, 0.0, -1.0]),
            )

    def test_parallel_up_vector_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not be parallel"):
            self._apply(
                position=np.array([0.0, 0.0, 1.0]),
                focal_point=np.array([0.0, 0.0, 0.0]),
                up=np.array([0.0, 0.0, 2.0]),
            )
