"""This module contains classes to test the convergence analysis functions."""

import unittest

import numpy as np

import pterasoftware as ps
from pterasoftware import convergence
from tests.unit.fixtures import (
    geometry_fixtures,
    movement_fixtures,
    problem_fixtures,
)


class TestConvergedParameterId(unittest.TestCase):
    """This class contains methods for testing
    convergence._converged_parameter_id.
    """

    def test_single_returns_this_id(self) -> None:
        """Test that a single tested value returns this iteration's own index."""
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=0, single=True, converged=False
            ),
            0,
        )

    def test_single_takes_precedence_over_converged(self) -> None:
        """Test that a single tested value returns this iteration's own index even
        when the converged flag is set.
        """
        self.assertEqual(
            convergence._converged_parameter_id(this_id=3, single=True, converged=True),
            3,
        )

    def test_converged_returns_coarser_id(self) -> None:
        """Test that a converged iteration returns the incrementally coarser
        index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=4, single=False, converged=True
            ),
            3,
        )

    def test_saturated_returns_this_id(self) -> None:
        """Test that an iteration that passed without converging returns this
        iteration's own index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=5, single=False, converged=False
            ),
            5,
        )

    def test_converged_id_depends_on_this_id(self) -> None:
        """Test that the coarser index tracks this iteration's index rather than a
        fixed value.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=7, single=False, converged=True
            ),
            6,
        )

    def test_returns_int(self) -> None:
        """Test that the converged index is returned as an int."""
        result = convergence._converged_parameter_id(
            this_id=2, single=False, converged=True
        )
        self.assertIsInstance(result, int)


class TestValidateCoefficientMask(unittest.TestCase):
    """This class contains methods for testing
    convergence._validate_coefficient_mask.
    """

    def test_none_returns_all_true(self) -> None:
        """Test that None returns a (6,) mask of all True."""
        result = convergence._validate_coefficient_mask(None)
        self.assertTrue(np.array_equal(result, np.ones(6, dtype=bool)))

    def test_tuple_is_returned_as_bool_array(self) -> None:
        """Test that a valid tuple is returned as an equivalent (6,) bool ndarray."""
        result = convergence._validate_coefficient_mask(
            (True, False, True, False, True, False)
        )
        self.assertTrue(
            np.array_equal(result, np.array([1, 0, 1, 0, 1, 0], dtype=bool))
        )

    def test_non_tuple_raises_type_error(self) -> None:
        """Test that a non-tuple, non-None mask raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_coefficient_mask([True] * 6)  # type: ignore[arg-type]

    def test_wrong_length_raises_value_error(self) -> None:
        """Test that a mask without exactly six elements raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_coefficient_mask(
                (True, True, True)  # type: ignore[arg-type]
            )

    def test_non_bool_element_raises_type_error(self) -> None:
        """Test that a mask with a non-bool element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_coefficient_mask(
                (True, 1, True, True, True, True)  # type: ignore[arg-type]
            )

    def test_all_false_raises_value_error(self) -> None:
        """Test that a mask with no True element raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_coefficient_mask((False,) * 6)


class TestCheckCoefficientConvergence(unittest.TestCase):
    """This class contains methods for testing
    convergence._check_coefficient_convergence.
    """

    def setUp(self) -> None:
        """Set up a full mask and the tolerances used across the tests.

        :return: None
        """
        self.mask = np.ones(6, dtype=bool)
        self.rtol = 0.05
        self.atol = 0.001

    def test_identical_coefficients_converge(self) -> None:
        """Test that identical coefficients converge with a perfect metric."""
        these = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        converged, metric, _ = convergence._check_coefficient_convergence(
            these, these.copy(), self.rtol, self.atol, self.mask
        )
        self.assertTrue(converged)
        self.assertEqual(metric, 100.0)

    def test_large_relative_change_does_not_converge(self) -> None:
        """Test that a coefficient changing by more than the relative tolerance does not
        converge and is reported as the limiting coefficient.
        """
        coarser = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        these = coarser.copy()
        these[0, 2] = 2.5
        converged, _, limiting_id = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertFalse(converged)
        self.assertEqual(limiting_id, 2)

    def test_masked_out_coefficient_is_ignored(self) -> None:
        """Test that masking out the only offending coefficient makes the check
        converge.
        """
        coarser = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        these = coarser.copy()
        these[0, 2] = 2.5
        mask = np.array([True, True, False, True, True, True], dtype=bool)
        converged, _, _ = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, mask
        )
        self.assertTrue(converged)

    def test_absolute_tolerance_floors_near_zero(self) -> None:
        """Test that a coefficient near zero converges via the absolute tolerance floor
        even though its relative change is large.
        """
        these = np.zeros((1, 6), dtype=float)
        coarser = np.zeros((1, 6), dtype=float)
        coarser[0, 0] = 0.5 * self.atol
        converged, _, _ = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertTrue(converged)

    def test_all_airplanes_must_converge(self) -> None:
        """Test that the check fails when any one Airplane has an unconverged
        coefficient.
        """
        coarser = np.array(
            [
                [1.0, 0.0, 2.0, 0.1, 0.0, 0.05],
                [1.0, 0.0, 2.0, 0.1, 0.0, 0.05],
            ]
        )
        these = coarser.copy()
        these[1, 0] = 2.0
        converged, _, limiting_id = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertFalse(converged)
        self.assertEqual(limiting_id, 0)

    def test_returns_bool_float_int(self) -> None:
        """Test that the result is a bool, a float, and an int."""
        these = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        converged, metric, limiting_id = convergence._check_coefficient_convergence(
            these, these.copy(), self.rtol, self.atol, self.mask
        )
        self.assertIsInstance(converged, bool)
        self.assertIsInstance(metric, float)
        self.assertIsInstance(limiting_id, int)


class TestValidatePanelAspectRatioBounds(unittest.TestCase):
    """This class contains methods for testing
    convergence._validate_panel_aspect_ratio_bounds.
    """

    def test_valid_descending_bounds_pass(self) -> None:
        """Test that a valid descending tuple of ints is accepted without raising."""
        # A valid descending tuple of ints does not raise.
        convergence._validate_panel_aspect_ratio_bounds((4, 1))

    def test_equal_bounds_pass(self) -> None:
        """Test that equal bounds, a single Panel aspect ratio, are accepted without
        raising.
        """
        # Equal bounds, a single Panel aspect ratio, do not raise.
        convergence._validate_panel_aspect_ratio_bounds((2, 2))

    def test_non_tuple_raises_type_error(self) -> None:
        """Test that a non-tuple of bounds raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_panel_aspect_ratio_bounds([4, 1])  # type: ignore[arg-type]

    def test_wrong_length_raises_type_error(self) -> None:
        """Test that a tuple without exactly two elements raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_panel_aspect_ratio_bounds((4, 2, 1))  # type: ignore[arg-type]

    def test_non_int_element_raises_type_error(self) -> None:
        """Test that a tuple with a non-int element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_panel_aspect_ratio_bounds((4.0, 1.0))  # type: ignore[arg-type]

    def test_ascending_bounds_raise_value_error(self) -> None:
        """Test that a tuple whose first value is less than its second raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            convergence._validate_panel_aspect_ratio_bounds((1, 4))

    def test_non_positive_second_value_raises_value_error(self) -> None:
        """Test that a second value at or below zero raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_panel_aspect_ratio_bounds((4, 0))


class TestValidateNumChordwisePanelsBounds(unittest.TestCase):
    """This class contains methods for testing
    convergence._validate_num_chordwise_panels_bounds.
    """

    def test_valid_ascending_bounds_pass(self) -> None:
        """Test that a valid ascending tuple of ints is accepted without raising."""
        # A valid ascending tuple of ints does not raise.
        convergence._validate_num_chordwise_panels_bounds((3, 12))

    def test_equal_bounds_pass(self) -> None:
        """Test that equal bounds, a single number of chordwise Panels, are accepted
        without raising.
        """
        # Equal bounds, a single number of chordwise Panels, do not raise.
        convergence._validate_num_chordwise_panels_bounds((5, 5))

    def test_non_tuple_raises_type_error(self) -> None:
        """Test that a non-tuple of bounds raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_num_chordwise_panels_bounds([3, 12])  # type: ignore[arg-type]

    def test_wrong_length_raises_type_error(self) -> None:
        """Test that a tuple without exactly two elements raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_num_chordwise_panels_bounds((3, 6, 12))  # type: ignore[arg-type]

    def test_non_int_element_raises_type_error(self) -> None:
        """Test that a tuple with a non-int element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_num_chordwise_panels_bounds((3.0, 12.0))  # type: ignore[arg-type]

    def test_descending_bounds_raise_value_error(self) -> None:
        """Test that a tuple whose second value is less than its first raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            convergence._validate_num_chordwise_panels_bounds((12, 3))

    def test_non_positive_first_value_raises_value_error(self) -> None:
        """Test that a first value at or below zero raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_num_chordwise_panels_bounds((0, 12))


class TestRejectUnrefinableWings(unittest.TestCase):
    """This class contains methods for testing
    convergence._reject_unrefinable_wings.
    """

    @staticmethod
    def _make_edge_defined_airplane() -> ps.geometry.airplane.Airplane:
        """Builds an Airplane holding a single edge-defined Wing.

        The Wing is built with Wing.from_edge_points, so its spanwise mesh marker is
        "edge_defined", which the convergence functions can refine.
        """
        ys = np.linspace(0.0, 2.0, 20)
        zeros = np.zeros_like(ys)
        leadingEdgePoints_Wn_Ler = np.column_stack((0.25 * ys, ys, zeros))
        trailingEdgePoints_Wn_Ler = np.column_stack((np.ones_like(ys), ys, zeros))
        return ps.geometry.airplane.Airplane(
            wings=[
                ps.geometry.wing.Wing.from_edge_points(
                    leadingEdgePoints_Wn_Ler=leadingEdgePoints_Wn_Ler,
                    trailingEdgePoints_Wn_Ler=trailingEdgePoints_Wn_Ler,
                    num_wing_cross_sections=5,
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                    name="Edge Wing",
                    symmetric=False,
                    num_chordwise_panels=4,
                )
            ],
            name="Edge Defined Airplane",
        )

    @staticmethod
    def _make_exploded_airplane() -> ps.geometry.airplane.Airplane:
        """Builds an Airplane holding a single exploded Wing named "Exploded Wing".

        The Wing is built with explode_into_strips=True, so its spanwise mesh marker is
        "exploded", which the convergence functions cannot refine.
        """
        return ps.geometry.airplane.Airplane(
            wings=[
                ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        ps.geometry.wing_cross_section.WingCrossSection(
                            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                            num_spanwise_panels=1,
                            chord=1.0,
                            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                            spanwise_spacing="uniform",
                        ),
                        ps.geometry.wing_cross_section.WingCrossSection(
                            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                            num_spanwise_panels=None,
                            chord=1.0,
                            Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                            spanwise_spacing=None,
                        ),
                    ],
                    name="Exploded Wing",
                    explode_into_strips=True,
                    num_chordwise_panels=2,
                    chordwise_spacing="uniform",
                )
            ],
            name="Exploded Airplane",
        )

    def test_trapezoidal_wing_is_accepted(self) -> None:
        """Test that a trapezoidal Wing, built from WingCrossSections, is accepted."""
        trapezoidal_airplane = geometry_fixtures.make_first_airplane_fixture()
        # A trapezoidal Wing does not raise.
        convergence._reject_unrefinable_wings(
            (trapezoidal_airplane,), "analyze_steady_convergence"
        )

    def test_edge_defined_wing_is_accepted(self) -> None:
        """Test that an edge-defined Wing, built from edge curves, is accepted."""
        # An edge-defined Wing does not raise.
        convergence._reject_unrefinable_wings(
            (self._make_edge_defined_airplane(),), "analyze_steady_convergence"
        )

    def test_exploded_wing_raises_value_error(self) -> None:
        """Test that an exploded Wing, which carries no edge curves, raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            convergence._reject_unrefinable_wings(
                (self._make_exploded_airplane(),), "analyze_steady_convergence"
            )

    def test_error_message_names_wing_and_function(self) -> None:
        """Test that the error message names the offending Wing and the calling
        function.
        """
        with self.assertRaisesRegex(ValueError, "Exploded Wing"):
            convergence._reject_unrefinable_wings(
                (self._make_exploded_airplane(),), "analyze_steady_convergence"
            )
        with self.assertRaisesRegex(ValueError, "analyze_steady_convergence"):
            convergence._reject_unrefinable_wings(
                (self._make_exploded_airplane(),), "analyze_steady_convergence"
            )

    def test_rejects_when_any_airplane_has_unrefinable_wing(self) -> None:
        """Test that a refinable Airplane paired with an unrefinable one still raises."""
        trapezoidal_airplane = geometry_fixtures.make_first_airplane_fixture()
        with self.assertRaises(ValueError):
            convergence._reject_unrefinable_wings(
                (trapezoidal_airplane, self._make_exploded_airplane()),
                "analyze_unsteady_convergence",
            )

    def test_empty_airplanes_is_accepted(self) -> None:
        """Test that an empty tuple of Airplanes is vacuously accepted."""
        # An empty tuple of Airplanes does not raise.
        convergence._reject_unrefinable_wings((), "analyze_steady_convergence")


class TestAnalyzeSteadyConvergenceValidation(unittest.TestCase):
    """This class contains methods for testing the input validation of
    convergence.analyze_steady_convergence that raises before any solving.
    """

    def setUp(self) -> None:
        """Set up a valid SteadyProblem shared by the tests.

        :return: None
        """
        self.steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

    def test_non_steady_problem_raises_type_error(self) -> None:
        """Test that a ref_problem that is not a SteadyProblem raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence.analyze_steady_convergence(
                ref_problem=problem_fixtures.make_basic_unsteady_problem_fixture(),
                solver_type="steady ring vortex lattice method",
            )

    def test_invalid_solver_type_raises_value_error(self) -> None:
        """Test that an unrecognized solver_type raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence.analyze_steady_convergence(
                ref_problem=self.steady_problem,
                solver_type="unsteady ring vortex lattice method",
            )


class TestAnalyzeUnsteadyConvergenceValidation(unittest.TestCase):
    """This class contains methods for testing the input validation of
    convergence.analyze_unsteady_convergence that raises before any solving.
    """

    def setUp(self) -> None:
        """Set up a static and a variable UnsteadyProblem shared by the tests.

        :return: None
        """
        self.variable_problem = problem_fixtures.make_basic_unsteady_problem_fixture()
        self.static_problem = ps.problems.UnsteadyProblem(
            movement=movement_fixtures.make_static_movement_fixture()
        )

    def test_non_unsteady_problem_raises_type_error(self) -> None:
        """Test that a ref_problem that is not an UnsteadyProblem raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence.analyze_unsteady_convergence(
                ref_problem=problem_fixtures.make_basic_steady_problem_fixture(),
                num_cycles_bounds=(1, 2),
            )

    def test_no_wake_state_raises_value_error(self) -> None:
        """Test that setting both prescribed_wake and free_wake to False raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.variable_problem,
                prescribed_wake=False,
                free_wake=False,
            )

    def test_static_geometry_rejects_num_cycles_bounds(self) -> None:
        """Test that supplying num_cycles_bounds for a static-geometry problem raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.static_problem,
                num_cycles_bounds=(1, 2),
            )

    def test_static_geometry_requires_num_chords_bounds(self) -> None:
        """Test that omitting num_chords_bounds for a static-geometry problem raises a
        TypeError.
        """
        with self.assertRaises(TypeError):
            convergence.analyze_unsteady_convergence(ref_problem=self.static_problem)

    def test_num_chords_bounds_wrong_length_raises_type_error(self) -> None:
        """Test that a num_chords_bounds without exactly two elements raises a
        TypeError.
        """
        with self.assertRaises(TypeError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.static_problem,
                num_chords_bounds=(1, 2, 3),  # type: ignore[arg-type]
            )

    def test_num_chords_bounds_non_int_raises_type_error(self) -> None:
        """Test that a num_chords_bounds with a non-int element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.static_problem,
                num_chords_bounds=(1.0, 4.0),  # type: ignore[arg-type]
            )

    def test_descending_num_chords_bounds_raises_value_error(self) -> None:
        """Test that a num_chords_bounds in descending order raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.static_problem,
                num_chords_bounds=(4, 1),
            )

    def test_non_positive_num_chords_bounds_raises_value_error(self) -> None:
        """Test that a num_chords_bounds at or below zero raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.static_problem,
                num_chords_bounds=(0, 0),
            )

    def test_variable_geometry_rejects_num_chords_bounds(self) -> None:
        """Test that supplying num_chords_bounds for a variable-geometry problem raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.variable_problem,
                num_chords_bounds=(1, 4),
            )

    def test_variable_geometry_requires_num_cycles_bounds(self) -> None:
        """Test that omitting num_cycles_bounds for a variable-geometry problem raises a
        TypeError.
        """
        with self.assertRaises(TypeError):
            convergence.analyze_unsteady_convergence(ref_problem=self.variable_problem)

    def test_num_cycles_bounds_non_int_raises_type_error(self) -> None:
        """Test that a num_cycles_bounds with a non-int element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.variable_problem,
                num_cycles_bounds=(1.0, 2.0),  # type: ignore[arg-type]
            )

    def test_descending_num_cycles_bounds_raises_value_error(self) -> None:
        """Test that a num_cycles_bounds in descending order raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.variable_problem,
                num_cycles_bounds=(2, 1),
            )

    def test_non_positive_num_cycles_bounds_raises_value_error(self) -> None:
        """Test that a num_cycles_bounds at or below zero raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence.analyze_unsteady_convergence(
                ref_problem=self.variable_problem,
                num_cycles_bounds=(0, 0),
            )
