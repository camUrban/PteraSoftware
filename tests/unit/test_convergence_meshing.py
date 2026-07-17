"""This module contains classes to test the convergence meshing functions."""

import inspect
import unittest
from typing import Any

import numpy as np

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _convergence_meshing
from tests.unit.fixtures import geometry_fixtures, operating_point_fixtures

# Both convergence builders rebuild each Airplane, Wing, and WingCrossSection by naming
# its parameters one by one, so these sets partition each class's public constructor
# parameters into the ones a build copies from its reference, the ones a build
# deliberately changes, and the ones a build deliberately omits.
# TestGeometryCopyParameterCoverage asserts that each partition covers every public
# parameter, so a parameter added to one of these classes fails there until it is
# classified here and passed at every site that constructs the class. Airfoil needs no
# partition, because a build shares each reference Airfoil rather than rebuilding one.
_AIRPLANE_COPIED_PARAMETERS = frozenset({"name", "Cg_GP1_CgP1", "weight"})

# A build gives each Airplane copy this iteration's Wings, and passes None for the
# reference dimensions so that the copy recalculates them for its own refined geometry.
_AIRPLANE_CHANGED_PARAMETERS = frozenset({"wings", "s_ref", "c_ref", "b_ref"})
_AIRPLANE_OMITTED_PARAMETERS: frozenset[str] = frozenset()

_WING_COPIED_PARAMETERS = frozenset(
    {
        "name",
        "Ler_Gs_Cgs",
        "angles_Gs_to_Wn_ixyz",
        "symmetric",
        "mirror_only",
        "symmetryNormal_G",
        "symmetryPoint_G_Cg",
        "chordwise_spacing",
    }
)
_WING_CHANGED_PARAMETERS = frozenset({"wing_cross_sections", "num_chordwise_panels"})

# A Wing never stores explode_into_strips, recording only its effect on spanwise_mesh,
# so a copy could not carry it even if a build tried. Omitting it is safe only because
# both analysis functions in convergence.py reject a Wing whose spanwise mesh is neither
# trapezoidal nor edge-defined, which leaves False as the only value a Wing reaching a
# build can have been built with.
_WING_OMITTED_PARAMETERS = frozenset({"explode_into_strips"})

# A build shares each reference WingCrossSection's Airfoil rather than rebuilding one,
# because refining a mesh changes nothing about an Airfoil.
_WING_CROSS_SECTION_COPIED_PARAMETERS = frozenset(
    {
        "chord",
        "Lp_Wcsp_Lpp",
        "angles_Wcsp_to_Wcs_ixyz",
        "control_surface_symmetry_type",
        "control_surface_hinge_point",
        "control_surface_deflection",
        "spanwise_spacing",
        "airfoil",
    }
)

# A build resolves each WingCrossSection copy's number of spanwise Panels for this
# iteration's mesh.
_WING_CROSS_SECTION_CHANGED_PARAMETERS = frozenset({"num_spanwise_panels"})
_WING_CROSS_SECTION_OMITTED_PARAMETERS: frozenset[str] = frozenset()


def _public_constructor_parameters(this_class) -> set[str]:
    """Returns the names of a class's public constructor parameters.

    A private parameter is internal machinery rather than part of the class's public
    construction contract, so it is excluded along with self.

    :param this_class: The class whose initialization method is read.
    :return: The set of the class's public constructor parameter names.
    """
    return {
        name
        for name in inspect.signature(this_class.__init__).parameters
        if name != "self" and not name.startswith("_")
    }


class TestGetWingSectionNumSpanwisePanels(unittest.TestCase):
    """This class contains methods for testing
    _convergence_meshing._get_wing_section_num_spanwise_panels, the search that picks a
    non-edge-defined Wing section's number of spanwise Panels.
    """

    def setUp(self) -> None:
        """Set up a root and tip WingCrossSection defining one wing section.

        :return: None
        """
        self.root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        self.tip_wing_cross_section = (
            geometry_fixtures.make_tip_wing_cross_section_fixture()
        )
        self.num_chordwise_panels = 4
        self.chordwise_spacing = "uniform"

    def _average_panel_aspect_ratio(self, num_spanwise_panels) -> float:
        """Meshes the wing section at a number of spanwise Panels and returns its average
        Panel aspect ratio.
        """
        return _convergence_meshing._get_wing_section_average_panel_aspect_ratio(
            self.num_chordwise_panels,
            self.chordwise_spacing,
            self.root_wing_cross_section,
            self.tip_wing_cross_section,
            num_spanwise_panels,
        )

    def _num_spanwise_panels_for(self, target) -> int:
        """Searches for the number of spanwise Panels that hits a target average Panel
        aspect ratio, starting from the smallest valid count.
        """
        return _convergence_meshing._get_wing_section_num_spanwise_panels(
            desired_average_panel_aspect_ratio=target,
            num_chordwise_panels=self.num_chordwise_panels,
            chordwise_spacing=self.chordwise_spacing,
            ref_root_wing_cross_section=self.root_wing_cross_section,
            ref_tip_wing_cross_section=self.tip_wing_cross_section,
            start_val=1,
        )

    def test_returns_int(self) -> None:
        """Test that the number of spanwise Panels is returned as an int."""
        self.assertIsInstance(self._num_spanwise_panels_for(4), int)

    def test_result_is_closest_to_target(self) -> None:
        """Test that the chosen number of spanwise Panels gives an average Panel aspect
        ratio at least as close to the target as either neighboring count.
        """
        target = 4
        result = self._num_spanwise_panels_for(target)

        chosen_difference = abs(self._average_panel_aspect_ratio(result) - target)
        upper_difference = abs(self._average_panel_aspect_ratio(result + 1) - target)
        self.assertLessEqual(chosen_difference, upper_difference)

        # The smallest valid number of spanwise Panels is one, so the lower neighbor
        # only exists above that.
        if result > 1:
            lower_difference = abs(
                self._average_panel_aspect_ratio(result - 1) - target
            )
            self.assertLessEqual(chosen_difference, lower_difference)

    def test_smaller_target_needs_at_least_as_many_panels(self) -> None:
        """Test that a smaller target Panel aspect ratio, a finer mesh, needs at least as
        many spanwise Panels as a larger one.
        """
        self.assertGreaterEqual(
            self._num_spanwise_panels_for(2), self._num_spanwise_panels_for(4)
        )


class TestGetNumWingCrossSectionsForPanelAr(unittest.TestCase):
    """This class contains methods for testing
    _convergence_meshing._get_num_wing_cross_sections_for_panel_ar, the search that picks an
    edge-defined Wing's number of WingCrossSections.
    """

    @staticmethod
    def _tapered_edge_points() -> tuple[np.ndarray, np.ndarray]:
        """Builds straight, tapered leading and trailing edge curves.

        The leading edge sweeps back from the origin to (0.5, 2.0, 0.0) and the trailing
        edge is unswept at x = 1.0, so the chord tapers from a unit root chord to a 0.5
        tip chord over a two meter half span.
        """
        ys = np.linspace(0.0, 2.0, 20)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.25 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        return leading, trailing

    def _make_edge_wing(self, symmetric=False) -> ps.geometry.wing.Wing:
        """Builds an edge-defined Wing from the tapered edge curves."""
        leading, trailing = self._tapered_edge_points()
        return ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            name="Edge Wing",
            symmetric=symmetric,
            symmetryNormal_G=(0.0, 1.0, 0.0) if symmetric else None,
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0) if symmetric else None,
            num_chordwise_panels=4,
        )

    def setUp(self) -> None:
        """Set up a reference edge-defined Wing.

        :return: None
        """
        self.ref_wing = self._make_edge_wing()
        self.num_chordwise_panels = 4

    def _average_panel_aspect_ratio(self, num_wing_cross_sections) -> float:
        """Rebuilds and meshes the edge-defined Wing at a number of WingCrossSections and
        returns its average Panel aspect ratio.
        """
        refined_wing = _convergence_meshing._build_edge_defined_wing(
            self.ref_wing, self.num_chordwise_panels, num_wing_cross_sections
        )
        airplane = ps.geometry.airplane.Airplane(wings=[refined_wing])
        average_panel_aspect_ratio = airplane.wings[0].average_panel_aspect_ratio
        assert average_panel_aspect_ratio is not None
        return average_panel_aspect_ratio

    def _num_wing_cross_sections_for(self, target, ref_wing=None) -> int:
        """Searches for the number of WingCrossSections that hits a target average Panel
        aspect ratio, starting from the smallest valid count.
        """
        return _convergence_meshing._get_num_wing_cross_sections_for_panel_ar(
            desired_average_panel_aspect_ratio=target,
            num_chordwise_panels=self.num_chordwise_panels,
            ref_wing=self.ref_wing if ref_wing is None else ref_wing,
            start_val=2,
        )

    def test_returns_int(self) -> None:
        """Test that the number of WingCrossSections is returned as an int."""
        self.assertIsInstance(self._num_wing_cross_sections_for(4), int)

    def test_result_is_closest_to_target(self) -> None:
        """Test that the chosen number of WingCrossSections gives an average Panel aspect
        ratio at least as close to the target as either neighboring count.
        """
        target = 4
        result = self._num_wing_cross_sections_for(target)

        chosen_difference = abs(self._average_panel_aspect_ratio(result) - target)
        upper_difference = abs(self._average_panel_aspect_ratio(result + 1) - target)
        self.assertLessEqual(chosen_difference, upper_difference)

        # An edge-defined Wing needs at least two WingCrossSections, so the lower
        # neighbor only exists above that.
        if result > 2:
            lower_difference = abs(
                self._average_panel_aspect_ratio(result - 1) - target
            )
            self.assertLessEqual(chosen_difference, lower_difference)

    def test_smaller_target_needs_at_least_as_many_wing_cross_sections(self) -> None:
        """Test that a smaller target Panel aspect ratio, a finer mesh, needs at least as
        many WingCrossSections as a larger one.
        """
        self.assertGreaterEqual(
            self._num_wing_cross_sections_for(2), self._num_wing_cross_sections_for(4)
        )

    def test_symmetric_matches_asymmetric(self) -> None:
        """Test that the count is measured on the half span, so a symmetric Wing and an
        asymmetric Wing built from the same half-span curves need the same count.
        """
        asymmetric_result = self._num_wing_cross_sections_for(
            4, ref_wing=self._make_edge_wing(symmetric=False)
        )
        symmetric_result = self._num_wing_cross_sections_for(
            4, ref_wing=self._make_edge_wing(symmetric=True)
        )
        self.assertEqual(asymmetric_result, symmetric_result)


class TestMemosComplete(unittest.TestCase):
    """This class contains methods for testing _convergence_meshing.memos_complete, the
    check that reports whether one mesh's build would only re-record already cached
    memos.
    """

    def setUp(self) -> None:
        """Set up reference Airplanes (one with a trapezoidal Wing and one with an
        edge-defined Wing) and memo caches holding every memo their build would record
        at the first mesh.

        :return: None
        """
        trapezoidal_wing = geometry_fixtures.make_three_section_wing_fixture()

        ys = np.linspace(0.0, 2.0, 20)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.25 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        edge_defined_wing = ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            name="Edge Wing",
            num_chordwise_panels=4,
        )

        self.ref_airplanes = [
            ps.geometry.airplane.Airplane(wings=[trapezoidal_wing]),
            ps.geometry.airplane.Airplane(wings=[edge_defined_wing]),
        ]

        # The trapezoidal Wing has three WingCrossSections, so its build records one
        # number of spanwise Panels memo per non-tip WingCrossSection. The edge-defined
        # Wing's build records one number of WingCrossSections memo.
        self.num_spanwise_panels_cache = {
            (0, 0, 0, 0, 0): 8,
            (0, 0, 0, 0, 1): 8,
        }
        self.num_wing_cross_sections_cache = {(0, 0, 1, 0): 5}

    def _memos_complete(self, ar_id=0, chord_id=0, delta_time_cache=None) -> bool:
        """Calls memos_complete against the reference Airplanes and memo caches."""
        return _convergence_meshing.memos_complete(
            ar_id,
            chord_id,
            self.ref_airplanes,
            self.num_spanwise_panels_cache,
            self.num_wing_cross_sections_cache,
            delta_time_cache,
        )

    def test_true_when_all_memos_cached(self) -> None:
        """Test that the check passes when every memo the build would record is
        cached.
        """
        self.assertTrue(self._memos_complete())

    def test_false_when_spanwise_panels_memo_missing(self) -> None:
        """Test that the check fails when a trapezoidal Wing's number of spanwise
        Panels memo is missing.
        """
        del self.num_spanwise_panels_cache[(0, 0, 0, 0, 1)]
        self.assertFalse(self._memos_complete())

    def test_false_when_wing_cross_sections_memo_missing(self) -> None:
        """Test that the check fails when an edge-defined Wing's number of
        WingCrossSections memo is missing.
        """
        del self.num_wing_cross_sections_cache[(0, 0, 1, 0)]
        self.assertFalse(self._memos_complete())

    def test_false_when_memos_are_for_another_mesh(self) -> None:
        """Test that memos cached for one mesh do not count toward a different mesh."""
        self.assertFalse(self._memos_complete(chord_id=1))

    def test_delta_time_checked_only_when_cache_given(self) -> None:
        """Test that the delta_time memo is required when a delta_time cache is given
        and ignored when it is None.
        """
        self.assertTrue(self._memos_complete(delta_time_cache=None))
        self.assertFalse(self._memos_complete(delta_time_cache={}))
        self.assertTrue(self._memos_complete(delta_time_cache={(0, 0): 0.01}))


class TestGeometryCopyParameterCoverage(unittest.TestCase):
    """This class contains methods for testing that this module's partitions of each
    geometry class's public constructor parameters stay complete.

    _convergence_meshing.build_steady_problem and
    _convergence_meshing.build_unsteady_problem each rebuild every Airplane, Wing,
    WingCrossSection, and Airfoil by naming its parameters one by one, so a parameter that
    a construction site forgets to pass silently falls back to its default and the
    analysis refines geometry that its reference problem never described. That is not
    hypothetical: it is what happened to WingMovement's rotationPointOffset_Gs_Ler, and
    these classes gain parameters at a similar rate. Wing alone has gained mirror_only,
    symmetryNormal_G, explode_into_strips, and tip_trim_fraction.

    These tests therefore fail whenever one of those classes gains a public constructor
    parameter that this module has not classified as copied, changed, or omitted. Making
    a failure here go away means passing the new parameter wherever both builders
    construct the class, not only classifying it here.
    """

    @staticmethod
    def _coverage_failure_message(class_name) -> str:
        """Builds the message shown when a geometry class's public constructor parameters
        are no longer the ones this module classifies.
        """
        return (
            f"This module's partition of {class_name}'s public constructor parameters "
            f"is no longer complete. If a parameter was added to {class_name}, "
            f"classifying it here is only half of the fix: "
            f"_convergence_meshing.build_steady_problem and "
            f"_convergence_meshing.build_unsteady_problem each rebuild {class_name} by "
            f"naming its parameters one by one, so a copied parameter must also be "
            f"passed at every site where either function constructs {class_name}. "
            f"Otherwise it falls back to its default and every convergence iteration "
            f"silently refines geometry that the reference problem never described. "
            f"Classify it as copied only if a build should carry it from the reference, "
            f"and add it to the copy tests below so that a build that drops it fails."
        )

    def test_airplane_parameters_are_classified(self) -> None:
        """Test that every one of Airplane's public constructor parameters is classified
        as copied, changed, or omitted.
        """
        self.assertEqual(
            _AIRPLANE_COPIED_PARAMETERS
            | _AIRPLANE_CHANGED_PARAMETERS
            | _AIRPLANE_OMITTED_PARAMETERS,
            _public_constructor_parameters(ps.geometry.airplane.Airplane),
            msg=self._coverage_failure_message("Airplane"),
        )

    def test_wing_parameters_are_classified(self) -> None:
        """Test that every one of Wing's public constructor parameters is classified as
        copied, changed, or omitted.
        """
        self.assertEqual(
            _WING_COPIED_PARAMETERS
            | _WING_CHANGED_PARAMETERS
            | _WING_OMITTED_PARAMETERS,
            _public_constructor_parameters(ps.geometry.wing.Wing),
            msg=self._coverage_failure_message("Wing"),
        )

    def test_wing_cross_section_parameters_are_classified(self) -> None:
        """Test that every one of WingCrossSection's public constructor parameters is
        classified as copied, changed, or omitted.
        """
        self.assertEqual(
            _WING_CROSS_SECTION_COPIED_PARAMETERS
            | _WING_CROSS_SECTION_CHANGED_PARAMETERS
            | _WING_CROSS_SECTION_OMITTED_PARAMETERS,
            _public_constructor_parameters(
                ps.geometry.wing_cross_section.WingCrossSection
            ),
            msg=self._coverage_failure_message("WingCrossSection"),
        )

    def test_partitions_do_not_overlap(self) -> None:
        """Test that no parameter is classified twice, since a parameter that a build both
        copies and changes would describe two different behaviors.
        """
        for copied, changed, omitted in (
            (
                _AIRPLANE_COPIED_PARAMETERS,
                _AIRPLANE_CHANGED_PARAMETERS,
                _AIRPLANE_OMITTED_PARAMETERS,
            ),
            (
                _WING_COPIED_PARAMETERS,
                _WING_CHANGED_PARAMETERS,
                _WING_OMITTED_PARAMETERS,
            ),
            (
                _WING_CROSS_SECTION_COPIED_PARAMETERS,
                _WING_CROSS_SECTION_CHANGED_PARAMETERS,
                _WING_CROSS_SECTION_OMITTED_PARAMETERS,
            ),
        ):
            self.assertEqual(len(copied & changed), 0)
            self.assertEqual(len(copied & omitted), 0)
            self.assertEqual(len(changed & omitted), 0)


class TestBuildSteadyProblem(unittest.TestCase):
    """This class contains methods for testing
    _convergence_meshing.build_steady_problem, the builder that copies a reference
    SteadyProblem's Airplanes at one mesh.

    The reference holds two Airplanes so that one build covers both refinement branches:
    the first Airplane's Wing is trapezoidal, which is refined by resolving each non tip
    WingCrossSection's number of spanwise Panels, and the second Airplane's Wing is
    edge-defined, which is refined by resampling its stored edge curves into a new number
    of WingCrossSections.
    """

    # The reference Airplanes' indices within the reference SteadyProblem, which are
    # also the indices their copies take and part of each cache key.
    trapezoidal_airplane_id = 0
    edge_defined_airplane_id = 1

    def setUp(self) -> None:
        """Set up a reference SteadyProblem holding a trapezoidal Wing and an
        edge-defined Wing, and the empty caches that a build fills.

        :return: None
        """
        self.trapezoidal_wing = geometry_fixtures.make_three_section_wing_fixture()

        ys = np.linspace(0.0, 2.0, 20)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.25 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        self.edge_defined_wing = ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            name="Edge Wing",
            num_chordwise_panels=4,
        )

        # The first Airplane in a simulation must keep a zero Cg_GP1_CgP1, so only the
        # second carries the nonzero value that the copy is checked for.
        self.ref_problem = ps.problems.SteadyProblem(
            airplanes=[
                ps.geometry.airplane.Airplane(
                    wings=[self.trapezoidal_wing],
                    name="Trapezoidal Airplane",
                    weight=10.0,
                ),
                ps.geometry.airplane.Airplane(
                    wings=[self.edge_defined_wing],
                    name="Edge Airplane",
                    Cg_GP1_CgP1=(1.0, 2.0, 3.0),
                    weight=20.0,
                ),
            ],
            operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
        )

        self.num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int] = {}
        self.num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int] = {}

    def _build(self, num_chordwise_panels=4) -> ps.problems.SteadyProblem:
        """Builds the SteadyProblem for the first mesh against the reference problem and
        the caches.
        """
        return _convergence_meshing.build_steady_problem(
            ar_id=0,
            chord_id=0,
            panel_aspect_ratio=4,
            num_chordwise_panels=num_chordwise_panels,
            ref_problem=self.ref_problem,
            num_spanwise_panels_cache=self.num_spanwise_panels_cache,
            num_wing_cross_sections_cache=self.num_wing_cross_sections_cache,
        )

    def test_returns_steady_problem(self) -> None:
        """Test that the build returns a SteadyProblem holding one copy of each reference
        Airplane.
        """
        this_problem = self._build()

        self.assertIsInstance(this_problem, ps.problems.SteadyProblem)
        self.assertEqual(len(this_problem.airplanes), 2)

    def test_operating_point_is_deep_copied(self) -> None:
        """Test that the built problem holds its own equal OperatingPoint rather than the
        reference's, so that solving it cannot populate the reference OperatingPoint's
        lazy caches and change the reference problem's content hash.
        """
        this_operating_point = self._build().operating_point
        ref_operating_point = self.ref_problem.operating_point

        self.assertIsNot(this_operating_point, ref_operating_point)
        self.assertEqual(this_operating_point.rho, ref_operating_point.rho)
        self.assertEqual(this_operating_point.vCg__E, ref_operating_point.vCg__E)
        self.assertEqual(this_operating_point.alpha, ref_operating_point.alpha)
        self.assertEqual(this_operating_point.beta, ref_operating_point.beta)

    def test_reference_problem_is_not_modified(self) -> None:
        """Test that the build leaves the reference problem's Wings at their own mesh, so
        that a later iteration still refines from the reference rather than from a
        previous iteration's result.
        """
        self._build(num_chordwise_panels=4)

        self.assertEqual(self.trapezoidal_wing.num_chordwise_panels, 8)
        self.assertEqual(len(self.edge_defined_wing.wing_cross_sections), 5)

    def test_num_chordwise_panels_is_applied_to_every_wing(self) -> None:
        """Test that every copied Wing takes this iteration's number of chordwise Panels
        rather than its reference's.
        """
        this_problem = self._build(num_chordwise_panels=6)

        for this_airplane in this_problem.airplanes:
            for this_wing in this_airplane.wings:
                self.assertEqual(this_wing.num_chordwise_panels, 6)

    def test_spanwise_mesh_types_are_preserved(self) -> None:
        """Test that a trapezoidal Wing's copy stays trapezoidal and an edge-defined
        Wing's copy stays edge-defined, so that a later iteration refines each by the same
        branch.
        """
        this_problem = self._build()

        self.assertEqual(
            this_problem.airplanes[self.trapezoidal_airplane_id].wings[0].spanwise_mesh,
            "trapezoidal",
        )
        self.assertEqual(
            this_problem.airplanes[self.edge_defined_airplane_id]
            .wings[0]
            .spanwise_mesh,
            "edge_defined",
        )

    def test_trapezoidal_wing_records_a_memo_per_non_tip_wing_cross_section(
        self,
    ) -> None:
        """Test that refining the trapezoidal Wing records one number of spanwise Panels
        memo for each of its non tip WingCrossSections, keyed by the mesh, Airplane, Wing,
        and WingCrossSection indices.
        """
        self._build()

        self.assertEqual(
            set(self.num_spanwise_panels_cache),
            {
                (0, 0, self.trapezoidal_airplane_id, 0, 0),
                (0, 0, self.trapezoidal_airplane_id, 0, 1),
            },
        )

    def test_edge_defined_wing_records_one_memo(self) -> None:
        """Test that refining the edge-defined Wing records one number of
        WingCrossSections memo, keyed by the mesh, Airplane, and Wing indices.
        """
        self._build()

        self.assertEqual(
            set(self.num_wing_cross_sections_cache),
            {(0, 0, self.edge_defined_airplane_id, 0)},
        )

    def test_cached_num_spanwise_panels_is_reused(self) -> None:
        """Test that a cached number of spanwise Panels is used as is rather than
        resolved again, so that one mesh's Wing sections resolve only once across a
        sweep.
        """
        self.num_spanwise_panels_cache[(0, 0, self.trapezoidal_airplane_id, 0, 0)] = 7

        this_wing = self._build().airplanes[self.trapezoidal_airplane_id].wings[0]

        self.assertEqual(this_wing.wing_cross_sections[0].num_spanwise_panels, 7)

    def test_cached_num_wing_cross_sections_is_reused(self) -> None:
        """Test that a cached number of WingCrossSections is used as is rather than
        resolved again, so that one mesh's edge-defined Wing resamples only once across a
        sweep.
        """
        self.num_wing_cross_sections_cache[(0, 0, self.edge_defined_airplane_id, 0)] = 3

        this_wing = self._build().airplanes[self.edge_defined_airplane_id].wings[0]

        self.assertEqual(len(this_wing.wing_cross_sections), 3)

    def _assert_copied(self, this_object, ref_object, parameter_names) -> None:
        """Asserts that a copied geometry object carries each named parameter from its
        reference.

        A parameter's value can be an array, a string, a bool, a number, or None, so array
        values are compared for closeness and the rest for equality.
        """
        for name in sorted(parameter_names):
            this_value = getattr(this_object, name)
            ref_value = getattr(ref_object, name)
            message = (
                f"The build's copy of {type(ref_object).__name__} did not carry {name}, "
                f"so a convergence analysis would refine geometry that its reference "
                f"problem never described. Pass {name} at every site where "
                f"build_steady_problem and build_unsteady_problem construct "
                f"{type(ref_object).__name__}."
            )
            if isinstance(ref_value, np.ndarray):
                np.testing.assert_allclose(this_value, ref_value, err_msg=message)
            else:
                self.assertEqual(this_value, ref_value, msg=message)

    def test_airplane_parameters_are_copied(self) -> None:
        """Test that an Airplane copy carries every one of its reference's copied
        parameters.
        """
        self._assert_copied(
            self._build().airplanes[self.edge_defined_airplane_id],
            self.ref_problem.airplanes[self.edge_defined_airplane_id],
            _AIRPLANE_COPIED_PARAMETERS,
        )

    def test_airplane_reference_values_are_recalculated(self) -> None:
        """Test that each Airplane copy calculates its own reference dimensions, because
        refining the mesh changes the geometry they describe.
        """
        this_airplane = self._build().airplanes[self.trapezoidal_airplane_id]

        self.assertIsNotNone(this_airplane.s_ref)
        self.assertIsNotNone(this_airplane.c_ref)
        self.assertIsNotNone(this_airplane.b_ref)

    def test_wing_parameters_are_copied(self) -> None:
        """Test that a trapezoidal Wing's copy carries every one of its reference's copied
        parameters.
        """
        self._assert_copied(
            self._build().airplanes[self.trapezoidal_airplane_id].wings[0],
            self.trapezoidal_wing,
            _WING_COPIED_PARAMETERS,
        )

    def test_wing_cross_section_parameters_are_copied(self) -> None:
        """Test that a trapezoidal Wing's copied WingCrossSections each carry every one of
        their reference's copied parameters.
        """
        these_wing_cross_sections = (
            self._build()
            .airplanes[self.trapezoidal_airplane_id]
            .wings[0]
            .wing_cross_sections
        )
        ref_wing_cross_sections = self.trapezoidal_wing.wing_cross_sections

        self.assertEqual(len(these_wing_cross_sections), len(ref_wing_cross_sections))
        for this_wing_cross_section, ref_wing_cross_section in zip(
            these_wing_cross_sections, ref_wing_cross_sections
        ):
            self._assert_copied(
                this_wing_cross_section,
                ref_wing_cross_section,
                _WING_CROSS_SECTION_COPIED_PARAMETERS,
            )

    def test_airfoil_is_shared_rather_than_rebuilt(self) -> None:
        """Test that each copied WingCrossSection holds its reference's own Airfoil, so
        that its outline is reproduced exactly.

        Rebuilding an Airfoil from a reference's outline revalidates and renormalizes it,
        which perturbs a cambered outline by a small amount rather than reproducing it, so
        each copy must share the reference's Airfoil instead. Sharing is safe because
        Airfoils are immutable.
        """
        these_wing_cross_sections = (
            self._build()
            .airplanes[self.trapezoidal_airplane_id]
            .wings[0]
            .wing_cross_sections
        )
        ref_wing_cross_sections = self.trapezoidal_wing.wing_cross_sections

        for this_wing_cross_section, ref_wing_cross_section in zip(
            these_wing_cross_sections, ref_wing_cross_sections
        ):
            self.assertIs(
                this_wing_cross_section.airfoil, ref_wing_cross_section.airfoil
            )
            np.testing.assert_array_equal(
                this_wing_cross_section.airfoil.outline_A_lp,
                ref_wing_cross_section.airfoil.outline_A_lp,
            )


class TestBuildUnsteadyProblemCopiesMotion(unittest.TestCase):
    """This class contains methods for testing that
    _convergence_meshing.build_unsteady_problem's AirplaneMovement, WingMovement, and
    WingCrossSectionMovement copies each carry every one of their reference's motion
    parameters.

    build_unsteady_problem rebuilds each of these movement classes by naming their
    parameters one by one, so a parameter that a construction site forgets to pass
    silently falls back to its default and every convergence iteration then solves motion
    that the reference problem never described. That is not hypothetical: it is what
    happened to WingMovement's rotationPointOffset_Gs_Ler, which changed a flapping
    case's mean thrust by 47 percent before it was caught.

    These tests therefore exist to fail whenever a motion parameter is added to one of
    those classes without also being passed at every one of build_unsteady_problem's
    construction sites for it. Making a failure here go away means adding the parameter
    to that function, not only to this class.
    """

    # Each dict names every motion parameter of one movement class, and gives it a value
    # that differs from its default so that a copy that dropped it would fall back to
    # that default and fail the comparison. A value equal to the parameter's default
    # would make these tests pass while the copy silently drops it, so the difference is
    # what gives them their teeth.
    #
    # A test below checks each dict against its class's signature. If one of those tests
    # fails because a parameter was added to a movement class, adding it here is only
    # half of the fix: it must also be passed wherever
    # _convergence_meshing.build_unsteady_problem constructs that class, or the
    # convergence analysis will quietly ignore it.
    #
    # Each dict's values are annotated as Any because they hold both float tuples and
    # spacing string tuples, so a narrower type could not be unpacked into any one
    # parameter.
    airplane_motion_parameters: dict[str, Any] = {
        "ampCg_GP1_CgP1": (0.1, 0.2, 0.3),
        "periodCg_GP1_CgP1": (1.0, 2.0, 3.0),
        "spacingCg_GP1_CgP1": ("uniform", "sine", "sine"),
        "phaseCg_GP1_CgP1": (10.0, 20.0, 30.0),
    }
    wing_motion_parameters: dict[str, Any] = {
        "ampLer_Gs_Cgs": (0.1, 0.2, 0.3),
        "periodLer_Gs_Cgs": (1.0, 2.0, 3.0),
        "spacingLer_Gs_Cgs": ("uniform", "sine", "sine"),
        "phaseLer_Gs_Cgs": (10.0, 20.0, 30.0),
        "ampAngles_Gs_to_Wn_ixyz": (1.0, 2.0, 3.0),
        "periodAngles_Gs_to_Wn_ixyz": (1.0, 2.0, 3.0),
        "spacingAngles_Gs_to_Wn_ixyz": ("sine", "uniform", "sine"),
        "phaseAngles_Gs_to_Wn_ixyz": (40.0, 50.0, 60.0),
        "rotationPointOffset_Gs_Ler": (0.01, 0.02, 0.03),
    }
    wing_cross_section_motion_parameters: dict[str, Any] = {
        "ampLp_Wcsp_Lpp": (0.01, 0.02, 0.03),
        "periodLp_Wcsp_Lpp": (1.0, 2.0, 3.0),
        "spacingLp_Wcsp_Lpp": ("uniform", "sine", "sine"),
        "phaseLp_Wcsp_Lpp": (10.0, 20.0, 30.0),
        "ampAngles_Wcsp_to_Wcs_ixyz": (1.0, 2.0, 3.0),
        "periodAngles_Wcsp_to_Wcs_ixyz": (1.0, 2.0, 3.0),
        "spacingAngles_Wcsp_to_Wcs_ixyz": ("sine", "uniform", "sine"),
        "phaseAngles_Wcsp_to_Wcs_ixyz": (40.0, 50.0, 60.0),
    }

    def setUp(self) -> None:
        """Set up a reference UnsteadyProblem whose Airplanes hold a trapezoidal Wing and
        an edge-defined Wing, each wrapped in an AirplaneMovement and a WingMovement
        carrying every motion parameter.

        :return: None
        """
        trapezoidal_wing = geometry_fixtures.make_three_section_wing_fixture()

        ys = np.linspace(0.0, 2.0, 20)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.25 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        edge_defined_wing = ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            name="Edge Wing",
            num_chordwise_panels=4,
        )

        # Each motion is carried by the one Airplane that can hold it, and the copies of
        # each are checked separately below. An edge-defined Wing's
        # WingCrossSectionMovements must all be static, because resampling changes the
        # WingCrossSection count and so cannot preserve per WingCrossSection motion, so
        # only the trapezoidal Wing's carry motion. The first Airplane in a simulation
        # must keep a zero Cg_GP1_CgP1, so only the second Airplane's AirplaneMovement
        # carries motion. Both Airplanes' WingMovements carry motion.
        airplane_movements = []
        for wing, wing_cross_section_motion, airplane_motion in (
            (trapezoidal_wing, self.wing_cross_section_motion_parameters, {}),
            (edge_defined_wing, {}, self.airplane_motion_parameters),
        ):
            airplane_movements.append(
                ps.movements.airplane_movement.AirplaneMovement(
                    base_airplane=ps.geometry.airplane.Airplane(wings=[wing]),
                    wing_movements=[
                        ps.movements.wing_movement.WingMovement(
                            base_wing=wing,
                            wing_cross_section_movements=[
                                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                                    base_wing_cross_section=wing_cross_section,
                                    # A root WingCrossSection must keep a zero
                                    # Lp_Wcsp_Lpp and zero angles_Wcsp_to_Wcs_ixyz, so
                                    # only the WingCrossSections outboard of it can
                                    # carry motion.
                                    **(
                                        {}
                                        if wing_cross_section_id == 0
                                        else wing_cross_section_motion
                                    ),
                                )
                                for wing_cross_section_id, wing_cross_section in (
                                    enumerate(wing.wing_cross_sections)
                                )
                            ],
                            **self.wing_motion_parameters,
                        )
                    ],
                    **airplane_motion,
                )
            )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        self.ref_problem = ps.problems.UnsteadyProblem(
            movement=ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=1,
                # An explicit, coarse time step matters here even though these tests
                # never solve anything: Movement generates its meshed geometry at every
                # time step during construction, and without this it first runs its
                # delta_time optimizer, which resolves a step over a thousand times
                # finer. This setUp runs once per test in this class.
                delta_time=0.5,
            )
        )

    def _built_airplane_movements(self) -> list:
        """Builds the UnsteadyProblem for the first mesh and returns its
        AirplaneMovements, one per reference Airplane, with the trapezoidal Wing's first
        and the edge-defined Wing's second.

        The delta_time cache is seeded so that the build reuses a time step instead of
        running Movement's iterative optimizer, and the seeded value is coarse so that
        the built Movement generates only a handful of time steps of geometry.
        """
        this_problem = _convergence_meshing.build_unsteady_problem(
            ar_id=0,
            chord_id=0,
            panel_aspect_ratio=4,
            num_chordwise_panels=4,
            wake_length=1,
            ref_problem=self.ref_problem,
            num_spanwise_panels_cache={},
            num_wing_cross_sections_cache={},
            delta_time_cache={(0, 0): 0.5},
        )
        return list(this_problem.movement.airplane_movements)

    def _built_wing_movements(self) -> list:
        """Builds the UnsteadyProblem for the first mesh and returns its WingMovements,
        one per reference Airplane.
        """
        return [
            this_airplane_movement.wing_movements[0]
            for this_airplane_movement in self._built_airplane_movements()
        ]

    @staticmethod
    def _coverage_failure_message(class_name) -> str:
        """Builds the message shown when a movement class's motion parameters are no
        longer the set that this class tests.
        """
        return (
            f"This class's motion parameters no longer match {class_name}'s. If a "
            f"parameter was added to {class_name}, listing it here is only half of the "
            f"fix: _convergence_meshing.build_unsteady_problem rebuilds each "
            f"{class_name} by naming its parameters one by one, so the new parameter "
            f"must also be passed at every site where that function constructs "
            f"{class_name}. Otherwise it falls back to its default and every "
            f"convergence iteration silently solves motion that the reference problem "
            f"never described. Give it a value here that differs from its default, so "
            f"that the copy tests below would catch that."
        )

    def _assert_parameters_match(self, actual_movement, expected_parameters) -> None:
        """Asserts that a copied movement class carries each expected motion parameter.

        The spacing parameters hold strings rather than numbers, so they are compared for
        equality while the rest are compared for closeness.
        """
        for name, expected in expected_parameters.items():
            message = (
                f"build_unsteady_problem's copy of "
                f"{type(actual_movement).__name__} did not carry {name}, so a "
                f"convergence analysis would solve motion that its reference problem "
                f"never described. Pass {name} at every site where "
                f"build_unsteady_problem constructs "
                f"{type(actual_movement).__name__}."
            )
            actual = getattr(actual_movement, name)
            if "spacing" in name:
                self.assertEqual(tuple(actual), expected, msg=message)
            else:
                np.testing.assert_allclose(actual, expected, err_msg=message)

    def test_airplane_motion_parameters_covers_every_constructor_parameter(
        self,
    ) -> None:
        """Test that airplane_motion_parameters names every AirplaneMovement constructor
        parameter except the two that the build is meant to change, so that a parameter
        added to AirplaneMovement later fails here until it is covered by the tests
        below.
        """
        parameters = inspect.signature(
            ps.movements.airplane_movement.AirplaneMovement.__init__
        ).parameters

        self.assertEqual(
            set(self.airplane_motion_parameters),
            set(parameters) - {"self", "base_airplane", "wing_movements"},
            msg=self._coverage_failure_message("AirplaneMovement"),
        )

    def test_wing_motion_parameters_covers_every_constructor_parameter(self) -> None:
        """Test that wing_motion_parameters names every WingMovement constructor
        parameter except the two that the build is meant to change, so that a parameter
        added to WingMovement later fails here until it is covered by the tests below.
        """
        parameters = inspect.signature(
            ps.movements.wing_movement.WingMovement.__init__
        ).parameters

        self.assertEqual(
            set(self.wing_motion_parameters),
            set(parameters) - {"self", "base_wing", "wing_cross_section_movements"},
            msg=self._coverage_failure_message("WingMovement"),
        )

    def test_wing_cross_section_motion_parameters_covers_every_constructor_parameter(
        self,
    ) -> None:
        """Test that wing_cross_section_motion_parameters names every
        WingCrossSectionMovement constructor parameter except the one that the build is
        meant to change, so that a parameter added to WingCrossSectionMovement later
        fails here until it is covered by the tests below.
        """
        parameters = inspect.signature(
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement.__init__
        ).parameters

        self.assertEqual(
            set(self.wing_cross_section_motion_parameters),
            set(parameters) - {"self", "base_wing_cross_section"},
            msg=self._coverage_failure_message("WingCrossSectionMovement"),
        )

    def test_rotation_point_offset_is_copied(self) -> None:
        """Test that both the trapezoidal and the edge-defined branch copy the reference
        WingMovement's rotationPointOffset_Gs_Ler rather than defaulting it to the Wing
        root leading edge.
        """
        for this_wing_movement in self._built_wing_movements():
            np.testing.assert_allclose(
                this_wing_movement.rotationPointOffset_Gs_Ler,
                self.wing_motion_parameters["rotationPointOffset_Gs_Ler"],
            )

    def test_every_airplane_motion_parameter_is_copied(self) -> None:
        """Test that an AirplaneMovement copy carries every one of its reference's motion
        parameters, so that a parameter added later cannot be silently dropped.

        Only the second Airplane is checked, because the first Airplane in a simulation
        must keep a zero Cg_GP1_CgP1 and so cannot carry this motion. The build copies
        every AirplaneMovement with the same code, so the second covers it.
        """
        self._assert_parameters_match(
            self._built_airplane_movements()[1], self.airplane_motion_parameters
        )

    def test_every_wing_cross_section_motion_parameter_is_copied(self) -> None:
        """Test that the trapezoidal branch's WingCrossSectionMovement copies carry every
        one of their reference's motion parameters, so that a parameter added later
        cannot be silently dropped.

        Only the trapezoidal branch is checked, because the edge-defined branch resamples
        its WingCrossSections and so deliberately builds motion free
        WingCrossSectionMovements. The root WingCrossSectionMovement is skipped because a
        root WingCrossSection cannot move.
        """
        trapezoidal_wing_movement = self._built_wing_movements()[0]
        these_wing_cross_section_movements = (
            trapezoidal_wing_movement.wing_cross_section_movements[1:]
        )

        self.assertEqual(len(these_wing_cross_section_movements), 2)
        for this_wing_cross_section_movement in these_wing_cross_section_movements:
            self._assert_parameters_match(
                this_wing_cross_section_movement,
                self.wing_cross_section_motion_parameters,
            )

    def test_every_wing_motion_parameter_is_copied(self) -> None:
        """Test that both branches copy every one of the WingMovement's motion
        parameters, so that a parameter added later cannot be silently dropped.
        """
        for this_wing_movement in self._built_wing_movements():
            self._assert_parameters_match(
                this_wing_movement, self.wing_motion_parameters
            )
