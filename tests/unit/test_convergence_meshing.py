"""This module contains classes to test the convergence meshing functions."""

import unittest

import numpy as np

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _convergence_meshing
from tests.unit.fixtures import geometry_fixtures


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

        # The smallest valid number of spanwise Panels is one, so the lower neighbor only
        # exists above that.
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

        # An edge-defined Wing needs at least two WingCrossSections, so the lower neighbor
        # only exists above that.
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
