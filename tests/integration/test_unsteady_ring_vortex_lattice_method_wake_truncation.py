"""This is a testing case for wake truncation in the unsteady ring vortex lattice method
solver.

It runs a static geometry unsteady simulation with and without wake truncation, and
verifies that the truncated simulation caps wake size at max_wake_rows and produces
loads that are close to the non truncated simulation.
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import airplane_fixtures, operating_point_fixtures


def _make_static_movement_components():
    """Creates the shared AirplaneMovement and OperatingPointMovement components for the
    wake truncation tests.

    :return: A tuple of (airplane_movement, operating_point_movement).
    """
    airplane = airplane_fixtures.make_symmetric_unsteady_validation_airplane()
    operating_point = operating_point_fixtures.make_validation_operating_point()

    root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=airplane.wings[0].wing_cross_sections[0]
        )
    )

    tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=airplane.wings[0].wing_cross_sections[1]
        )
    )

    wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=[
            root_wing_cross_section_movement,
            tip_wing_cross_section_movement,
        ],
    )

    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[wing_movement],
    )

    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point
        )
    )

    return airplane_movement, operating_point_movement


class TestWakeTruncation(unittest.TestCase):
    """This is a class for testing wake truncation in the
    UnsteadyRingVortexLatticeMethodSolver."""

    @classmethod
    def setUpClass(cls):
        """Set up both truncated and non truncated solvers."""
        airplane_movement, operating_point_movement = _make_static_movement_components()

        # Non truncated: 6 chord lengths of wake.
        non_truncated_movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            num_chords=6,
        )
        non_truncated_problem = ps.problems.UnsteadyProblem(
            movement=non_truncated_movement
        )
        cls.non_truncated_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            non_truncated_problem
        )
        cls.non_truncated_solver.run(
            prescribed_wake=True,
            calculate_streamlines=False,
            show_progress=False,
        )

        # Truncated: uses the same num_steps and delta_time but caps the wake at 3 chord
        # lengths.
        truncated_movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=non_truncated_movement.delta_time,
            num_steps=non_truncated_movement.num_steps,
            max_wake_chords=3,
        )
        truncated_problem = ps.problems.UnsteadyProblem(movement=truncated_movement)
        cls.truncated_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            truncated_problem
        )
        cls.truncated_solver.run(
            prescribed_wake=True,
            calculate_streamlines=False,
            show_progress=False,
        )

        cls.max_wake_rows = truncated_movement.max_wake_rows

    def test_wake_size_is_capped(self):
        """Test that the truncated simulation's wake has exactly max_wake_rows chordwise
        rows at the final time step."""
        solver = self.truncated_solver
        final_step = solver.num_steps - 1

        # Total spanwise panel count summed over all wings of all airplanes.
        total_spanwise_panels = sum(
            wing.num_spanwise_panels
            for airplane in solver.current_airplanes
            for wing in airplane.wings
        )
        expected_num_wake_vortices = self.max_wake_rows * total_spanwise_panels

        self.assertEqual(
            solver.list_num_wake_vortices[final_step], expected_num_wake_vortices
        )

    def test_non_truncated_wake_is_larger(self):
        """Test that the non truncated simulation's wake has more chordwise rows than
        the truncated simulation's wake at the final time step."""
        non_truncated_final_step = self.non_truncated_solver.num_steps - 1
        truncated_final_step = self.truncated_solver.num_steps - 1

        non_truncated_num_wake_vortices = (
            self.non_truncated_solver.list_num_wake_vortices[non_truncated_final_step]
        )
        truncated_num_wake_vortices = self.truncated_solver.list_num_wake_vortices[
            truncated_final_step
        ]

        self.assertGreater(non_truncated_num_wake_vortices, truncated_num_wake_vortices)

    def test_truncated_loads_close_to_non_truncated(self):
        """Test that the truncated simulation's loads are close to the non truncated
        simulation's loads."""
        non_truncated_airplane = self.non_truncated_solver.current_airplanes[0]
        truncated_airplane = self.truncated_solver.current_airplanes[0]

        # Compare lift coefficients (CL = -forceCoefficients_W[2]).
        non_truncated_cl = -non_truncated_airplane.forceCoefficients_W[2]
        truncated_cl = -truncated_airplane.forceCoefficients_W[2]

        # The truncated loads should be within 10% of the non truncated loads.
        cl_error = abs(truncated_cl - non_truncated_cl) / abs(non_truncated_cl)
        self.assertLess(cl_error, 0.10)

    def test_pre_allocated_wake_arrays_capped(self):
        """Test that the truncated solver's pre-allocated wake arrays are correctly
        sized with truncation."""
        num_steps = self.truncated_solver.num_steps

        # After the wake reaches max_wake_rows, the number of wake vortices should
        # plateau.
        for step in range(self.max_wake_rows + 1, num_steps):
            self.assertEqual(
                self.truncated_solver.list_num_wake_vortices[step],
                self.truncated_solver.list_num_wake_vortices[self.max_wake_rows],
            )

    def test_wake_point_grid_is_capped(self):
        """Test that the truncated simulation's wake point grid has the correct number
        of rows."""
        final_airplane = self.truncated_solver.current_airplanes[0]

        for wing in final_airplane.wings:
            grid = wing.gridWrvp_GP1_CgP1
            self.assertIsNotNone(grid)

            # The point grid has one more row than the vortex grid.
            self.assertEqual(grid.shape[0], self.max_wake_rows + 1)


if __name__ == "__main__":
    unittest.main()
