"""This module contains classes to test FreeFlightUnsteadyProblems."""

import unittest
from unittest.mock import MagicMock

import numpy as np

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _transformations
from tests.unit.fixtures import problem_fixtures


def _movement_and_mass():
    """Return a fresh FreeFlightMovement and the mass consistent with its Airplane's
    weight and gravitational field (weight == mass * |g_E|).

    The mass is taken from a basic fixture problem, so it stays consistent with whatever
    weight and gravity the fixture uses.

    :return: A tuple of the FreeFlightMovement and the consistent mass in kilograms.
    """
    fixture_problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
    # noinspection PyProtectedMember
    return fixture_problem._free_flight_movement, fixture_problem.mass


class TestFreeFlightUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test FreeFlightUnsteadyProblems."""

    def setUp(self):
        """Set up a fresh FreeFlightUnsteadyProblem for each test."""
        self.problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test FreeFlightUnsteadyProblem initialization with valid parameters."""
        self.assertIsInstance(
            self.problem,
            ps.problems.FreeFlightUnsteadyProblem,
        )
        self.assertIsInstance(
            self.problem,
            ps.problems._CoupledUnsteadyProblem,
        )
        self.assertIsInstance(
            self.problem.movement,
            ps.movements.free_flight_movement.FreeFlightMovement,
        )

    def test_movement_type_validation(self):
        """Test that movement must be a FreeFlightMovement."""
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement="not_a_movement",
                mass=1.0,
                I_BP1_CgP1=np.eye(3, dtype=float),
            )

    def test_single_airplane_movement_validation(self):
        """Test that the FreeFlightMovement has exactly one FreeFlightAirplaneMovement."""
        movement, mass = _movement_and_mass()
        two_airplane_movement = ps.movements.free_flight_movement.FreeFlightMovement(
            airplane_movements=[
                movement.airplane_movements[0],
                movement.airplane_movements[0],
            ],
            operating_point_movement=movement.operating_point_movement,
            delta_time=movement.delta_time,
            prescribed_num_steps=3,
            free_num_steps=2,
        )
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=two_airplane_movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
            )

    def test_I_BP1_CgP1_symmetry_validation(self):
        """Test that I_BP1_CgP1 must be symmetric."""
        movement, mass = _movement_and_mass()
        asymmetric_inertia = np.array(
            [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
        )
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=asymmetric_inertia,
            )

    def test_external_loads_fn_validation(self):
        """Test that external_loads_fn must be callable or None."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                external_loads_fn="not_callable",
            )

    def test_extra_xml_type_validation(self):
        """Test that extra_xml must be a dict or None."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                extra_xml="invalid",
            )

    def test_extra_xml_key_validation(self):
        """Test that an extra_xml key must be a permitted injection point."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                extra_xml={"wordbody": "<geom/>"},
            )

    def test_extra_xml_value_validation(self):
        """Test that an extra_xml value must be a str."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                extra_xml={"visual": 123},
            )

    def test_mujoco_assets_type_validation(self):
        """Test that mujoco_assets must be a dict or None."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                mujoco_assets="invalid",
            )

    def test_mujoco_assets_key_validation(self):
        """Test that a mujoco_assets key must be a str."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                mujoco_assets={123: b"data"},
            )

    def test_mujoco_assets_value_validation(self):
        """Test that a mujoco_assets value must be bytes."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                mujoco_assets={"dummy.stl": "not bytes"},
            )

    def test_integrator_type_validation(self):
        """Test that integrator must be a str."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                integrator=4,
            )

    def test_integrator_value_validation(self):
        """Test that integrator must be a supported MuJoCo integrator."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                integrator="rk4",
            )

    def test_externalFX_W_validation(self):
        """Test that a nonzero externalFX_W on the initial OperatingPoint raises."""
        base_operating_point = ps.operating_point.OperatingPoint(externalFX_W=10.0)
        with self.assertRaises(ValueError):
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture(
                base_operating_point=base_operating_point
            )

    def test_integrator_forwarded_to_mujoco_model(self):
        """Test that the integrator choice reaches the generated MuJoCo XML."""
        movement, mass = _movement_and_mass()
        problem = ps.problems.FreeFlightUnsteadyProblem(
            movement=movement,
            mass=mass,
            I_BP1_CgP1=np.eye(3, dtype=float),
            integrator="implicitfast",
        )
        self.assertIn('integrator="implicitfast"', problem._mujoco_model.xml_str)

    def test_k_max_type_validation(self):
        """Test that k_max must be an int."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                k_max=20.0,
            )

    def test_k_max_value_validation(self):
        """Test that k_max must be greater than zero."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(3, dtype=float),
                k_max=0,
            )

    def test_k_max_default(self):
        """Test that k_max defaults to 20."""
        self.assertEqual(self.problem.k_max, 20)

    def test_k_max_stored(self):
        """Test that a valid k_max is stored and returned."""
        movement, mass = _movement_and_mass()
        problem = ps.problems.FreeFlightUnsteadyProblem(
            movement=movement,
            mass=mass,
            I_BP1_CgP1=np.eye(3, dtype=float),
            k_max=7,
        )
        self.assertEqual(problem.k_max, 7)

    def test_external_loads_fn_default_none(self):
        """Test that external_loads_fn defaults to None."""
        self.assertIsNone(self.problem.external_loads_fn)

    def test_external_loads_fn_stored_when_callable(self):
        """Test that a callable external_loads_fn is stored and returned."""

        # noinspection PyUnusedLocal
        def external_loads_fn(operating_point, airplane):
            return np.zeros(3, dtype=float), np.zeros(3, dtype=float)

        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture(
            external_loads_fn=external_loads_fn
        )
        self.assertIs(problem.external_loads_fn, external_loads_fn)

    def test_only_final_results_forced_false(self):
        """Test that only_final_results is always False for coupled problems."""
        self.assertFalse(self.problem.only_final_results)

    def test_num_steps_from_movement(self):
        """Test that num_steps is taken from the movement."""
        self.assertEqual(self.problem.num_steps, self.problem.movement.num_steps)

    def test_delta_time_from_movement(self):
        """Test that delta_time is taken from the movement."""
        self.assertEqual(self.problem.delta_time, self.problem.movement.delta_time)

    def test_steady_problems_initialization(self):
        """Test that steady_problems is initialized with one SteadyProblem."""
        self.assertIsInstance(self.problem.steady_problems, tuple)
        self.assertEqual(len(self.problem.steady_problems), 1)
        self.assertIsInstance(
            self.problem.steady_problems[0],
            ps.problems.SteadyProblem,
        )

    def test_initialization_of_load_lists(self):
        """Test that load lists are initialized as empty."""
        self.assertIsInstance(self.problem.forces_W, list)
        self.assertEqual(len(self.problem.forces_W), 0)

        self.assertIsInstance(self.problem.forceCoefficients_W, list)
        self.assertEqual(len(self.problem.forceCoefficients_W), 0)

        self.assertIsInstance(self.problem.moments_W_Cg, list)
        self.assertEqual(len(self.problem.moments_W_Cg), 0)

        self.assertIsInstance(self.problem.momentCoefficients_W_Cg, list)
        self.assertEqual(len(self.problem.momentCoefficients_W_Cg), 0)

    def test_I_BP1_CgP1_attribute(self):
        """Test that I_BP1_CgP1 is stored correctly."""
        self.assertIsInstance(self.problem.I_BP1_CgP1, np.ndarray)
        self.assertEqual(self.problem.I_BP1_CgP1.shape, (3, 3))
        np.testing.assert_array_equal(self.problem.I_BP1_CgP1, np.diag([1.0, 1.0, 1.0]))

    def test_I_BP1_CgP1_accepts_array_like(self):
        """Test that I_BP1_CgP1 accepts a nested list and converts it to a float
        ndarray.
        """
        movement, mass = _movement_and_mass()
        inertia_list = [
            [2, 0, 0],
            [0, 3, 0],
            [0, 0, 4],
        ]
        problem = ps.problems.FreeFlightUnsteadyProblem(
            movement=movement,
            mass=mass,
            I_BP1_CgP1=inertia_list,
        )
        self.assertIsInstance(problem.I_BP1_CgP1, np.ndarray)
        self.assertEqual(problem.I_BP1_CgP1.dtype, np.dtype(float))
        np.testing.assert_array_equal(problem.I_BP1_CgP1, np.diag([2.0, 3.0, 4.0]))

    def test_I_BP1_CgP1_rejects_wrong_shape(self):
        """Test that I_BP1_CgP1 must have shape (3, 3)."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass,
                I_BP1_CgP1=np.eye(2, dtype=float),
            )

    def test_free_flight_movement_property_narrows_movement(self):
        """Test that the _free_flight_movement property returns the same object as
        movement, narrowed to FreeFlightMovement.
        """
        self.assertIs(self.problem._free_flight_movement, self.problem.movement)
        self.assertIsInstance(
            self.problem._free_flight_movement,
            ps.movements.free_flight_movement.FreeFlightMovement,
        )

    def test_mujoco_model_attribute(self):
        """Test that the MuJoCoModel is constructed during initialization."""
        # noinspection PyProtectedMember
        from pterasoftware import _mujoco_model

        self.assertIsInstance(self.problem._mujoco_model, _mujoco_model.MuJoCoModel)

    def test_mass_attribute(self):
        """Test that mass is stored as a float."""
        movement, mass = _movement_and_mass()
        problem = ps.problems.FreeFlightUnsteadyProblem(
            movement=movement,
            mass=mass,
            I_BP1_CgP1=np.eye(3, dtype=float),
        )
        self.assertIsInstance(problem.mass, float)
        self.assertAlmostEqual(problem.mass, mass)

    def test_mass_rejects_non_positive(self):
        """Test that mass must be greater than zero."""
        movement, mass = _movement_and_mass()
        for invalid_mass in (0.0, -1.0):
            with self.assertRaises(ValueError):
                ps.problems.FreeFlightUnsteadyProblem(
                    movement=movement,
                    mass=invalid_mass,
                    I_BP1_CgP1=np.eye(3, dtype=float),
                )

    def test_weight_mass_gravity_consistency_validation(self):
        """Test that the Airplane's weight must equal mass * |g_E| within tolerance."""
        movement, mass = _movement_and_mass()
        with self.assertRaises(ValueError):
            ps.problems.FreeFlightUnsteadyProblem(
                movement=movement,
                mass=mass * 2.0,
                I_BP1_CgP1=np.eye(3, dtype=float),
            )


class TestFreeFlightUnsteadyProblemInitializeNextProblem(unittest.TestCase):
    """Tests for FreeFlightUnsteadyProblem.initialize_next_problem."""

    def setUp(self):
        """Set up a fresh FreeFlightUnsteadyProblem and a primed mock solver."""
        self.problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )
        self.final_step = self.problem.num_steps - 1

        # The solver populates the current Airplane's load attributes via
        # _calculate_loads before invoking initialize_next_problem, so emulate that here.
        self.current_airplane = self.problem.steady_problems[0].airplanes[0]
        self.current_airplane.forces_W = np.array([1.0, 2.0, 3.0], dtype=float)
        self.current_airplane.forceCoefficients_W = np.array(
            [0.1, 0.2, 0.3], dtype=float
        )
        self.current_airplane.moments_W_CgP1 = np.array([4.0, 5.0, 6.0], dtype=float)
        self.current_airplane.momentCoefficients_W_CgP1 = np.array(
            [0.4, 0.5, 0.6], dtype=float
        )
        self.current_operating_point = self.problem.steady_problems[0].operating_point

        self.solver = MagicMock()
        self.solver.current_airplanes = [self.current_airplane]
        self.solver.current_operating_point = self.current_operating_point

    def _mock_mujoco_model(self):
        """Replace the problem's MuJoCoModel with a MagicMock returning a fixed state.

        This keeps the rigid body dynamics deterministic and avoids stepping the real
        physics engine during the unit test.

        :return: The MagicMock that replaced the problem's MuJoCoModel.
        """
        mock_model = MagicMock()
        mock_model.get_state.return_value = {
            "position_E_Eo": np.array([1.0, 0.0, 0.0], dtype=float),
            "R_pas_E_to_BP1": np.eye(3, dtype=float),
            "velocity_E__E": np.array([10.0, 0.0, 0.0], dtype=float),
            "omegas_BP1__E": np.zeros(3, dtype=float),
            "time": self.problem.delta_time,
        }
        self.problem._mujoco_model = mock_model
        return mock_model

    @staticmethod
    def _primed_problem_and_solver(external_loads_fn):
        """Build a problem carrying the given external_loads_fn and a primed mock solver.

        Mirrors setUp's load priming and _mock_mujoco_model, but for a fresh problem that
        carries an external_loads_fn so its first-invocation return validation can be
        exercised on a non final step.

        :param external_loads_fn: The external_loads_fn to attach to the problem.
        :return: A tuple of the FreeFlightUnsteadyProblem and the primed mock solver.
        """
        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture(
            external_loads_fn=external_loads_fn
        )
        airplane = problem.steady_problems[0].airplanes[0]
        airplane.forces_W = np.array([1.0, 2.0, 3.0], dtype=float)
        airplane.forceCoefficients_W = np.array([0.1, 0.2, 0.3], dtype=float)
        airplane.moments_W_CgP1 = np.array([4.0, 5.0, 6.0], dtype=float)
        airplane.momentCoefficients_W_CgP1 = np.array([0.4, 0.5, 0.6], dtype=float)

        solver = MagicMock()
        solver.current_airplanes = [airplane]
        solver.current_operating_point = problem.steady_problems[0].operating_point

        mock_model = MagicMock()
        mock_model.get_state.return_value = {
            "position_E_Eo": np.array([1.0, 0.0, 0.0], dtype=float),
            "R_pas_E_to_BP1": np.eye(3, dtype=float),
            "velocity_E__E": np.array([10.0, 0.0, 0.0], dtype=float),
            "omegas_BP1__E": np.zeros(3, dtype=float),
            "time": problem.delta_time,
        }
        problem._mujoco_model = mock_model
        return problem, solver

    def test_appends_next_steady_problem_on_non_final_step(self):
        """Test that a new SteadyProblem is appended on a non final step."""
        self._mock_mujoco_model()
        before = len(self.problem.steady_problems)

        self.problem.initialize_next_problem(self.solver, step=0)

        self.assertEqual(len(self.problem.steady_problems), before + 1)
        self.assertIsInstance(
            self.problem.steady_problems[-1], ps.problems.SteadyProblem
        )

    def test_appends_next_operating_point_on_non_final_step(self):
        """Test that a new OperatingPoint is appended to the movement on a non final
        step.
        """
        self._mock_mujoco_model()
        operating_points = (
            self.problem._free_flight_movement.operating_point_movement.operating_points
        )
        before = len(operating_points)

        self.problem.initialize_next_problem(self.solver, step=0)

        self.assertEqual(len(operating_points), before + 1)

    def test_steps_dynamics_but_withholds_loads_on_prescribed_step(self):
        """Test that on a prescribed phase step the MuJoCo model is stepped but no loads
        are applied, so the body coasts at its trimmed condition.
        """
        mock_model = self._mock_mujoco_model()

        # Step 0 is always in the prescribed phase, since prescribed_num_steps is at
        # least one.
        self.problem.initialize_next_problem(self.solver, step=0)

        mock_model.apply_loads.assert_not_called()
        mock_model.step.assert_called_once()

    def test_records_loads_on_non_final_step(self):
        """Test that the current Airplane's loads are recorded on a non final step."""
        self._mock_mujoco_model()

        self.problem.initialize_next_problem(self.solver, step=0)

        self.assertEqual(len(self.problem.forces_W), 1)
        np.testing.assert_array_equal(
            self.problem.forces_W[0], self.current_airplane.forces_W
        )
        np.testing.assert_array_equal(
            self.problem.forceCoefficients_W[0],
            self.current_airplane.forceCoefficients_W,
        )
        np.testing.assert_array_equal(
            self.problem.moments_W_Cg[0], self.current_airplane.moments_W_CgP1
        )
        np.testing.assert_array_equal(
            self.problem.momentCoefficients_W_Cg[0],
            self.current_airplane.momentCoefficients_W_CgP1,
        )

    def test_no_steady_problem_appended_on_final_step(self):
        """Test that no SteadyProblem is appended on the final step."""
        before = len(self.problem.steady_problems)

        self.problem.initialize_next_problem(self.solver, step=self.final_step)

        self.assertEqual(len(self.problem.steady_problems), before)

    def test_records_loads_on_final_step(self):
        """Test that the current Airplane's loads are still recorded on the final
        step.
        """
        self.problem.initialize_next_problem(self.solver, step=self.final_step)

        self.assertEqual(len(self.problem.forces_W), 1)
        np.testing.assert_array_equal(
            self.problem.forces_W[0], self.current_airplane.forces_W
        )

    def test_does_not_step_dynamics_on_final_step(self):
        """Test that the MuJoCo model is not stepped on the final step."""
        mock_model = self._mock_mujoco_model()

        self.problem.initialize_next_problem(self.solver, step=self.final_step)

        mock_model.apply_loads.assert_not_called()
        mock_model.step.assert_not_called()

    def test_recorded_loads_are_copies(self):
        """Test that the recorded loads are copies, so later mutation of the Airplane's
        loads does not change the recorded history.
        """
        # Use the final step so no rigid body dynamics integration is required.
        self.problem.initialize_next_problem(self.solver, step=self.final_step)

        self.current_airplane.forces_W[0] = 999.0

        self.assertNotEqual(self.problem.forces_W[0][0], 999.0)

    def test_external_loads_fn_invoked_on_non_final_step(self):
        """Test that a set external_loads_fn is invoked with the current OperatingPoint
        and Airplane on a non final step.
        """
        external_loads_fn = MagicMock(
            return_value=(np.zeros(3, dtype=float), np.zeros(3, dtype=float))
        )
        problem, solver = self._primed_problem_and_solver(external_loads_fn)
        airplane = solver.current_airplanes[0]
        operating_point = solver.current_operating_point

        problem.initialize_next_problem(solver, step=0)

        external_loads_fn.assert_called_once_with(operating_point, airplane)

    def test_external_loads_fn_not_invoked_on_later_prescribed_step(self):
        """Test that the external_loads_fn is not invoked on a prescribed phase step
        after the first, since its loads would be withheld there in any case.
        """
        external_loads_fn = MagicMock(
            return_value=(np.zeros(3, dtype=float), np.zeros(3, dtype=float))
        )
        problem, solver = self._primed_problem_and_solver(external_loads_fn)

        # The basic fixture has more than one prescribed step, so step 1 is in the
        # prescribed phase but is not the first step.
        self.assertGreater(problem._free_flight_movement.prescribed_num_steps, 1)
        problem.initialize_next_problem(solver, step=1)

        external_loads_fn.assert_not_called()

    def test_external_loads_fn_valid_return_passes(self):
        """Test that a well formed external_loads_fn return passes validation and marks
        the return as validated.
        """
        external_loads_fn = MagicMock(
            return_value=(
                np.array([1.0, 0.0, 0.0], dtype=float),
                np.array([0.0, 1.0, 0.0], dtype=float),
            )
        )
        problem, solver = self._primed_problem_and_solver(external_loads_fn)

        problem.initialize_next_problem(solver, step=0)

        self.assertTrue(problem._external_loads_validated)

    def test_external_loads_fn_wrong_arity_raises(self):
        """Test that an external_loads_fn returning other than two items raises."""
        external_loads_fn = MagicMock(return_value=(np.zeros(3, dtype=float),))
        problem, solver = self._primed_problem_and_solver(external_loads_fn)

        with self.assertRaises(ValueError):
            problem.initialize_next_problem(solver, step=0)

    def test_external_loads_fn_wrong_shape_raises(self):
        """Test that an external_loads_fn returning a non (3,) component raises."""
        external_loads_fn = MagicMock(
            return_value=(np.zeros(2, dtype=float), np.zeros(3, dtype=float))
        )
        problem, solver = self._primed_problem_and_solver(external_loads_fn)

        with self.assertRaises(ValueError):
            problem.initialize_next_problem(solver, step=0)

    def test_external_loads_fn_non_finite_raises(self):
        """Test that an external_loads_fn returning a non finite value raises."""
        external_loads_fn = MagicMock(
            return_value=(
                np.array([np.inf, 0.0, 0.0], dtype=float),
                np.zeros(3, dtype=float),
            )
        )
        problem, solver = self._primed_problem_and_solver(external_loads_fn)

        with self.assertRaises(ValueError):
            problem.initialize_next_problem(solver, step=0)


class TestFreeFlightUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for FreeFlightUnsteadyProblem attribute immutability."""

    def setUp(self):
        """Set up a fresh FreeFlightUnsteadyProblem for each immutability test."""
        self.problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )

    def test_immutable_movement_property(self):
        """Test that the movement property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.movement = None

    def test_immutable_I_BP1_CgP1_property(self):
        """Test that the I_BP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.I_BP1_CgP1 = np.eye(3, dtype=float)

    def test_I_BP1_CgP1_array_not_writeable(self):
        """Test that the I_BP1_CgP1 numpy array is not writeable."""
        with self.assertRaises(ValueError):
            self.problem.I_BP1_CgP1[0, 0] = 999.0

    def test_immutable_mass_property(self):
        """Test that the mass property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.mass = 1.0

    def test_immutable_external_loads_fn_property(self):
        """Test that the external_loads_fn property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.external_loads_fn = None

    def test_immutable_num_steps_property(self):
        """Test that the num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.num_steps = 100

    def test_immutable_delta_time_property(self):
        """Test that the delta_time property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.delta_time = 0.5

    def test_mutable_load_lists(self):
        """Test that load lists remain mutable for solver population."""
        self.problem.forces_W.append(np.array([1.0, 2.0, 3.0]))
        self.assertEqual(len(self.problem.forces_W), 1)


class TestFreeFlightUnsteadyProblemStateConversions(unittest.TestCase):
    """Tests for the rigid body state conversions the sub-iteration relies on."""

    def setUp(self):
        """Set up a fresh FreeFlightUnsteadyProblem and its reference OperatingPoint."""
        self.problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )
        self.reference_operating_point = self.problem.steady_problems[0].operating_point

    def test_state_to_vector_block_layout_and_units(self):
        """Test that _state_to_vector concatenates the position, attitude, velocity, and
        angular rate blocks in order, with the angular quantities converted to radians.

        With an identity orientation the attitude block is zero and the body angular
        rate carries into the Earth frame unchanged, so the only conversion to check is
        degrees per second to radians per second on the angular rate.
        """
        state = {
            "position_E_Eo": np.array([1.0, 2.0, 3.0], dtype=float),
            "R_pas_E_to_BP1": np.eye(3, dtype=float),
            "velocity_E__E": np.array([4.0, 5.0, 6.0], dtype=float),
            "omegas_BP1__E": np.array([7.0, 8.0, 9.0], dtype=float),
            "time": self.problem.delta_time,
        }

        x = self.problem._state_to_vector(state)

        expected_x = np.concatenate(
            [
                np.array([1.0, 2.0, 3.0]),
                np.zeros(3),
                np.array([4.0, 5.0, 6.0]),
                np.deg2rad(np.array([7.0, 8.0, 9.0])),
            ]
        )
        self.assertEqual(x.shape, (12,))
        np.testing.assert_allclose(x, expected_x)

    def test_state_vector_round_trips_through_operating_point(self):
        """Test that a state converted to the 12-vector and back into an OperatingPoint
        reproduces the position, speed, attitude, and body angular rate.

        The round trip pins the attitude and angular rate conversions in
        _state_to_vector against their inverses in _operating_point_from_vector, which a
        single one-way check cannot, since the angular rate is rotated into the Earth
        frame one way and back the other.
        """
        angles_E_to_BP1_izyx = np.array([10.0, 20.0, 30.0], dtype=float)
        R_pas_E_to_BP1 = _transformations.generate_rot_T(
            angles=angles_E_to_BP1_izyx,
            passive=True,
            intrinsic=True,
            order="zyx",
        )[:3, :3]
        position_E_Eo = np.array([1.0, -2.0, 3.0], dtype=float)
        velocity_E__E = np.array([12.0, -3.0, 4.0], dtype=float)
        omegas_BP1__E = np.array([5.0, -6.0, 7.0], dtype=float)
        state = {
            "position_E_Eo": position_E_Eo,
            "R_pas_E_to_BP1": R_pas_E_to_BP1,
            "velocity_E__E": velocity_E__E,
            "omegas_BP1__E": omegas_BP1__E,
            "time": self.problem.delta_time,
        }

        x = self.problem._state_to_vector(state)
        operating_point = self.problem._operating_point_from_vector(
            x, self.reference_operating_point
        )

        np.testing.assert_allclose(operating_point.CgP1_E_Eo, position_E_Eo)
        self.assertAlmostEqual(
            operating_point.vCg__E, float(np.linalg.norm(velocity_E__E))
        )
        np.testing.assert_allclose(
            operating_point.angles_E_to_BP1_izyx, angles_E_to_BP1_izyx
        )
        np.testing.assert_allclose(operating_point.omegas_BP1__E, omegas_BP1__E)

    def test_operating_point_from_vector_carries_environment_from_reference(self):
        """Test that _operating_point_from_vector carries the environmental quantities
        across from the reference OperatingPoint unchanged.

        Only the rigid body state comes from the 12-vector; the fluid density, kinematic
        viscosity, gravity, surface geometry, and external force are inherited.
        """
        x = np.zeros(12, dtype=float)
        x[6] = 10.0

        operating_point = self.problem._operating_point_from_vector(
            x, self.reference_operating_point
        )

        reference = self.reference_operating_point
        self.assertEqual(operating_point.rho, reference.rho)
        self.assertEqual(operating_point.nu, reference.nu)
        np.testing.assert_array_equal(operating_point.g_E, reference.g_E)
        np.testing.assert_array_equal(
            operating_point.surfaceNormal_E, reference.surfaceNormal_E
        )
        np.testing.assert_array_equal(
            operating_point.surfacePoint_E_Eo, reference.surfacePoint_E_Eo
        )
        self.assertEqual(operating_point.externalFX_W, reference.externalFX_W)


class TestFreeFlightUnsteadyProblemRelaxationWeights(unittest.TestCase):
    """Tests for the sub-iteration weighting matrix diagonal."""

    def test_build_relaxation_weights_matches_formula(self):
        """Test that _build_relaxation_weights builds the weighting diagonal from the
        reference chord, mass, inertia, and time step per its definition.

        A non-identity inertia is used so the mean of the principal moments (one third of
        the trace) is exercised rather than collapsing to one. The velocity and angular
        rate blocks must equal the position and attitude blocks scaled by the time step.
        """
        movement, mass = _movement_and_mass()
        problem = ps.problems.FreeFlightUnsteadyProblem(
            movement=movement,
            mass=mass,
            I_BP1_CgP1=np.diag([2.0, 3.0, 4.0]),
        )

        weights = problem._build_relaxation_weights()

        L_ref = problem._free_flight_movement.airplanes[0][0].c_ref
        mean_inertia = (2.0 + 3.0 + 4.0) / 3.0
        theta_ref = mass * L_ref**2 / mean_inertia
        delta_time = problem.delta_time
        expected_weights = np.concatenate(
            [
                np.full(3, 1.0 / L_ref),
                np.full(3, 1.0 / theta_ref),
                np.full(3, delta_time / L_ref),
                np.full(3, delta_time / theta_ref),
            ]
        )

        self.assertEqual(weights.shape, (12,))
        np.testing.assert_allclose(weights, expected_weights)

        # The rate blocks are the configuration blocks scaled by the time step.
        np.testing.assert_allclose(weights[6:9], delta_time * weights[0:3])
        np.testing.assert_allclose(weights[9:12], delta_time * weights[3:6])
