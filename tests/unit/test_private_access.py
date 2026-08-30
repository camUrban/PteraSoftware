"""This module contains classes to test the private access registration functions."""

import unittest
from typing import Any

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _mujoco_model, _private_access
from tests.unit.fixtures import problem_fixtures


class TestRegisterMuJoCoModelGetter(unittest.TestCase):
    """This class contains methods for testing
    _private_access.register_mujoco_model_getter."""

    def test_importing_problems_registers_a_getter(self) -> None:
        """Test that importing problems.py registers a getter at import time."""
        # noinspection PyProtectedMember
        self.assertIsNotNone(_private_access._mujoco_model_getter)

    def test_register_replaces_the_getter(self) -> None:
        """Test that registering a new getter replaces the previously registered one."""
        # Save the registered getter and restore it after the test, since the
        # replacement below would otherwise leak into the rest of the suite.
        # noinspection PyProtectedMember
        original_getter = _private_access._mujoco_model_getter
        assert original_getter is not None
        self.addCleanup(_private_access.register_mujoco_model_getter, original_getter)

        sentinel_model: Any = object()
        recorded_problems: list[ps.problems.FreeFlightUnsteadyProblem] = []

        def stub_getter(
            problem: ps.problems.FreeFlightUnsteadyProblem,
        ) -> Any:
            """Records the FreeFlightUnsteadyProblem and returns the sentinel."""
            recorded_problems.append(problem)
            return sentinel_model

        _private_access.register_mujoco_model_getter(stub_getter)

        sentinel_problem: Any = object()
        self.assertIs(
            _private_access.get_mujoco_model(sentinel_problem), sentinel_model
        )
        self.assertEqual(recorded_problems, [sentinel_problem])


class TestGetMuJoCoModel(unittest.TestCase):
    """This class contains methods for testing _private_access.get_mujoco_model."""

    basic_free_flight_unsteady_problem: ps.problems.FreeFlightUnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up the shared test fixtures."""
        cls.basic_free_flight_unsteady_problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )

    def test_returns_the_problems_mujoco_model(self) -> None:
        """Test that the accessor returns the FreeFlightUnsteadyProblem's own
        MuJoCoModel."""
        model = _private_access.get_mujoco_model(
            self.basic_free_flight_unsteady_problem
        )
        self.assertIsInstance(model, _mujoco_model.MuJoCoModel)
        # noinspection PyProtectedMember
        self.assertIs(model, self.basic_free_flight_unsteady_problem._mujoco_model)
