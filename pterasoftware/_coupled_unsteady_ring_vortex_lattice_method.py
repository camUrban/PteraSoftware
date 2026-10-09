"""Contains the CoupledUnsteadyRingVortexLatticeMethodSolver class."""

from __future__ import annotations

from typing import cast

from . import _logging, problems
from .unsteady_ring_vortex_lattice_method import UnsteadyRingVortexLatticeMethodSolver

_logger = _logging.get_logger("_coupled_unsteady_ring_vortex_lattice_method")


class CoupledUnsteadyRingVortexLatticeMethodSolver(
    UnsteadyRingVortexLatticeMethodSolver
):
    """A subclass of UnsteadyRingVortexLatticeMethodSolver that solves
    CoupledUnsteadyProblems.

    Geometry in a CoupledUnsteadyProblem is determined step by step from the solver's
    results at the previous step, so bound vortices cannot be initialized upfront. This
    class inherits the parent's run() and initialize_step_geometry() unchanged and
    overrides three hooks: _initialize_step_vortices (per step bound vortex init),
    _update_next_step_hook (calls CoupledUnsteadyProblem.initialize_next_problem between
    steps), and _get_steady_problem_at (dynamic dispatch through the problem's
    get_steady_problem accessor).
    """

    __slots__ = ()

    def __init__(self, unsteady_problem: problems.CoupledUnsteadyProblem) -> None:
        """The initialization method.

        :param unsteady_problem: The CoupledUnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(unsteady_problem, problems.CoupledUnsteadyProblem):
            raise TypeError("unsteady_problem must be a CoupledUnsteadyProblem.")
        super().__init__(unsteady_problem)

    @property
    def _coupled_unsteady_problem(self) -> problems.CoupledUnsteadyProblem:
        """Type narrowed view of the inherited unsteady_problem attribute.

        The parent stores unsteady_problem as a CoreUnsteadyProblem (widened to let
        subclasses pass their own variants). __init__ validates that this subclass
        always receives a CoupledUnsteadyProblem, so the cast here is safe.

        :return: The unsteady_problem narrowed to CoupledUnsteadyProblem.
        """
        return cast(problems.CoupledUnsteadyProblem, self.unsteady_problem)

    def _initialize_step_vortices(self, step: int) -> None:
        _logger.debug(
            _logging.indent() + f"Initializing step {step}'s bound ring vortices"
        )
        self._initialize_panel_vortices_at(step)

    def _update_next_step_hook(self, step: int) -> None:
        self._coupled_unsteady_problem.initialize_next_problem(self, step)
        if step < self.num_steps - 1:
            self._initialize_panel_vortices_at(step + 1)

    def _get_steady_problem_at(self, step: int) -> problems.SteadyProblem:
        return self._coupled_unsteady_problem.get_steady_problem(step)
