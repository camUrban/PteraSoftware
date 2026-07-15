"""Contains the CoupledUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

CoupledUnsteadyRingVortexLatticeMethodSolver: A subclass of
UnsteadyRingVortexLatticeMethodSolver that solves _CoupledUnsteadyProblems, whose
geometry is initialized and updated step by step rather than being fully precomputed.

**Contains the following functions:**

None
"""

from __future__ import annotations

from typing import cast

from . import _logging, problems
from .unsteady_ring_vortex_lattice_method import UnsteadyRingVortexLatticeMethodSolver

_logger = _logging.get_logger("_coupled_unsteady_ring_vortex_lattice_method")


class CoupledUnsteadyRingVortexLatticeMethodSolver(
    UnsteadyRingVortexLatticeMethodSolver
):
    """A subclass of UnsteadyRingVortexLatticeMethodSolver that solves
    _CoupledUnsteadyProblems.

    Geometry in a _CoupledUnsteadyProblem is determined step by step from the solver's
    results at the previous step, so bound vortices cannot be initialized upfront. This
    class inherits the parent's run() and initialize_step_geometry() unchanged and
    overrides three hooks: _initialize_step_vortices (per step bound vortex init),
    _update_next_step_hook (calls _CoupledUnsteadyProblem.initialize_next_problem
    between steps), and _get_steady_problem_at (dynamic dispatch through the problem's
    get_steady_problem accessor).

    **Contains the following methods:**

    None
    """

    __slots__ = ()

    # _CoupledUnsteadyProblem is deliberately given a private class name inside the
    # public problems module so it stays out of the user-facing API. Referencing it from
    # this in-package solver is the intended design, so silence PyCharm's cross-module
    # protected-member inspection on the signature and the isinstance guard below.

    # noinspection PyProtectedMember
    def __init__(self, unsteady_problem: problems._CoupledUnsteadyProblem) -> None:
        """The initialization method.

        :param unsteady_problem: The _CoupledUnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(unsteady_problem, problems._CoupledUnsteadyProblem):
            raise TypeError("unsteady_problem must be a _CoupledUnsteadyProblem.")
        super().__init__(unsteady_problem)

    # As in __init__, the reference to the deliberately private _CoupledUnsteadyProblem
    # is by design; silence the cross-module protected-member inspection.

    # noinspection PyProtectedMember
    @property
    def _coupled_unsteady_problem(self) -> problems._CoupledUnsteadyProblem:
        """Type narrowed view of the inherited unsteady_problem attribute.

        The parent stores unsteady_problem as a CoreUnsteadyProblem (widened to let
        subclasses pass their own variants). __init__ validates that this subclass
        always receives a _CoupledUnsteadyProblem, so the cast here is safe.

        :return: The unsteady_problem narrowed to _CoupledUnsteadyProblem.
        """
        return cast(problems._CoupledUnsteadyProblem, self.unsteady_problem)

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
