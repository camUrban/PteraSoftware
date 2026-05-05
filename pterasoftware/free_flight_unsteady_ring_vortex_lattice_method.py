"""Contains the FreeFlightUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

FreeFlightUnsteadyRingVortexLatticeMethodSolver: A stub subclass of
CoupledUnsteadyRingVortexLatticeMethodSolver for free flight simulations. Not yet
implemented.

**Contains the following functions:**

None
"""

from __future__ import annotations

from . import problems
from ._coupled_unsteady_ring_vortex_lattice_method import (
    CoupledUnsteadyRingVortexLatticeMethodSolver,
)


class FreeFlightUnsteadyRingVortexLatticeMethodSolver(
    CoupledUnsteadyRingVortexLatticeMethodSolver
):
    """A stub subclass of CoupledUnsteadyRingVortexLatticeMethodSolver for free flight
    simulations.

    This solver will handle FreeFlightUnsteadyProblems where the operating point is
    updated dynamically based on computed aerodynamic forces and moments.

    Not yet implemented.
    """

    def __init__(
        self,
        free_flight_unsteady_problem: problems.FreeFlightUnsteadyProblem,
    ) -> None:
        """Initialize the solver for a FreeFlightUnsteadyProblem.

        :param free_flight_unsteady_problem: The FreeFlightUnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(
            free_flight_unsteady_problem, problems.FreeFlightUnsteadyProblem
        ):
            raise TypeError(
                "free_flight_unsteady_problem must be a " "FreeFlightUnsteadyProblem."
            )

        super().__init__(free_flight_unsteady_problem)
