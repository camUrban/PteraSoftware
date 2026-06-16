"""Contains the FreeFlightUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

FreeFlightUnsteadyRingVortexLatticeMethodSolver: A subclass of
CoupledUnsteadyRingVortexLatticeMethodSolver that solves FreeFlightUnsteadyProblems,
contributing the body angular rate (omega cross r) to the apparent velocity at every
evaluation point so that the inherited unsteady ring vortex lattice method models the
six-degree-of-freedom motion that the coupled MuJoCo dynamics produce.

**Contains the following functions:**

None
"""

from __future__ import annotations

from typing import cast

import numpy as np

from . import problems
from ._coupled_unsteady_ring_vortex_lattice_method import (
    CoupledUnsteadyRingVortexLatticeMethodSolver,
)

# Body and geometry axes differ by a 180-degree rotation about y, so transforming a free
# vector (such as an angular velocity) from the first Airplane's body axes to its geometry
# axes negates the x and z components.
_BP1_TO_GP1_FLIP = np.array([-1.0, 1.0, -1.0], dtype=float)
_BP1_TO_GP1_FLIP.flags.writeable = False


class FreeFlightUnsteadyRingVortexLatticeMethodSolver(
    CoupledUnsteadyRingVortexLatticeMethodSolver
):
    """A subclass of CoupledUnsteadyRingVortexLatticeMethodSolver that solves
    FreeFlightUnsteadyProblems.

    In a FreeFlightUnsteadyProblem the body state at each time step comes from the
    coupled MuJoCo rigid-body dynamics (carried out in FreeFlightUnsteadyProblem's
    initialize_next_problem), so each step's OperatingPoint carries a body angular rate
    (omegas_BP1__E) that the standard solver assumes is zero. This solver contributes
    the apparent velocity from that body rate (omega cross r) at every collocation point
    and bound line vortex leg center by overriding _currentOmegasRad_GP1__E, which the
    inherited velocity calculations feed through _apply_body_rate.

    **Key additions over parent CoupledUnsteadyRingVortexLatticeMethodSolver:** sets
    _models_body_rates to True so the inherited constructor permits non-zero body rates,
    and overrides _currentOmegasRad_GP1__E and _convectionOmegasRad_GP1__E to supply the
    current and next OperatingPoints' body angular rates (in the first Airplane's
    geometry axes, in radians per second), the latter used when convecting the wake to
    the next time step.
    """

    __slots__ = ()

    _models_body_rates = True

    def __init__(
        self,
        free_flight_unsteady_problem: problems.FreeFlightUnsteadyProblem,
    ) -> None:
        """The initialization method.

        :param free_flight_unsteady_problem: The FreeFlightUnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(
            free_flight_unsteady_problem, problems.FreeFlightUnsteadyProblem
        ):
            raise TypeError(
                "free_flight_unsteady_problem must be a FreeFlightUnsteadyProblem."
            )
        super().__init__(free_flight_unsteady_problem)

    @property
    def _free_flight_unsteady_problem(self) -> problems.FreeFlightUnsteadyProblem:
        """The solver's FreeFlightUnsteadyProblem, narrowed from the inherited
        unsteady_problem.

        The inherited unsteady_problem slot is typed as the base CoreUnsteadyProblem so
        the parent solver can hold any coupled problem. This solver's constructor only
        accepts a FreeFlightUnsteadyProblem, so the cast here is safe.

        :return: This solver's FreeFlightUnsteadyProblem.
        """
        return cast(problems.FreeFlightUnsteadyProblem, self.unsteady_problem)

    def _currentOmegasRad_GP1__E(self) -> np.ndarray:
        """Finds the current time step's body angular velocity (in the first Airplane's
        geometry axes, observed from the Earth frame).

        **Notes:**

        The current OperatingPoint stores the body angular rate in the first Airplane's
        body axes, in degrees per second. This method transforms it to the first
        Airplane's geometry axes and converts it to radians per second for the omega
        cross r calculation in _apply_body_rate.

        :return: A (3,) ndarray of floats representing the body angular velocity (in the
            first Airplane's geometry axes, observed from the Earth frame). Its units
            are in radians per second.
        """
        omegas_GP1__E = self.current_operating_point.omegas_BP1__E * _BP1_TO_GP1_FLIP
        return cast(np.ndarray, np.deg2rad(omegas_GP1__E))

    def _convectionOmegasRad_GP1__E(self) -> np.ndarray:
        """Finds the next time step's body angular velocity (in the first Airplane's
        geometry axes, observed from the Earth frame) used to convect the wake.

        **Notes:**

        The wake convects over the interval ending at the next time step, so its
        apparent velocity from body rotation uses the next OperatingPoint's body rate
        rather than the current step's. The next OperatingPoint stores the body angular
        rate in the first Airplane's body axes, in degrees per second. This method
        transforms it to the first Airplane's geometry axes and converts it to radians
        per second for the omega cross r calculation in _apply_body_rate.

        :return: A (3,) ndarray of floats representing the next time step's body angular
            velocity (in the first Airplane's geometry axes, observed from the Earth
            frame). Its units are in radians per second.
        """
        next_operating_point = self._get_steady_problem_at(
            self._current_step + 1
        ).operating_point
        omegas_GP1__E = next_operating_point.omegas_BP1__E * _BP1_TO_GP1_FLIP
        return cast(np.ndarray, np.deg2rad(omegas_GP1__E))
