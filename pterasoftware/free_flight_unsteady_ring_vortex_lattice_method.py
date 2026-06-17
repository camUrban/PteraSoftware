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

from . import geometry, problems
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

    **Strongly coupled solve interface:** freeze_substep, evaluate_trial_aero_loads, and
    restore_substep are the public methods the FreeFlightUnsteadyProblem's step solve
    drives to evaluate the aerodynamics at successive trial body states within a single
    time step. They share transient working state (the frozen wake induced velocities
    and the snapshot bound ring vortex strengths) that freeze_substep establishes and
    restore_substep clears.
    """

    __slots__ = (
        "_substep_stackVIndGridWrvp_GP1__E",
        "_substep_gamma_n",
        "_substep_gamma_n_minus_1",
    )

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

        # Transient working state for the strongly coupled sub-iteration, established by
        # freeze_substep and cleared by restore_substep. Each is None between substeps.
        # The frozen induced velocities are the iterate-independent part of the wake
        # transport; the two strength snapshots are the current and previous steps' solved
        # bound ring vortex strengths.
        self._substep_stackVIndGridWrvp_GP1__E: list[list[np.ndarray]] | None = None
        self._substep_gamma_n: np.ndarray | None = None
        self._substep_gamma_n_minus_1: np.ndarray | None = None

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

    def freeze_substep(self) -> None:
        """Freezes the snapshot data the sub-iteration reuses across its trials.

        Called once at the start of a strongly coupled free-flight step solve, before
        any trial. Precomputes the induced (Biot-Savart) part of the wake transport
        velocity, which depends only on the current step's bound geometry, strengths,
        and wake and so is identical across the step's trials (skipped for a prescribed
        wake, which has no induced part), and snapshots the current and previous steps'
        solved bound ring vortex strengths. evaluate_trial_aero_loads reuses these, and
        restore_substep clears them. The solver's current step, Airplanes, and strengths
        must be those of the step being solved, as they are when the run loop invokes
        the step solve.

        The wake-grid induced velocity is recomputed by the official wake build after
        the solve accepts, so its singularity events are logged there; the throwaway
        counters here are not logged.

        :return: None
        """
        if self._prescribed_wake:
            self._substep_stackVIndGridWrvp_GP1__E = None
        else:
            bound_singularity_counts = np.zeros(4, dtype=np.int64)
            wake_singularity_counts = np.zeros(4, dtype=np.int64)
            self._substep_stackVIndGridWrvp_GP1__E = (
                self._calculate_wake_grid_induced_velocities(
                    bound_singularity_counts, wake_singularity_counts
                )
            )

        self._substep_gamma_n = self._current_bound_vortex_strengths.copy()
        self._substep_gamma_n_minus_1 = self._last_bound_vortex_strengths.copy()

    def evaluate_trial_aero_loads(self, step: int) -> geometry.airplane.Airplane:
        """Evaluates the aerodynamic loads at a trial body state for the next time step.

        Called once per sub-iteration, with the trial OperatingPoint already installed
        in the next step's SteadyProblem by the caller. Builds the trial bound geometry
        at the next step, ages the current step's wake into the next step using the
        frozen induced transport from freeze_substep plus the trial frame terms, then
        solves the circulation and loads at the next step. Returns the first Airplane
        carrying the trial aerodynamic loads, which the caller assembles into the
        interval load.

        The wake build runs with the current step as the reference step: its wake is the
        one aged, and its snapshot bound ring vortex strengths set the newly shed row,
        so the current step's context and snapshot strengths are restored before it. The
        circulation and load solve then runs at the next step against the trial geometry
        and the freshly aged wake, with the snapshot strengths supplied as the previous
        strengths for the unsteady load term.

        :param step: The current time step index (zero indexed). The trial body state is
            the next step's.
        :return: The first Airplane, carrying the trial aerodynamic forces and moments
            the solver calculated at the next step.
        """
        assert self._substep_gamma_n is not None

        next_step = step + 1

        # G: build the trial bound ring vortex geometry at the next step from the trial
        # OperatingPoint the caller installed in the next step's SteadyProblem.
        self._initialize_panel_vortices_at(next_step)

        # P: age the current step's wake into the next step. The wake builder reads the
        # current step as the reference step, the next step's OperatingPoint for the trial
        # frame terms, and the snapshot bound strengths for the newly shed row, so set the
        # current step's context and restore the snapshot strengths first. The induced
        # transport is the frozen precompute from freeze_substep.
        current_problem = self._get_steady_problem_at(step)
        self._current_step = step
        self.current_airplanes = current_problem.airplanes
        self.current_operating_point = current_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E
        self._current_bound_vortex_strengths = self._substep_gamma_n
        self._populate_next_airplanes_wake_vortex_points(
            self._substep_stackVIndGridWrvp_GP1__E
        )
        self._populate_next_airplanes_wake_vortices()

        # C and L: solve the circulation and loads at the next step against the trial
        # geometry and the freshly aged wake. The unsteady load term differences the
        # solved strengths against the snapshot strengths.
        next_problem = self._get_steady_problem_at(next_step)
        self._current_step = next_step
        self.current_airplanes = next_problem.airplanes
        self.current_operating_point = next_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E
        self._last_bound_vortex_strengths = self._substep_gamma_n
        self._evaluate_step_aerodynamics()

        return self.current_airplanes[0]

    def restore_substep(self, step: int) -> None:
        """Restores the solver's current-step state after a strongly coupled step solve.

        Called once after the sub-iteration accepts, before the inherited run loop
        builds the official wake for the accepted next step. The trials overwrote the
        current step's solver working state (the geometry stacks, influences, and
        strengths) while evaluating the next step, so this re-runs the current step's
        aerodynamic evaluation to reconstruct that state exactly from the current step's
        intact per step geometry and wake, leaving the solver as if the step had just
        been solved. The previous strengths are restored first so the re-evaluation
        reproduces the original solve. Clears the frozen sub-iteration state.

        :param step: The current time step index (zero indexed).
        :return: None
        """
        assert self._substep_gamma_n_minus_1 is not None

        current_problem = self._get_steady_problem_at(step)
        self._current_step = step
        self.current_airplanes = current_problem.airplanes
        self.current_operating_point = current_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E
        self._last_bound_vortex_strengths = self._substep_gamma_n_minus_1
        self._evaluate_step_aerodynamics()

        self._substep_stackVIndGridWrvp_GP1__E = None
        self._substep_gamma_n = None
        self._substep_gamma_n_minus_1 = None
