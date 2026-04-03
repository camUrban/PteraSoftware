"""Contains the SteadyProblem, UnsteadyProblem, and CoupledUnsteadyProblem classes.

**Contains the following classes:**

SteadyProblem: A class used to contain steady aerodynamics problems.

UnsteadyProblem: A class used to contain unsteady aerodynamics problems.

CoupledUnsteadyProblem: A class used to contain coupled unsteady aerodynamics problems
where SteadyProblems are generated dynamically at each time step.

**Contains the following functions:**

None
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import _core, _parameter_validation, _transformations, geometry, movements
from . import operating_point as operating_point_mod

if TYPE_CHECKING:
    from .unsteady_ring_vortex_lattice_method import (
        UnsteadyRingVortexLatticeMethodSolver,
    )


class SteadyProblem:
    """A class used to contain steady aerodynamics problems.

    **Contains the following methods:**

    reynolds_numbers: A tuple of Reynolds numbers, one for each Airplane in the
    SteadyProblem.
    """

    __slots__ = (
        "_airplanes",
        "_operating_point",
        "_reynolds_numbers",
    )

    def __init__(
        self,
        airplanes: list[geometry.airplane.Airplane],
        operating_point: operating_point_mod.OperatingPoint,
    ) -> None:
        """The initialization method.

        :param airplanes: The list of the Airplanes for this SteadyProblem.
        :param operating_point: The OperatingPoint for this SteadyProblem.
        :return: None
        """
        # Validate and store immutable attributes.
        if not isinstance(airplanes, list):
            raise TypeError("airplanes must be a list.")
        if len(airplanes) < 1:
            raise ValueError("airplanes must have at least one element.")
        for airplane in airplanes:
            if not isinstance(airplane, geometry.airplane.Airplane):
                raise TypeError("Every element in airplanes must be an Airplane.")
        # Store as tuple to prevent external mutation via .append(), .pop(), etc.
        self._airplanes: tuple[geometry.airplane.Airplane, ...] = tuple(airplanes)

        if not isinstance(operating_point, operating_point_mod.OperatingPoint):
            raise TypeError("operating_point must be an OperatingPoint.")
        self._operating_point = operating_point

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._reynolds_numbers: tuple[float, ...] | None = None

        # Validate that the first Airplane has Cg_GP1_CgP1 set to zeros.
        self._airplanes[0].validate_first_airplane_constraints()

        # Populate GP1_CgP1 coordinates for all Airplanes' Panels. This finds the
        # Panels' positions in the first Airplane's geometry axes, relative to the
        # first Airplane's CG based on their locally defined positions.
        for airplane in self._airplanes:
            # Compute the passive transformation matrix from this Airplane's local
            # geometry axes, relative to its CG, to the first Airplane's geometry axes,
            # relative to the first Airplane's CG.
            T_pas_G_Cg_to_GP1_CgP1 = airplane.T_pas_G_Cg_to_GP1_CgP1

            for wing in airplane.wings:
                assert wing.panels is not None

                for panel in np.ravel(wing.panels):
                    panel.Frpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Frpp_G_Cg, has_point=True
                    )
                    panel.Flpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Flpp_G_Cg, has_point=True
                    )
                    panel.Blpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Blpp_G_Cg, has_point=True
                    )
                    panel.Brpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Brpp_G_Cg, has_point=True
                    )

    # --- Immutable: read only properties ---
    @property
    def airplanes(self) -> tuple[geometry.airplane.Airplane, ...]:
        return self._airplanes

    @property
    def operating_point(self) -> operating_point_mod.OperatingPoint:
        return self._operating_point

    # --- Immutable derived: manual lazy caching ---
    @property
    def reynolds_numbers(self) -> tuple[float, ...]:
        """A tuple of Reynolds numbers, one for each Airplane in the SteadyProblem.

        **Notes:**

        The Reynolds number is calculated as: Re = (V x L) / nu, where V is the
        freestream speed, observed from the Earth frame (vCg__E from OperatingPoint,
        m/s), L is the characteristic length (c_ref from Airplane, m), and nu is the
        kinematic viscosity (nu from OperatingPoint, m^2/s).

        These Reynolds numbers only consider the freestream speed, not any apparent
        velocity due to prescribed motion, so be careful interpreting it for cases where
        this SteadyProblem corresponds to one time step in an UnsteadyProblem.

        :return: A tuple of Reynolds numbers, one for each Airplane.
        """
        if self._reynolds_numbers is None:
            v = self._operating_point.vCg__E
            nu = self._operating_point.nu

            reynolds_list = []
            for airplane in self._airplanes:
                c_ref = airplane.c_ref
                assert c_ref is not None, "Airplane c_ref must be set to calculate Re"
                re = (v * c_ref) / nu
                reynolds_list.append(re)

            # Store as tuple to prevent external mutation.
            self._reynolds_numbers = tuple(reynolds_list)
        return self._reynolds_numbers


class UnsteadyProblem(_core.CoreUnsteadyProblem):
    """A class used to contain unsteady aerodynamics problems.

    **Contains the following methods:**

    only_final_results: Determines whether the solver will only calculate loads for the
    final time step or final cycle.

    num_steps: The number of time steps.

    delta_time: The time step size in seconds.

    first_averaging_step: The first time step included in cycle averaging.

    first_results_step: The first time step for which loads are calculated.

    max_wake_rows: The maximum chordwise wake rows per Wing.

    movement: The Movement that contains this UnsteadyProblem's OperatingPointMovement
    and AirplaneMovements.

    steady_problems: A tuple of SteadyProblems, one for each time step.
    """

    __slots__ = (
        "_movement",
        "_steady_problems",
    )

    def __init__(
        self,
        movement: movements.movement.Movement,
        only_final_results: bool | np.bool_ = False,
    ) -> None:
        """The initialization method.

        :param movement: The Movement that contains this UnsteadyProblem's
            OperatingPointMovement and AirplaneMovements.
        :param only_final_results: Determines whether the Solver will only calculate
            loads for the final time step (for static Movements) or (for non static
            Movements) for will only calculate loads for the time steps in the final
            complete motion cycle (of the Movement's sub Movement with the longest
            period), which increases simulation speed. Can be a bool or a numpy bool and
            will be converted internally to a bool. The default is False.
        :return: None
        """
        # Validate and store the Movement before calling super().__init__() because
        # the Movement provides the parameters that the core class needs.
        if not isinstance(movement, movements.movement.Movement):
            raise TypeError("movement must be a Movement.")
        self._movement = movement

        # Delegate shared initialization (validation, first_averaging_step computation,
        # load list initialization) to the core class.
        super().__init__(
            only_final_results=only_final_results,
            delta_time=self._movement.delta_time,
            num_steps=self._movement.num_steps,
            max_wake_rows=self._movement.max_wake_rows,
            lcm_period=self._movement.lcm_period,
        )

        # Initialize an empty list to hold the SteadyProblems as they are generated.
        steady_problems_temp: list[SteadyProblem] = []

        # Iterate through the UnsteadyProblem's time steps.
        for step_id in range(self._num_steps):

            # Get the Airplanes and the OperatingPoint associated with this time step.
            these_airplanes = []
            for this_base_airplane in movement.airplanes:
                these_airplanes.append(this_base_airplane[step_id])
            this_operating_point = movement.operating_points[step_id]

            # Initialize the SteadyProblem at this time step.
            this_steady_problem = SteadyProblem(
                airplanes=these_airplanes, operating_point=this_operating_point
            )

            # Append this SteadyProblem to the temporary list.
            steady_problems_temp.append(this_steady_problem)

        # Store as tuple to prevent external mutation via .append(), .pop(), etc.
        self._steady_problems: tuple[SteadyProblem, ...] = tuple(steady_problems_temp)

    # --- Immutable: read only properties ---
    @property
    def movement(self) -> movements.movement.Movement:
        return self._movement

    @property
    def steady_problems(self) -> tuple[SteadyProblem, ...]:
        return self._steady_problems


class CoupledUnsteadyProblem(_core.CoreUnsteadyProblem):
    """A class used to contain coupled unsteady aerodynamics problems.

    Unlike UnsteadyProblem, which pre-computes all SteadyProblems during initialization,
    CoupledUnsteadyProblem generates SteadyProblems dynamically at each time step. This
    enables two-way coupling between the aerodynamic solver and external models (e.g.,
    structural deformation, rigid body dynamics) where the geometry or operating
    conditions at step N + 1 depend on the aerodynamic solution at step N.

    The base implementation generates each step's SteadyProblem from the movement's
    prescribed geometry and operating conditions. Subclasses
    (AeroelasticUnsteadyProblem, FreeFlightUnsteadyProblem) override
    initialize_next_problem to incorporate feedback from the solver.

    **Contains the following methods:**

    movement: The CoreMovement that contains this CoupledUnsteadyProblem's
    OperatingPointMovement and AirplaneMovements.

    steady_problems: A tuple of the SteadyProblems generated so far.

    get_steady_problem: Retrieves the SteadyProblem at a given time step.

    initialize_next_problem: Generates and appends the next time step's SteadyProblem.
    """

    __slots__ = (
        "_movement",
        "_coupled_steady_problems",
    )

    def __init__(
        self,
        movement: _core.CoreMovement,
        only_final_results: bool | np.bool_ = False,
    ) -> None:
        """The initialization method.

        :param movement: The CoreMovement (or subclass) that contains this
            CoupledUnsteadyProblem's OperatingPointMovement and AirplaneMovements. It
            must be an instance of CoreMovement.
        :param only_final_results: Determines whether the solver will only calculate
            loads for the final time step (for static movements) or will only calculate
            loads for the time steps in the final complete motion cycle (for non-static
            movements), which increases simulation speed. Can be a bool or a numpy bool
            and will be converted internally to a bool. The default is False.
        :return: None
        """
        if not isinstance(movement, _core.CoreMovement):
            raise TypeError("movement must be a CoreMovement or one of its subclasses.")
        self._movement = movement

        # Delegate shared initialization (validation, first_averaging_step computation,
        # load list initialization) to the core class.
        super().__init__(
            only_final_results=only_final_results,
            delta_time=self._movement.delta_time,
            num_steps=self._movement.num_steps,
            max_wake_rows=self._movement.max_wake_rows,
            lcm_period=self._movement.lcm_period,
        )

        # Generate the initial SteadyProblem for step 0 from the movement's prescribed
        # geometry and operating conditions.
        initial_airplanes: list[geometry.airplane.Airplane] = []
        for airplane_movement in self._movement.airplane_movements:
            initial_airplanes.append(
                airplane_movement.generate_airplane_at_time_step(
                    0, self._movement.delta_time
                )
            )
        initial_operating_point = self._movement.operating_point_movement.generate_operating_point_at_time_step(
            0, self._movement.delta_time
        )
        initial_problem = SteadyProblem(
            airplanes=initial_airplanes,
            operating_point=initial_operating_point,
        )

        # Initialize the mutable list of SteadyProblems with the initial problem.
        self._coupled_steady_problems: list[SteadyProblem] = [initial_problem]

    # --- Immutable: read only properties ---
    @property
    def movement(self) -> _core.CoreMovement:
        return self._movement

    @property
    def steady_problems(self) -> tuple[SteadyProblem, ...]:
        """A tuple of the SteadyProblems generated so far.

        During simulation, this tuple grows as initialize_next_problem is called at each
        time step. After simulation completes, it contains one SteadyProblem per time
        step.

        :return: A tuple of SteadyProblems.
        """
        return tuple(self._coupled_steady_problems)

    def get_steady_problem(self, step: int) -> SteadyProblem:
        """Retrieves the SteadyProblem at a given time step.

        :param step: An int representing the time step index. It must be non-negative
            and less than the number of SteadyProblems generated so far.
        :return: The SteadyProblem at the given time step.
        """
        if step < 0 or step >= len(self._coupled_steady_problems):
            raise ValueError(
                f"Step index {step} is out of range. Only"
                f" {len(self._coupled_steady_problems)} problems have been"
                " initialized so far."
            )
        return self._coupled_steady_problems[step]

    def initialize_next_problem(
        self, solver: UnsteadyRingVortexLatticeMethodSolver
    ) -> None:
        """Generates and appends the next time step's SteadyProblem.

        The base implementation generates the next SteadyProblem from the movement's
        prescribed geometry and operating conditions. Subclasses override this method to
        incorporate feedback from the solver (e.g., structural deformation angles,
        dynamics-integrated operating conditions).

        :param solver: The solver instance, which provides access to the current
            aerodynamic solution state. Subclasses use this to extract forces, moments,
            or other quantities needed to compute the next step's geometry or operating
            conditions.
        :return: None
        """
        next_step = len(self._coupled_steady_problems)

        next_airplanes: list[geometry.airplane.Airplane] = []
        for airplane_movement in self._movement.airplane_movements:
            next_airplanes.append(
                airplane_movement.generate_airplane_at_time_step(
                    next_step, self._movement.delta_time
                )
            )
        next_operating_point = self._movement.operating_point_movement.generate_operating_point_at_time_step(
            next_step, self._movement.delta_time
        )
        next_problem = SteadyProblem(
            airplanes=next_airplanes,
            operating_point=next_operating_point,
        )
        self._coupled_steady_problems.append(next_problem)
