"""Contains the SteadyProblem and UnsteadyProblem classes.

**Contains the following classes:**

SteadyProblem: A class used to contain steady aerodynamics problems.

UnsteadyProblem: A class used to contain unsteady aerodynamics problems.

**Contains the following functions:**

None
"""

from __future__ import annotations

import math

import numpy as np

from . import geometry
from . import movements

from . import _parameter_validation
from . import _transformations
from . import operating_point as operating_point_mod


class SteadyProblem:
    """A class used to contain steady aerodynamics problems.

    **Contains the following methods:**

    None
    """

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
        if not isinstance(airplanes, list):
            raise TypeError("airplanes must be a list.")
        if len(airplanes) < 1:
            raise ValueError("airplanes must have at least one element.")
        for airplane in airplanes:
            if not isinstance(airplane, geometry.airplane.Airplane):
                raise TypeError("Every element in airplanes must be an Airplane.")
        self.airplanes = airplanes
        if not isinstance(operating_point, operating_point_mod.OperatingPoint):
            raise TypeError("operating_point must be an OperatingPoint.")
        self.operating_point = operating_point

        # Validate that the first Airplane has Cg_GP1_CgP1 set to zeros.
        self.airplanes[0].validate_first_airplane_constraints()

        # Populate GP1_CgP1 coordinates for all Airplanes' Panels This finds the Panels'
        # positions in the first Airplanes' geometry axes, relative to the first
        # Airplanes' CG based on their locally defined positions.
        for airplane_id, airplane in enumerate(self.airplanes):
            # Compute the passive transformation matrix from this Airplane's local
            # geometry axes, relative to its CG, to the first Airplanes' geometry axes,
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


class UnsteadyProblem:
    """A class used to contain unsteady aerodynamics problems.

    **Contains the following methods:**

    None
    """

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
        if not isinstance(movement, movements.movement.Movement):
            raise TypeError("movement must be a Movement.")
        self.movement = movement
        self.only_final_results = _parameter_validation.boolLike_return_bool(
            only_final_results, "only_final_results"
        )

        self.num_steps: int = self.movement.num_steps
        self.delta_time: float = self.movement.delta_time

        # For UnsteadyProblems with a static Movement, we are typically interested in
        # the final time step's forces and moments, which, assuming convergence, will be
        # the most accurate. For UnsteadyProblems with cyclic movement, (e.g. flapping
        # wings) we are typically interested in the forces and moments averaged over the
        # last cycle simulated. Therefore, determine which time step will be the first
        # with relevant results based on if the Movement is static or cyclic.
        _movement_max_period = self.movement.max_period
        self.first_averaging_step: int
        if _movement_max_period == 0:
            self.first_averaging_step = self.num_steps - 1
        else:
            self.first_averaging_step = max(
                0,
                math.floor(self.num_steps - (_movement_max_period / self.delta_time)),
            )

        # If we only wants to calculate forces and moments for the final cycle (for a
        # cyclic Movement) or for the final time step (for a static Movement) set the
        # first step to calculate results to the first averaging step. Otherwise, set it
        # to the zero, which is the first time step.
        self.first_results_step: int
        if self.only_final_results:
            self.first_results_step = self.first_averaging_step
        else:
            self.first_results_step = 0

        # Initialize empty lists to hold the final loads and load coefficients each
        # Airplane experiences. These will only be populated if this UnsteadyProblem's
        # Movement is static.
        self.finalForces_W: list[np.ndarray] = []
        self.finalForceCoefficients_W: list[np.ndarray] = []
        self.finalMoments_W_CgP1: list[np.ndarray] = []
        self.finalMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize empty lists to hold the final cycle-averaged loads and load
        # coefficients each Airplane experiences. These will only be populated if this
        # UnsteadyProblem's Movement is cyclic.
        self.finalMeanForces_W: list[np.ndarray] = []
        self.finalMeanForceCoefficients_W: list[np.ndarray] = []
        self.finalMeanMoments_W_CgP1: list[np.ndarray] = []
        self.finalMeanMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize empty lists to hold the final cycle-root-mean-squared loads and
        # load coefficients each airplane object experiences. These will only be
        # populated for variable geometry problems.
        self.finalRmsForces_W: list[np.ndarray] = []
        self.finalRmsForceCoefficients_W: list[np.ndarray] = []
        self.finalRmsMoments_W_CgP1: list[np.ndarray] = []
        self.finalRmsMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize an empty list to hold the SteadyProblems.
        self.steady_problems: list[SteadyProblem] = []

        # Iterate through the UnsteadyProblem's time steps.
        for step_id in range(self.num_steps):

            # Get the Airplanes and the OperatingPoint associated with this time step.
            these_airplanes = []
            for this_base_airplane in movement.airplanes:
                these_airplanes.append(this_base_airplane[step_id])
            this_operating_point = movement.operating_points[step_id]

            # Initialize the SteadyProblem at this time step.
            this_steady_problem = SteadyProblem(
                airplanes=these_airplanes, operating_point=this_operating_point
            )

            # Append this SteadyProblem to the list of SteadyProblems.
            self.steady_problems.append(this_steady_problem)
