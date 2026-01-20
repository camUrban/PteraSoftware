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
import scipy.signal as sp_sig


from copy import deepcopy
from scipy.integrate import quad, solve_ivp
import matplotlib.pyplot as plt
from .movements.single_step.single_step_movement import SingleStepMovement
from . import _parameter_validation, _transformations, geometry, movements
from . import operating_point as operating_point_mod


class SteadyProblem:
    """A class used to contain steady aerodynamics problems.

    **Contains the following methods:**

    reynolds_numbers: A list of Reynolds numbers, one for each Airplane in the
    SteadyProblem.
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

    @property
    def reynolds_numbers(self) -> list[float]:
        """A list of Reynolds numbers, one for each Airplane in the SteadyProblem.

        **Notes:**

        The Reynolds number is calculated as: Re = (V x L) / nu, where V is the
        freestream speed, observed from the Earth frame (vCg__E from OperatingPoint,
        m/s), L is the characteristic length (c_ref from Airplane, m), and nu is the
        kinematic viscosity (nu from OperatingPoint, m^2/s).

        These Reynolds numbers only consider the freestream speed, not any apparent
        velocity due to prescribed motion, so be careful interpreting it for cases where
        this SteadyProblem corresponds to one time step in an UnsteadyProblem.

        :return: A list of Reynolds numbers, one for each Airplane.
        """
        v = self.operating_point.vCg__E
        nu = self.operating_point.nu

        reynolds_list = []
        for airplane in self.airplanes:
            c_ref = airplane.c_ref
            assert c_ref is not None, "Airplane c_ref must be set to calculate Re"
            re = (v * c_ref) / nu
            reynolds_list.append(re)

        return reynolds_list


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
        # last cycle simulated. Use the LCM of all motion periods to ensure we average
        # over a complete cycle of all motions.
        _movement_lcm_period = self.movement.lcm_period
        self.first_averaging_step: int
        if _movement_lcm_period == 0:
            self.first_averaging_step = self.num_steps - 1
        else:
            self.first_averaging_step = max(
                0,
                math.floor(self.num_steps - (_movement_lcm_period / self.delta_time)),
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

class CoupledUnsteadyProblem():
    """This is a class for coupled unsteady problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None
    """

    def __init__(self, movement, only_final_results=False):
        """This is the initialization method.

        :param movement: Movement

            This is the Movement that contains this UnsteadyProblem's
            OperatingPointMovement and AirplaneMovements.

        :param only_final_results: boolLike, optional

            If set to True, the Solver will only calculate forces, moments,
            and pressures for the final complete cycle (of the Movement's
            sub-Movement with the longest period), which increases simulation speed.
            The default value is False.
        """
        if not isinstance(movement, movements.movement.Movement):
            raise TypeError("movement must be a Movement.")
        self.movement = movement
        self.only_final_results = _parameter_validation.boolLike_return_bool(
            only_final_results, "only_final_results"
        )

        self.num_steps = self.movement.num_steps
        self.delta_time = self.movement.delta_time

        # For UnsteadyProblems with a static Movement, users are typically interested
        # in the final time step's forces and moments, which, assuming convergence,
        # will be the most accurate. For UnsteadyProblems with cyclic movement,
        # (e.g. flapping wings) users are typically interested in the forces and
        # moments averaged over the last cycle simulated. Therefore, determine which
        # time step will be the first with relevant results based on if the Movement
        # is static or cyclic.
        _movement_max_period = self.movement.max_period
        if _movement_max_period == 0:
            self.first_averaging_step = self.num_steps - 1
        else:
            self.first_averaging_step = max(
                0,
                math.floor(self.num_steps - (_movement_max_period / self.delta_time)),
            )

        # If the user only wants to calculate forces and moments for the final cycle
        # (for a cyclic Movement) or for the final time step (for a static Movement)
        # set the first step to calculate results to the first averaging step.
        # Otherwise, set it to the zero, which is the first time step.
        if self.only_final_results:
            self.first_results_step = self.first_averaging_step
        else:
            self.first_results_step = 0

        # Initialize empty lists to hold the final loads and load coefficients each
        # Airplane experiences. These will only be populated if this
        # UnsteadyProblem's Movement is static.
        self.finalForces_W = []
        self.finalForceCoefficients_W = []
        self.finalMoments_W_CgP1 = []
        self.finalMomentCoefficients_W_CgP1 = []

        # Initialize empty lists to hold the final cycle-averaged loads and load
        # coefficients each Airplane experiences. These will only be populated if
        # this UnsteadyProblem's Movement is cyclic.
        self.finalMeanForces_W = []
        self.finalMeanForceCoefficients_W = []
        self.finalMeanMoments_W_CgP1 = []
        self.finalMeanMomentCoefficients_W_CgP1 = []

        # Initialize empty lists to hold the final cycle-root-mean-squared loads and
        # load coefficients each airplane object experiences. These will only be
        # populated for variable geometry problems.
        self.finalRmsForces_W = []
        self.finalRmsForceCoefficients_W = []
        self.finalRmsMoments_W_CgP1 = []
        self.finalRmsMomentCoefficients_W_CgP1 = []

        # this set of steady problems should essnetially be treated as private
        # and the getter method should be used to obtain it
        self._steady_problems = [
            SteadyProblem(
                [self.movement.airplane_movements[0].base_airplane],
                self.movement.operating_point_movement.base_operating_point,
            )
        ]

        # Initialize an empty list to hold the SteadyProblems.
        self.steady_problems = []

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

    def get_steady_problem(self, step):
        """
        Return the steady-state problem associated with the given step index.

        Parameters
        ----------
        step : int
            Index of the steady problem to retrieve.

        Returns
        -------
        Any
            The steady-state problem object stored at the specified index.

        Raises
        ------
        Exception
            If `step` is greater than or equal to the number of initialized
            steady problems.
        """
        # Ensure the requested step index is valid.
        if step >= len(self._steady_problems):
            raise Exception(
                f"Step index {step} is out of range of the number of initialized problems"
            )

        # Return the corresponding steady-state problem.
        return self._steady_problems[step]

    def initialize_next_problem(self, solver):
        self._steady_problems.append(self.steady_problems[len(self._steady_problems)])

class AeroelasticUnsteadyProblem(CoupledUnsteadyProblem):
    def __init__(self, single_step_movement: SingleStepMovement, movement, only_final_results=False):
        # TODO: fix this constructor to properly inherit from CoupledUnsteadyProblem
        super().__init__(movement, only_final_results)
        self.prev_velocities = []
        self.single_step_movement = single_step_movement
        self.G = 1e8
        self.I_area = 1e-14
        self.curr_airplanes = [movement.airplane_movements[0].base_airplane]
        self.curr_operating_point = movement.operating_point_movement.base_operating_point
        self.last_torsion_angles = None
        self.error_exceeded_air = False
        self.storage_steady_problem = None
        self.net_torsion_angles = np.zeros(9)

    def define_mass_matrix(self, M_wing, airplane):
        """
        Currently treats all panels as having equal mass. This will

        :param M_wing: float
            This parameter is the total mass of one wing in (kg).
        :param _geometry.airplane.Airplane:
            The current

        :return numpy.ndarray
            A 3D array of shape (num_spanwise_panels, num_chordwise_panels)
            containing a float value for the mass of each panel
        """
        # yes it's bad practice to have this in both functions, but I intend to update
        # this with more complex methods
        num_spanwise_panels = airplane.wings[0].num_spanwise_panels
        num_chordwise_panels = airplane.wings[0].num_chordwise_panels
        point_mass = M_wing / (num_spanwise_panels * num_chordwise_panels)

        return (
            np.ones((num_chordwise_panels, num_spanwise_panels, 3), dtype=float)
            * point_mass
        )

    def calculate_wing_panel_accelerations(self, solver, num_panels):
        dt = self.movement.delta_time

        curr_wing_panel_veloctiy = solver.calculate_solution_velocity(
            solver.stackCpp_GP1_CgP1
        )
        if len(self.prev_velocities) < 1:
            # Set the flapping velocities to be zero for all points. Then, return the
            # flapping velocities.
            return np.zeros((num_panels, 3))

        return (curr_wing_panel_veloctiy - self.prev_velocities[-1])[:num_panels] / dt

    def get_steady_problem(self, step):
        if self.error_exceeded_air:
            return self.storage_steady_problem
        else:
            return super().get_steady_problem(step)

    def initialize_next_problem(self, solver):
        deformation_matrices = self.calculate_wing_deformation(solver, len(self._steady_problems))
        self.curr_airplanes, self.curr_operating_point = (
            self.single_step_movement.generate_next_movement(
                base_airplanes=self.curr_airplanes,
                base_operating_point=self.curr_operating_point,
                step=len(self._steady_problems),
                deformation_matrices=deformation_matrices,
            )
        )
        if self.error_exceeded_air:
            self.storage_steady_problem = SteadyProblem(
                airplanes=self.curr_airplanes, operating_point=self.curr_operating_point
            )
        else:
            self._steady_problems.append(
                SteadyProblem(
                    airplanes=self.curr_airplanes, operating_point=self.curr_operating_point
                )
            )

    def calculate_wing_deformation(self, solver, step):
        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]

        wing = airplane.wings[0]
        mass_matrix = self.define_mass_matrix(0.12, airplane)

        # Panel number definitions
        num_chordwise_panels = wing.num_chordwise_panels
        num_spanwise_panels = wing.num_spanwise_panels
        num_panels = num_chordwise_panels * num_spanwise_panels

        current_torsion_aero = np.zeros(num_chordwise_panels)
        current_torsion_inertia = np.zeros(num_chordwise_panels)
        span_torsion_angles = np.zeros(int(num_spanwise_panels))
        chord_torsion_angles = np.zeros(int(num_spanwise_panels))
        torsion_matrix = np.zeros((num_chordwise_panels, num_spanwise_panels))

        points = np.array(solver.stackCpp_GP1_CgP1)[:num_panels, :]
        x_values = points.reshape((num_chordwise_panels, num_spanwise_panels, 3))[
            :num_panels, :, 0
        ]
        panelAeroForces_G = (
            np.stack([o.forces_GP1 for o in np.ravel(wing.panels)]).reshape(
                (num_chordwise_panels, num_spanwise_panels, 3)
            )
        ) / 10000

        panelInertialForces = self.calculate_wing_panel_accelerations(solver, num_panels).reshape(num_chordwise_panels, num_spanwise_panels, 3) * mass_matrix

        print("\nMaximums", max(panelAeroForces_G.flatten()), max(panelInertialForces.flatten()))
        print("\nMinimums", min(panelAeroForces_G.flatten()), min(panelInertialForces.flatten()))
        # Iterate over spanwise and chordwise panels to find cumulative torsion due to force on each mesh element
        # Force across spanwise panel is distinct
        for span_panel in range(num_spanwise_panels):
            # Force on each chordwise panel from LE to TE
            # from each spanwise point is added to produce torsion at LE
            for chord_panel in range(num_chordwise_panels):
                # Torsion due to UVLM aero forces on LE

                current_torsion_aero[chord_panel] = (
                    panelAeroForces_G[chord_panel][span_panel][2]
                ) * (x_values[chord_panel][span_panel] - x_values[0][span_panel])
                # Torsion due to panel inertia
                current_torsion_inertia[chord_panel] = (
                    panelInertialForces[chord_panel][span_panel][2]
                ) * (x_values[chord_panel][span_panel] - x_values[0][span_panel])
                # Total torsion on LE due to chordwise panel. Aero and Inertial Torsion are oriented in opposite sense
                ct_angle = quad(
                    self.d_alpha_dy_air_static,
                    0,
                    0.01,
                    args=(
                        -(current_torsion_aero[chord_panel])
                        + (current_torsion_inertia[chord_panel]),
                        self.G * self.I_area,
                    ),
                )[0]

                chord_torsion_angles[span_panel] += ct_angle
                torsion_matrix[chord_panel][span_panel] = ct_angle
            # Torsion on span-wise collection of panels
            span_torsion_angles[span_panel] = (
                span_torsion_angles[span_panel - 1] + chord_torsion_angles[span_panel]
                if span_panel > 0
                else chord_torsion_angles[span_panel]
            )

        # Inserting torsion of static span-wise collection of panels at wing root
        span_torsion_angles = np.insert(span_torsion_angles, 0, 0)

        # ### ********************** Convergence of torsion angle ***************************** ###
        # Error in torsion angle (radians)
        # if self.last_torsion_angles is None:
        #     self.last_torsion_angles = np.zeros(num_spanwise_panels + 1)
        # error_torsion = abs(span_torsion_angles - self.last_torsion_angles)

        # # if step < 0.5 * self.num_steps / 3:
        # #     # Error threshold for 1st half cycle is high to account for initial spike in forces in UVLM
        # #     self.error_exceeded_air = False
        # # else:
        # #     # Error threshold in subsequeny timesteps is set to 0.01
        # #     # TODO : Determine appropriate error threshold by running test cases on changing wing twist
        # error_threshold = 0.01 * np.pi / 180  #
        # # Boolean to determine if change in torsion on any panel exceeds error threshold
        # self.error_exceeded_air = np.any(error_torsion[1:] > error_threshold)

        # self.last_torsion_angles = span_torsion_angles
        # TODO: add logic for when to append vs overwrite
        if (True):
            self.prev_velocities.append(
                solver.calculate_solution_velocity(solver.stackCpp_GP1_CgP1)
            )
        span_torsion_angles = span_torsion_angles / -40
        if step > 30:
            self.net_torsion_angles += span_torsion_angles
            print("Net torsion angles (deg): ", self.net_torsion_angles)
        else:
            span_torsion_angles = np.zeros_like(span_torsion_angles)
        return span_torsion_angles

    def d_alpha_dy_air_static(self, y, tau_torsion, GI):
        return (tau_torsion * y) / GI

    def rotational_inertia(self, m, x, theta):
        return (1 / 3) * m * (x / np.cos(theta)) ** 3


class BetterAeroelasticUnsteadyProblem(CoupledUnsteadyProblem):
    def __init__(
        self,
        single_step_movement: SingleStepMovement,
        movement,
        custom_spacing_second_derivative=None,
        only_final_results=False,
    ):
        # TODO: fix this constructor to properly inherit from CoupledUnsteadyProblem
        super().__init__(movement, only_final_results)
        self.prev_velocities = []
        self.single_step_movement = single_step_movement
        self.curr_airplanes = [movement.airplane_movements[0].base_airplane]
        self.curr_operating_point = (
            movement.operating_point_movement.base_operating_point
        )
        self.positions = []
        self.net_deformation = None
        self.angluar_velocities = None

        # Tunable Parameters
        self.wing_density = 0.024  # per unit height kg/m^2
        self.moment_scaling_factor = 1
        self.spring_constant = 0
        self.damping_constant = 1
        self.aero_scaling = 3.0gi
        self.numerical_integration = True # use numerical integration or closed form solution
        self.damping_eps = 1e-3  # critical damping tolerance

        # self.wing_density = 0.012  # per unit height kg/m^2
        # self.moment_scaling_factor = 5
        # self.spring_constant = 1
        # self.aero_scaling = 0.0
        # self.new_integrand = True

        self.per_step_data = []
        self.net_data = []
        self.angluar_velocity_data = []
        self.per_step_inertial = []
        self.per_step_aero = []
        self.per_step_spring = []
        self.base_wing_positions = None
        self.flap_points = []

        # For custom spacing defined in movement.
        self.custom_spacing_second_derivative = custom_spacing_second_derivative

    def calculate_wing_panel_accelerations(self):
        if len(self.positions) <= 2:
            return np.zeros_like(self.positions[0])
        dt = self.movement.delta_time
        return (self.positions[-1] - 2 * self.positions[-2] + self.positions[-3]) / (dt * dt)

    def calculate_mass_matrix(self, wing):
        areas = np.array([[panel.area for panel in row] for row in wing.panels])
        return (
            np.repeat(areas[:, :, None], 3, axis=2)
            * self.wing_density
        )

    def initialize_next_problem(self, solver):

        deformation_matrices = self.calculate_wing_deformation(
            solver, len(self._steady_problems)
        )
        self.curr_airplanes, self.curr_operating_point = (
            self.single_step_movement.generate_next_movement(
                base_airplanes=self.curr_airplanes,
                base_operating_point=self.curr_operating_point,
                step=len(self._steady_problems),
                deformation_matrices=deformation_matrices,
            )
        )
        self._steady_problems.append(
            SteadyProblem(
                airplanes=self.curr_airplanes,
                operating_point=self.curr_operating_point,
            )
        )

    def calculate_wing_deformation(self, solver, step):
        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]

        wing: geometry.wing.Wing = airplane.wings[0]

        # Panel number definitions
        num_chordwise_panels = wing.num_chordwise_panels
        num_spanwise_panels = wing.num_spanwise_panels
        num_panels = num_chordwise_panels * num_spanwise_panels
        if self.net_deformation is None:
            self.net_deformation = np.zeros((num_spanwise_panels + 1, 3))
            self.angluar_velocities = np.zeros((num_spanwise_panels + 1, 3))

        aeroMoments_GP1_Slep = np.array(solver.moments_GP1_Slep[:num_panels]).reshape(
            num_chordwise_panels, num_spanwise_panels, 3
        ) * self.aero_scaling

        self.positions.append(np.array(
            [[panel.Cpp_GP1_CgP1 for panel in row] for row in wing.panels]
        ))

        mass_matrix = self.calculate_mass_matrix(wing)

        inertial_forces = (
            self.calculate_wing_panel_accelerations()
            * mass_matrix
        )

        inertial_moments = np.cross(
            self.positions[-1] - solver.stack_leading_edge_points[:num_panels].reshape((num_chordwise_panels, num_spanwise_panels, 3)), inertial_forces, axis=2
        )

        undeforemed_wing = self.steady_problems[step].airplanes[0].wings[0]
        undeformed_postions = np.array(
            [[panel.Cpp_GP1_CgP1 for panel in row] for row in undeforemed_wing.panels]
        )

        thetas, omegas, spring_moments = self.calculate_spring_moments(
            num_spanwise_panels=num_spanwise_panels,
            wing=wing,
            mass_matrix=mass_matrix,
            aero_moments=aeroMoments_GP1_Slep,
            step=step,
        )
        self.angluar_velocities[:, 1] = omegas
        if self.base_wing_positions is None:
            self.base_wing_positions = np.array(undeformed_postions)

        self.flap_points.append(np.array(undeformed_postions) - self.base_wing_positions)
        self.per_step_inertial.append(inertial_moments.copy())
        self.per_step_aero.append(aeroMoments_GP1_Slep.copy())
        self.per_step_spring.append(spring_moments.copy())

        total_moments = aeroMoments_GP1_Slep - inertial_moments #+ spring_moments
        deformation_moments = total_moments[:, :, 2]  # Z-axis moments

        step_deformation = np.array(
            [
                np.array(
                    [
                        0,
                        # (np.sum(deformation_moments[:, i]) - self.net_deformation[i + 1][1] * self.spring_constant) * self.moment_scaling_factor,
                        thetas[i + 1] * self.moment_scaling_factor,
                        0,
                    ]
                )
                for i in range(num_spanwise_panels)
            ]
        )

        step_deformation = np.insert(step_deformation, 0, np.array([0,0,0]), axis=0)
        if step > 5:
            self.net_deformation = step_deformation

        self.per_step_data.append(step_deformation)
        self.net_data.append(self.net_deformation.copy())
        self.angluar_velocity_data.append(self.angluar_velocities.copy())

        if step % 10  == 3:
            # print("Net deformation: ", self.net_deformation)
            # print("step deformation: ", step_deformation)
            aeroMoments_GP1_Slep = np.array(solver.moments_GP1_Slep[:num_panels]).reshape(num_chordwise_panels, num_spanwise_panels, 3)

            # print("Aero Moments Slep", self.net_deformation[:, 1])
            print("Thetas: ", thetas)

        if step == self.num_steps - 1:
            zero_curve = np.zeros((1, np.array(self.per_step_inertial).shape[0]))
            print(np.array(self.per_step_inertial).shape)
            print(np.array(self.per_step_aero).shape)
            print(np.array(self.per_step_spring).shape)
            print(np.array(self.per_step_data).shape)
            print(np.array(self.net_data).shape)
            plot_curves(np.array(self.per_step_data)[:, :, 1].T.tolist(), "Per Step Deformation")
            plot_curves(np.array(self.net_data)[:, :, 1].T.tolist(), "Net Deformation")
            plot_curves(np.vstack((zero_curve, np.array(self.per_step_inertial)[:, :, :, 2].sum(axis=1).T)).tolist(), "Per Step Inertial Moments")
            plot_curves(np.vstack((zero_curve, np.array(self.per_step_aero)[:, :, :, 2].sum(axis=1).T)).tolist(), "Per Step Aero Moments")
            plot_curves(np.vstack((zero_curve, np.array(self.per_step_spring)[:, :, 2].sum(axis=1).T)).tolist(), "Per Step Spring Moments")
            plot_curves(
                np.vstack(
                    (
                        zero_curve,
                        np.array(self.flap_points)[:, :, :, 2].sum(axis=1).T,
                    )
                ).tolist(),
                "Flap Points Z",
            )

        return self.net_deformation

    def calculate_spring_moments(self, num_spanwise_panels, wing, mass_matrix, aero_moments, step):
        spring_moments = np.zeros((num_spanwise_panels, 3))
        thetas = np.zeros(num_spanwise_panels + 1)
        omegas = np.zeros(num_spanwise_panels + 1)
        d = 0.0 # distance from flapping axis to panel centroid
        for span_panel in range(num_spanwise_panels):
            aero_span_moment = np.sum(aero_moments[:, span_panel, 2])
            if span_panel == 0:
                theta0 = 0.0
                omega0 = 0.0
            else:
                theta0 = self.net_deformation[span_panel][1] 
                omega0 = self.angluar_velocities[span_panel][1]  

            dt = self.movement.delta_time
            mass = mass_matrix[:, span_panel, :].sum()
            # Equation for rotational inertia of rectangular prism about flapping axis
            # Considers two factors, the first is the rotational inertial of a rectangular
            # prism about its centroid, the second is the parallel axis theorem to
            # account for distance from flapping axis to the panel centroid
            L = (
                wing.wing_cross_sections[span_panel].chord
                + wing.wing_cross_sections[span_panel + 1].chord
            ) / 2
            W = np.linalg.norm(wing.panels[0][span_panel].frontLeg_G)
            d += W / 2
            span_I = 1/12 * mass * (L ** 2 + W ** 2)  + mass * (d ** 2) 
            # span_I = d * mass * L
            theta, omega, moment = self.calculate_torsional_spring_moment(
                dt,
                # 1/2 * M * L^2
                # I=mass * (wing.wing_cross_sections[span_panel].chord ** 2) / 2,
                # I= 4/3 * mass * (L ** 2),
                I=span_I,
                theta0=theta0,
                omega0=omega0,
                aero_span_moment=aero_span_moment,
                step=step,
                span_I=span_I,
            )
            d += W / 2
            print("Theta", theta, "Omega", omega, "Moment", moment)
            thetas[span_panel + 1] = theta
            omegas[span_panel + 1] = omega
            spring_moments[span_panel] = np.array([0, moment, 0])

        return thetas, omegas, spring_moments

    def calculate_torsional_spring_moment(
        self, dt, I, theta0, omega0, aero_span_moment, step, span_I, num_steps=2, 
    ):
        k = self.spring_constant
        c = self.damping_constant

        t = np.linspace(dt * (step - 1), dt * step, num_steps)

        if self.numerical_integration:
            # ---- Forced numerical integration ----
            theta, omega = self.spring_numerical_ode(t, k, c, I, theta0, omega0, aero_span_moment, self.generate_inertial_torque_function(span_I))

        else:
            # ---- Closed-form with constant torque ----
            const_tau = aero_span_moment + self.generate_inertial_torque_function(span_I)(t[0])
            theta, omega = self.closed_form_spring_ode(
                t, k, c, I, theta0, omega0, const_tau=const_tau
            )

        # ---- Internal spring-damper moment ----
        spring_moment = -k * theta - c * omega

        # ---- Net moment (optional, depending on sign convention) ----
        # net_moment = spring_moment + tau(t)

        return theta, omega, spring_moment

    def generate_inertial_torque_function(self, span_I):
        """
        Docstring for generate_inertial_torque_function
        
        :param span_I: float
            The rotational inertia of the wing span section about alpha (the flapping axis).
        """
        spacing = (
            self.single_step_movement.airplane_movements[0]
            .wing_movements[0]
            .spacingAngles_Gs_to_Wn_ixyz[0]
        )
        wing_movement = self.single_step_movement.airplane_movements[0].wing_movements[0]
        amp = wing_movement.ampAngles_Gs_to_Wn_ixyz[0]
        b = 2 * np.pi / wing_movement.periodAngles_Gs_to_Wn_ixyz[0]
        h = np.deg2rad(wing_movement.phaseAngles_Gs_to_Wn_ixyz[0])
        if spacing == "sine":
            torque_func = lambda time: -1 * (b ** 2) * np.sin(b * time + h) * amp * span_I
        elif spacing == "uniform":
            raise ValueError("Sawtooth function (uniform spacing) is not differentiable, cannot be used for inertial torque function.")
        elif callable(spacing):
            if self.custom_spacing_second_derivative is not None:
                torque_func = lambda time: self.custom_spacing_second_derivative(time) * span_I
            else:
                raise ValueError("Custom spacing function provided without second derivative function for inertial torque calculation.")

        return torque_func 

    def spring_numerical_ode(self, t, k, c, I, theta0, omega0, aero_torque, inertial_torque_func):
        """ Numerical solution to the spring-damper ODE with arbitrary torque input.
        t: numpy.ndarray    
            Time array.
        k: float
            Spring constant.
        c: float
            Damping constant.
        I: float
            Rotational inertia.
        theta0: float
            Initial angular displacement.
        omega0: float
            Initial angular velocity."""
        def tau(time):
            return aero_torque + inertial_torque_func(time) 
        def ode(time, y):
            theta, omega = y
            return [omega, (tau(time) - c * omega - k * theta) / I]

        sol = solve_ivp(
            ode, (t[0], t[-1]), [theta0, omega0], t_eval=t, rtol=1e-9, atol=1e-12
        )

        theta = sol.y[0][-1]
        omega = sol.y[1][-1]

        return theta, omega

    def closed_form_spring_ode(
        self,
        t,
        k,
        c,
        I,
        theta0,
        omega0,
        const_tau=0.0,
    ):
        """ Closed-form solution to the spring-damper ODE with constant torque input. 
        t: numpy.ndarray    
            Time array.
        k: float
            Spring constant.
        c: float
            Damping constant.
        I: float
            Rotational inertia.
        theta0: float
            Initial angular displacement.
        omega0: float
            Initial angular velocity.
        const_tau: float, optional
            Constant torque input. Default is 0.0.
        """
        # equilibrium shift
        theta_eq = const_tau / k

        theta0s = theta0 - theta_eq
        omega0s = omega0

        omega_n = np.sqrt(k / I)
        zeta = c / (2 * np.sqrt(k * I))

        # ---- Underdamped ----
        if zeta < 1.0 - self.damping_eps:
            omega_d = omega_n * np.sqrt(1 - zeta**2)

            A = theta0s
            B = (omega0s + zeta * omega_n * theta0s) / omega_d

            exp_term = np.exp(-zeta * omega_n * t)

            phi = exp_term * (
                A * np.cos(omega_d * t) +
                B * np.sin(omega_d * t)
            )

            omega = exp_term * (
                -zeta * omega_n * (A * np.cos(omega_d * t) + B * np.sin(omega_d * t))
                - A * omega_d * np.sin(omega_d * t)
                + B * omega_d * np.cos(omega_d * t)
            )

        elif zeta > 1.0 + self.damping_eps:
            # ---- Overdamped ----
            s = np.sqrt(zeta**2 - 1)
            r1 = -omega_n * (zeta - s)
            r2 = -omega_n * (zeta + s)

            A = (omega0s - r2 * theta0s) / (r1 - r2)
            B = theta0s - A

            phi = A * np.exp(r1 * t) + B * np.exp(r2 * t)
            omega = A * r1 * np.exp(r1 * t) + B * r2 * np.exp(r2 * t)

        else:
            # ---- Critically damped (limit form) ----
            A = theta0s
            B = omega0s + omega_n * theta0s

            exp_term = np.exp(-omega_n * t)

            phi = (A + B * t) * exp_term
            omega = (B - omega_n * (A + B * t)) * exp_term

        # shift back to physical angle
        theta = phi + theta_eq

        return theta[-1], omega[-1]


def plot_curves(data, title, flap_cycle=None):
    """
    data: list of lists
          each inner list is a curve
    """
    plt.figure(figsize=(12, 6), dpi=200)

    for i, curve in enumerate(data):
        x = range(len(curve))
        plt.plot(x, curve, label=f"Curve {i}")
    if flap_cycle is not None:
        plt.plot(range(len(flap_cycle)), flap_cycle, label=f"Flap Cycle", color="black")
    plt.xlabel("Step")
    plt.ylabel("Value")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{title.replace(' ', '_')}.png")
    plt.show()
