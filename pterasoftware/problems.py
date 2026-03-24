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

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from .movements.single_step.single_step_movement import SingleStepMovement
from . import _parameter_validation, _transformations, geometry, movements
from . import operating_point as operating_point_mod


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


class UnsteadyProblem:
    """A class used to contain unsteady aerodynamics problems.

    **Contains the following methods:**

    None
    """

    __slots__ = (
        "_movement",
        "_only_final_results",
        "_num_steps",
        "_delta_time",
        "_max_wake_rows",
        "_first_averaging_step",
        "_first_results_step",
        "finalForces_W",
        "finalForceCoefficients_W",
        "finalMoments_W_CgP1",
        "finalMomentCoefficients_W_CgP1",
        "finalMeanForces_W",
        "finalMeanForceCoefficients_W",
        "finalMeanMoments_W_CgP1",
        "finalMeanMomentCoefficients_W_CgP1",
        "finalRmsForces_W",
        "finalRmsForceCoefficients_W",
        "finalRmsMoments_W_CgP1",
        "finalRmsMomentCoefficients_W_CgP1",
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
        # Validate and store immutable attributes.
        if not isinstance(movement, movements.movement.Movement):
            raise TypeError("movement must be a Movement.")
        self._movement = movement
        self._only_final_results = _parameter_validation.boolLike_return_bool(
            only_final_results, "only_final_results"
        )

        self._num_steps: int = self._movement.num_steps
        self._delta_time: float = self._movement.delta_time
        self._max_wake_rows: int | None = self._movement.max_wake_rows

        # For UnsteadyProblems with a static Movement, we are typically interested in
        # the final time step's forces and moments, which, assuming convergence, will be
        # the most accurate. For UnsteadyProblems with cyclic movement, (e.g. flapping
        # wings) we are typically interested in the forces and moments averaged over the
        # last cycle simulated. Use the LCM of all motion periods to ensure we average
        # over a complete cycle of all motions.
        _movement_lcm_period = self._movement.lcm_period
        self._first_averaging_step: int
        if _movement_lcm_period == 0:
            self._first_averaging_step = self._num_steps - 1
        else:
            self._first_averaging_step = max(
                0,
                math.floor(self._num_steps - (_movement_lcm_period / self._delta_time)),
            )

        # If we only wants to calculate forces and moments for the final cycle (for a
        # cyclic Movement) or for the final time step (for a static Movement) set the
        # first step to calculate results to the first averaging step. Otherwise, set it
        # to the zero, which is the first time step.
        self._first_results_step: int
        if self._only_final_results:
            self._first_results_step = self._first_averaging_step
        else:
            self._first_results_step = 0

        # Initialize empty lists to hold the final loads and load coefficients each
        # Airplane experiences. These will only be populated if this UnsteadyProblem's
        # Movement is static. These are mutable and populated by the solver.
        self.finalForces_W: list[np.ndarray] = []
        self.finalForceCoefficients_W: list[np.ndarray] = []
        self.finalMoments_W_CgP1: list[np.ndarray] = []
        self.finalMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize empty lists to hold the final cycle-averaged loads and load
        # coefficients each Airplane experiences. These will only be populated if this
        # UnsteadyProblem's Movement is cyclic. These are mutable and populated by the
        # solver.
        self.finalMeanForces_W: list[np.ndarray] = []
        self.finalMeanForceCoefficients_W: list[np.ndarray] = []
        self.finalMeanMoments_W_CgP1: list[np.ndarray] = []
        self.finalMeanMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize empty lists to hold the final cycle-root-mean-squared loads and
        # load coefficients each airplane object experiences. These will only be
        # populated for variable geometry problems. These are mutable and populated by
        # the solver.
        self.finalRmsForces_W: list[np.ndarray] = []
        self.finalRmsForceCoefficients_W: list[np.ndarray] = []
        self.finalRmsMoments_W_CgP1: list[np.ndarray] = []
        self.finalRmsMomentCoefficients_W_CgP1: list[np.ndarray] = []

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

            # Append this SteadyProblem to the list of SteadyProblems.
            # Append this SteadyProblem to the temporary list.
            steady_problems_temp.append(this_steady_problem)

        # Store as tuple to prevent external mutation via .append(), .pop(), etc.
        self._steady_problems: tuple[SteadyProblem, ...] = tuple(steady_problems_temp)

        # --- Immutable: read only properties ---
    @property
    def movement(self) -> movements.movement.Movement:
        return self._movement

    @property
    def only_final_results(self) -> bool:
        return self._only_final_results

    @property
    def num_steps(self) -> int:
        return self._num_steps

    @property
    def delta_time(self) -> float:
        return self._delta_time

    @property
    def first_averaging_step(self) -> int:
        return self._first_averaging_step

    @property
    def first_results_step(self) -> int:
        return self._first_results_step

    @property
    def max_wake_rows(self) -> int | None:
        return self._max_wake_rows

    @property
    def steady_problems(self) -> tuple[SteadyProblem, ...]:
        return self._steady_problems

class CoupledUnsteadyProblem(UnsteadyProblem):
    """This is a class for coupled unsteady problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None
    """

    def __init__(self, single_step_movement, only_final_results=False):
        """This is the initialization method.

        :param single_step_movement: SingleStepMovement

            This is the SingleStepMovement that contains this CoupledUnsteadyProblem's
            SingleStepOperatingPointMovement and SingleStepAirplaneMovements.
            OperatingPointMovement and AirplaneMovements.

        :param only_final_results: boolLike, optional

            If set to True, the Solver will only calculate forces, moments,
            and pressures for the final complete cycle (of the Movement's
            sub-Movement with the longest period), which increases simulation speed.
            The default value is False.
        """
        if not isinstance(single_step_movement, movements.single_step.single_step_movement.SingleStepMovement):
            raise TypeError("single_step_movement must be a SingleStepMovement.")
        self.single_step_movement = single_step_movement
        self.movement = single_step_movement.corresponding_movement
        self.only_final_results = _parameter_validation.boolLike_return_bool(
            only_final_results, "only_final_results"
        )

        super().__init__(movement=self.movement, only_final_results=self.only_final_results)

        # this set of steady problems should essnetially be treated as private
        # and the getter method should be used to obtain it
        self.coupled_steady_problems = [
            SteadyProblem(
                [self.movement.airplane_movements[0].base_airplane],
                self.movement.operating_point_movement.base_operating_point,
            )
        ]

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
        if step >= len(self.coupled_steady_problems):
            raise Exception(
                f"Step index {step} is out of range of the number of initialized problems"
            )

        # Return the corresponding steady-state problem.
        return self.coupled_steady_problems[step]

    def initialize_next_problem(self, solver):
        self.coupled_steady_problems.append(self.steady_problems[len(self.coupled_steady_problems)])

class AeroelasticUnsteadyProblem(CoupledUnsteadyProblem):

    def __init__(
        self,
        single_step_movement: SingleStepMovement,
        wing_density,
        spring_constant,
        damping_constant,
        aero_scaling=1.0,
        moment_scaling_factor=1.0,
        damping_eps=1e-3,
        plot_flap_cycle=False,
        custom_spacing_second_derivative=None,
        only_final_results=False,
    ):
        super().__init__(single_step_movement=single_step_movement, only_final_results=only_final_results)
        self.plot_flap_cycle = plot_flap_cycle
        self.prev_velocities = []
        self.curr_airplanes = [self.movement.airplane_movements[0].base_airplane]
        self.curr_operating_point = (
            self.movement.operating_point_movement.base_operating_point
        )
        self.positions = []
        self.net_deformation = None
        self.angluar_velocities = None

        # Tunable Parameters
        self.wing_density = wing_density  # per unit height kg/m^2
        self.moment_scaling_factor = moment_scaling_factor
        self.spring_constant = spring_constant
        self.damping_constant = damping_constant
        self.aero_scaling = aero_scaling
        self.numerical_integration = True # use numerical integration or closed form solution
        self.damping_eps = damping_eps  # critical damping tolerance

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
            solver, len(self.coupled_steady_problems)
        )
        self.curr_airplanes, self.curr_operating_point = (
            self.single_step_movement.generate_next_movement(
                base_airplanes=self.curr_airplanes,
                base_operating_point=self.curr_operating_point,
                step=len(self.coupled_steady_problems),
                deformation_matrices=deformation_matrices,
            )
        )
        self.coupled_steady_problems.append(
            SteadyProblem(
                airplanes=self.curr_airplanes,
                operating_point=self.curr_operating_point,
            )
        )

    def calculate_wing_deformation(self, solver, step):
        curr_problem: SteadyProblem = self.coupled_steady_problems[-1]
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

        step_deformation = np.array(
            [
                np.array(
                    [
                        0,
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

        if self.plot_flap_cycle and step == self.num_steps - 1:
            zero_curve = np.zeros((1, np.array(self.per_step_inertial).shape[0]))
            self.plot_flap_cycle_curves(np.array(self.per_step_data)[:, :, 1].T.tolist(), "Per Step Deformation")
            self.plot_flap_cycle_curves(np.array(self.net_data)[:, :, 1].T.tolist(), "Net Deformation")
            self.plot_flap_cycle_curves(np.vstack((zero_curve, np.array(self.per_step_inertial)[:, :, :, 2].sum(axis=1).T)).tolist(), "Per Step Inertial Moments")
            self.plot_flap_cycle_curves(np.vstack((zero_curve, np.array(self.per_step_aero)[:, :, :, 2].sum(axis=1).T)).tolist(), "Per Step Aero Moments")
            self.plot_flap_cycle_curves(np.vstack((zero_curve, np.array(self.per_step_spring)[:, :, 2].sum(axis=1).T)).tolist(), "Per Step Spring Moments")
            self.plot_flap_cycle_curves(
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

        # ---- Forced numerical integration ----
        theta, omega = self.spring_numerical_ode(t, k, c, I, theta0, omega0, aero_span_moment, self.generate_inertial_torque_function(span_I))

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

    def plot_flap_cycle_curves(self, data, title, flap_cycle=None):
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
