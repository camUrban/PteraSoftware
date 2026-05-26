"""Contains the SteadyProblem and UnsteadyProblem classes.

**Contains the following classes:**

SteadyProblem: A class used to contain steady aerodynamics problems.

UnsteadyProblem: A class used to contain unsteady aerodynamics problems.

**Contains the following functions:**

None
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from . import _core, _parameter_validation, _transformations, geometry, movements
from . import operating_point as operating_point_mod
from .movements import aeroelastic_movement as aeroelastic_movement_mod

if TYPE_CHECKING:
    from ._coupled_unsteady_ring_vortex_lattice_method import (
        CoupledUnsteadyRingVortexLatticeMethodSolver,
    )
    from .aeroelastic_unsteady_ring_vortex_lattice_method import (
        AeroelasticUnsteadyRingVortexLatticeMethodSolver,
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


class _CoupledUnsteadyProblem(_core.CoreUnsteadyProblem):
    """A class for coupled unsteady aerodynamics problems.

    This class extends CoreUnsteadyProblem to manage SteadyProblems for coupled
    simulations where the geometry at each time step depends on the solver's results
    from previous time steps.

    **Contains the following methods:**

    movement: The CoreMovement that defines the motion parameters for this problem.

    steady_problems: A tuple of SteadyProblems, one for each time step that has been
    initialized so far.

    get_steady_problem: Gets the SteadyProblem at a specified time step.

    initialize_next_problem: Initializes the next time step's SteadyProblem. Must be
    overridden by subclasses.
    """

    __slots__ = (
        "_movement",
        "_steady_problems",
    )

    def __init__(
        self,
        movement: _core.CoreMovement,
        initial_airplanes: list[geometry.airplane.Airplane],
        initial_operating_point: operating_point_mod.OperatingPoint,
    ) -> None:
        """The initialization method.

        Initializes the coupled unsteady problem with the first time step's geometry and
        the motion parameters from the provided CoreMovement.

        :param movement: A CoreMovement object that defines the motion parameters
            (delta_time, num_steps, max_wake_rows, lcm_period) for this problem.
        :param initial_airplanes: The list of Airplanes at the first time step.
        :param initial_operating_point: The OperatingPoint at the first time step.
        :return: None
        """
        self._movement = movement

        # Delegate shared initialization (validation, first_averaging_step computation,
        # load list initialization) to the core class. _CoupledUnsteadyProblems require
        # per step results to feed the coupling hook, so only_final_results is always
        # False.
        super().__init__(
            only_final_results=False,
            delta_time=self._movement.delta_time,
            num_steps=self._movement.num_steps,
            max_wake_rows=self._movement.max_wake_rows,
            lcm_period=self._movement.lcm_period,
        )

        # Coupled-specific state: a mutable list of SteadyProblems that grows as the
        # solver advances. Subclass initialize_next_problem overrides append to this
        # list; external code reads through the steady_problems tuple-view property to
        # preserve the read-only contract inherited from UnsteadyProblem. Seed with a
        # SteadyProblem built from the initial geometry so step zero is always ready.
        self._steady_problems: list[SteadyProblem] = [
            SteadyProblem(
                airplanes=initial_airplanes,
                operating_point=initial_operating_point,
            )
        ]

    # --- Immutable: read only properties ---
    @property
    def movement(self) -> _core.CoreMovement:
        return self._movement

    @property
    def steady_problems(self) -> tuple[SteadyProblem, ...]:
        return tuple(self._steady_problems)

    def get_steady_problem(self, step: int) -> SteadyProblem:
        """Get the SteadyProblem at a given time step.

        :param step: The time step index (zero indexed). Must be greater than or equal
            to zero and less than the total number of time steps.
        :return: The SteadyProblem at the specified time step.
        """
        step = _parameter_validation.int_in_range_return_int(
            step, "step", 0, True, len(self._steady_problems), False
        )

        return self._steady_problems[step]

    def initialize_next_problem(
        self, solver: CoupledUnsteadyRingVortexLatticeMethodSolver
    ) -> None:
        """Initialize the next time step's SteadyProblem.

        Must be overridden by subclasses to compute the geometry for the next time step
        based on the solver's results.

        :param solver: The solver instance providing aerodynamic data from the current
            time step.
        :return: None
        :raises NotImplementedError: Always. Subclasses must override this method.
        """
        raise NotImplementedError("Subclasses must implement initialize_next_problem.")


class AeroelasticUnsteadyProblem(_CoupledUnsteadyProblem):
    """A subclass of _CoupledUnsteadyProblem used to couple aeroelastic wing
    deformations with unsteady aerodynamics.

    This class couples aerodynamic loads with wing structural dynamics (spring-mass-
    damper system) to simulate aeroelastic deformation. Each time step, wing
    deformations are calculated based on the combined effects of aerodynamic moments,
    inertial forces, and spring-damper restoring forces.

    **Contains the following methods:**

    calculate_wing_panel_accelerations: Computes panel accelerations from finite
    difference of positions.

    calculate_mass_matrix: Generates the mass distribution matrix for wing panels.

    calculate_wing_deformation: Computes cumulative wing deformation for the current
    step.

    calculate_spring_moments: Calculates spring-damper moments acting on each spanwise
    section.

    calculate_torsional_spring_moment: Solves the torsional spring-damper ODE for a
    single span section.

    generate_inertial_torque_function: Creates a torque function from prescribed wing
    motion.

    spring_numerical_ode: Numerically integrates the spring-damper differential
    equation.

    plot_flap_cycle_curves: Visualizes moment and deformation time histories.

    **Notes:**

    The aeroelastic coupling assumes a torsional spring-mass-damper model for each
    spanwise section. Wing motion is prescribed through wing flapping, and aerodynamic
    moments from the solver are combined with inertial and spring restoring forces via
    ODE integration to produce structural deformations.
    """

    def __init__(
        self,
        movement: aeroelastic_movement_mod.AeroelasticMovement,
        wing_density: float,
        spring_constant: float,
        damping_constant: float,
        aero_scaling: float = 1.0,
        moment_scaling_factor: float = 1.0,
        damping_eps: float = 1e-3,
        plot_flap_cycle: bool = False,
        custom_spacing_second_derivative=None,
    ) -> None:
        """The initialization method.

        Sets up the aeroelastic problem with structural parameters for the torsional
        spring-mass-damper model applied to each wing spanwise section. Initializes
        storage for aerodynamic loads, deformations, moments, and solver state.

        See _CoupledUnsteadyProblem's initialization method for descriptions of
        inherited parameters.

        :param movement: An AeroelasticMovement object containing the prescribed motion
            and aerodynamic setup for the aeroelastic simulation.
        :param wing_density: The mass per unit span area of the wing (kg/m^2). Used to
            distribute wing mass across panels for inertial calculations.
        :param spring_constant: The torsional spring stiffness for the spring-mass-
            damper model (N*m/rad). Controls the restoring torque opposing deformation.
        :param damping_constant: The torsional damping coefficient (N*m*s/rad). Controls
            the viscous damping in the spring-mass-damper system.
        :param aero_scaling: A scaling factor applied to aerodynamic moments (unitless).
            The default is 1.0. Use values less than 1 to reduce aerodynamic influence.
        :param moment_scaling_factor: A scaling factor applied to the computed wing
            deformation angles (unitless). The default is 1.0. Useful for adjusting the
            magnitude of structural response.
        :param damping_eps: The critical damping tolerance used for diagnostics
            (unitless). The default is 1e-3. This parameter is not currently used in the
            solver.
        :param plot_flap_cycle: If True, plots time histories of moments and
            deformations at the end of the simulation. The default is False.
        :param custom_spacing_second_derivative: An optional callable function of time
            that returns the second time derivative of a custom wing motion spacing
            function. Required if custom (non-sinusoidal) wing motion spacing is used.
            The default is None.
        :return: None
        """
        if not isinstance(movement, aeroelastic_movement_mod.AeroelasticMovement):
            raise TypeError("movement must be an AeroelasticMovement.")

        # Generate the initial airplane at step 0 with no deformation.
        initial_airplane = movement.generate_airplane_at_time_step(
            airplane_movement_index=0, step=0
        )

        super().__init__(
            movement=movement,
            initial_airplanes=[initial_airplane],
            initial_operating_point=movement.operating_points[0],
        )

        # Store a typed reference for aeroelastic-specific operations.
        self._aeroelastic_movement: aeroelastic_movement_mod.AeroelasticMovement = (
            movement
        )

        self.plot_flap_cycle = plot_flap_cycle
        self.prev_velocities: list[np.ndarray] = []
        self.positions: list[np.ndarray] = []
        self.net_deformation: np.ndarray = np.zeros((0, 3))
        self.angular_velocities: np.ndarray = np.zeros((0, 3))

        # Tunable Parameters
        self.wing_density = wing_density  # per unit height kg/m^2
        self.moment_scaling_factor = moment_scaling_factor
        self.spring_constant = spring_constant
        self.damping_constant = damping_constant
        self.aero_scaling = aero_scaling
        self.damping_eps = damping_eps  # critical damping tolerance

        # Permanent parameters
        self.step_discards = (
            5  # number of initial steps to discard for numerical stability
        )
        self.spacing = (
            self._aeroelastic_movement.airplane_movements[0]
            .wing_movements[0]
            .spacingAngles_Gs_to_Wn_ixyz[0]
        )
        self.wing_movement = self._aeroelastic_movement.airplane_movements[
            0
        ].wing_movements[0]

        self.net_data: list[np.ndarray] = []
        self.angular_velocity_data: list[np.ndarray] = []
        self.per_step_inertial: list[np.ndarray] = []
        self.per_step_aero: list[np.ndarray] = []
        self.base_wing_positions: np.ndarray = np.zeros(0)
        self.flap_points: list[np.ndarray] = []

        # For custom spacing defined in movement.
        self.custom_spacing_second_derivative = custom_spacing_second_derivative

    def calculate_wing_panel_accelerations(self) -> np.ndarray:
        """Compute panel accelerations using finite difference of stored positions.

        Calculates second-order accelerations using the finite difference formula: a =
        (p[n] - 2*p[n-1] + p[n-2]) / dt^2.

        :return: An (N_chordwise, N_spanwise, 3) ndarray of floats representing panel
            center accelerations in the global frame. Returns zeros if fewer than 3
            position snapshots are available.
        """
        if len(self.positions) <= 2:
            if len(self.positions) == 0:
                return np.zeros(1)
            return np.zeros_like(self.positions[0])
        dt = self.movement.delta_time
        # If given a relatively large dt value, the finite difference calculation can produce
        # very large accelerations that cause numerical instability in the spring ODE integration.
        # A higher order model may be useful if this is the case.
        pos_m1: np.ndarray = self.positions[-1]
        pos_m2: np.ndarray = self.positions[-2]
        pos_m3: np.ndarray = self.positions[-3]
        return np.array((pos_m1 - 2 * pos_m2 + pos_m3) / (dt * dt))

    def calculate_mass_matrix(self, wing: geometry.wing.Wing) -> np.ndarray:
        """Generate the mass distribution matrix for all wing panels.

        Distributes the total spanwise mass (wing_density) across panel areas to form a
        panel-by-panel mass matrix. Each panel's mass is proportional to its area times
        the specified wing_density.

        :param wing: A Wing object whose panels define the mass distribution.
        :return: An (N_chordwise, N_spanwise, 3) ndarray of floats representing the mass
            at each panel. The three components are identical (mass scalar replicated
            for x, y, z axes).
        """
        assert wing.panels is not None
        areas = np.array([[panel.area for panel in row] for row in wing.panels])
        return np.repeat(areas[:, :, None], 3, axis=2) * self.wing_density

    def initialize_next_problem(
        self, solver: CoupledUnsteadyRingVortexLatticeMethodSolver
    ) -> None:
        # Circular at module level: the aeroelastic solver imports problems.py.
        # Needed at runtime for cast().
        from .aeroelastic_unsteady_ring_vortex_lattice_method import (
            AeroelasticUnsteadyRingVortexLatticeMethodSolver,
        )

        aeroelastic_solver = cast(
            AeroelasticUnsteadyRingVortexLatticeMethodSolver, solver
        )

        step = len(self._steady_problems)
        deformation_matrices = self.calculate_wing_deformation(aeroelastic_solver, step)

        # Build the per-wing deformation list. The main wing (index 0) and its
        # symmetric reflection (index 1) both receive the same deformation angles.
        # The mirror mesh generation in the Airplane constructor ensures that
        # applying the same angles produces physically symmetric deformation.
        # All other wings (e.g. the v-tail) get None (no deformation applied).
        num_wings = len(self._aeroelastic_movement.airplane_movements[0].wing_movements)
        wing_deformation_angles_ixyz: list[np.ndarray | None] = [None] * num_wings
        wing_deformation_angles_ixyz[0] = deformation_matrices
        wing_deformation_angles_ixyz[1] = deformation_matrices

        # Generate the deformed airplane at this step.
        airplane = self._aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0,
            step=step,
            wing_deformation_angles_ixyz=wing_deformation_angles_ixyz,
        )
        operating_point = self._aeroelastic_movement.operating_points[step]

        self._steady_problems.append(
            SteadyProblem(
                airplanes=[airplane],
                operating_point=operating_point,
            )
        )

    def calculate_wing_deformation(
        self,
        solver: AeroelasticUnsteadyRingVortexLatticeMethodSolver,
        step: int,
    ) -> np.ndarray:
        """Compute cumulative wing deformation for the current time step.

        Orchestrates the calculation of inertial moments, spring moments, and cumulative
        deformation. Updates internal state and optionally generates plots.

        :param solver: The solver instance providing aerodynamic moment data
            (moments_GP1_Slep).
        :param step: The current time step index (0-indexed).
        :return: An (N_spanwise+1, 3) ndarray of floats representing cumulative
            deformation angles at each spanwise station. The y-component (index 1)
            contains torsional angles in radians; x and z components are zero.
        """
        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]
        wing: geometry.wing.Wing = airplane.wings[0]

        # Compute panel parameters and mass matrix once
        num_chordwise_panels = wing.num_chordwise_panels
        num_spanwise_panels = wing.num_spanwise_panels
        assert num_spanwise_panels is not None, "num_spanwise_panels must not be None"
        num_panels = num_chordwise_panels * num_spanwise_panels
        mass_matrix = self.calculate_mass_matrix(wing)

        # Initialize deformation state if needed
        if self.net_deformation.size == 0:
            self.net_deformation = np.zeros((num_spanwise_panels + 1, 3))
            self.angular_velocities = np.zeros((num_spanwise_panels + 1, 3))

        # Extract aerodynamic and inertial moments
        aeroMoments_GP1_Slep = self._extract_aero_moments(
            solver, num_chordwise_panels, num_spanwise_panels, num_panels
        )
        inertial_moments = self._calculate_inertial_moments(
            solver,
            wing,
            mass_matrix,
            num_chordwise_panels,
            num_spanwise_panels,
            num_panels,
        )

        # Calculate spring moments and deformation via ODE integration
        thetas, omegas, spring_moments = self.calculate_spring_moments(
            num_spanwise_panels=num_spanwise_panels,
            wing=wing,
            mass_matrix=mass_matrix,
            aero_moments=aeroMoments_GP1_Slep,
            step=step,
        )

        # Build deformation vector and update state
        step_deformation = self._build_deformation_vector(thetas, num_spanwise_panels)
        self._apply_moment_updates(
            step=step,
            step_deformation=step_deformation,
            omegas=omegas,
            inertial_moments=inertial_moments,
            aeroMoments_GP1_Slep=aeroMoments_GP1_Slep,
            spring_moments=spring_moments,
        )

        # Plot results at end of simulation if enabled
        if self.plot_flap_cycle and step == self.num_steps - 1:
            self._plot_aeroelastic_results()

        return self.net_deformation

    def _extract_aero_moments(
        self,
        solver: AeroelasticUnsteadyRingVortexLatticeMethodSolver,
        num_chordwise_panels: int,
        num_spanwise_panels: int,
        num_panels: int,
    ) -> np.ndarray:
        """Extract and scale aerodynamic moments from the solver output.

        Uses the strip leading edge points as the reference point for moment
        calculations, consistent with the assumption of a torsional spring at the
        leading edge.

        :param solver: The solver instance with moments_GP1_Slep data.
        :param num_chordwise_panels: Number of chordwise panel rows.
        :param num_spanwise_panels: Number of spanwise panel rows.
        :param num_panels: Total number of panels (num_chordwise * num_spanwise).
        :return: An (N_chordwise, N_spanwise, 3) ndarray of scaled aerodynamic moments
            in the global panel frame.
        """
        aeroMoments_GP1_Slep = (
            np.array(solver.moments_GP1_Slep[:num_panels]).reshape(
                num_chordwise_panels, num_spanwise_panels, 3
            )
            * self.aero_scaling
        )
        return aeroMoments_GP1_Slep

    def _calculate_inertial_moments(
        self,
        solver: AeroelasticUnsteadyRingVortexLatticeMethodSolver,
        wing: geometry.wing.Wing,
        mass_matrix: np.ndarray,
        num_chordwise_panels: int,
        num_spanwise_panels: int,
        num_panels: int,
    ) -> np.ndarray:
        """Calculate inertial moments from panel accelerations and mass distribution.

        Computes panel accelerations via finite difference, multiplies by mass to get
        forces, then calculates moments about the leading edge reference point using
        cross products.

        :param solver: The solver instance providing leading edge point positions.
        :param wing: The Wing object containing panel definitions.
        :param mass_matrix: An (N_chordwise, N_spanwise, 3) ndarray of panel masses.
        :param num_chordwise_panels: Number of chordwise panel rows.
        :param num_spanwise_panels: Number of spanwise panel rows.
        :param num_panels: Total number of panels (num_chordwise * num_spanwise).
        :return: An (N_chordwise, N_spanwise, 3) ndarray of inertial moment vectors.
        """
        # Store current panel center positions
        assert wing.panels is not None
        self.positions.append(
            np.array([[panel.Cpp_GP1_CgP1 for panel in row] for row in wing.panels])
        )

        # Calculate panel accelerations and inertial forces
        inertial_forces = self.calculate_wing_panel_accelerations() * mass_matrix

        # Calculate moments about leading edge points via cross product
        inertial_moments = np.cross(
            self.positions[-1]
            - solver.stack_leading_edge_points[:num_panels].reshape(
                (num_chordwise_panels, num_spanwise_panels, 3)
            ),
            inertial_forces,
            axis=2,
        )
        return np.array(inertial_moments)

    def _build_deformation_vector(
        self, thetas: np.ndarray, num_spanwise_panels: int
    ) -> np.ndarray:
        """Construct the step deformation vector from torsional angles.

        Converts the torsional angles output from the spring-damper ODE (one per
        spanwise section) into a full (N_spanwise+1, 3) deformation vector with scaling
        applied to the y-component (torsional angle).

        :param thetas: An (N_spanwise+1,) ndarray of torsional angles in radians.
        :param num_spanwise_panels: Number of spanwise panel rows.
        :return: An (N_spanwise+1, 3) ndarray with zero-valued x and z components and
            scaled torsional angles in the y component.
        """
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
        step_deformation = np.insert(step_deformation, 0, np.array([0, 0, 0]), axis=0)
        return step_deformation

    def _apply_moment_updates(
        self,
        step: int,
        step_deformation: np.ndarray,
        omegas: np.ndarray,
        inertial_moments: np.ndarray,
        aeroMoments_GP1_Slep: np.ndarray,
        spring_moments: np.ndarray,
    ) -> None:
        """Update internal moment and deformation state arrays.

        Stores per-step moment and deformation data, updates the cumulative net
        deformation (with discarding of early unstable steps), and tracks wing
        deflection points relative to the undeformed baseline.

        :param step: The current time step index.
        :param step_deformation: The (N_spanwise+1, 3) deformation vector for this step.
        :param omegas: An (N_spanwise+1,) ndarray of angular velocities.
        :param inertial_moments: An (N_chordwise, N_spanwise, 3) ndarray of inertial
            moments.
        :param aeroMoments_GP1_Slep: An (N_chordwise, N_spanwise, 3) ndarray of aero
            moments.
        :param spring_moments: An (N_spanwise, 3) ndarray of spring-damper moments.
        :return: None
        """
        # Update angular velocity state
        self.angular_velocities[:, 1] = omegas

        # Generate the reference (undeformed) airplane at this step to get the
        # baseline panel positions for tracking wing deflection.
        ref_airplane = self._aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0, step=step
        )
        ref_problem = SteadyProblem(
            [ref_airplane], self._aeroelastic_movement.operating_points[step]
        )
        undeformed_wing = ref_problem.airplanes[0].wings[0]
        assert undeformed_wing.panels is not None
        undeformed_positions = np.array(
            [[panel.Cpp_GP1_CgP1 for panel in row] for row in undeformed_wing.panels]
        )
        if self.base_wing_positions.size == 0:
            self.base_wing_positions = np.array(undeformed_positions)

        # Track wing deflection relative to undeformed baseline
        self.flap_points.append(
            np.array(undeformed_positions) - self.base_wing_positions
        )

        # Store per-step moment components for later analysis/plotting
        self.per_step_inertial.append(inertial_moments.copy())
        self.per_step_aero.append(aeroMoments_GP1_Slep.copy())

        # Update cumulative deformation (with numerical stability discarding)
        # Accounts for numerical instability causing large aerodynamic forces in initial steps
        if step > self.step_discards:
            self.net_deformation = step_deformation

        # Store deformation and angular velocity history
        self.net_data.append(self.net_deformation.copy())
        self.angular_velocity_data.append(self.angular_velocities.copy())

    def _plot_aeroelastic_results(self) -> None:
        """Generate and display time-history plots of aeroelastic results.

        Creates plots of per-step and cumulative deformations, moment components
        (inertial, aerodynamic, spring), and wing deflection points. Useful for
        visualizing the aeroelastic coupling behavior.

        :return: None
        """
        zero_curve = np.zeros((1, np.array(self.per_step_inertial).shape[0]))

        # Deformation time histories
        self.plot_flap_cycle_curves(
            np.array(self.net_data)[:, :, 1].T.tolist(), "Net Deformation"
        )

        # Moment component time histories
        self.plot_flap_cycle_curves(
            np.vstack(
                (
                    zero_curve,
                    np.array(self.per_step_inertial)[:, :, :, 2].sum(axis=1).T,
                )
            ).tolist(),
            "Per Step Inertial Moments",
        )
        self.plot_flap_cycle_curves(
            np.vstack(
                (zero_curve, np.array(self.per_step_aero)[:, :, :, 2].sum(axis=1).T)
            ).tolist(),
            "Per Step Aero Moments",
        )

        # Wing deflection tracking
        self.plot_flap_cycle_curves(
            np.vstack(
                (
                    zero_curve,
                    np.array(self.flap_points)[:, :, :, 2].sum(axis=1).T,
                )
            ).tolist(),
            "Flap Points Z",
        )

    def calculate_spring_moments(
        self,
        num_spanwise_panels: int,
        wing: geometry.wing.Wing,
        mass_matrix: np.ndarray,
        aero_moments: np.ndarray,
        step: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate spring-damper moments and angular states for each spanwise section.

        Solves the torsional spring-damper ODE independently for each spanwise section,
        accounting for aerodynamic moments, inertial forces, and structural properties.
        Uses the parallel axis theorem to compute rotational inertia about the flapping
        axis.

        :param num_spanwise_panels: Number of spanwise panel rows in the wing.
        :param wing: The Wing object containing geometric and structural definitions.
        :param mass_matrix: An (N_chordwise, N_spanwise, 3) ndarray of panel masses.
        :param aero_moments: An (N_chordwise, N_spanwise, 3) ndarray of aerodynamic
            moments from the aerodynamic solver.
        :param step: The current time step index.
        :return: A tuple of three ndarrays: - thetas: (N_spanwise+1,) ndarray of
            torsional angles (radians) at each station. - omegas: (N_spanwise+1,)
            ndarray of angular velocities (rad/s) at each station. - spring_moments:
            (N_spanwise, 3) ndarray of spring-damper moment vectors. **Notes:** The
            rotational inertia is computed as: I = (1/12)*M*(L^2 + W^2) + M*d^2, where M
            is panel mass, L is chord, W is span width, and d is distance from the
            flapping axis (computed cumulatively using the parallel axis theorem).
        """
        spring_moments = np.zeros((num_spanwise_panels, 3))
        thetas = np.zeros(num_spanwise_panels + 1)
        omegas = np.zeros(num_spanwise_panels + 1)
        d = 0.0  # distance from flapping axis to panel centroid (computed in half-span increments)
        for span_panel in range(num_spanwise_panels):
            aero_span_moment = np.sum(aero_moments[:, span_panel, 2])
            theta0: float = 0.0
            omega0: float = 0.0
            if span_panel != 0:
                theta0 = self.net_deformation[span_panel][1]
                omega0 = self.angular_velocities[span_panel][1]

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
            assert wing.panels is not None
            W: float = float(np.linalg.norm(wing.panels[0][span_panel].frontLeg_G))
            d += W / 2
            span_I = 1 / 12 * mass * (L**2 + W**2) + mass * (d**2)
            theta, omega, moment = self.calculate_torsional_spring_moment(
                dt,
                # A potential knob to tweak in representation of the torsional inertia
                # I=mass * (wing.wing_cross_sections[span_panel].chord ** 2) / 2,
                I=1 / 2 * mass * (L**2),
                # I=span_I,
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
        self,
        dt: float,
        I: float,
        theta0: float,
        omega0: float,
        aero_span_moment: float,
        step: int,
        span_I: float,
        num_steps: int = 2,
    ) -> tuple[float, float, float]:
        """Solve the torsional spring-damper ODE for a single wing section.

        Integrates the forced torsional damped harmonic oscillator equation: I*d(omega)/dt =
        tau_aero + tau_inertial - k*theta - c*omega

        Returns the angular displacement and velocity at the end of the time step, along
        with the spring-damper restoring moment.

        :param dt: The time step duration (seconds).
        :param I: The rotational inertia about the flapping axis (kg*m^2).
        :param theta0: Initial torsional angle at the start of the time step (radians).
        :param omega0: Initial angular velocity at the start of the time step (rad/s).
        :param aero_span_moment: The z-component aerodynamic moment summed over
            chordwise panels for this spanwise section (N*m).
        :param step: The current time step index (used for inertial torque evaluation).
        :param span_I: The rotational inertia including parallel axis theorem (kg*m^2).
            This is the actual inertia used in the ODE solver.
        :param num_steps: Number of time sub-steps for numerical integration. The
            default is 2.
        :return: A tuple of (theta, omega, spring_moment) where: - theta: Final
            torsional angle (radians). - omega: Final angular velocity (rad/s). -
            spring_moment: The z-component spring-damper moment tau = -k*theta - c*omega (N*m).
        """
        k = self.spring_constant
        c = self.damping_constant
        t = np.linspace(dt * (step - 1), dt * step, num_steps)

        # Forced numerical integration of the spring-damper ODE
        theta, omega = self.spring_numerical_ode(
            t,
            k,
            c,
            I,
            theta0,
            omega0,
            aero_span_moment,
            self.generate_inertial_torque_function(span_I),
        )

        # Internal spring-damper moment (restoring force from structural springs/dampers)
        spring_moment = -k * theta - c * omega

        return theta, omega, spring_moment

    def generate_inertial_torque_function(self, span_I: float):
        """Generate the prescribed wing motion inertial torque function.

        Extracts the prescribed flapping motion from the wing_movement definition and
        creates a callable inertial torque function tau_inertial = I * d^2(theta_prescribed)/dt^2.
        Supports sinusoidal and custom spacing functions.

        :param span_I: The rotational inertia of the wing span section about the
            flapping axis (kg*m^2).
        :return: A callable function that accepts time and returns the inertial torque
            (N*m) due to the prescribed wing motion acceleration. **Notes:** For
            sinusoidal spacing: tau = -I * b^2 * sin(b*t + h) * A, where b = 2*pi/period, h
            = phase, A = amplitude. For custom spacing, requires
            custom_spacing_second_derivative to be defined.
        """
        amp = self.wing_movement.ampAngles_Gs_to_Wn_ixyz[0]
        b = 2 * np.pi / self.wing_movement.periodAngles_Gs_to_Wn_ixyz[0]
        h = np.deg2rad(self.wing_movement.phaseAngles_Gs_to_Wn_ixyz[0])
        if self.spacing == "sine":
            torque_func = lambda time: -1 * (b**2) * np.sin(b * time + h) * amp * span_I
        elif self.spacing == "uniform":
            raise ValueError(
                "Sawtooth function (uniform spacing) is not differentiable, "
                "cannot be used for inertial torque function."
            )
        elif callable(self.spacing):
            if self.custom_spacing_second_derivative is not None:
                torque_func = (
                    lambda time: self.custom_spacing_second_derivative(time) * span_I
                )
            else:
                raise ValueError(
                    "Custom spacing function provided without second derivative function "
                    "for inertial torque calculation."
                )

        return torque_func

    def spring_numerical_ode(
        self,
        t: np.ndarray,
        k: float,
        c: float,
        I: float,
        theta0: float,
        omega0: float,
        aero_torque: float,
        inertial_torque_func,
    ) -> tuple[float, float]:
        """Numerically integrate the torsional spring-damper ODE.

        Solves the second-order forced ODE: I * d^2(theta)/dt^2 = tau_aero + tau_inertial(t) - k*theta -
        c*d(theta)/dt

        using scipy.integrate.solve_ivp with strict tolerances.

        :param t: A (N,) ndarray of time points for integration evaluation.
        :param k: Spring constant (N*m/rad).
        :param c: Damping constant (N*m*s/rad).
        :param I: Rotational inertia (kg*m^2). This parameter is present for potential
            alternative models of inertia.
        :param theta0: Initial angular displacement (radians).
        :param omega0: Initial angular velocity (rad/s).
        :param aero_torque: Constant aerodynamic torque acting on the section (N*m).
        :param inertial_torque_func: A callable function of time that returns the
            inertial torque from prescribed motion acceleration (N*m).
        :return: A tuple of (theta, omega) representing the final angle and angular
            velocity at the last time point in t.
        """

        def tau(time: float) -> float:
            """Total external torque (aerodynamic + inertial from prescribed motion)."""
            return float(aero_torque + inertial_torque_func(time))

        def ode(time: float, y: np.ndarray) -> np.ndarray:
            """ODE system: d(theta)/dt = omega, d(omega)/dt = (tau - c*omega - k*theta)/I."""
            theta, omega = y
            return np.array([omega, (tau(time) - c * omega - k * theta) / I])

        sol = solve_ivp(
            ode,
            (t[0], t[-1]),
            np.array([theta0, omega0]),
            t_eval=t,
            rtol=1e-9,
            atol=1e-12,
        )

        theta = float(sol.y[0][-1])
        omega = float(sol.y[1][-1])

        return theta, omega

    def plot_flap_cycle_curves(
        self,
        data: list,
        title: str,
        flap_cycle=None,
    ) -> None:
        """Visualize time histories of moments, deformations, or forces.

        Creates a multi-curve line plot showing moment or deformation values across all
        time steps, with optional overlay of a reference flap cycle.

        :param data: A list of lists where each inner list represents a curve to plot.
            Values in each curve are plotted against step number.
        :param title: The title for the plot and the output PNG filename (spaces
            replaced with underscores).
        :param flap_cycle: Optional reference curve to overlay on the plot. If provided,
            should be a list of values to plot with label "Flap Cycle" in black. The
            default is None.
        :return: None **Notes:** The plot is saved as a PNG file with the title as the
            filename. The plot window is displayed to the user. Figure size is 12x6
            inches at 200 DPI.
        """
        plt.figure(figsize=(12, 6), dpi=200)

        for i, curve in enumerate(data):
            x = range(len(curve))
            plt.plot(x, curve, label=f"Curve {i}")
        if flap_cycle is not None:
            plt.plot(
                range(len(flap_cycle)), flap_cycle, label=f"Flap Cycle", color="black"
            )
        plt.xlabel("Step")
        plt.ylabel("Value")
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{title.replace(' ', '_')}.png")
        plt.show()
