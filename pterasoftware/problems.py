"""This module contains the class definitions for different types of problems.

This module contains the following classes:
    SteadyProblem: This is a class for steady aerodynamics problems.
    UnsteadyProblem: This is a class for unsteady aerodynamics problems.

This module contains the following functions:
    None
"""

import math

import numpy as np


from . import geometry
from . import movements
from . import _parameter_validation
from . import _transformations
from . import operating_point as op
from copy import deepcopy
from scipy.integrate import quad
from .movements.single_step.single_step_movement import SingleStepMovement


class SteadyProblem:
    """This is a class for steady aerodynamics problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, airplanes, operating_point):
        """This is the initialization method.

        :param airplanes: list of Airplanes

            This is the list of the Airplanes for this SteadyProblem.

        :param operating_point: OperatingPoint

            This is the OperatingPoint for this SteadyProblem.
        """
        if not isinstance(airplanes, list):
            raise TypeError("airplanes must be a list.")
        if len(airplanes) < 1:
            raise ValueError("airplanes must have at least one element.")
        for airplane in airplanes:
            if not isinstance(airplane, geometry.airplane.Airplane):
                raise TypeError("Every element in airplanes must be an Airplane.")
        self.airplanes = airplanes
        if not isinstance(operating_point, op.OperatingPoint):
            raise TypeError("operating_point must be an OperatingPoint.")
        self.operating_point = operating_point

        # Validate that the first Airplane has Cg_E_CgP1 set to zeros
        self.airplanes[0].validate_first_airplane_constraints()

        # Populate GP1_CgP1 coordinates for all Airplanes' Panels This finds the
        # Panels' positions in the first Airplanes' geometry axes, relative to the
        # first Airplanes' CG based on their locally defined positions.
        for airplane_id, airplane in enumerate(self.airplanes):
            if airplane_id == 0:
                # First Airplane: use identity transformation (G_Cg == GP1_CgP1)
                T_pas_G_Cg_to_GP1_CgP1 = np.eye(4, dtype=float)
            else:
                # Other Airplanes: compute the passive transformation matrix from
                # this Airplane's local geometry axes, relative to its CG,
                # to the first Airplanes' geometry axes, relative to the first
                # Airplane's CG.
                T_pas_G_Cg_to_GP1_CgP1 = airplane.compute_T_pas_G_Cg_to_GP1_CgP1(
                    first_airplane=self.airplanes[0]
                )

            for wing in airplane.wings:
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
    """This is a class for unsteady aerodynamics problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
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

        # Tunable Parameters
        self.wing_density = 0.012  # per unit height kg/m^2
        self.moment_scaling_factor = 1e-2
        self.spring_constant = 1e-1
        self.aero_scaling = 1.0
        self.new_integrand = False

        # self.wing_density = 0.012  # per unit height kg/m^2
        # self.moment_scaling_factor = 5
        # self.spring_constant = 1
        # self.aero_scaling = 0.0
        # self.new_integrand = True

        self.integration_method = getattr(self, "integration_method", "substep")  # "substep" or "newmark"
        self.max_delta_theta_per_substep = getattr(self, "max_delta_theta_per_substep", 0.005)
        self.max_theta_abs = getattr(self, "max_theta_abs", np.deg2rad(30.0))
        self.max_integration_substeps = getattr(self, "max_integration_substeps", 200)
        self.theta_under_relax = getattr(self, "theta_under_relax", 1.0)
        # Newmark params
        self.newmark_beta = getattr(self, "newmark_beta", 0.25)
        self.newmark_gamma = getattr(self, "newmark_gamma", 0.5)

    def calculate_wing_panel_accelerations(self):
        if len(self.positions) <= 2:
            return np.zeros_like(self.positions[0])
        dt = self.movement.delta_time
        return (self.positions[-1] - 2 * self.positions[-2] + self.positions[-3]) / (dt * dt)

    def initialize_next_problem(self, solver):
        if self.new_integrand:
            deformation_matrices = self.calculate_wing_deformation_new(
                solver, len(self._steady_problems)
            )
        else:
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

    def calculate_wing_deformation_new(self, solver, step):
        """
        Improved torsional-strip update — minimal but correct structural reference axis,
        inertial moment calculation using panel accelerations, and twist measured as
        chord angle. Returns self.net_deformation (same shape as before).
        """

        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]
        wing: geometry.wing.Wing = airplane.wings[0]

        # panel counts
        nc = wing.num_chordwise_panels
        ns = wing.num_spanwise_panels

        # Gather arrays (vectorized)
        aeroMoments = np.array([[p.moments_GP1_CgP1 for p in row] for row in wing.panels])  # (nc,ns,3)
        # append current collocation positions to self.positions (same as you did)
        self.positions.append(np.array([[p.Cpp_GP1_CgP1 for p in row] for row in wing.panels]))  # (nc,ns,3)
        areas = np.array([[p.area for p in row] for row in wing.panels])  # (nc,ns)

        # --- Build structural reference axis points (undeformed mid-chord) ----------
        undeformed_wing = self.steady_problems[-1].airplanes[0].wings[0]
        # we'll compute LE and TE from panel corner properties for each undeformed panel
        span_axis_points = np.zeros((ns, 3))
        span_y = np.zeros(ns)
        for i in range(ns):
            LE_pts = []
            TE_pts = []
            for j in range(nc):
                p = undeformed_wing.panels[j][i]
                LE = np.asarray(p.Flpp_GP1_CgP1)               # front-left point
                TE = np.asarray(p.Brpp_GP1_CgP1) + np.asarray(p.rightLeg_GP1)  # back-right + rightLeg -> TE approx
                LE_pts.append(LE)
                TE_pts.append(TE)
            LE_mean = np.mean(LE_pts, axis=0)
            TE_mean = np.mean(TE_pts, axis=0)
            span_axis_points[i, :] = 0.5 * (LE_mean + TE_mean)  # mid-chord
            span_y[i] = np.mean([pt[1] for pt in LE_pts])

        # approximate dy
        dy = np.mean(np.diff(span_y)) if ns > 1 else 1.0
        dt = float(self.movement.delta_time)

        # --- Panel accelerations & inertial forces --------------------------------
        panel_accels = self.calculate_wing_panel_accelerations()  # (nc,ns,3)
        panel_masses = areas * self.wing_density                     # (nc,ns)
        F_inertial = panel_accels * panel_masses[:, :, np.newaxis]   # (nc,ns,3)

        # --- reference points expanded to panel shape -----------------------------
        ref_points = np.repeat(span_axis_points[np.newaxis, :, :], nc, axis=0)  # (nc,ns,3)

        # --- inertial moment per panel: r x F (vector) ---------------------------
        curr_pos = self.positions[-1]  # (nc,ns,3)
        r = curr_pos - ref_points      # vector from structural axis to panel CG
        panel_inertial_moments = np.cross(r, F_inertial, axis=2)  # (nc,ns,3)

        # --- Choose the correct moment component for torsion ---------------------
        # GP1: +y is span direction. Torsion about span axis => use the y-component (index 1).
        aero_M_y = aeroMoments[:, :, 1]          # (nc,ns)
        inertial_M_y = panel_inertial_moments[:, :, 1]  # (nc,ns)

        # sum chordwise to get per-span driving moment
        M_aero_span = np.sum(aero_M_y, axis=0)      # (ns,)
        M_inertial_span = np.sum(inertial_M_y, axis=0)  # (ns,)
        M_net = M_aero_span - M_inertial_span          # (ns,)

        # --- Rotational inertia about span axis per strip -----------------------
        r_perp_sq = np.sum(r[..., [0, 2]] ** 2, axis=2)  # (nc,ns) using x,z dist
        I_theta = np.sum(panel_masses * r_perp_sq, axis=0)  # (ns,)

        # --- Compute twist angle per span station from chord orientation -------
        # use undeformed LE/TE to get theta_ref if not already set
        if not hasattr(self, "_theta_ref") or self._theta_ref is None:
            self._theta_ref = np.zeros(ns)
            for i in range(ns):
                p_le = np.asarray(undeformed_wing.panels[0][i].Flpp_GP1_CgP1)
                p_te = np.asarray(undeformed_wing.panels[-1][i].Brpp_GP1_CgP1) + np.asarray(undeformed_wing.panels[-1][i].rightLeg_GP1)
                v0 = p_te[[0, 2]] - p_le[[0, 2]]  # projected into x-z
                self._theta_ref[i] = np.arctan2(v0[1], v0[0])
        theta_ref = self._theta_ref

        # initialize dynamic state if missing
        if not hasattr(self, "_theta"):
            self._theta = theta_ref.copy()
            self._theta_dot = np.zeros_like(self._theta)

        # measure current chord-based theta (from current geometry)
        current_theta = np.zeros(ns)
        for i in range(ns):
            # define LE/TE from current panels if corner data available, otherwise use Cpp proxies
            try:
                p_le = np.asarray(wing.panels[0][i].Flpp_GP1_CgP1)
                p_te = np.asarray(wing.panels[-1][i].Brpp_GP1_CgP1) + np.asarray(wing.panels[-1][i].rightLeg_GP1)
            except Exception:
                p_le = curr_pos[0, i]
                p_te = curr_pos[-1, i]
            v = p_te[[0, 2]] - p_le[[0, 2]]
            current_theta[i] = np.arctan2(v[1], v[0])

        # structural parameters (tunable on self)
        GJ = getattr(self, "GJ", 1e4)                 # torsional rigidity
        c_theta = getattr(self, "c_theta", 1e2)       # damping
        k_theta = getattr(self, "k_theta", 0.0)       # optional torsion spring per station (N·m/rad)

        # --- discrete laplacian (spanwise torsion coupling) ---------------------
        lap = np.zeros(ns)
        if ns > 1:
            lap[1:-1] = (self._theta[2:] - 2.0 * self._theta[1:-1] + self._theta[0:-2]) / (dy * dy)
            # boundaries: mirror (free-tip) or use clamped flags
            if getattr(self, "root_clamped", False):
                lap[0] = (self._theta[1] - 2.0 * self._theta[0] + theta_ref[0]) / (dy * dy)
            else:
                lap[0] = (self._theta[1] - 2.0 * self._theta[0] + self._theta[1]) / (dy * dy)
            if getattr(self, "tip_clamped", False):
                lap[-1] = (theta_ref[-1] - 2.0 * self._theta[-1] + self._theta[-2]) / (dy * dy)
            else:
                lap[-1] = (self._theta[-2] - 2.0 * self._theta[-1] + self._theta[-2]) / (dy * dy)

        # --- equation: I*theta_ddot + c*theta_dot + k_theta*(theta-theta_ref) - GJ*lap = M_net
        eps = 1e-12
        I_safe = np.where(I_theta > eps, I_theta, eps)
        torque_spring = k_theta * (self._theta - theta_ref)
        theta_ddot = (M_net - torque_spring + GJ * lap - c_theta * self._theta_dot) / I_safe

        # ---------- SUBSTEPPING + CLIPPING (explicit, conservative) ----------
        if getattr(self, "integration_method", "substep") == "substep":
            # discover / ensure arrays
            dt = float(self.movement.delta_time)
            max_delta = float(self.max_delta_theta_per_substep)
            max_theta_abs = float(self.max_theta_abs)
            under_relax = float(self.theta_under_relax)
            max_substeps_limit = int(self.max_integration_substeps)
            c_theta = float(getattr(self, "c_theta", c_theta))
            k_theta = np.asarray(getattr(self, "k_theta", k_theta))

            # initial theta_ddot (from earlier formula)
            # if you computed theta_ddot earlier as array use it, otherwise compute quickly:
            theta_ddot_curr = (M_net - k_theta * (self._theta - theta_ref) + GJ * lap - c_theta * self._theta_dot) / I_safe

            # heuristic estimate of delta per macro step
            delta_est = np.abs(theta_ddot_curr) * (dt * dt)
            required_substeps = np.ceil(np.maximum(1.0, delta_est / (max_delta + 1e-12))).astype(int)
            substeps = int(min(max(required_substeps.max(), 1), max_substeps_limit))
            dt_sub = dt / float(substeps)

            theta_local = self._theta.copy()
            theta_dot_local = self._theta_dot.copy()

            for s in range(substeps):
                # recompute lap using theta_local
                if ns > 1:
                    lap_local = np.zeros(ns)
                    lap_local[1:-1] = (theta_local[2:] - 2.0 * theta_local[1:-1] + theta_local[0:-2]) / (dy * dy)
                    if getattr(self, "root_clamped", False):
                        lap_local[0] = (theta_local[1] - 2.0 * theta_local[0] + theta_ref[0]) / (dy * dy)
                    else:
                        lap_local[0] = (theta_local[1] - 2.0 * theta_local[0] + theta_local[1]) / (dy * dy)
                    if getattr(self, "tip_clamped", False):
                        lap_local[-1] = (theta_ref[-1] - 2.0 * theta_local[-1] + theta_local[-2]) / (dy * dy)
                    else:
                        lap_local[-1] = (theta_local[-2] - 2.0 * theta_local[-1] + theta_local[-2]) / (dy * dy)
                else:
                    lap_local = np.zeros_like(theta_local)

                torque_spring_local = k_theta * (theta_local - theta_ref)
                theta_ddot_local = (M_net - torque_spring_local + GJ * lap_local - c_theta * theta_dot_local) / I_safe

                # semi-implicit Euler for substep
                theta_dot_new_local = theta_dot_local + dt_sub * theta_ddot_local
                theta_new_local = theta_local + dt_sub * theta_dot_new_local

                # limit per-substep delta
                delta = theta_new_local - theta_local
                too_big = np.abs(delta) > max_delta
                if np.any(too_big):
                    delta[too_big] = np.sign(delta[too_big]) * max_delta

                # under-relaxation
                delta *= under_relax

                theta_local = theta_local + delta
                # update local angular velocity consistent with delta
                theta_dot_local = np.where(dt_sub > 0.0, delta / dt_sub, theta_dot_new_local)

                # absolute clipping
                theta_local = np.clip(theta_local, -max_theta_abs, max_theta_abs)
                clipped = (np.abs(theta_local) >= max_theta_abs - 1e-15)
                if np.any(clipped):
                    theta_dot_local[clipped] = 0.0

            # accept
            theta_new = theta_local
            theta_dot_new = theta_dot_local

            # update object
            self._theta = theta_new
            self._theta_dot = theta_dot_new

        # ---------- NEWMARK-BETA (implicit solver) ----------
        if getattr(self, "integration_method", "substep") == "newmark":
            dt = float(self.movement.delta_time)
            beta = float(self.newmark_beta)
            gamma = float(self.newmark_gamma)
            # diagonal mass and damping
            M_diag = I_safe.copy()                # shape (ns,)
            C_diag = (getattr(self, "c_theta", c_theta) * np.ones(ns)).copy()
            # build discrete second-derivative operator L_matrix (size ns x ns)
            # L_matrix @ theta = (theta_{i+1} - 2 theta_i + theta_{i-1}) / dy^2
            L = np.zeros((ns, ns))
            if ns == 1:
                L[0, 0] = 0.0
            else:
                invdy2 = 1.0 / (dy * dy)
                for i in range(ns):
                    if i == 0:
                        if getattr(self, "root_clamped", False):
                            # clamped: second derivative uses theta_ref later
                            L[i, i] = -2.0 * invdy2
                            L[i, i + 1] = 1.0 * invdy2
                        else:
                            L[i, i] = -2.0 * invdy2
                            L[i, i + 1] = 2.0 * invdy2  # mirrored
                    elif i == ns - 1:
                        if getattr(self, "tip_clamped", False):
                            L[i, i - 1] = 1.0 * invdy2
                            L[i, i] = -2.0 * invdy2
                        else:
                            L[i, i - 1] = 2.0 * invdy2
                            L[i, i] = -2.0 * invdy2
                    else:
                        L[i, i - 1] = 1.0 * invdy2
                        L[i, i] = -2.0 * invdy2
                        L[i, i + 1] = 1.0 * invdy2

            # static stiffness matrix: K_static = diag(k_theta) - GJ * L
            k_theta_arr = np.asarray(getattr(self, "k_theta", k_theta * np.ones(ns)))
            K_static = np.diag(k_theta_arr) - GJ * L  # shape (ns, ns)

            # Effective matrices for Newmark
            # K_eff = M/(beta dt^2) + C*(gamma/(beta dt)) + K_static
            M_over = np.diag(M_diag / (beta * dt * dt))
            C_term = np.diag(C_diag * (gamma / (beta * dt)))
            K_eff = M_over + C_term + K_static

            # right-hand side: F_eff = F_{n+1} + M * a_hat + C * v_hat
            # where F_{n+1} = M_net + k_theta * theta_ref (note k_theta*theta_ref moves to RHS)
            F_ext = M_net + k_theta_arr * theta_ref

            # compute acceleration / velocity predictors from current state
            # assume self._theta_ddot exists or compute current accel:
            if not hasattr(self, "_theta_ddot"):
                # compute current theta_ddot using original formula (elementwise)
                torque_spring = k_theta_arr * (self._theta - theta_ref)
                self._theta_ddot = (M_net - torque_spring + GJ * lap - C_diag * self._theta_dot) / np.where(M_diag > 1e-12, M_diag, 1e-12)

            a_n = self._theta_ddot
            v_n = self._theta_dot
            u_n = self._theta

            a_hat = ( (1.0/(beta*dt*dt)) * u_n +
                    (1.0/(beta*dt)) * v_n +
                    (1.0/(2.0*beta) - 1.0) * a_n )
            v_hat = ( (gamma/(beta*dt)) * u_n +
                    (gamma/beta - 1.0) * v_n +
                    dt * (gamma/(2.0*beta) - 1.0) * a_n )

            F_eff = F_ext + M_diag * a_hat + C_diag * v_hat

            # Handle Dirichlet BCs (clamped nodes) by modifying K_eff and F_eff:
            # for clamped nodes, we enforce theta = theta_ref (Dirichlet)
            clamp_nodes = []
            if getattr(self, "root_clamped", False):
                clamp_nodes.append(0)
            if getattr(self, "tip_clamped", False):
                clamp_nodes.append(ns - 1)
            # implement simple row/col replacement for clamp nodes
            for node in clamp_nodes:
                K_eff[node, :] = 0.0
                K_eff[node, node] = 1.0
                F_eff[node] = theta_ref[node]

            # Solve linear system for theta_{n+1}
            theta_new = np.linalg.solve(K_eff, F_eff)

            # compute acceleration and velocity at n+1
            a_new = (1.0 / (beta * dt * dt)) * (theta_new - u_n) - (1.0 / (beta * dt)) * v_n - (1.0 / (2.0 * beta) - 1.0) * a_n
            v_new = v_n + dt * ( (1.0 - gamma) * a_n + gamma * a_new )

            # update state
            self._theta = theta_new
            self._theta_dot = v_new
            self._theta_ddot = a_new

            # absolute clipping as a final guard
            self._theta = np.clip(self._theta, -self.max_theta_abs, self.max_theta_abs)
            self._theta_dot[np.abs(self._theta) >= self.max_theta_abs - 1e-15] = 0.0

        # --- convert theta_new to net_deformation (backwards-compatible) ----------
        # small-angle vertical displacement approx at half-chord
        half_chords = np.zeros(ns)
        for i in range(ns):
            p_le = np.asarray(undeformed_wing.panels[0][i].Flpp_GP1_CgP1)
            p_te = np.asarray(undeformed_wing.panels[-1][i].Brpp_GP1_CgP1) + np.asarray(undeformed_wing.panels[-1][i].rightLeg_GP1)
            half_chords[i] = 0.5 * abs(p_te[0] - p_le[0])

        delta_theta = theta_new - theta_ref
        step_deformation = np.zeros((ns, 3))
        # put small-angle z (vertical) displacement into index 1 as your original code did
        step_deformation[:, 1] = half_chords * delta_theta * self.moment_scaling_factor

        step_deformation_full = np.insert(step_deformation, 0, np.zeros(3), axis=0)
        if self.net_deformation is None:
            self.net_deformation = np.zeros((ns + 1, 3))
        self.net_deformation += step_deformation_full

        self.net_deformation[:, 1] = np.clip(
            self.net_deformation[:, 1], -90, 90)

        if step % 10 == 3:
            print("Net deformation: ", self.net_deformation)
            print("step deformation: ", step_deformation_full)

        return self.net_deformation

    def calculate_wing_deformation(self, solver, step):
        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]

        wing: geometry.wing.Wing = airplane.wings[0]

        # Panel number definitions
        num_chordwise_panels = wing.num_chordwise_panels
        num_spanwise_panels = wing.num_spanwise_panels
        num_panels = num_chordwise_panels * num_spanwise_panels

        aeroMoments_GP1_CgP1 = np.array(
            [[panel.moments_GP1_CgP1 for panel in row] for row in wing.panels]
        ) * self.aero_scaling
        self.positions.append(np.array(
            [[panel.Cpp_GP1_CgP1 for panel in row] for row in wing.panels]
        ))
        areas = np.array(
            [[panel.area for panel in row] for row in wing.panels]
        )

        inertial_forces = (
            self.calculate_wing_panel_accelerations()
            * np.repeat(areas[:, :, None], 3, axis=2)
            * self.wing_density
        )
        inertial_moments = np.cross(
            self.positions[-1] - self.positions[0],
            inertial_forces,
            axis=2
        )
        total_moments = aeroMoments_GP1_CgP1 - inertial_moments

        deformation_moments = total_moments[:, :, 2]  # Z-axis moments

        if self.net_deformation is None:
            self.net_deformation = np.zeros((num_spanwise_panels + 1, 3))

        undeforemed_wing = self.steady_problems[-1].airplanes[0].wings[0]
        undeformed_postions = np.array(
            [[panel.Cpp_GP1_CgP1 for panel in row] for row in undeforemed_wing.panels]
        )
        step_deformation = np.array(
            [
                np.array(
                    [
                        0,
                        np.sum(deformation_moments[:, i]) * self.moment_scaling_factor
                        - self.spring_constant
                        * np.sum(
                            self.positions[-1][:, i, 2] - undeformed_postions[:, i, 2]
                        ),
                        0,
                    ]
                )
                for i in range(num_spanwise_panels)
            ]
        )
        step_deformation = np.insert(step_deformation, 0, np.array([0,0,0]), axis=0)
        self.net_deformation -= step_deformation

        if step % 10  == 3:
            print("Net deformation: ", self.net_deformation)
            print("step deformation: ", step_deformation)

        return self.net_deformation
