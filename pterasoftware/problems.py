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
from coupled_unsteady_ring_vortex_lattice_method import CoupledUnsteadyRingVortexLatticeMethodSolver
from copy import deepcopy
from scipy.integrate import quad


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
    def __init__(self, movement, only_final_results=False):
        super().__init__(movement, only_final_results)
        self.prev_velocities = None
        self.G = 1e8
        self.I_area = 1e-14

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

    def calculate_wing_panel_accelerations(self, solver: CoupledUnsteadyRingVortexLatticeMethodSolver, num_panels):
        dt = self.movement.delta_time
        if len(self._steady_problems) - len():
            self.forces_per_timestep.append(solver.stackCpp_GP1_CgP1)

        wing_panel_veloctiy = solver.calculate_solution_velocity()
        return 

    def initialize_next_problem(self, solver: CoupledUnsteadyRingVortexLatticeMethodSolver):
        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]

        wing = airplane.wings[0]
        mass_matrix = self.define_mass_matrix(0.02, airplane)

        # Panel number definitions
        num_chordwise_panels = wing.num_chordwise_panels
        num_spanwise_panels = wing.num_spanwise_panels
        num_panels = num_chordwise_panels * num_spanwise_panels

        current_torsion_aero = np.zeros(num_chordwise_panels)
        current_torsion_inertia = np.zeros(num_chordwise_panels)
        span_torsion_angles = np.zeros(int(num_spanwise_panels))
        chord_torsion_angles = np.zeros(int(num_spanwise_panels))
        torsion_matrix = np.zeros((num_chordwise_panels, num_spanwise_panels))

        points = np.array(self.stackCpp_G_Cg)[:num_panels, :]
        x_values = points.reshape((num_chordwise_panels, num_spanwise_panels, 3))[
            :num_panels, :, 0
        ]
        panelAeroForces_G = np.stack(
            [o.forces_G for o in np.ravel(wing.panels)]
        ).reshape((num_chordwise_panels, num_spanwise_panels, 3))

        panelInertialForces = self.calculate_wing_panel_accelerations(solver, num_panels) * mass_matrix
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
        # # Error in torsion angle (radians)
        # if self.last_torsion_angles is None:
        #     self.last_torsion_angles = np.zeros(num_spanwise_panels + 1)
        # error_torsion = abs(span_torsion_angles - self.last_torsion_angles)

        # if step < 0.5 * self.num_steps / 3:
        #     # Error threshold for 1st half cycle is high to account for initial spike in forces in UVLM
        #     error_exceeded_air = False
        # else:
        #     # Error threshold in subsequeny timesteps is set to 0.01
        #     # TODO : Determine appropriate error threshold by running test cases on changing wing twist
        #     error_threshold = 0.01 * np.pi / 180  #
        #     # Boolean to determine if change in torsion on any panel exceeds error threshold
        #     error_exceeded_air = np.any(error_torsion[1:] > error_threshold)

        # # If error exceeds the given threshold, the current step is re-run
        # if error_exceeded_air:
        #     return True
        # else:
        #     # Move to next timestep
        #     return False

        # Create new wing based on calculated aero-elastic response
        self.create_new_wing(
            wing=wing,
            airplane=airplane,
            deformation_matrices=span_torsion_angles * 180 / np.pi,
            freestream_velocity=op.vInf_G__E,
            num_chordwise_panels=num_chordwise_panels,
            num_spanwise_panels=num_spanwise_panels,
        )

        # self.last_torsion_angles = span_torsion_angles

        self._steady_problems.append()


    def create_new_wing(
        self,
        wing,
        airplane,
        deformation_matrices,
        freestream_velocity,
        num_chordwise_panels,
        num_spanwise_panels,
    ):
        """This method redefines the current airplane by defining :
        1. new wing cross-section objects, each cross-section's twist = calculated torsion angle
        2. new wing object

        :param wing: Wing object
        :param airplane : Airplane object
        :param torsion_angle : A map from wing cross section to deformations
        :param freestream_velocity : float, current freestream velocity

        :return: None
        """
        # Initialize variable to hold new cross-section objects
        these_cross_sections = []

        # Normalizes the deformation matrix
        if np.max((np.abs(deformation_matrices))) > 1:
            deformation_matrices = deformation_matrices / np.max(
                (np.abs(deformation_matrices))
            )

        # # Create new cross-section objects with calculated torsion angle as wing twist
        for i in range(len(deformation_matrices)):
            this_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                # Every wing cross section has an airfoil object.
                airfoil=geometry.airfoil.Airfoil(
                    name="naca0012",
                    outline_A_lp=None,
                    resample=True,
                    n_points_per_side=400,
                ),
                num_spanwise_panels=wing.wing_cross_sections[i].num_spanwise_panels,
                chord=wing.wing_cross_sections[i].chord,
                Lp_Wcsp_Lpp=wing.wing_cross_sections[i].Lp_Wcsp_Lpp,
                # Every cross-section's twist is due to aero-elastic response
                angles_Wcsp_to_Wcs_ixyz=(
                    np.array([0.0, 0.0, 0.0])
                    if i == 0
                    else np.array([0.0, deformation_matrices[i], 0.0])
                ),
                control_surface_symmetry_type="symmetric",
                control_surface_hinge_point=0.75,
                control_surface_deflection=0.0,
                spanwise_spacing=wing.wing_cross_sections[i].spanwise_spacing,
            )
            these_cross_sections.append(this_wing_cross_section)

        # Define new Wing object with determined cross-sections
        this_wing = geometry.wing.Wing(
            name="Main Wing",
            # Wing root position remains same
            Ler_Gs_Cgs=wing.Ler_Gs_Cgs,
            angles_Gs_to_Wn_ixyz=wing.angles_Gs_to_Wn_ixyz,
            # Wing cross section is redefined
            wing_cross_sections=these_cross_sections,
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=wing.num_chordwise_panels,
            chordwise_spacing=wing.chordwise_spacing,
        )

        # Create new Airplane object from new Wing objects
        these_wings = [this_wing]
        this_airplane = geometry.airplane.Airplane(
            name="Example Airplane",
            wings=these_wings,
            Cg_E_CgP1=airplane.Cg_E_CgP1,
            angles_E_to_B_izyx=airplane.angles_E_to_B_izyx,
            weight=airplane.weight,
        )
        print("sup", airplane.wings, this_airplane.wings)
        # Redefine airplane at current timestep with newly created Airplane object
        self.current_airplanes[0] = this_airplane
        self.steady_problems[self._current_step] = SteadyProblem(
            airplanes=self.current_airplanes, operating_point=self.current_operating_point
        )

    
    def d_alpha_dy_air_static(self, y, tau_torsion, GI):
        return (tau_torsion * y) / GI

    def rotational_inertia(self, m, x, theta):
        return (1 / 3) * m * (x / np.cos(theta)) ** 3
