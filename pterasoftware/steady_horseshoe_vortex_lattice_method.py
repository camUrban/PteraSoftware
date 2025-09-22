# NOTE: I've started refactoring this module.
"""This module contains the class definition of this package's steady horseshoe
vortex lattice solver.

This module contains the following classes:
    SteadyHorseshoeVortexLatticeMethodSolver: This is an aerodynamics solver that
    uses a steady horseshoe vortex lattice method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import logging

import numpy as np

from . import aerodynamics
from . import functions
from . import geometry
from . import problems


# NOTE: I've started refactoring this class.
class SteadyHorseshoeVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses a steady horseshoe vortex lattice
    method.

    Citation:
        Adapted from:         aerodynamics.vlm3.py in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    This class contains the following public methods:
        run: Run the solver on the SteadyProblem.

        calculate_solution_velocity: This function takes in a group of points. At
        every point, it finds the induced velocity due to every vortex and the
        freestream velocity.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, steady_problem):
        """This is the initialization method.

        :param steady_problem: SteadyProblem
            This is the SteadyProblem to be solved.
        :return: None
        """
        if not isinstance(steady_problem, problems.SteadyProblem):
            raise TypeError("steady_problem must be a SteadyProblem.")
        self.steady_problem = steady_problem

        self.airplanes = self.steady_problem.airplanes
        self.operating_point = self.steady_problem.operating_point
        self.num_airplanes = len(self.airplanes)

        # Calculate the total number of Panels for all of this SteadyProblem's
        # Airplanes.
        self.num_panels = 0
        for airplane in self.airplanes:
            self.num_panels += airplane.num_panels

        # Initialize attributes to hold aerodynamic data that pertains to this
        # SteadyProblem.
        self.stackWingInfluCpp_G__E = np.zeros(
            (self.num_panels, self.num_panels), dtype=float
        )
        self.vInf_G__E = self.operating_point.vInf_G__E
        self.stackVInfInfluCpp_G__E = np.zeros(self.num_panels, dtype=float)

        self.vortex_strengths = np.zeros(self.num_panels, dtype=float)

        self.stackUnitNormal_G = np.zeros((self.num_panels, 3), dtype=float)
        self.panel_areas = np.zeros(self.num_panels, dtype=float)
        self.stackCpp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)

        self.stackBrhvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFrhvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFlhvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
        self.stackBlhvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)

        self.panels = np.empty(self.num_panels, dtype=object)
        self.stackPanelBoundVortexCenter_G_Cg = np.zeros(
            (self.num_panels, 3), dtype=float
        )
        self.stackPanelBoundVortexVector_G = np.zeros((self.num_panels, 3), dtype=float)
        self.stackPanelMomentReference_G_Cg = np.zeros(
            (self.num_panels, 3), dtype=float
        )

        self.stackSeedPoint_G_Cg = np.empty((0, 3), dtype=float)
        self.stackStreamlinePoint_G_Cg = None

    def run(self, logging_level="Warning"):
        """Run the solver on the SteadyProblem.

        :param logging_level: str, optional

            This parameter determines the detail of information that the solver's
            logger will output while running. The options are, in order of detail and
            severity, "Debug", "Info", "Warning", "Error", "Critical". The default
            value is "Warning".

        :return: None
        """
        # Configure the SteadyProblem's logger.
        logging_level_value = functions.convert_logging_level_name_to_value(
            logging_level
        )
        logging.basicConfig(level=logging_level_value)

        # Initialize this Panels to have vortices congruent with this solver type.
        logging.info("Initializing the Panels' HorseshoeVortices.")
        self._initialize_panel_vortices()

        # Collapse the geometry matrices into 1D ndarrays of attributes.
        logging.info("Collapsing the geometry.")
        self._collapse_geometry()

        # Find the matrix of aerodynamic influence coefficients associated with this
        # SteadyProblem's geometry.
        logging.info("Calculating the wing-wing influences.")
        self._calculate_wing_wing_influences()

        # Find the normal velocity (in geometry axes, observed from the Earth frame)
        # at every collocation point due solely to the freestream.
        logging.info("Calculating the freestream-wing influences.")
        functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each Panel's HorseshoeVortex's strength.
        logging.info("Calculating the HorseshoeVortex strengths.")
        self._calculate_vortex_strengths()

        # Solve for the forces (in geometry axes) and moments (in geometry axes,
        # relative to the CG) on each Panel.
        logging.info("Calculating the forces and moments.")
        self._calculate_forces_and_moments()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges.
        logging.info("Calculating streamlines.")
        functions.calculate_streamlines(self)

    def _initialize_panel_vortices(self):
        """This method calculates the locations of the HorseshoeVortex vertices,
        and then initializes the HorseshoeVortices.

        Every Panel has a HorseshoeVortex. The HorseshoeVortices' front legs runs along
        their Panel's quarter chord from right to left. Their quasi-infinite legs point
        backward in the direction of the freestream.

        :return: None
        """
        # Find the freestream direction (in geometry axes, observed from the Earth
        # frame).
        vInfHat_G__E = self.operating_point.vInfHat_G__E

        # Iterate through each Airplane's Wings.
        for airplane in self.airplanes:
            for wing in airplane.wings:

                # Find a suitable length for the quasi-infinite legs of the
                # HorseshoeVortices on this Wing. At twenty-times the Wing's span,
                # these legs are essentially infinite.
                infinite_leg_length = wing.span * 20

                # Iterate through the chordwise and spanwise positions of this Wing's
                # Panels.
                for chordwise_position in range(wing.num_chordwise_panels):
                    for spanwise_position in range(wing.num_spanwise_panels):
                        # Pull the Panel out of the Wing's list of Panels.
                        panel = wing.panels[chordwise_position, spanwise_position]

                        # Initialize this Panel's HorseshoeVortex.
                        panel.horseshoe_vortex = aerodynamics.HorseshoeVortex(
                            Frhvp_G_Cg=panel.Frbvp_G_Cg,
                            Flhvp_G_Cg=panel.Flbvp_G_Cg,
                            leftLegVector_G=vInfHat_G__E,
                            left_right_leg_lengths=infinite_leg_length,
                            strength=None,
                        )

    # NOTE: I've started refactoring this method.
    def _collapse_geometry(self):
        """This method converts attributes of the SteadyProblem's geometry into 1D
        ndarrays. This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize a variable to hold the global position of the current Panel as we
        # iterate through them.
        global_panel_position = 0

        # Iterate through each Airplane's Wings.
        for airplane in self.airplanes:
            for wing in airplane.wings:

                # NOTE: I've refactored up to here.
                # Convert this wing's 2D array of panels into a 1D array.
                panels = np.ravel(wing.panels)

                # Iterate through the 1D array of this wing's panels.
                for panel in panels:

                    # Update the solver's list of attributes with this panel's
                    # attributes.
                    self.panels[global_panel_position] = panel
                    self.stackUnitNormal_G[global_panel_position, :] = (
                        panel.unitNormal_G
                    )
                    self.panel_areas[global_panel_position] = panel.area
                    self.stackCpp_G_Cg[global_panel_position, :] = panel.Cpp_G_Cg
                    self.stackBrhvp_G_Cg[global_panel_position, :] = (
                        panel.horseshoe_vortex.right_leg.Slvp_G_Cg
                    )
                    self.stackFrhvp_G_Cg[global_panel_position, :] = (
                        panel.horseshoe_vortex.right_leg.Elvp_G_Cg
                    )
                    self.stackFlhvp_G_Cg[global_panel_position, :] = (
                        panel.horseshoe_vortex.left_leg.Slvp_G_Cg
                    )
                    self.stackBlhvp_G_Cg[global_panel_position, :] = (
                        panel.horseshoe_vortex.left_leg.Elvp_G_Cg
                    )
                    self.stackPanelBoundVortexCenter_G_Cg[global_panel_position, :] = (
                        panel.horseshoe_vortex.finite_leg.Clvp_G_Cg
                    )
                    self.stackPanelBoundVortexVector_G[global_panel_position, :] = (
                        panel.horseshoe_vortex.finite_leg.vector_G
                    )
                    self.stackPanelMomentReference_G_Cg[global_panel_position, :] = (
                        airplane.Cgi_E_I
                    )

                    # Check if this panel is on the trailing edge.
                    if panel.is_trailing_edge:
                        # If it is, calculate it's streamline seed point and add it
                        # to the solver's array of seed points.
                        self.stackSeedPoint_G_Cg = np.vstack(
                            (
                                self.stackSeedPoint_G_Cg,
                                panel.Blpp_G_Cg
                                + 0.5 * (panel.Brpp_G_Cg - panel.Blpp_G_Cg),
                            )
                        )

                    # Increment the global panel position.
                    global_panel_position += 1

    # NOTE: I haven't yet started refactoring this method.
    def _calculate_wing_wing_influences(self):
        """This method finds the matrix of wing-wing influence coefficients
        associated with the airplanes' geometry.

        :return: None
        """
        # Find the matrix of normalized velocities induced at every panel's
        # collocation point by every panel's horseshoe vortex.
        induced_velocities = aerodynamics.expanded_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.stackCpp_G_Cg,
            stackBrhvp_G_Cg=self.stackBrhvp_G_Cg,
            stackFrhvp_G_Cg=self.stackFrhvp_G_Cg,
            stackFlhvp_G_Cg=self.stackFlhvp_G_Cg,
            stackBlhvp_G_Cg=self.stackBlhvp_G_Cg,
            strengths=np.ones(self.num_panels),
        )

        # Take the batch dot product of the normalized velocities with each panel's
        # normal direction. This is now the problem's matrix of wing-wing influence
        # coefficients.
        self.stackWingInfluCpp_G__E = np.einsum(
            "...k,...k->...",
            induced_velocities,
            np.expand_dims(self.stackUnitNormal_G, axis=1),
        )

    # NOTE: I haven't yet started refactoring this method.
    def _calculate_vortex_strengths(self):
        """Solve for each panel's vortex strengths.

        :return: None
        """
        # Solve for the strength of each panel's vortex.
        self.vortex_strengths = np.linalg.solve(
            self.stackWingInfluCpp_G__E, -self.stackVInfInfluCpp_G__E
        )

        # Iterate through the panels and update their vortex strengths.
        for panel_num in range(len(self.panels)):
            this_panel: geometry.panel.Panel = self.panels[panel_num]

            # Update this panel's horseshoe vortex strength.
            this_panel.horseshoe_vortex.update_strength(
                self.vortex_strengths[panel_num]
            )

    # NOTE: I haven't yet started refactoring this method.
    def calculate_solution_velocity(self, points):
        """This function takes in a group of points. At every point, it finds the
        induced velocity due to every vortex and the freestream velocity.

        Note: The velocity calculated by this method is in geometry axes. Also,
        this method assumes that the correct vortex strengths have already been
        calculated.

        This method uses vectorization, and therefore is much faster for batch
        operations than using the vortex objects' class methods for calculating
        induced velocity.

        :param points: 2D array of floats
            This variable is an array of shape (N x 3), where N is the number of
            points. Each row contains the x, y, and z float coordinates of that
            point's position in meters.
        :return: 2D array of floats
            The output is the summed effects from every vortex, and from the
            freestream on a given point. The result will be of shape (N x 3),
            where each row identifies the velocity at a point. The results units are
            meters per second.
        """
        induced_velocities = aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_G_Cg=points,
            stackBrhvp_G_Cg=self.stackBrhvp_G_Cg,
            stackFrhvp_G_Cg=self.stackFrhvp_G_Cg,
            stackFlhvp_G_Cg=self.stackFlhvp_G_Cg,
            stackBlhvp_G_Cg=self.stackBlhvp_G_Cg,
            strengths=self.vortex_strengths,
        )

        return induced_velocities + self.vInf_G__E

    # NOTE: I haven't yet started refactoring this method.
    def _calculate_forces_and_moments(self):
        """Calculate the forces and moments.

        Note: The forces and moments calculated are in geometry axes. The moment is
        about each airplane's reference point, which should be at the center of
        gravity. The units are Newtons and Newton-meters.

        :return: None
        """
        # Calculate the total velocity at every panel's bound vortex center.
        total_velocities = self.calculate_solution_velocity(
            points=self.stackPanelBoundVortexCenter_G_Cg
        )

        # Calculate the force, in geometry axes, on each panel's bound
        # vortex.
        forces_G = (
            self.operating_point.density
            * np.expand_dims(self.vortex_strengths, axis=1)
            * np.cross(total_velocities, self.stackPanelBoundVortexVector_G, axis=-1)
        )

        # Calculate the moments, in geometry axes, on each panel's bound
        # vortex.
        moments_G_Cg = np.cross(
            self.stackPanelBoundVortexCenter_G_Cg - self.stackPanelMomentReference_G_Cg,
            forces_G,
            axis=-1,
        )

        functions.process_steady_solver_loads(
            steady_solver=self, forces_G=forces_G, moments_G_Cg=moments_G_Cg
        )
