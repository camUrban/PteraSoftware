"""Contains the class definition of this package's steady horseshoe vortex lattice
solver.

This module contains the following classes:
    SteadyHorseshoeVortexLatticeMethodSolver: This is an aerodynamics solver that
    uses a steady horseshoe vortex lattice method.

This module contains the following functions:
    None
"""

import logging

import numpy as np

from . import _aerodynamics, geometry, operating_point
from . import _functions
from . import _panel
from . import _parameter_validation
from . import problems


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
class SteadyHorseshoeVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses a steady horseshoe vortex lattice
    method.

    Citation:
        Adapted from:         aerodynamics.vlm3.py in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    This class contains the following public methods:
        run: Run the solver on the SteadyProblem.

        calculate_solution_velocity: This function takes in a group of points (in the
        first Airplane's geometry axes, relative to the first Airplane's CG). At
        every point, it finds the fluid velocity (in the first Airplane's geometry
        axes, observed from the Earth frame) at that point due to the freestream
        velocity and the induced velocity from every HorseshoeVortex.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, steady_problem: problems.SteadyProblem):
        """This is the initialization method.

        :param steady_problem: SteadyProblem
            This is the SteadyProblem to be solved.
        :return: None
        """
        if not isinstance(steady_problem, problems.SteadyProblem):
            raise TypeError("steady_problem must be a SteadyProblem.")
        self._steady_problem: problems.SteadyProblem = steady_problem

        self.airplanes = self._steady_problem.airplanes
        self.operating_point: operating_point.OperatingPoint = (
            self._steady_problem.operating_point
        )
        self.num_airplanes = len(self.airplanes)

        # Calculate the total number of Panels for all of this SteadyProblem's
        # Airplanes.
        self.num_panels = 0
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            self.num_panels += airplane.num_panels

        # Initialize attributes to hold aerodynamic data that pertains to this
        # SteadyProblem.
        self._gridWingWingInfluences__E = np.zeros(
            (self.num_panels, self.num_panels), dtype=float
        )
        self.vInf_GP1__E = self.operating_point.vInf_GP1__E
        self.stackFreestreamWingInfluences__E = np.zeros(self.num_panels, dtype=float)

        # TODO: The steady ring vortex lattice initializes the strengths to ones,
        #  which makes more sense because then they can be passed in to find the
        #  normalized wing-wing influence coefficients. Switch to that declaration if
        #  doing so doesn't change the results (it shouldn't).
        self._vortex_strengths = np.zeros(self.num_panels, dtype=float)

        self.stackUnitNormals_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._panel_areas = np.zeros(self.num_panels, dtype=float)
        self._stackCpp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        self._stackBrhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._stackFrhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._stackFlhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._stackBlhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        self.panels = np.empty(self.num_panels, dtype=object)
        self._stackBoundVortexCenters_GP1_CgP1 = np.zeros(
            (self.num_panels, 3), dtype=float
        )
        self._stackBoundVortexVectors_GP1 = np.zeros((self.num_panels, 3), dtype=float)

        self.stackSeedPoints_GP1_CgP1 = np.empty((0, 3), dtype=float)
        self.gridStreamlinePoints_GP1_CgP1 = None

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
        logging_level_value = _functions.convert_logging_level_name_to_value(
            logging_level
        )
        logging.basicConfig(level=logging_level_value)

        # Initialize the Panels' HorseshoeVortices.
        logging.info("Initializing the Panels' HorseshoeVortices.")
        self._initialize_panel_vortices()

        # Collapse the geometry matrices into 1D ndarrays of attributes.
        logging.info("Collapsing the geometry.")
        self._collapse_geometry()

        # Find the matrix of Wing-wing influence coefficients associated with this
        # SteadyProblem's geometry.
        logging.info("Calculating the Wing-Wing influences.")
        self._calculate_wing_wing_influences()

        # Find the normal velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at every collocation point due solely to the freestream.
        logging.info("Calculating the freestream-Wing influences.")
        _functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each Panel's HorseshoeVortex's strength.
        logging.info("Calculating the HorseshoeVortex strengths.")
        self._calculate_vortex_strengths()

        # Solve for the forces (in the first Airplane's geometry axes) and moments (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG)
        # on each Panel.
        logging.info("Calculating the forces and moments.")
        self._calculate_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges.
        logging.info("Calculating streamlines.")
        _functions.calculate_streamlines(self)

    def _initialize_panel_vortices(self):
        """This method calculates the locations of the HorseshoeVortex vertices,
        and then initializes the HorseshoeVortices.

        Every Panel has a HorseshoeVortex. The HorseshoeVortices' front legs runs along
        their Panel's quarter chord from right to left. Their quasi-infinite legs point
        backward in the direction of the freestream.

        :return: None
        """
        # Find the freestream direction (in the first Airplane's geometry axes,
        # observed from the Earth frame).
        vInfHat_GP1__E = self.operating_point.vInfHat_GP1__E

        # Iterate through each Airplane's Wings.
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            wing: geometry.wing.Wing
            for wing in airplane.wings:
                # Find a suitable length for the quasi-infinite legs of the
                # HorseshoeVortices on this Wing. At twenty-times the Wing's span,
                # these legs are essentially infinite.
                infinite_leg_length = wing.span * 20

                # Iterate through the chordwise and spanwise positions of this Wing's
                # Panels.
                for chordwise_position in range(wing.num_chordwise_panels):
                    for spanwise_position in range(wing.num_spanwise_panels):
                        # Pull the Panel out of the Wing's 2D ndarray of Panels.
                        panel: _panel.Panel = wing.panels[
                            chordwise_position, spanwise_position
                        ]

                        # Initialize this Panel's HorseshoeVortex's location (in the
                        # first Airplane's geometry axes, relative to the first
                        # Airplane's CG).
                        panel.horseshoe_vortex = _aerodynamics.HorseshoeVortex(
                            Frhvp_GP1_CgP1=panel.Frbvp_GP1_CgP1,
                            Flhvp_GP1_CgP1=panel.Flbvp_GP1_CgP1,
                            leftLegVector_GP1=vInfHat_GP1__E,
                            left_right_leg_lengths=infinite_leg_length,
                            strength=None,
                        )

    def _collapse_geometry(self):
        """This method converts attributes of the SteadyProblem's geometry into 1D
        ndarrays. This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize a variable to hold the global position of the current Panel as we
        # iterate through them.
        global_panel_position = 0

        # Iterate through each Airplane's Wings.
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            wing: geometry.wing.Wing
            for wing in airplane.wings:

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(wing.panels)

                # Iterate through the 1D ndarray of this Wing's Panels.
                panel: _panel.Panel
                for panel in panels:

                    horseshoe_vortex: _aerodynamics.HorseshoeVortex = (
                        panel.horseshoe_vortex
                    )

                    # Update the solver's list of attributes with this Panel's
                    # attributes (in the first Airplane's geometry axes, relative to
                    # the first Airplane's CG).
                    self.panels[global_panel_position] = panel
                    self.stackUnitNormals_GP1[global_panel_position, :] = (
                        panel.unitNormal_GP1
                    )
                    self._panel_areas[global_panel_position] = panel.area
                    self._stackCpp_GP1_CgP1[global_panel_position, :] = (
                        panel.Cpp_GP1_CgP1
                    )
                    self._stackBrhvp_GP1_CgP1[global_panel_position, :] = (
                        horseshoe_vortex.Brhvp_GP1_CgP1
                    )
                    self._stackFrhvp_GP1_CgP1[global_panel_position, :] = (
                        horseshoe_vortex.Frhvp_GP1_CgP1
                    )
                    self._stackFlhvp_GP1_CgP1[global_panel_position, :] = (
                        horseshoe_vortex.Flhvp_GP1_CgP1
                    )
                    self._stackBlhvp_GP1_CgP1[global_panel_position, :] = (
                        horseshoe_vortex.Blhvp_GP1_CgP1
                    )
                    self._stackBoundVortexCenters_GP1_CgP1[global_panel_position, :] = (
                        horseshoe_vortex.finite_leg.Clvp_GP1_CgP1
                    )
                    self._stackBoundVortexVectors_GP1[global_panel_position, :] = (
                        horseshoe_vortex.finite_leg.vector_GP1
                    )

                    if panel.is_trailing_edge:
                        # Calculate this Panel's streamline seed point (in the first
                        # Airplane's geometry axes, relative to the first Airplane's
                        # CG). Add it to the solver's 1D ndarray of seed points.
                        self.stackSeedPoints_GP1_CgP1 = np.vstack(
                            (
                                self.stackSeedPoints_GP1_CgP1,
                                panel.Blpp_GP1_CgP1
                                + 0.5 * (panel.Brpp_GP1_CgP1 - panel.Blpp_GP1_CgP1),
                            )
                        )

                    # Increment the global Panel position variable.
                    global_panel_position += 1

    def _calculate_wing_wing_influences(self):
        """This method finds this SteadyProblem's 2D ndarray of Wing-Wing influence
        coefficients (observed from the Earth frame).

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at each Panel's
        # collocation point by each HorseshoeVortex.
        gridNormVIndCpp_GP1__E = (
            _aerodynamics.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self._stackCpp_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=np.ones(self.num_panels, dtype=float),
                ages=None,
                nu=self.operating_point.nu,
            )
        )

        # Take the batch dot product of the normalized induced velocities (in the
        # first Airplane's geometry axes, observed from the Earth frame) with each
        # Panel's unit normal direction (in the first Airplane's geometry axes). This
        # is now the Problem's 2D ndarray of Wing-Wing influence coefficients (
        # observed from the Earth frame).
        self._gridWingWingInfluences__E = np.einsum(
            "...k,...k->...",
            gridNormVIndCpp_GP1__E,
            np.expand_dims(self.stackUnitNormals_GP1, axis=1),
        )

    def _calculate_vortex_strengths(self):
        """Solve for the strength of each Panel's HorseshoeVortex.

        :return: None
        """
        self._vortex_strengths = np.linalg.solve(
            self._gridWingWingInfluences__E, -self.stackFreestreamWingInfluences__E
        )

        # Update the HorseshoeVortices' strengths.
        this_panel: _panel.Panel
        for panel_num, this_panel in enumerate(self.panels):
            horseshoe_vortex: _aerodynamics.HorseshoeVortex = (
                this_panel.horseshoe_vortex
            )
            horseshoe_vortex.update_strength(self._vortex_strengths[panel_num])

    def calculate_solution_velocity(self, stackP_GP1_CgP1):
        """This function takes in a group of points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). At every point, it finds the
        fluid velocity (in the first Airplane's geometry axes, observed from the
        Earth frame) at that point due to the freestream velocity and the induced
        velocity from every HorseshoeVortex.

        Note: This method assumes that the correct strengths for the
        HorseshoeVortices have already been calculated and set.

        :param stackP_GP1_CgP1: (N,3) array-like of numbers

            Positions of the evaluation points (in the first Airplane's geometry
            axes, relative to the first Airplane's CG). Can be any array-like object
            (tuple, list, or ndarray) with size (N, 3) that has numeric elements (int
            or float). Values are converted to floats internally. The units are in
            meters.

        :return: (N,3) ndarray of floats

            The velocity (in the first Airplane's geometry axes, observed from the
            Earth frame) at every evaluation point due to the summed effects of the
            freestream velocity and the induced velocity from every HorseshoeVortex.
            The units are in meters per second.
        """
        stackP_GP1_CgP1 = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_GP1_CgP1, "stackP_GP1_CgP1"
            )
        )

        stackVInd_GP1__E = _aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
            stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
            stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
            stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
            strengths=self._vortex_strengths,
            ages=None,
            nu=self.operating_point.nu,
        )

        return stackVInd_GP1__E + self.vInf_GP1__E

    def _calculate_loads(self):
        """Calculate the forces (in the first Airplane's geometry axes) and moments (
        in the first Airplane's geometry axes, relative to the first Airplane's CG)
        on every Panel.

        Note: This method assumes that the correct strengths for the
        HorseshoeVortices have already been calculated and set.

        :return: None
        """
        # Calculate the velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at the center of every Panel's HorseshoeVortex's
        # finite leg.
        stackVelocityBoundVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self._stackBoundVortexCenters_GP1_CgP1
        )

        # Calculate the force (in the first Airplane's geometry axes) on each Panel's
        # HorseshoeVortex's finite leg using the Kutta-Joukowski theorem.
        forces_GP1 = (
            self.operating_point.rho
            * np.expand_dims(self._vortex_strengths, axis=1)
            * np.cross(
                stackVelocityBoundVortexCenters_GP1__E,
                self._stackBoundVortexVectors_GP1,
                axis=-1,
            )
        )

        # TODO: Determine if we get any performance gains by switching to the
        #  functions.numba1d_explicit_cross function here.
        # Calculate the moment (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) on each Panel's HorseshoeVortex's finite leg.
        moments_GP1_CgP1 = np.cross(
            self._stackBoundVortexCenters_GP1_CgP1,
            forces_GP1,
            axis=-1,
        )

        # TODO: Transform forces_GP1 and moments_GP1_CgP1 to each Airplane's local
        #  geometry axes before passing to process_solver_loads.
        _functions.process_solver_loads(self, forces_GP1, moments_GP1_CgP1)
