"""Contains the class definition of this package's steady ring vortex lattice solver.

This module contains the following classes:
    SteadyRingVortexLatticeMethodSolver: This is an aerodynamics solver that uses a
    steady ring vortex lattice method.

This module contains the following functions:
    None
"""

from __future__ import annotations

import logging

import numpy as np

from . import _aerodynamics, operating_point, geometry
from . import _functions
from . import _panel
from . import _parameter_validation
from . import problems


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
class SteadyRingVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses a steady ring vortex lattice method.

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
        velocity and the induced velocity from every RingVortex and every
        HorseshoeVortex.

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
        self.num_panels = 0
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            self.num_panels += airplane.num_panels

        # Initialize attributes to hold aerodynamic data that pertains to this
        # simulation.
        self.vInf_GP1__E = self.operating_point.vInf_GP1__E
        self.stackFreestreamWingInfluences__E = np.zeros(self.num_panels, dtype=float)
        self._gridWingWingInfluences__E = np.zeros(
            (self.num_panels, self.num_panels), dtype=float
        )
        self._vortex_strengths = np.ones(self.num_panels, dtype=float)

        self.panels = np.empty(self.num_panels, dtype=object)
        self.stackUnitNormals_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.panel_areas = np.zeros(self.num_panels, dtype=float)

        # Collocation panel points (in the first Airplane's geometry axes, relative
        # to the first Airplane's CG)
        self.stackCpp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Back-right, front-right, front-left, and back-left bound RingVortex points
        # (in the first Airplane's geometry axes, relative to the first Airplane's CG).
        self.stackBrbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFrbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFlbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackBlbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Center bound LineVortex points for the right, front, left, and back legs (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG).
        self.stackCblvpr_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpf_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpl_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpb_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Right, front, left, and back bound RingVortex vectors (in the first
        # Airplane's geometry axes).
        self.stackRbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackLbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackBbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Initialize variables that will hold data which characterizes this Panels'
        # HorseshoeVortex. If the Panel does not have a HorseshoeVortex, these values
        # will not be updated. However, this does not adversely affect the results,
        # because the default HorseshoeVortex strength is zero. The default
        # coordinates will also be updated by the collapse geometry method for Panels
        # that have a HorseshoeVortex.
        self._stackBrhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._stackFrhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._stackFlhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._stackBlhvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self._horseshoe_vortex_strengths = np.zeros(self.num_panels, dtype=float)

        # Initialize variables to hold details about each Panels' location on its Wing.
        self.panel_is_trailing_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_leading_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_left_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_right_edge = np.zeros(self.num_panels, dtype=bool)

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
        logging_level = _parameter_validation.string_return_string(
            logging_level, "logging_level"
        )
        logging_level_value = _functions.convert_logging_level_name_to_value(
            logging_level
        )
        logging.basicConfig(level=logging_level_value)

        # Initialize the Panels' RingVortices and HorseshoeVortices.
        logging.info("Initializing Panels' RingVortices and HorseshoeVortices.")
        self._initialize_panel_vortices()

        # Collapse the geometry matrices into 1D ndarrays of attributes.
        logging.info("Collapsing geometry.")
        self._collapse_geometry()

        # Find the matrix of Wing-Wing influence coefficients associated with this
        # SteadyProblem's geometry.
        logging.info("Calculating the Wing-Wing influences.")
        self._calculate_wing_wing_influences()

        # Find the normal fluid speed (observed from the Earth frame) at every
        # collocation point due solely to the freestream.
        logging.info("Calculating the freestream-Wing influences.")
        _functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each Panel's RingVortex's and HorseshoeVortex's strength.
        logging.info("Calculating RingVortex and HorseshoeVortex strengths.")
        self._calculate_vortex_strengths()

        # Solve for the forces (in the first Airplane's geometry axes) and moments (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG)
        # on each Panel.
        logging.info("Calculating forces and moments.")
        self._calculate_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges.
        logging.info("Calculating streamlines.")
        _functions.calculate_streamlines(self)

    def _initialize_panel_vortices(self):
        """This method calculates the locations of the RingVortex and HorseshoeVortex
        vertices, and then initializes them.

        Every Panel has a RingVortex, which is a quadrangle whose front leg is a
        LineVortex at the Panel's quarter chord. The left and right legs are
        LineVortices running along the Panel's left and right legs. If the Panel is
        not along the trailing edge, they extend backwards and meet the back
        LineVortex, at the rear Panel's quarter chord. Otherwise, they extend
        backwards and meet the back LineVortex one quarter chord back from the
        Panel's back leg.

        Panels that are at the trailing edge of a Wing have a HorseshoeVortex in
        addition to their RingVortex. The HorseshoeVortex's finite leg runs along the
        RingVortex's back leg but in the opposite direction. Its infinite legs point
        backwards in the direction of the freestream. The RingVortex and
        HorseshoeVortex have the same strength, so the effects of the RingVortex's
        back leg's LineVortex and the HorseshoeVortex's finite leg's LineVortex
        cancel each other out.

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
                # HorseshoeVortices on this wing. At twenty-times the Wing's span,
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

                        # Find the location of this Panel's front-left and
                        # front-right RingVortex points (in the first Airplane's
                        # geometry axes, relative to the first Airplane's CG).
                        Flrvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1
                        Frrvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1

                        # Define the location of the back-left and back-right
                        # RingVortex points based on whether the Panel is along the
                        # trailing edge or not.
                        if not panel.is_trailing_edge:
                            next_chordwise_panel = wing.panels[
                                chordwise_position + 1, spanwise_position
                            ]
                            Blrvp_GP1_CgP1 = next_chordwise_panel.Flbvp_GP1_CgP1
                            Brrvp_GP1_CgP1 = next_chordwise_panel.Frbvp_GP1_CgP1
                        else:
                            Blrvp_GP1_CgP1 = Flrvp_GP1_CgP1 + (
                                panel.Blpp_GP1_CgP1 - panel.Flpp_GP1_CgP1
                            )
                            Brrvp_GP1_CgP1 = Frrvp_GP1_CgP1 + (
                                panel.Brpp_GP1_CgP1 - panel.Frpp_GP1_CgP1
                            )

                            # If the Panel is along the trailing edge, initialize its
                            # HorseshoeVortex.
                            panel.horseshoe_vortex = _aerodynamics.HorseshoeVortex(
                                Frhvp_GP1_CgP1=Brrvp_GP1_CgP1,
                                Flhvp_GP1_CgP1=Blrvp_GP1_CgP1,
                                leftLegVector_GP1=vInfHat_GP1__E,
                                left_right_leg_lengths=infinite_leg_length,
                                strength=None,
                            )

                        # Initialize the Panel's RingVortex.
                        panel.ring_vortex = _aerodynamics.RingVortex(
                            Flrvp_GP1_CgP1=Flrvp_GP1_CgP1,
                            Frrvp_GP1_CgP1=Frrvp_GP1_CgP1,
                            Blrvp_GP1_CgP1=Blrvp_GP1_CgP1,
                            Brrvp_GP1_CgP1=Brrvp_GP1_CgP1,
                            strength=None,
                        )

    def _collapse_geometry(self):
        """This method converts attributes of the SteadyProblem's geometry into 1D
        ndarrays. This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize a variable to hold the global position of the Panel as we
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
                    # Update the solver's list of attributes with this Panel's
                    # attributes.
                    _functions.update_ring_vortex_solvers_panel_attributes(
                        ring_vortex_solver=self,
                        global_panel_position=global_panel_position,
                        panel=panel,
                    )

                    if panel.is_trailing_edge:
                        horseshoe_vortex: _aerodynamics.HorseshoeVortex = (
                            panel.horseshoe_vortex
                        )

                        # Update the attribute lists' HorseshoeVortex attributes at
                        # this position with this Panel's HorseshoeVortex attributes
                        self._stackBrhvp_GP1_CgP1[global_panel_position] = (
                            horseshoe_vortex.right_leg.Slvp_GP1_CgP1
                        )
                        self._stackFrhvp_GP1_CgP1[global_panel_position] = (
                            horseshoe_vortex.right_leg.Elvp_GP1_CgP1
                        )
                        self._stackFlhvp_GP1_CgP1[global_panel_position] = (
                            horseshoe_vortex.left_leg.Slvp_GP1_CgP1
                        )
                        self._stackBlhvp_GP1_CgP1[global_panel_position] = (
                            horseshoe_vortex.left_leg.Elvp_GP1_CgP1
                        )

                        # Set the HorseshoeVortex strength at this position to 1.0.
                        # This will be updated after the correct strengths are
                        # calculated.
                        self._horseshoe_vortex_strengths[global_panel_position] = 1.0

                    # Increment the global Panel position variable.
                    global_panel_position += 1

    def _calculate_wing_wing_influences(self):
        """This method finds this SteadyProblem's 2D ndarray of Wing-Wing influence
        coefficients (observed from the Earth frame).

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at each Panel's
        # collocation point by each RingVortex. The answer is normalized because the
        # solver's list of RingVortex strengths was initialized to all be 1.0. This
        # will be updated once the correct strengths are calculated.
        gridRingNormVIndCpp_GP1__E = (
            _aerodynamics.expanded_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                strengths=self._vortex_strengths,
                ages=None,
                nu=self.operating_point.nu,
            )
        )

        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at every Panel's
        # collocation point by every HorseshoeVortex. The answer is normalized
        # because the solver's list of HorseshoeVortex strengths was initialized to
        # 1.0 for locations which have a HorseshoeVortex, and zeros everywhere else.
        # The strengths at the positions with a HorseshoeVortex will be updated once
        # the correct vortex strengths are calculated. The positions elsewhere will
        # remain zero.
        gridHorseshoeNormVIndCpp_GP1__E = (
            _aerodynamics.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=self._horseshoe_vortex_strengths,
                ages=None,
                nu=self.operating_point.nu,
            )
        )

        gridNormVIndCpp_GP1__E = (
            gridRingNormVIndCpp_GP1__E + gridHorseshoeNormVIndCpp_GP1__E
        )

        # Take the batch dot product of the normalized induced velocities (in the
        # first Airplane's geometry axes, observed from the Earth frame) with each
        # Panel's unit normal direction (in the first Airplane's geometry axes). This
        # is now the 2D ndarray of Wing-Wing influence coefficients (observed from
        # the Earth frame).
        self._gridWingWingInfluences__E = np.einsum(
            "...k,...k->...",
            gridNormVIndCpp_GP1__E,
            np.expand_dims(self.stackUnitNormals_GP1, axis=1),
        )

    def _calculate_vortex_strengths(self):
        """Solve for the strength of each Panel's RingVortex and HorseshoeVortex.

        :return: None
        """
        self._vortex_strengths = np.linalg.solve(
            self._gridWingWingInfluences__E, -self.stackFreestreamWingInfluences__E
        )

        # Update the RingVortices' and HorseshoeVortices' strengths.
        panel: _panel.Panel
        for panel_num, panel in enumerate(self.panels):
            ring_vortex: _aerodynamics.RingVortex = panel.ring_vortex
            ring_vortex.update_strength(self._vortex_strengths[panel_num])

            horseshoe_vortex: _aerodynamics.HorseshoeVortex = panel.horseshoe_vortex
            if horseshoe_vortex is not None:
                horseshoe_vortex.update_strength(self._vortex_strengths[panel_num])

                # Also update 1D ndarray of HorseshoeVortex strengths at Panel's
                # location.
                self._horseshoe_vortex_strengths[panel_num] = self._vortex_strengths[
                    panel_num
                ]

    def calculate_solution_velocity(self, stackP_GP1_CgP1):
        """This function takes in a group of points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). At every point, it finds the
        fluid velocity (in the first Airplane's geometry axes, observed from the
        Earth frame) at that point due to the freestream velocity and the induced
        velocity from every RingVortex and every HorseshoeVortex.

        Note: This method assumes that the correct strengths for the RingVortices and
        HorseshoeVortices have already been calculated and set.

        :param stackP_GP1_CgP1: (N,3) array-like of numbers

            Positions of the evaluation points (in the first Airplane's geometry
            axes, relative to the first Airplane's CG). Can be any array-like object
            (tuple, list, or ndarray) with size (N,3) that has numeric elements (int
            or float). Values are converted to floats internally. The units are in
            meters.

        :return: (N,3) ndarray of floats

            The velocity (in the first Airplane's geometry axes, observed from the
            Earth frame) at every evaluation point due to the summed effects of the
            freestream velocity and the induced velocity from every RingVortex and
            HorseshoeVortex. The units are in meters per second.
        """
        stackP_GP1_CgP1 = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_GP1_CgP1, "stackP_GP1_CgP1"
            )
        )

        stackRingVInd_GP1__E = _aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
            stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
            stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
            stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
            strengths=self._vortex_strengths,
            ages=None,
            nu=self.operating_point.nu,
        )
        stackHorseshoeVInd_GP1__E = (
            _aerodynamics.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=self._horseshoe_vortex_strengths,
                ages=None,
                nu=self.operating_point.nu,
            )
        )

        return stackRingVInd_GP1__E + stackHorseshoeVInd_GP1__E + self.vInf_GP1__E

    def _calculate_loads(self):
        """Calculate the forces (in the first Airplane's geometry axes) and moments (
        in the first Airplane's geometry axes, relative to the first Airplane's CG)
        on every Panel.

        Citation: This method uses logic described on pages 9-11 of "Modeling of
        aerodynamic forces in flapping flight with the Unsteady Vortex Lattice
        Method" by Thomas Lambert.

        Note: This method assumes that the correct strengths for the RingVortices and
        HorseshoeVortices have already been calculated and set.

        :return: None
        """
        # Initialize a variable to hold the global Panel position as we iterate
        # through them.
        global_panel_position = 0

        # Initialize three 1D ndarrays to hold the effective strength of the Panels'
        # RingVortices' LineVortices.
        effective_right_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_front_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_left_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)

        # Iterate through the Airplanes' Wings.
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            wing: geometry.wing.Wing
            for wing in airplane.wings:
                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(wing.panels)

                # Iterate through this Wing's 1D ndarray of Panels.
                panel: _panel.Panel
                for panel in panels:

                    # FIXME: After rereading pages 9-10 of "Modeling of aerodynamic
                    #  forces in flapping flight with the Unsteady Vortex Lattice
                    #  Method" by Thomas Lambert, I think our implementation here is
                    #  critically wrong. Consider we have a wing with a (1,2) ndarray
                    #  of Panels. Let's call them Panel A and Panel B. With our
                    #  current method, we calculate the force on Panel A's right
                    #  LineVortex as though it had a strength of Gamma_a - Gamma_b,
                    #  and the force on Panel B's left LineVortex as Gamma_b -
                    #  Gamma_a. I think these forces will precisely cancel-out!

                    if panel.is_right_edge:
                        # Set the effective right LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            self._vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_right: _panel.Panel = wing.panels[
                            panel.local_chordwise_position,
                            panel.local_spanwise_position + 1,
                        ]

                        # Set the effective right LineVortex strength to the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the right.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            self._vortex_strengths[global_panel_position]
                            - panel_to_right.ring_vortex.strength
                        )

                    if panel.is_leading_edge:
                        # Set the effective front LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            self._vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_front: _panel.Panel = wing.panels[
                            panel.local_chordwise_position - 1,
                            panel.local_spanwise_position,
                        ]

                        # Set the effective front LineVortex strength to the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel in front of it.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            self._vortex_strengths[global_panel_position]
                            - panel_to_front.ring_vortex.strength
                        )

                    if panel.is_left_edge:
                        # Set the effective left LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            self._vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_left: _panel.Panel = wing.panels[
                            panel.local_chordwise_position,
                            panel.local_spanwise_position - 1,
                        ]

                        # Set the effective left LineVortex strength to the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the left.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            self._vortex_strengths[global_panel_position]
                            - panel_to_left.ring_vortex.strength
                        )

                    # Increment the global Panel position variable.
                    global_panel_position += 1

        # Calculate the velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at the center of every Panels' RingVortex's right
        # LineVortex, front LineVortex, and left LineVortex.
        stackVelocityRightLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpr_GP1_CgP1
        )
        stackVelocityFrontLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpf_GP1_CgP1
        )
        stackVelocityLeftLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpl_GP1_CgP1
        )

        # Using the effective LineVortex strengths and the Kutta-Joukowski theorem,
        # find the forces (in the first Airplane's geometry axes) on the Panels'
        # RingVortex's right LineVortex, front LineVortex, and left LineVortex using
        # the effective vortex strengths.
        rightLegForces_GP1 = (
            self.operating_point.rho
            * np.expand_dims(effective_right_line_vortex_strengths, axis=1)
            * np.cross(
                stackVelocityRightLineVortexCenters_GP1__E,
                self.stackRbrv_GP1,
                axis=-1,
            )
        )
        frontLegForces_GP1 = (
            self.operating_point.rho
            * np.expand_dims(effective_front_line_vortex_strengths, axis=1)
            * np.cross(
                stackVelocityFrontLineVortexCenters_GP1__E,
                self.stackFbrv_GP1,
                axis=-1,
            )
        )
        leftLegForces_GP1 = (
            self.operating_point.rho
            * np.expand_dims(effective_left_line_vortex_strengths, axis=1)
            * np.cross(
                stackVelocityLeftLineVortexCenters_GP1__E,
                self.stackLbrv_GP1,
                axis=-1,
            )
        )

        forces_GP1 = rightLegForces_GP1 + frontLegForces_GP1 + leftLegForces_GP1

        # TODO: Determine if we get any performance gains by switching to the
        #  functions.numba_1d_explicit_cross function here.
        # Find the moments (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) on the Panels' RingVortex's right LineVortex,
        # front LineVortex, and left LineVortex.
        rightLegMoments_GP1_CgP1 = np.cross(
            self.stackCblvpr_GP1_CgP1,
            rightLegForces_GP1,
            axis=-1,
        )
        frontLegMoments_GP1_CgP1 = np.cross(
            self.stackCblvpf_GP1_CgP1,
            frontLegForces_GP1,
            axis=-1,
        )
        leftLegMoments_GP1_CgP1 = np.cross(
            self.stackCblvpl_GP1_CgP1,
            leftLegForces_GP1,
            axis=-1,
        )

        moments_GP1_CgP1 = (
            rightLegMoments_GP1_CgP1
            + frontLegMoments_GP1_CgP1
            + leftLegMoments_GP1_CgP1
        )

        # TODO: Transform forces_GP1 and moments_GP1_CgP1 to each Airplane's local
        #  geometry axes before passing to process_solver_loads.
        _functions.process_solver_loads(self, forces_GP1, moments_GP1_CgP1)
