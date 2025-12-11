"""Contains the SteadyHorseshoeVortexLatticeMethodSolver class.

**Contains the following classes:**

SteadyHorseshoeVortexLatticeMethodSolver: A class used to solve SteadyProblems with the
horseshoe vortex lattice method.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import numpy as np

from . import _aerodynamics, geometry, operating_point
from . import _functions
from . import _logging
from . import _panel
from . import _parameter_validation
from . import problems

_logger = _logging.get_logger("steady_horseshoe_vortex_lattice_method")


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
class SteadyHorseshoeVortexLatticeMethodSolver:
    """A class used to solve SteadyProblems with the horseshoe vortex lattice method.

    **Contains the following methods:**

    run: Runs the solver on the SteadyProblem.

    calculate_solution_velocity: Finds the fluid velocity (in the first Airplane's
    geometry axes, observed from the Earth frame) at one or more points (in the first
    Airplane's geometry axes, relative to the first Airplane's CG) due to the freestream
    velocity and the induced velocity from every HorseshoeVortex.

    **Citation:**

    Adapted from: aerodynamics.vlm3.py in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/28/2020
    """

    def __init__(self, steady_problem: problems.SteadyProblem) -> None:
        """The initialization method.

        :param steady_problem: The SteadyProblem to be solved.
        :return: None
        """
        if not isinstance(steady_problem, problems.SteadyProblem):
            raise TypeError("steady_problem must be a SteadyProblem.")
        self._steady_problem: problems.SteadyProblem = steady_problem

        self.airplanes = self._steady_problem.airplanes
        self.operating_point: operating_point.OperatingPoint = (
            self._steady_problem.operating_point
        )
        self.reynolds_numbers = self._steady_problem.reynolds_numbers
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
        self.gridStreamlinePoints_GP1_CgP1 = np.empty((0, 3), dtype=float)

        self.ran = False

    def run(self) -> None:
        """Runs the solver on the SteadyProblem.

        :return: None
        """
        # Initialize the Panels' HorseshoeVortices.
        _logger.info("Initializing the Panels' HorseshoeVortices.")
        self._initialize_panel_vortices()

        # Collapse the geometry matrices into 1D ndarrays of attributes.
        _logger.info("Collapsing the geometry.")
        self._collapse_geometry()

        # Find the matrix of Wing-wing influence coefficients associated with this
        # SteadyProblem's geometry.
        _logger.info("Calculating the Wing Wing influences.")
        self._calculate_wing_wing_influences()

        # Find the normal velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at every collocation point due solely to the freestream.
        _logger.info("Calculating the freestream Wing influences.")
        _functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each Panel's HorseshoeVortex's strength.
        _logger.info("Calculating the HorseshoeVortex strengths.")
        self._calculate_vortex_strengths()

        # Solve for the forces (in the first Airplane's geometry axes) and moments (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG)
        # on each Panel.
        _logger.info("Calculating the forces and moments.")
        self._calculate_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges.
        _logger.info("Calculating streamlines.")
        _functions.calculate_streamlines(self)

        # Mark that the solver has run.
        self.ran = True

    def _initialize_panel_vortices(self) -> None:
        """Calculates the locations of the HorseshoeVortex vertices, and then
        initializes the HorseshoeVortices.

        Every Panel has a HorseshoeVortex. The HorseshoeVortices' front legs runs along
        their Panel's quarter chord from right to left. Their quasi infinite legs point
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
                _span = wing.span
                assert _span is not None

                # Find a suitable length for the quasi infinite legs of the
                # HorseshoeVortices on this Wing. At twenty-times the Wing's span,
                # these legs are essentially infinite.
                infinite_leg_length = _span * 20

                _num_spanwise_panels = wing.num_spanwise_panels
                assert _num_spanwise_panels is not None

                # Iterate through the chordwise and spanwise positions of this Wing's
                # Panels.
                for chordwise_position in range(wing.num_chordwise_panels):
                    for spanwise_position in range(_num_spanwise_panels):
                        _panels = wing.panels
                        assert _panels is not None

                        # Pull the Panel out of the Wing's 2D ndarray of Panels.
                        panel: _panel.Panel = _panels[
                            chordwise_position, spanwise_position
                        ]

                        _Frbvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1
                        assert _Frbvp_GP1_CgP1 is not None

                        _Flbvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1
                        assert _Flbvp_GP1_CgP1 is not None

                        # Initialize this Panel's HorseshoeVortex's location (in the
                        # first Airplane's geometry axes, relative to the first
                        # Airplane's CG).
                        panel.horseshoe_vortex = _aerodynamics.HorseshoeVortex(
                            Frhvp_GP1_CgP1=_Frbvp_GP1_CgP1,
                            Flhvp_GP1_CgP1=_Flbvp_GP1_CgP1,
                            leftLegVector_GP1=vInfHat_GP1__E,
                            left_right_leg_lengths=infinite_leg_length,
                            strength=0.0,
                        )

    def _collapse_geometry(self) -> None:
        """Converts attributes of the SteadyProblem's geometry into 1D ndarrays.

        This facilitates vectorization, which speeds up the solver.

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
                _panels = wing.panels
                assert _panels is not None

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(_panels)

                # Iterate through the 1D ndarray of this Wing's Panels.
                panel: _panel.Panel
                for panel in panels:
                    _horseshoe_vortex = panel.horseshoe_vortex
                    assert _horseshoe_vortex is not None

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
                        _horseshoe_vortex.Brhvp_GP1_CgP1
                    )
                    self._stackFrhvp_GP1_CgP1[global_panel_position, :] = (
                        _horseshoe_vortex.Frhvp_GP1_CgP1
                    )
                    self._stackFlhvp_GP1_CgP1[global_panel_position, :] = (
                        _horseshoe_vortex.Flhvp_GP1_CgP1
                    )
                    self._stackBlhvp_GP1_CgP1[global_panel_position, :] = (
                        _horseshoe_vortex.Blhvp_GP1_CgP1
                    )
                    self._stackBoundVortexCenters_GP1_CgP1[global_panel_position, :] = (
                        _horseshoe_vortex.finite_leg.Clvp_GP1_CgP1
                    )
                    self._stackBoundVortexVectors_GP1[global_panel_position, :] = (
                        _horseshoe_vortex.finite_leg.vector_GP1
                    )

                    if panel.is_trailing_edge:
                        _Blpp_GP1_CgP1 = panel.Blpp_GP1_CgP1
                        assert _Blpp_GP1_CgP1 is not None

                        _Brpp_GP1_CgP1 = panel.Brpp_GP1_CgP1
                        assert _Brpp_GP1_CgP1 is not None

                        # Calculate this Panel's streamline seed point (in the first
                        # Airplane's geometry axes, relative to the first Airplane's
                        # CG). Add it to the solver's 1D ndarray of seed points.
                        self.stackSeedPoints_GP1_CgP1 = np.vstack(
                            (
                                self.stackSeedPoints_GP1_CgP1,
                                _Blpp_GP1_CgP1
                                + 0.5 * (_Brpp_GP1_CgP1 - _Blpp_GP1_CgP1),
                            )
                        )

                    # Increment the global Panel position variable.
                    global_panel_position += 1

    def _calculate_wing_wing_influences(self) -> None:
        """Finds this SteadyProblem's 2D ndarray of Wing Wing influence coefficients
        (observed from the Earth frame).

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
        # is now the Problem's 2D ndarray of Wing Wing influence coefficients (observed
        # from the Earth frame).
        self._gridWingWingInfluences__E = np.einsum(
            "...k,...k->...",
            gridNormVIndCpp_GP1__E,
            np.expand_dims(self.stackUnitNormals_GP1, axis=1),
        )

    def _calculate_vortex_strengths(self) -> None:
        """Solves for the strength of each Panel's HorseshoeVortex.

        :return: None
        """
        self._vortex_strengths = np.linalg.solve(
            self._gridWingWingInfluences__E, -self.stackFreestreamWingInfluences__E
        )

        # Update the HorseshoeVortices' strengths.
        panel: _panel.Panel
        for panel_num, panel in enumerate(self.panels):
            horseshoe_vortex = panel.horseshoe_vortex
            assert horseshoe_vortex is not None

            horseshoe_vortex.update_strength(self._vortex_strengths[panel_num])

    def calculate_solution_velocity(
        self, stackP_GP1_CgP1: np.ndarray | Sequence[Sequence[float | int]]
    ) -> np.ndarray:
        """Finds the fluid velocity (in the first Airplane's geometry axes, observed
        from the Earth frame) at one or more points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG) due to the freestream velocity and
        the induced velocity from every HorseshoeVortex.

        **Notes:**

        This method assumes that the correct strengths for the HorseshoeVortices have
        already been calculated and set.

        :param stackP_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (N,3) representing the positions of the evaluation points (in the
            first Airplane's geometry axes, relative to the first Airplane's CG). Can be
            a tuple, list,or ndarray. Values are converted to floats internally. The
            units are in meters.
        :return: A (N,3) ndarray of floats representing the velocity (in the first
            Airplane's geometry axes, observed from the Earth frame) at each evaluation
            point due to the summed effects of the freestream velocity and the induced
            velocity from every HorseshoeVortex. The units are in meters per second.
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

        return cast(np.ndarray, stackVInd_GP1__E + self.vInf_GP1__E)

    def _calculate_loads(self) -> None:
        """Calculates the forces (in the first Airplane's geometry axes) and moments (in
        the first Airplane's geometry axes, relative to the first Airplane's CG) on
        every Panel.

        **Notes:**

        This method assumes that the correct strengths for the HorseshoeVortices have
        already been calculated and set.

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
