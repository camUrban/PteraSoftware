"""Contains the SteadyHorseshoeVortexLatticeMethodSolver class.

**Contains the following classes:**

SteadyHorseshoeVortexLatticeMethodSolver: A class used to solve SteadyProblems with the
horseshoe vortex lattice method.

**Contains the following functions:**

None
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import cast

import numpy as np

from . import (
    _aerodynamics_functions,
    _functions,
    _logging,
    _panel,
    _parameter_validation,
    _transformations,
    geometry,
    operating_point,
    problems,
)

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
    velocity and the induced velocity from every horseshoe vortex.

    **Citation:**

    Adapted from: aerodynamics.vlm3.py in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/28/2020
    """

    __slots__ = (
        "_steady_problem",
        "airplanes",
        "operating_point",
        "reynolds_numbers",
        "num_airplanes",
        "num_panels",
        "_gridWingWingInfluences__E",
        "vInf_GP1__E",
        "stackFreestreamWingInfluences__E",
        "_vortex_strengths",
        "stackUnitNormals_GP1",
        "_panel_areas",
        "_stackCpp_GP1_CgP1",
        "_stackBrhvp_GP1_CgP1",
        "_stackFrhvp_GP1_CgP1",
        "_stackFlhvp_GP1_CgP1",
        "_stackBlhvp_GP1_CgP1",
        "panels",
        "_stackBoundVortexCenters_GP1_CgP1",
        "_stackBoundVortexVectors_GP1",
        "_stackRc0s",
        "stackSeedPoints_GP1_CgP1",
        "gridStreamlinePoints_GP1_CgP1",
        "ran",
    )

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

        # Initialize the vortex strengths to ones so that they can be passed in to find
        # the normalized wing wing influence coefficients.
        self._vortex_strengths = np.ones(self.num_panels, dtype=float)

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

        self._stackRc0s = np.zeros(self.num_panels, dtype=float)

        self.stackSeedPoints_GP1_CgP1 = np.empty((0, 3), dtype=float)
        self.gridStreamlinePoints_GP1_CgP1 = np.empty((0, 3), dtype=float)

        self.ran = False

    def run(self, calculate_streamlines: bool | np.bool_ = True) -> None:
        """Runs the solver on the SteadyProblem.

        :param calculate_streamlines: Determines whether to calculate the streamlines
            emanating from the back of the wing after running the solver. Can be a bool
            or a numpy bool and will be converted internally to a bool. The default is
            True.
        :return: None
        """
        calculate_streamlines = _parameter_validation.boolLike_return_bool(
            calculate_streamlines, "calculate_streamlines"
        )

        # Compute the horseshoe vortex geometries and collapse them, along with each
        # Panel's per panel scalars, into 1D ndarrays of attributes.
        _logger.debug("Collapsing the geometry.")
        self._collapse_geometry()

        # Find the matrix of Wing-wing influence coefficients associated with this
        # SteadyProblem's geometry.
        _logger.debug("Calculating the Wing Wing influences.")
        self._calculate_wing_wing_influences()

        # Find the normal velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at every collocation point due solely to the freestream.
        _logger.debug("Calculating the freestream Wing influences.")
        _functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each Panel's horseshoe vortex's strength.
        _logger.debug("Calculating the horseshoe vortex strengths.")
        self._calculate_vortex_strengths()

        # Solve for the forces (in the first Airplane's geometry axes) and moments (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG)
        # on each Panel.
        _logger.debug("Calculating the forces and moments.")
        self._calculate_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges, if requested.
        if calculate_streamlines:
            _logger.debug("Calculating streamlines.")
            _functions.calculate_streamlines(self)

        # Mark that the solver has run.
        self.ran = True

    def _collapse_geometry(self) -> None:
        """Computes the horseshoe vortex geometries and collapses them, along with each
        Panel's per panel scalars, into the solver's 1D ndarrays of attributes.

        Every Panel carries a horseshoe vortex. The finite leg runs along the Panel's
        quarter chord from right to left. The semi infinite legs extend downstream from
        the front corners along the freestream direction.

        :return: None
        """
        # Find the freestream direction (in the first Airplane's geometry axes,
        # observed from the Earth frame).
        vInfHat_GP1__E = self.operating_point.vInfHat_GP1__E

        # Initialize a variable to hold the global position of the current Panel as
        # we iterate through them.
        global_panel_position = 0

        # Iterate through each Airplane's Wings.
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            wing: geometry.wing.Wing
            for wing in airplane.wings:
                _span = wing.span
                assert _span is not None
                # At twenty times the Wing's span, the horseshoe vortex legs are
                # essentially infinite.
                infinite_leg_offset_GP1 = vInfHat_GP1__E * (_span * 20)

                _panels = wing.panels
                assert _panels is not None

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(_panels)

                # Iterate through the 1D ndarray of this Wing's Panels.
                panel: _panel.Panel
                for panel in panels:
                    Frhvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1
                    assert Frhvp_GP1_CgP1 is not None

                    Flhvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1
                    assert Flhvp_GP1_CgP1 is not None

                    # The semi infinite legs trail downstream from the front
                    # corners along the freestream direction.
                    Brhvp_GP1_CgP1 = Frhvp_GP1_CgP1 + infinite_leg_offset_GP1
                    Blhvp_GP1_CgP1 = Flhvp_GP1_CgP1 + infinite_leg_offset_GP1

                    # Update the solver's list of attributes with this Panel's
                    # attributes (in the first Airplane's geometry axes, relative
                    # to the first Airplane's CG).
                    self.panels[global_panel_position] = panel
                    self.stackUnitNormals_GP1[global_panel_position, :] = (
                        panel.unitNormal_GP1
                    )
                    self._panel_areas[global_panel_position] = panel.area
                    self._stackCpp_GP1_CgP1[global_panel_position, :] = (
                        panel.Cpp_GP1_CgP1
                    )

                    self._stackBrhvp_GP1_CgP1[global_panel_position, :] = Brhvp_GP1_CgP1
                    self._stackFrhvp_GP1_CgP1[global_panel_position, :] = Frhvp_GP1_CgP1
                    self._stackFlhvp_GP1_CgP1[global_panel_position, :] = Flhvp_GP1_CgP1
                    self._stackBlhvp_GP1_CgP1[global_panel_position, :] = Blhvp_GP1_CgP1

                    # The finite leg runs from the front right to the front left.
                    self._stackBoundVortexCenters_GP1_CgP1[global_panel_position, :] = (
                        0.5 * (Frhvp_GP1_CgP1 + Flhvp_GP1_CgP1)
                    )
                    self._stackBoundVortexVectors_GP1[global_panel_position, :] = (
                        Flhvp_GP1_CgP1 - Frhvp_GP1_CgP1
                    )

                    if panel.is_trailing_edge:
                        _Blpp_GP1_CgP1 = panel.Blpp_GP1_CgP1
                        assert _Blpp_GP1_CgP1 is not None

                        _Brpp_GP1_CgP1 = panel.Brpp_GP1_CgP1
                        assert _Brpp_GP1_CgP1 is not None

                        # Calculate this Panel's streamline seed point (in the
                        # first Airplane's geometry axes, relative to the first
                        # Airplane's CG). Add it to the solver's 1D ndarray of
                        # seed points.
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

        When an image surface is defined on the OperatingPoint, the influence
        coefficients also include the contributions from image horseshoe vortices
        reflected across that surface.

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at each Panel's
        # collocation point by each horseshoe vortex.
        singularity_counts = np.zeros(4, dtype=np.int64)
        gridNormVIndCpp_GP1__E = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self._stackCpp_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=self._vortex_strengths,
                r_c0s=self._stackRc0s,
                singularity_counts=singularity_counts,
                nu=self.operating_point.nu,
            )
        )

        # Add the image contribution if an image surface is defined.
        surfaceReflect_T_act_GP1_CgP1 = (
            self.operating_point.surfaceReflect_T_act_GP1_CgP1
        )
        if surfaceReflect_T_act_GP1_CgP1 is not None:
            stackReflectedCpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                self._stackCpp_GP1_CgP1,
                has_point=True,
            )
            gridImageVIndCpp_GP1__E = (
                _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                    stackP_GP1_CgP1=stackReflectedCpp_GP1_CgP1,
                    stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                    stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                    stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                    stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                    strengths=self._vortex_strengths,
                    r_c0s=self._stackRc0s,
                    singularity_counts=singularity_counts,
                    nu=self.operating_point.nu,
                )
            )
            gridNormVIndCpp_GP1__E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                gridImageVIndCpp_GP1__E,
                has_point=False,
            )

        unexpected_singularity_counts = np.copy(singularity_counts)

        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.ERROR,
            "_calculate_wing_wing_influences",
            unexpected_singularity_counts,
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
        """Solves for the strength of each Panel's horseshoe vortex.

        :return: None
        """
        self._vortex_strengths = np.linalg.solve(
            self._gridWingWingInfluences__E, -self.stackFreestreamWingInfluences__E
        )

    def calculate_solution_velocity(
        self,
        stackP_GP1_CgP1: np.ndarray | Sequence[Sequence[float | int]],
        bound_singularity_counts: np.ndarray | None = None,
    ) -> np.ndarray:
        """Finds the fluid velocity (in the first Airplane's geometry axes, observed
        from the Earth frame) at one or more points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG) due to the freestream velocity and
        the induced velocity from every horseshoe vortex.

        When an image surface is defined on the OperatingPoint, the returned velocity
        also includes the induced velocity from image horseshoe vortices reflected
        across that surface.

        **Notes:**

        This method assumes that the correct strengths for the horseshoe vortices have
        already been calculated and set.

        :param stackP_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (N,3) representing the positions of the evaluation points (in the
            first Airplane's geometry axes, relative to the first Airplane's CG). Can be
            a tuple, list,or ndarray. Values are converted to floats internally. The
            units are in meters.
        :param bound_singularity_counts: An optional (4,) ndarray of int64 for
            accumulating singularity event counts from bound horseshoe vortices. If
            None, counts are discarded.
        :return: A (N,3) ndarray of floats representing the velocity (in the first
            Airplane's geometry axes, observed from the Earth frame) at each evaluation
            point due to the summed effects of the freestream velocity and the induced
            velocity from every horseshoe vortex. The units are in meters per second.
        """
        stackP_GP1_CgP1 = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_GP1_CgP1, "stackP_GP1_CgP1"
            )
        )

        if bound_singularity_counts is None:
            bound_singularity_counts = np.zeros(4, dtype=np.int64)

        stackVInd_GP1__E = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=self._vortex_strengths,
                r_c0s=self._stackRc0s,
                singularity_counts=bound_singularity_counts,
                nu=self.operating_point.nu,
            )
        )

        # Add the image contribution if an image surface is defined.
        surfaceReflect_T_act_GP1_CgP1 = (
            self.operating_point.surfaceReflect_T_act_GP1_CgP1
        )
        if surfaceReflect_T_act_GP1_CgP1 is not None:
            stackReflectedP_GP1_CgP1 = _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackP_GP1_CgP1,
                has_point=True,
            )
            stackImageVInd_GP1__E = (
                _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                    stackP_GP1_CgP1=stackReflectedP_GP1_CgP1,
                    stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                    stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                    stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                    stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                    strengths=self._vortex_strengths,
                    r_c0s=self._stackRc0s,
                    singularity_counts=bound_singularity_counts,
                    nu=self.operating_point.nu,
                )
            )
            stackVInd_GP1__E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackImageVInd_GP1__E,
                has_point=False,
            )

        return cast(np.ndarray, stackVInd_GP1__E + self.vInf_GP1__E)

    def _calculate_loads(self) -> None:
        """Calculates the forces (in the first Airplane's geometry axes) and moments (in
        the first Airplane's geometry axes, relative to the first Airplane's CG) on
        every Panel.

        **Notes:**

        This method assumes that the correct strengths for the horseshoe vortices have
        already been calculated and set.

        :return: None
        """
        # Calculate the velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at the center of every Panel's horseshoe vortex's
        # finite leg.
        bound_singularity_counts = np.zeros(4, dtype=np.int64)
        stackVelocityBoundVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self._stackBoundVortexCenters_GP1_CgP1,
            bound_singularity_counts=bound_singularity_counts,
        )

        unexpected_bound_singularity_counts = np.copy(bound_singularity_counts)

        # Subtract the expected structural collinearity before logging. Each
        # bound vortex center is collinear with its own finite leg, producing
        # exactly one collinearity singularity per Panel.
        unexpected_bound_singularity_counts[3] -= self.num_panels
        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.ERROR,
            "_calculate_loads (bound)",
            unexpected_bound_singularity_counts,
        )

        # Calculate the force (in the first Airplane's geometry axes) on each Panel's
        # horseshoe vortex's finite leg using the Kutta-Joukowski theorem.
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
        # first Airplane's CG) on each Panel's horseshoe vortex's finite leg.
        moments_GP1_CgP1 = np.cross(
            self._stackBoundVortexCenters_GP1_CgP1,
            forces_GP1,
            axis=-1,
        )

        # TODO: Transform forces_GP1 and moments_GP1_CgP1 to each Airplane's local
        #  geometry axes before passing to process_solver_loads.
        _functions.process_solver_loads(self, forces_GP1, moments_GP1_CgP1)
