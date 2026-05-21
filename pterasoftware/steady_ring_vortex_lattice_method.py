"""Contains the SteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

SteadyRingVortexLatticeMethodSolver: A class used to solve SteadyProblems with the ring
vortex lattice method.

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

_logger = _logging.get_logger("steady_ring_vortex_lattice_method")


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
class SteadyRingVortexLatticeMethodSolver:
    """A class used to solve SteadyProblems with the ring vortex lattice method.

    **Contains the following methods:**

    run: Runs the solver on the SteadyProblem.

    calculate_solution_velocity: Finds the fluid velocity (in the first Airplane's
    geometry axes, observed from the Earth frame) at one or more points (in the first
    Airplane's geometry axes, relative to the first Airplane's CG) due to the freestream
    velocity and the induced velocity from every ring vortex and horseshoe vortex.

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
        "num_trailing_edge_panels",
        "vInf_GP1__E",
        "stackFreestreamWingInfluences__E",
        "_gridWingWingInfluences__E",
        "_vortex_strengths",
        "panels",
        "stackUnitNormals_GP1",
        "panel_areas",
        "stackCpp_GP1_CgP1",
        "stackBrbrvp_GP1_CgP1",
        "stackFrbrvp_GP1_CgP1",
        "stackFlbrvp_GP1_CgP1",
        "stackBlbrvp_GP1_CgP1",
        "stackCblvpr_GP1_CgP1",
        "stackCblvpf_GP1_CgP1",
        "stackCblvpl_GP1_CgP1",
        "stackCblvpb_GP1_CgP1",
        "stackRbrv_GP1",
        "stackFbrv_GP1",
        "stackLbrv_GP1",
        "stackBbrv_GP1",
        "_stackBrhvp_GP1_CgP1",
        "_stackFrhvp_GP1_CgP1",
        "_stackFlhvp_GP1_CgP1",
        "_stackBlhvp_GP1_CgP1",
        "_horseshoe_vortex_strengths",
        "_stackRc0s",
        "panel_is_trailing_edge",
        "panel_is_leading_edge",
        "panel_is_left_edge",
        "panel_is_right_edge",
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
        self.num_panels = 0
        self.num_trailing_edge_panels = 0
        airplane: geometry.airplane.Airplane
        for airplane in self.airplanes:
            self.num_panels += airplane.num_panels
            for wing in airplane.wings:
                assert wing.num_spanwise_panels is not None
                self.num_trailing_edge_panels += wing.num_spanwise_panels

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

        # Back-right, front right, front left, and back left bound ring vortex points
        # (in the first Airplane's geometry axes, relative to the first Airplane's CG).
        self.stackBrbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFrbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFlbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackBlbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Center bound line vortex points for the right, front, left, and back legs (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG).
        self.stackCblvpr_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpf_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpl_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpb_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Right, front, left, and back bound ring vortex vectors (in the first
        # Airplane's geometry axes).
        self.stackRbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackLbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
        self.stackBbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)

        # Initialize variables that will hold data which characterizes the trailing
        # edge Panels' horseshoe vortices. These arrays are sized by the number of
        # trailing edge Panels, in trailing edge order, since non trailing edge Panels
        # do not carry a horseshoe vortex.
        self._stackBrhvp_GP1_CgP1 = np.zeros(
            (self.num_trailing_edge_panels, 3), dtype=float
        )
        self._stackFrhvp_GP1_CgP1 = np.zeros(
            (self.num_trailing_edge_panels, 3), dtype=float
        )
        self._stackFlhvp_GP1_CgP1 = np.zeros(
            (self.num_trailing_edge_panels, 3), dtype=float
        )
        self._stackBlhvp_GP1_CgP1 = np.zeros(
            (self.num_trailing_edge_panels, 3), dtype=float
        )
        self._horseshoe_vortex_strengths = np.zeros(
            self.num_trailing_edge_panels, dtype=float
        )

        self._stackRc0s = np.zeros(self.num_panels, dtype=float)

        # Initialize variables to hold details about each Panels' location on its Wing.
        self.panel_is_trailing_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_leading_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_left_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_right_edge = np.zeros(self.num_panels, dtype=bool)

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

        # Compute the bound ring vortex and trailing edge horseshoe vortex geometries
        # and collapse them, along with each Panel's per panel scalars, into 1D
        # ndarrays of attributes.
        _logger.debug("Collapsing geometry.")
        self._collapse_geometry()

        # Find the matrix of Wing Wing influence coefficients associated with this
        # SteadyProblem's geometry.
        _logger.debug("Calculating the Wing Wing influences.")
        self._calculate_wing_wing_influences()

        # Find the normal fluid speed (observed from the Earth frame) at every
        # collocation point due solely to the freestream.
        _logger.debug("Calculating the freestream Wing influences.")
        _functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each Panel's ring vortex's and horseshoe vortex's strength.
        _logger.debug("Calculating ring vortex and horseshoe vortex strengths.")
        self._calculate_vortex_strengths()

        # Solve for the forces (in the first Airplane's geometry axes) and moments (
        # in the first Airplane's geometry axes, relative to the first Airplane's CG)
        # on each Panel.
        _logger.debug("Calculating forces and moments.")
        self._calculate_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges, if requested.
        if calculate_streamlines:
            _logger.debug("Calculating streamlines.")
            _functions.calculate_streamlines(self)

        # Mark that the solver has run.
        self.ran = True

    def _collapse_geometry(self) -> None:
        """Computes the bound ring vortex and trailing edge horseshoe vortex geometries
        and collapses them, along with each Panel's per panel scalars, into the solver's
        1D ndarrays of attributes.

        Every Panel has a bound ring vortex, a quadrangle whose front leg lies on the
        Panel's quarter chord. The left and right legs run along the Panel's left and
        right legs. If the Panel is not along the trailing edge, the back leg meets the
        next Panel's quarter chord. Otherwise, it is offset one quarter of the Panel's
        chord behind the trailing edge.

        Trailing edge Panels also carry a horseshoe vortex whose finite leg overlaps the
        bound ring vortex's back leg in the opposite direction, and whose semi infinite
        legs extend downstream along the freestream. The ring vortex and horseshoe
        vortex carry the same strength, so the back leg contributions cancel.

        :return: None
        """
        # Find the freestream direction (in the first Airplane's geometry axes,
        # observed from the Earth frame).
        vInfHat_GP1__E = self.operating_point.vInfHat_GP1__E

        # Initialize a variable to hold the global position of the Panel as we
        # iterate through them.
        global_panel_position = 0

        # Initialize a variable to hold the trailing edge position as we encounter
        # trailing edge Panels. This indexes the trailing edge sized horseshoe
        # vortex arrays.
        te_position = 0

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

                _num_spanwise_panels = wing.num_spanwise_panels
                assert _num_spanwise_panels is not None

                _panels = wing.panels
                assert _panels is not None

                # Iterate through the chordwise and spanwise positions of this
                # Wing's Panels.
                for chordwise_position in range(wing.num_chordwise_panels):
                    for spanwise_position in range(_num_spanwise_panels):
                        panel: _panel.Panel = _panels[
                            chordwise_position, spanwise_position
                        ]

                        # Find the location of this Panel's front right and front
                        # left bound ring vortex points (in the first Airplane's
                        # geometry axes, relative to the first Airplane's CG).
                        Frrvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1
                        assert Frrvp_GP1_CgP1 is not None

                        Flrvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1
                        assert Flrvp_GP1_CgP1 is not None

                        # Define the location of the back left and back right
                        # bound ring vortex points based on whether the Panel is
                        # along the trailing edge or not.
                        if not panel.is_trailing_edge:
                            next_chordwise_panel = _panels[
                                chordwise_position + 1, spanwise_position
                            ]

                            Blrvp_GP1_CgP1 = next_chordwise_panel.Flbvp_GP1_CgP1
                            assert Blrvp_GP1_CgP1 is not None

                            Brrvp_GP1_CgP1 = next_chordwise_panel.Frbvp_GP1_CgP1
                            assert Brrvp_GP1_CgP1 is not None
                        else:
                            _Blpp_GP1_CgP1 = panel.Blpp_GP1_CgP1
                            assert _Blpp_GP1_CgP1 is not None

                            _Flpp_GP1_CgP1 = panel.Flpp_GP1_CgP1
                            assert _Flpp_GP1_CgP1 is not None

                            _Brpp_GP1_CgP1 = panel.Brpp_GP1_CgP1
                            assert _Brpp_GP1_CgP1 is not None

                            _Frpp_GP1_CgP1 = panel.Frpp_GP1_CgP1
                            assert _Frpp_GP1_CgP1 is not None

                            Blrvp_GP1_CgP1 = Flrvp_GP1_CgP1 + (
                                _Blpp_GP1_CgP1 - _Flpp_GP1_CgP1
                            )
                            Brrvp_GP1_CgP1 = Frrvp_GP1_CgP1 + (
                                _Brpp_GP1_CgP1 - _Frpp_GP1_CgP1
                            )

                        # Update the solver's bound ring vortex stack arrays and
                        # per Panel scalars.
                        _functions.update_ring_vortex_solvers_panel_attributes(
                            ring_vortex_solver=self,
                            global_panel_position=global_panel_position,
                            panel=panel,
                            Frrvp_GP1_CgP1=Frrvp_GP1_CgP1,
                            Flrvp_GP1_CgP1=Flrvp_GP1_CgP1,
                            Blrvp_GP1_CgP1=Blrvp_GP1_CgP1,
                            Brrvp_GP1_CgP1=Brrvp_GP1_CgP1,
                        )

                        if panel.is_trailing_edge:
                            # Populate the horseshoe vortex stack arrays. The
                            # horseshoe vortex shares its front corners with the
                            # bound ring vortex's back corners and its semi infinite
                            # legs extend downstream along the freestream. The
                            # strength is initialized to 1.0 and will be replaced
                            # after the strength solve.
                            self._stackBrhvp_GP1_CgP1[te_position] = (
                                Brrvp_GP1_CgP1 + infinite_leg_offset_GP1
                            )
                            self._stackFrhvp_GP1_CgP1[te_position] = Brrvp_GP1_CgP1
                            self._stackFlhvp_GP1_CgP1[te_position] = Blrvp_GP1_CgP1
                            self._stackBlhvp_GP1_CgP1[te_position] = (
                                Blrvp_GP1_CgP1 + infinite_leg_offset_GP1
                            )
                            self._horseshoe_vortex_strengths[te_position] = 1.0
                            te_position += 1

                        global_panel_position += 1

    def _calculate_wing_wing_influences(self) -> None:
        """Finds this SteadyProblem's 2D ndarray of Wing Wing influence coefficients
        (observed from the Earth frame).

        When an image surface is defined on the OperatingPoint, the influence
        coefficients also include the contributions from image ring vortices and image
        horseshoe vortices reflected across that surface.

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at each Panel's
        # collocation point by each ring vortex. The answer is normalized because the
        # solver's list of ring vortex strengths was initialized to all be 1.0. This
        # will be updated once the correct strengths are calculated.
        singularity_counts = np.zeros(4, dtype=np.int64)
        gridRingNormVIndCpp_GP1__E = (
            _aerodynamics_functions.expanded_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                strengths=self._vortex_strengths,
                r_c0s=self._stackRc0s,
                singularity_counts=singularity_counts,
                ages=None,
                nu=self.operating_point.nu,
            )
        )

        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at every Panel's
        # collocation point by every trailing edge Panel's horseshoe vortex. The
        # answer is normalized because the horseshoe vortex strengths were
        # initialized to 1.0 and will be updated once the correct vortex strengths
        # are calculated. The result has shape (num_panels, num_trailing_edge_panels,
        # 3) and is later scatter added into the ring contribution at the trailing
        # edge columns.
        gridHorseshoeNormVIndCpp_GP1__E = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=self._horseshoe_vortex_strengths,
                r_c0s=self._stackRc0s,
                singularity_counts=singularity_counts,
                nu=self.operating_point.nu,
            )
        )

        # Add the image contributions if an image surface is defined.
        surfaceReflect_T_act_GP1_CgP1 = (
            self.operating_point.surfaceReflect_T_act_GP1_CgP1
        )
        if surfaceReflect_T_act_GP1_CgP1 is not None:
            stackReflectedCpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                self.stackCpp_GP1_CgP1,
                has_point=True,
            )
            gridImageRingVIndCpp_GP1__E = (
                _aerodynamics_functions.expanded_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=stackReflectedCpp_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                    strengths=self._vortex_strengths,
                    r_c0s=self._stackRc0s,
                    singularity_counts=singularity_counts,
                    ages=None,
                    nu=self.operating_point.nu,
                )
            )
            gridRingNormVIndCpp_GP1__E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                gridImageRingVIndCpp_GP1__E,
                has_point=False,
            )
            gridImageHorseshoeVIndCpp_GP1__E = (
                _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                    stackP_GP1_CgP1=stackReflectedCpp_GP1_CgP1,
                    stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                    stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                    stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                    stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                    strengths=self._horseshoe_vortex_strengths,
                    r_c0s=self._stackRc0s,
                    singularity_counts=singularity_counts,
                    nu=self.operating_point.nu,
                )
            )
            gridHorseshoeNormVIndCpp_GP1__E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                gridImageHorseshoeVIndCpp_GP1__E,
                has_point=False,
            )

        unexpected_singularity_counts = np.copy(singularity_counts)
        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.ERROR,
            "_calculate_wing_wing_influences",
            unexpected_singularity_counts,
        )

        # Fold the trailing edge sized horseshoe contribution into the ring grid at
        # the trailing edge columns. For non trailing edge columns the ring grid is
        # untouched.
        gridRingNormVIndCpp_GP1__E[
            :, self.panel_is_trailing_edge, :
        ] += gridHorseshoeNormVIndCpp_GP1__E

        # Take the batch dot product of the normalized induced velocities (in the
        # first Airplane's geometry axes, observed from the Earth frame) with each
        # Panel's unit normal direction (in the first Airplane's geometry axes). This
        # is now the 2D ndarray of Wing Wing influence coefficients (observed from
        # the Earth frame).
        self._gridWingWingInfluences__E = np.einsum(
            "...k,...k->...",
            gridRingNormVIndCpp_GP1__E,
            np.expand_dims(self.stackUnitNormals_GP1, axis=1),
        )

    def _calculate_vortex_strengths(self) -> None:
        """Solves for the bound ring vortex and trailing edge horseshoe vortex
        strengths.

        :return: None
        """
        self._vortex_strengths = np.linalg.solve(
            self._gridWingWingInfluences__E, -self.stackFreestreamWingInfluences__E
        )

        # Mirror the trailing edge entries of the solved strengths into the
        # trailing edge sized horseshoe vortex strength array.
        self._horseshoe_vortex_strengths = self._vortex_strengths[
            self.panel_is_trailing_edge
        ]

    def calculate_solution_velocity(
        self,
        stackP_GP1_CgP1: np.ndarray | Sequence[Sequence[float | int]],
        bound_singularity_counts: np.ndarray | None = None,
    ) -> np.ndarray:
        """Finds the fluid velocity (in the first Airplane's geometry axes, observed
        from the Earth frame) at one or more points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG) due to the freestream velocity and
        the induced velocity from every ring vortex and horseshoe vortex.

        When an image surface is defined on the OperatingPoint, the returned velocity
        also includes the induced velocity from image ring vortices and image horseshoe
        vortices reflected across that surface.

        **Notes:**

        This method assumes that the correct strengths for the ring vortices and
        horseshoe vortices have already been calculated and set.

        :param stackP_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (N,3) representing the positions of the evaluation points (in the
            first Airplane's geometry axes, relative to the first Airplane's CG). Can be
            a tuple, list, or ndarray. Values are converted to floats internally. The
            units are in meters.
        :param bound_singularity_counts: An optional (4,) ndarray of int64 for
            accumulating singularity event counts from bound ring vortices and horseshoe
            vortices. If None, counts are discarded.
        :return: A (N,3) ndarray of floats representing the velocity (in the first
            Airplane's geometry axes, observed from the Earth frame) at each evaluation
            point due to the summed effects of the freestream velocity and the induced
            velocity from every ring vortex and horseshoe vortex. The units are in
            meters per second.
        """
        stackP_GP1_CgP1 = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_GP1_CgP1, "stackP_GP1_CgP1"
            )
        )

        if bound_singularity_counts is None:
            bound_singularity_counts = np.zeros(4, dtype=np.int64)

        stackRingVInd_GP1__E = (
            _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                strengths=self._vortex_strengths,
                r_c0s=self._stackRc0s,
                singularity_counts=bound_singularity_counts,
                ages=None,
                nu=self.operating_point.nu,
            )
        )
        stackHorseshoeVInd_GP1__E = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                strengths=self._horseshoe_vortex_strengths,
                r_c0s=self._stackRc0s,
                singularity_counts=bound_singularity_counts,
                nu=self.operating_point.nu,
            )
        )

        # Add the image contributions if an image surface is defined.
        surfaceReflect_T_act_GP1_CgP1 = (
            self.operating_point.surfaceReflect_T_act_GP1_CgP1
        )
        if surfaceReflect_T_act_GP1_CgP1 is not None:
            stackReflectedP_GP1_CgP1 = _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackP_GP1_CgP1,
                has_point=True,
            )
            stackImageRingVInd_GP1__E = (
                _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=stackReflectedP_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                    strengths=self._vortex_strengths,
                    r_c0s=self._stackRc0s,
                    singularity_counts=bound_singularity_counts,
                    ages=None,
                    nu=self.operating_point.nu,
                )
            )
            stackRingVInd_GP1__E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackImageRingVInd_GP1__E,
                has_point=False,
            )
            stackImageHorseshoeVInd_GP1__E = (
                _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                    stackP_GP1_CgP1=stackReflectedP_GP1_CgP1,
                    stackBrhvp_GP1_CgP1=self._stackBrhvp_GP1_CgP1,
                    stackFrhvp_GP1_CgP1=self._stackFrhvp_GP1_CgP1,
                    stackFlhvp_GP1_CgP1=self._stackFlhvp_GP1_CgP1,
                    stackBlhvp_GP1_CgP1=self._stackBlhvp_GP1_CgP1,
                    strengths=self._horseshoe_vortex_strengths,
                    r_c0s=self._stackRc0s,
                    singularity_counts=bound_singularity_counts,
                    nu=self.operating_point.nu,
                )
            )
            stackHorseshoeVInd_GP1__E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackImageHorseshoeVInd_GP1__E,
                has_point=False,
            )

        return cast(
            np.ndarray,
            stackRingVInd_GP1__E + stackHorseshoeVInd_GP1__E + self.vInf_GP1__E,
        )

    def _calculate_loads(self) -> None:
        """Calculates the forces (in the first Airplane's geometry axes) and moments (in
        the first Airplane's geometry axes, relative to the first Airplane's CG) on
        every Panel.

        **Notes:**

        This method assumes that the correct strengths for the ring vortices and
        horseshoe vortices have already been calculated and set.

        This method used to accidentally double-count the load on each Panel due to the
        left and right line vortex legs. Additionally, it didn't include contributions
        to the load on each Panel from their back line vortex legs. Thankfully, these
        issues only introduced small errors in most typical simulations. They have both
        now been fixed by (1) using a 1/2 factor for each "effective" vortex strength
        shared between two Panels, and (2) including the effects each Panel's back line
        vortex with its own effective strength.

        :return: None
        """
        # Initialize a variable to hold the global Panel position as we iterate
        # through them.
        global_panel_position = 0

        # Initialize three 1D ndarrays to hold the effective strength of the Panels'
        # ring vortices' line vortices.
        effective_right_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_front_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_left_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_back_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)

        # Iterate through the Airplanes' Wings. Within a Wing, Panels are laid out
        # in row major (chordwise outer, spanwise inner) order, so neighbouring
        # Panel strengths can be found at fixed offsets from the current global
        # Panel position: +1 right, -1 left, +num_spanwise back, -num_spanwise
        # front.
        for airplane in self.airplanes:
            for wing in airplane.wings:
                _panels = wing.panels
                assert _panels is not None

                num_spanwise = wing.num_spanwise_panels
                assert num_spanwise is not None

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(_panels)

                # Iterate through this Wing's 1D ndarray of Panels.
                panel: _panel.Panel
                for panel in panels:
                    this_strength = self._vortex_strengths[global_panel_position]

                    if panel.is_right_edge:
                        # Set the effective right line vortex strength to this Panel's
                        # bound ring vortex strength.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            this_strength
                        )
                    else:
                        # Set the effective right line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex strength
                        # and that of the Panel to the right.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._vortex_strengths[global_panel_position + 1]
                        ) / 2

                    if panel.is_leading_edge:
                        # Set the effective front line vortex strength to this Panel's
                        # bound ring vortex strength.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            this_strength
                        )
                    else:
                        # Set the effective front line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex strength
                        # and that of the Panel in front of it.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._vortex_strengths[
                                global_panel_position - num_spanwise
                            ]
                        ) / 2

                    if panel.is_left_edge:
                        # Set the effective left line vortex strength to this Panel's
                        # bound ring vortex strength.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            this_strength
                        )
                    else:
                        # Set the effective left line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex strength
                        # and that of the Panel to the left.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._vortex_strengths[global_panel_position - 1]
                        ) / 2

                    if panel.is_trailing_edge:
                        # Set the effective back line vortex strength to zero, as it
                        # is perfectly canceled by the wake horseshoe vortex's finite
                        # leg line vortex.
                        effective_back_line_vortex_strengths[global_panel_position] = (
                            0.0
                        )
                    else:
                        # Set the effective back line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex strength
                        # and that of the Panel to the back.
                        effective_back_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._vortex_strengths[
                                global_panel_position + num_spanwise
                            ]
                        ) / 2

                    # Increment the global Panel position variable.
                    global_panel_position += 1

        # Calculate the velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at the center of every Panels' ring vortex's right
        # line vortex, front line vortex, left line vortex, and back line vortex.
        bound_singularity_counts = np.zeros(4, dtype=np.int64)
        stackVelocityRightLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpr_GP1_CgP1,
            bound_singularity_counts=bound_singularity_counts,
        )
        stackVelocityFrontLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpf_GP1_CgP1,
            bound_singularity_counts=bound_singularity_counts,
        )
        stackVelocityLeftLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpl_GP1_CgP1,
            bound_singularity_counts=bound_singularity_counts,
        )
        stackVelocityBackLineVortexCenters_GP1__E = self.calculate_solution_velocity(
            stackP_GP1_CgP1=self.stackCblvpb_GP1_CgP1,
            bound_singularity_counts=bound_singularity_counts,
        )

        unexpected_bound_singularity_counts = np.copy(bound_singularity_counts)

        # Subtract expected structural collinearity counts before logging. For
        # each Wing with C chordwise and S spanwise Panels, the four leg center
        # evaluations produce (8 * C * S - 2 * C - 2 * S) collinearity
        # singularities from ring vortex self and adjacent shared edge pairs.
        # Each trailing edge Panel's back leg center is also collinear with the
        # corresponding wake horseshoe vortex's finite leg, adding S per Wing.
        expected_collinearity = 0
        for airplane in self.airplanes:
            for wing in airplane.wings:
                num_chordwise = wing.num_chordwise_panels
                num_spanwise = wing.num_spanwise_panels
                assert num_spanwise is not None
                n = num_chordwise * num_spanwise
                expected_collinearity += 8 * n - 2 * num_chordwise - 2 * num_spanwise
                expected_collinearity += num_spanwise
        unexpected_bound_singularity_counts[3] -= expected_collinearity
        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.ERROR,
            "_calculate_loads (bound)",
            unexpected_bound_singularity_counts,
        )

        # Using the effective line vortex strengths and the Kutta-Joukowski theorem,
        # find the forces (in the first Airplane's geometry axes) on the Panels'
        # ring vortex's right line vortex, front line vortex, left line vortex, and back
        # line vortex using the effective vortex strengths.
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
        backLegForces_GP1 = (
            self.operating_point.rho
            * np.expand_dims(effective_back_line_vortex_strengths, axis=1)
            * np.cross(
                stackVelocityBackLineVortexCenters_GP1__E,
                self.stackBbrv_GP1,
                axis=-1,
            )
        )

        forces_GP1 = (
            rightLegForces_GP1
            + frontLegForces_GP1
            + leftLegForces_GP1
            + backLegForces_GP1
        )

        # TODO: Determine if we get any performance gains by switching to the
        #  functions.numba_1d_explicit_cross function here.
        # Find the moments (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) on the Panels' ring vortex's right line vortex,
        # front line vortex, left line vortex, and back line vortex.
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
        backLegMoments_GP1_CgP1 = np.cross(
            self.stackCblvpb_GP1_CgP1,
            backLegForces_GP1,
            axis=-1,
        )

        moments_GP1_CgP1 = (
            rightLegMoments_GP1_CgP1
            + frontLegMoments_GP1_CgP1
            + leftLegMoments_GP1_CgP1
            + backLegMoments_GP1_CgP1
        )

        # TODO: Transform forces_GP1 and moments_GP1_CgP1 to each Airplane's local
        #  geometry axes before passing to process_solver_loads.
        _functions.process_solver_loads(self, forces_GP1, moments_GP1_CgP1)
