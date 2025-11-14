"""Contains shared utility functions."""

from __future__ import annotations

import logging
from typing import cast

import numpy as np
from numba import njit

from . import _panel
from . import _transformations
from . import steady_horseshoe_vortex_lattice_method
from . import steady_ring_vortex_lattice_method
from . import unsteady_ring_vortex_lattice_method


# TEST: Consider adding unit tests for this function.
def cosspace(
    minimum: float, maximum: float, n_points: int = 50, endpoint: bool = True
) -> np.ndarray:
    """Creates a ndarray containing a specified number of values between a minimum and
    maximum value, spaced via a cosine function.

    **Citation:**

    Adapted from: geometry.cosspace in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/28/2020

    :param minimum: The minimum value of the range of numbers you would like spaced.
    :param maximum: The maximum value of the range of numbers you would like spaced.
    :param n_points: The number of points to space. The default is 50.
    :param endpoint: Determines if the maximum value will be included in the output. The
        default is True.
    :return: A (N,) ndarray of floats ranging from the minimum to the maximum value
        (inclusive), spaced via a cosine function.
    """
    # Find the mean and the amplitude of the cosine function.
    mean = (maximum + minimum) / 2
    amp = (maximum - minimum) / 2

    # Space the points by applying cosine to the output of linspace and return them.
    return cast(
        np.ndarray,
        mean + amp * np.cos(np.linspace(np.pi, 0, n_points, endpoint=endpoint)),
    )


# TEST: Consider adding unit tests for this function.
@njit(cache=True, fastmath=False)
def numba_centroid_of_quadrilateral(
    frontLeftPoint_A_a: np.ndarray,
    frontRightPoint_A_a: np.ndarray,
    backLeftPoint_A_a: np.ndarray,
    backRightPoint_A_a: np.ndarray,
) -> np.ndarray:
    """Finds the centroid of a quadrilateral.

    This function has been optimized for JIT compilation using Numba.

    :param frontLeftPoint_A_a: A (3,) ndarray of floats containing the x, y, and z
        components of the front left point of the quadrilateral (in A axes, relative to
        point a).
    :param frontRightPoint_A_a: A (3,) ndarray of floats containing the x, y, and z
        components of the front right point of the quadrilateral (in A axes, relative to
        point a).
    :param backLeftPoint_A_a: A (3,) ndarray of floats containing the x, y, and z
        components of the back left point of the quadrilateral (in A axes, relative to
        point a).
    :param backRightPoint_A_a: A (3,) ndarray of floats containing the x, y, and z
        components of the back right point of the quadrilateral (in A axes, relative to
        point a).
    :return: A (3,) ndarray of floats containing the x, y, and z components of the
        centroid of the quadrilateral (in A axes, relative to point a).
    """
    x_average = (
        frontLeftPoint_A_a[0]
        + frontRightPoint_A_a[0]
        + backLeftPoint_A_a[0]
        + backRightPoint_A_a[0]
    ) / 4
    y_average = (
        frontLeftPoint_A_a[1]
        + frontRightPoint_A_a[1]
        + backLeftPoint_A_a[1]
        + backRightPoint_A_a[1]
    ) / 4
    z_average = (
        frontLeftPoint_A_a[2]
        + frontRightPoint_A_a[2]
        + backLeftPoint_A_a[2]
        + backRightPoint_A_a[2]
    ) / 4

    return np.array([x_average, y_average, z_average])


# TEST: Consider adding unit tests for this function.
def calculate_streamlines(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    num_steps: int = 25,
    delta_time: float = 0.02,
) -> None:
    """Calculates the location of the streamlines coming off the back of the Wings in a
    solver.

    :param solver: The solver for which to calculate the streamlines.
    :param num_steps: The number of points along each streamline (not including the
        initial points). It can be increased for higher fidelity visuals. The default is
        25.
    :param delta_time: The time in seconds between each time step It can be decreased
        for higher fidelity visuals or to make the streamlines shorter. Its default is
        0.02 seconds.
    :return: None
    """
    # Initialize a ndarray to hold this solver's grid of streamline points (in the first
    # Airplane's geometry axes, relative to the first Airplane's CG).
    assert solver.stackSeedPoints_GP1_CgP1 is not None
    solver.gridStreamlinePoints_GP1_CgP1 = np.expand_dims(
        solver.stackSeedPoints_GP1_CgP1, axis=0
    )

    # Iterate through the streamline time steps.
    for step in range(num_steps):
        # Get the previous row of streamline points (in the first Airplane's geometry
        # axes, relative to the first Airplane's CG).
        lastRowStackStreamlinePoints_GP1_CgP1 = solver.gridStreamlinePoints_GP1_CgP1[
            -1, :, :
        ]

        # Finds the fluid velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at each streamline point in the previous row due to
        # the freestream velocity and the induced velocity from the vortices.
        stackVLastRowStreamlinePoints_GP1__E = solver.calculate_solution_velocity(
            stackP_GP1_CgP1=lastRowStackStreamlinePoints_GP1_CgP1
        )

        # Interpolate to find the new row of streamline points (in the first
        # Airplane's geometry axes, relative to the first Airplane's CG).
        newRowStackStreamlinePoints_GP1_CgP1 = (
            lastRowStackStreamlinePoints_GP1_CgP1
            + stackVLastRowStreamlinePoints_GP1__E * delta_time
        )

        # Stack the new row of streamline points to the bottom of the ndarray of
        # streamline points.
        solver.gridStreamlinePoints_GP1_CgP1 = np.vstack(
            (
                solver.gridStreamlinePoints_GP1_CgP1,
                np.expand_dims(newRowStackStreamlinePoints_GP1_CgP1, axis=0),
            )
        )


# TEST: Consider adding unit tests for this function.
def convert_logging_level_name_to_value(name: str) -> int:
    """Takes in a str, checks that it represents a valid logging level, and converts it
    to the int representation of that level.

    :param name: The string representation of the logging level. The options are
        "Debug", "Info", "Warning", "Error", and "Critical".
    :return: The int that can used to set the appropriate logging level.
    """
    logging_levels = {
        "Debug": logging.DEBUG,
        "Info": logging.INFO,
        "Warning": logging.WARNING,
        "Error": logging.ERROR,
        "Critical": logging.CRITICAL,
    }
    try:
        return logging_levels[name]
    except KeyError:
        raise ValueError(f"{name} is not a valid value of name.")


# TEST: Consider adding unit tests for this function.
def process_solver_loads(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    stackPanelForces_GP1: np.ndarray,
    stackPanelMoments_GP1_CgP1: np.ndarray,
) -> None:
    """Uses the loads a solver has found on its Panels to find and set the net loads and
    associated coefficients for each Airplane.

    :param solver: The solver whose loads will be processed.
    :param stackPanelForces_GP1: A (N,3) ndarray of floats representing the forces (in
        the first Airplane's geometry axes) on each of the solver's Panels. The units
        are in Newtons.
    :param stackPanelMoments_GP1_CgP1: A (N,3) ndarray of floats representing the
        moments (in the first Airplane's geometry axes, relative to the first Airplane's
        CG) on each of the solver's Panels. The units are in Newton-meters.
    :return: None
    """
    if isinstance(
        solver,
        (
            steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
        ),
    ):
        assert solver.airplanes is not None
        these_airplanes = solver.airplanes
        assert solver.operating_point is not None
        this_operating_point = solver.operating_point
    elif isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        assert solver.current_airplanes is not None
        these_airplanes = solver.current_airplanes
        assert solver.current_operating_point is not None
        this_operating_point = solver.current_operating_point
    else:
        raise ValueError(
            f"solver must be a SteadyHorseshoeVortexLatticeMethodSolver, "
            f"a SteadyRingVortexLatticeMethodSolver, or an "
            f"UnsteadyRingVortexLatticeMethodSolver."
        )

    if stackPanelForces_GP1.shape[0] != solver.num_panels:
        raise ValueError(
            f"The first dimension of stackPanelForces_GP1 must equal "
            f"solver.num_panels ({solver.num_panels}), "
            f"got {stackPanelForces_GP1.shape[0]}."
        )

    if stackPanelMoments_GP1_CgP1.shape[0] != solver.num_panels:
        raise ValueError(
            f"The first dimension of stackPanelMoments_GP1_CgP1 must equal "
            f"solver.num_panels ({solver.num_panels}), "
            f"got {stackPanelMoments_GP1_CgP1.shape[0]}."
        )

    qInf__E = this_operating_point.qInf__E
    T_pas_GP1_CgP1_to_W_CgP1 = this_operating_point.T_pas_GP1_CgP1_to_W_CgP1

    # Iterate through this solver's Panels.
    assert solver.panels is not None
    for panel_num, panel in enumerate(solver.panels):
        theseForces_GP1 = stackPanelForces_GP1[panel_num, :]
        theseMoments_GP1_CgP1 = stackPanelMoments_GP1_CgP1[panel_num, :]

        theseForces_W = _transformations.apply_T_to_vectors(
            T_pas_GP1_CgP1_to_W_CgP1, theseForces_GP1, has_point=False
        )
        theseMoments_W_CgP1 = _transformations.apply_T_to_vectors(
            T_pas_GP1_CgP1_to_W_CgP1, theseMoments_GP1_CgP1, has_point=True
        )

        # Update this Panel's loads.
        panel.forces_GP1 = theseForces_GP1
        panel.moments_GP1_CgP1 = theseMoments_GP1_CgP1
        panel.forces_W = theseForces_W
        panel.moments_W_CgP1 = theseMoments_W_CgP1

    # Initialize ndarrays to hold each Airplane's loads.
    stackAirplaneForces_GP1 = np.zeros((solver.num_airplanes, 3), dtype=float)
    stackAirplaneMoments_GP1_CgP1 = np.zeros((solver.num_airplanes, 3), dtype=float)

    # Iterate through each Airplane and find the total loads on each by summing up
    # the contributions from its Panels.
    for airplane_num, airplane in enumerate(these_airplanes):
        for wing in airplane.wings:
            assert wing.panels is not None
            for panel in np.ravel(wing.panels):
                stackAirplaneForces_GP1[airplane_num, :] += panel.forces_GP1
                stackAirplaneMoments_GP1_CgP1[airplane_num, :] += panel.moments_GP1_CgP1

    for airplane_num, airplane in enumerate(these_airplanes):
        airplane.forces_W = _transformations.apply_T_to_vectors(
            T_pas_GP1_CgP1_to_W_CgP1,
            stackAirplaneForces_GP1[airplane_num],
            has_point=False,
        )
        airplane.moments_W_CgP1 = _transformations.apply_T_to_vectors(
            T_pas_GP1_CgP1_to_W_CgP1,
            stackAirplaneMoments_GP1_CgP1[airplane_num],
            has_point=True,
        )

    # Iterate through the Airplanes and calculate each one's load coefficients.
    for airplane in these_airplanes:
        assert airplane.forces_W is not None
        cFX_W = airplane.forces_W[0] / qInf__E / airplane.s_ref
        cFY_W = airplane.forces_W[1] / qInf__E / airplane.s_ref
        cFZ_W = airplane.forces_W[2] / qInf__E / airplane.s_ref

        assert airplane.moments_W_CgP1 is not None
        cMX_W_CgP1 = (
            airplane.moments_W_CgP1[0] / qInf__E / airplane.s_ref / airplane.b_ref
        )
        cMY_W_CgP1 = (
            airplane.moments_W_CgP1[1] / qInf__E / airplane.s_ref / airplane.c_ref
        )
        cMZ_W_CgP1 = (
            airplane.moments_W_CgP1[2] / qInf__E / airplane.s_ref / airplane.b_ref
        )

        # Populate this Airplane's load coefficients.
        airplane.forceCoefficients_W = np.array([cFX_W, cFY_W, cFZ_W])
        airplane.momentCoefficients_W_CgP1 = np.array(
            [cMX_W_CgP1, cMY_W_CgP1, cMZ_W_CgP1]
        )


# TEST: Consider adding unit tests for this function.
def update_ring_vortex_solvers_panel_attributes(
    ring_vortex_solver: (
        steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    global_panel_position: int,
    panel: _panel.Panel,
) -> None:
    """Populates a ring vortex solver's attributes with the attributes of a given Panel.

    :param ring_vortex_solver: The solver whose attributes are to be updated.
    :param global_panel_position: The position of the Panel with respect to the global
        array of Panels.
    :param panel: The Panel whose attributes will be used to update the solver's
        attributes.
    :return: None
    """

    assert ring_vortex_solver.panels is not None
    ring_vortex_solver.panels[global_panel_position] = panel
    assert ring_vortex_solver.stackUnitNormals_GP1 is not None
    ring_vortex_solver.stackUnitNormals_GP1[global_panel_position, :] = (
        panel.unitNormal_GP1
    )
    assert ring_vortex_solver.panel_areas is not None
    ring_vortex_solver.panel_areas[global_panel_position] = panel.area
    assert ring_vortex_solver.stackCpp_GP1_CgP1 is not None
    ring_vortex_solver.stackCpp_GP1_CgP1[global_panel_position, :] = panel.Cpp_GP1_CgP1

    assert panel.ring_vortex is not None
    ring_vortex = panel.ring_vortex

    assert ring_vortex.right_leg is not None
    right_leg = ring_vortex.right_leg
    assert ring_vortex.front_leg is not None
    front_leg = ring_vortex.front_leg
    assert ring_vortex.left_leg is not None
    left_leg = ring_vortex.left_leg
    assert ring_vortex.back_leg is not None
    back_leg = ring_vortex.back_leg

    assert ring_vortex_solver.stackBrbrvp_GP1_CgP1 is not None
    ring_vortex_solver.stackBrbrvp_GP1_CgP1[global_panel_position, :] = (
        right_leg.Slvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackFrbrvp_GP1_CgP1 is not None
    ring_vortex_solver.stackFrbrvp_GP1_CgP1[global_panel_position, :] = (
        right_leg.Elvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackFlbrvp_GP1_CgP1 is not None
    ring_vortex_solver.stackFlbrvp_GP1_CgP1[global_panel_position, :] = (
        left_leg.Slvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackBlbrvp_GP1_CgP1 is not None
    ring_vortex_solver.stackBlbrvp_GP1_CgP1[global_panel_position, :] = (
        left_leg.Elvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackCblvpr_GP1_CgP1 is not None
    ring_vortex_solver.stackCblvpr_GP1_CgP1[global_panel_position, :] = (
        right_leg.Clvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackRbrv_GP1 is not None
    ring_vortex_solver.stackRbrv_GP1[global_panel_position, :] = right_leg.vector_GP1
    assert ring_vortex_solver.stackCblvpf_GP1_CgP1 is not None
    ring_vortex_solver.stackCblvpf_GP1_CgP1[global_panel_position, :] = (
        front_leg.Clvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackFbrv_GP1 is not None
    ring_vortex_solver.stackFbrv_GP1[global_panel_position, :] = front_leg.vector_GP1
    assert ring_vortex_solver.stackCblvpl_GP1_CgP1 is not None
    ring_vortex_solver.stackCblvpl_GP1_CgP1[global_panel_position, :] = (
        left_leg.Clvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackLbrv_GP1 is not None
    ring_vortex_solver.stackLbrv_GP1[global_panel_position, :] = left_leg.vector_GP1
    assert ring_vortex_solver.stackCblvpb_GP1_CgP1 is not None
    ring_vortex_solver.stackCblvpb_GP1_CgP1[global_panel_position, :] = (
        back_leg.Clvp_GP1_CgP1
    )
    assert ring_vortex_solver.stackBbrv_GP1 is not None
    ring_vortex_solver.stackBbrv_GP1[global_panel_position, :] = back_leg.vector_GP1
    assert ring_vortex_solver.panel_is_trailing_edge is not None
    ring_vortex_solver.panel_is_trailing_edge[global_panel_position] = (
        panel.is_trailing_edge
    )
    assert ring_vortex_solver.panel_is_leading_edge is not None
    ring_vortex_solver.panel_is_leading_edge[global_panel_position] = (
        panel.is_leading_edge
    )
    assert ring_vortex_solver.panel_is_right_edge is not None
    ring_vortex_solver.panel_is_right_edge[global_panel_position] = panel.is_right_edge
    assert ring_vortex_solver.panel_is_left_edge is not None
    ring_vortex_solver.panel_is_left_edge[global_panel_position] = panel.is_left_edge

    # Check if this Panel is on the trailing edge. If so, calculate its streamline
    # seed point (in the first Airplane's geometry axes, relative to the first
    # Airplane's CG) and add it to the solver's ndarray of streamline seed points.
    if panel.is_trailing_edge:
        assert ring_vortex_solver.stackSeedPoints_GP1_CgP1 is not None
        assert panel.Brpp_GP1_CgP1 is not None
        assert panel.Blpp_GP1_CgP1 is not None
        ring_vortex_solver.stackSeedPoints_GP1_CgP1 = np.vstack(
            (
                ring_vortex_solver.stackSeedPoints_GP1_CgP1,
                panel.Blpp_GP1_CgP1 + 0.5 * (panel.Brpp_GP1_CgP1 - panel.Blpp_GP1_CgP1),
            )
        )


# TEST: Consider adding unit tests for this function.
def calculate_steady_freestream_wing_influences(
    steady_solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
    ),
) -> None:
    """Finds and sets the vector of freestream-Wing influence coefficients associated
    with a steady solver. These coefficients are the normal velocity (in the first
    Airplane's geometry axes, observed from the Earth frame) at every collocation point
    due solely to the freestream.

    :param steady_solver: The steady solver for which to calculate the freestream-Wing
    influences.
    :return: None
    """
    # Take the batch dot product of the freestream velocity (in the first Airplane's
    # geometry axes, observed from the Earth frame) with each of the N Panel's unit
    # normal vectors (in the first Airplane's geometry axes). This is now the solver's
    # (N,) ndarray of freestream-Wing influence coefficients.
    steady_solver.stackFreestreamWingInfluences__E = np.einsum(
        "ij,j->i",
        steady_solver.stackUnitNormals_GP1,
        steady_solver.vInf_GP1__E,
    )


# TEST: Consider adding unit tests for this function.
@njit(cache=True, fastmath=False)
def numba_1d_explicit_cross(
    stackVectors1_A: np.ndarray, stackVectors2_A: np.ndarray
) -> np.ndarray:
    """Takes in two ndarrays, each of which contain N (3,) vectors and returns the cross
    products of the two vectors at each of the N positions.

    This function has been optimized for JIT compilation using Numba.

    **Citation:**

    Adapted from: https://stackoverflow.com/a/66757029/13240504

    Author: Jérôme Richard

    Date of retrieval: 03/23/2021

    :param stackVectors1_A: The (N,3) ndarray of floats representing the first N vectors
        (in A axes).
    :param stackVectors2_A: The (N,3) ndarray of floats representing the second N
        vectors (in A axes).
    :return: The (N,3) ndarray of floats representing the cross products (in A axes) of
        the two vectors at each of the N positions.
    """
    stackCrossProducts = np.zeros(stackVectors1_A.shape, dtype=float)
    for i in range(stackCrossProducts.shape[0]):
        stackCrossProducts[i, 0] = (
            stackVectors1_A[i, 1] * stackVectors2_A[i, 2]
            - stackVectors1_A[i, 2] * stackVectors2_A[i, 1]
        )
        stackCrossProducts[i, 1] = (
            stackVectors1_A[i, 2] * stackVectors2_A[i, 0]
            - stackVectors1_A[i, 0] * stackVectors2_A[i, 2]
        )
        stackCrossProducts[i, 2] = (
            stackVectors1_A[i, 0] * stackVectors2_A[i, 1]
            - stackVectors1_A[i, 1] * stackVectors2_A[i, 0]
        )
    return stackCrossProducts


# TEST: Consider adding unit tests for this function.
@njit(cache=True, fastmath=False)
def interp_between_points(
    stackStartPoints_A_a: np.ndarray,
    stackEndPoints_A_a: np.ndarray,
    norm_spacings: np.ndarray,
) -> np.ndarray:
    """Returns the (M,N,3) grid of points between M pairs of (3,) points given a ndarray
    of N normalized spacings.

    :param stackStartPoints_A_a: The (M,3) ndarray of floats representing the positions
        of the M start points (in A axes, relative to point a). The units are in meters.
    :param stackEndPoints_A_a: The (M,3) ndarray of floats representing the positions
        of the M end points (in A axes, relative to point a). The units are in meters.
    :param norm_spacings: The (N,) ndarray of floats representing the N normalized
        spacing values between the start and end points. The values are unitless and
        must be normalized to lie in the range from 0.0 to 1.0.
    :return: The (M,N,3) ndarray of floats representing the positions of the MxN
        interpolated points (in A axes, relative to point a). The units are in meters.
    """
    m = stackStartPoints_A_a.shape[0]
    n = norm_spacings.size

    gridInterpolatedPoints_A_a = np.zeros((m, n, 3), dtype=float)

    for i in range(m):
        startPoint_A_a = stackStartPoints_A_a[i, :]
        endPoint_A_a = stackEndPoints_A_a[i, :]

        vector_A = endPoint_A_a - startPoint_A_a

        for j in range(n):
            norm_spacing = norm_spacings[j]

            scaledVector_A = norm_spacing * vector_A

            gridInterpolatedPoints_A_a[i, j, :] = startPoint_A_a + scaledVector_A

    return gridInterpolatedPoints_A_a
