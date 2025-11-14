"""This module contains functions shared by other modules in the pterasoftware
package."""

from __future__ import annotations

import logging

import numpy as np
from numba import njit

from . import _panel
from . import _parameter_validation
from . import _transformations
from . import geometry
from . import steady_horseshoe_vortex_lattice_method
from . import steady_ring_vortex_lattice_method
from . import unsteady_ring_vortex_lattice_method


# TEST: Consider adding unit tests for this function.
def cosspace(minimum, maximum, n_points=50, endpoint=True):
    """This function is used to create an array containing a specified number of
    values between a specified minimum and maximum value that are spaced via a cosine
    function.

    Citation:
        Adapted from:         geometry.cosspace in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    :param minimum: float
        This is the minimum value of the range of numbers you would like spaced.
    :param maximum: float
        This is the maximum value of the range of numbers you would like spaced.
    :param n_points: int, optional
        This is the number of points to space. The default is 50.
    :param endpoint: bool, optional
        This determines if the maximum value will be included in the output. The
        default is True.
    :return cosine_spaced_points: (N,) ndarray of floats
        This is a ndarray of the N points, ranging from the minimum to the maximum
        value (inclusive), spaced via a cosine function.
    """

    # Find the mean and the amplitude of the cosine function.
    mean = (maximum + minimum) / 2
    amp = (maximum - minimum) / 2

    # Space the points by applying cosine to the output of linspace and return them
    return mean + amp * np.cos(np.linspace(np.pi, 0, n_points, endpoint=endpoint))


# TEST: Consider adding unit tests for this function.
@njit(cache=True, fastmath=False)
def numba_centroid_of_quadrilateral(
    frontLeftPoint_A_a,
    frontRightPoint_A_a,
    backLeftPoint_A_a,
    backRightPoint_A_a,
):
    """This function is used to find the centroid of a quadrilateral. It has been
    optimized for JIT compilation using Numba.

    :param frontLeftPoint_A_a: (3,) ndarray of floats
        This is a ndarray containing the x, y, and z components of the front left
        point of the quadrilateral (in A axes, relative to point a).
    :param frontRightPoint_A_a: (3,) ndarray of floats
        This is a ndarray containing the x, y, and z components of the front right
        point of the quadrilateral (in A axes, relative to point a).
    :param backLeftPoint_A_a: (3,) ndarray of floats
        This is a ndarray containing the x, y, and z components of the back left
        point of the quadrilateral (in A axes, relative to point a).
    :param backRightPoint_A_a: (3,) ndarray of floats
        This is a ndarray containing the x, y, and z components of the back right
        point of the quadrilateral (in A axes, relative to point a).
    :return: (3,) ndarray of floats
        This is a ndarray containing the x, y, and z components of the centroid of the
        quadrilateral (in A axes, relative to point a).
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
    num_steps=25,
    delta_time=0.02,
):
    """This function calculates the location of the streamlines coming off the back
    of the Wings.

    :param solver: steady_horseshoe_vortex_lattice_method
    .SteadyHorseshoeVortexLatticeMethodSolver or
    steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver or
    unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver

        This is the solver for which to calculate the streamlines.

    :param num_steps: int, optional

        This is the integer number of points along each streamline (not including the
        initial points). It can be increased for higher fidelity visuals. The default
        value is 25.

    :param delta_time: float, optional

        This is the time in seconds between each time step It can be decreased for
        higher fidelity visuals or to make the streamlines shorter. Its default value
        is 0.02 seconds.

    :return: None
    """
    # Initialize a ndarray to hold this solver's grid of streamline points (in the
    # first Airplane's geometry axes, relative to the first Airplane's CG).
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
def convert_logging_level_name_to_value(name):
    """This function takes in a string that represents the logging level and returns
    the integer that can be used to set the logger to this level.

    :param name: str
        This is the string representation of the logging level. The options are
        "Debug", "Info", "Warning", "Error", and "Critical".
    :return: int
        This is the integer value that can used to set the appropriate logging level.
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
    stackPanelForces_GP1,
    stackPanelMoments_GP1_CgP1,
):
    """This function uses the forces (in the first Airplane's geometry axes) and
    moments (in the first Airplane's geometry axes, relative to the first Airplane's
    CG) a solver has found on its Panels to find the net loads and associated
    coefficients for each Airplane.

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver

        This is the solver whose loads will be processed.

    :param stackPanelForces_GP1: (N,3) array-like of numbers

        This is an array of the forces (in the first Airplane's geometry axes) on
        each of the solver's Panels. Can be any array-like object (tuple, list,
        or ndarray) with size (N,3) that has numeric elements (int or float). N must
        be equal to the solver's number of Panels. Values are converted to floats
        internally. The units are Newtons.

    :param stackPanelMoments_GP1_CgP1: (N,3) array-like of numbers

        This is an array of the moments (in the first Airplane's geometry axes,
        relative to the first Airplane's CG) on each of the solver's Panels. Can be
        any array-like object (tuple, list, or ndarray) with size (N,3) that has
        numeric elements (int or float). N must be equal to the solver's number of
        Panels. Values are converted to floats internally. The units are Newton-meters.

    :return: None
    """
    if isinstance(
        solver,
        (
            steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
            steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
        ),
    ):
        these_airplanes = solver.airplanes
        this_operating_point = solver.operating_point
    elif isinstance(
        solver,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
    ):
        these_airplanes = solver.current_airplanes
        this_operating_point = solver.current_operating_point
    else:
        raise ValueError(
            f"solver must be a SteadyHorseshoeVortexLatticeMethodSolver, a SteadyRingVortexLatticeMethodSolver, or an UnsteadyRingVortexLatticeMethodSolver."
        )

    stackPanelForces_GP1 = (
        _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
            stackPanelForces_GP1, "stackPanelForces_GP1"
        )
    )
    if stackPanelForces_GP1.shape[0] != solver.num_panels:
        raise ValueError(
            f"The first dimension of stackPanelForces_GP1 must equal solver.num_panels ({solver.num_panels}), got {stackPanelForces_GP1.shape[0]}."
        )

    stackPanelMoments_GP1_CgP1 = (
        _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
            stackPanelMoments_GP1_CgP1, "stackPanelMoments_GP1_CgP1"
        )
    )
    if stackPanelMoments_GP1_CgP1.shape[0] != solver.num_panels:
        raise ValueError(
            f"The first dimension of stackPanelMoments_GP1_CgP1 must equal solver.num_panels ({solver.num_panels}), got {stackPanelMoments_GP1_CgP1.shape[0]}."
        )

    qInf__E = this_operating_point.qInf__E
    T_pas_GP1_CgP1_to_W_CgP1 = this_operating_point.T_pas_GP1_CgP1_to_W_CgP1

    # Iterate through this solver's Panels.
    panel: _panel.Panel
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
    airplane: geometry.airplane.Airplane
    for airplane_num, airplane in enumerate(these_airplanes):
        wing: geometry.wing.Wing
        for wing in airplane.wings:
            panel: _panel.Panel
            for panel in np.ravel(wing.panels):
                stackAirplaneForces_GP1[airplane_num, :] += panel.forces_GP1
                stackAirplaneMoments_GP1_CgP1[airplane_num, :] += panel.moments_GP1_CgP1

    airplane: geometry.airplane.Airplane
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
    airplane: geometry.airplane.Airplane
    for airplane in these_airplanes:
        cFX_W = airplane.forces_W[0] / qInf__E / airplane.s_ref
        cFY_W = airplane.forces_W[1] / qInf__E / airplane.s_ref
        cFZ_W = airplane.forces_W[2] / qInf__E / airplane.s_ref

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
):
    """This function populates a ring vortex solver's attributes with the attributes
    of a given Panel.

    Note: This function doesn't perform any parameter validation.

    :param ring_vortex_solver: SteadyRingVortexLatticeMethodSolver or
    UnsteadyRingVortexLatticeMethodSolver

        This is the solver whose attributes are to be updated.

    :param global_panel_position: int

        This is the position of the Panel with respect to the global array of Panels.

    :param panel: Panel

        This is the Panel whose attributes will be used to update the solver's
        attributes.

    :return: None
    """
    ring_vortex_solver.panels[global_panel_position] = panel
    ring_vortex_solver.stackUnitNormals_GP1[global_panel_position, :] = (
        panel.unitNormal_GP1
    )
    ring_vortex_solver.panel_areas[global_panel_position] = panel.area
    ring_vortex_solver.stackCpp_GP1_CgP1[global_panel_position, :] = panel.Cpp_GP1_CgP1
    ring_vortex_solver.stackBrbrvp_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.right_leg.Slvp_GP1_CgP1
    )
    ring_vortex_solver.stackFrbrvp_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.right_leg.Elvp_GP1_CgP1
    )
    ring_vortex_solver.stackFlbrvp_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.left_leg.Slvp_GP1_CgP1
    )
    ring_vortex_solver.stackBlbrvp_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.left_leg.Elvp_GP1_CgP1
    )
    ring_vortex_solver.stackCblvpr_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.right_leg.Clvp_GP1_CgP1
    )
    ring_vortex_solver.stackRbrv_GP1[global_panel_position, :] = (
        panel.ring_vortex.right_leg.vector_GP1
    )
    ring_vortex_solver.stackCblvpf_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.front_leg.Clvp_GP1_CgP1
    )
    ring_vortex_solver.stackFbrv_GP1[global_panel_position, :] = (
        panel.ring_vortex.front_leg.vector_GP1
    )
    ring_vortex_solver.stackCblvpl_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.left_leg.Clvp_GP1_CgP1
    )
    ring_vortex_solver.stackLbrv_GP1[global_panel_position, :] = (
        panel.ring_vortex.left_leg.vector_GP1
    )
    ring_vortex_solver.stackCblvpb_GP1_CgP1[global_panel_position, :] = (
        panel.ring_vortex.back_leg.Clvp_GP1_CgP1
    )
    ring_vortex_solver.stackBbrv_GP1[global_panel_position, :] = (
        panel.ring_vortex.back_leg.vector_GP1
    )
    ring_vortex_solver.panel_is_trailing_edge[global_panel_position] = (
        panel.is_trailing_edge
    )
    ring_vortex_solver.panel_is_leading_edge[global_panel_position] = (
        panel.is_leading_edge
    )
    ring_vortex_solver.panel_is_right_edge[global_panel_position] = panel.is_right_edge
    ring_vortex_solver.panel_is_left_edge[global_panel_position] = panel.is_left_edge

    # Check if this Panel is on the trailing edge. If so, calculate its streamline
    # seed point (in the first Airplane's geometry axes, relative to the first
    # Airplane's CG) and add it to the solver's ndarray of streamline seed points.
    if panel.is_trailing_edge:
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
):
    """This function finds the vector of freestream-wing influence coefficients
    associated with this solver. These coefficients are the normal velocity (in
    the first Airplane's geometry axes, observed from the Earth frame) at every
    collocation point due solely to the freestream.

    :param steady_solver: steady_horseshoe_vortex_lattice_method
    .SteadyHorseshoeVortexLatticeMethodSolver or
    steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver

        This is the steady solver for which to calculate the freestream-wing influences.

    :return: None
    """
    # Take the batch dot product of the freestream velocity with each of the N
    # panel's normal direction. This is now the problem's (N,) ndarray of
    # freestream-wing influence coefficients.
    steady_solver.stackFreestreamWingInfluences__E = np.einsum(
        "ij,j->i",
        steady_solver.stackUnitNormals_GP1,
        steady_solver.vInf_GP1__E,
    )


# TEST: Consider adding unit tests for this function.
@njit(cache=True, fastmath=False)
def numba_1d_explicit_cross(stackVectors1, stackVectors2):
    """This function takes in two ndarrays, each of which contain N vectors with 3
    components. The function then calculates and returns the cross product of the two
    vectors at each of the N positions.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba. It also doesn't perform any parameter validation.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param stackVectors1: (N,3) ndarray of numbers
        This is the first ndarray of N vectors.
    :param stackVectors2: (N,3) ndarray of numbers
        This is the second ndarray of N vectors.
    :return stackCrossProducts: (N,3) ndarray of numbers
        This is the cross product of the two vectors at each of the N positions.
    """
    stackCrossProducts = np.zeros(stackVectors1.shape)
    for i in range(stackCrossProducts.shape[0]):
        stackCrossProducts[i, 0] = (
            stackVectors1[i, 1] * stackVectors2[i, 2]
            - stackVectors1[i, 2] * stackVectors2[i, 1]
        )
        stackCrossProducts[i, 1] = (
            stackVectors1[i, 2] * stackVectors2[i, 0]
            - stackVectors1[i, 0] * stackVectors2[i, 2]
        )
        stackCrossProducts[i, 2] = (
            stackVectors1[i, 0] * stackVectors2[i, 1]
            - stackVectors1[i, 1] * stackVectors2[i, 0]
        )
    return stackCrossProducts


# TEST: Consider adding unit tests for this function.
@njit(cache=True, fastmath=False)
def interp_between_points(stackStartPoints_A_a, stackEndPoints_A_a, norm_spacings):
    """This function finds the (M,N) grid of points between M pairs of points given
    an array of N normalized spacings.

    :param stackStartPoints_A_a: (M,3) ndarray of floats

        This is the (M,3) ndarray containing the positions of the M start points (in
        A axes, relative to point a). The units are meters.

    :param stackEndPoints_A_a: (M,3) ndarray of floats

        This is the (M,3) ndarray containing the positions of the M end points (in A
        axes, relative to point a). The units are meters.

    :param norm_spacings: (N,) ndarray of floats

        This is the (N,) ndarray of the N normalized spacing values between the start
        and end points. The values are unitless and must be normalized to lie in the
        range from 0.0 to 1.0.

    :return gridInterpolatedPoints_A_a: (M,N,3) ndarray of floats

        This is the (M,N,3) ndarray of the positions of the M*N interpolated points (
        in A axes, relative to point a). The units are meters.
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
