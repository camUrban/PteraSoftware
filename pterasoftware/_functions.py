# REFACTOR: I've started refactoring this module.
"""This module contains functions shared by other modules in the pterasoftware
package."""

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


# REFACTOR: I haven't yet started refactoring this function.
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
    :return cosine_spaced_points: 1D array
        This is a 1D array of the points, ranging from the minimum to the maximum
        value (inclusive), spaced via a cosine function.
    """

    # Find the mean and the amplitude of the cosine function.
    mean = (maximum + minimum) / 2
    amp = (maximum - minimum) / 2

    # Space the points by applying cosine to the linspace function. Then return the
    # points.
    cosine_spaced_points = mean + amp * np.cos(
        np.linspace(np.pi, 0, n_points, endpoint=endpoint)
    )
    return cosine_spaced_points


# REFACTOR: I haven't yet started refactoring this function.
@njit(cache=True, fastmath=False)
def numba_centroid_of_quadrilateral(
    front_left_vertex, front_right_vertex, back_left_vertex, back_right_vertex
):
    """This function is used to find the centroid of a quadrilateral. It has been
    optimized for JIT compilation using Numba.

    :param front_left_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the front left
        vertex of the quadrilateral.
    :param front_right_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the front right
        vertex of the quadrilateral.
    :param back_left_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the back left
        vertex of the quadrilateral.
    :param back_right_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the back right
        vertex of the quadrilateral.
    :return: 1D array of floats
        This is an array containing the x, y, and z components of the centroid of the
        quadrilateral.
    """

    x_average = (
        front_left_vertex[0]
        + front_right_vertex[0]
        + back_left_vertex[0]
        + back_right_vertex[0]
    ) / 4
    y_average = (
        front_left_vertex[1]
        + front_right_vertex[1]
        + back_left_vertex[1]
        + back_right_vertex[1]
    ) / 4
    z_average = (
        front_left_vertex[2]
        + front_right_vertex[2]
        + back_left_vertex[2]
        + back_right_vertex[2]
    ) / 4

    return np.array([x_average, y_average, z_average])


# REFACTOR: I haven't yet started refactoring this function.
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
    of the wings.

    This method is vectorized to increase performance.

    :param solver: Solver
        This is the solver object for which to calculate the streamlines.
    :param num_steps: int, optional
        This is the integer number of points along each streamline (not including the
        initial points). It can be increased for higher fidelity visuals. The default
        value is 25.
    :param delta_time: float, optional
        This is the time in seconds between each time step It can be
        decreased for higher fidelity visuals or to make the streamlines shorter.
        Its default value is 0.02 seconds.
    :return: None
    """
    # Initialize an array to hold this solver's matrix of streamline points.
    solver.stackStreamlinePoints_G_Cg = np.expand_dims(
        solver.stackSeedPoints_G_Cg, axis=0
    )

    # Iterate through the streamline steps.
    for step in range(num_steps):
        # Get the last row of streamline points.
        last_row_streamline_points = solver.stackStreamlinePoints_G_Cg[-1, :, :]

        # Add the freestream velocity to the induced velocity to get the total
        # velocity at each of the last row of streamline points.
        total_velocities = solver.calculate_solution_velocity(
            stackP_G_Cg=last_row_streamline_points
        )

        # Interpolate the positions on a new row of streamline points.
        new_row_streamline_points = (
            last_row_streamline_points + total_velocities * delta_time
        )

        # Stack the new row of streamline points to the bottom of the matrix of
        # streamline points.
        solver.stackStreamlinePoints_G_Cg = np.vstack(
            (
                solver.stackStreamlinePoints_G_Cg,
                np.expand_dims(new_row_streamline_points, axis=0),
            )
        )


# REFACTOR: I haven't yet started refactoring this function.
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
        raise Exception("The name of the logging level provided is not a valid option.")


# TEST: Add unit tests for this method.
def process_solver_loads(
    solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    stackPanelForces_G,
    stackPanelMoments_G_Cg,
):
    """This function uses the forces (in geometry axes) and moments (in geometry
    axes, relative to the CG) a solver has found on its Panels to find the net loads
    and associated coefficients for each Airplane.

    :param solver: SteadyHorseshoeVortexLatticeMethodSolver or
    SteadyRingVortexLatticeMethodSolver or UnsteadyRingVortexLatticeMethodSolver

        This is the solver whose loads will be processed.

    :param stackPanelForces_G: (N,3) array-like of numbers

        This is an array of the forces (in geometry axes) on each of the solver's
        Panels. Can be any array-like object (tuple, list, or ndarray) with size (N,
        3) that has numeric elements (int or float). N must be equal to the solver's
        number of Panels. Values are converted to floats internally. The units are
        Newtons.

    :param stackPanelMoments_G_Cg: (N,3) array-like of numbers

        This is an array of the moments (in geometry axes, relative to the CG) on
        each of the solver's Panels. Can be any array-like object (tuple, list,
        or ndarray) with size (N, 3) that has numeric elements (int or float). N must
        be equal to the solver's number of Panels. Values are converted to floats
        internally. The units are Newton-meters.

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

    stackPanelForces_G = (
        _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
            stackPanelForces_G, "stackPanelForces_G"
        )
    )
    if stackPanelForces_G.shape[0] != solver.num_panels:
        raise ValueError(
            f"The first dimension of stackPanelForces_G must equal solver.num_panels ({solver.num_panels}), got {stackPanelForces_G.shape[0]}."
        )

    stackPanelMoments_G_Cg = (
        _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
            stackPanelMoments_G_Cg, "stackPanelMoments_G_Cg"
        )
    )
    if stackPanelMoments_G_Cg.shape[0] != solver.num_panels:
        raise ValueError(
            f"The first dimension of stackPanelMoments_G_Cg must equal solver.num_panels ({solver.num_panels}), got {stackPanelMoments_G_Cg.shape[0]}."
        )

    qInf__E = this_operating_point.qInf__E
    T_pas_G_Cg_to_W_Cg = this_operating_point.T_pas_G_Cg_to_W_Cg

    # Iterate through this solver's Panels.
    panel: _panel.Panel
    for panel_num, panel in enumerate(solver.panels):
        theseForces_G = stackPanelForces_G[panel_num, :]
        theseMoments_G_Cg = stackPanelMoments_G_Cg[panel_num, :]

        theseForces_W = _transformations.apply_T_to_vectors(
            T_pas_G_Cg_to_W_Cg, theseForces_G, has_point=False
        )
        theseMoments_W_Cg = _transformations.apply_T_to_vectors(
            T_pas_G_Cg_to_W_Cg, theseMoments_G_Cg, has_point=True
        )

        # Update this Panel's loads.
        panel.forces_G = theseForces_G
        panel.moments_G_Cg = theseMoments_G_Cg
        panel.forces_W = theseForces_W
        panel.moments_W_Cg = theseMoments_W_Cg

        # Update this Panel's load coefficients.
        panel.update_coefficients(qInf__E)

    # Initialize ndarrays to hold each Airplane's loads.
    stackAirplaneForces_G = np.zeros((solver.num_airplanes, 3), dtype=float)
    stackAirplaneMoments_G_Cg = np.zeros((solver.num_airplanes, 3), dtype=float)

    # Iterate through each Airplane and find the total loads on each by summing up
    # the contributions from its Panels.
    airplane: geometry.airplane.Airplane
    for airplane_num, airplane in enumerate(these_airplanes):
        wing: geometry.wing.Wing
        for wing in airplane.wings:
            panel: _panel.Panel
            for panel in np.ravel(wing.panels):
                stackAirplaneForces_G[airplane_num, :] += panel.forces_G
                stackAirplaneMoments_G_Cg[airplane_num, :] += panel.moments_G_Cg

    airplane: geometry.airplane.Airplane
    for airplane_num, airplane in enumerate(these_airplanes):
        airplane.forces_W = _transformations.apply_T_to_vectors(
            T_pas_G_Cg_to_W_Cg,
            stackAirplaneForces_G[airplane_num],
            has_point=False,
        )
        airplane.moments_W_Cg = _transformations.apply_T_to_vectors(
            T_pas_G_Cg_to_W_Cg,
            stackAirplaneMoments_G_Cg[airplane_num],
            has_point=True,
        )

    # Iterate through the Airplanes and calculate each one's load coefficients.
    airplane: geometry.airplane.Airplane
    for airplane in these_airplanes:
        cFX_W = airplane.forces_W[0] / qInf__E / airplane.s_ref
        cFY_W = airplane.forces_W[1] / qInf__E / airplane.s_ref
        cFZ_W = airplane.forces_W[2] / qInf__E / airplane.s_ref

        cMX_W_Cg = airplane.moments_W_Cg[0] / qInf__E / airplane.s_ref / airplane.b_ref
        cMY_W_Cg = airplane.moments_W_Cg[1] / qInf__E / airplane.s_ref / airplane.c_ref
        cMZ_W_Cg = airplane.moments_W_Cg[2] / qInf__E / airplane.s_ref / airplane.b_ref

        # Populate this Airplane's load coefficients.
        airplane.forceCoefficients_W = np.array([cFX_W, cFY_W, cFZ_W])
        airplane.momentCoefficients_W_Cg = np.array([cMX_W_Cg, cMY_W_Cg, cMZ_W_Cg])


# REFACTOR: I haven't yet started refactoring this function.
def update_ring_vortex_solvers_panel_attributes(
    ring_vortex_solver: (
        steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    ),
    global_panel_position: int,
    panel: _panel.Panel,
):
    """This function populates a ring vortex solver's attributes with the attributes
    of a given panel.

    :param ring_vortex_solver: SteadySolver or UnsteadySolver
        This is the solver object whose attributes are to be updated. It should be a
        solver that uses ring vortices.
    :param global_panel_position: int
        This is the position of the panel with respect to the global array of all
        panels.
    :param panel: Panel
        This is the panel object whose attributes will be used to update the solver's
        attributes.
    :return: None
    """

    # Update the solver's list of attributes with this panel's attributes.
    ring_vortex_solver.panels[global_panel_position] = panel
    ring_vortex_solver.stackUnitNormals_G[global_panel_position, :] = panel.unitNormal_G
    ring_vortex_solver.panel_areas[global_panel_position] = panel.area
    ring_vortex_solver.stackCpp_G_Cg[global_panel_position, :] = panel.Cpp_G_Cg
    ring_vortex_solver.stackBrbrvp_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.right_leg.Slvp_G_Cg
    )
    ring_vortex_solver.stackFrbrvp_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.right_leg.Elvp_G_Cg
    )
    ring_vortex_solver.stackFlbrvp_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.left_leg.Slvp_G_Cg
    )
    ring_vortex_solver.stackBlbrvp_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.left_leg.Elvp_G_Cg
    )
    ring_vortex_solver.stackCblvpr_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.right_leg.Clvp_G_Cg
    )
    ring_vortex_solver.stackRbrv_G[global_panel_position, :] = (
        panel.ring_vortex.right_leg.vector_G
    )
    ring_vortex_solver.stackCblvpf_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.front_leg.Clvp_G_Cg
    )
    ring_vortex_solver.stackFbrv_G[global_panel_position, :] = (
        panel.ring_vortex.front_leg.vector_G
    )
    ring_vortex_solver.stackCblvpl_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.left_leg.Clvp_G_Cg
    )
    ring_vortex_solver.stackLbrv_G[global_panel_position, :] = (
        panel.ring_vortex.left_leg.vector_G
    )
    ring_vortex_solver.stackCblvpb_G_Cg[global_panel_position, :] = (
        panel.ring_vortex.back_leg.Clvp_G_Cg
    )
    ring_vortex_solver.stackBbrv_G[global_panel_position, :] = (
        panel.ring_vortex.back_leg.vector_G
    )
    ring_vortex_solver.panel_is_trailing_edge[global_panel_position] = (
        panel.is_trailing_edge
    )
    ring_vortex_solver.panel_is_leading_edge[global_panel_position] = (
        panel.is_leading_edge
    )
    ring_vortex_solver.panel_is_right_edge[global_panel_position] = panel.is_right_edge
    ring_vortex_solver.panel_is_left_edge[global_panel_position] = panel.is_left_edge

    # Check if this panel is on the trailing edge. If it is, calculate its
    # streamline seed point and add it to the solver's # array of seed points.
    if panel.is_trailing_edge:
        ring_vortex_solver.stackSeedPoints_G_Cg = np.vstack(
            (
                ring_vortex_solver.stackSeedPoints_G_Cg,
                panel.Blpp_G_Cg + 0.5 * (panel.Brpp_G_Cg - panel.Blpp_G_Cg),
            )
        )


# REFACTOR: I haven't yet started refactoring this function.
def calculate_steady_freestream_wing_influences(
    steady_solver: (
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
    ),
):
    """This function finds the vector of freestream-wing influence coefficients
    associated with this problem.

    :param steady_solver:
    :return:
    """
    # Take the batch dot product of the freestream velocity with each panel's
    # normal direction. This is now the problem's 1D array of freestream-wing
    # influence coefficients.
    steady_solver.stackFreestreamWingInfluences__E = np.einsum(
        "ij,j->i",
        steady_solver.stackUnitNormals_G,
        steady_solver.vInf_G__E,
    )


# REFACTOR: I haven't yet started refactoring this function.
@njit(cache=True, fastmath=False)
def numba_1d_explicit_cross(vectors_1, vectors_2):
    """This function takes in two arrays, each of which contain N vectors of 3
    components. The function then calculates and returns the cross product of the two
    vectors at each position.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param vectors_1: array of floats of size (N x 3)
        This is the first array of N vectors.
    :param vectors_2: array of floats of size (N x 3)
        This is the second array of N vectors.
    :return crosses: array of floats of size (N x 3)
        This is the cross product of the two inputted vectors at each of the N
        positions.
    """
    crosses = np.zeros(vectors_1.shape)
    for i in range(crosses.shape[0]):
        crosses[i, 0] = (
            vectors_1[i, 1] * vectors_2[i, 2] - vectors_1[i, 2] * vectors_2[i, 1]
        )
        crosses[i, 1] = (
            vectors_1[i, 2] * vectors_2[i, 0] - vectors_1[i, 0] * vectors_2[i, 2]
        )
        crosses[i, 2] = (
            vectors_1[i, 0] * vectors_2[i, 1] - vectors_1[i, 1] * vectors_2[i, 0]
        )
    return crosses


# REFACTOR: I haven't yet started refactoring this function.
@njit(cache=True, fastmath=False)
def interp_between_points(start_points, end_points, norm_spacings):
    """This function finds the MxN points between M pairs of points in 3D space given
    an array of N normalized spacings.

    :param start_points: (M, 3) array of floats
        This is the (M, 3) array containing the coordinates of the M starting points.
        The units are meters.
    :param end_points: (M, 3) array of floats
        This is the (M, 3) array containing the coordinates of the M ending points.
        The units are meters.
    :param norm_spacings: (N,) array of floats
        This is the (N,) array of the N spacings between the starting points and
        ending points. The values are unitless and must be normalized from 0 to 1.
    :return points: (M, N, 3) array of floats
        This is the (M, N, 3) array of the coordinates of the MxN interpolated
        points. The units are meters.
    """
    m = start_points.shape[0]
    n = norm_spacings.size

    points = np.zeros((m, n, 3))

    for i in range(m):
        start_point = start_points[i, :]
        end_point = end_points[i, :]

        vector = end_point - start_point

        for j in range(n):
            norm_spacing = norm_spacings[j]

            spacing = norm_spacing * vector

            points[i, j, :] = start_point + spacing

    return points
