"""This module contains functions to create fixtures for the output rendering tests."""

import numpy as np
import pyvista as pv
import webp

import pterasoftware as ps
from pterasoftware import _transformations

from . import (
    airplane_movement_fixtures,
    geometry_fixtures,
    operating_point_fixtures,
    problem_fixtures,
    wing_cross_section_movement_fixtures,
)


def _make_playback_solver(
    airplane_movement: ps.AirplaneMovement,
    base_operating_point: ps.OperatingPoint,
    delta_time: float,
) -> ps.UnsteadyRingVortexLatticeMethodSolver:
    """Makes an unrun solver whose time step characteristics are set explicitly.

    The playback arithmetic reads only the solver's delta_time, its num_steps, and its
    Movement, so every playback fixture is one of these built over a different Movement.
    Both time step characteristics are passed rather than derived, since the resolved
    stride and frame rate are quotients of the time step.

    :param airplane_movement: The AirplaneMovement the solver's Movement is built over.
    :param base_operating_point: The OperatingPoint the solver's OperatingPointMovement
        is built over.
    :param delta_time: The time step's length in seconds.
    :return: The unrun UnsteadyRingVortexLatticeMethodSolver.
    """
    operating_point_movement = ps.OperatingPointMovement(
        base_operating_point=base_operating_point
    )
    movement = ps.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        delta_time=delta_time,
        num_steps=11,
    )
    unsteady_problem = ps.UnsteadyProblem(movement=movement, only_final_results=False)
    return ps.UnsteadyRingVortexLatticeMethodSolver(unsteady_problem)


def make_playback_solver_fixture() -> ps.UnsteadyRingVortexLatticeMethodSolver:
    """Makes a fixture that is a solver whose playback arithmetic comes out even.

    A time step of 0.01 seconds carries 0.5 seconds of simulation per second of playback
    at the maximum frame rate, so the resolved stride and frame rate are whole numbers
    at the speeds the tests ask for. The motion's shortest period is 2.0 seconds, which
    is slow enough that dropping frames never trips the aliasing warning.

    :return: The unrun UnsteadyRingVortexLatticeMethodSolver of 11 time steps.
    """
    return _make_playback_solver(
        airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
        operating_point_fixtures.make_basic_operating_point_fixture(),
        0.01,
    )


def make_long_step_playback_solver_fixture() -> (
    ps.UnsteadyRingVortexLatticeMethodSolver
):
    """Makes a fixture that is a solver the maximum frame rate can play at true speed.

    A time step of 0.05 seconds needs only 20 frames per second of playback to run at
    true speed, which is where the default speed stops being held down by the maximum
    frame rate.

    :return: The unrun UnsteadyRingVortexLatticeMethodSolver of 11 time steps.
    """
    return _make_playback_solver(
        airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
        operating_point_fixtures.make_basic_operating_point_fixture(),
        0.05,
    )


def make_fast_motion_playback_solver_fixture() -> (
    ps.UnsteadyRingVortexLatticeMethodSolver
):
    """Makes a fixture that is a solver whose motion aliases as soon as frames are
    dropped.

    The Wing's rotation has a period of 0.1 seconds, so a time step of 0.01 seconds
    shows 10 frames per cycle while every frame is kept and 5 once every other one is
    dropped, which is below the floor the aliasing warning is held to.

    :return: The unrun UnsteadyRingVortexLatticeMethodSolver of 11 time steps.
    """
    # Build the base Airplane first, then build the movements around its own Wing and
    # WingCrossSections.
    base_airplane = geometry_fixtures.make_origin_airplane_fixture()
    base_wing = base_airplane.wings[0]
    fast_wing_movement = ps.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=[
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(
                base_wing.wing_cross_sections[0]
            ),
            wing_cross_section_movement_fixtures.make_basic_wing_cross_section_movement_fixture(
                base_wing.wing_cross_sections[1]
            ),
        ],
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(5.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.1, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )
    fast_airplane_movement = ps.AirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=[fast_wing_movement],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )
    return _make_playback_solver(
        fast_airplane_movement,
        operating_point_fixtures.make_basic_operating_point_fixture(),
        0.01,
    )


def make_static_playback_solver_fixture() -> ps.UnsteadyRingVortexLatticeMethodSolver:
    """Makes a fixture that is a solver whose geometry never moves.

    A static geometry has no motion to alias, which is one of the three cases the
    aliasing warning is skipped for.

    :return: The unrun UnsteadyRingVortexLatticeMethodSolver of 11 time steps.
    """
    return _make_playback_solver(
        airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        operating_point_fixtures.make_basic_operating_point_fixture(),
        0.01,
    )


def make_image_surface_solver_fixture() -> ps.UnsteadyRingVortexLatticeMethodSolver:
    """Makes a fixture that is a solver whose OperatingPoint defines an image surface.

    The surface is the ground, 10.0 meters below the first Airplane's CG, so the
    geometry and its reflected copy sit well apart and a bounding box spanning both is
    easy to tell from one spanning the geometry alone.

    :return: The unrun UnsteadyRingVortexLatticeMethodSolver of 11 time steps.
    """
    return _make_playback_solver(
        airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
        operating_point_fixtures.make_with_ground_surface_operating_point_fixture(),
        0.01,
    )


def make_loaded_airplanes_fixture() -> tuple[ps.Airplane, ...]:
    """Makes a fixture that is a tuple of one Airplane whose Panels carry known loads.

    A solver sets each Panel's loads while it runs, which no unit test does, so they are
    set here instead. Each Panel's load rises with its position in the Wing's unraveled
    ndarray of Panels, and the three components differ from one another, so a test can
    tell one Panel's scalar from another's and one component from another.

    :return: A tuple of one Airplane whose every Panel has its forces (in wind axes)
        set.
    """
    airplane = geometry_fixtures.make_basic_airplane_fixture()
    for wing in airplane.wings:
        assert wing.panels is not None
        panels = np.ravel(wing.panels)
        for panel_num, panel in enumerate(panels):
            panel.forces_W = np.array(
                [
                    -(panel_num + 1.0),
                    2.0 * (panel_num + 1.0),
                    -3.0 * (panel_num + 1.0),
                ],
                dtype=float,
            )
    return (airplane,)


def make_placed_airplanes_fixture() -> tuple[ps.Airplane, ...]:
    """Makes a fixture that is a tuple of one Airplane placed into a problem.

    A Panel's corner positions in the first Airplane's geometry axes are what the Panel
    surfaces are built from, and those are set when an Airplane is placed into a
    problem. A bare Airplane straight from the geometry fixtures carries only its own
    geometry axis positions, so it cannot stand in here.

    :return: A tuple of one placed Airplane.
    """
    return problem_fixtures.make_basic_steady_problem_fixture().airplanes


def make_formation_airplanes_fixture() -> tuple[ps.Airplane, ...]:
    """Makes a fixture that is a tuple of two Airplanes placed into one problem.

    Both are placed in the first Airplane's geometry axes, which is what lets a
    formation's Airplanes share one mesh.

    :return: A tuple of two placed Airplanes.
    """
    return problem_fixtures.make_multi_airplane_steady_problem_fixture().airplanes


def make_streamline_points_fixture() -> np.ndarray:
    """Makes a fixture that is the points along each of three streamlines.

    Every point holds a different value, so a test can tell which streamline a point
    belongs to and where along that streamline it sits.

    :return: A (4,3,3) ndarray of floats representing the points along each streamline
        (in the first Airplane's geometry axes, relative to the first Airplane's CG).
    """
    return np.arange(4 * 3 * 3, dtype=float).reshape(4, 3, 3)


def make_cube_mesh_fixture() -> pv.PolyData:
    """Makes a fixture that is a cube of side 2.0 meters centered on the origin.

    Its bounding box corners sit 1.0 meter from the origin along each axis, so the
    extent a camera has to frame is 1.0 whichever pair of axes it measures.

    :return: The cube's PolyData mesh.
    """
    return pv.Cube(center=(0.0, 0.0, 0.0), x_length=2.0, y_length=2.0, z_length=2.0)


def make_deep_box_mesh_fixture() -> pv.PolyData:
    """Makes a fixture that is a box stretched along the z axis.

    It spans the same 2.0 meters as the cube across the x and y axes and 20.0 meters
    along the z axis, so a camera looking down the z axis frames it exactly as it frames
    the cube.

    :return: The box's PolyData mesh.
    """
    return pv.Cube(center=(0.0, 0.0, 0.0), x_length=2.0, y_length=2.0, z_length=20.0)


def make_geometry_bounds_fixture() -> tuple[float, float, float, float, float, float]:
    """Makes a fixture that is a bounding box for sizing an image surface plane.

    The box is 2.0 by 4.0 by 6.0 meters, so its diagonal is sqrt(56.0) meters, and it is
    centered on the origin, so the plane's center is the origin's projection onto the
    surface.

    :return: The (xmin, xmax, ymin, ymax, zmin, zmax) bounding box in meters.
    """
    return -1.0, 1.0, -2.0, 2.0, -3.0, 3.0


def make_outlier_scalars_fixture() -> np.ndarray:
    """Makes a fixture that is a set of scalars with one far outlier.

    A hundred Panels carry 1.0 and one carries 101.0, which puts the largest value far
    more than three standard deviations above the mean. Without the sigma bound, that
    one Panel would take the whole upper end of the color map.

    :return: A (101,) ndarray of floats representing the scalar value at each Panel.
    """
    scalars = np.ones(101, dtype=float)
    scalars[-1] = 101.0
    return scalars


def make_animation_frames_fixture(
    num_frames: int, width: int = 16, height: int = 12
) -> list[webp.Image.Image]:
    """Makes a fixture that is a list of frames for an animation.

    Each frame is filled with random colors from a fixed seed, so the frames differ from
    one another and the encoder cannot merge any of them, while the list is the same on
    every call.

    :param num_frames: The number of frames to make.
    :param width: The width of each frame in pixels. The default is 16.
    :param height: The height of each frame in pixels. The default is 12.
    :return: A list of num_frames Images with transparent backgrounds.
    """
    rng = np.random.default_rng(0)
    return [
        webp.Image.fromarray(
            rng.integers(0, 255, size=(height, width, 4), dtype=np.uint8)
        )
        for _ in range(num_frames)
    ]


def make_no_vortex_points_fixture() -> np.ndarray:
    """Makes a fixture that is an empty stack of vortex points.

    It stands in for every corner point stack and unit normal stack of a kind of vortex
    that a solver doesn't place, as the solvers themselves pass when drawing their
    diagrams.

    :return: A (0,3) ndarray of floats.
    """
    return np.empty((0, 3), dtype=float)


def make_square_ring_vortex_fixture() -> (
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
):
    """Makes a fixture that is a single ring vortex with sides of 1.0 meter.

    The ring vortex lies in the xy plane, with its front leg along the y axis and its
    back leg 1.0 meter downstream of it, and its Panel's unit normal points along the
    positive z direction. Its center is (0.5, 0.5, 0.0), so the simplified ring vortex's
    corners sit 0.45 meters from its center along the x and y axes.

    :return: A tuple of five (1,3) ndarrays of floats. In order, they hold the ring
        vortex's front right, front left, back left, and back right points (in diagram
        axes, relative to the diagram origin), and its Panel's unit normal (in diagram
        axes). The units of the points are in meters.
    """
    return (
        np.array([[0.0, 1.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 0.0]], dtype=float),
        np.array([[1.0, 0.0, 0.0]], dtype=float),
        np.array([[1.0, 1.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 1.0]], dtype=float),
    )


def make_neighboring_ring_vortices_fixture() -> (
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
):
    """Makes a fixture that is two ring vortices with sides of 1.0 meter that share a
    leg.

    Both ring vortices lie in the xy plane, side by side along the y axis, so the first
    one's right leg is the second one's left leg. Their Panels' unit normals point along
    the positive z direction.

    :return: A tuple of five (2,3) ndarrays of floats. In order, they hold each ring
        vortex's front right, front left, back left, and back right points (in diagram
        axes, relative to the diagram origin), and each one's Panel's unit normal (in
        diagram axes). The units of the points are in meters.
    """
    return (
        np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float),
        np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], dtype=float),
        np.array([[1.0, 1.0, 0.0], [1.0, 2.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=float),
    )


def make_straight_horseshoe_vortex_fixture() -> (
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
):
    """Makes a fixture that is a single horseshoe vortex whose trailing legs run along
    the positive x direction.

    The finite leg runs along the y axis from (0.0, 1.0, 0.0) to the origin, so the
    vortices' bounding box is flat in x and every trailing leg ends exactly the overhang
    downstream of the finite leg. The back points sit 20.0 meters downstream, as they
    would for a solver's wing of span 1.0 meter, and the Panel's unit normal points
    along the positive z direction.

    :return: A tuple of five (1,3) ndarrays of floats. In order, they hold the horseshoe
        vortex's front right, front left, back left, and back right points (in diagram
        axes, relative to the diagram origin), and its Panel's unit normal (in diagram
        axes). The units of the points are in meters.
    """
    return (
        np.array([[0.0, 1.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 0.0]], dtype=float),
        np.array([[20.0, 0.0, 0.0]], dtype=float),
        np.array([[20.0, 1.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 1.0]], dtype=float),
    )


def make_staggered_horseshoe_vortices_fixture() -> (
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
):
    """Makes a fixture that is two horseshoe vortices whose neighboring trailing legs
    lie on top of each other but start at different points.

    The first horseshoe vortex's finite leg runs along the y axis from (0.0, 1.0, 0.0)
    to the origin. The second one's runs from (0.3, 2.0, 0.0) to (0.3, 1.0, 0.0), so its
    left trailing leg starts 0.3 meters downstream along the first one's right trailing
    leg. Every trailing leg runs along the positive x direction, and both Panels' unit
    normals point along the positive z direction.

    :return: A tuple of five (2,3) ndarrays of floats. In order, they hold each
        horseshoe vortex's front right, front left, back left, and back right points (in
        diagram axes, relative to the diagram origin), and each one's Panel's unit
        normal (in diagram axes). The units of the points are in meters.
    """
    return (
        np.array([[0.0, 1.0, 0.0], [0.3, 2.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 0.0], [0.3, 1.0, 0.0]], dtype=float),
        np.array([[20.0, 0.0, 0.0], [20.3, 1.0, 0.0]], dtype=float),
        np.array([[20.0, 1.0, 0.0], [20.3, 2.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=float),
    )


def make_slanted_horseshoe_vortex_fixture() -> (
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
):
    """Makes a fixture that is a single horseshoe vortex whose trailing legs run
    diagonally in the xz plane.

    The finite leg runs along the y axis from (0.0, 1.0, 0.0) to the origin, like the
    straight horseshoe vortex's, but the back points sit 20.0 meters along both the x
    and z axes, so the trailing legs run along the unit vector (1.0, 0.0, 1.0) /
    sqrt(2.0). The vortices' bounding box is flat along that direction, so every
    trailing leg still ends exactly the overhang past the finite leg. The Panel's unit
    normal points along the positive z direction.

    :return: A tuple of five (1,3) ndarrays of floats. In order, they hold the horseshoe
        vortex's front right, front left, back left, and back right points (in diagram
        axes, relative to the diagram origin), and its Panel's unit normal (in diagram
        axes). The units of the points are in meters.
    """
    return (
        np.array([[0.0, 1.0, 0.0]], dtype=float),
        np.array([[0.0, 0.0, 0.0]], dtype=float),
        np.array([[20.0, 0.0, 20.0]], dtype=float),
        np.array([[20.0, 1.0, 20.0]], dtype=float),
        np.array([[0.0, 0.0, 1.0]], dtype=float),
    )


def make_offset_axes_transformation_fixture() -> np.ndarray:
    """Makes a fixture that is the transformation of an axes set turned 90 degrees about
    the z axis and placed away from the diagram origin.

    The axes set's x basis direction points along the positive y direction, its y basis
    direction points along the negative x direction, and its z basis direction points
    along the positive z direction (all in diagram axes). Its point sits at (1.0, 2.0,
    3.0) (in diagram axes, relative to the diagram origin).

    :return: A (4,4) ndarray of floats representing the passive transformation matrix
        which maps in homogeneous coordinates from the axes set, relative to its point,
        to diagram axes, relative to the diagram origin.
    """
    return np.array(
        [
            [0.0, -1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def make_x_turned_axes_transformation_fixture() -> np.ndarray:
    """Makes a fixture that is the transformation of an axes set turned 90 degrees about
    the x axis, at the diagram origin.

    The axes set's x basis direction points along the positive x direction, its y basis
    direction points along the positive z direction, and its z basis direction points
    along the negative y direction (all in diagram axes). Drawn alongside an axes set
    aligned with diagram axes at the same point, its x arrow coincides with that axes
    set's x arrow and its y arrow coincides with that axes set's z arrow, while its z
    arrow coincides with none of them.

    :return: A (4,4) ndarray of floats representing the passive transformation matrix
        which maps in homogeneous coordinates from the axes set, relative to its point,
        to diagram axes, relative to the diagram origin.
    """
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def make_reversed_axes_transformation_fixture() -> np.ndarray:
    """Makes a fixture that is the transformation of an axes set turned 180 degrees
    about the z axis, at the diagram origin.

    The axes set's x and y basis directions point along the negative x and negative y
    directions, and its z basis direction points along the positive z direction (all in
    diagram axes), so its x and y arrows point left and down on screen when viewed down
    the negative z direction.

    :return: A (4,4) ndarray of floats representing the passive transformation matrix
        which maps in homogeneous coordinates from the axes set, relative to its point,
        to diagram axes, relative to the diagram origin.
    """
    return np.array(
        [
            [-1.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def make_oblique_axes_transformation_fixture() -> np.ndarray:
    """Makes a fixture that is the transformation of an axes set turned obliquely and
    placed away from the diagram origin.

    The axes set is turned by intrinsic xyz angles of (30.0, -50.0, 120.0) degrees, so
    none of its basis directions lie along a diagram axis, and each of them has a
    negative x component (in diagram axes). Its point sits at (0.5, -1.0, 2.0) (in
    diagram axes, relative to the diagram origin).

    :return: A (4,4) ndarray of floats representing the passive transformation matrix
        which maps in homogeneous coordinates from the axes set, relative to its point,
        to diagram axes, relative to the diagram origin.
    """
    # Turning first and then translating within diagram axes leaves the turned basis
    # directions in the first three columns and the translation in the last column.
    rot_T_act = _transformations.generate_rot_T(
        np.array([30.0, -50.0, 120.0], dtype=float),
        passive=False,
        intrinsic=True,
        order="xyz",
    )
    trans_T_act = _transformations.generate_trans_T(
        np.array([0.5, -1.0, 2.0], dtype=float), passive=False
    )
    return _transformations.compose_T_act(rot_T_act, trans_T_act)


def make_cross_directions_fixture() -> np.ndarray:
    """Makes a fixture that is the directions of an extra point's cross's two arms.

    The arms lie along the diagonals of the xy plane, so neither one lies along an axis,
    and a test can tell a cross oriented by these directions from one oriented by any
    axis.

    :return: A (2,3) ndarray of floats whose rows hold the unit vectors (in diagram
        axes) along which the cross's two arms lie.
    """
    crossDirections_D: np.ndarray = np.array(
        [[1.0, 1.0, 0.0], [1.0, -1.0, 0.0]],
        dtype=float,
    ) / np.sqrt(2.0)
    return crossDirections_D
