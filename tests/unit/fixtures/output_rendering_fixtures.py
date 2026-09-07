"""This module contains functions to create fixtures for the output rendering tests."""

import numpy as np
import pyvista as pv
import webp

import pterasoftware as ps

from . import (
    airplane_movement_fixtures,
    geometry_fixtures,
    operating_point_fixtures,
    problem_fixtures,
    wing_cross_section_movement_fixtures,
)


def _make_playback_solver(
    airplane_movement: ps.movements.airplane_movement.AirplaneMovement,
    base_operating_point: ps.operating_point.OperatingPoint,
    delta_time: float,
) -> ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver:
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
    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point
        )
    )
    movement = ps.movements.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        delta_time=delta_time,
        num_steps=11,
    )
    unsteady_problem = ps.problems.UnsteadyProblem(
        movement=movement, only_final_results=False
    )
    return ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem
    )


def make_playback_solver_fixture() -> (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
):
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
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
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
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
):
    """Makes a fixture that is a solver whose motion aliases as soon as frames are
    dropped.

    The Wing's rotation has a period of 0.1 seconds, so a time step of 0.01 seconds
    shows 10 frames per cycle while every frame is kept and 5 once every other one is
    dropped, which is below the floor the aliasing warning is held to.

    :return: The unrun UnsteadyRingVortexLatticeMethodSolver of 11 time steps.
    """
    fast_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=geometry_fixtures.make_origin_wing_fixture(),
        wing_cross_section_movements=[
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
            wing_cross_section_movement_fixtures.make_basic_wing_cross_section_movement_fixture(),
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
    fast_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=geometry_fixtures.make_first_airplane_fixture(),
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


def make_static_playback_solver_fixture() -> (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
):
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


def make_image_surface_solver_fixture() -> (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
):
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


def make_loaded_airplanes_fixture() -> tuple[ps.geometry.airplane.Airplane, ...]:
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


def make_placed_airplanes_fixture() -> tuple[ps.geometry.airplane.Airplane, ...]:
    """Makes a fixture that is a tuple of one Airplane placed into a problem.

    A Panel's corner positions in the first Airplane's geometry axes are what the Panel
    surfaces are built from, and those are set when an Airplane is placed into a
    problem. A bare Airplane straight from the geometry fixtures carries only its own
    geometry axis positions, so it cannot stand in here.

    :return: A tuple of one placed Airplane.
    """
    return problem_fixtures.make_basic_steady_problem_fixture().airplanes


def make_formation_airplanes_fixture() -> tuple[ps.geometry.airplane.Airplane, ...]:
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
