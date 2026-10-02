"""This module contains functions to create fixtures for the vector export tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware import _vector_export


def make_top_down_camera_fixture() -> _vector_export.VectorCamera:
    """Makes a VectorCamera that looks down the negative z direction at the diagram
    origin, with the y direction up, into a 200 by 100 pixel window that spans 10.0 m
    vertically.

    Through this camera, the diagram axes' x and y directions run right and up on screen
    at 10.0 pixels per meter, the diagram origin is shown at the window's center, and a
    point's depth is -10.0 pixels per meter of its z position, so a larger z is nearer.

    :return: The VectorCamera.
    """
    return _vector_export.VectorCamera(
        position_D_Do=np.array([0.0, 0.0, 10.0], dtype=float),
        focalPoint_D_Do=np.zeros(3, dtype=float),
        viewUp_D=np.array([0.0, 1.0, 0.0], dtype=float),
        parallel_scale=5.0,
        window_width=200,
        window_height=100,
    )


def make_crossing_triangles_fixture() -> list[tuple[np.ndarray, int]]:
    """Makes two triangles that share a footprint on screen and cross each other halfway
    across it, so that neither is wholly in front of the other.

    The first triangle's depth rises with x, and the second's falls with it, so the
    first is nearer on the left half of the footprint and the second is nearer on the
    right half.

    :return: A list of two polygons, each a tuple of a (3,3) ndarray of floats holding
        its vertices (in display coordinates) and its int ID, which is its index.
    """
    return [
        (
            np.array(
                [[0.0, 0.0, 0.0], [100.0, 0.0, 100.0], [0.0, 100.0, 0.0]],
                dtype=float,
            ),
            0,
        ),
        (
            np.array(
                [[0.0, 0.0, 100.0], [100.0, 0.0, 0.0], [0.0, 100.0, 100.0]],
                dtype=float,
            ),
            1,
        ),
    ]


def make_random_triangles_fixture() -> list[tuple[np.ndarray, int]]:
    """Makes 40 seeded random triangles, which overlap and cross one another throughout
    a 100 by 100 pixel region of the screen.

    :return: A list of 40 polygons, each a tuple of a (3,3) ndarray of floats holding
        its vertices (in display coordinates) and its int ID, which is its index.
    """
    rng = np.random.default_rng(0)
    return [
        (rng.uniform(0.0, 100.0, size=(3, 3)), triangle_id) for triangle_id in range(40)
    ]


def make_flat_triangle_fixture() -> np.ndarray:
    """Makes a triangle at a constant depth of 50.0 pixels, whose footprint is the lower
    left half of the square from (0.0, 0.0) to (100.0, 100.0) on screen.

    :return: A (1,3,3) ndarray of floats holding the triangle's vertices (in display
        coordinates).
    """
    return np.array(
        [[[0.0, 0.0, 50.0], [100.0, 0.0, 50.0], [0.0, 100.0, 50.0]]], dtype=float
    )


def make_cube_points_fixture() -> np.ndarray:
    """Makes the corners of a cube with 2.0 m sides, centered on the diagram origin.

    :return: A (8,3) ndarray of floats holding the cube's corners (in diagram axes,
        relative to the diagram origin). Corner i has the x, y, and z positions -1.0 or
        1.0 according to bits 0, 1, and 2 of i. The units are in meters.
    """
    return np.array(
        [
            [
                1.0 if corner_id & 1 else -1.0,
                1.0 if corner_id & 2 else -1.0,
                1.0 if corner_id & 4 else -1.0,
            ]
            for corner_id in range(8)
        ],
        dtype=float,
    )


def make_cube_faces_fixture() -> list[np.ndarray]:
    """Makes the faces of the cube make_cube_points_fixture makes.

    :return: A list of six (4,) ndarrays of ints, each holding the indices of a face's
        corners in order around it. The faces are the ones at z = -1.0, z = 1.0, y =
        -1.0, y = 1.0, x = -1.0, and x = 1.0, in that order.
    """
    return [
        np.array(face, dtype=int)
        for face in (
            [0, 1, 3, 2],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 3, 7, 6],
            [0, 2, 6, 4],
            [1, 3, 7, 5],
        )
    ]
