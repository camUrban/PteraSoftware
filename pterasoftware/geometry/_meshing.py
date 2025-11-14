"""Contains the function for meshing Wings."""

from __future__ import annotations

import numpy as np

from . import airfoil as airfoil_mod
from . import wing_cross_section as wing_cross_section_mod
from . import wing as wing_mod
from .. import _functions
from .. import _panel
from .. import _transformations


def mesh_wing(wing: wing_mod.Wing) -> None:
    """Takes in a Wing, creates a quadrilateral mesh of its geometry, and then populates
    its array of Panels with the mesh data.

    **Citation:**

    Adapted from: vlm3.make_panels in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 05/01/2020

    :param wing: The Wing to be meshed.
    :return: None
    """
    # Gather this Wing's attributes
    wing_cross_sections = wing.wing_cross_sections
    symmetry_type = wing.symmetry_type
    symmetryNormal_G = wing.symmetryNormal_G
    symmetryPoint_G_Cg = wing.symmetryPoint_G_Cg
    num_chordwise_panels = wing.num_chordwise_panels
    chordwise_spacing = wing.chordwise_spacing
    T_pas_Wn_Ler_to_G_Cg = wing.T_pas_Wn_Ler_to_G_Cg
    children_T_pas_Wcs_Lp_to_Wn_Ler = wing.children_T_pas_Wcs_Lp_to_Wn_Ler

    # Define the number of points.
    num_chordwise_coordinates = num_chordwise_panels + 1

    # Get the chordwise coordinates.
    if chordwise_spacing == "uniform":
        chordwise_coordinates = np.linspace(0, 1, num_chordwise_coordinates)
    else:
        chordwise_coordinates = _functions.cosspace(0, 1, num_chordwise_coordinates)

    # Get the number of WingCrossSections and wing sections.
    num_wing_cross_sections = len(wing_cross_sections)
    num_wing_sections = num_wing_cross_sections - 1

    # Initialize an empty array that will hold the Panels of this Wing. It currently
    # has 0 columns and M rows, where M is the number of the Wing's chordwise Panels.
    wing_panels: np.ndarray = np.empty((num_chordwise_panels, 0), dtype=object)

    # Make the Panels for each wing section.
    for wing_section_num in range(num_wing_sections):
        # Define variables to hold the indices of this wing section's inner
        # WingCrossSection.
        inner_wing_cross_section_num = wing_section_num

        # Define the relevant WingCrossSections.
        inner_wing_cross_section = wing_cross_sections[inner_wing_cross_section_num]
        outer_wing_cross_section = wing_cross_sections[inner_wing_cross_section_num + 1]

        # Get the homogeneous passive transformation matrices from the inner and
        # outer WingCrossSection's axes relative to their respective leading points
        # to wing axes relative to the leading edge root point.
        T_pas_Wcsi_Lpi_to_Wn_Ler = children_T_pas_Wcs_Lp_to_Wn_Ler[
            inner_wing_cross_section_num
        ]
        T_pas_Wcso_Lpo_to_Wn_Ler = children_T_pas_Wcs_Lp_to_Wn_Ler[
            inner_wing_cross_section_num + 1
        ]

        # Define the Airfoils at each WingCrossSection and modify them with their
        # control surfaces.
        inner_airfoil = inner_wing_cross_section.airfoil.add_control_surface(
            deflection=inner_wing_cross_section.control_surface_deflection,
            hinge_point=inner_wing_cross_section.control_surface_hinge_point,
        )
        outer_airfoil = outer_wing_cross_section.airfoil.add_control_surface(
            deflection=outer_wing_cross_section.control_surface_deflection,
            hinge_point=outer_wing_cross_section.control_surface_hinge_point,
        )

        # Get the MCL point components for the inner and outer Airfoils.
        [
            inner_mcl_pointsY_Ai_lpAi,
            inner_mcl_pointsX_Ai_lpAi,
            outer_mcl_pointsY_Ao_lpAo,
            outer_mcl_pointsX_Ao_lpAo,
        ] = _get_mcl_points(inner_airfoil, outer_airfoil, chordwise_coordinates)

        # Define number of spanwise points and Panels. This is based on the inner
        # WingCrossSection.
        num_spanwise_panels = inner_wing_cross_section.num_spanwise_panels
        assert num_spanwise_panels is not None
        num_spanwise_coordinates = num_spanwise_panels + 1

        # Get the spanwise coordinates.
        if inner_wing_cross_section.spanwise_spacing == "uniform":
            spanwise_coordinates = np.linspace(0, 1, num_spanwise_coordinates)
        else:
            spanwise_coordinates = _functions.cosspace(0, 1, num_spanwise_coordinates)

        # Get this wing section's MCS points (in wing axes, relative to the leading
        # edge root point).
        [
            Fipp_Wn_Ler,
            Fopp_Wn_Ler,
            Bipp_Wn_Ler,
            Bopp_Wn_Ler,
        ] = _get_mcs_points(
            T_pas_Wcsi_Lpi_to_Wn_Ler,
            T_pas_Wcso_Lpo_to_Wn_Ler,
            inner_wing_cross_section,
            outer_wing_cross_section,
            inner_mcl_pointsY_Ai_lpAi,
            inner_mcl_pointsX_Ai_lpAi,
            outer_mcl_pointsY_Ao_lpAo,
            outer_mcl_pointsX_Ao_lpAo,
            spanwise_coordinates,
        )

        # Find the MCS points expressed in geometry axes, relative to the CG.
        Fipp_G_Cg = _transformations.apply_T_to_vectors(
            T_pas_Wn_Ler_to_G_Cg, Fipp_Wn_Ler, has_point=True
        )
        Fopp_G_Cg = _transformations.apply_T_to_vectors(
            T_pas_Wn_Ler_to_G_Cg, Fopp_Wn_Ler, has_point=True
        )
        Bipp_G_Cg = _transformations.apply_T_to_vectors(
            T_pas_Wn_Ler_to_G_Cg, Bipp_Wn_Ler, has_point=True
        )
        Bopp_G_Cg = _transformations.apply_T_to_vectors(
            T_pas_Wn_Ler_to_G_Cg, Bopp_Wn_Ler, has_point=True
        )

        # Compute a matrix that is (M,N), where M and N are the number of chordwise
        # and spanwise Panels. The values are either 1 if the Panel at that location
        # is a trailing edge, or 0 if not.
        wing_section_is_trailing_edge = np.vstack(
            [
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
                np.ones((1, num_spanwise_panels), dtype=bool),
            ]
        )

        # Compute a matrix that is (M,N), where M and N are the number of chordwise
        # and spanwise Panels. The values are either 1 if the Panel at that location
        # is a leading edge, or 0 if not.
        wing_section_is_leading_edge = np.vstack(
            [
                np.ones((1, num_spanwise_panels), dtype=bool),
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
            ]
        )

        Flpp_G_Cg = Fipp_G_Cg
        Frpp_G_Cg = Fopp_G_Cg
        Blpp_G_Cg = Bipp_G_Cg
        Brpp_G_Cg = Bopp_G_Cg

        if symmetry_type in (2, 3):
            Flpp_G_Cg = Fopp_G_Cg
            Frpp_G_Cg = Fipp_G_Cg
            Blpp_G_Cg = Bopp_G_Cg
            Brpp_G_Cg = Bipp_G_Cg

        # Get this wing section's Panels.
        wing_section_panels = _get_panels(
            Flpp_G_Cg=Flpp_G_Cg,
            Frpp_G_Cg=Frpp_G_Cg,
            Blpp_G_Cg=Blpp_G_Cg,
            Brpp_G_Cg=Brpp_G_Cg,
            is_trailing_edge=wing_section_is_trailing_edge,
            is_leading_edge=wing_section_is_leading_edge,
        )

        if symmetry_type in (2, 3):
            wing_panels = np.hstack([np.flip(wing_section_panels, axis=1), wing_panels])
        else:
            wing_panels = np.hstack([wing_panels, wing_section_panels])

        # Handle type 4 symmetry.
        if symmetry_type == 4:

            # Handle control surface symmetry.
            if inner_wing_cross_section.control_surface_symmetry_type == "symmetric":
                # Define the Airfoil at the inner WingCrossSection with flap-like
                # control surface symmetry.
                inner_airfoil = inner_wing_cross_section.airfoil.add_control_surface(
                    deflection=inner_wing_cross_section.control_surface_deflection,
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )
            else:  # If not "symmetric", this must be "asymmetric" for type 4 symmetry.
                # Define the Airfoil at the inner WingCrossSection with aileron-like
                # control surface symmetry.
                inner_airfoil = inner_wing_cross_section.airfoil.add_control_surface(
                    deflection=-inner_wing_cross_section.control_surface_deflection,
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )
            if outer_wing_cross_section.control_surface_symmetry_type == "symmetric":
                # Define the Airfoil at the outer WingCrossSection with flap-like
                # control surface symmetry.
                outer_airfoil = outer_wing_cross_section.airfoil.add_control_surface(
                    deflection=outer_wing_cross_section.control_surface_deflection,
                    hinge_point=outer_wing_cross_section.control_surface_hinge_point,
                )
            else:  # If not "symmetric", this must be "asymmetric" for type 4 symmetry.
                # Define the Airfoil at the outer WingCrossSection with aileron-like
                # control surface symmetry.
                outer_airfoil = outer_wing_cross_section.airfoil.add_control_surface(
                    deflection=-outer_wing_cross_section.control_surface_deflection,
                    hinge_point=outer_wing_cross_section.control_surface_hinge_point,
                )

            # Get the MCL point components for the inner and outer Airfoils.
            [
                inner_mcl_pointsY_Ai_lpAi,
                inner_mcl_pointsX_Ai_lpAi,
                outer_mcl_pointsY_Ao_lpAo,
                outer_mcl_pointsX_Ao_lpAo,
            ] = _get_mcl_points(inner_airfoil, outer_airfoil, chordwise_coordinates)

            # Get this wing section's preliminary panel points.
            [
                Fipp_Wn_Ler,
                Fopp_Wn_Ler,
                Bipp_Wn_Ler,
                Bopp_Wn_Ler,
            ] = _get_mcs_points(
                T_pas_Wcsi_Lpi_to_Wn_Ler,
                T_pas_Wcso_Lpo_to_Wn_Ler,
                inner_wing_cross_section,
                outer_wing_cross_section,
                inner_mcl_pointsY_Ai_lpAi,
                inner_mcl_pointsX_Ai_lpAi,
                outer_mcl_pointsY_Ao_lpAo,
                outer_mcl_pointsX_Ao_lpAo,
                spanwise_coordinates,
            )

            # Compute a matrix that is (M,N), where M and N are the number of
            # chordwise and spanwise Panels. The values are either 1 if the Panel at
            # that location is a trailing edge, or 0 if not.
            wing_section_is_trailing_edge = np.vstack(
                [
                    np.zeros(
                        (num_chordwise_panels - 1, num_spanwise_panels), dtype=bool
                    ),
                    np.ones((1, num_spanwise_panels), dtype=bool),
                ]
            )

            # Compute a matrix that is (M,N), where M and N are the number of
            # chordwise and spanwise Panels. The values are either 1 if the Panel at
            # that location is a leading edge, or 0 if not.
            wing_section_is_leading_edge = np.vstack(
                [
                    np.ones((1, num_spanwise_panels), dtype=bool),
                    np.zeros(
                        (num_chordwise_panels - 1, num_spanwise_panels), dtype=bool
                    ),
                ]
            )

            # TODO: Understand how this block of code works. We're reflecting about a
            #  plane defined by a normal vector in geometry axes and a point in
            #  geometry axes relative to the CG. But that's okay I guess and we end
            #  up with reflected points in wing axes relative to the leading edge
            #  root point? Also then we perform an passive transformation to find
            #  them in geometry axes relative to the CG? Why do we even need to
            #  reflect them if we are staying in wing axes? If we take a reflected a
            #  non-reflected Wing that are otherwise identical, they should have the
            #  same coordinates in their respective wing axes relative to their
            #  respective leading edge root points. So why is the active
            #  transformation necessary?
            # DOCUMENT: Document the logic in this block of code.
            reflect_T_act = _transformations.generate_reflect_T(
                plane_point_A_a=symmetryPoint_G_Cg,
                plane_normal_A=symmetryNormal_G,
                passive=False,
            )
            reflected_Fipp_Wn_Ler = _transformations.apply_T_to_vectors(
                reflect_T_act, Fipp_Wn_Ler, has_point=True
            )
            reflected_Fopp_Wn_Ler = _transformations.apply_T_to_vectors(
                reflect_T_act, Fopp_Wn_Ler, has_point=True
            )
            reflected_Bipp_Wn_Ler = _transformations.apply_T_to_vectors(
                reflect_T_act, Bipp_Wn_Ler, has_point=True
            )
            reflected_Bopp_Wn_Ler = _transformations.apply_T_to_vectors(
                reflect_T_act, Bopp_Wn_Ler, has_point=True
            )
            reflected_Fipp_G_Cg = _transformations.apply_T_to_vectors(
                T_pas_Wn_Ler_to_G_Cg, reflected_Fipp_Wn_Ler, has_point=True
            )
            reflected_Fopp_G_Cg = _transformations.apply_T_to_vectors(
                T_pas_Wn_Ler_to_G_Cg, reflected_Fopp_Wn_Ler, has_point=True
            )
            reflected_Bipp_G_Cg = _transformations.apply_T_to_vectors(
                T_pas_Wn_Ler_to_G_Cg, reflected_Bipp_Wn_Ler, has_point=True
            )
            reflected_Bopp_G_Cg = _transformations.apply_T_to_vectors(
                T_pas_Wn_Ler_to_G_Cg, reflected_Bopp_Wn_Ler, has_point=True
            )

            # Get the reflected wing section's Panels.
            wing_section_panels = _get_panels(
                Flpp_G_Cg=reflected_Fopp_G_Cg,
                Frpp_G_Cg=reflected_Fipp_G_Cg,
                Blpp_G_Cg=reflected_Bopp_G_Cg,
                Brpp_G_Cg=reflected_Bipp_G_Cg,
                is_trailing_edge=wing_section_is_trailing_edge,
                is_leading_edge=wing_section_is_leading_edge,
            )

            # This wing section's Panels array is stacked horizontally, to the left
            # of the Wing's Panel array.
            wing_panels = np.hstack([np.flip(wing_section_panels, axis=1), wing_panels])

    # Iterate through the Panels and populate their local position attributes.
    for chordwise_position in range(wing_panels.shape[0]):
        for spanwise_position in range(wing_panels.shape[1]):
            this_panel = wing_panels[chordwise_position, spanwise_position]
            this_panel.local_chordwise_position = chordwise_position
            this_panel.local_spanwise_position = spanwise_position
            if spanwise_position == 0:
                this_panel.is_left_edge = True
            else:
                this_panel.is_left_edge = False
            if spanwise_position == wing_panels.shape[1] - 1:
                this_panel.is_right_edge = True
            else:
                this_panel.is_right_edge = False

    # Populate the Wing's Panels attribute.
    wing.panels = wing_panels


def _get_mcl_points(
    inner_airfoil: airfoil_mod.Airfoil,
    outer_airfoil: airfoil_mod.Airfoil,
    chordwise_coordinates: np.ndarray,
) -> list[np.ndarray]:
    """Takes in the inner and outer Airfoils of a wing section and its normalized
    chordwise coordinates. It returns a list of four column vectors containing the
    normalized components of the positions of points along the mean camber line (MCL)
    (in each Airfoil's axes, relative to each Airfoil's leading point).

    :param inner_airfoil: The wing section's inner Airfoil.
    :param outer_airfoil: The wing section's outer Airfoil.
    :param chordwise_coordinates: A (N,) ndarray of floats for the normalized chordwise
        coordinates where we'd like to sample each Airfoil's MCL. The values are
        normalized from 0.0 to 1.0 and are unitless.
    :return: A list of four (N,1) ndarrays of floats, where N is the number of points at
        which we'd like to sample each Airfoil's MCL. The ndarrays contain components of
        the positions of points along each Airfoil's MCL. In order, the ndarrays
        returned are, (1) the inner Airfoil's MCL points' y-components, (2) the inner
        Airfoil's MCL points' x-components (3) the outer Airfoil's MCL points'
        y-components, and (4) the outer Airfoil's MCL points' x-components. The values
        are normalized from 0.0 to 1.0 and are unitless.
    """

    # Make the MCLs for each Airfoil. First index is point number, second index is
    # the coordinates of that point on the MCL (in each Airfoil's axes, relative to
    # each Airfoil's leading point).
    inner_mcl_points_Ai_lpAi = inner_airfoil.get_resampled_mcl(chordwise_coordinates)
    outer_mcl_points_Ao_lpAo = outer_airfoil.get_resampled_mcl(chordwise_coordinates)

    # Extract the y-components of the inner Airfoil's MCL points (in the inner
    # Airfoil's axes, relative to the inner Airfoil's leading point) and put them in
    # a column vector.
    inner_mcl_pointsY_Ai_lpAi = np.expand_dims(inner_mcl_points_Ai_lpAi[:, 1], 1)

    # Extract the x-components of the inner Airfoil's MCL points (in the inner
    # Airfoil's axes, relative to the inner Airfoil's leading point) and put them in
    # a column vector.
    inner_mcl_pointsX_Ai_lpAi = np.expand_dims(inner_mcl_points_Ai_lpAi[:, 0], 1)

    # Extract the y-components of the outer Airfoil's MCL points (in the outer
    # Airfoil's axes, relative to the outer Airfoil's leading point) and put them in
    # a column vector.
    outer_mcl_pointsY_Ao_lpAo = np.expand_dims(outer_mcl_points_Ao_lpAo[:, 1], 1)

    # Extract the x-components of the outer Airfoil's MCL points (in the outer
    # Airfoil's axes, relative to the outer Airfoil's leading point) and put them in
    # a column vector.
    outer_mcl_pointsX_Ao_lpAo = np.expand_dims(outer_mcl_points_Ao_lpAo[:, 0], 1)

    return [
        inner_mcl_pointsY_Ai_lpAi,
        inner_mcl_pointsX_Ai_lpAi,
        outer_mcl_pointsY_Ao_lpAo,
        outer_mcl_pointsX_Ao_lpAo,
    ]


def _get_mcs_points(
    T_pas_Wcsi_Lpi_Wn_Ler: np.ndarray,
    T_pas_Wcso_Lpo_Wn_Ler: np.ndarray,
    inner_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    outer_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    inner_mcl_pointsY_Ai_lpAi: np.ndarray,
    inner_mcl_pointsX_Ai_lpAi: np.ndarray,
    outer_mcl_pointsY_Ao_lpAo: np.ndarray,
    outer_mcl_pointsX_Ao_lpAo: np.ndarray,
    spanwise_coordinates: np.ndarray,
) -> list[np.ndarray]:
    """Calculates the points on a wing section's mean camber surface (MCS) (in wing
    axes, relative to the leading edge root point).

    :param T_pas_Wcsi_Lpi_Wn_Ler: A (4,4) ndarray of floats representing a passive
        transformation matrix which maps in homogeneous coordinates from the inner
        WingCrossSection's axes, relative to its leading point to wing axes relative to
        the leading edge root point.
    :param T_pas_Wcso_Lpo_Wn_Ler: A (4,4) ndarray of floats representing a passive
        transformation matrix which maps in homogeneous coordinates from the outer
        WingCrossSection's axes, relative to its leading point to wing axes relative to
        the leading edge root point.
    :param inner_wing_cross_section: The wing section's inner WingCrossSection.
    :param outer_wing_cross_section: The wing section's outer WingCrossSection.
    :param inner_mcl_pointsY_Ai_lpAi: A (M,1) ndarray of floats, where M is the number
        of chordwise points in the mesh. Each element represents the y-component of the
        inner Airfoil's MCL points (in the inner Airfoil's axes, relative to the inner
        Airfoil's leading point). The values are normalized from 0.0 to 1.0 and are
        unitless.
    :param inner_mcl_pointsX_Ai_lpAi: A (M,1) ndarray of floats, where M is the number
        of chordwise points in the mesh. Each element represents the x-component of the
        inner Airfoil's MCL points (in the inner Airfoil's axes, relative to the inner
        Airfoil's leading point). The values are normalized from 0.0 to 1.0 and are
        unitless.
    :param outer_mcl_pointsY_Ao_lpAo: A (M,1) ndarray of floats, where M is the number
        of chordwise points in the mesh. Each element represents the y-component of the
        outer Airfoil's MCL points (in the outer Airfoil's axes, relative to the outer
        Airfoil's leading point). The values are normalized from 0.0 to 1.0 and are
        unitless.
    :param outer_mcl_pointsX_Ao_lpAo: A (M,1) ndarray of floats, where M is the number
        of chordwise points in the mesh. Each element represents the x-component of the
        outer Airfoil's MCL points (in the outer Airfoil's axes, relative to the outer
        Airfoil's leading point). The values are normalized from 0.0 to 1.0 and are
        unitless.
    :param spanwise_coordinates: A (N,1) ndarray of floats, where N is the number of
        spanwise points. It holds the distances of each spanwise point along the wing
        section. The values are normalized from 0.0 to 1.0 and are unitless.
    :return: A list of four (M,N,3) ndarrays of floats, where M is the number of
        chordwise points and N is the number of spanwise points. The four ndarrays are,
        in order, this wing section's Panel's (1) forward inner, (2) forward outer, (3)
        backward inner, and (4) backward outer panel points (in wing axes, relative to
        the leading edge root point). The units are in meters.
    """
    inner_mcl_pointsX_Wcsi_Lpi = (
        inner_wing_cross_section.chord * inner_mcl_pointsX_Ai_lpAi
    )

    inner_mcl_pointsZ_Wcsi_Lpi = (
        inner_wing_cross_section.chord * inner_mcl_pointsY_Ai_lpAi
    )

    outer_mcl_pointsX_Wcso_Lpo = (
        outer_wing_cross_section.chord * outer_mcl_pointsX_Ao_lpAo
    )

    outer_mcl_pointsZ_Wcso_Lpo = (
        outer_wing_cross_section.chord * outer_mcl_pointsY_Ao_lpAo
    )

    inner_mcl_points_Wcsi_Lpi = np.hstack(
        [
            inner_mcl_pointsX_Wcsi_Lpi,
            np.zeros_like(inner_mcl_pointsX_Wcsi_Lpi),
            inner_mcl_pointsZ_Wcsi_Lpi,
        ]
    )
    outer_mcl_points_Wcso_Lpo = np.hstack(
        [
            outer_mcl_pointsX_Wcso_Lpo,
            np.zeros_like(outer_mcl_pointsX_Wcso_Lpo),
            outer_mcl_pointsZ_Wcso_Lpo,
        ]
    )

    inner_mcl_points_Wn_Ler = _transformations.apply_T_to_vectors(
        T_pas_Wcsi_Lpi_Wn_Ler, inner_mcl_points_Wcsi_Lpi, has_point=True
    )
    outer_mcl_points_Wn_Ler = _transformations.apply_T_to_vectors(
        T_pas_Wcso_Lpo_Wn_Ler, outer_mcl_points_Wcso_Lpo, has_point=True
    )

    # Find the vertices of the points on this wing section's (MCS) (in wing axes,
    # relative to the leading edge root point) with interpolation. This returns an (
    # M,N,3) array, where M and N are the number of chordwise points and spanwise
    # points.
    wing_section_mcl_vertices = _functions.interp_between_points(
        inner_mcl_points_Wn_Ler, outer_mcl_points_Wn_Ler, spanwise_coordinates
    )

    # Extract the coordinates for corners of each panel point.
    Fipp_Wn_Ler = wing_section_mcl_vertices[:-1, :-1, :]
    Fopp_Wn_Ler = wing_section_mcl_vertices[:-1, 1:, :]
    Bipp_Wn_Ler = wing_section_mcl_vertices[1:, :-1, :]
    Bopp_Wn_Ler = wing_section_mcl_vertices[1:, 1:, :]

    return [
        Fipp_Wn_Ler,
        Fopp_Wn_Ler,
        Bipp_Wn_Ler,
        Bopp_Wn_Ler,
    ]


def _get_panels(
    Flpp_G_Cg: np.ndarray,
    Frpp_G_Cg: np.ndarray,
    Blpp_G_Cg: np.ndarray,
    Brpp_G_Cg: np.ndarray,
    is_trailing_edge: np.ndarray,
    is_leading_edge: np.ndarray,
) -> np.ndarray:
    """Takes in arrays of Panel attributes and returns a 2D ndarray of Panels.

    :param Flpp_G_Cg: A (M,N,3) ndarray of floats, where M is the number of chordwise
        Panels, N is the number of spanwise Panels, and the last dimension contains the
        position vector of each Panel's front left vertex (in geometry axes, relative to
        the CG). The values are in meters.
    :param Frpp_G_Cg: A (M,N,3) ndarray of floats, where M is the number of chordwise
        Panels, N is the number of spanwise Panels, and the last dimension contains the
        position vector of each Panel's front right vertex (in geometry axes, relative
        to the CG). The values are in meters.
    :param Blpp_G_Cg: A (M,N,3) ndarray of floats, where M is the number of chordwise
        Panels, N is the number of spanwise Panels, and the last dimension contains the
        position vector of each Panel's back left vertex (in geometry axes, relative to
        the CG). The values are in meters.
    :param Brpp_G_Cg: A (M,N,3) ndarray of floats, where M is the number of chordwise
        Panels, N is the number of spanwise Panels, and the last dimension contains the
        position vector of each Panel's back right vertex (in geometry axes, relative to
        the CG). The values are in meters.
    :param is_trailing_edge: A (M,N) ndarray of bools that denote if the Panel in each
        location is on the trailing edge of the Wing.
    :param is_leading_edge: A (M,N) ndarray of bools that denote if the Panel in each
        location is on the leading edge of the Wing.
    :return panel_array: A (M,N) ndarray of Panels constructed using the given
        parameters.
    """
    num_chordwise_panels = Flpp_G_Cg.shape[0]
    num_spanwise_panels = Flpp_G_Cg.shape[1]

    # Initialize an empty ndarray to hold the wing section's Panels. The ndarray is
    # size (M,N), where M and N are the number of chordwise and spanwise Panels.
    panels = np.empty((num_chordwise_panels, num_spanwise_panels), dtype=object)

    # Loop through the empty Panels array and create a new Panel in each position.
    for chordwise_position in range(num_chordwise_panels):
        for spanwise_position in range(num_spanwise_panels):
            panels[chordwise_position, spanwise_position] = _panel.Panel(
                Frpp_G_Cg=Frpp_G_Cg[chordwise_position, spanwise_position],
                Flpp_G_Cg=Flpp_G_Cg[chordwise_position, spanwise_position],
                Blpp_G_Cg=Blpp_G_Cg[chordwise_position, spanwise_position],
                Brpp_G_Cg=Brpp_G_Cg[chordwise_position, spanwise_position],
                is_leading_edge=is_leading_edge[chordwise_position, spanwise_position],
                is_trailing_edge=is_trailing_edge[
                    chordwise_position, spanwise_position
                ],
            )

    return panels
