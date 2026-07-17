"""Contains functions for building and refining the meshes of convergence iterations."""

from __future__ import annotations

import copy
from collections.abc import Callable, Sequence

import numpy as np

from . import _logging, geometry, movements, problems

_logger = _logging.get_logger("_convergence_meshing")


def build_steady_problem(
    ar_id: int,
    chord_id: int,
    panel_aspect_ratio: int,
    num_chordwise_panels: int,
    ref_problem: problems.SteadyProblem,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
) -> problems.SteadyProblem:
    """Builds the SteadyProblem for one convergence iteration.

    Each of the reference SteadyProblem's Airplanes is copied with the given Panel
    aspect ratio and number of chordwise Panels. Each Wing is refined by its spanwise
    mesh: a trapezoidal Wing has each non-tip WingCrossSection's number of spanwise
    Panels resolved (and cached) to achieve the desired Panel aspect ratio, while an
    edge-defined Wing has its stored edge curves resampled into the number of
    WingCrossSections (also resolved and cached) that achieves it. The copied Airplanes
    are wrapped in a new SteadyProblem with the reference OperatingPoint.

    :param ar_id: The index of the current Panel aspect ratio within the list of Panel
        aspect ratios being tested.
    :param chord_id: The index of the current number of chordwise Panels within the list
        of numbers of chordwise Panels being tested.
    :param panel_aspect_ratio: The Panel aspect ratio to target for this iteration.
    :param num_chordwise_panels: The number of chordwise Panels to use for this
        iteration.
    :param ref_problem: The reference SteadyProblem whose Airplanes and OperatingPoint
        are copied.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels, used for
        trapezoidal Wings. It is read and updated in place.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections, used for edge-defined Wings.
        It is read and updated in place.
    :return: The SteadyProblem for this iteration.
    """
    ref_airplanes = ref_problem.airplanes
    ref_operating_point = ref_problem.operating_point

    # Initialize an empty list to hold this iteration's Airplanes. Then, fill the list
    # by making new copies of each of the Airplanes with modified values for Panel
    # aspect ratio and number of chordwise Panels.
    these_airplanes = []
    for ref_airplane_id, ref_airplane in enumerate(ref_airplanes):
        ref_wings = ref_airplane.wings
        these_wings = []

        for ref_wing_id, ref_wing in enumerate(ref_wings):

            # An edge-defined Wing is refined by resampling its stored edge curves into
            # the number of WingCrossSections that achieves the desired Panel aspect
            # ratio, rather than by sweeping the number of spanwise Panels.
            if ref_wing.spanwise_mesh == "edge_defined":
                # The resolution's log messages nest one level under the build's.
                with _logging.nested():
                    this_num_wing_cross_sections = _resolve_num_wing_cross_sections(
                        ar_id,
                        chord_id,
                        ref_airplane_id,
                        ref_wing_id,
                        ref_airplane.name,
                        ref_wing.name,
                        num_wing_cross_sections_cache,
                        lambda start: _get_num_wing_cross_sections_for_panel_ar(
                            panel_aspect_ratio,
                            num_chordwise_panels,
                            ref_wing,
                            start,
                        ),
                    )
                these_wings.append(
                    _build_edge_defined_wing(
                        ref_wing, num_chordwise_panels, this_num_wing_cross_sections
                    )
                )
                continue

            ref_wing_cross_sections = ref_wing.wing_cross_sections
            these_wing_cross_sections = []

            for ref_wing_cross_section_id, ref_wing_cross_section in enumerate(
                ref_wing_cross_sections
            ):

                # If this is not the last WingCrossSection, find the number of spanwise
                # Panels to use for this section of the Wing, based on the desired Panel
                # aspect ratio and number of chordwise Panels.
                if ref_wing_cross_section_id < (len(ref_wing_cross_sections) - 1):
                    # The resolution's log messages nest one level under the build's.
                    with _logging.nested():
                        this_num_spanwise_panels = _resolve_num_spanwise_panels(
                            ar_id,
                            chord_id,
                            ref_airplane_id,
                            ref_wing_id,
                            ref_wing_cross_section_id,
                            ref_airplane.name,
                            ref_wing.name,
                            num_spanwise_panels_cache,
                            lambda start: _get_wing_section_num_spanwise_panels(
                                panel_aspect_ratio,
                                num_chordwise_panels,
                                ref_wing.chordwise_spacing,
                                ref_wing_cross_section,
                                ref_wing_cross_sections[ref_wing_cross_section_id + 1],
                                start,
                            ),
                        )
                else:
                    this_num_spanwise_panels = None

                these_wing_cross_sections.append(
                    geometry.wing_cross_section.WingCrossSection(
                        # These values are copied from the reference WingCrossSection.
                        chord=ref_wing_cross_section.chord,
                        Lp_Wcsp_Lpp=ref_wing_cross_section.Lp_Wcsp_Lpp,
                        angles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                        control_surface_symmetry_type=ref_wing_cross_section.control_surface_symmetry_type,
                        control_surface_hinge_point=ref_wing_cross_section.control_surface_hinge_point,
                        control_surface_deflection=ref_wing_cross_section.control_surface_deflection,
                        spanwise_spacing=ref_wing_cross_section.spanwise_spacing,
                        # These values change.
                        num_spanwise_panels=this_num_spanwise_panels,
                        airfoil=geometry.airfoil.Airfoil(
                            name=ref_wing_cross_section.airfoil.name,
                            outline_A_lp=ref_wing_cross_section.airfoil.outline_A_lp,
                            resample=ref_wing_cross_section.airfoil.resample,
                            n_points_per_side=ref_wing_cross_section.airfoil.n_points_per_side,
                        ),
                    )
                )

            these_wings.append(
                geometry.wing.Wing(
                    # These values are copied from the reference Wing.
                    name=ref_wing.name,
                    Ler_Gs_Cgs=ref_wing.Ler_Gs_Cgs,
                    angles_Gs_to_Wn_ixyz=ref_wing.angles_Gs_to_Wn_ixyz,
                    symmetric=ref_wing.symmetric,
                    mirror_only=ref_wing.mirror_only,
                    symmetryNormal_G=ref_wing.symmetryNormal_G,
                    symmetryPoint_G_Cg=ref_wing.symmetryPoint_G_Cg,
                    chordwise_spacing=ref_wing.chordwise_spacing,
                    # These values change.
                    wing_cross_sections=these_wing_cross_sections,
                    num_chordwise_panels=num_chordwise_panels,
                )
            )

        these_airplanes.append(
            geometry.airplane.Airplane(
                # These values are copied from the reference Airplane.
                name=ref_airplane.name,
                Cg_GP1_CgP1=ref_airplane.Cg_GP1_CgP1,
                weight=ref_airplane.weight,
                # These values change.
                wings=these_wings,
                s_ref=None,
                c_ref=None,
                b_ref=None,
            )
        )

    # Create a new SteadyProblem for this iteration with its own deep-copied
    # OperatingPoint. Solving populates the OperatingPoint's lazy caches, so sharing the
    # reference OperatingPoint would mutate the reference problem and change its content
    # hash between an uncached run and a later cached one.
    return problems.SteadyProblem(
        airplanes=these_airplanes,
        operating_point=copy.deepcopy(ref_operating_point),
    )


def build_unsteady_problem(
    ar_id: int,
    chord_id: int,
    panel_aspect_ratio: int,
    num_chordwise_panels: int,
    wake_length: int,
    ref_problem: problems.UnsteadyProblem,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
    delta_time_cache: dict[tuple[int, int], float],
) -> problems.UnsteadyProblem:
    """Builds the UnsteadyProblem for one convergence iteration.

    Each of the reference Movement's AirplaneMovements (and its nested WingMovements and
    WingCrossSectionMovements) is copied with the given Panel aspect ratio and number of
    chordwise Panels, preserving every motion parameter. Each Wing is refined by its
    spanwise mesh: a trapezoidal Wing has each non-tip WingCrossSection's number of
    spanwise Panels resolved (and cached) to achieve the desired Panel aspect ratio,
    while an edge-defined Wing has its stored edge curves resampled into the number of
    WingCrossSections (also resolved and cached) that achieves it. An edge-defined
    Wing's WingCrossSectionMovements must all be static, because resampling changes the
    WingCrossSection count and so cannot preserve per-WingCrossSection motion; the
    refined WingCrossSections are wrapped in motion free WingCrossSectionMovements and
    the WingMovement's own motion parameters are copied. The copied AirplaneMovements
    are wrapped in a new Movement (with the given wake length) and UnsteadyProblem.

    The Movement's time step is optimized for this iteration's mesh with Movement's
    iterative delta_time optimizer. The optimum depends only on the mesh, so it is
    resolved once per (Panel aspect ratio index, number of chordwise Panels index) and
    cached, then reused across every wake-state and wake-length combination at that
    mesh.

    :param ar_id: The index of the current Panel aspect ratio within the list of Panel
        aspect ratios being tested.
    :param chord_id: The index of the current number of chordwise Panels within the list
        of numbers of chordwise Panels being tested.
    :param panel_aspect_ratio: The Panel aspect ratio to target for this iteration.
    :param num_chordwise_panels: The number of chordwise Panels to use for this
        iteration.
    :param wake_length: The wake length (number of chords if static, number of cycles if
        variable) to use for this iteration.
    :param ref_problem: The reference UnsteadyProblem whose Movement is copied.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels, used for
        trapezoidal Wings. It is read and updated in place.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections, used for edge-defined Wings.
        It is read and updated in place.
    :param delta_time_cache: The cache mapping a (Panel aspect ratio index, number of
        chordwise Panels index) tuple to a previously optimized delta_time for that
        mesh. It is read and updated in place.
    :return: The UnsteadyProblem for this iteration.
    """
    ref_movement = ref_problem.movement
    static = ref_movement.static
    ref_airplane_movements = ref_movement.airplane_movements
    ref_operating_point_movement = ref_movement.operating_point_movement

    # Initialize an empty list for this iteration's base AirplaneMovements.
    these_base_airplanes = []

    # Create an empty list for the AirplaneMovement copies.
    these_airplane_movements = []

    # Now, we will begin iterating through this iteration's reference AirplaneMovements,
    # WingMovements, and WingCrossSectionMovements, and creating copies of them. These
    # copies will have identical parameters to their respective reference movements
    # except for the number of spanwise Panels (which is based on the Panel aspect
    # ratio), and the number of chordwise Panels.
    #
    # To do this, we iterate over the AirplaneMovements and perform a several step
    # procedure:
    # 1. Reference the AirplaneMovement's base Airplane.
    # 2. Reference the AirplaneMovement's list of WingMovements.
    # 3. Create an empty list for the WingMovements' base Wing copies.
    # 4. Create an empty list for the WingMovement copies.
    # 5. Iterate over the WingMovements.
    #   5.1. Reference the WingMovement's base Wing.
    #   5.2. Reference the WingMovement's list of WingCrossSectionMovements.
    #   5.3. Create an empty list for the WingCrossSectionMovements' base
    #        WingCrossSection copies.
    #   5.4. Create an empty list for the WingCrossSectionMovement copies.
    #   5.5. Iterate over the WingCrossSectionMovements.
    #     5.5.1. Reference the WingCrossSectionMovement's base WingCrossSection.
    #     5.5.2. Calculate the number of spanwise Panels that corresponds to the desired
    #            combination of Panel aspect ratio and number of chordwise Panels.
    #     5.5.3. Create a copy of the base WingCrossSection.
    #     5.5.4. Create a copy of the WingCrossSectionMovement.
    #     5.5.5. Append the base WingCrossSection copy to the list of base
    #            WingCrossSection copies.
    #     5.5.6. Append the WingCrossSectionMovement copy to the list of
    #            WingCrossSectionMovement copies.
    #   5.6. Create a copy of the base Wing.
    #   5.7. Create a copy of the WingMovement.
    #   5.8. Append the base Wing copy to the list of base Wing copies.
    #   5.9. Append the WingMovement copy to the list of WingMovement copies.
    # 6. Create a copy of the base Airplane.
    # 7. Create a copy of the AirplaneMovement.
    # 8. Append the base Airplane copy to the list of base Airplane copies.
    # 9. Append the AirplaneMovement copy to the list of AirplaneMovement copies.
    ref_airplane_movement: movements.airplane_movement.AirplaneMovement
    for ref_airplane_movement_id, ref_airplane_movement in enumerate(
        ref_airplane_movements
    ):
        # 1. Reference the AirplaneMovement's base Airplane.
        ref_base_airplane = ref_airplane_movement.base_airplane

        # 2. Reference the AirplaneMovement's list of WingMovements.
        ref_wing_movements = ref_airplane_movement.wing_movements

        # 3. Create an empty list for the WingMovements' base Wing copies.
        these_base_wings = []

        # 4. Create an empty list for the WingMovement copies.
        these_wing_movements = []

        # 5. Iterate over the WingMovements.
        for ref_wing_movement_id, ref_wing_movement in enumerate(ref_wing_movements):
            # 5.1. Reference the WingMovement's base Wing.
            ref_base_wing = ref_wing_movement.base_wing

            # An edge-defined Wing is refined by resampling its stored edge curves into
            # the number of WingCrossSections that achieves the desired Panel aspect
            # ratio. Resampling changes the WingCrossSection count, so it cannot
            # preserve per-WingCrossSection motion. Every WingCrossSectionMovement must
            # therefore be static; the refined WingCrossSections are wrapped in motion
            # free WingCrossSectionMovements and only the WingMovement's own motion is
            # carried over.
            if ref_base_wing.spanwise_mesh == "edge_defined":
                for (
                    ref_wing_cross_section_movement_id,
                    ref_wing_cross_section_movement,
                ) in enumerate(ref_wing_movement.wing_cross_section_movements):
                    if np.any(ref_wing_cross_section_movement.ampLp_Wcsp_Lpp) or np.any(
                        ref_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz
                    ):
                        raise ValueError(
                            "analyze_unsteady_convergence cannot refine an edge-defined "
                            "Wing whose WingCrossSectionMovements are not all static, "
                            "because resampling the Wing changes its number of "
                            "WingCrossSections. The WingCrossSectionMovement at index "
                            f"{ref_wing_cross_section_movement_id} of the WingMovement "
                            f'for the Wing named "{ref_base_wing.name}" is not static.'
                        )

                # The resolution's log messages nest one level under the build's.
                with _logging.nested():
                    this_num_wing_cross_sections = _resolve_num_wing_cross_sections(
                        ar_id,
                        chord_id,
                        ref_airplane_movement_id,
                        ref_wing_movement_id,
                        ref_base_airplane.name,
                        ref_base_wing.name,
                        num_wing_cross_sections_cache,
                        lambda start: _get_num_wing_cross_sections_for_panel_ar(
                            panel_aspect_ratio,
                            num_chordwise_panels,
                            ref_base_wing,
                            start,
                        ),
                    )

                this_base_wing = _build_edge_defined_wing(
                    ref_base_wing, num_chordwise_panels, this_num_wing_cross_sections
                )

                these_wing_cross_section_movements = [
                    movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=this_base_wing_cross_section
                    )
                    for this_base_wing_cross_section in (
                        this_base_wing.wing_cross_sections
                    )
                ]

                this_wing_movement = movements.wing_movement.WingMovement(
                    # These values are copied from the reference WingMovement.
                    ampLer_Gs_Cgs=ref_wing_movement.ampLer_Gs_Cgs,
                    periodLer_Gs_Cgs=ref_wing_movement.periodLer_Gs_Cgs,
                    spacingLer_Gs_Cgs=ref_wing_movement.spacingLer_Gs_Cgs,
                    phaseLer_Gs_Cgs=ref_wing_movement.phaseLer_Gs_Cgs,
                    ampAngles_Gs_to_Wn_ixyz=ref_wing_movement.ampAngles_Gs_to_Wn_ixyz,
                    periodAngles_Gs_to_Wn_ixyz=(
                        ref_wing_movement.periodAngles_Gs_to_Wn_ixyz
                    ),
                    spacingAngles_Gs_to_Wn_ixyz=(
                        ref_wing_movement.spacingAngles_Gs_to_Wn_ixyz
                    ),
                    phaseAngles_Gs_to_Wn_ixyz=(
                        ref_wing_movement.phaseAngles_Gs_to_Wn_ixyz
                    ),
                    rotationPointOffset_Gs_Ler=(
                        ref_wing_movement.rotationPointOffset_Gs_Ler
                    ),
                    # These values change.
                    base_wing=this_base_wing,
                    wing_cross_section_movements=these_wing_cross_section_movements,
                )

                these_base_wings.append(this_base_wing)
                these_wing_movements.append(this_wing_movement)
                continue

            # 5.2. Reference the WingMovement's list of WingCrossSectionMovements.
            ref_wing_cross_section_movements = (
                ref_wing_movement.wing_cross_section_movements
            )

            # 5.3. Create an empty list for the WingCrossSectionMovements' base
            # WingCrossSection copies.
            these_base_wing_cross_sections = []

            # 5.4. Create an empty list for the WingCrossSectionMovement copies.
            these_wing_cross_section_movements = []

            # 5.5. Iterate over the WingCrossSectionMovements.
            for (
                ref_wing_cross_section_movement_id,
                ref_wing_cross_section_movement,
            ) in enumerate(ref_wing_cross_section_movements):
                # 5.5.1. Reference the WingCrossSectionMovement's base WingCrossSection.
                ref_base_wing_cross_section = (
                    ref_wing_cross_section_movement.base_wing_cross_section
                )

                # 5.5.2. Calculate the number of spanwise Panels that corresponds to the
                # desired combination of Panel aspect ratio and number of chordwise
                # Panels.
                if ref_wing_cross_section_movement_id < (
                    len(ref_wing_cross_section_movements) - 1
                ):
                    # The resolution's log messages nest one level under the build's.
                    with _logging.nested():
                        this_num_spanwise_panels = _resolve_num_spanwise_panels(
                            ar_id,
                            chord_id,
                            ref_airplane_movement_id,
                            ref_wing_movement_id,
                            ref_wing_cross_section_movement_id,
                            ref_base_airplane.name,
                            ref_base_wing.name,
                            num_spanwise_panels_cache,
                            lambda start: _get_wing_section_movement_num_spanwise_panels(
                                panel_aspect_ratio,
                                num_chordwise_panels,
                                ref_base_wing.chordwise_spacing,
                                ref_movement.airplanes[ref_airplane_movement_id],
                                ref_wing_movement_id,
                                ref_wing_cross_section_movement_id,
                                ref_wing_cross_section_movement_id + 1,
                                start,
                                ref_problem.first_averaging_step,
                            ),
                        )
                else:
                    this_num_spanwise_panels = None

                # 5.5.3. Create a copy of the base WingCrossSection.
                this_base_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                    # These values are copied from the reference base WingCrossSection.
                    chord=ref_base_wing_cross_section.chord,
                    Lp_Wcsp_Lpp=ref_base_wing_cross_section.Lp_Wcsp_Lpp,
                    angles_Wcsp_to_Wcs_ixyz=ref_base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                    control_surface_symmetry_type=ref_base_wing_cross_section.control_surface_symmetry_type,
                    control_surface_hinge_point=ref_base_wing_cross_section.control_surface_hinge_point,
                    control_surface_deflection=ref_base_wing_cross_section.control_surface_deflection,
                    spanwise_spacing=ref_base_wing_cross_section.spanwise_spacing,
                    # These values change.
                    num_spanwise_panels=this_num_spanwise_panels,
                    airfoil=geometry.airfoil.Airfoil(
                        name=ref_base_wing_cross_section.airfoil.name,
                        outline_A_lp=ref_base_wing_cross_section.airfoil.outline_A_lp,
                        resample=ref_base_wing_cross_section.airfoil.resample,
                        n_points_per_side=ref_base_wing_cross_section.airfoil.n_points_per_side,
                    ),
                )

                # 5.5.4. Create a copy of the WingCrossSectionMovement.
                this_wing_cross_section_movement = movements.wing_cross_section_movement.WingCrossSectionMovement(
                    # These values are copied from the reference
                    # WingCrossSectionMovement.
                    ampLp_Wcsp_Lpp=ref_wing_cross_section_movement.ampLp_Wcsp_Lpp,
                    periodLp_Wcsp_Lpp=ref_wing_cross_section_movement.periodLp_Wcsp_Lpp,
                    spacingLp_Wcsp_Lpp=ref_wing_cross_section_movement.spacingLp_Wcsp_Lpp,
                    phaseLp_Wcsp_Lpp=ref_wing_cross_section_movement.phaseLp_Wcsp_Lpp,
                    ampAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz,
                    periodAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.periodAngles_Wcsp_to_Wcs_ixyz,
                    spacingAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.spacingAngles_Wcsp_to_Wcs_ixyz,
                    phaseAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.phaseAngles_Wcsp_to_Wcs_ixyz,
                    # These values change.
                    base_wing_cross_section=this_base_wing_cross_section,
                )

                # 5.5.5. Append the base WingCrossSection copy to the list of base
                # WingCrossSection copies.
                these_base_wing_cross_sections.append(this_base_wing_cross_section)

                # 5.5.6. Append the WingCrossSectionMovement copy to the list of
                # WingCrossSectionMovement copies.
                these_wing_cross_section_movements.append(
                    this_wing_cross_section_movement
                )

            # 5.6. Create a copy of base Wing.
            this_base_wing = geometry.wing.Wing(
                # These values are copied from the reference Wing.
                name=ref_base_wing.name,
                Ler_Gs_Cgs=ref_base_wing.Ler_Gs_Cgs,
                angles_Gs_to_Wn_ixyz=ref_base_wing.angles_Gs_to_Wn_ixyz,
                symmetric=ref_base_wing.symmetric,
                mirror_only=ref_base_wing.mirror_only,
                symmetryNormal_G=ref_base_wing.symmetryNormal_G,
                symmetryPoint_G_Cg=ref_base_wing.symmetryPoint_G_Cg,
                chordwise_spacing=ref_base_wing.chordwise_spacing,
                # These values change.
                wing_cross_sections=these_base_wing_cross_sections,
                num_chordwise_panels=num_chordwise_panels,
            )

            # 5.7. Create a copy of the WingMovement.
            this_wing_movement = movements.wing_movement.WingMovement(
                # These values are copied from the reference WingMovement.
                ampLer_Gs_Cgs=ref_wing_movement.ampLer_Gs_Cgs,
                periodLer_Gs_Cgs=ref_wing_movement.periodLer_Gs_Cgs,
                spacingLer_Gs_Cgs=ref_wing_movement.spacingLer_Gs_Cgs,
                phaseLer_Gs_Cgs=ref_wing_movement.phaseLer_Gs_Cgs,
                ampAngles_Gs_to_Wn_ixyz=ref_wing_movement.ampAngles_Gs_to_Wn_ixyz,
                periodAngles_Gs_to_Wn_ixyz=ref_wing_movement.periodAngles_Gs_to_Wn_ixyz,
                spacingAngles_Gs_to_Wn_ixyz=ref_wing_movement.spacingAngles_Gs_to_Wn_ixyz,
                phaseAngles_Gs_to_Wn_ixyz=ref_wing_movement.phaseAngles_Gs_to_Wn_ixyz,
                rotationPointOffset_Gs_Ler=ref_wing_movement.rotationPointOffset_Gs_Ler,
                # These values change.
                base_wing=this_base_wing,
                wing_cross_section_movements=these_wing_cross_section_movements,
            )

            # 5.8. Append the base Wing copy to the list of base Wing copies.
            these_base_wings.append(this_base_wing)

            # 5.9. Append the WingMovement copy to the list of WingMovement copies.
            these_wing_movements.append(this_wing_movement)

        # 6. Create a copy of the base Airplane.
        this_base_airplane = geometry.airplane.Airplane(
            # These values are copied from the reference Airplane.
            name=ref_base_airplane.name,
            Cg_GP1_CgP1=ref_base_airplane.Cg_GP1_CgP1,
            weight=ref_base_airplane.weight,
            # These values change.
            wings=these_base_wings,
            s_ref=None,
            c_ref=None,
            b_ref=None,
        )

        # 7. Create a copy of the AirplaneMovement.
        this_airplane_movement = movements.airplane_movement.AirplaneMovement(
            # These values are copied from the reference AirplaneMovement.
            ampCg_GP1_CgP1=ref_airplane_movement.ampCg_GP1_CgP1,
            periodCg_GP1_CgP1=ref_airplane_movement.periodCg_GP1_CgP1,
            spacingCg_GP1_CgP1=ref_airplane_movement.spacingCg_GP1_CgP1,
            phaseCg_GP1_CgP1=ref_airplane_movement.phaseCg_GP1_CgP1,
            # These values change.
            base_airplane=this_base_airplane,
            wing_movements=these_wing_movements,
        )

        # 8. Append the base Airplane copy to the list of base Airplane copies.
        these_base_airplanes.append(this_base_airplane)

        # 9. Append the AirplaneMovement copy to the list of AirplaneMovement copies.
        these_airplane_movements.append(this_airplane_movement)

    # Resolve the time step for this mesh. The optimum depends only on the mesh, so on a
    # cache miss the iterative optimizer finds it (delta_time="optimize"), and on a hit
    # the cached value is reused. The Movement resolves the string to a float, which is
    # read back and cached below.
    delta_time_key = (ar_id, chord_id)
    this_delta_time: str | float = delta_time_cache.get(delta_time_key, "optimize")

    # Create a new Movement for this iteration.
    if static:
        this_movement = movements.movement.Movement(
            airplane_movements=these_airplane_movements,
            operating_point_movement=ref_operating_point_movement,
            num_chords=wake_length,
            delta_time=this_delta_time,
        )
    else:
        this_movement = movements.movement.Movement(
            airplane_movements=these_airplane_movements,
            operating_point_movement=ref_operating_point_movement,
            num_cycles=wake_length,
            delta_time=this_delta_time,
        )

    # Cache the optimized delta_time on a miss so it is reused across the other
    # wake-state and wake-length combinations at this mesh. Only a hit is logged, since
    # on a miss the optimizer has just logged the value itself.
    if delta_time_key not in delta_time_cache:
        delta_time_cache[delta_time_key] = this_movement.delta_time
    else:
        _logger.info(
            _logging.indent()
            + "Cached delta_time: "
            + f"{this_movement.delta_time:#.3G}"
            + " s"
        )

    # Create a new UnsteadyProblem for this iteration.
    return problems.UnsteadyProblem(
        movement=this_movement,
        only_final_results=True,
    )


def memos_complete(
    ar_id: int,
    chord_id: int,
    ref_airplanes: Sequence[geometry.airplane.Airplane],
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
    delta_time_cache: dict[tuple[int, int], float] | None = None,
) -> bool:
    """Reports whether every memo that building one convergence iteration's problem
    would record is already cached.

    Building an iteration's problem is what resolves and records its mesh memos: the
    number of spanwise Panels for each trapezoidal Wing's non-tip WingCrossSections, the
    number of WingCrossSections for each edge-defined Wing, and, for an unsteady
    analysis, the mesh's optimized delta_time. When an iteration's solve is served from
    the solve cache and this function returns True, nothing needs the iteration's
    problem, so the caller can skip building it. Only the reference geometry's structure
    is walked, so this check never meshes anything.

    :param ar_id: The index of the current Panel aspect ratio within the list of Panel
        aspect ratios being tested.
    :param chord_id: The index of the current number of chordwise Panels within the list
        of numbers of chordwise Panels being tested.
    :param ref_airplanes: The reference Airplanes whose structure determines which memo
        keys the build would record. For an unsteady analysis, these are the reference
        AirplaneMovements' base Airplanes.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels, used for
        trapezoidal Wings.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections, used for edge-defined Wings.
    :param delta_time_cache: The cache mapping a (Panel aspect ratio index, number of
        chordwise Panels index) tuple to a previously optimized delta_time, or None for
        a steady analysis, which has no delta_time. The default is None.
    :return: True if every memo the build would record is already cached, otherwise
        False.
    """
    if delta_time_cache is not None and (ar_id, chord_id) not in delta_time_cache:
        return False

    for airplane_id, ref_airplane in enumerate(ref_airplanes):
        for wing_id, ref_wing in enumerate(ref_airplane.wings):

            # An edge-defined Wing records one memo: its number of WingCrossSections.
            if ref_wing.spanwise_mesh == "edge_defined":
                if (
                    ar_id,
                    chord_id,
                    airplane_id,
                    wing_id,
                ) not in num_wing_cross_sections_cache:
                    return False
                continue

            # A trapezoidal Wing records one memo per non-tip WingCrossSection: its
            # number of spanwise Panels.
            for wing_cross_section_id in range(len(ref_wing.wing_cross_sections) - 1):
                if (
                    ar_id,
                    chord_id,
                    airplane_id,
                    wing_id,
                    wing_cross_section_id,
                ) not in num_spanwise_panels_cache:
                    return False

    return True


def _get_wing_section_movement_num_spanwise_panels(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_airplanes: tuple[geometry.airplane.Airplane, ...],
    ref_wing_id: int,
    ref_root_wing_cross_section_id: int,
    ref_tip_wing_cross_section_id: int,
    start_val: int,
    first_applicable_time_step_id: int,
) -> int:
    """Calculates the number of spanwise Panels to use for the wing section of a
    WingMovement based on a desired average Panel aspect ratio.

    :param desired_average_panel_aspect_ratio: The target average Panel aspect ratio to
        achieve. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes). It must be a positive int.
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param chordwise_spacing: The type of spacing between the chordwise Panels. Can be
        "cosine" or "uniform".
    :param ref_airplanes: A tuple of the Airplanes at each time step.
    :param ref_wing_id: The index of the Wing within each Airplane in ref_airplanes. It
        must be a non negative int.
    :param ref_root_wing_cross_section_id: The index of the root WingCrossSection of the
        wing section within the Wing. It must be a non negative int.
    :param ref_tip_wing_cross_section_id: The index of the tip WingCrossSection of the
        wing section within the Wing. It must be a non negative int and greater than
        ``ref_root_wing_cross_section_id``.
    :param start_val: The initial number of spanwise Panels to start the search from. It
        must be a positive int. Using a higher value can speed up the search if a lower
        bound is already known.
    :param first_applicable_time_step_id: The index within ref_airplanes of the first
        time step to consider. It must be a non negative int. All Airplanes from this
        index onward will be analyzed.
    :return: The maximum number of spanwise Panels needed across all applicable time
        steps to achieve the desired average Panel aspect ratio.
    """
    # Slice the tuple of Airplanes to only the applicable ones. For cases with static
    # geometry, this is just the last time step's Airplane. For cases with variable
    # geometry, this is the last max-period cycle's time steps' Airplanes
    ref_airplanes = ref_airplanes[first_applicable_time_step_id:]

    num_time_steps = len(ref_airplanes)
    these_num_spanwise_panels = np.zeros_like(ref_airplanes, dtype=int)

    for time_step_id, ref_airplane_at_time_step in enumerate(ref_airplanes):
        ref_wing_at_time_step = ref_airplane_at_time_step.wings[ref_wing_id]
        ref_root_wing_cross_section_at_time_step = (
            ref_wing_at_time_step.wing_cross_sections[ref_root_wing_cross_section_id]
        )
        ref_tip_wing_cross_section_at_time_step = (
            ref_wing_at_time_step.wing_cross_sections[ref_tip_wing_cross_section_id]
        )

        _logger.debug(
            _logging.indent()
            + "Calculating the number of spanwise Panels for time step "
            f"{time_step_id+1}/{num_time_steps}"
        )

        num_spanwise_panels_at_step = _get_wing_section_num_spanwise_panels(
            desired_average_panel_aspect_ratio=desired_average_panel_aspect_ratio,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
            ref_root_wing_cross_section=ref_root_wing_cross_section_at_time_step,
            ref_tip_wing_cross_section=ref_tip_wing_cross_section_at_time_step,
            start_val=start_val,
        )

        these_num_spanwise_panels[time_step_id] = num_spanwise_panels_at_step

        _logger.debug(
            _logging.indent()
            + f"Number of spanwise Panels: {num_spanwise_panels_at_step}"
        )

    return int(max(these_num_spanwise_panels))


def _get_wing_section_num_spanwise_panels(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_root_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    ref_tip_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    start_val: int,
) -> int:
    """Calculates the number of spanwise Panels to use for the wing section of a Wing
    based on a desired average Panel aspect ratio.

    :param desired_average_panel_aspect_ratio: The target average Panel aspect ratio to
        achieve. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes). It must be a positive int.
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param chordwise_spacing: The type of spacing between the chordwise Panels. Can be
        "cosine" or "uniform".
    :param ref_root_wing_cross_section: The root WingCrossSection of the wing section.
    :param ref_tip_wing_cross_section: The tip WingCrossSection of the wing section.
    :param start_val: The initial number of spanwise Panels to start the search from. It
        must be a positive int. Using a higher value can speed up the search if a lower
        bound is already known.
    :return: The number of spanwise Panels that results in an average Panel aspect ratio
        closest to the desired value.
    """
    this_num_spanwise_panels = start_val
    average_panel_aspect_ratios = []

    while True:
        this_average_panel_aspect_ratio = _get_wing_section_average_panel_aspect_ratio(
            num_chordwise_panels,
            chordwise_spacing,
            ref_root_wing_cross_section,
            ref_tip_wing_cross_section,
            num_spanwise_panels=this_num_spanwise_panels,
        )
        average_panel_aspect_ratios.append(this_average_panel_aspect_ratio)

        if this_average_panel_aspect_ratio <= desired_average_panel_aspect_ratio:
            break

        this_num_spanwise_panels += 1

    if len(average_panel_aspect_ratios) < 2:
        return this_num_spanwise_panels

    this_aspect_ratio_difference = abs(
        average_panel_aspect_ratios[-1] - desired_average_panel_aspect_ratio
    )
    last_aspect_ratio_difference = abs(
        average_panel_aspect_ratios[-2] - desired_average_panel_aspect_ratio
    )

    if last_aspect_ratio_difference < this_aspect_ratio_difference:
        this_num_spanwise_panels -= 1
    return this_num_spanwise_panels


def _get_wing_section_average_panel_aspect_ratio(
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_root_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    ref_tip_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    num_spanwise_panels: int,
) -> float:
    """Calculates the average aspect ratio of Panels in a wing section with a particular
    number of chordwise and spanwise Panels.

    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param chordwise_spacing: The type of spacing between the chordwise Panels. Can be
        "cosine" or "uniform".
    :param ref_root_wing_cross_section: The root WingCrossSection of the wing section.
    :param ref_tip_wing_cross_section: The tip WingCrossSection of the wing section.
    :param num_spanwise_panels: The number of spanwise Panels to use. It must be a
        positive int.
    :return: The average Panel aspect ratio for the wing section with the given Panel
        counts. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes).
    """
    this_airplane = geometry.airplane.Airplane(
        wings=[
            geometry.wing.Wing(
                wing_cross_sections=[
                    geometry.wing_cross_section.WingCrossSection(
                        airfoil=geometry.airfoil.Airfoil(
                            name=ref_root_wing_cross_section.airfoil.name,
                            outline_A_lp=ref_root_wing_cross_section.airfoil.outline_A_lp,
                            resample=ref_root_wing_cross_section.airfoil.resample,
                            n_points_per_side=ref_root_wing_cross_section.airfoil.n_points_per_side,
                        ),
                        num_spanwise_panels=num_spanwise_panels,
                        chord=ref_root_wing_cross_section.chord,
                        spanwise_spacing=ref_root_wing_cross_section.spanwise_spacing,
                    ),
                    geometry.wing_cross_section.WingCrossSection(
                        airfoil=geometry.airfoil.Airfoil(
                            name=ref_tip_wing_cross_section.airfoil.name,
                            outline_A_lp=ref_tip_wing_cross_section.airfoil.outline_A_lp,
                            resample=ref_tip_wing_cross_section.airfoil.resample,
                            n_points_per_side=ref_tip_wing_cross_section.airfoil.n_points_per_side,
                        ),
                        num_spanwise_panels=None,
                        chord=ref_tip_wing_cross_section.chord,
                        Lp_Wcsp_Lpp=ref_tip_wing_cross_section.Lp_Wcsp_Lpp,
                        angles_Wcsp_to_Wcs_ixyz=ref_tip_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                    ),
                ],
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing,
            )
        ]
    )

    _average_panel_aspect_ratio = this_airplane.wings[0].average_panel_aspect_ratio
    assert _average_panel_aspect_ratio is not None

    return _average_panel_aspect_ratio


def _build_edge_defined_wing(
    ref_wing: geometry.wing.Wing,
    num_chordwise_panels: int,
    num_wing_cross_sections: int,
) -> geometry.wing.Wing:
    """Rebuilds an edge-defined Wing from its stored edge curves at a new resolution.

    The reference Wing's stored leading edge and trailing edge curves are resampled into
    the requested number of WingCrossSections with Wing.from_edge_points, reusing every
    other geometric parameter (position, orientation, symmetry, chordwise spacing, and
    tip trim) from the reference Wing. The single Airfoil is read from the reference
    Wing's first WingCrossSection and shared across the rebuilt WingCrossSections, which
    is safe because Airfoils are immutable.

    :param ref_wing: The reference edge-defined Wing whose stored edge curves and
        parameters are reused. Its spanwise_mesh must be "edge_defined".
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param num_wing_cross_sections: The number of WingCrossSections to resample the edge
        curves into. It must be an int of at least 2.
    :return: The rebuilt edge-defined Wing.
    """
    # An edge-defined Wing always stores its edge curves and tip trim fraction, so these
    # are never None here. The assertions narrow the accessors' optional types.
    leadingEdgePoints_Wn_Ler = ref_wing.leadingEdgePoints_Wn_Ler
    trailingEdgePoints_Wn_Ler = ref_wing.trailingEdgePoints_Wn_Ler
    tip_trim_fraction = ref_wing.tip_trim_fraction
    assert leadingEdgePoints_Wn_Ler is not None
    assert trailingEdgePoints_Wn_Ler is not None
    assert tip_trim_fraction is not None

    return geometry.wing.Wing.from_edge_points(
        leadingEdgePoints_Wn_Ler=leadingEdgePoints_Wn_Ler,
        trailingEdgePoints_Wn_Ler=trailingEdgePoints_Wn_Ler,
        num_wing_cross_sections=num_wing_cross_sections,
        airfoil=ref_wing.wing_cross_sections[0].airfoil,
        name=ref_wing.name,
        Ler_Gs_Cgs=ref_wing.Ler_Gs_Cgs,
        angles_Gs_to_Wn_ixyz=ref_wing.angles_Gs_to_Wn_ixyz,
        symmetric=ref_wing.symmetric,
        mirror_only=ref_wing.mirror_only,
        symmetryNormal_G=ref_wing.symmetryNormal_G,
        symmetryPoint_G_Cg=ref_wing.symmetryPoint_G_Cg,
        num_chordwise_panels=num_chordwise_panels,
        chordwise_spacing=ref_wing.chordwise_spacing,
        tip_trim_fraction=tip_trim_fraction,
    )


def _get_num_wing_cross_sections_for_panel_ar(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    ref_wing: geometry.wing.Wing,
    start_val: int,
) -> int:
    """Calculates the number of WingCrossSections that gives an edge-defined Wing a
    target average Panel aspect ratio.

    An edge-defined Wing is refined by resampling its stored edge curves into more
    WingCrossSections, not by adding spanwise Panels to fixed WingCrossSections, so the
    number of WingCrossSections is the refinement knob. Increasing it monotonically
    lowers the meshed average Panel aspect ratio.

    Rather than estimate the count from a closed form, this searches for it by meshing.
    A closed form built from the mean chord would undershoot the count on a tapered
    planform, because average_panel_aspect_ratio is the mean of each Panel's individual
    aspect ratio and is dominated by the small-chord tip Panels. Instead, the reference
    Wing is rebuilt at trial WingCrossSection counts (with _build_edge_defined_wing),
    each trial is meshed, and its average Panel aspect ratio is read. The search
    proportionally jumps to bracket the target, then bisects to the crossing, and
    returns whichever of the two counts straddling the target gives the closer average
    Panel aspect ratio. This matches how the trapezoidal
    _get_wing_section_num_spanwise_panels searches, so a target Panel aspect ratio means
    the same physical thing for an edge-defined Wing as for a trapezoidal one. Only
    meshing is performed here, never solving.

    :param desired_average_panel_aspect_ratio: The target average Panel aspect ratio to
        achieve. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes). It must be a positive int.
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param ref_wing: The reference edge-defined Wing whose stored edge curves are
        resampled. Its spanwise_mesh must be "edge_defined".
    :param start_val: The initial number of WingCrossSections to start the search from.
        It must be a positive int and a lower bound on the result. Using a higher value
        can speed up the search if a lower bound is already known.
    :return: The number of WingCrossSections that gives an average Panel aspect ratio
        closest to the desired value.
    """
    meshed_average_panel_aspect_ratios: dict[int, float] = {}

    def average_panel_aspect_ratio_at(num_wing_cross_sections: int) -> float:
        if num_wing_cross_sections not in meshed_average_panel_aspect_ratios:
            refined_wing = _build_edge_defined_wing(
                ref_wing, num_chordwise_panels, num_wing_cross_sections
            )
            refined_airplane = geometry.airplane.Airplane(wings=[refined_wing])
            this_average_panel_aspect_ratio = refined_airplane.wings[
                0
            ].average_panel_aspect_ratio
            assert this_average_panel_aspect_ratio is not None
            meshed_average_panel_aspect_ratios[num_wing_cross_sections] = (
                this_average_panel_aspect_ratio
            )
        return meshed_average_panel_aspect_ratios[num_wing_cross_sections]

    # from_edge_points requires at least two WingCrossSections.
    lower_num = max(int(start_val), 2)
    lower_average_panel_aspect_ratio = average_panel_aspect_ratio_at(lower_num)

    # If the starting count already meets the target, it is the smallest count that
    # does, so return it.
    if lower_average_panel_aspect_ratio <= desired_average_panel_aspect_ratio:
        return lower_num

    # Proportionally jump to bracket the target from above. The average Panel aspect
    # ratio scales roughly as 1.0 / (num_wing_cross_sections - 1), so this estimate
    # lands near the crossing in a few steps. Each jump strictly increases the count,
    # and the aspect ratio falls toward zero as the count grows, so the loop terminates.
    upper_num = lower_num
    upper_average_panel_aspect_ratio = lower_average_panel_aspect_ratio
    while upper_average_panel_aspect_ratio > desired_average_panel_aspect_ratio:
        projected_num = 1 + round(
            upper_average_panel_aspect_ratio
            * (upper_num - 1)
            / desired_average_panel_aspect_ratio
        )
        upper_num = max(projected_num, upper_num + 1)
        upper_average_panel_aspect_ratio = average_panel_aspect_ratio_at(upper_num)

    # Bisect between the count above the target (lower_num) and the count at or below it
    # (upper_num) for the integer crossing.
    while upper_num - lower_num > 1:
        middle_num = (lower_num + upper_num) // 2
        if (
            average_panel_aspect_ratio_at(middle_num)
            > desired_average_panel_aspect_ratio
        ):
            lower_num = middle_num
        else:
            upper_num = middle_num

    # Return whichever of the two straddling counts gives the closer average Panel
    # aspect ratio. On a tie, prefer the finer mesh (upper_num), matching the
    # trapezoidal search.
    lower_difference = abs(
        average_panel_aspect_ratio_at(lower_num) - desired_average_panel_aspect_ratio
    )
    upper_difference = abs(
        average_panel_aspect_ratio_at(upper_num) - desired_average_panel_aspect_ratio
    )
    if lower_difference < upper_difference:
        return lower_num
    return upper_num


def _resolve_num_wing_cross_sections(
    panel_aspect_ratio_id: int,
    num_chordwise_panels_id: int,
    airplane_id: int,
    wing_id: int,
    airplane_name: str,
    wing_name: str,
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
    compute_num_wing_cross_sections: Callable[[int], int],
) -> int:
    """Resolves the number of WingCrossSections for one edge-defined Wing, using and
    updating the shared cache.

    The result for a given (Panel aspect ratio, number of chordwise Panels, Airplane,
    Wing) combination is returned from the cache if present. Otherwise, the search
    starts from a conservative lower bound (the smallest number of WingCrossSections
    already found for this Wing at an incrementally coarser mesh, since the current
    finer mesh must need at least that many), ``compute_num_wing_cross_sections`` is
    called with that starting value to find the count, and the result is cached.

    :param panel_aspect_ratio_id: The index of the current Panel aspect ratio within the
        list of Panel aspect ratios being tested.
    :param num_chordwise_panels_id: The index of the current number of chordwise Panels
        within the list of numbers of chordwise Panels being tested.
    :param airplane_id: The index of the current Airplane.
    :param wing_id: The index of the current Wing within the Airplane.
    :param airplane_name: The name of the current Airplane, used in the log messages.
    :param wing_name: The name of the current Wing, used in the log messages.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections. It is read and updated in
        place.
    :param compute_num_wing_cross_sections: A callable that takes a starting number of
        WingCrossSections and returns the number of WingCrossSections that achieves the
        desired Panel aspect ratio. It is called only on a cache miss.
    :return: The number of WingCrossSections for the Wing.
    """
    num_wing_cross_sections_key = (
        panel_aspect_ratio_id,
        num_chordwise_panels_id,
        airplane_id,
        wing_id,
    )

    _logger.debug(_logging.indent() + f"{airplane_name}'s {wing_name}:")

    if num_wing_cross_sections_key in num_wing_cross_sections_cache:
        _logger.debug(
            _logging.indent(1) + "Getting the cached number of WingCrossSections"
        )
        this_num_wing_cross_sections = num_wing_cross_sections_cache[
            num_wing_cross_sections_key
        ]
    else:
        # Start the search from a conservative lower bound: the smallest number of
        # WingCrossSections already found for this Wing at an incrementally coarser mesh
        # (in Panel aspect ratio, number of chordwise Panels, or both), since the
        # current finer mesh must use at least that many. The floor is two, the fewest
        # WingCrossSections an edge-defined Wing can have.
        starting_num_wing_cross_sections = 2
        last_ar_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id,
            airplane_id,
            wing_id,
        )
        last_chord_key = (
            panel_aspect_ratio_id,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
        )
        last_ar_and_chord_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
        )
        last_cache_val = min(
            num_wing_cross_sections_cache.get(last_ar_key, np.inf),
            num_wing_cross_sections_cache.get(last_chord_key, np.inf),
            num_wing_cross_sections_cache.get(last_ar_and_chord_key, np.inf),
        )
        if last_cache_val != np.inf:
            starting_num_wing_cross_sections = int(last_cache_val)

        _logger.debug(
            _logging.indent(1) + "Calculating the number of WingCrossSections"
        )
        _logger.debug(
            _logging.indent(1)
            + f"Starting the search at {starting_num_wing_cross_sections}"
        )

        # The search's log messages nest one level under the detail lines above.
        with _logging.nested(2):
            this_num_wing_cross_sections = compute_num_wing_cross_sections(
                starting_num_wing_cross_sections
            )

        num_wing_cross_sections_cache[num_wing_cross_sections_key] = (
            this_num_wing_cross_sections
        )

    _logger.debug(
        _logging.indent(1)
        + f"Number of WingCrossSections: {this_num_wing_cross_sections}"
    )

    return this_num_wing_cross_sections


def _resolve_num_spanwise_panels(
    panel_aspect_ratio_id: int,
    num_chordwise_panels_id: int,
    airplane_id: int,
    wing_id: int,
    wing_cross_section_id: int,
    airplane_name: str,
    wing_name: str,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    compute_num_spanwise_panels: Callable[[int], int],
) -> int:
    """Resolves the number of spanwise Panels for one Wing section, using and updating
    the shared cache.

    The result for a given (Panel aspect ratio, number of chordwise Panels, Airplane,
    Wing, WingCrossSection) combination is returned from the cache if present.
    Otherwise, the search starts from a conservative lower bound (the smallest number of
    spanwise Panels already found for this Wing section at an incrementally coarser
    mesh, since the current finer mesh must need at least that many),
    ``compute_num_spanwise_panels`` is called with that starting value to find the
    count, and the result is cached.

    :param panel_aspect_ratio_id: The index of the current Panel aspect ratio within the
        list of Panel aspect ratios being tested.
    :param num_chordwise_panels_id: The index of the current number of chordwise Panels
        within the list of numbers of chordwise Panels being tested.
    :param airplane_id: The index of the current Airplane.
    :param wing_id: The index of the current Wing within the Airplane.
    :param wing_cross_section_id: The index of the current WingCrossSection within the
        Wing.
    :param airplane_name: The name of the current Airplane, used in the log messages.
    :param wing_name: The name of the current Wing, used in the log messages.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels. It is read
        and updated in place.
    :param compute_num_spanwise_panels: A callable that takes a starting number of
        spanwise Panels and returns the number of spanwise Panels that achieves the
        desired Panel aspect ratio. It is called only on a cache miss.
    :return: The number of spanwise Panels for the Wing section.
    """
    num_spanwise_panels_key = (
        panel_aspect_ratio_id,
        num_chordwise_panels_id,
        airplane_id,
        wing_id,
        wing_cross_section_id,
    )

    _logger.debug(
        _logging.indent() + f"{airplane_name}'s {wing_name}, "
        f"WingCrossSection #{wing_cross_section_id + 1}:"
    )

    if num_spanwise_panels_key in num_spanwise_panels_cache:
        _logger.debug(
            _logging.indent(1) + "Getting the cached number of spanwise Panels"
        )
        this_num_spanwise_panels = num_spanwise_panels_cache[num_spanwise_panels_key]
    else:
        # Start the search from a conservative lower bound: the smallest number of
        # spanwise Panels already found for this Wing section at an incrementally
        # coarser mesh (in Panel aspect ratio, number of chordwise Panels, or both),
        # since the current finer mesh must use at least that many.
        starting_num_spanwise_panels = 1
        last_ar_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id,
            airplane_id,
            wing_id,
            wing_cross_section_id,
        )
        last_chord_key = (
            panel_aspect_ratio_id,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
            wing_cross_section_id,
        )
        last_ar_and_chord_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
            wing_cross_section_id,
        )
        last_cache_val = min(
            num_spanwise_panels_cache.get(last_ar_key, np.inf),
            num_spanwise_panels_cache.get(last_chord_key, np.inf),
            num_spanwise_panels_cache.get(last_ar_and_chord_key, np.inf),
        )
        if last_cache_val != np.inf:
            starting_num_spanwise_panels = int(last_cache_val)

        _logger.debug(_logging.indent(1) + "Calculating the number of spanwise Panels")
        _logger.debug(
            _logging.indent(1)
            + f"Starting the search at {starting_num_spanwise_panels}"
        )

        # The search's log messages nest one level under the detail lines above.
        with _logging.nested(2):
            this_num_spanwise_panels = compute_num_spanwise_panels(
                starting_num_spanwise_panels
            )

        num_spanwise_panels_cache[num_spanwise_panels_key] = this_num_spanwise_panels

    _logger.debug(
        _logging.indent(1) + f"Number of spanwise Panels: {this_num_spanwise_panels}"
    )

    return this_num_spanwise_panels
