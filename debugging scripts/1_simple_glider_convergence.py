"""This script is for finding the converged parameters for the simple glider."""

from turtledemo.penrose import inflatedart

import numpy as np

import pterasoftware as ps

simple_glider_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=5,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            chordwise_spacing="uniform",
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=5,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=5,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 0.5),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=5,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=5,
                    control_surface_symmetry_type=None,
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 2.0, 0.0),
                    control_surface_symmetry_type=None,
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 1.0),
            angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=5,
        ),
    ],
    weight=420,
)

simple_glider_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=simple_glider_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[2],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

del simple_glider_airplane

simple_glider_operating_point = ps.operating_point.OperatingPoint(
    externalFX_W=7.5,
)

simple_glider_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=simple_glider_operating_point
    )
)

del simple_glider_operating_point

simple_glider_movement = ps.movements.movement.Movement(
    airplane_movements=[simple_glider_airplane_movement],
    operating_point_movement=simple_glider_operating_point_movement,
    num_chords=5,
)

del simple_glider_airplane_movement
del simple_glider_operating_point_movement

simple_glider_problem = ps.problems.UnsteadyProblem(
    movement=simple_glider_movement,
    only_final_results=True,
)

# simple_glider_solver = (
#     ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
#         unsteady_problem=simple_glider_problem
#     )
# )
#
# ps.output.draw(simple_glider_solver)
#
# del simple_glider_solver

ps.convergence.analyze_unsteady_convergence(
    ref_problem=simple_glider_problem,
    prescribed_wake=True,
    free_wake=False,
    num_chords_bounds=(13 - 2, 13 + 2),
    panel_aspect_ratio_bounds=(1 + 2, 1),
    num_chordwise_panels_bounds=(6 - 2, 6 + 2),
    convergence_criteria=1.0,
)

# prescribed_wake = True
# num_chords = 13
# panel_aspect_ratio = 3
# num_chordwise_panels = 7
# num_spanwise_panels:
# Wing 1:
#     WingCrossSection 1: 12
#     WingCrossSection 2: None
# Wing 2:
#     WingCrossSection 1: 2
#     WingCrossSection 2: None
# Wing 3:
#     WingCrossSection 1: 5
#     WingCrossSection 2: None
