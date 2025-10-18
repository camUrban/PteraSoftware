# REFACTOR: I haven't yet started refactoring this module.
import pterasoftware as ps

flapping_frequency = 1
num_chordwise_panels = 5
num_spanwise_panels = 20

example_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    chord=1.75,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    control_surface_symmetry_type="symmetric",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    spanwise_spacing=None,
                    Lp_Wcsp_Lpp=(0.625, 5.0, 0.0),
                    chord=0.5,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        ),
    ],
    name="Example Airplane",
)

upper_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0]
    )
)

upper_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1]
    )
)

upper_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_section_movements=[
        upper_wing_root_wing_cross_section_movement,
        upper_wing_tip_wing_cross_section_movement,
    ],
)

del upper_wing_root_wing_cross_section_movement
del upper_wing_tip_wing_cross_section_movement

airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=example_airplane, wing_movements=[upper_wing_movement]
)

del example_airplane
del upper_wing_movement

example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0
)

operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point
)

del example_operating_point

movement = ps.movements.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    num_chords=5,
)

unsteady_problem = ps.problems.UnsteadyProblem(
    movement=movement, only_final_results=True
)

del airplane_movement
del operating_point_movement

ps.convergence.analyze_unsteady_convergence(
    ref_problem=unsteady_problem,
    prescribed_wake=True,
    free_wake=True,
    num_cycles_bounds=(1, 3),
    num_chords_bounds=(1, 5),
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(5, 10),
    convergence_criteria=5.0,
)

# Converged result:
# wake=false
# cycles=2
# ar=4
# panels=9
