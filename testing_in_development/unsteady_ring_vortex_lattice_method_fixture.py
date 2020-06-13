import aviansoftwareminimumviableproduct as asmvp

unsteady_solver_validation_airplane = asmvp.geometry.Airplane(
    name="Unsteady Solver Testing Airplane",
    # Define a list of the current_airplane's wings.
    wings=[

        # Initialize the wing object.
        asmvp.geometry.Wing(

            # Name the wing.
            name="Wing",

            symmetric=False,
            num_chordwise_panels=6,
            chordwise_spacing='uniform',

            # Define a list of the wing's cross sections.
            wing_cross_sections=[

                # Initialize the root cross section object.
                asmvp.geometry.WingCrossSection(
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    num_spanwise_panels=24
                ),

                # Initialize the tip cross section object.
                asmvp.geometry.WingCrossSection(
                    y_le=4.0,
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    num_spanwise_panels=24
                )
            ]
        )
    ]
)

unsteady_solver_validation_operating_point = asmvp.operating_point.OperatingPoint()

unsteady_solver_validation_root_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
    base_wing_cross_section=unsteady_solver_validation_airplane.wings[0].wing_cross_sections[0]
)

unsteady_solver_validation_tip_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
    base_wing_cross_section=unsteady_solver_validation_airplane.wings[0].wing_cross_sections[1]
)

unsteady_solver_validation_wing_movement = asmvp.movement.WingMovement(
    base_wing=unsteady_solver_validation_airplane.wings[0],
    wing_cross_sections_movements=[unsteady_solver_validation_root_wing_cross_section_movement,
                                   unsteady_solver_validation_tip_wing_cross_section_movement]
)

del unsteady_solver_validation_root_wing_cross_section_movement
del unsteady_solver_validation_tip_wing_cross_section_movement

unsteady_solver_validation_airplane_movement = asmvp.movement.AirplaneMovement(
    base_airplane=unsteady_solver_validation_airplane,
    wing_movements=[unsteady_solver_validation_wing_movement]
)

del unsteady_solver_validation_wing_movement

unsteady_solver_validation_operating_point_movement = asmvp.movement.OperatingPointMovement(
    base_operating_point=unsteady_solver_validation_operating_point
)

unsteady_solver_validation_movement = asmvp.movement.Movement(
    airplane_movement=unsteady_solver_validation_airplane_movement,
    operating_point_movement=unsteady_solver_validation_operating_point_movement,
    num_steps=60,
    delta_time=0.0166
)

del unsteady_solver_validation_airplane_movement
del unsteady_solver_validation_operating_point_movement

# asmvp.output.make_flapping_gif(unsteady_solver_validation_movement)

unsteady_solver_validation_problem = asmvp.problems.UnsteadyProblem(
    airplane=unsteady_solver_validation_airplane,
    operating_point=unsteady_solver_validation_operating_point,
    movement=unsteady_solver_validation_movement
)

del unsteady_solver_validation_airplane
del unsteady_solver_validation_operating_point
del unsteady_solver_validation_movement

unsteady_solver_validation_solver = asmvp.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    unsteady_solver_validation_problem)

del unsteady_solver_validation_problem

unsteady_solver_validation_solver.run()
