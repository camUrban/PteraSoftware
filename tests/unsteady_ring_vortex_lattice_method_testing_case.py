"""This is a testing case for the unsteady ring vortex lattice method solver.

"""

import aviansoftwareminimumviableproduct as asmvp

# Initialize the problem's geometry.
unsteady_solver_validation_airplane = asmvp.geometry.Airplane(

    # Name the airplane.
    name="Unsteady Solver Testing Airplane",

    # Define a list of the airplane's wings.
    wings=[

        # Initialize the wing object.
        asmvp.geometry.Wing(

            # Name the wing.
            name="Wing",

            # This will be a symmetrical wing.
            symmetric=True,

            # Define a list of the wing's cross sections.
            cross_sections=[

                # Initialize the root cross section object.
                asmvp.geometry.WingCrossSection(
                    # Define the cross section's twist chord.
                    chord=1.0,

                    # Initialize this cross section's airfoil object.
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                ),

                # Initialize the tip cross section object.
                asmvp.geometry.WingCrossSection(

                    # Define the cross section's leading edge placement.
                    x_le=1.0,
                    y_le=5.0,
                    z_le=0.0,

                    # Define the cross section's twist and chord.
                    twist=5.0,
                    chord=0.75,

                    # Initialize this cross section's airfoil object.
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                )
            ]
        )
    ]
)

# Initialize the problem's operating point.
unsteady_solver_validation_operating_point = asmvp.performance.OperatingPoint()

# Initialize the problem's operating point.
unsteady_solver_validation_movement = None

# Initialize the problem.
unsteady_solver_validation_problem = asmvp.problems.UnsteadyProblem(
    airplane=unsteady_solver_validation_airplane,
    operating_point=unsteady_solver_validation_operating_point,
    movement=unsteady_solver_validation_movement,
    simulation_duration=1.0,
    simulation_time_step=0.1
)

# Delete placeholder objects that are now stored in the problem object.
del unsteady_solver_validation_airplane
del unsteady_solver_validation_operating_point

# Initialize the solver.
unsteady_solver_validation_solver = asmvp.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    unsteady_solver_validation_problem)

# Delete placeholder object that are now stored in the simulation object.
del unsteady_solver_validation_problem

# Run the solver and draw it's output.
unsteady_solver_validation_solver.run()
asmvp.output.draw(unsteady_solver_validation_solver.airplane, show_delta_pressures=True)
