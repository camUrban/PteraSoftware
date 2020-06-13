import aviansoftwareminimumviableproduct as asmvp

steady_solver_validation_airplane = asmvp.geometry.Airplane(
    name="Unsteady Solver Testing Airplane",
    wings=[
        asmvp.geometry.Wing(
            name="Wing",
            symmetric=True,
            num_chordwise_panels=3,
            wing_cross_sections=[
                asmvp.geometry.WingCrossSection(
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    num_spanwise_panels=5
                ),
                asmvp.geometry.WingCrossSection(
                    y_le=1.0,
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    num_spanwise_panels=5
                )
            ]
        )
    ]
)

steady_solver_validation_operating_point = asmvp.operating_point.OperatingPoint()

steady_solver_validation_problem = asmvp.problems.SteadyProblem(
    airplane=steady_solver_validation_airplane,
    operating_point=steady_solver_validation_operating_point
)

del steady_solver_validation_operating_point

steady_solver_validation_solver = asmvp.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
    steady_solver_validation_problem)

del steady_solver_validation_problem

steady_solver_validation_solver.run()

asmvp.output.draw(airplane=steady_solver_validation_airplane, show_delta_pressures=False)
