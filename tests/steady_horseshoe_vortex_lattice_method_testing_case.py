# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp

steady_solver_validation_airplane = asmvp.geometry.Airplane(
    name="Steady Solver Testing Airplane",
    x_ref=0,
    y_ref=0,
    z_ref=0,
    wings=[
        asmvp.geometry.Wing(
            name="Wing",
            x_le=0,
            y_le=0,
            z_le=0,
            symmetric=True,
            chordwise_spacing="cosine",
            cross_sections=[
                asmvp.geometry.WingCrossSection(
                    x_le=0,
                    y_le=0,
                    z_le=0,
                    twist=0,
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    spanwise_spacing="cosine"
                ),
                asmvp.geometry.WingCrossSection(
                    x_le=1.0,
                    y_le=5,
                    z_le=0.5,
                    twist=5.0,
                    chord=0.75,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    spanwise_spacing="cosine"
                )
            ]
        )
    ]
)

steady_solver_validation_operating_point = asmvp.performance.OperatingPoint()
steady_solver_validation_problem = asmvp.problems.SteadyProblem(
    airplane=steady_solver_validation_airplane, operating_point=steady_solver_validation_operating_point)

steady_solver_validation_simulation = (asmvp.steady_horseshoe_vortex_lattice_method.
                                       SteadyHorseshoeVortexLatticeMethodSolver(steady_solver_validation_problem))

steady_solver_validation_simulation.run()
asmvp.output.draw(steady_solver_validation_airplane, show_delta_pressures=True)
