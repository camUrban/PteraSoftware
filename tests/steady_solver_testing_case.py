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
            chordwise_spacing="uniform",
            cross_sections=[
                asmvp.geometry.WingCrossSection(
                    x_le=0,
                    y_le=0,
                    z_le=0,
                    twist=0,
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca0010"),
                    spanwise_spacing="uniform"
                ),
                asmvp.geometry.WingCrossSection(
                    x_le=0,
                    y_le=5,
                    z_le=0,
                    twist=0,
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca0010"),
                    spanwise_spacing="uniform"
                )
            ]
        )
    ]
)
# steady_solver_validation_airplane.set_paneling_everywhere(10, 10)

steady_solver_validation_operating_point = asmvp.performance.OperatingPoint()

steady_solver_validation_simulation = asmvp.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
    steady_solver_validation_airplane, steady_solver_validation_operating_point)

steady_solver_validation_simulation.run()
asmvp.output.draw_airplane(steady_solver_validation_airplane)
