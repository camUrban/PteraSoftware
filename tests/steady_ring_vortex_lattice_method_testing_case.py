"""This is a testing case for the steady ring vortex lattice method solver.

    Based on an identical XFLR5 testing case, the expected output for this case is:
        CL:     0.788
        CDi:    0.019
        Cl:     -0.000
        Cm:     -0.687
        Cn:     -0.000

    Note: The expected output was created using XFLR5's VLM2 analysis type, which is a ring vortex lattice method
    solver.
"""

import aviansoftwareminimumviableproduct as asmvp

# Initialize the problem's geometry.
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
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    twist=0.0,
                    chord=1.0,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    spanwise_spacing="cosine"
                ),
                asmvp.geometry.WingCrossSection(
                    x_le=1.0,
                    y_le=5.0,
                    z_le=0.0,
                    twist=5.0,
                    chord=0.75,
                    airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    spanwise_spacing="cosine"
                )
            ]
        )
    ]
)

# Initialize the problem's operating point.
steady_solver_validation_operating_point = asmvp.performance.OperatingPoint()

# Initialize the problem.
steady_solver_validation_problem = asmvp.problems.SteadyProblem(
    airplane=steady_solver_validation_airplane, operating_point=steady_solver_validation_operating_point)

# Delete placeholder objects that are now stored in the problem object.
del steady_solver_validation_airplane
del steady_solver_validation_operating_point

# Initialize the solver.
steady_solver_validation_solver = asmvp.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
    steady_solver_validation_problem)

# Delete placeholder object that are now stored in the simulation object.
del steady_solver_validation_problem

# Run the solver and draw it's output.
steady_solver_validation_solver.run()
asmvp.output.draw(steady_solver_validation_solver.airplane, show_delta_pressures=True)
