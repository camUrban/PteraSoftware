import aviansoftwareminimumviableproduct as asmvp


this_airplane = asmvp.geometry.Airplane(
    # Name the current_airplane.
    name="Steady Solver Testing Airplane",
    # Define a list of the current_airplane's wings.
    wings=[
        # Initialize the wing object.
        asmvp.geometry.Wing(
            # Name the wing.
            name="Wing",
            # This will be a symmetrical wing.
            symmetric=True,
            # Define a list of the wing's cross sections.
            wing_cross_sections=[
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
                ),
            ],
        )
    ],
)

this_operating_point = asmvp.operating_point.OperatingPoint()

this_steady_problem = asmvp.problems.SteadyProblem(
    airplane=this_airplane, operating_point=this_operating_point
)

this_steady_solver = asmvp.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
    steady_problem=this_steady_problem
)

asmvp.output.draw(
    airplane=this_airplane,
    show_delta_pressures=False,
    show_streamlines=False,
    show_wake_vortices=False,
)
