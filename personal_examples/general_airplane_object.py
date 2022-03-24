import pterasoftware as ps

default_airplane = ps.geometry.Airplane(
    wings=[
        ps.geometry.Wing(
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    airfoil=ps.geometry.Airfoil(
                        name="naca2412",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
        ps.geometry.Wing(
            x_le=7.50,
            z_le=0.25,
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    chord=0.5,
                    twist=-5.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=1.0,
                    chord=0.5,
                    twist=-5.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

default_operating_point = ps.operating_point.OperatingPoint()

default_problem = ps.problems.SteadyProblem(
    airplanes=[default_airplane], operating_point=default_operating_point
)

default_solver = (
    ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=default_problem
    )
)

default_solver.run()

ps.output.draw(solver=default_solver, scalar_type="lift", show_streamlines=True)
ps.output.print_steady_results(steady_solver=default_solver)

ps.trim.analyze_steady_trim(problem=default_problem)

print("done!")
