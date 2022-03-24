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

trim_conditions = ps.trim.analyze_steady_trim(problem=default_problem)

print("Trim Velocity:\t%.2f m/s" % trim_conditions[0])
print("Trim Alpha:\t\t%.2f deg" % trim_conditions[1])
print("Trim Beta:\t\t%.2f deg" % trim_conditions[2])
print("Trim Thrust:\t%.2f N" % trim_conditions[3])
