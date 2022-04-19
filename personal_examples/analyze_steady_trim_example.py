# ToDo: Document this script.
import logging

import pterasoftware as ps

example_logger = logging.getLogger("example")
example_logger.setLevel(logging.DEBUG)

default_airplane = ps.geometry.Airplane(
    weight=250,
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

default_operating_point = ps.operating_point.OperatingPoint(external_thrust=5)

default_problem = ps.problems.SteadyProblem(
    airplanes=[default_airplane], operating_point=default_operating_point
)

default_solver = (
    ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=default_problem
    )
)
default_solver.run()

trim_conditions = ps.trim.analyze_steady_trim(
    problem=default_problem,
    velocity_bounds=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-1, 1),
    external_thrust_bounds=(0, 10),
)

example_logger.info("Trim Velocity:\t\t\t%.2f m/s" % trim_conditions[0])
example_logger.info("Trim Alpha:\t\t\t%.2f deg" % trim_conditions[1])
example_logger.info("Trim Beta:\t\t\t\t%.2f deg" % trim_conditions[2])
example_logger.info("Trim External Thrust:\t%.2f N" % trim_conditions[3])
