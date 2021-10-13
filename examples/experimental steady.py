"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """

import src

experimental_airplane = src.geometry.Airplane(
    name="Experimental Airplane",
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0875,
                    y_le=0.350,
                    chord=0.0219,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
    ],
)

experimental_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
    velocity=1.0,
    alpha=5.0,
    nu=15.06e-6,
)

experimental_problem = src.problems.SteadyProblem(
    airplanes=[experimental_airplane], operating_point=experimental_operating_point
)

del experimental_airplane
del experimental_operating_point

experimental_solver = (
    src.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=experimental_problem,
    )
)

del experimental_problem

experimental_solver.run()

src.output.print_steady_results(steady_solver=experimental_solver)

src.output.draw(
    solver=experimental_solver, show_delta_pressures=True, show_streamlines=True
)
