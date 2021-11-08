"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on three airplanes (with geometry similar to the NMT
experimental setup) flying in formation."""
import src

horseshoe_vortex_method = False

alpha = 0
x_spacing = 0.5
y_spacing = 0.5

experimental_lead_airplane = src.geometry.Airplane(
    name="Experimental Lead Airplane",
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    twist=alpha,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.350,
                    chord=0.0219,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
    ],
)

experimental_right_airplane = src.geometry.Airplane(
    name="Experimental Right Airplane",
    x_ref=x_spacing,
    y_ref=y_spacing,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            x_le=x_spacing,
            y_le=y_spacing,
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    twist=alpha,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.350,
                    chord=0.0219,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
    ],
)

experimental_left_airplane = src.geometry.Airplane(
    name="Experimental Left Airplane",
    x_ref=x_spacing,
    y_ref=-y_spacing,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            x_le=x_spacing,
            y_le=-y_spacing,
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    twist=alpha,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
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
    alpha=0.0,
    nu=15.06e-6,
)


experimental_problem = src.problems.SteadyProblem(
    airplanes=[
        experimental_lead_airplane,
        experimental_right_airplane,
        experimental_left_airplane,
    ],
    operating_point=experimental_operating_point,
)

del experimental_lead_airplane
del experimental_right_airplane
del experimental_left_airplane
del experimental_operating_point

if horseshoe_vortex_method:
    experimental_solver = src.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=experimental_problem,
    )
else:
    experimental_solver = (
        src.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
            steady_problem=experimental_problem,
        )
    )

del experimental_problem

experimental_solver.run()

src.output.print_steady_results(steady_solver=experimental_solver)

src.output.draw(
    solver=experimental_solver, show_delta_pressures=True, show_streamlines=True
)
