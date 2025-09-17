"""This script is an example of analyzing the steady convergence of a problem with
multiple airplanes. It should take a few minutes to run. It will display the
convergence progress and results in the console."""

import pterasoftware as ps

# Create two airplane objects. Read through the solver and formation examples for
# more details on creating these objects.
leading_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.25,
                    y_le=10.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
        ps.geometry.wing.Wing(
            x_le=5.0,
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)
trailing_airplane = ps.geometry.airplane.Airplane(
    x_ref=10,
    y_ref=-5,
    wings=[
        ps.geometry.wing.Wing(
            x_le=10,
            y_le=-5,
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.25,
                    y_le=10.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
        ps.geometry.wing.Wing(
            x_le=15.0,
            y_le=-5,
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

# Create an operating point object.
operating_point = ps.operating_point.OperatingPoint()

# Create a steady problem. We will pass this into the convergence function.
problem = ps.problems.SteadyProblem(
    airplanes=[leading_airplane, trailing_airplane],
    operating_point=operating_point,
)

del leading_airplane
del trailing_airplane
del operating_point

# Run the steady convergence analysis. This will run the problem several times,
# modifying average panel aspect ratio, and number of chordwise panels with each
# iteration. Once it detects that the net load coefficients haven't change by more
# than the convergence criteria (measured as an absolute percent error), it will
# return the parameters it found to result in a converged solution. See the
# analyze_steady_convergence function docstring for more details. The progress and
# results are displayed to the console.
ps.convergence.analyze_steady_convergence(
    ref_problem=problem,
    solver_type="steady ring vortex lattice method",
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 8),
    convergence_criteria=1.0,
)

# Check the console that the convergence analysis found that the solution converged
# with the following parameters:
# Panel aspect ratio: 4
# Chordwise panels: 4
