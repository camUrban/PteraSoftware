"""This script is an example of analyzing the convergence of a SteadyProblem with
multiple Airplanes. It should take a few minutes to run. It will display the
convergence progress and results in the console."""

import pterasoftware as ps

# Configure logging to display info level messages. This is important for seeing the
# output from the convergence function.
ps.set_up_logging(level="Info")

# Create two Airplanes. Read through the solver and formation examples for
# more details on creating these Airplanes.
leading_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=8,
                    control_surface_symmetry_type="asymmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=8,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.25, 10.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            name="Main Wing",
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            chordwise_spacing="uniform",
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=8,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            name="Tail",
            Ler_Gs_Cgs=(10, 0, 0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            chordwise_spacing="uniform",
        ),
    ],
    name="Leading Airplane",
)

trailing_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=8,
                    control_surface_symmetry_type="asymmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=8,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.25, 10.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            name="Main Wing",
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            chordwise_spacing="uniform",
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=8,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            name="Tail",
            Ler_Gs_Cgs=(10, 0, 0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            chordwise_spacing="uniform",
        ),
    ],
    name="Trailing Airplane",
    Cg_GP1_CgP1=(-20, 5, 0),
)

# Create an OperatingPoint.
operating_point = ps.operating_point.OperatingPoint()

# Create a SteadyProblem. We will pass this into the convergence function.
steady_problem = ps.problems.SteadyProblem(
    airplanes=[leading_airplane, trailing_airplane],
    operating_point=operating_point,
)

del leading_airplane
del trailing_airplane
del operating_point

# Run the steady convergence analysis. This will run several simulations, modifying
# average Panel aspect ratio and number of chordwise Panels with each iteration. Once
# it detects that the net load coefficients haven't change by more than the
# convergence criteria (measured as an absolute percent error), it will return the
# parameters it found to result in a converged solution. See the
# analyze_steady_convergence function docstring for more details. The progress and
# results are displayed to the console.
ps.convergence.analyze_steady_convergence(
    ref_problem=steady_problem,
    solver_type="steady ring vortex lattice method",
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 8),
    convergence_criteria=1.0,
)

# Check the console that the convergence analysis found the following converged
# parameters:
#   Panel aspect ratio: 4
#   Chordwise Panels: 5
