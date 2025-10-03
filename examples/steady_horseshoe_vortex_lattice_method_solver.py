"""This is script is an example of how to run Ptera Software's
SteadyHorseshoeVortexLatticeMethodSolver with a custom Airplane."""

# First, import the software's main package. Note that if you wished to import this
# software into another package, you would first install it by running "pip install
# pterasoftware" in your terminal.
import pterasoftware as ps

# Create an Airplane with our custom geometry. I am going to declare every parameter
# for Airplane, even though most of them have usable default values. This is for
# educational purposes, but keep in mind that it makes the code much longer than it
# needs to be. For details about each parameter, read the detailed class docstring.
# The same caveats apply to the other classes, methods, and functions I call in this
# script.
example_airplane = ps.geometry.airplane.Airplane(
    name="Example Airplane",
    Cgi_E_I=(0.0, 0.0, 0.0),
    angles_E_to_B_izyx=(0.0, 0.0, 0.0),
    weight=0.0,
    s_ref=None,
    b_ref=None,
    c_ref=None,
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=1.75,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.75, 6.0, 1.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="Main Wing",
            prelimLer_G_Cg=(0.0, 0.0, 0.0),
            angles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetry_normal_G=(0.0, 1.0, 0.0),
            symmetry_point_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=6,
            chordwise_spacing="cosine",
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.5, 2.0, 1.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="V-Tail",
            prelimLer_G_Cg=(6.75, 0.0, 0.25),
            angles_G_to_prelimWn_izyx=(0.0, 5.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetry_normal_G=(0.0, 1.0, 0.0),
            symmetry_point_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=8,
            chordwise_spacing="cosine",
        ),
    ],
)

# Define a new OperatingPoint, which we'll pass into the SteadyProblem.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=5.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define a new SteadyProblem, which contains the OperatingPoint and a list of one or
# more Airplanes.
example_problem = ps.problems.SteadyProblem(
    airplanes=[example_airplane],
    operating_point=example_operating_point,
)

# Now that the Airplane and OperatingPoints exist within the SteadyProblem, I like to
# delete the external pointers to these objects to ease debugging.
del example_airplane
del example_operating_point

# Define a new solver. The available solver classes are
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver,
# and UnsteadyRingVortexLatticeMethodSolver. We'll create a
# SteadyHorseshoeVortexLatticeMethodSolver, which requires a SteadyProblem.
example_solver = (
    ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=example_problem
    )
)

del example_problem

# Run the solver.
example_solver.run(
    logging_level="Warning",
)

# Call this function from the output module to print the results.
ps.output.print_results(example_solver)

# Call the output module's draw function on the solver.
ps.output.draw(
    solver=example_solver,
    scalar_type="lift",
    show_streamlines=True,
    show_wake_vortices=False,
    save=False,
    testing=False,
)
