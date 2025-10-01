# NOTE: I've started refactoring this module.
"""This is script is an example of how to run Ptera Software's steady horseshoe
vortex lattice method solver on a custom airplane."""

# First, import the software's main package. Note that if you wished to import this
# software into another package, you would first install the software by running "pip
# install pterasoftware" in your terminal.
import pterasoftware as ps

# Create an airplane object. Note, I am going to declare every attribute for each
# class, even most of them have usable default values. This is simply for educational
# purposes, even though it makes the code much longer than what it needs to be.
example_airplane = ps.geometry.airplane.Airplane(
    # Give the airplane object a name. This value defaults to "Untitled".
    name="Example Airplane",
    # Specify the location of the Airplane's starting point (the location of its
    # center of gravity at t=0). This is the point about which the solver will
    # calculate moments. This value is in Earth axes, relative to the simulation's
    # starting point. For the first Airplane in a simulation, this must be (0.0,
    # 0.0, 0.0) since the simulation's starting point is defined as the first
    # Airplane's starting point. The default is (0.0, 0.0, 0.0). This and every
    # input and output of this program is in SI units. Units are meters.
    Cgi_E_I=(0.0, 0.0, 0.0),
    # Define the orientation of the Airplane's body axes relative to Earth axes
    # using rotation angles with an intrinsic z-y'-x" sequence. For a standard
    # orientation aligned with Earth axes, use (0.0, 0.0, 0.0). The default is
    # (0.0, 0.0, 0.0). Units are degrees.
    angles_E_to_B_izyx=(0.0, 0.0, 0.0),
    # Specify the weight of the Airplane in Newtons. This is used by trim
    # functions. The default is 0.0.
    weight=0.0,
    # Give the reference dimensions of this aircraft. "s_ref" is the reference area
    # in meters squared, "b_ref" is the reference span in meters, and "c_ref" is the
    # reference chord in meters. I set these values to None, which is their default,
    # so that they will be populated by the first wing object's calculated
    # characteristics. Note that the reference area used in this program is the
    # wetted area of the wing's mean-camberline surface.
    s_ref=None,
    b_ref=None,
    c_ref=None,
    wings=[
        ps.geometry.wing.Wing(
            name="Main Wing",
            # Define the position of this Wing's leading edge root point (in geometry
            # axes, relative to the CG point). Geometry axes are defined with +x
            # pointing aft, +y pointing right, and +z pointing up. The default is (0.0,
            # 0.0, 0.0). Units are meters.
            prelimLer_G_Cg=(0.0, 0.0, 0.0),
            # Define the orientation of this Wing's axes relative to geometry axes using
            # rotation angles with an intrinsic z-y'-x" sequence. For a Wing aligned
            # with geometry axes, use (0.0, 0.0, 0.0). The default is (0.0, 0.0, 0.0).
            # Units are degrees.
            angles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            # Declare that this Wing is symmetric. Set to True to reflect the geometry
            # across a symmetry plane while retaining the non-reflected geometry. The
            # default is False.
            symmetric=True,
            # Set to True to reflect the geometry across a symmetry plane without
            # retaining the non-reflected geometry. If symmetric is True, mirror_only
            # must be False. The default is False.
            mirror_only=False,
            # Define the unit normal vector to the symmetry plane (in wing axes). For a
            # standard symmetric wing across the xz-plane, use (0.0, 1.0, 0.0). Must be
            # specified if symmetric or mirror_only is True.
            symmetry_normal_Wn=(0.0, 1.0, 0.0),
            # Define a point on the symmetry plane (in wing axes, relative to the
            # leading edge root point). For a standard symmetric wing, use (0.0, 0.0,
            # 0.0). Must be specified if symmetric or mirror_only is True. Units are
            # meters.
            symmetry_point_Wn_Ler=(0.0, 0.0, 0.0),
            # Define the number of chordwise panels on the wing, and the spacing
            # between them. The number of chordwise panels defaults to 8 panels. The
            # spacing defaults to "cosine", which makes the panels relatively finer,
            # in the chordwise direction, near the leading and trailing edges. The
            # other option is "uniform".
            num_chordwise_panels=6,
            chordwise_spacing="cosine",
            # Every wing has a list of wing cross sections. In order for the geometry
            # output to be sensible, each wing must have at least two wing cross
            # sections.
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    # Define the number of spanwise panels on the wing cross section.
                    # The number of spanwise panels
                    # defaults to 8 panels.
                    num_spanwise_panels=8,
                    # Set the chord of this cross section to be 1.75 meters. This
                    # value defaults to 1.0 meter.
                    chord=1.75,
                    # Define the location of the leading edge of the wing cross
                    # section relative to the wing's leading edge. These values all
                    # default to 0.0 meters.
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    # Define the orientation of this WingCrossSection's axes relative
                    # to the wing cross section parent axes using rotation angles with
                    # an intrinsic z-y'-x" sequence. The angle vector has the form
                    # (angleX, angleY, angleZ) where angleX is rotation about the
                    # x-axis, angleY is rotation about the y-axis, and angleZ is
                    # rotation about the z-axis. For the root WingCrossSection, this
                    # must be (0.0, 0.0, 0.0). The default is (0.0, 0.0, 0.0). Units
                    # are degrees.
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    # Define the type of control surface. The options are "symmetric"
                    # and "asymmetric". This is only applicable if your wing is also
                    # symmetric. If so, symmetric control surfaces will deflect in
                    # the same direction, like flaps, while asymmetric control
                    # surfaces will deflect in opposite directions, like ailerons.
                    # The default value is "symmetric".
                    control_surface_symmetry_type="symmetric",
                    # Define the point on the airfoil where the control surface
                    # hinges. This is expressed as a fraction of the chord length,
                    # back from the leading edge. The default value is 0.75.
                    control_surface_hinge_point=0.75,
                    # Define the deflection of the control surface in degrees. The
                    # default is 0.0 degrees.
                    control_surface_deflection=0.0,
                    # Define the spacing between the spanwise panels. The spacing
                    # defaults to "cosine", which makes the panels relatively finer,
                    # in the spanwise direction, near the cross section ends. The
                    # other option is "uniform".
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        # Give the airfoil a name. This defaults to "Untitled
                        # Airfoil". This name should correspond to a name in the
                        # airfoil directory or a NACA four series airfoil, unless you
                        # are passing in your own coordinates.
                        name="naca2412",
                        # If you wish to pass in coordinates, set this to an N x 2
                        # array of the airfoil's coordinates, where N is the number
                        # of coordinates. Treat this as an immutable, don't edit
                        # directly after initialization. If you wish to load
                        # coordinates from the airfoil directory, leave this as None.
                        # The default is None. Make sure that any airfoil coordinates
                        # used range in x from 0 to 1.
                        outline_A_lp=None,
                        # This is the variable that determines whether you would like
                        # to repanel the airfoil coordinates. This applies to
                        # coordinates passed in by the user or to the directory
                        # coordinates. I highly recommended setting this to True. The
                        # default is True.
                        resample=True,
                        # This is number of points to use if repaneling the airfoil.
                        # It is ignored if the repanel is False. The default is 400.
                        n_points_per_side=400,
                    ),
                ),
                # Define the tip WingCrossSection. This cross section defines the
                # wing's tip geometry. The num_spanwise_panels is None because this
                # is the tip cross section.
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.75, 6.0, 1.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
        # Define the next wing.
        ps.geometry.wing.Wing(
            name="V-Tail",
            prelimLer_G_Cg=(6.75, 0.0, 0.25),
            angles_G_to_prelimWn_izyx=(0.0, 5.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetry_normal_Wn=(0.0, 1.0, 0.0),
            symmetry_point_Wn_Ler=(0.0, 0.0, 0.0),
            num_chordwise_panels=8,
            chordwise_spacing="cosine",
            # Define this wing's root wing cross section.
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    # Give the root wing cross section an airfoil.
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                # Define the wing's tip wing cross section.
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.5, 2.0, 1.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

# Define a new operating point object. This defines the state at which the airplane
# object is operating.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=5.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define a new steady problem. A steady problem contains an airplane object and an
# operating point object.
example_problem = ps.problems.SteadyProblem(
    # Set this steady problem's airplane object to be the one we just created.
    airplanes=[example_airplane],
    # Set this steady problem's operating point object to be the one we just created.
    operating_point=example_operating_point,
)

# Now, the airplane and operating point object exist within the steady problem
# object. I like to delete the external pointers to these objects to ease debugging.
del example_airplane
del example_operating_point

# Define a new solver. The available solver objects are the steady horseshoe vortex
# lattice method solver, the steady ring vortex lattice method solver, and the
# unsteady ring vortex lattice method solver.
example_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
    # Solvers just take in one attribute: the problem they are going to solve.
    steady_problem=example_problem
)

# Delete the extraneous pointer to the problem as it is now contained within the
# solver. Again, this is unnecessary, I just like to do this to ease debugging.
del example_problem

# Run the example solver.
example_solver.run(
    # This parameter determines the detail of information that the solver's logger
    # will output while running. The options are, in order of detail and severity,
    # "Debug", "Info", "Warning", "Error", "Critical". The default value is "Warning".
    logging_level="Warning",
)

# Call this function from the output module to print the results.
ps.output.print_steady_results(steady_solver=example_solver)

# Call the software's draw function on the solver.
ps.output.draw(
    solver=example_solver,
    # Tell the draw function to color the aircraft's wing panels with the local lift
    # coefficient. The valid arguments for this parameter are None, "induced drag",
    # "side force", or "lift".
    scalar_type="lift",
    # Tell the draw function to show the calculated streamlines. This value defaults
    # to false.
    show_streamlines=True,
    # Tell the draw function to not show any wake vortices. As this is a steady
    # solver, no vortices have been shed into the wake. This value defaults to false.
    show_wake_vortices=False,
    # Tell the draw function to not save the drawing as an image file. This way,
    # the drawing will still be displayed but not saved. This value defaults to false.
    save=False,
)

# Compare the output you see with the expected outputs saved in the "docs/examples
# expected output" directory.
