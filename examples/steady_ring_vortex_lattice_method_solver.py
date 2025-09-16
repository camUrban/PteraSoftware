"""This is script is an example of how to run Ptera Software's steady ring vortex
lattice method solver on a custom airplane.
"""

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
    # Specify the location of the airplane's center of gravity. This is the point
    # around about which the solver will calculate the moments on the airplane. These
    # three values default to 0.0 meters. This and every input and output of this
    # program is in SI units. Note: these values are relative to the global
    # coordinate system fixed front left corner of the first airplane's first wing's
    # root wing cross section.
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
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
            # Define the location of the leading edge of the wing relative to the
            # global coordinate system fixed front left corner of the first
            # airplane's first wing's root wing cross section. These values all
            # default to 0.0 meters.
            x_le=0.0,
            y_le=0.0,
            z_le=0.0,
            # Declare that this wing is symmetric. This means that the geometry will
            # be reflected across plane of this wing's root wing cross section. Note
            # that the geometry coordinates are defined as such: If you were riding
            # in the airplane, the positive x direction would point behind you,
            # the positive y direction would point out of your right wing, and the
            # positive z direction would point upwards, out of your chair. These
            # directions form a right-handed coordinate system. The default value of
            # "symmetric" is false.
            symmetric=True,
            # Define the number of chordwise panels on the wing, and the spacing
            # between them. The number of chordwise panels defaults to 8 panels. The
            # spacing defaults to "cosine", which makes the panels relatively finer,
            # in the chordwise direction, near the leading and trailing edges. The
            # other option is "uniform".
            num_chordwise_panels=8,
            chordwise_spacing="cosine",
            # Every wing has a list of wing cross sections. In order for the geometry
            # output to be sensible, each wing must have at least two wing cross
            # sections.
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    # Define the location of the leading edge of the wing cross
                    # section relative to the wing's leading edge. These values all
                    # default to 0.0 meters.
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    # Define the twist of the wing cross section in degrees. This is
                    # equivalent to incidence angle of cross section. The twist is
                    # about the leading edge. Note that the twist is only stable up
                    # to 45.0 degrees. Values above that produce unexpected results.
                    # This will be fixed in a future release. The default value is
                    # 0.0 degrees. Positive twist corresponds to positive rotation
                    # about the y axis, as defined by the right-hand rule.
                    twist=0.0,  # Define the type of control surface. The options are "symmetric"
                    # and "asymmetric". This is only applicable if your wing is also
                    # symmetric. If so, symmetric control surfaces will deflect in
                    # the same direction, like flaps, while asymmetric control
                    # surfaces will deflect in opposite directions, like ailerons.
                    # The default value is "symmetric".
                    control_surface_type="asymmetric",
                    # Define the point on the airfoil where the control surface
                    # hinges. This is expressed as a fraction of the chord length,
                    # back from the leading edge. The default value is 0.75.
                    control_surface_hinge_point=0.75,
                    # Define the deflection of the control surface in degrees. The
                    # default is 0.0 degrees. We'll set it to 10.0 degrees to show an
                    # example of an aileron deflection.
                    control_surface_deflection=10.0,
                    # Define the number of spanwise panels on the wing cross section,
                    # and the spacing between them. The number of spanwise panels
                    # defaults to 8 panels. The spacing defaults to "cosine",
                    # which makes the panels relatively finer, in the spanwise
                    # direction, near the cross section ends. The other option is
                    # "uniform".
                    num_spanwise_panels=8,
                    spanwise_spacing="cosine",
                    # Set the chord of this cross section to be 1.75 meters. This
                    # value defaults to 1.0 meter.
                    chord=1.5,
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
                        coordinates=None,
                        # This is the variable that determines whether you would like
                        # to repanel the airfoil coordinates. This applies to
                        # coordinates passed in by the user or to the directory
                        # coordinates. I highly recommended setting this to True. The
                        # default is True.
                        repanel=True,
                        # This is number of points to use if repaneling the airfoil.
                        # It is ignored if the repanel is False. The default is 400.
                        n_points_per_side=400,
                    ),
                ),
                # Define the next wing cross section. From here on out,
                # the declarations will not be as commented as the previous. See the
                # above comments if you have questions.
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=1.5,
                    y_le=6.0,
                    z_le=0.5,
                    chord=0.75,
                    control_surface_type="asymmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=10.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
        # Define the next wing.
        ps.geometry.wing.Wing(
            name="Horizontal Stabilizer",
            x_le=6.75,
            z_le=0.25,
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    chord=1.5,
                    # Give the root wing cross section an airfoil.
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    twist=-5.0,
                ),
                # Define the wing's tip wing cross section.
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.5,
                    y_le=2.0,
                    chord=1.0,
                    twist=-5.0,
                    # Give the tip wing cross section an airfoil.
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
        # Define the next wing.
        ps.geometry.wing.Wing(
            name="Vertical Stabilizer",
            x_le=6.75,
            z_le=0.5,
            symmetric=False,
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    chord=1.5,
                    # Give the root wing cross section an airfoil.
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                # Define the wing's tip wing cross section.
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.5,
                    z_le=2.0,
                    chord=1.0,
                    # Give the tip wing cross section an airfoil.
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
    # Define the density of the fluid the airplane is flying in. This defaults to
    # 1.225 kilograms per meters cubed.
    density=1.225,
    # Define the angle of sideslip the airplane is experiencing. This defaults to 0.0
    # degrees.
    beta=0.0,
    # Define the freestream velocity at which the airplane is flying. This defaults
    # to 10.0 meters per second.
    velocity=10.0,
    # Define the angle of attack the airplane is experiencing. This defaults to 5.0
    # degrees.
    alpha=1.0,
)

# Define a new steady problem. A steady problem contains an airplane object and an
# operating point object.
example_problem = ps.problems.SteadyProblem(
    # Set this steady problem's list of airplane objects to be the one we just created.
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
example_solver = ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
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
    # Tell the draw function to color the aircraft's wing panels with the local
    # lift coefficient. The valid arguments for this parameter are None, "induced drag",
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
