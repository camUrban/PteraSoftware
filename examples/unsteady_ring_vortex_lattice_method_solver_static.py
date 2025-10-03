# NOTE: I haven't yet started refactoring this module.
"""This is script is an example of how to run Ptera Software's unsteady ring vortex
lattice method solver on a custom airplane with static geometry."""

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
    c_ref=None,  # All airplane objects have a list of wings.
    wings=[  # Create the first wing object in this airplane.
        ps.geometry.wing.Wing(
            wing_cross_sections=[  # Create a new wing cross section object.
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
                    twist=0.0,
                    # Define the type of control surface. The options are "symmetric"
                    # and "asymmetric". This is only applicable if your wing is also
                    # symmetric. If so, symmetric control surfaces will deflect in
                    # the same direction, like flaps, while asymmetric control
                    # surfaces will deflect in opposite directions, like ailerons.
                    # The default value is "symmetric".
                    control_surface_type="symmetric",
                    # Define the point on the airfoil where the control surface
                    # hinges. This is expressed as a fraction of the chord length,
                    # back from the leading edge. The default value is 0.75.
                    control_surface_hinge_point=0.75,
                    # Define the deflection of the control surface in degrees. The
                    # default is 0.0 degrees.
                    control_surface_deflection=0.0,
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
                    chord=1.75,  # Every wing cross section has an airfoil object.
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
                    x_le=0.75,
                    y_le=6.0,
                    z_le=1.0,
                    chord=1.5,
                    twist=5.0,
                    # Give this wing cross section an airfoil.
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
        ),
        # Define the next wing.
        ps.geometry.wing.Wing(
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
                    z_le=1.0,
                    chord=1.0,
                    twist=-5.0,
                    # Give the tip wing cross section an airfoil.
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
            name="V-Tail",
            symmetric=True,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
        ),
    ],
)

# Now define the main wing's root wing cross section's movement. Cross sections can
# move in three ways: sweeping, pitching, and heaving. Sweeping is defined as the
# relative rotation of this wing cross section's leading edge to its preceding wing
# cross section's leading edge about the airplane's body x axis. Pitching is defined
# as the relative rotation of this wing cross section's leading edge to the preceding
# wing cross section's leading edge about the body y axis. Heaving is defined as the
# relative rotation of this wing cross section's leading edge to the preceding wing
# cross section's leading edge about the body z axis. The sign of all rotations is
# determined via the right-hand-rule.
main_wing_root_wing_cross_section_movement = ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
    # Provide the base cross section.
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
    # Define the sweeping amplitude. This value is in degrees. As this is the first
    # wing cross section, this must be 0.0 degrees, which is the default value.
    sweeping_amplitude=0.0,
    # Define the sweeping period. This value is in seconds. As this is the first wing
    # cross section, this must be 0.0 seconds, which is the default value.
    sweeping_period=0.0,
    # Define the time step spacing of the sweeping. This is "sine" by default. The
    # options are "sine" and "uniform".
    sweeping_spacing="sine",
    # Define the pitching amplitude. This value is in degrees. As this is the first
    # wing cross section, this must be 0.0 degrees, which is the default value.
    pitching_amplitude=0.0,
    # Define the pitching period. This value is in seconds. As this is the first wing
    # cross section, this must be 0.0 seconds, which is the default value.
    pitching_period=0.0,
    # Define the time step spacing of the pitching. This is "sine" by default. The
    # options are "sine" and "uniform".
    pitching_spacing="sine",
    # Define the heaving amplitude. This value is in degrees. As this is the first
    # wing cross section, this must be 0.0 degrees, which is the default value.
    heaving_amplitude=0.0,
    # Define the heaving period. This value is in seconds. As this is the first wing
    # cross section, this must be 0.0 seconds, which is the default value.
    heaving_period=0.0,
    # Define the time step spacing of the heaving. This is "sine" by default. The
    # options are "sine" and "uniform".
    heaving_spacing="sine",
)

# Define the main wing's tip wing cross section's movement. As the example has static
# geometry, the movement attributes can be excluded, and the default values will
# suffice.
main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
    )
)

# Define the v-tail's root wing cross section's movement. As the example has static
# geometry, the movement attributes can be excluded, and the default values will
# suffice.
v_tail_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[0],
    )
)

# Define the v-tail's tip wing cross section's movement. As the example has static
# geometry, the movement attributes can be excluded, and the default values will
# suffice.
v_tail_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[1],
    )
)

# Now define the main wing's movement. In addition to their wing cross sections'
# relative movements, wings' leading edge positions can move as well.
main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_section_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    x_le_amplitude=0.0,
    x_le_period=0.0,
    x_le_spacing="sine",
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing="sine",
    z_le_amplitude=0.0,
    z_le_period=0.0,
    z_le_spacing="sine",
)

# Delete the extraneous wing cross section movement objects, as these are now
# contained within the wing movement object. This is unnecessary, but it can make
# debugging easier.
del main_wing_root_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement

# Make the v-tail's wing movement object.
v_tail_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[1],
    wing_cross_section_movements=[
        v_tail_root_wing_cross_section_movement,
        v_tail_tip_wing_cross_section_movement,
    ],
)

# Delete the extraneous wing cross section movement objects, as these are now
# contained within the wing movement object. This is unnecessary, but it can make
# debugging easier.
del v_tail_root_wing_cross_section_movement
del v_tail_tip_wing_cross_section_movement

# Now define the airplane's movement object. In addition to their wing's and wing
# cross sections' relative movements, airplane's reference positions can move as well.
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(  # Define the base airplane object.
    base_airplane=example_airplane,  # Add the list of wing movement objects.
    wing_movements=[main_wing_movement, v_tail_movement],
    # Define the amplitude of the reference position's change in x position. This
    # value is in meters. This is set to 0.0 meters, which is the default value.
    x_ref_amplitude=0.0,
    # Define the period of the reference position's change in x position. This value
    # is in seconds. This is set to 0.0 seconds, which is the default value.
    x_ref_period=0.0,
    # Define the time step spacing of the reference position's change in x position.
    # This is "sine" by default. The options are "sine" and "uniform".
    x_ref_spacing="sine",
    # Define the amplitude of the reference position's change in y position. This
    # value is in meters. This is set to 0.0 meters, which is the default value.
    y_ref_amplitude=0.0,
    # Define the period of the reference position's change in y position. This value
    # is in seconds. This is set to 0.0 seconds, which is the default value.
    y_ref_period=0.0,
    # Define the time step spacing of the reference position's change in y position.
    # This is "sine" by default. The options are "sine" and "uniform".
    y_ref_spacing="sine",
    # Define the amplitude of the reference position's change in z position. This
    # value is in meters. This is set to 0.0 meters, which is the default value.
    z_ref_amplitude=0.0,
    # Define the period of the reference position's change in z position. This value
    # is in seconds. This is set to 0.0 seconds, which is the default value.
    z_ref_period=0.0,
    # Define the time step spacing of the reference position's change in z position.
    # This is "sine" by default. The options are "sine" and "uniform".
    z_ref_spacing="sine",
)

# Delete the extraneous wing movement objects, as these are now contained within the
# airplane movement object.
del main_wing_movement
del v_tail_movement

# Define a new operating point object. This defines the state at which the airplane
# object is operating.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=1.0, beta=0.0, nu=15.06e-6
)

# Define the operating point's movement. The operating point's velocity can change
# with respect to time.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point, periodVCg__E=0.0, spacingVCg__E="sine"
)

# Define the movement object. This contains the airplane movement and the operating
# point movement.
movement = ps.movement.Movement(  # Add the airplane movement.
    airplane_movements=[airplane_movement],  # Add the operating point movement.
    operating_point_movement=operating_point_movement,
    # Leave the number of time steps and the length of each time step unspecified.
    # The solver will automatically set the length of the time steps so that the wake
    # ring vortices and the bound ring vortices have approximately the same area. The
    # solver will also determine if the geometry is static or not. If it is static,
    # the number of steps will be set such that the wake extends ten chord lengths
    # back from the main wing. If the geometry isn't static, the number of steps will
    # be set such that three periods of the slowest movement oscillation complete.
    num_steps=None,
    delta_time=None,
)

# Delete the extraneous airplane and operating point movement objects, as these are
# now contained within the movement object.
del airplane_movement
del operating_point_movement

# Define the unsteady example problem.
example_problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Define a new solver. The available solver objects are the steady horseshoe vortex
# lattice method solver, the steady ring vortex lattice method solver, and the
# unsteady ring vortex lattice method solver.
example_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    # Solvers just take in one attribute: the problem they are going to solve.
    unsteady_problem=example_problem,
)

# Delete the extraneous pointer to the problem as it is now contained within the
# solver.
del example_problem

# Run the example solver.
example_solver.run(
    # This parameter determines the detail of information that the solver's logger
    # will output while running. The options are, in order of detail and severity,
    # "Debug", "Info", "Warning", "Error", "Critical". The default value is "Warning".
    logging_level="Warning",
    # Use a prescribed wake model. This is faster, but may be slightly less accurate.
    prescribed_wake=True,
)

# Call the software's draw function on the solver. Press "q" to close the plotter
# after it draws the output.
ps.output.draw(  # Set the solver to the one we just ran.
    solver=example_solver,
    # Tell the draw function to color the aircraft's wing panels with the local lift
    # coefficient. The valid arguments for this parameter are None, "induced drag",
    # "side force", or "lift".
    scalar_type="lift",
    # Tell the draw function to show the calculated streamlines. This value defaults
    # to False.
    show_streamlines=True,
    # Tell the draw function to not show the wake vortices. This value defaults to
    # False.
    show_wake_vortices=False,
    # Tell the draw function to not save the drawing as an image file. This way,
    # the drawing will still be displayed but not saved. This value defaults to False.
    save=False,
)

# Call the software's animate function on the solver. This produces a GIF of the wake
# being shed. The GIF is saved in the same directory as this script. Press "q",
# after orienting the view, to begin the animation.
ps.output.animate(  # Set the unsteady solver to the one we just ran.
    unsteady_solver=example_solver,
    # Tell the animate function to color the aircraft's wing panels with the local
    # lift coefficient. The valid arguments for this parameter are None, "induced drag",
    # "side force", or "lift".
    scalar_type="lift",
    # Tell the animate function to show the wake vortices. This value defaults to
    # False.
    show_wake_vortices=True,
    # Tell the animate function to not save the animation as file. This way,
    # the animation will still be displayed but not saved. This value defaults to
    # False.
    save=False,
)

# Call the software's plotting function on the solver. This produces graphs of the
# output forces and moments with respect to time.
ps.output.plot_results_versus_time(  # Set the unsteady solver to the one we just ran.
    unsteady_solver=example_solver,
    # Set the show attribute to True, which is the default value. With this set to
    # show, some IDEs (such as PyCharm in "Scientific Mode") will display the plots
    # in a sidebar. Other IDEs may not display the plots, in which case you should
    # set the save attribute to True, and open the files after they've been saved to
    # the current directory.
    show=True,
    save=False,
)

# Compare the output you see with the expected outputs saved in the "docs/examples
# expected output" directory.
