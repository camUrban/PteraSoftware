# NOTE: I haven't yet started refactoring this module.
"""This is script is an example of how to run Ptera Software's unsteady ring vortex
lattice method solver on three airplanes, flying in formation, each with custom
geometry and motion. Note, I will comment this example less rigorously than the
single-airplane examples for readability. I recommend you read and understand those
examples before reading this example."""

import pterasoftware as ps

x_spacing = 13
y_spacing = 13

# Create the lead airplane object.
lead_airplane = ps.geometry.airplane.Airplane(
    name="Lead Airplane",
    # Specify the location of the lead airplane's center of gravity. This is the
    # point around about which the solver will calculate the moments on the airplane.
    # These three values default to 0.0 meters. Note: these values are relative to
    # the global coordinate system fixed front left corner of the first airplane's
    # first wing's root wing cross section.
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    # Define the location of the leading edge of the wing cross
                    # section relative to the wing's leading edge. These values all
                    # default to 0.0 meters.
                    x_le=0.0,
                    y_le=0.0,
                    # Assign the twist of this wing cross section. Note: when
                    # assigning angles of attack to multiple airplanes, it is better
                    # to set the operating point's angle of attack to zero, and then
                    # use offset the twist values of all the wing cross sections to
                    # simulate each aircraft having an angle of attack.
                    twist=5.0,
                    chord=1.75,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.75,
                    y_le=6.0,
                    chord=1.5,
                    twist=5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
)

# Now define the lead airplane's movement object.
lead_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=lead_airplane,
    wing_movements=[  # Define the main wing's movement.
        ps.movements.wing_movement.WingMovement(
            base_wing=lead_airplane.wings[0],
            wing_cross_section_movements=[
                # Define the root wing cross section's movement object.
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                # Define the tip wing cross section's movement object.
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[
                        1
                    ],
                    sweeping_amplitude=15.0,
                    sweeping_period=1.5,
                    sweeping_spacing="sine",
                ),
            ],
        ),
    ],
)

# Create the trailing right airplane object.
right_airplane = ps.geometry.airplane.Airplane(
    name="Right Airplane",
    # Specify the location of the right airplane's center of gravity. This is the
    # point around about which the solver will calculate the moments on the airplane.
    # These three values default to 0.0 meters. Note: these values are relative to
    # the global coordinate system fixed front left corner of the first airplane's
    # first wing's root wing cross section.
    x_ref=x_spacing,
    y_ref=y_spacing,
    z_ref=0.0,
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    twist=5.0,
                    chord=1.75,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.75,
                    y_le=6.0,
                    chord=1.5,
                    twist=5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
)

# Now define the trailing right airplane's movement object.
right_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=right_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=right_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=right_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=right_airplane.wings[0].wing_cross_sections[
                        1
                    ],
                    sweeping_amplitude=15.0,
                    sweeping_period=1.5,
                    sweeping_spacing="sine",
                ),
            ],
        ),
    ],
)

# Create the trailing left airplane object.
left_airplane = ps.geometry.airplane.Airplane(
    name="Left Airplane",
    # Specify the location of the left airplane's center of gravity. This is the
    # point around about which the solver will calculate the moments on the airplane.
    # These three values default to 0.0 meters. Note: these values are relative to
    # the global coordinate system fixed front left corner of the first airplane's
    # first wing's root wing cross section.
    x_ref=x_spacing,
    y_ref=-y_spacing,
    z_ref=0.0,
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    twist=5.0,
                    chord=1.75,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.75,
                    y_le=6.0,
                    chord=1.5,
                    twist=5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
)

# Now define the trailing left airplane's movement object.
left_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=left_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=left_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=left_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=left_airplane.wings[0].wing_cross_sections[
                        1
                    ],
                    sweeping_amplitude=15.0,
                    sweeping_period=1.5,
                    sweeping_spacing="sine",
                ),
            ],
        ),
    ],
)

# Define a new operating point object. This defines the state at which all the
# airplanes objects are operating. Note: when assigning angles of attack to multiple
# airplanes, it is better to set the operating point's angle of attack to zero,
# and then use offset the twist values of all the wing cross sections to simulate
# each aircraft having an angle of attack.
operating_point = ps.operating_point.OperatingPoint(vCg__E=10.0, alpha=0.0)

# Define the operating point's movement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=operating_point
)

# Delete the extraneous airplane and operating point objects, as these are now
# contained within their respective movement objects.
del lead_airplane
del right_airplane
del left_airplane
del operating_point

# Define the movement object. This contains each airplane's movement and the operating
# point movement.
movement = ps.movement.Movement(
    airplane_movements=[
        lead_airplane_movement,
        right_airplane_movement,
        left_airplane_movement,
    ],
    operating_point_movement=operating_point_movement,
    num_cycles=2,
)

# Delete the extraneous airplane and operating point movement objects, as these are
# now contained within the movement object.
del lead_airplane_movement
del right_airplane_movement
del left_airplane_movement
del operating_point_movement

# Define the unsteady example problem.
problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Define a new solver. The available solver objects are the steady horseshoe vortex
# lattice method solver, the steady ring vortex lattice method solver, and the
# unsteady ring vortex lattice method solver.
solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    # Solvers just take in one attribute: the problem they are going to solve.
    unsteady_problem=problem,
)

# Delete the extraneous pointer to the problem as it is now contained within the
# solver.
del problem

# Run the example solver.
solver.run(
    prescribed_wake=False,
)

# Call the software's animate function on the solver. This produces a WebP animation
# of the wake being shed. The animation is saved in the same directory as this
# script. Press "q", after orienting the view, to begin the animation.
ps.output.animate(
    unsteady_solver=solver,
    scalar_type="lift",
    show_wake_vortices=True,
    # Tell the animate function to not save the animation as file. This way,
    # the animation will still be displayed but not saved. This value defaults to
    # false.
    save=True,
)

# Compare the output you see with the expected outputs saved in the "docs/examples
# expected output" directory.
