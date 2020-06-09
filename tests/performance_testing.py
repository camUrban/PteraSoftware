"""

"""

import aviansoftwareminimumviableproduct as asmvp

# Initialize the airplane.
performance_testing_airplane = asmvp.geometry.Airplane(

    # Name the airplane.
    name="Performance Testing Airplane",

    # Define a list of the airplane's wings.
    wings=[

        # Initialize the wing object.
        asmvp.geometry.Wing(

            # Name the wing.
            name="Wing",

            # This will be a symmetrical wing.
            symmetric=True,

            # Define a list of the wing's cross sections.
            wing_cross_sections=[

                # Initialize the root cross section object.
                asmvp.geometry.WingCrossSection(
                    # Define the cross section's twist chord.
                    chord=1.0,
                    twist=0,
                    # Initialize this cross section's airfoil object.
                    airfoil=asmvp.geometry.Airfoil(name="naca5512"),
                ),

                # Initialize the tip cross section object.
                asmvp.geometry.WingCrossSection(

                    # Define the cross section's leading edge placement.
                    x_le=1.0,
                    y_le=5.0,
                    z_le=0.0,

                    # Define the cross section's twist and chord.
                    twist=5,
                    chord=0.75,

                    # Initialize this cross section's airfoil object.
                    airfoil=asmvp.geometry.Airfoil(name="naca5512"),
                )
            ]
        )
    ]
)

root_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
    base_wing_cross_section=performance_testing_airplane.wings[0].wing_cross_sections[0]
)

tip_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
    base_wing_cross_section=performance_testing_airplane.wings[0].wing_cross_sections[1],
    x_le_amplitude=2.5,
    x_le_period=5.0,
    x_le_spacing='sine',
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing='sine',
    z_le_amplitude=2.0,
    z_le_period=2.5,
    z_le_spacing='sine',
    twist_amplitude=30.0,
    twist_period=5.0,
    twist_spacing='sine'
)

wing_movement = asmvp.movement.WingMovement(
    base_wing=performance_testing_airplane.wings[0],
    wing_cross_sections_movements=[root_wing_cross_section_movement, tip_wing_cross_section_movement]
)

airplane_movement = asmvp.movement.AirplaneMovement(
    base_airplane=performance_testing_airplane,
    wing_movements=[wing_movement]
)

performance_testing_operating_point = asmvp.operating_point.OperatingPoint()

operating_point_movement = asmvp.movement.OperatingPointMovement(
    base_operating_point=performance_testing_operating_point
)

movement = asmvp.movement.Movement(
    airplane_movement=airplane_movement,
    operating_point_movement=operating_point_movement,
)

asmvp.output.make_flapping_gif(movement)
