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
                    z_le=0.5,

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

root_wing_cross_section_movement = asmvp.performance.WingCrossSectionMovement(
    base_wing_cross_section=performance_testing_airplane.wings[0].wing_cross_sections[0],
    x_le_amplitude=0.0,
    x_le_period=0.0,
    x_le_spacing='sine',
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing='sine',
    z_le_amplitude=0.0,
    z_le_period=0.0,
    z_le_spacing='sine',
    twist_amplitude=0.0,
    twist_period=0.0,
    twist_spacing='sine',
    control_surface_deflection_amplitude=0.0,
    control_surface_deflection_period=0.0,
    control_surface_deflection_spacing='sine'
)

tip_wing_cross_section_movement = asmvp.performance.WingCrossSectionMovement(
    base_wing_cross_section=performance_testing_airplane.wings[0].wing_cross_sections[1],
    x_le_amplitude=1.0,
    x_le_period=2.5,
    x_le_spacing='sine',
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing='sine',
    z_le_amplitude=2.0,
    z_le_period=2.5,
    z_le_spacing='sine',
    twist_amplitude=15.0,
    twist_period=2.5,
    twist_spacing='sine',
    control_surface_deflection_amplitude=0.0,
    control_surface_deflection_period=0.0,
    control_surface_deflection_spacing='sine'
)

wing_movement = asmvp.performance.WingMovement(
    base_wing=performance_testing_airplane.wings[0],
    wing_cross_sections_movements=[root_wing_cross_section_movement, tip_wing_cross_section_movement],
    x_le_amplitude=0.0,
    x_le_period=0.0,
    x_le_spacing='sine',
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing='sine',
    z_le_amplitude=0.0,
    z_le_period=0.0,
    z_le_spacing='sine'
)

airplane_movement = asmvp.performance.AirplaneMovement(
    base_airplane=performance_testing_airplane,
    wing_movements=[wing_movement],
    x_ref_amplitude=0.0,
    x_ref_period=0.0,
    x_ref_spacing='sine',
    y_ref_amplitude=0.0,
    y_ref_period=0.0,
    y_ref_spacing='sine',
    z_ref_amplitude=0.0,
    z_ref_period=0.0,
    z_ref_spacing='sine'
)

# Initialize the problem's operating point.
performance_testing_operating_point = asmvp.performance.OperatingPoint()

operating_point_movement = asmvp.performance.OperatingPointMovement(
    base_operating_point=performance_testing_operating_point,
    velocity_amplitude=0.0,
    velocity_period=0.0,
    velocity_spacing='sine',
    alpha_amplitude=0.0,
    alpha_period=0.0,
    alpha_spacing='sine',
    beta_amplitude=0.0,
    beta_period=0.0,
    beta_spacing='sine'
)

movement = asmvp.performance.Movement(
    airplane_movement=airplane_movement,
    operating_point_movement=operating_point_movement,
    num_steps=100,
    delta_time=0.1
)

asmvp.output.make_flapping_gif(movement)
