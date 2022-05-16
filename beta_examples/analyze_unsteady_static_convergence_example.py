# ToDo: Document this script.

import logging

import pterasoftware as ps

# Configure a logger for this example.
example_logger = logging.getLogger("example")
example_logger.setLevel(logging.DEBUG)

# Create an airplane object. Read through the solver examples for more details on
# creating this object.
leading_airplane = ps.geometry.Airplane(
    weight=250,
    wings=[
        ps.geometry.Wing(
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca2412",
                    ),
                    spanwise_spacing="uniform",
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=3.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
    ],
)

leading_airplane_movement = ps.movement.AirplaneMovement(
    base_airplane=leading_airplane,
    wing_movements=[
        ps.movement.WingMovement(
            base_wing=leading_airplane.wings[0],
            wing_cross_sections_movements=[
                ps.movement.WingCrossSectionMovement(
                    base_wing_cross_section=leading_airplane.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                ps.movement.WingCrossSectionMovement(
                    base_wing_cross_section=leading_airplane.wings[
                        0
                    ].wing_cross_sections[1],
                ),
            ],
        )
    ],
)

operating_point = ps.operating_point.OperatingPoint()
operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=operating_point
)

movement = ps.movement.Movement(
    airplane_movements=[leading_airplane_movement],
    operating_point_movement=operating_point_movement,
)

del leading_airplane_movement
del operating_point_movement

converged_attributes = ps.convergence.analyze_unsteady_convergence(
    ref_movement=movement,
)
