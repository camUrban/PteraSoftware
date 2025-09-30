# NOTE: I haven't yet started refactoring this module.
"""This example script demonstrates how to automatically find the trim condition for
an unsteady simulation. It is not as well documented as some solver example scripts,
as it assumes you have read and understood those first."""

import logging

import pterasoftware as ps

# Configure a logger for this example.
example_logger = logging.getLogger("example")
example_logger.setLevel(logging.DEBUG)

# Create an airplane object. Read through the solver examples for more details on
# creating this object.
example_airplane = ps.geometry.airplane.Airplane(
    x_ref=0.14,
    weight=420,
    wings=[
        ps.geometry.wing.Wing(
            symmetric=True,
            num_chordwise_panels=5,
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=5,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
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
            ],
        ),
        ps.geometry.wing.Wing(
            x_le=5,
            symmetric=True,
            num_chordwise_panels=5,
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=5,
                    twist=-5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=1.0,
                    chord=1.0,
                    twist=-5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

# Create a movement for this example's airplane. Read through the unsteady solver
# examples for more details on this type of object.
example_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=example_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=example_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=example_airplane.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_period=1.0,
                    sweeping_amplitude=5.0,
                    heaving_period=1.0,
                    heaving_amplitude=10.0,
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=example_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=example_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=example_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

# Create an operating point object for this example's problem using the default values.
example_operating_point = ps.operating_point.OperatingPoint()

# Create an operating point movement object. Read through the unsteady solver
# examples for more details on this type of object.
example_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=example_operating_point
    )
)

# Construct this example's movement and problem object. Only calculate final results
# to speed up the solver.
example_movement = ps.movement.Movement(
    airplane_movements=[example_airplane_movement],
    operating_point_movement=example_operating_point_movement,
)
example_problem = ps.problems.UnsteadyProblem(
    movement=example_movement,
    only_final_results=False,
)

# Call the analyze_unsteady_trim function to search for a trim condition (thrust
# balances drag, weight balances lift, and all moments are close to zero) within a
# certain set of bounds.
trim_conditions = ps.trim.analyze_unsteady_trim(
    airplane_movement=example_airplane_movement,
    operating_point=example_operating_point,
    velocity_bounds=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-0.1, 0.1),
)

# Log the trim conditions. If these display "nan", then the trim function couldn't
# find a trimmed state.
example_logger.info("Trim Velocity:\t%.2f m/s" % trim_conditions[0])
example_logger.info("Trim Alpha:\t%.2f deg" % trim_conditions[1])
example_logger.info("Trim Beta:\t\t%.2f deg" % trim_conditions[2])

# The expected results are:
# Trim Velocity: 10.02 m/s
# Trim Alpha: 4.94 deg
# Trim Beta: 0.00 deg
