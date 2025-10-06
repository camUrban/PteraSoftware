# NOTE: I haven't yet started refactoring this module.
"""This example script demonstrates how to automatically find the trim condition for
a steady simulation. It is not as well documented as some solver example scripts,
as it assumes you have read and understood those first."""

import logging

import pterasoftware as ps

# Configure a logger for this example.
example_logger = logging.getLogger("example")
example_logger.setLevel(logging.DEBUG)

# Create an airplane object. Read through the solver examples for more details on
# creating this object.
default_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
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
            symmetric=True,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    chord=0.5,
                    twist=-5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=1.0,
                    chord=0.5,
                    twist=-5.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
            symmetric=True,
        ),
    ],
    weight=250,
)

# Create an operating point object for this example's problem. Be sure to specify an
# external thrust because this aircraft is not flapping, and will therefore generate
# no thrust of its own. This external thrust could be due to a propeller or other
# type of engine.
default_operating_point = ps.operating_point.OperatingPoint(externalFX_W=5)

# Construct this example's problem object.
default_problem = ps.problems.SteadyProblem(
    airplanes=[default_airplane], operating_point=default_operating_point
)

# Call the analyze_steady_trim function to search for a trim condition (thrust
# balances drag, weight balances lift, and all moments are close to zero) within a
# certain set of bounds.
trim_conditions = ps.trim.analyze_steady_trim(
    problem=default_problem,
    velocity_bounds=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-1, 1),
    external_thrust_bounds=(0, 10),
)

# Log the trim conditions. If these display "nan", then the trim function couldn't
# find a trimmed state.
example_logger.info("Trim Velocity:\t\t\t%.2f m/s" % trim_conditions[0])
example_logger.info("Trim Alpha:\t\t\t%.2f deg" % trim_conditions[1])
example_logger.info("Trim Beta:\t\t\t\t%.2f deg" % trim_conditions[2])
example_logger.info("Trim External Thrust:\t%.2f N" % trim_conditions[3])

# The expected results are:
# Trim Velocity: 10.02 m/s
# Trim Alpha: 3.00 deg
# Trim Beta: -0.00 deg
# Trim External Thrust: 4.98 N
