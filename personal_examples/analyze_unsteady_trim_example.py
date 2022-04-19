# ToDo: Document this script.
import logging

import pterasoftware as ps

example_logger = logging.getLogger("example")
example_logger.setLevel(logging.DEBUG)

default_airplane = ps.geometry.Airplane(
    x_ref=0.14,
    weight=420,
    wings=[
        ps.geometry.Wing(
            symmetric=True,
            num_chordwise_panels=5,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    num_spanwise_panels=5,
                    airfoil=ps.geometry.Airfoil(
                        name="naca2412",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
        ),
        ps.geometry.Wing(
            x_le=5,
            symmetric=True,
            num_chordwise_panels=5,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    num_spanwise_panels=5,
                    twist=-5.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=1.0,
                    chord=1.0,
                    twist=-5.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

default_airplane_movement = ps.movement.AirplaneMovement(
    base_airplane=default_airplane,
    wing_movements=[
        ps.movement.WingMovement(
            base_wing=default_airplane.wings[0],
            wing_cross_sections_movements=[
                ps.movement.WingCrossSectionMovement(
                    base_wing_cross_section=default_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movement.WingCrossSectionMovement(
                    base_wing_cross_section=default_airplane.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_period=1.0,
                    sweeping_amplitude=5.0,
                    heaving_period=1.0,
                    heaving_amplitude=10.0,
                ),
            ],
        ),
        ps.movement.WingMovement(
            base_wing=default_airplane.wings[1],
            wing_cross_sections_movements=[
                ps.movement.WingCrossSectionMovement(
                    base_wing_cross_section=default_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movement.WingCrossSectionMovement(
                    base_wing_cross_section=default_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

default_operating_point = ps.operating_point.OperatingPoint()

default_operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=default_operating_point
)

default_movement = ps.movement.Movement(
    airplane_movements=[default_airplane_movement],
    operating_point_movement=default_operating_point_movement,
)

default_problem = ps.problems.UnsteadyProblem(
    movement=default_movement,
    only_final_results=False,
)

default_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=default_problem
    )
)
default_solver.run(logging_level="Critical")

trim_conditions = ps.trim.analyze_unsteady_trim(
    airplane_movement=default_airplane_movement,
    operating_point=default_operating_point,
    velocity_bounds=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-0.1, 0.1),
)

example_logger.info("Trim Velocity:\t%.2f m/s" % trim_conditions[0])
example_logger.info("Trim Alpha:\t%.2f deg" % trim_conditions[1])
example_logger.info("Trim Beta:\t\t%.2f deg" % trim_conditions[2])
