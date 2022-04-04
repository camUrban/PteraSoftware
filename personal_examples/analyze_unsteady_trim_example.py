# ToDo: Document this script.
import pterasoftware as ps

default_airplane = ps.geometry.Airplane(
    x_ref=0.25,
    weight=500,
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
                    sweeping_amplitude=10.0,
                    heaving_period=1.0,
                    heaving_amplitude=10.0,
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
default_solver.run()

ps.output.animate(
    unsteady_solver=default_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=False,
)
ps.output.print_unsteady_results(unsteady_solver=default_solver)

print("Untrimmed Results:")

trim_conditions = ps.trim.analyze_unsteady_trim(
    airplane_movement=default_airplane_movement,
    operating_point=default_operating_point,
    velocity_bounds=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-1, 1),
    objective_cut_off=1,
    num_calls=100,
)
#
# print("\nTrim Values:")
# print("Velocity:\t%.2f m/s" % trim_conditions[0])
# print("Alpha:\t\t%.2f deg" % trim_conditions[1])
# print("Beta:\t\t%.2f deg" % trim_conditions[2])
# print("External Thrust:\t%.2f N" % trim_conditions[3])
#
# trim_operating_point = ps.operating_point.OperatingPoint(
#     velocity=trim_conditions[0],
#     alpha=trim_conditions[1],
#     beta=trim_conditions[2],
#     external_thrust=trim_conditions[3],
# )
# trim_problem = ps.problems.SteadyProblem(
#     operating_point=trim_operating_point,
#     airplanes=[default_airplane],
# )
# trim_solver = (
#     ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
#         steady_problem=trim_problem
#     )
# )
# trim_solver.run()
#
# print("\nTrimmed Results:")
# ps.output.print_steady_results(steady_solver=trim_solver)
