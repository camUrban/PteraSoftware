"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import src

this_airplane = src.geometry.Airplane(
    name="Experimental Airplane",
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0875,
                    y_le=0.350,
                    chord=0.0219,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
    ],
)

this_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
    velocity=1.0,
    alpha=5.0,
    nu=15.06e-6,
)

this_airplane_movement = src.movement.AirplaneMovement(
    base_airplane=this_airplane,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=this_airplane.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=this_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=this_airplane.wings[0].wing_cross_sections[
                        1
                    ],
                    sweeping_amplitude=15.0,
                    sweeping_period=1,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

this_operating_point_movement = src.movement.OperatingPointMovement(
    base_operating_point=this_operating_point,
)

this_movement = src.movement.Movement(
    airplane_movements=[this_airplane_movement],
    operating_point_movement=this_operating_point_movement,
)

del this_airplane_movement
del this_operating_point_movement

this_problem = src.problems.UnsteadyProblem(
    movement=this_movement, only_final_results=True
)

this_solver = (
    src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=this_problem,
    )
)

del this_problem

this_solver.run(
    prescribed_wake=True,
    calculate_streamlines=False,
)

src.output.draw(
    solver=this_solver,
    show_delta_pressures=True,
    show_wake_vortices=True,
)

# src.output.animate(
#     unsteady_solver=this_solver,
#     show_delta_pressures=True,
#     show_wake_vortices=True,
#     keep_file=False,
# )
