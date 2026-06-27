"""Demonstrates analyzing the convergence of an unsteady simulation with non-static
geometry.

The script will likely take several minutes to run, and will log convergence progress
and results in a log file.
"""

import logging

import pterasoftware as ps

# Configure logging to write info level messages to a file. To display log messages on
# the console instead, omit the handler argument.
ps.set_up_logging(level="Info", handler=logging.FileHandler("example_convergence.log"))

# Create an Airplane and AirplaneMovement
example_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=8,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="asymmetric",
                    spanwise_spacing="uniform",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 3.0, 0.0),
                    control_surface_symmetry_type="asymmetric",
                    spanwise_spacing=None,
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.5, 0.0),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            chordwise_spacing="uniform",
        ),
    ],
    name="Example Airplane",
)
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
                    ].wing_cross_sections[1]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
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
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
    ],
)

# Create an OperatingPoint and an OperatingPointMovement.
example_operating_point = ps.operating_point.OperatingPoint()
example_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=example_operating_point
    )
)

# Create a Movement using the AirplaneMovement and OperatingPointMovement
# objects.
example_movement = ps.movements.movement.Movement(
    airplane_movements=[example_airplane_movement],
    operating_point_movement=example_operating_point_movement,
    num_cycles=1,
)

del example_airplane_movement
del example_operating_point_movement

# Create an UnsteadyProblem. We will pass this into the convergence function.
example_problem = ps.problems.UnsteadyProblem(
    movement=example_movement, only_final_results=True
)

del example_movement

# Run the unsteady convergence analysis. This will run several simulations, modifying
# the wake state, wake length, average Panel aspect ratio, and number of chordwise
# Panels with each iteration. Once it detects that the net load coefficients haven't
# change by more than the convergence criteria (measured as an absolute percent
# error), it will return the parameters it found to result in a converged solution.
# See the analyze_unsteady_convergence function docstring for more details.
ps.convergence.analyze_unsteady_convergence(
    ref_problem=example_problem,
    prescribed_wake=True,
    free_wake=True,
    num_cycles_bounds=(1, 4),
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 5),
    convergence_criteria=2.5,
    show_solver_progress=False,
)
