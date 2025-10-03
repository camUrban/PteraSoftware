# NOTE: I haven't yet started refactoring this module.
"""This script is an example of analyzing the unsteady convergence of a problem with
variable geometry. It should take about 30 minutes to run. It will display the
convergence progress and results in the console."""

import pterasoftware as ps

# Create an airplane and airplane movement object. Read through the unsteady solver
# examples for more details on creating this object.
airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    spanwise_spacing="uniform",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    x_le=0.0,
                    y_le=3.0,
                    z_le=0.0,
                    chord=1.0,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                ),
            ],
            symmetric=True,
            chordwise_spacing="uniform",
        ),
    ],
)
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane.wings[0].wing_cross_sections[0],
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane.wings[0].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=1.0,
                ),
            ],
        )
    ],
)

# Create an operating point and an operating point movement object.
operating_point = ps.operating_point.OperatingPoint()
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=operating_point
)

# Create a movement object from the airplane movement and operating point movement
# objects.
movement = ps.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
)

del airplane_movement
del operating_point_movement

# Create an unsteady problem. We will pass this into the convergence function.
problem = ps.problems.UnsteadyProblem(movement=movement, only_final_results=True)

del movement

# Run the unsteady convergence analysis. This will run the problem several times,
# modifying the wake state, wake length, average panel aspect ratio, and number of
# chordwise panels with each iteration. Once it detects that the net load
# coefficients haven't change by more than the convergence criteria (measured as an
# absolute percent error), it will return the parameters it found to result in a
# converged solution. See the analyze_unsteady_convergence function docstring for
# more details. The progress and results are displayed to the console.
ps.convergence.analyze_unsteady_convergence(
    ref_problem=problem,
    prescribed_wake=True,
    free_wake=True,
    num_cycles_bounds=(1, 4),
    panel_aspect_ratio_bounds=(2, 1),
    num_chordwise_panels_bounds=(4, 7),
    convergence_criteria=1.0,
)

# Check the console that the convergence analysis found that the solution converged
# with the following parameters:
# Wake state: free
# Wake length: 2 cycles
# Panel aspect ratio: 1
# Chordwise panels: 5
