# REFACTOR: I haven't yet started refactoring this module.
# DOCUMENT: Document this script.
import pterasoftware as ps

# Known Converged Values (0.5% convergence, 0 degrees angle of attack):
#   1 Airplane (flaps, chordwise panels, prescribed wake):
#       2, 7, False
#   3 Airplanes (flaps, chordwise panels, prescribed wake):
#       3, 7, False
#   5 Airplanes (flaps, chordwise panels, prescribed wake):
#       3, 7, False
num_airplanes = 1
num_flaps = 2
num_chord = 7
prescribed_wake = False

aspect_ratio = 5.0
speed = 1.0
alpha = 0.0
x_spacing = 0.5
y_spacing = 0.5
root_to_mid_span = 0.2275
root_chord = 0.1094
mid_to_tip_span = 0.350 - 0.2275
tip_chord = 0.0219
flapping_amplitude = 15.0

period = x_spacing / speed
root_to_mid_chord = root_chord
mid_to_tip_chord = (root_chord + tip_chord) / 2

this_operating_point = ps.operating_point.OperatingPoint(vCg__E=speed, alpha=0.0)
this_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=this_operating_point
    )
)
del this_operating_point

print("Prescribed Wake:", prescribed_wake)
print("Number of flaps:", num_flaps)
print("Number of chordwise panels:", num_chord)

root_to_mid_panel_chord = root_to_mid_chord / num_chord
mid_to_tip_panel_chord = mid_to_tip_chord / num_chord

root_to_mid_num_span = round(
    root_to_mid_span / (aspect_ratio * root_to_mid_panel_chord)
)
mid_to_tip_num_span = round(mid_to_tip_span / (aspect_ratio * mid_to_tip_panel_chord))

these_airplane_movements = []
row = None
position = None
offset_sign = None
for airplane_id in range(num_airplanes):
    if airplane_id == 0:
        row = 1
        position = ""
        offset_sign = 1
    elif airplane_id % 2 != 0:
        row += 1
        position = "Right "
        offset_sign = 1
    else:
        position = "Left "
        offset_sign = -1

    this_name = "Airplane " + position + str(row)

    offset = row - 1

    this_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        twist=alpha,
                        chord=root_chord,
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=root_to_mid_num_span,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        twist=alpha,
                        y_le=root_to_mid_span,
                        chord=root_chord,
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=mid_to_tip_num_span,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        twist=alpha,
                        y_le=root_to_mid_span + mid_to_tip_span,
                        chord=tip_chord,
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                    ),
                ],
                name="Main Wing",
                symmetric=True,
                num_chordwise_panels=num_chord,
                chordwise_spacing="uniform",
            ),
        ],
        name=this_name,
    )

    this_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=this_airplane,
        wing_movements=[
            ps.movements.wing_movement.WingMovement(
                base_wing=this_airplane.wings[0],
                wing_cross_section_movements=[
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=this_airplane.wings[
                            0
                        ].wing_cross_sections[0]
                    ),
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=this_airplane.wings[
                            0
                        ].wing_cross_sections[1]
                    ),
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=this_airplane.wings[
                            0
                        ].wing_cross_sections[2]
                    ),
                ],
            )
        ],
    )

    these_airplane_movements.append(this_airplane_movement)

    del this_airplane
    del this_airplane_movement

this_movement = ps.movement.Movement(
    airplane_movements=these_airplane_movements,
    operating_point_movement=this_operating_point_movement,
    num_steps=None,
    num_cycles=num_flaps,
    delta_time=None,
)

del these_airplane_movements

this_problem = ps.problems.UnsteadyProblem(
    movement=this_movement,
    only_final_results=False,
)

del this_movement

this_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=this_problem,
    )
)

del this_problem

this_solver.run(
    prescribed_wake=prescribed_wake,
    calculate_streamlines=False,
)

ps.output.print_unsteady_results(unsteady_solver=this_solver)

ps.output.plot_results_versus_time(unsteady_solver=this_solver, save=True)

ps.output.draw(
    solver=this_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
)

ps.output.animate(
    unsteady_solver=this_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
)
