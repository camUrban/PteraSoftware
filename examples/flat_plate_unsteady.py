import pterasoftware as ps

chord = 1.0
aspect_ratio = 4.0
num_chordwise_panels = 5

num_spanwise_panels = round(num_chordwise_panels * aspect_ratio)
tip_y_le = chord * aspect_ratio

flat_plate_airplane = ps.geometry.Airplane(
    name="Flat Plate Airplane",
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
    wings=[
        ps.geometry.Wing(
            name="Flat Plate Wing",
            x_le=0.0,
            y_le=0.0,
            z_le=0.0,
            symmetric=False,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    chord=chord,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0000",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=tip_y_le,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0000",
                    ),
                ),
            ],
        ),
    ],
)

del chord
del aspect_ratio
del num_chordwise_panels
del num_spanwise_panels
del tip_y_le

flat_plate_wing_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=flat_plate_airplane.wings[0].wing_cross_sections[0],
    sweeping_amplitude=0.0,
    sweeping_period=0.0,
    sweeping_spacing="sine",
    pitching_amplitude=0.0,
    pitching_period=0.0,
    pitching_spacing="sine",
    heaving_amplitude=0.0,
    heaving_period=0.0,
    heaving_spacing="sine",
)

flat_plate_wing_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=flat_plate_airplane.wings[0].wing_cross_sections[1],
)

flat_plate_wing_movement = ps.movement.WingMovement(
    base_wing=flat_plate_airplane.wings[0],
    wing_cross_sections_movements=[
        flat_plate_wing_root_wing_cross_section_movement,
        flat_plate_wing_tip_wing_cross_section_movement,
    ],
    x_le_amplitude=0.0,
    x_le_period=0.0,
    x_le_spacing="sine",
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing="sine",
    z_le_amplitude=0.0,
    z_le_period=0.0,
    z_le_spacing="sine",
)

del flat_plate_wing_root_wing_cross_section_movement
del flat_plate_wing_tip_wing_cross_section_movement

flat_plate_airplane_movement = ps.movement.AirplaneMovement(
    base_airplane=flat_plate_airplane,
    wing_movements=[flat_plate_wing_movement],
    x_ref_amplitude=0.0,
    x_ref_period=0.0,
    x_ref_spacing="sine",
    y_ref_amplitude=0.0,
    y_ref_period=0.0,
    y_ref_spacing="sine",
    z_ref_amplitude=0.0,
    z_ref_period=0.0,
    z_ref_spacing="sine",
)

del flat_plate_airplane
del flat_plate_wing_movement

flat_plate_operating_point = ps.operating_point.OperatingPoint(
    density=1.225,
    beta=0.0,
    velocity=1.0,
    alpha=15.0,
    nu=15.06e-6,
)

flat_plate_operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=flat_plate_operating_point,
    velocity_amplitude=0.0,
    velocity_period=0.0,
    velocity_spacing="sine",
)

del flat_plate_operating_point

flat_plate_movement = ps.movement.Movement(
    airplane_movements=[flat_plate_airplane_movement],
    operating_point_movement=flat_plate_operating_point_movement,
    num_steps=None,
    delta_time=None,
)

del flat_plate_airplane_movement
del flat_plate_operating_point_movement

flat_plate_problem = ps.problems.UnsteadyProblem(
    movement=flat_plate_movement,
    only_final_results=False,
)

del flat_plate_movement

flat_plate_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=flat_plate_problem,
    )
)

del flat_plate_problem

flat_plate_solver.run(
    logging_level="Warning",
    prescribed_wake=False,
    calculate_streamlines=False,
)

ps.output.draw(
    solver=flat_plate_solver,
    scalar_type="lift",
    show_streamlines=False,
    show_wake_vortices=True,
    save=False,
)

# ps.output.animate(
#     unsteady_solver=flat_plate_solver,
#     scalar_type="lift",
#     show_wake_vortices=True,
#     save=False,
# )

ps.output.plot_results_versus_time(
    unsteady_solver=flat_plate_solver,
    show=False,
    save=True,
)
