"""This is script is an example of how to run Ptera Software's
UnsteadyRingVortexLatticeMethodSolver on three Airplanes, flying in formation,
each with custom geometry and non static motion."""

import pterasoftware as ps

x_spacing = 13
y_spacing = 13

# Create the lead Airplane.
lead_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=10,
                    chord=1.75,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.75, 6.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.25, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
    name="Lead Airplane",
    Cg_E_CgP1=(0.0, 0.0, 0.0),
    angles_E_to_B_izyx=(0.0, 0.0, 0.0),
)

# Now define the lead Airplane's AirplaneMovement.
lead_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=lead_airplane,
    wing_movements=[
        # Define the main Wing's WingMovement.
        ps.movements.wing_movement.WingMovement(
            base_wing=lead_airplane.wings[0],
            wing_cross_section_movements=[
                # Define the root WingCrossSection's WingCrossSectionMovement.
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[
                        0
                    ]
                ),
                # Define the tip WingCrossSection's WingCrossSectionMovement.
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[
                        1
                    ]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(25.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
        # Define the reflected main Wing's WingMovement.
        ps.movements.wing_movement.WingMovement(
            base_wing=lead_airplane.wings[1],
            wing_cross_section_movements=[
                # Define the root WingCrossSection's WingCrossSectionMovement.
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[1].wing_cross_sections[
                        0
                    ]
                ),
                # Define the tip WingCrossSection's WingCrossSectionMovement.
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[1].wing_cross_sections[
                        1
                    ]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(25.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
    ],
)

# Create the trailing right Airplane.
trailing_right_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=10,
                    chord=1.75,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.75, 6.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.25, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
    name="Trailing Right Airplane",
    Cg_E_CgP1=(x_spacing, y_spacing, 0.0),
    angles_E_to_B_izyx=(0.0, 0.0, 0.0),
)

# Create the trailing right Airplane's AirplaneMovement.
trailing_right_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=trailing_right_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=trailing_right_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_right_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_right_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(25.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=trailing_right_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_right_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_right_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(25.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
    ],
)

# Create the trailing left Airplane.
trailing_left_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=10,
                    chord=1.75,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.75, 6.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.25, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
    name="Trailing Left Airplane",
    Cg_E_CgP1=(x_spacing, -y_spacing, 0.0),
    angles_E_to_B_izyx=(0.0, 0.0, 0.0),
)

# Create the trailing left Airplane's AirplaneMovement.
left_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=trailing_left_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=trailing_left_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_left_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_left_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(25.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=trailing_left_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_left_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trailing_left_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
            ampAngles_Gs_to_Wn_ixyz=(25.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        ),
    ],
)

# Define an OperatingPoint. This defines the state at which all the Airplanes are
# operating.
operating_point = ps.operating_point.OperatingPoint(vCg__E=15.0, alpha=0.0)

# Define the OperatingPoint's OperatingPointMovement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=operating_point
)

# Delete the extraneous pointers to the Airplanes and the OperatingPoint, as these
# are now accessible within their respective movement objects.
del lead_airplane
del trailing_right_airplane
del trailing_left_airplane
del operating_point

# Define the Movement. This contains each AirplaneMovement and the
# OperatingPointMovement.
movement = ps.movements.movement.Movement(
    airplane_movements=[
        lead_airplane_movement,
        trailing_right_airplane_movement,
        left_airplane_movement,
    ],
    operating_point_movement=operating_point_movement,
    num_cycles=3,
)

# Delete the extraneous pointers to the AirplaneMovements and the
# OperatingPointMovement, as these are now accessible within the Movement.
del lead_airplane_movement
del trailing_right_airplane_movement
del left_airplane_movement
del operating_point_movement

# Using the Movement, create an UnsteadyProblem.
unsteady_problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Define a new UnsteadyRingVortexLatticeMethodSolver.
solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    unsteady_problem=unsteady_problem,
)

# Delete the extraneous pointer to the UnsteadyProblem.
del unsteady_problem

# Run the UnsteadyRingVortexLatticeMethodSolver.
solver.run(
    prescribed_wake=True,
)

# Now that we have run the solver, we can create an animation of the results.
ps.output.animate(
    unsteady_solver=solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=False,
)
