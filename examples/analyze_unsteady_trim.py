"""This script is an example of how to automatically find the trim condition for an
unsteady simulation."""

import pterasoftware as ps

# Configure logging to display info level messages. This is important for seeing the
# output from the trim function.
ps.set_up_logging(level="Info")

# Create an Airplane. We must specify a weight (in Newtons) for the Airplane. We will
# later find a trim condition where the weight is exactly balanced by lift.
trim_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=5,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=5,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=5,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            Ler_Gs_Cgs=(10, 0, 0),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=5,
        ),
    ],
    weight=420,
)

# Create an AirplaneMovement for this example's Airplane.
trim_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=trim_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=trim_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trim_airplane.wings[0].wing_cross_sections[
                        0
                    ]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trim_airplane.wings[0].wing_cross_sections[
                        1
                    ]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=trim_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trim_airplane.wings[1].wing_cross_sections[
                        0
                    ]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=trim_airplane.wings[1].wing_cross_sections[
                        1
                    ]
                ),
            ],
        ),
    ],
)

# Create an OperatingPoint using default values.
trim_operating_point = ps.operating_point.OperatingPoint(
    externalFX_W=7.5,
)

# Create an OperatingPointMovement using default values.
trim_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=trim_operating_point
    )
)

trim_movement = ps.movements.movement.Movement(
    airplane_movements=[trim_airplane_movement],
    operating_point_movement=trim_operating_point_movement,
    num_chords=5,
)

trim_problem = ps.problems.UnsteadyProblem(
    movement=trim_movement,
    only_final_results=True,
)

# Call the analyze_unsteady_trim function to search for a trim condition (thrust
# balances drag, weight balances lift, and all moments are close to zero) within a
# certain set of bounds.
trim_conditions = ps.trim.analyze_unsteady_trim(
    problem=trim_problem,
    boundsVCg__E=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-0.1, 0.1),
    boundsExternalFX_W=(5, 15),
    objective_cut_off=0.01,
    num_calls=100,
    show_solver_progress=True,
)
