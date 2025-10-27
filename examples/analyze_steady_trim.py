"""This script is an example of how to automatically find the trim condition for a
steady simulation."""

import pterasoftware as ps

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
                    num_spanwise_panels=8,
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
                ),
            ],
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=8,
                    chord=0.5,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    chord=0.5,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            Ler_Gs_Cgs=(5.0, 0.0, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        ),
    ],
    weight=250,
)

# Create an OperatingPoint. We must specify an external thrust because this Airplane
# is not flapping, so it won't generate thrust via its Wings. Therefore, to balance
# induced drag, we need an external thrust force which could be due to a propeller or
# other type of engine.
trim_operating_point = ps.operating_point.OperatingPoint(externalFX_W=5)

# Construct a SteadyProblem containing the Airplane and OperatingPoint
trim_problem = ps.problems.SteadyProblem(
    airplanes=[trim_airplane], operating_point=trim_operating_point
)

# Call the analyze_steady_trim function to search for a trim condition (thrust
# balances drag, weight balances lift, and all moments are close to zero) within a
# certain set of bounds.
trim_conditions = ps.trim.analyze_steady_trim(
    problem=trim_problem,
    solver_type="steady horseshoe vortex lattice method",
    boundsVCg__E=(5, 15),
    alpha_bounds=(-10, 10),
    beta_bounds=(-1, 1),
    boundsExternalFX_W=(0, 10),
    objective_cut_off=0.01,
    num_calls=1000,
)
