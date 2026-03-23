import pterasoftware as ps
import numpy as np
from scipy.spatial.transform import Rotation as R

def get_relative_transform(Point1, Unit1, Point2, Unit2):

    Point1 = np.array(Point1, dtype=float)
    Unit1  = np.array(Unit1, dtype=float)
    Point2 = np.array(Point2, dtype=float)
    Unit2  = np.array(Unit2, dtype=float)

    Xp = Unit1 / np.linalg.norm(Unit1)
    Zp = np.array([0, 0, 1])
    Yp = np.cross(Zp, Xp)
    Yp /= np.linalg.norm(Yp)
    Zp = np.cross(Xp, Yp)
    Zp /= np.linalg.norm(Zp)
    Rp = np.column_stack((Xp, Yp, Zp))

    Xe = Unit2 / np.linalg.norm(Unit2)
    Ze = np.array([0, 0, 1])
    Ye = np.cross(Ze, Xe)
    Ye /= np.linalg.norm(Ye)
    Ze = np.cross(Xe, Ye)
    Ze /= np.linalg.norm(Ze)
    Re = np.column_stack((Xe, Ye, Ze))

    Rrel = Rp.T @ Re

    rot = R.from_matrix(Rrel)
    angles_Wcsp_to_Wcs_ixyz = rot.as_euler('xyz', degrees=True)

    Lp_Wcsp_Lpp = Rp.T @ (Point2 - Point1)

    return Lp_Wcsp_Lpp, angles_Wcsp_to_Wcs_ixyz


wing_cross_section_1 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=6,
                    chord=0.25,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="uniform",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                )

wing_cross_section_2 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=6,
                    chord=np.linalg.norm((0.1559,0,-0.0931)),
                    Lp_Wcsp_Lpp=get_relative_transform((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0889,0.2249, 0.0955), (0.1559,0,-0.0931))[0],
                    angles_Wcsp_to_Wcs_ixyz=get_relative_transform((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0889,0.2249, 0.0955), (0.1559,0,-0.0931))[1],
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="uniform",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                )

wing_cross_section_3 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=6,
                    chord=np.linalg.norm((0.2864,0,-0.1878)),
                    Lp_Wcsp_Lpp=get_relative_transform((0.0889,0.2249, 0.0955), (0.1559,0,-0.0931), (-0.0521, 0.5749, 0.1940), (0.2864,0,-0.1878))[0],
                    angles_Wcsp_to_Wcs_ixyz=get_relative_transform((0.0889,0.2249, 0.0955), (0.1559,0,-0.0931), (-0.0521, 0.5749, 0.1940), (0.2864,0,-0.1878))[1],
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="uniform",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                )

wing_cross_section_4 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=6,
                    chord=np.linalg.norm((0.322,0,-0.2256)),
                    Lp_Wcsp_Lpp=get_relative_transform((-0.0521, 0.5749, 0.1940), (0.2864,0,-0.1878), (-0.0946, 0.8282, 0.2345), (0.322,0,-0.2256))[0],
                    angles_Wcsp_to_Wcs_ixyz=get_relative_transform((-0.0521, 0.5749, 0.1940), (0.2864,0,-0.1878), (-0.0946, 0.8282, 0.2345), (0.322,0,-0.2256))[1],
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="uniform",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                )

wing_cross_section_5 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=1,
                    chord=np.linalg.norm((0.005,0,-0.0003)),
                    Lp_Wcsp_Lpp=get_relative_transform((-0.0946, 0.8282, 0.2345), (0.322,0,-0.2256), (0.1829, 2.4373, 0.0266), (0.005,0,-0.0003))[0],
                    angles_Wcsp_to_Wcs_ixyz=get_relative_transform((-0.0946, 0.8282, 0.2345), (0.322,0,-0.2256), (0.1829, 2.4373, 0.0266), (0.005,0,-0.0003))[1],
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                )
wing_cross_section_6 = ps.geometry.wing_cross_section.WingCrossSection(
    num_spanwise_panels=None,
    chord=np.linalg.norm((0.005, 0, -0.0003)),
    Lp_Wcsp_Lpp=get_relative_transform(
        (0.322, 0, -0.2256),
        (0.323, 0, -0.2251),
        (0.005, 0, -0.0003),
        (0.006, 0, -0.0002),
    )[0],
    angles_Wcsp_to_Wcs_ixyz=get_relative_transform(
        (-0.0946, 0.8282, 0.2345),
        (0.322, 0, -0.2256),
        (0.1829, 2.4373, 0.0266),
        (0.005, 0, -0.0003),
    )[1],
    control_surface_symmetry_type="symmetric",
    control_surface_hinge_point=0.75,
    control_surface_deflection=0.0,
    spanwise_spacing=None,
    airfoil=ps.geometry.airfoil.Airfoil(
        name="naca0012",
        outline_A_lp=None,
        resample=True,
        n_points_per_side=400,
    ),
)


pterasaure = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[wing_cross_section_1, wing_cross_section_2, wing_cross_section_3, wing_cross_section_4, wing_cross_section_5, wing_cross_section_6],
            name="Main Wing",
            Ler_Gs_Cgs= [0.0, 0.025, 0.0],
            angles_Gs_to_Wn_ixyz= [4, 0.0, 0.0],
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=5,
            chordwise_spacing="uniform",
        ),
    ],
    name="Pterosaure",
    Cg_GP1_CgP1=(0.0, 0.0, 0.0),
    weight=0,
    s_ref=None,
    c_ref=None,
    b_ref=None,
)


# Define the airplane's AirplaneMovement.
wing_movements = []
for wing in pterasaure.wings:
    wing_cross_section_movements = []
    for wing_cross_section in wing.wing_cross_sections:
        wing_cross_section_movements.append(ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section)
        )
    wing_movements.append(wing_cross_section_movements)

main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=pterasaure.wings[0],
    wing_cross_section_movements= wing_movements[0],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(30.0, 0.0, 0.0), 
    periodAngles_Gs_to_Wn_ixyz=(1/3, 0.0, 0.0), 
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)

reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=pterasaure.wings[1],
    wing_cross_section_movements=wing_movements[1],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(30.0, 0.0, 0.0),  
    periodAngles_Gs_to_Wn_ixyz=(1/3, 0.0, 0.0), 
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)

# Now define the example airplane's AirplaneMovement. For now, no movement of the airplane is possible.
pterasaure_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=pterasaure,
    wing_movements=[main_wing_movement, reflected_main_wing_movement],
    ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
    periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
    spacingCg_GP1_CgP1=("sine", "sine", "sine"),
    phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
)

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=30.0, alpha=5.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define the operating point's OperatingPointMovement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point)

# Define the Movement. This contains the AirplaneMovement and the
# OperatingPointMovement.
movement = ps.movements.movement.Movement(
    airplane_movements=[pterasaure_movement],
    operating_point_movement=operating_point_movement,
    delta_time=None,
    num_cycles=1,
    num_chords=None,
    num_steps=None,
)

# Define the UnsteadyProblem.
example_problem = ps.problems.UnsteadyProblem(
    movement=movement, 
)

# Define a new solver. The available solver classes are
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver,
# and UnsteadyRingVortexLatticeMethodSolver. We'll create an
# UnsteadyRingVortexLatticeMethodSolver, which requires a UnsteadyProblem.
example_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)


# Run the solver.
example_solver.run(
    prescribed_wake=False,
    show_progress=True,
)

# Call the animate function on the solver. This produces a GIF of the wake being
# shed. The GIF is saved in the same directory as this script. Press "q",
# after orienting the view, to begin the animation.

ps.output.animate(
    unsteady_solver=example_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
)

# ps.output.print_results(example_solver)

# ps.output.plot_wing_loads_versus_time(unsteady_solver=example_solver, show=True)
