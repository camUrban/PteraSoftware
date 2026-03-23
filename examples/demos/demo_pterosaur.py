from email.mime import base

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


def interpolate_between_wing_cross_sections(wcs1, wcs2, first_wcs):
    """
    Wing cross section panels are between the line of wcs1 and wcs2.
    When exploding a wing to 1 spanwise panel per cross section,
    we need to interpolate the intermediate cross sections.
    """

    interpolated = []

    if first_wcs:
        interpolated.append(
            ps.geometry.wing_cross_section.WingCrossSection(
                num_spanwise_panels=1,
                chord=wcs1.chord,
                Lp_Wcsp_Lpp=wcs1.Lp_Wcsp_Lpp,
                angles_Wcsp_to_Wcs_ixyz=wcs1.angles_Wcsp_to_Wcs_ixyz,
                control_surface_symmetry_type=wcs1.control_surface_symmetry_type,
                control_surface_hinge_point=wcs1.control_surface_hinge_point,
                control_surface_deflection=wcs1.control_surface_deflection,
                spanwise_spacing="uniform",
                airfoil=wcs1.airfoil,
            )
        )

    N = wcs1.num_spanwise_panels

    for i in range(N):
        t = (i + 1) / N  # interpolation parameter between 0 and 1

        chord = (1 - t) * wcs1.chord + t * wcs2.chord
        Lp_Wcsp_Lpp = tuple(np.array(wcs2.Lp_Wcsp_Lpp) / N)
        # angles_Wcsp_to_Wcs_ixyz = tuple((1 - t) * np.array(wcs1.angles_Wcsp_to_Wcs_ixyz) + t * np.array(wcs2.angles_Wcsp_to_Wcs_ixyz))
        angles_Wcsp_to_Wcs_ixyz = wcs2.angles_Wcsp_to_Wcs_ixyz / N
        is_final_section = wcs2.num_spanwise_panels is None and i == N - 1 

        interpolated.append(
            ps.geometry.wing_cross_section.WingCrossSection(
                num_spanwise_panels=None if is_final_section else 1,
                chord=chord,
                Lp_Wcsp_Lpp=Lp_Wcsp_Lpp,
                angles_Wcsp_to_Wcs_ixyz=angles_Wcsp_to_Wcs_ixyz,
                control_surface_symmetry_type=wcs1.control_surface_symmetry_type,
                control_surface_hinge_point=wcs1.control_surface_hinge_point,
                control_surface_deflection=wcs1.control_surface_deflection,
                spanwise_spacing=None if is_final_section else "uniform",
                airfoil=wcs1.airfoil,
            )
        )
    return interpolated

def explode_wing(wing):
    """
    Takes a ps.geometry.wing.Wing and returns a NEW Wing
    where all cross sections have num_spanwise_panels = 1.
    """

    new_cross_sections = []

    for i in range(len(wing.wing_cross_sections) - 1):
        new_cross_sections.extend(
            interpolate_between_wing_cross_sections(
                wing.wing_cross_sections[i], wing.wing_cross_sections[i + 1], i == 0
            )
        )

    # Rebuild the wing (copying everything else verbatim)
    return ps.geometry.wing.Wing(
        wing_cross_sections=new_cross_sections,
        name=wing.name,
        Ler_Gs_Cgs=wing.Ler_Gs_Cgs,
        angles_Gs_to_Wn_ixyz=wing.angles_Gs_to_Wn_ixyz,
        symmetric=wing.symmetric,
        mirror_only=wing.mirror_only,
        symmetryNormal_G=wing.symmetryNormal_G,
        symmetryPoint_G_Cg=wing.symmetryPoint_G_Cg,
        num_chordwise_panels=wing.num_chordwise_panels,
        chordwise_spacing=wing.chordwise_spacing,
    )


wing_cross_section_1 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=4,
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

# points = [(0.0, 0.0, 0.0),
#           (0.0889,0.2249, 0.0955),
#           (-0.0521, 0.5749, 0.1940),
#           (-0.0946, 0.8282, 0.2345),
#           (0.1829, 2.4373, 0.0266)]

# vectors = [(1.0, 0.0, 0.0),
#            (0.1559,0,-0.0931),
#            (0.2864,0,-0.1878),
#            (0.3291,0,-0.2154),
#            (0.1829, 2.4373, 0.0266)]

wing_cross_section_2 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=1,
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
                    num_spanwise_panels=1,
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
                    num_spanwise_panels=1,
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
wing_cross_section_4_5 = ps.geometry.wing_cross_section.WingCrossSection(
    num_spanwise_panels=1,
    chord=(0.39316581743584983 + 0.005008991914547277) / 2,
    Lp_Wcsp_Lpp= 0.0 * np.array((-0.05774876, 0.2533, 0.01056319)) + 0.5 * np.array((0.34656431, 1.6091, -0.01103809)),
    angles_Wcsp_to_Wcs_ixyz=get_relative_transform(
        (-0.0946, 0.8282, 0.2345),
        (0.322, 0, -0.2256),
        (0.1829, 2.4373, 0.0266),
        (0.005, 0, -0.0003),
    )[1]  * 0,
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

wing_cross_section_new_5 = ps.geometry.wing_cross_section.WingCrossSection(
    num_spanwise_panels=None,
    chord=np.linalg.norm((0.005, 0, -0.0003)),
    Lp_Wcsp_Lpp=np.array([0.34656431, 1.6091, -0.01103809]) / 2,
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

wing_cross_section_5 = ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
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

original_wing = ps.geometry.wing.Wing(
        wing_cross_sections=[
            wing_cross_section_1,
            wing_cross_section_2,
            wing_cross_section_3,
            wing_cross_section_4,
            wing_cross_section_5,
        ],
        name="Main Wing",
        Ler_Gs_Cgs=[0.0, 0.025, 0.0],
        angles_Gs_to_Wn_ixyz=[4, 0.0, 0.0],
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        single_step_wing = True,
        num_chordwise_panels=5,
        chordwise_spacing="uniform",
    )

pterasaure = ps.geometry.airplane.Airplane(
    wings=[original_wing],
    name="Pterosaur",
    Cg_GP1_CgP1=(0.0, 0.0, 0.0),
    weight=0,
    s_ref=None,
    c_ref=None,
    b_ref=None,
)

dephase_x = 0.0
period_x = 0.0
amplitude_x = 0.0

dephase_y = 0.0
period_y = 0.0
amplitude_y = 0.0

dephase_z = 0.0
period_z = 0.0
amplitude_z = 0.0

# A list of movements for the main wing
main_single_step_movements_list = []

# A list of movements for the reflected wing
reflected_single_step_movements_list = []

for i in range(len(pterasaure.wings[0].wing_cross_sections)):
    if i == 0:
        single_step_movement = ps.movements.single_step.single_step_wing_cross_section_movement.SingleStepWingCrossSectionMovement(
            base_wing_cross_section=pterasaure.wings[0].wing_cross_sections[i],
        )

        main_single_step_movements_list.append(single_step_movement)
        reflected_single_step_movements_list.append(single_step_movement)

    else:
        single_step_movement = ps.movements.single_step.single_step_wing_cross_section_movement.SingleStepWingCrossSectionMovement(
            base_wing_cross_section=pterasaure.wings[0].wing_cross_sections[i],
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(amplitude_x, amplitude_y, amplitude_z),
            periodAngles_Wcsp_to_Wcs_ixyz=(period_x, period_y, period_z),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(dephase_x, dephase_y, dephase_z),
        )

        main_single_step_movements_list.append(single_step_movement)
        reflected_single_step_movements_list.append(single_step_movement)

single_step_main_wing_movement = (
    ps.movements.single_step.single_step_wing_movement.SingleStepWingMovement(
        base_wing=pterasaure.wings[0],
        single_step_wing_cross_section_movements=main_single_step_movements_list,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(30.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1 / 3, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )
)

single_step_reflected_main_wing_movement = (
    ps.movements.single_step.single_step_wing_movement.SingleStepWingMovement(
        base_wing=pterasaure.wings[1],
        single_step_wing_cross_section_movements=reflected_single_step_movements_list,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(30.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1 / 3, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )
)

# Now define the example airplane's AirplaneMovement. For now, no movement of the airplane is possible.
pterasaure_single_step_movement = (
    ps.movements.single_step.single_step_airplane_movement.SingleStepAirplaneMovement(
        base_airplane=pterasaure,
        single_step_wing_movements=[
            single_step_main_wing_movement,
            single_step_reflected_main_wing_movement,
        ],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )
)

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=30.0, alpha=5.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define the operating point's OperatingPointMovement.
single_step_operating_point_movement = ps.movements.single_step.single_step_operating_point_movement.SingleStepOperatingPointMovement(
    base_operating_point=example_operating_point, ampVCg__E=0.0, periodVCg__E=0.0, spacingVCg__E="sine"
)

# Define the Movement. This contains the AirplaneMovement and the
# OperatingPointMovement.
single_step_movement = ps.movements.single_step.single_step_movement.SingleStepMovement(
    single_step_airplane_movements=[pterasaure_single_step_movement],
    single_step_operating_point_movement=single_step_operating_point_movement,
    delta_time=None,
    num_cycles=2,
)

# Define the UnsteadyProblem.
example_problem = ps.problems.AeroelasticUnsteadyProblem(
    single_step_movement=single_step_movement,
    wing_density=1,
    spring_constant=0.1,
    plot_flap_cycle=False,
    damping_constant=1.0,
)

# Define a new solver. The available solver classes are
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver,
# and UnsteadyRingVortexLatticeMethodSolver. We'll create an
# UnsteadyRingVortexLatticeMethodSolver, which requires a UnsteadyProblem.
example_solver = (
    ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
        coupled_unsteady_problem=example_problem,
    )
)

# Run the solver.
example_solver.run(
    prescribed_wake=True,
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
