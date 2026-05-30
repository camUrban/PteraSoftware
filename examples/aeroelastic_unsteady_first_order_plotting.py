"""This script runs the coupled aeroelastic UVLM solver sweeping one parameter
at a time (spring constant k, damping constant b, or wing density) and overlays
the Curve 16 Net Deformation from each run."""

import matplotlib.pyplot as plt
import numpy as np

import pterasoftware as ps

# The curve index to extract from each Net Deformation result. Curve 16 corresponds
# to the wing-tip spanwise station.
CURVE_INDEX = 16

# Default values used when a parameter is not being swept
DEFAULT_K = 1.0
DEFAULT_B = 1.0
DEFAULT_DENSITY = 0.012

# Populate exactly ONE of these lists to sweep that parameter while holding the
# others at their defaults. Leave the other two as empty lists.
K_VALUES: list[float] = [1.0, 4.0, 10.0, 40.0]
B_VALUES: list[float] = []
DENSITY_VALUES: list[float] = []


def run_aeroelastic(
    spring_constant: float = DEFAULT_K,
    damping_constant: float = DEFAULT_B,
    wing_density: float = DEFAULT_DENSITY,
) -> tuple[list, object]:
    """Run the aeroelastic solver and return the net deformation data.

    :param spring_constant: The torsional spring stiffness value.
    :param damping_constant: The damping constant value.
    :param wing_density: The wing density per unit height (kg/m^2).
    :return: A tuple of (net_data list, solved problem object).
    """

    # Wing cross section initialization
    num_spanwise_panels = 2
    Lp_Wcsp_Lpp_Offsets = (0.1, 0.5, 0.0)
    cross_section_chords = [1.75, 1.75, 1.75, 1.75, 1.65, 1.55, 1.4, 1.2, 1.0]
    wing_cross_sections = []

    for i in range(len(cross_section_chords)):
        wing_cross_sections.append(
            ps.geometry.wing_cross_section.WingCrossSection(
                num_spanwise_panels=(
                    num_spanwise_panels if i < len(cross_section_chords) - 1 else None
                ),
                chord=cross_section_chords[i],
                Lp_Wcsp_Lpp=Lp_Wcsp_Lpp_Offsets if i > 0 else (0.0, 0.0, 0.0),
                angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                control_surface_symmetry_type="symmetric",
                control_surface_hinge_point=0.75,
                control_surface_deflection=0.0,
                spanwise_spacing=(
                    "uniform" if i < len(cross_section_chords) - 1 else None
                ),
                airfoil=ps.geometry.airfoil.Airfoil(
                    name="naca2412",
                    outline_A_lp=None,
                    resample=True,
                    n_points_per_side=400,
                ),
            )
        )

    wing_1 = ps.geometry.wing.Wing(
        wing_cross_sections=wing_cross_sections,
        name="Main Wing",
        Ler_Gs_Cgs=(0.0, 0.5, 0.0),
        angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        explode_into_strips=True,
        num_chordwise_panels=6,
        chordwise_spacing="uniform",
    )

    example_airplane = ps.geometry.airplane.Airplane(
        wings=[
            wing_1,
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        num_spanwise_panels=8,
                        chord=1.5,
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
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        num_spanwise_panels=None,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.5, 2.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
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
                    ),
                ],
                name="V-Tail",
                Ler_Gs_Cgs=(5.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                explode_into_strips=False,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
            ),
        ],
        name="Example Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        c_ref=None,
        b_ref=None,
    )

    # Wing cross section movement parameters
    dephase_x = 0.0
    period_x = 0.0
    amplitude_x = 0.0

    dephase_y = 0.0
    period_y = 0.0
    amplitude_y = 0.0

    dephase_z = 0.0
    period_z = 0.0
    amplitude_z = 0.0

    main_wcs_movements_list = []
    reflected_wcs_movements_list = []

    for i in range(len(example_airplane.wings[0].wing_cross_sections)):
        if i == 0:
            wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
                base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[
                    i
                ],
            )
            main_wcs_movements_list.append(wcs_movement)
            reflected_wcs_movements_list.append(wcs_movement)

        else:
            wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
                base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[
                    i
                ],
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(
                    amplitude_x,
                    amplitude_y,
                    amplitude_z,
                ),
                periodAngles_Wcsp_to_Wcs_ixyz=(period_x, period_y, period_z),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(dephase_x, dephase_y, dephase_z),
            )
            main_wcs_movements_list.append(wcs_movement)
            reflected_wcs_movements_list.append(wcs_movement)

    v_tail_root_wcs_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=example_airplane.wings[2].wing_cross_sections[0],
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )
    v_tail_tip_wcs_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=example_airplane.wings[2].wing_cross_sections[1],
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )
    dephase = 169.0

    main_wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
        base_wing=example_airplane.wings[0],
        wing_cross_section_movements=main_wcs_movements_list,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(dephase, 0.0, 0.0),
    )

    reflected_main_wing_movement = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=example_airplane.wings[1],
            wing_cross_section_movements=reflected_wcs_movements_list,
            ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
            periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
            spacingLer_Gs_Cgs=("sine", "sine", "sine"),
            phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(dephase, 0.0, 0.0),
        )
    )

    v_tail_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=example_airplane.wings[2],
        wing_cross_section_movements=[
            v_tail_root_wcs_movement,
            v_tail_tip_wcs_movement,
        ],
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    example_airplane_movement = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=example_airplane,
            wing_movements=[
                main_wing_movement,
                reflected_main_wing_movement,
                v_tail_wing_movement,
            ],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    example_operating_point = ps.operating_point.OperatingPoint(
        rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
    )

    example_operating_point_movement = ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
        base_operating_point=example_operating_point,
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
    )

    example_movement = ps.movements.aeroelastic_movement.AeroelasticMovement(
        airplane_movements=[example_airplane_movement],
        operating_point_movement=example_operating_point_movement,
        delta_time=0.03,
        num_steps=100,
    )

    example_problem = ps.problems.AeroelasticUnsteadyProblem(
        movement=example_movement,
        wing_density=wing_density,
        spring_constant=spring_constant,
        damping_constant=damping_constant,
        aero_scaling=1.0,
        step_discards=5,
        moment_scaling_factor=1.0,
        plot_flap_cycle=False,
    )

    example_solver = ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver(
        aeroelastic_unsteady_problem=example_problem,
    )

    example_solver.run(
        prescribed_wake=True,
    )

    problem = example_solver.unsteady_problem

    # ps.output.animate(
    #     unsteady_solver=example_solver,
    #     scalar_type="lift",
    #     show_wake_vortices=True,
    #     save=True,
    # )
    return problem.net_data_per_wing[0], problem


# Determine which parameter is being swept
active_sweeps = sum(1 for v in (K_VALUES, B_VALUES, DENSITY_VALUES) if v)
if active_sweeps > 1:
    raise ValueError(
        "Only one of K_VALUES, B_VALUES, or DENSITY_VALUES should be non-empty."
    )
if active_sweeps == 0:
    raise ValueError(
        "At least one of K_VALUES, B_VALUES, or DENSITY_VALUES must be non-empty."
    )

if K_VALUES:
    sweep_values = K_VALUES
    sweep_name = "Spring Constant"
    sweep_symbol = "k"
    sweep_kwarg = "spring_constant"
elif B_VALUES:
    sweep_values = B_VALUES
    sweep_name = "Damping Constant"
    sweep_symbol = "b"
    sweep_kwarg = "damping_constant"
else:
    sweep_values = DENSITY_VALUES
    sweep_name = "Wing Density"
    sweep_symbol = "rho"
    sweep_kwarg = "wing_density"

# Run for each swept value and collect Curve 16 of the Net Deformation
results = {}
flap_angle = None
for val in sweep_values:
    print(f"Running with {sweep_symbol}={val}...")
    net_data, problem = run_aeroelastic(**{sweep_kwarg: val})
    # Extract y-component (torsional angle) for Curve 16 across all time steps
    curve_16 = np.array(net_data)[:, CURVE_INDEX, 1].tolist()
    results[val] = curve_16
    print(f"  Completed {sweep_symbol}={val}, {len(curve_16)} steps")

    # Compute the prescribed flap angle once (it is the same for all runs)
    if flap_angle is None:
        wm = problem.wing_movement
        amp = wm.ampAngles_Gs_to_Wn_ixyz[0]
        period = wm.periodAngles_Gs_to_Wn_ixyz[0]
        phase = np.deg2rad(wm.phaseAngles_Gs_to_Wn_ixyz[0])
        dt = problem.movement.delta_time
        num_steps = len(curve_16)
        t = np.arange(num_steps) * dt
        flap_angle = (amp * np.sin((2 * np.pi / period) * t + phase)).tolist()

# Overlay plot of Curve 16 Net Deformation for all swept values
plt.figure(figsize=(12, 6), dpi=200)
for val, curve in results.items():
    plt.plot(range(len(curve)), curve, label=f"{sweep_symbol}={val}")
plt.plot(
    range(len(flap_angle)),
    flap_angle,
    label="Flap Position",
    color="black",
    linestyle="--",
)
plt.xlabel("Step")
plt.ylabel("Angle (degrees)")
plt.title(f"Net Deformation (Curve {CURVE_INDEX}) - Varying {sweep_name}")
plt.legend()
plt.grid(True)
filename = f"Net_Deformation_Curve_{CURVE_INDEX}_{sweep_name.replace(' ', '_')}.png"
plt.savefig(filename)
plt.show()
