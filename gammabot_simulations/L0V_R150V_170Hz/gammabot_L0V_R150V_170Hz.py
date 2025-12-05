from pathlib import Path

import pterasoftware as ps
from gammabot_simulations import dxf_to_csv

fine_mesh = False

velocity = 0.9
alpha = 0.0
flapping_frequency = 170.0
wing_spacing = 0.02746

# Left Wing
phi_max_left = 0.0
phi_v_shift_left = 0.0
psi_max_left = 0.0
psi_v_shift_left = 0.0
delta_left = 0.0

# Right Wing
phi_max_right = 15.63
phi_v_shift_right = 3.58
psi_max_right = 37.07
psi_v_shift_right = 2.72
delta_right = 15.40

x_offset = 1.215e-6
y_offset = -0.00215

flapping_period = 1.0 / flapping_frequency

# TODO: Find the converged values for these parameters.
# Set the number of flap cycles to run the simulation for. The converged result is X
# flaps. Set the number of chordwise panels. The converged result is X panels. Set
# the number of sections to map on each Wing. There will be this number +1
# WingCrossSections per Wing. The converged result is X spanwise sections. Set the
# chordwise spacing scheme for the panels. This is set to uniform, as is standard for
# UVLM simulations.
prescribed_wake = True
chordwise_spacing = "uniform"

if fine_mesh:
    num_flaps = 3
    num_chordwise_panels = 8
    num_spanwise_sections = 16
    delta_time = 0.000130
else:
    num_flaps = 2
    num_chordwise_panels = 4
    num_spanwise_sections = 8
    delta_time = 0.000272

flapping_amplitude_angleX_left = phi_max_left
flapping_amplitude_angleY_left = psi_max_left

flapping_period_angleX_left = flapping_period
if flapping_amplitude_angleX_left == 0.0:
    flapping_period_angleX_left = 0.0
flapping_period_angleY_left = flapping_period
if flapping_amplitude_angleY_left == 0.0:
    flapping_period_angleY_left = 0.0

flapping_phase_angleY_left = 90.0 + delta_left
if flapping_amplitude_angleY_left == 0.0:
    flapping_phase_angleY_left = 0.0

flapping_amplitude_angleX_right = phi_max_right
flapping_amplitude_angleY_right = psi_max_right

flapping_period_angleX_right = flapping_period
if flapping_amplitude_angleX_right == 0.0:
    flapping_period_angleX_right = 0.0
flapping_period_angleY_right = flapping_period
if flapping_amplitude_angleY_right == 0.0:
    flapping_period_angleY_right = 0.0

flapping_phase_angleY_right = 90.0 + delta_right
if flapping_amplitude_angleY_right == 0.0:
    flapping_phase_angleY_right = 0.0

num_wing_cross_sections = num_spanwise_sections + 1

dxf_filepath = Path(__file__).parent.parent / "gammabot_approximate_wing.dxf"
wing_section_data = dxf_to_csv.process_dxf_to_wing_section_data(
    str(dxf_filepath), num_spanwise_sections
)

left_wing_cross_sections = []
right_wing_cross_sections = []

# Iterate through the wing section data to create the WingCrossSections.
for i in range(num_wing_cross_sections):
    if i < (num_wing_cross_sections - 1):
        this_num_spanwise_panels = 1
    else:
        this_num_spanwise_panels = None

    left_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        Lp_Wcsp_Lpp=wing_section_data[i, :3],
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        chord=wing_section_data[i, 3],
        airfoil=ps.geometry.airfoil.Airfoil(
            name="naca0012",
        ),
        num_spanwise_panels=this_num_spanwise_panels,
        control_surface_symmetry_type=None,
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
    )
    right_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        Lp_Wcsp_Lpp=wing_section_data[i, :3],
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        chord=wing_section_data[i, 3],
        airfoil=ps.geometry.airfoil.Airfoil(
            name="naca0012",
        ),
        num_spanwise_panels=this_num_spanwise_panels,
        control_surface_symmetry_type=None,
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
    )

    # Append this WingCrossSection to the lists of WingCrossSections.
    left_wing_cross_sections.append(left_wing_cross_section)
    right_wing_cross_sections.append(right_wing_cross_section)

# Define the GammaBot Airplane.
airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=left_wing_cross_sections,
            Ler_Gs_Cgs=(0.0, wing_spacing / 2, 0.0),
            angles_Gs_to_Wn_ixyz=(phi_v_shift_left, psi_v_shift_left, 0.0),
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=right_wing_cross_sections,
            Ler_Gs_Cgs=(0.0, wing_spacing / 2, 0.0),
            angles_Gs_to_Wn_ixyz=(phi_v_shift_right, psi_v_shift_right, 0.0),
            symmetric=False,
            mirror_only=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
        ),
    ],
    name="GammaBot",
)

# Delete the extraneous pointers.
del left_wing_cross_sections
del right_wing_cross_sections

left_wing_cross_section_movements = []
for j in range(num_wing_cross_sections):
    movement = ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=airplane.wings[0].wing_cross_sections[j]
    )
    left_wing_cross_section_movements.append(movement)

right_wing_cross_section_movements = []
for j in range(num_wing_cross_sections):
    movement = ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=airplane.wings[1].wing_cross_sections[j]
    )
    right_wing_cross_section_movements.append(movement)

# Define the WingMovement that contains the WingCrossSectionMovements.
left_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=airplane.wings[0],
    wing_cross_section_movements=left_wing_cross_section_movements,
    rotationPointOffset_Gs_Ler=(x_offset, y_offset, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(
        flapping_amplitude_angleX_left,
        flapping_amplitude_angleY_left,
        0.0,
    ),
    periodAngles_Gs_to_Wn_ixyz=(
        flapping_period_angleX_left,
        flapping_period_angleY_left,
        0.0,
    ),
    spacingAngles_Gs_to_Wn_ixyz=(
        "sine",
        "sine",
        "sine",
    ),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, flapping_phase_angleY_left, 0.0),
)

# Define the WingMovement for the right Wing.
right_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=airplane.wings[1],
    wing_cross_section_movements=right_wing_cross_section_movements,
    rotationPointOffset_Gs_Ler=(x_offset, y_offset, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(
        flapping_amplitude_angleX_right,
        flapping_amplitude_angleY_right,
        0.0,
    ),
    periodAngles_Gs_to_Wn_ixyz=(
        flapping_period_angleX_right,
        flapping_period_angleY_right,
        0.0,
    ),
    spacingAngles_Gs_to_Wn_ixyz=(
        "sine",
        "sine",
        "sine",
    ),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, flapping_phase_angleY_right, 0.0),
)

# Delete the extraneous pointers.
del left_wing_cross_section_movements
del right_wing_cross_section_movements

# Define the AirplaneMovement that contains the WingMovement.
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=airplane,
    wing_movements=[
        left_wing_movement,
        right_wing_movement,
    ],
)

# Delete the extraneous pointers.
del airplane
del left_wing_movement
del right_wing_movement

# Define an OperatingPoint corresponding to the conditions of the GammaBot study.
operating_point = ps.operating_point.OperatingPoint(vCg__E=velocity, alpha=alpha)

# Define an OperatingPointMovement that contains the OperatingPoint.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=operating_point,
)

# Delete the extraneous pointer.
del operating_point

# Define the overall Movement.
movement = ps.movements.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    num_cycles=num_flaps,
    delta_time=delta_time,
)

# Delete the extraneous pointers.
del airplane_movement
del operating_point_movement

# Define the GammaBot UnsteadyProblem.
problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Delete the extraneous pointer.
del movement

# Define the GammaBot solver.
solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    unsteady_problem=problem,
)

solver.run(prescribed_wake=prescribed_wake)

ps.output.draw(
    solver=solver,
    show_wake_vortices=True,
    scalar_type="induced drag",
    save=True,
)

ps.output.plot_results_versus_time(
    unsteady_solver=solver,
    show=True,
    save=True,
)

ps.output.animate(
    unsteady_solver=solver,
    show_wake_vortices=True,
    scalar_type="induced drag",
    save=True,
)

ps.output.print_results(solver=solver)
