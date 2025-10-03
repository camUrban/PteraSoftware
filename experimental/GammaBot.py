# Import Python’s math and copy packages.
import copy

# Import Numpy.
import numpy as np

import pterasoftware as ps
import dxf_to_csv

gammabot_velocity = 0.9
gammabot_alpha = 0.0
gammabot_flapping_frequency = 165.0
wing_spacing = 0.02172011
gammabot_flapping_amplitude_angleX = 34.0
gammabot_flapping_amplitude_angleY = 43.0
# TODO: Get actual value for the phase
delta = -30
gammabot_flapping_phase_angleY = 0 - delta

gammabot_flapping_period = 1.0 / gammabot_flapping_frequency

# Set the number of flap cycles to run the simulation for. The converged result is X
# flaps. Set the number of chordwise panels. The converged result is X panels. Set
# the number of sections to map on each Wing. There will be this number +1
# WingCrossSections per Wing. The converged result is X spanwise sections. Set the
# chordwise spacing scheme for the panels. This is set to uniform, as is standard for
# UVLM simulations.
num_flaps = 5
num_chordwise_panels = 5
num_spanwise_sections = 20
chordwise_spacing = "uniform"

num_wing_cross_sections = num_spanwise_sections + 1

wing_section_data = dxf_to_csv.process_dxf_to_wing_section_data(
    "gammabot_approximate_wing.dxf", num_spanwise_sections
)

# Define an empty list to hold the WingCrossSections.
gammabot_airplane_wing_cross_sections = []

# Iterate through the wing section data to create the WingCrossSections.
for i in range(num_wing_cross_sections):
    if i < (num_wing_cross_sections - 1):
        this_num_spanwise_panels = 1
    else:
        this_num_spanwise_panels = None

    this_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        Lp_Wcsp_Lpp=wing_section_data[i, :3],
        angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        chord=wing_section_data[i, 3],
        airfoil=ps.geometry.airfoil.Airfoil(
            name="naca0012",
        ),
        num_spanwise_panels=this_num_spanwise_panels,
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
    )

    # Append this WingCrossSection to the list of WingCrossSections.
    gammabot_airplane_wing_cross_sections.append(this_wing_cross_section)

# Define the GammaBot Airplane.
gammabot_airplane = ps.geometry.airplane.Airplane(
    name="GammaBot Airplane",
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=gammabot_airplane_wing_cross_sections,
            prelimLer_G_Cg=(0.0, wing_spacing / 2, 0.0),
            symmetric=True,
            symmetry_normal_G=(0, 1, 0),
            symmetry_point_G_Cg=(0, -wing_spacing / 2, 0),
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
        ),
    ],
)

# Delete the extraneous pointer.
del gammabot_airplane_wing_cross_sections


# TODO: Update this with the actual custom function.
def gammabot_angleX_function(thetaRad):
    """GammaBot's custom angleX function."""
    return np.sin(thetaRad)


# TODO: Update this with the actual custom function.
def gammabot_angleY_function(thetaRad):
    """GammaBot's custom angleY function."""
    return np.sin(thetaRad)


# Initialize an empty list to hold each WingCrossSectionMovement.
gammabot_wing_cross_section_movements = []

# Iterate through each of the WingCrossSections.
for j in range(num_wing_cross_sections):
    # Define the WingCrossSectionMovement for this WingCrossSections. The
    # amplitude and period are both set to one because the true amplitude and period
    # are already accounted for in the custom sweep function.
    this_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=gammabot_airplane.wings[0].wing_cross_sections[j],
        )
    )

    # Append this WingCrossSectionMovement to the list of WingCrossSectionMovements.
    gammabot_wing_cross_section_movements.append(this_wing_cross_section_movement)

# Define the WingMovement that contains the WingCrossSectionMovements.
gammabot_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=gammabot_airplane.wings[0],
    wing_cross_section_movements=gammabot_wing_cross_section_movements,
    ampAngles_G_to_prelimWn_izyx=(
        gammabot_flapping_amplitude_angleX,
        gammabot_flapping_amplitude_angleY,
        0.0,
    ),
    periodAngles_G_to_prelimWn_izyx=(
        gammabot_flapping_period,
        gammabot_flapping_period,
        0.0,
    ),
    spacingAngles_G_to_prelimWn_izyx=(
        gammabot_angleX_function,
        gammabot_angleY_function,
        "sine",
    ),
    phaseAngles_G_to_prelimWn_izyx=(0.0, gammabot_flapping_phase_angleY, 0.0),
)

# Define the WingMovement for the mirrored Wing.
gammabot_mirrored_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=gammabot_airplane.wings[1],
    wing_cross_section_movements=copy.deepcopy(gammabot_wing_cross_section_movements),
    ampAngles_G_to_prelimWn_izyx=(
        gammabot_flapping_amplitude_angleX,
        gammabot_flapping_amplitude_angleY,
        0.0,
    ),
    periodAngles_G_to_prelimWn_izyx=(
        gammabot_flapping_period,
        gammabot_flapping_period,
        0.0,
    ),
    spacingAngles_G_to_prelimWn_izyx=(
        gammabot_angleX_function,
        gammabot_angleY_function,
        "sine",
    ),
    phaseAngles_G_to_prelimWn_izyx=(0.0, gammabot_flapping_phase_angleY, 0.0),
)

# Delete the extraneous pointer.
del gammabot_wing_cross_section_movements

# Define the AirplaneMovement that contains the WingMovement.
gammabot_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=gammabot_airplane,
    wing_movements=[
        gammabot_main_wing_movement,
        gammabot_mirrored_wing_movement,
    ],
)

# Delete the extraneous pointers.
del gammabot_airplane
del gammabot_main_wing_movement
del gammabot_mirrored_wing_movement

# Define an OperatingPoint corresponding to the conditions of the GammaBot study.
gammabot_operating_point = ps.operating_point.OperatingPoint(
    vCg__E=gammabot_velocity, alpha=gammabot_alpha
)

# Define an OperatingPointMovement that contains the OperatingPoint.
gammabot_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=gammabot_operating_point,
    )
)

# Delete the extraneous pointer.
del gammabot_operating_point

# Define the overall Movement.
gammabot_movement = ps.movements.movement.Movement(
    airplane_movements=[gammabot_airplane_movement],
    operating_point_movement=gammabot_operating_point_movement,
    num_cycles=num_flaps,
    delta_time=gammabot_flapping_period / 40.0,
)

# Delete the extraneous pointers.
del gammabot_airplane_movement
del gammabot_operating_point_movement

# Define the GammaBot UnsteadyProblem.
gammabot_problem = ps.problems.UnsteadyProblem(
    movement=gammabot_movement,
)

# Delete the extraneous pointer.
del gammabot_movement

# Define the GammaBot solver.
gammabot_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=gammabot_problem,
    )
)

# Run the GammaBot solver. This study was run using an X wake.
gammabot_solver.run(prescribed_wake=True)

ps.output.draw(
    solver=gammabot_solver,
    show_wake_vortices=True,
    scalar_type="lift",
    save=False,
)

# ps.output.plot_results_versus_time(
#     unsteady_solver=gammabot_solver,
#     show=False,
#     save=True,
# )

ps.output.animate(
    unsteady_solver=gammabot_solver,
    show_wake_vortices=True,
    scalar_type="lift",
    save=True,
)

ps.output.print_results(solver=gammabot_solver)
