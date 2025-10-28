# TODO: Redo the convergence analysis for this simulation and update the relevant
#  parameters if necessary.
"""This script runs a validation case of Ptera Software’s UVLM.

I first emulate the geometry and kinematics of a flapping robotic test stand from
"Experimental and Analytical Pressure Characterization of a Rigid Flapping Wing for
Ornithopter Development" by Derrick Yeo, Ella M. Atkins, and Wei Shyy. Then,
I run the UVLM simulation of an experiment from this paper. Finally, I compare the
simulated results to the published experimental results.

WebPlotDigitizer, by Ankit Rohatgi, was used to extract data from Yeo et al., 2011.

More information can be found in my accompanying report: "Validating an Open-Source
UVLM Solver for Analyzing Flapping Wing Flight: An Experimental Approach." """

# Import Python’s math package.
import math

# Import Numpy and MatPlotLib’s PyPlot package.
import matplotlib.pyplot as plt
import numpy as np

# Import the source package.
import pterasoftware as ps

# Set the given characteristics of the wing in meters.
half_span = 0.213
chord = 0.072

# Set the given forward flight velocity in meters per second.
validation_velocity = 2.9

# Set the given angle of attack in degrees. Note: If you analyze a different
# operating point where this is not zero, you need to modify the code to rotate the
# experimental lift into the wind axes.
validation_alpha = 0

# Set the given flapping frequency in Hertz.
validation_flapping_frequency = 3.3

# This wing planform has a rounded tip so the outermost WingCrossSection needs to
# be inset some amount. This value is in meters.
tip_inset = 0.005

# A similar constraint is that Ptera Software requires symmetric, flapping Wings have
# some small midline offset. This value is in meters.
wing_midline_offset = 0.005

# Import the extracted points from the paper’s diagram of the planform. The resulting
# array is of the form [spanwise coordinate, chordwise coordinate], and is ordered
# from the leading edge root, to the tip, to the trailing edge root. The origin is
# the trailing edge root point. The positive spanwise axis extends from root to tip
# and the positive chordwise axis from trailing edge to leading edge. The values
# are in millimeters. I'll call this the Yeo axis system.
stackPlanformPointsMm_Yeo_Ter = np.genfromtxt(
    "extracted_planform_coordinates.csv", delimiter=","
)

# Convert the points to SI units.
stackPlanformPoints_Yeo_Ter = stackPlanformPointsMm_Yeo_Ter / 1000

# Set the origin to the leading edge root point.
stackPlanformPoints_Yeo_Ler = stackPlanformPoints_Yeo_Ter - np.array(
    [0, chord], dtype=float
)

# Switch the sign of the points' chordwise components.
stackPlanformPoints_YeoXReversed_Ler = stackPlanformPoints_Yeo_Ler * np.array(
    [1, -1], dtype=float
)

# Swap the axes to the form [chordwise coordinate, spanwise coordinate]. The
# coordinates are now in wing axes projected onto its xy-plane, and relative to the
# leading edge root point.
stackPlanformPointsXY_Wn_Ler = stackPlanformPoints_YeoXReversed_Ler[:, [1, 0]]

# Find the index of the point where the planform point's x component equals the half
# span.
tip_index = np.where(stackPlanformPointsXY_Wn_Ler[:, 1] == half_span)[0][0]

# Using the tip index, split the points into two ndarrays of leading and trailing
# edge points (in wing axes projected onto its xy-plane, relative to the leading edge
# root point).
stackLeadingPointsXY_Wn_Ler = stackPlanformPointsXY_Wn_Ler[:tip_index, :]
stackTrailingPointsXY_Wn_Ler = np.flip(
    stackPlanformPointsXY_Wn_Ler[tip_index:, :], axis=0
)

# Set the number of flap cycles to run the simulation for. The converged result is 3
# flaps.
num_flaps = 3

# Set the number of chordwise Panels. The converged result is 5 Panels.
num_chordwise_panels = 5

# Set the number of sections to map on each Wing half. There will be this number +1
# WingCrossSections per Wing half. The converged result is 18 spanwise sections.
num_spanwise_sections = 18

# Set the chordwise spacing scheme for the Panels. This is set to uniform,
# as is standard for UVLM simulations.
chordwise_spacing = "uniform"

# Calculate the spanwise distance between the WingCrossSections.
spanwise_step = (half_span - tip_inset) / num_spanwise_sections

# Define four ndarrays to hold the leading and trailing points of each section’s left
# and right WingCrossSections (in wing axes projected onto its xy-plane, relative to
# the leading edge root point).
stackLeftLpsXY_Wn_Ler = np.zeros((num_spanwise_sections, 2), dtype=float)
stackRightLpsXY_Wn_Ler = np.zeros((num_spanwise_sections, 2), dtype=float)
stackLeftTpsXY_Wn_Ler = np.zeros((num_spanwise_sections, 2), dtype=float)
stackRightTpsXY_Wn_Ler = np.zeros((num_spanwise_sections, 2), dtype=float)

# Iterate through the locations of the future sections to populate the left and right
# WingCrossSection's leading and trailing points (in wing axes projected onto its
# xy-plane, relative to the leading edge root point).
for spanwise_loc in range(num_spanwise_sections):
    # Find the y component of the leading and trailing points (in wing axes projected
    # onto its xy-plane, relative to the leading edge root point).
    stackLeftLpsXY_Wn_Ler[spanwise_loc, 1] = spanwise_loc * spanwise_step
    stackLeftTpsXY_Wn_Ler[spanwise_loc, 1] = spanwise_loc * spanwise_step
    stackRightLpsXY_Wn_Ler[spanwise_loc, 1] = (spanwise_loc + 1) * spanwise_step
    stackRightTpsXY_Wn_Ler[spanwise_loc, 1] = (spanwise_loc + 1) * spanwise_step

    # Interpolate between the points to find their x components (in wing axes
    # projected onto its xy-plane, relative to the leading edge root point).
    stackLeftLpsXY_Wn_Ler[spanwise_loc, 0] = np.interp(
        spanwise_loc * spanwise_step,
        stackLeadingPointsXY_Wn_Ler[:, 1],
        stackLeadingPointsXY_Wn_Ler[:, 0],
    )
    stackLeftTpsXY_Wn_Ler[spanwise_loc, 0] = np.interp(
        spanwise_loc * spanwise_step,
        stackTrailingPointsXY_Wn_Ler[:, 1],
        stackTrailingPointsXY_Wn_Ler[:, 0],
    )
    stackRightLpsXY_Wn_Ler[spanwise_loc, 0] = np.interp(
        (spanwise_loc + 1) * spanwise_step,
        stackLeadingPointsXY_Wn_Ler[:, 1],
        stackLeadingPointsXY_Wn_Ler[:, 0],
    )
    stackRightTpsXY_Wn_Ler[spanwise_loc, 0] = np.interp(
        (spanwise_loc + 1) * spanwise_step,
        stackTrailingPointsXY_Wn_Ler[:, 1],
        stackTrailingPointsXY_Wn_Ler[:, 0],
    )

# Define an empty list to hold the WingCrossSections.
validation_airplane_wing_cross_sections = []

# Iterate through the leading and trailing point ndarrays to create the
# WingCrossSections.
for i in range(num_spanwise_sections):
    if i == 0:
        thisLpY_Wcsp_Lpp = 0.0
        thisLpX_Wcsp_Lpp = 0.0
    else:
        thisLpY_Wcsp_Lpp = spanwise_step
        thisLpX_Wcsp_Lpp = stackLeftLpsXY_Wn_Ler[i, 0] - stackLeftLpsXY_Wn_Ler[i - 1, 0]

    this_chord = stackLeftTpsXY_Wn_Ler[i, 0] - stackLeftLpsXY_Wn_Ler[i, 0]

    # Create this WingCrossSection.
    this_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=ps.geometry.airfoil.Airfoil(
            name="naca0012",
        ),
        num_spanwise_panels=1,
        chord=this_chord,
        Lp_Wcsp_Lpp=(thisLpX_Wcsp_Lpp, thisLpY_Wcsp_Lpp, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing="uniform",
    )

    # Append this WingCrossSection to the list of WingCrossSections.
    validation_airplane_wing_cross_sections.append(this_wing_cross_section)

    # If this is the last section, also create the right WingCrossSection and append
    # it to the list.
    if i == num_spanwise_sections - 1:
        thisLpY_Wcsp_Lpp = spanwise_step
        thisLpX_Wcsp_Lpp = (
            stackRightLpsXY_Wn_Ler[i, 0] - stackRightLpsXY_Wn_Ler[i - 1, 0]
        )

        this_chord = stackRightTpsXY_Wn_Ler[i, 0] - stackRightLpsXY_Wn_Ler[i, 0]

        this_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=ps.geometry.airfoil.Airfoil(
                name="naca0012",
            ),
            num_spanwise_panels=None,
            chord=this_chord,
            Lp_Wcsp_Lpp=(thisLpX_Wcsp_Lpp, thisLpY_Wcsp_Lpp, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            control_surface_symmetry_type="symmetric",
            control_surface_hinge_point=0.75,
            control_surface_deflection=0.0,
            spanwise_spacing=None,
        )

        validation_airplane_wing_cross_sections.append(this_wing_cross_section)

# Create the Airplane.
validation_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=validation_airplane_wing_cross_sections,
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, wing_midline_offset / 2, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
        ),
    ],
    name="Validation Airplane",
)

# Delete the extraneous pointer.
del validation_airplane_wing_cross_sections

# Initialize empty lists to hold the WingCrossSectionMovements for the main and
# reflected main Wings.
main_wing_cross_section_movements = []
reflected_main_wing_cross_section_movements = []

# Create static WingCrossSectionMovements for each WingCrossSection in the main and
# reflected main Wings.
for i in range(num_spanwise_sections + 1):
    this_main_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=validation_airplane.wings[0].wing_cross_sections[i]
        )
    )
    this_reflected_main_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=validation_airplane.wings[1].wing_cross_sections[i]
        )
    )

    main_wing_cross_section_movements.append(this_main_wing_cross_section_movement)
    reflected_main_wing_cross_section_movements.append(
        this_reflected_main_wing_cross_section_movement
    )

# Delete the extraneous pointers to make debugging easier.
del this_main_wing_cross_section_movement
del this_reflected_main_wing_cross_section_movement


def validation_geometry_sweep_function(time):
    """This function takes in the time during a flap cycle and returns the flap angle
    in degrees. It uses the flapping frequency defined in the encompassing script,
    and is based on a fourth-order Fourier series. The coefficients were calculated
    by Yeo et al., 2011.

    :param time: float or a (N,) ndarray of floats

        This is a single time or a ndarray of N times at which to calculate the flap
        angle. The units are seconds.

    :return flap_angle: float a (N,) ndarray of floats

        This is a single flap angle or a ndarray of N flap angles at the inputted
        time value or values. The units are degrees.
    """

    # Set the Fourier series coefficients and the flapping frequency.
    a_0 = 0.0354
    a_1 = 4.10e-5
    b_1 = 0.3793
    a_2 = -0.0322
    b_2 = -1.95e-6
    a_3 = -8.90e-7
    b_3 = -0.0035
    a_4 = 0.00046
    b_4 = -3.60e-6
    f = 2 * math.pi * validation_flapping_frequency

    # Calculate and return the flap angle(s).
    return (
        a_0
        + a_1 * np.cos(1 * f * time)
        + b_1 * np.sin(1 * f * time)
        + a_2 * np.cos(2 * f * time)
        + b_2 * np.sin(2 * f * time)
        + a_3 * np.cos(3 * f * time)
        + b_3 * np.sin(3 * f * time)
        + a_4 * np.cos(4 * f * time)
        + b_4 * np.sin(4 * f * time)
    ) / 0.0174533


def time_normalized_validation_geometry_sweep_function_rad(time):
    """This function takes in the time during a flap cycle and returns the flap angle
    in radians. It uses a normalized flapping frequency of 1 Hz, and is based on a
    fourth-order Fourier series. The coefficients were calculated by Yeo et al., 2011.

    :param time: float or a (N,) ndarray of floats

        This is a single time or a ndarray of N times at which to calculate the flap
        angle. The units are seconds.

    :return flap_angle: float or a (N,) ndarray of floats

        This is a single flap angle or a ndarray of N flap angles at the inputted
        time value or values. The units are radians.
    """

    # Set the Fourier series coefficients.
    a_0 = 0.0354
    a_1 = 4.10e-5
    b_1 = 0.3793
    a_2 = -0.0322
    b_2 = -1.95e-6
    a_3 = -8.90e-7
    b_3 = -0.0035
    a_4 = 0.00046
    b_4 = -3.60e-6

    # Calculate and return the flap angle(s).
    return -(
        a_0
        + a_1 * np.cos(1 * time)
        + b_1 * np.sin(1 * time)
        + a_2 * np.cos(2 * time)
        + b_2 * np.sin(2 * time)
        + a_3 * np.cos(3 * time)
        + b_3 * np.sin(3 * time)
        + a_4 * np.cos(4 * time)
        + b_4 * np.sin(4 * time)
    )


# Define the WingMovements for the main and reflected main Wings.
main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=validation_airplane.wings[0],
    wing_cross_section_movements=main_wing_cross_section_movements,
    # TODO: Replace with actual angle movement values.
    ampAngles_Gs_to_Wn_ixyz=(22.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1 / validation_flapping_frequency, 0.0, 0.0),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
)
reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=validation_airplane.wings[1],
    wing_cross_section_movements=main_wing_cross_section_movements,
    # TODO: Replace with actual angle movement values.
    ampAngles_Gs_to_Wn_ixyz=(22.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1 / validation_flapping_frequency, 0.0, 0.0),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
)

# Delete the extraneous pointer.
del main_wing_cross_section_movements

# Define the AirplaneMovement that contains the WingMovements.
validation_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=validation_airplane,
    wing_movements=[main_wing_movement, reflected_main_wing_movement],
)

# Delete the extraneous pointers.
del validation_airplane
del main_wing_movement
del reflected_main_wing_movement

# Define an OperatingPoint and OperatingPointMovement corresponding to the conditions
# of the validation study.
validation_operating_point = ps.operating_point.OperatingPoint(
    vCg__E=validation_velocity, alpha=validation_alpha
)
validation_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=validation_operating_point
    )
)

# Define the Movement.
validation_movement = ps.movements.movement.Movement(
    airplane_movements=[validation_airplane_movement],
    operating_point_movement=validation_operating_point_movement,
    num_cycles=num_flaps,
)

# Delete the extraneous pointers.
del validation_airplane_movement
del validation_operating_point_movement

# Define the UnsteadyProblem.
validation_problem = ps.problems.UnsteadyProblem(
    movement=validation_movement,
    only_final_results=False,
)

# Define the UnsteadyRingVortexLatticeMethodSolver.
validation_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=validation_problem,
    )
)

# Delete the extraneous pointer.
del validation_problem

# Define the position of the points of interest and the area of their rectangles.
# These values were extracted by digitizing the figures in Yeo et al., 2011.
blueTrailingPointsXY_Wn_Ler = [0.060, 0.036]
blue_trailing_area = 0.072 * 0.024
blueMiddlePointsXY_Wn_Ler = [0.036, 0.036]
blue_middle_area = 0.072 * 0.024
blueLeadingPointsXY_Wn_Ler = [0.012, 0.036]
blue_leading_area = 0.072 * 0.024
orangeTrailingPointsXY_Wn_Ler = [0.05532, 0.107]
orange_trailing_area = 0.07 * 0.02112
orangeMiddlePointsXY_Wn_Ler = [0.0342, 0.107]
orange_middle_area = 0.07 * 0.02112
orangeLeadingPointsXY_Wn_Ler = [0.01308, 0.107]
orange_leading_area = 0.07 * 0.02112
greenTrailingPointsXY_Wn_Ler = [0.04569, 0.162825]
green_trailing_area = 0.04165 * 0.015
greenMiddlePointsXY_Wn_Ler = [0.03069, 0.176]
green_middle_area = 0.06565 * 0.015
greenLeadingPointsXY_Wn_Ler = [0.01569, 0.1775]
green_leading_area = 0.071 * 0.015

# Run the validation solver using a prescribed wake.
validation_solver.run(prescribed_wake=True)

# Extract the Movement's num_steps and delta_time attributes.
validation_num_steps = validation_movement.num_steps
validation_delta_time = validation_movement.delta_time

# Create a variable to hold the time in seconds at each of the simulation’s time steps.
times = np.linspace(
    0,
    validation_num_steps * validation_delta_time,
    validation_num_steps,
    endpoint=False,
)

# Discretize the time period of the final flap analyzed into 100 steps. Store this to
# a ndarray.
final_flap_times = np.linspace(
    (num_flaps - 1) / validation_flapping_frequency,
    num_flaps / validation_flapping_frequency,
    100,
    endpoint=False,
)

# Discretize the normalized flap cycle times into 100 steps. Store this to a ndarray.
normalized_times = np.linspace(0, 1, 100, endpoint=False)

# Pull the experimental pressure vs. time histories from the digitized data. These
# data sets are stored in CSV files in the same directory as this script. The
# pressure units used are inAq and time units are normalized flap cycle times from 0
# to 1.
exp_blue_trailing_point_pressures = np.genfromtxt(
    "blue_trailing_point_experimental_pressures.csv", delimiter=","
)
exp_blue_middle_point_pressures = np.genfromtxt(
    "blue_middle_point_experimental_pressures.csv", delimiter=","
)
exp_blue_leading_point_pressures = np.genfromtxt(
    "blue_leading_point_experimental_pressures.csv", delimiter=","
)
exp_orange_trailing_point_pressures = np.genfromtxt(
    "orange_trailing_point_experimental_pressures.csv", delimiter=","
)
exp_orange_middle_point_pressures = np.genfromtxt(
    "orange_middle_point_experimental_pressures.csv", delimiter=","
)
exp_orange_leading_point_pressures = np.genfromtxt(
    "orange_leading_point_experimental_pressures.csv", delimiter=","
)
exp_green_trailing_point_pressures = np.genfromtxt(
    "green_trailing_point_experimental_pressures.csv", delimiter=","
)
exp_green_middle_point_pressures = np.genfromtxt(
    "green_middle_point_experimental_pressures.csv", delimiter=","
)
exp_green_leading_point_pressures = np.genfromtxt(
    "green_leading_point_experimental_pressures.csv", delimiter=","
)

# Interpolate the experimental pressure data to ensure that they all reference the
# same normalized timescale.
exp_blue_trailing_point_pressures_norm = np.interp(
    normalized_times,
    exp_blue_trailing_point_pressures[:, 0],
    exp_blue_trailing_point_pressures[:, 1],
)
exp_blue_middle_point_pressures_norm = np.interp(
    normalized_times,
    exp_blue_middle_point_pressures[:, 0],
    exp_blue_middle_point_pressures[:, 1],
)
exp_blue_leading_point_pressures_norm = np.interp(
    normalized_times,
    exp_blue_leading_point_pressures[:, 0],
    exp_blue_leading_point_pressures[:, 1],
)
exp_orange_trailing_point_pressures_norm = np.interp(
    normalized_times,
    exp_orange_trailing_point_pressures[:, 0],
    exp_orange_trailing_point_pressures[:, 1],
)
exp_orange_middle_point_pressures_norm = np.interp(
    normalized_times,
    exp_orange_middle_point_pressures[:, 0],
    exp_orange_middle_point_pressures[:, 1],
)
exp_orange_leading_point_pressures_norm = np.interp(
    normalized_times,
    exp_orange_leading_point_pressures[:, 0],
    exp_orange_leading_point_pressures[:, 1],
)
exp_green_trailing_point_pressures_norm = np.interp(
    normalized_times,
    exp_green_trailing_point_pressures[:, 0],
    exp_green_trailing_point_pressures[:, 1],
)
exp_green_middle_point_pressures_norm = np.interp(
    normalized_times,
    exp_green_middle_point_pressures[:, 0],
    exp_green_middle_point_pressures[:, 1],
)
exp_green_leading_point_pressures_norm = np.interp(
    normalized_times,
    exp_green_leading_point_pressures[:, 0],
    exp_green_leading_point_pressures[:, 1],
)

# Find the normal force time history on each of the experimental panels in Newtons.
exp_blue_trailing_normal_forces = (
    248.84 * exp_blue_trailing_point_pressures_norm * blue_trailing_area
)
exp_blue_middle_normal_forces = (
    248.84 * exp_blue_middle_point_pressures_norm * blue_middle_area
)
exp_blue_leading_normal_forces = (
    248.84 * exp_blue_leading_point_pressures_norm * blue_leading_area
)
exp_orange_trailing_normal_forces = (
    248.84 * exp_orange_trailing_point_pressures_norm * orange_trailing_area
)
exp_orange_middle_normal_forces = (
    248.84 * exp_orange_middle_point_pressures_norm * orange_middle_area
)
exp_orange_leading_normal_forces = (
    248.84 * exp_orange_leading_point_pressures_norm * orange_leading_area
)
exp_green_trailing_normal_forces = (
    248.84 * exp_green_trailing_point_pressures_norm * green_trailing_area
)
exp_green_middle_normal_forces = (
    248.84 * exp_green_middle_point_pressures_norm * green_middle_area
)
exp_green_leading_normal_forces = (
    248.84 * exp_green_leading_point_pressures_norm * green_leading_area
)

# Convert each experimental panel's normal force time history to a time history of
# the force's geometry axes' z-component.
stackExpBlueTrailingForcesZ_G = exp_blue_trailing_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpBlueMiddleForcesZ_G = exp_blue_middle_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpBlueLeadingForcesZ_G = exp_blue_leading_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpOrangeTrailingForcesZ_G = exp_orange_trailing_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpOrangeMiddleForcesZ_G = exp_orange_middle_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpOrangeLeadingForcesZ_G = exp_orange_leading_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpGreenTrailingForcesZ_G = exp_green_trailing_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpGreenMiddleForcesZ_G = exp_green_middle_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)
stackExpGreenLeadingForcesZ_G = exp_green_leading_normal_forces * np.cos(
    time_normalized_validation_geometry_sweep_function_rad(normalized_times)
)

# Calculate the net experimental force's geometry axes' z-component. This is
# multiplied by two because the experimental panels only cover one of the symmetric
# wing halves.
stackExpNetForcesZ_G = 2 * (
    stackExpBlueTrailingForcesZ_G
    + stackExpBlueMiddleForcesZ_G
    + stackExpBlueLeadingForcesZ_G
    + stackExpOrangeTrailingForcesZ_G
    + stackExpOrangeMiddleForcesZ_G
    + stackExpOrangeLeadingForcesZ_G
    + stackExpGreenTrailingForcesZ_G
    + stackExpGreenMiddleForcesZ_G
    + stackExpGreenLeadingForcesZ_G
)

# Initialize a ndarray to hold the net force's wind axes' z-component.
stackExpNetForcesZ_W = np.zeros(stackExpNetForcesZ_G.size, dtype=float)

# Get the passive transformation matrix which maps in homogeneous coordinates from
# the first Airplane's geometry axes relative to the first Airplane's CG to wind
# axes relative to the first Airplane's CG.
T_pas_GP1_CgP1_to_W_CgP1 = validation_operating_point.T_pas_GP1_CgP1_to_W_CgP1

# Delete the extraneous pointer.
del validation_operating_point

# Transform from the first Airplane's geometry axes to wind axes.
for force_id, expNetForceZ_GP1 in enumerate(stackExpNetForcesZ_G):
    expNetForceHomog_GP1 = np.array([0.0, 0.0, expNetForceZ_GP1, 0.0], dtype=float)
    expNetForceHomog_W = T_pas_GP1_CgP1_to_W_CgP1 @ expNetForceHomog_GP1
    expNetForceZ_W = expNetForceHomog_W[2]
    stackExpNetForcesZ_W[force_id] = expNetForceZ_W

# Get the experimental lift values. Lift is defined as the force's wind axes'
# z-component multiplied by negative one.
exp_lifts = -1 * stackExpNetForcesZ_W

# Get this solver’s SteadyProblems' Airplanes.
airplanes = []
for steady_problem in validation_solver.steady_problems:
    airplanes.append(steady_problem.airplanes[0])

# Initialize a ndarray to hold the force at each time step (in wind axes).
stackSimForces_W = np.zeros((3, validation_num_steps))

# Iterate through the time steps and populate the ndarray.
for step, airplane in enumerate(airplanes):
    stackSimForces_W[:, step] = airplane.forces_W

# Initialize the figure and axes of the experimental versus simulated lift plot.
lift_figure, lift_axes = plt.subplots(figsize=(5, 4))

# Get the simulated lift values. Lift is defined as the force's wind axes' z-component
# multiplied by negative one.
sim_lifts = -1 * stackSimForces_W[2, :]

# Interpolate the simulated lift values to find them with respect to the normalized
# final flap timescale.
final_flap_sim_lifts = np.interp(final_flap_times, times, sim_lifts[:])

sim_lift_color = "#D81E5B"
exp_lift_color = "#003F91"

num_markers = 6
marker_size = 8
text_color = "black"
figure_background_color = "None"

lift_axes.spines.right.set_visible(False)
lift_axes.spines.top.set_visible(False)
lift_axes.spines.bottom.set_color(text_color)
lift_axes.spines.left.set_color(text_color)
lift_axes.xaxis.label.set_color(text_color)
lift_axes.yaxis.label.set_color(text_color)
lift_axes.tick_params(axis="x", colors=text_color)
lift_axes.tick_params(axis="y", colors=text_color)
lift_figure.patch.set_facecolor(figure_background_color)
lift_axes.set_facecolor(figure_background_color)

marker_spacing = 1.0 / num_markers

# Plot the simulated lift values. The x-axis is set to the normalized times,
# which may seem odd because we just interpolated to get them in terms of the
# normalized final flap times. But, they are discretized in exactly the same way as
# the normalized times, just horizontally shifted.
lift_axes.plot(
    normalized_times,
    final_flap_sim_lifts,
    label="Simulated",
    color=sim_lift_color,
    marker=".",
    markevery=(marker_spacing * 0 / 2, marker_spacing),
    markersize=marker_size,
)

# Plot the experimental lift values.
lift_axes.plot(
    normalized_times,
    exp_lifts,
    label="Experimental",
    color=exp_lift_color,
    marker=".",
    markevery=(marker_spacing * 1 / 2, marker_spacing),
    markersize=marker_size,
)

# Add a gray box to signify which part of the graph is the downstroke.
plt.axvspan(0.25, 0.75, facecolor="darkgray", label="Downstroke")

# Label the axis, add a title, and add a legend.
lift_axes.set_xlabel(
    "Normalized Flap Cycle Time",
)
lift_axes.set_ylabel(
    "Lift (N)",
)
lift_axes.set_title(
    "Simulated and Experimental Lift Versus Time",
)
lift_axes.legend(
    loc="upper left",
    facecolor=figure_background_color,
    edgecolor=figure_background_color,
    labelcolor=text_color,
)

# Save the lift comparison figure.
lift_figure.savefig(
    fname="Lift comparison.jpg",
    dpi=300,
    bbox_inches="tight",
)

# Delete the extraneous pointers.
del airplanes
del stackSimForces_W
del step

# Calculate the lift mean absolute error (MAE). The experimental and simulated lift
# comparison here is valid because, due to the interpolation steps, the experimental
# and simulated lifts time histories are discretized so that they are with respect to
# the same timescale.
lift_absolute_errors = np.abs(final_flap_sim_lifts - exp_lifts)
lift_mean_absolute_error = np.mean(lift_absolute_errors)

sim_lift_rms = math.sqrt(np.mean(final_flap_sim_lifts**2))
exp_lift_rms = math.sqrt(np.mean(exp_lifts**2))
lift_rmsape = 100 * abs((sim_lift_rms - exp_lift_rms) / exp_lift_rms)
print("\nLift RMS Absolute Percent Error: " + str(np.round(lift_rmsape, 2)) + "%")
print("Simulated Lift RMS: " + str(np.round(sim_lift_rms, 4)) + " N")
print("Experimental Lift RMS: " + str(np.round(exp_lift_rms, 4)) + " N")

# Print the MAE.
print(
    "\nMean Absolute Error on Lift: " + str(np.round(lift_mean_absolute_error, 4)) + "N"
)

# Calculate the experimental root-mean-square (RMS) lift.
exp_rms_lift = np.sqrt(np.mean(np.power(exp_lifts, 2)))

# Print the experimental RMS lift.
print("Experimental RMS Lift: " + str(np.round(exp_rms_lift, 4)) + " N")

ps.output.draw(
    solver=validation_solver,
    show_wake_vortices=True,
    scalar_type="lift",
    save=True,
)

ps.output.plot_results_versus_time(
    unsteady_solver=validation_solver,
    show=False,
    save=True,
)

ps.output.animate(
    unsteady_solver=validation_solver,
    show_wake_vortices=True,
    scalar_type="lift",
    save=True,
)
