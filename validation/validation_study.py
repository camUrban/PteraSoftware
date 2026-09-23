"""This script runs a validation case of Ptera Software's UVLM.

I first emulate the geometry and kinematics of a flapping robotic test stand from
"Experimental and Analytical Pressure Characterization of a Rigid Flapping Wing for
Ornithopter Development" by Derrick Yeo, Ella M. Atkins, and Wei Shyy. Then, I run the
UVLM simulation of an experiment from this paper. Finally, I compare the simulated
results to the published experimental results.

WebPlotDigitizer, by Ankit Rohatgi, was used to extract data from Yeo et al., 2011.

More information can be found in my accompanying report: "Validating an Open-Source UVLM
Solver for Analyzing Flapping Wing Flight: An Experimental Approach."
"""

# Import Python's math and pathlib packages.
import math
from pathlib import Path

# Import NumPy and MatPlotLib's PyPlot package.
import matplotlib.pyplot as plt
import numpy as np

# Import the source package.
import pterasoftware as ps

# Find this script's directory so that the data files it reads and the figure it saves
# resolve correctly regardless of the current working directory. The experimental data
# extracted from the paper is stored in CSV files in a subdirectory.
VALIDATION_DIRECTORY = Path(__file__).resolve().parent
EXPERIMENTAL_DATA_DIRECTORY = VALIDATION_DIRECTORY / "experimental_data"

# Configure logging to display info level messages on the console alongside progress
# bars. Keep the configured logger so this script can log its own results alongside the
# package's messages.
validation_logger = ps.set_up_logging(level="Info")

# Set the given characteristics of the wing in meters.
HALF_SPAN = 0.213
CHORD = 0.072

# Set the given forward flight velocity in meters per second.
VALIDATION_VELOCITY = 2.9

# Set the given angle of attack in degrees. If you analyze a different operating point
# where this is not zero, you need to modify the code to rotate the experimental lift
# into the wind axes.
VALIDATION_ALPHA = 0

# Set the given flapping frequency in Hertz.
VALIDATION_FLAPPING_FREQUENCY = 3.3

# This wing planform has a rounded tip so the outermost WingCrossSection needs to be
# inset some amount. This value is in meters.
TIP_INSET = 0.005

# A similar constraint is that Ptera Software requires symmetric, flapping Wings have
# some small midline offset. This value is in meters.
WING_MIDLINE_OFFSET = 0.005

# Import the extracted points from the paper's diagram of the planform. The resulting
# array is of the form [spanwise coordinate, chordwise coordinate], and is ordered from
# the leading edge root, to the tip, to the trailing edge root. The origin is the
# trailing edge root point. The positive spanwise axis extends from root to tip and the
# positive chordwise axis from trailing edge to leading edge. The values are in
# millimeters. I'll call this the Yeo axis system.
stackPlanformPointsMm_Yeo_Ter = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "extracted_planform_coordinates.csv", delimiter=","
)

# Convert the points to SI units.
stackPlanformPoints_Yeo_Ter = stackPlanformPointsMm_Yeo_Ter / 1000

# Set the origin to the leading edge root point.
stackPlanformPoints_Yeo_Ler = stackPlanformPoints_Yeo_Ter - np.array(
    [0, CHORD], dtype=float
)

# Switch the sign of the points' chordwise components.
stackPlanformPoints_YeoXReversed_Ler = stackPlanformPoints_Yeo_Ler * np.array(
    [1, -1], dtype=float
)

# Swap the axes to the form [chordwise coordinate, spanwise coordinate]. The coordinates
# are now in wing axes projected onto its xy plane, and relative to the leading edge
# root point.
stackPlanformPointsXY_Wn_Ler = stackPlanformPoints_YeoXReversed_Ler[:, [1, 0]]

# Find the index of the point where the planform point's y component equals the half
# span.
tip_index = np.where(stackPlanformPointsXY_Wn_Ler[:, 1] == HALF_SPAN)[0][0]

# Using the tip index, split the points into two ndarrays of leading and trailing edge
# points (in wing axes projected onto its xy plane, relative to the leading edge root
# point). Both curves include the tip point so they span the same maximum y component.
stackLeadingPointsXY_Wn_Ler = stackPlanformPointsXY_Wn_Ler[: tip_index + 1, :]
stackTrailingPointsXY_Wn_Ler = np.flip(
    stackPlanformPointsXY_Wn_Ler[tip_index:, :], axis=0
)

# Add zero z components so the curves are 3D points (in wing axes, relative to the
# leading edge root point), which is the form that Wing.from_edge_points expects.
leadingEdgePoints_Wn_Ler = np.column_stack(
    (stackLeadingPointsXY_Wn_Ler, np.zeros(len(stackLeadingPointsXY_Wn_Ler)))
)
trailingEdgePoints_Wn_Ler = np.column_stack(
    (stackTrailingPointsXY_Wn_Ler, np.zeros(len(stackTrailingPointsXY_Wn_Ler)))
)

# Set the starting values for the number of flap cycles to run the simulation for, the
# number of chordwise Panels, and the number of sections to map on each Wing half (there
# will be this number + 1 WingCrossSections per Wing half). These only define the
# reference problem for the convergence analysis below, so they are set to the coarse
# end of its sweep (3 spanwise sections gives an average Panel aspect ratio of about 4
# at 3 chordwise Panels). The validation results come from the solver that analysis runs
# at whatever values it finds are converged.
NUM_FLAPS = 1
NUM_CHORDWISE_PANELS = 3
NUM_SPANWISE_SECTIONS = 3

# Set the chordwise spacing scheme for the Panels. This is set to uniform, as is
# standard for UVLM simulations.
CHORDWISE_SPACING = "uniform"


def validation_flap_angle_series(
    cycle_angle_rad: float | np.ndarray,
) -> float | np.ndarray:
    """Returns the flap angle in degrees at a position within one flap cycle.

    The position is measured in cycle radians, so 0.0 is the start of a flap and 2.0 *
    pi is the end of that same flap. This function contains no flapping frequency. The
    frequency enters later, when the solver converts time in seconds to cycle radians.
    The Fourier series is fourth order, and its coefficients were calculated by Yeo et
    al., 2011.

    :param cycle_angle_rad: A float or a (N,) ndarray of floats representing the
        position or positions within a flap cycle at which to evaluate the flap angle.
        The units are radians of flap cycle.
    :return: A float or a (N,) ndarray of floats representing the flap angle or angles
        at the given cycle positions. The units are degrees.
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
    return np.rad2deg(
        a_0
        + a_1 * np.cos(1 * cycle_angle_rad)
        + b_1 * np.sin(1 * cycle_angle_rad)
        + a_2 * np.cos(2 * cycle_angle_rad)
        + b_2 * np.sin(2 * cycle_angle_rad)
        + a_3 * np.cos(3 * cycle_angle_rad)
        + b_3 * np.sin(3 * cycle_angle_rad)
        + a_4 * np.cos(4 * cycle_angle_rad)
        + b_4 * np.sin(4 * cycle_angle_rad)
    )


# The custom spacing API expects the flap angle split into three parts: the value at the
# start of a flap, an amplitude, and a unit shape that starts at 0.0, returns to 0.0
# after 2.0 * pi, and has an amplitude of 1.0. Sample the series over one flap cycle to
# find the first two parts. The units are degrees.
flap_cycle_angles_rad = np.linspace(0.0, 2.0 * np.pi, 10001, dtype=float)
flap_angles = validation_flap_angle_series(flap_cycle_angles_rad)
validation_flap_angle_at_start = float(validation_flap_angle_series(0.0))
validation_flap_angle_amplitude = float(
    (np.max(flap_angles) - np.min(flap_angles)) / 2.0
)

# Delete the extraneous pointers.
del flap_cycle_angles_rad
del flap_angles


def validation_flap_angle_shape(cycle_angle_rad: float) -> float:
    """Returns the unit shape of the flap angle series at a position within one flap
    cycle.

    This is the third part of the split described above. It is the flap angle series
    with its value at the start of a flap subtracted and its amplitude divided out, so
    it meets the requirements for a custom spacing function.

    :param cycle_angle_rad: A float representing the position within a flap cycle at
        which to evaluate the shape. The units are radians of flap cycle.
    :return: A float representing the unit shape at the given cycle position. It is
        dimensionless.
    """
    return float(
        (validation_flap_angle_series(cycle_angle_rad) - validation_flap_angle_at_start)
        / validation_flap_angle_amplitude
    )


# Create the Airplane. The main Wing is built directly from the digitized leading and
# trailing edge curves, resampled at uniformly spaced WingCrossSections and trimmed at
# the tip by the inset. The x component of the main Wing's angles describing the
# orientation of the wing axes relative to the geometry axes (after accounting for
# symmetry) using an intrinsic xy'z" sequence is the flap angle at the start of a flap.
# The WingMovements oscillate it about that value.
validation_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leadingEdgePoints_Wn_Ler,
            trailingEdgePoints_Wn_Ler=trailingEdgePoints_Wn_Ler,
            num_wing_cross_sections=NUM_SPANWISE_SECTIONS + 1,
            airfoil=ps.geometry.airfoil.Airfoil(
                name="naca0012",
            ),
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, WING_MIDLINE_OFFSET / 2, 0.0),
            angles_Gs_to_Wn_ixyz=(validation_flap_angle_at_start, 0.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=NUM_CHORDWISE_PANELS,
            chordwise_spacing=CHORDWISE_SPACING,
            tip_trim_fraction=TIP_INSET / HALF_SPAN,
        ),
    ],
    name="Validation Airplane",
)

# Initialize empty lists to hold the WingCrossSectionMovements for the main and
# reflected main Wings.
main_wing_cross_section_movements = []
reflected_main_wing_cross_section_movements = []

# Create static WingCrossSectionMovements for each WingCrossSection in the main and
# reflected main Wings.
for i in range(NUM_SPANWISE_SECTIONS + 1):
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


def time_normalized_validation_geometry_sweep_function_rad(
    time: float | np.ndarray,
) -> float | np.ndarray:
    """This function takes in the time during a flap cycle and returns the flap angle in
    radians. It uses a normalized flapping frequency of 1 Hz, and is based on a fourth-
    order Fourier series. The coefficients were calculated by Yeo et al., 2011.

    :param time: float or a (N,) ndarray of floats This is a single time or a ndarray of
        N times at which to calculate the flap angle. The units are seconds.
    :return flap_angle: float or a (N,) ndarray of floats This is a single flap angle or
        a ndarray of N flap angles at the inputted time value or values. The units are
        radians.
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
    ampAngles_Gs_to_Wn_ixyz=(validation_flap_angle_amplitude, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1 / VALIDATION_FLAPPING_FREQUENCY, 0.0, 0.0),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=(validation_flap_angle_shape, "sine", "sine"),
)
reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=validation_airplane.wings[1],
    wing_cross_section_movements=reflected_main_wing_cross_section_movements,
    ampAngles_Gs_to_Wn_ixyz=(validation_flap_angle_amplitude, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1 / VALIDATION_FLAPPING_FREQUENCY, 0.0, 0.0),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=(validation_flap_angle_shape, "sine", "sine"),
)

# Delete the extraneous pointers.
del main_wing_cross_section_movements
del reflected_main_wing_cross_section_movements

# Define the AirplaneMovement that contains the WingMovements.
validation_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=validation_airplane,
    wing_movements=[main_wing_movement, reflected_main_wing_movement],
)

# Delete the extraneous pointers.
del validation_airplane
del main_wing_movement
del reflected_main_wing_movement

# Define an OperatingPoint and OperatingPointMovement corresponding to the conditions of
# the validation study.
validation_operating_point = ps.operating_point.OperatingPoint(
    vCg__E=VALIDATION_VELOCITY, alpha=VALIDATION_ALPHA
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
    num_cycles=NUM_FLAPS,
)

# Delete the extraneous pointers.
del validation_airplane_movement
del validation_operating_point_movement

# Define the reference UnsteadyProblem for the convergence analysis. It only needs each
# iteration's final results, which is all the analysis compares.
reference_problem = ps.problems.UnsteadyProblem(
    movement=validation_movement,
    only_final_results=True,
)

# Delete the extraneous pointer.
del validation_movement

# Run the convergence analysis. This will run several simulations, modifying the wake
# state, number of flaps, average Panel aspect ratio, and number of chordwise Panels
# with each iteration, until each of the Airplane's final-flap mean load coefficients
# stops changing by more than the relative tolerance (rtol) plus the absolute tolerance
# (atol) between successive iterations. The main Wing is edge-defined, so the analysis
# refines it by resampling its edge curves into more WingCrossSections. Each iteration's
# results are cached in a JSON file next to this script, so rerunning this study with an
# existing cache skips the simulations it has already done. The analysis also rebuilds
# and runs a solver at the converged parameters, and that solver's results are the ones
# compared to the experimental results below. See the analyze_unsteady_convergence
# function docstring for more details.
(
    converged_prescribed_wake,
    converged_num_flaps,
    converged_panel_aspect_ratio,
    converged_num_chordwise_panels,
    converged_solver,
) = ps.convergence.analyze_unsteady_convergence(
    ref_problem=reference_problem,
    prescribed_wake=True,
    free_wake=True,
    num_cycles_bounds=(1, 4),
    panel_aspect_ratio_bounds=(4, 2),
    num_chordwise_panels_bounds=(6, 10),
    rtol=0.025,
    atol=0.0025,
    show_solver_progress=True,
    resolve_converged_solver=True,
    cache_path=VALIDATION_DIRECTORY / "validation_convergence_cache.json",
)

# Delete the extraneous pointer.
del reference_problem

# The analysis returns Nones if it did not find a converged case within its bounds, in
# which case there is nothing to validate.
if (
    converged_prescribed_wake is None
    or converged_num_flaps is None
    or converged_num_chordwise_panels is None
    or converged_solver is None
):
    raise RuntimeError(
        "The convergence analysis did not find a converged case within its bounds."
    )

# Print and log the converged parameters.
convergence_message = (
    "Converged parameters: prescribed wake = "
    + str(converged_prescribed_wake)
    + ", flaps = "
    + str(converged_num_flaps)
    + ", Panel aspect ratio = "
    + str(converged_panel_aspect_ratio)
    + ", chordwise Panels = "
    + str(converged_num_chordwise_panels)
)
print("\n" + convergence_message)
validation_logger.info(convergence_message)

# Delete the extraneous pointers.
del leadingEdgePoints_Wn_Ler
del trailingEdgePoints_Wn_Ler

# Define the position of the points of interest and the area of their rectangles. These
# values were extracted by digitizing the figures in Yeo et al., 2011.
blueTrailingPointsXY_Wn_Ler = [0.060, 0.036]
BLUE_TRAILING_AREA = 0.072 * 0.024
blueMiddlePointsXY_Wn_Ler = [0.036, 0.036]
BLUE_MIDDLE_AREA = 0.072 * 0.024
blueLeadingPointsXY_Wn_Ler = [0.012, 0.036]
BLUE_LEADING_AREA = 0.072 * 0.024
orangeTrailingPointsXY_Wn_Ler = [0.05532, 0.107]
ORANGE_TRAILING_AREA = 0.07 * 0.02112
orangeMiddlePointsXY_Wn_Ler = [0.0342, 0.107]
ORANGE_MIDDLE_AREA = 0.07 * 0.02112
orangeLeadingPointsXY_Wn_Ler = [0.01308, 0.107]
ORANGE_LEADING_AREA = 0.07 * 0.02112
greenTrailingPointsXY_Wn_Ler = [0.04569, 0.162825]
GREEN_TRAILING_AREA = 0.04165 * 0.015
greenMiddlePointsXY_Wn_Ler = [0.03069, 0.176]
GREEN_MIDDLE_AREA = 0.06565 * 0.015
greenLeadingPointsXY_Wn_Ler = [0.01569, 0.1775]
GREEN_LEADING_AREA = 0.071 * 0.015

# The converged solver has already been run, and its results are compared to the
# experimental results directly. It was run with only final results, so it holds loads
# only for the time steps in the final flap, from its first results step onward. That is
# the only flap the comparison uses.
first_results_step = converged_solver.first_results_step

# Create a variable to hold the time in seconds at each of the time steps with results.
times = (
    np.arange(first_results_step, converged_solver.num_steps)
    * converged_solver.delta_time
)

# Discretize the time period of the final flap analyzed into 100 steps. Store this to a
# ndarray.
final_flap_times = np.linspace(
    (converged_num_flaps - 1) / VALIDATION_FLAPPING_FREQUENCY,
    converged_num_flaps / VALIDATION_FLAPPING_FREQUENCY,
    100,
    endpoint=False,
)

# Discretize the normalized flap cycle times into 100 steps. Store this to a ndarray.
normalized_times = np.linspace(0, 1, 100, endpoint=False)

# Pull the experimental pressure vs. time histories from the digitized data. These data
# sets are stored in CSV files in the experimental data subdirectory. The pressure units
# used are inAq and time units are normalized flap cycle times from 0 to 1.
exp_blue_trailing_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "blue_trailing_point_experimental_pressures.csv",
    delimiter=",",
)
exp_blue_middle_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "blue_middle_point_experimental_pressures.csv",
    delimiter=",",
)
exp_blue_leading_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "blue_leading_point_experimental_pressures.csv",
    delimiter=",",
)
exp_orange_trailing_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "orange_trailing_point_experimental_pressures.csv",
    delimiter=",",
)
exp_orange_middle_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "orange_middle_point_experimental_pressures.csv",
    delimiter=",",
)
exp_orange_leading_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "orange_leading_point_experimental_pressures.csv",
    delimiter=",",
)
exp_green_trailing_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "green_trailing_point_experimental_pressures.csv",
    delimiter=",",
)
exp_green_middle_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "green_middle_point_experimental_pressures.csv",
    delimiter=",",
)
exp_green_leading_point_pressures = np.genfromtxt(
    EXPERIMENTAL_DATA_DIRECTORY / "green_leading_point_experimental_pressures.csv",
    delimiter=",",
)

# Interpolate the experimental pressure data to ensure that they all reference the same
# normalized timescale.
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
    248.84 * exp_blue_trailing_point_pressures_norm * BLUE_TRAILING_AREA
)
exp_blue_middle_normal_forces = (
    248.84 * exp_blue_middle_point_pressures_norm * BLUE_MIDDLE_AREA
)
exp_blue_leading_normal_forces = (
    248.84 * exp_blue_leading_point_pressures_norm * BLUE_LEADING_AREA
)
exp_orange_trailing_normal_forces = (
    248.84 * exp_orange_trailing_point_pressures_norm * ORANGE_TRAILING_AREA
)
exp_orange_middle_normal_forces = (
    248.84 * exp_orange_middle_point_pressures_norm * ORANGE_MIDDLE_AREA
)
exp_orange_leading_normal_forces = (
    248.84 * exp_orange_leading_point_pressures_norm * ORANGE_LEADING_AREA
)
exp_green_trailing_normal_forces = (
    248.84 * exp_green_trailing_point_pressures_norm * GREEN_TRAILING_AREA
)
exp_green_middle_normal_forces = (
    248.84 * exp_green_middle_point_pressures_norm * GREEN_MIDDLE_AREA
)
exp_green_leading_normal_forces = (
    248.84 * exp_green_leading_point_pressures_norm * GREEN_LEADING_AREA
)

# Convert each experimental panel's normal force time history to a time history of the
# force's z component in geometry axes.
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

# Calculate the net experimental force's z component in geometry axes. This is
# multiplied by two because the experimental panels only cover one of the symmetric wing
# halves.
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

# Initialize a ndarray to hold the net force's wind axes' z component.
stackExpNetForcesZ_W = np.zeros(stackExpNetForcesZ_G.size, dtype=float)

# Get the passive transformation matrix which maps in homogeneous coordinates from the
# first Airplane's geometry axes relative to the first Airplane's CG to wind axes
# relative to the first Airplane's CG.
T_pas_GP1_CgP1_to_W_CgP1 = validation_operating_point.T_pas_GP1_CgP1_to_W_CgP1

# Delete the extraneous pointer.
del validation_operating_point

# Transform from the first Airplane's geometry axes to wind axes.
for force_id, expNetForceZ_GP1 in enumerate(stackExpNetForcesZ_G):
    expNetForceHomog_GP1 = np.array([0.0, 0.0, expNetForceZ_GP1, 0.0], dtype=float)
    expNetForceHomog_W = T_pas_GP1_CgP1_to_W_CgP1 @ expNetForceHomog_GP1
    expNetForceZ_W = expNetForceHomog_W[2]
    stackExpNetForcesZ_W[force_id] = expNetForceZ_W

# Get the experimental lift values. Lift is defined as the force's wind axes' z
# component multiplied by negative one.
exp_lifts = -1 * stackExpNetForcesZ_W

# Get the converged solver's SteadyProblems' Airplanes for the time steps with results.
airplanes = []
for steady_problem in converged_solver.steady_problems[first_results_step:]:
    airplanes.append(steady_problem.airplanes[0])

# Initialize a ndarray to hold the force at each time step with results (in wind axes).
stackSimForces_W = np.zeros((3, len(airplanes)), dtype=float)

# Iterate through the time steps and populate the ndarray.
for step, airplane in enumerate(airplanes):
    stackSimForces_W[:, step] = airplane.forces_W

# Initialize the figure and axes of the experimental versus simulated lift plot.
lift_figure, lift_axes = plt.subplots(figsize=(5, 4))

# Get the simulated lift values. Lift is defined as the force's wind axes' z component
# multiplied by negative one.
sim_lifts = -1 * stackSimForces_W[2, :]

# Interpolate the simulated lift values to find them with respect to the normalized
# final flap timescale.
final_flap_sim_lifts = np.interp(final_flap_times, times, sim_lifts[:])

SIM_LIFT_COLOR = "#D81E5B"
EXP_LIFT_COLOR = "#003F91"

NUM_MARKERS = 6
MARKER_SIZE = 8
TEXT_COLOR = "black"
FIGURE_BACKGROUND_COLOR = "None"

lift_axes.spines.right.set_visible(False)
lift_axes.spines.top.set_visible(False)
lift_axes.spines.bottom.set_color(TEXT_COLOR)
lift_axes.spines.left.set_color(TEXT_COLOR)
lift_axes.xaxis.label.set_color(TEXT_COLOR)
lift_axes.yaxis.label.set_color(TEXT_COLOR)
lift_axes.tick_params(axis="x", colors=TEXT_COLOR)
lift_axes.tick_params(axis="y", colors=TEXT_COLOR)
lift_figure.patch.set_facecolor(FIGURE_BACKGROUND_COLOR)
lift_axes.set_facecolor(FIGURE_BACKGROUND_COLOR)

marker_spacing = 1.0 / NUM_MARKERS

# Plot the simulated lift values. The x axis is set to the normalized times, which may
# seem odd because we just interpolated to get them in terms of the normalized final
# flap times. But, they are discretized in exactly the same way as the normalized times,
# just horizontally shifted.
lift_axes.plot(
    normalized_times,
    final_flap_sim_lifts,
    label="Simulated",
    color=SIM_LIFT_COLOR,
    marker=".",
    markevery=(marker_spacing * 0 / 2, marker_spacing),
    markersize=MARKER_SIZE,
)

# Plot the experimental lift values.
lift_axes.plot(
    normalized_times,
    exp_lifts,
    label="Experimental",
    color=EXP_LIFT_COLOR,
    marker=".",
    markevery=(marker_spacing * 1 / 2, marker_spacing),
    markersize=MARKER_SIZE,
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
    facecolor=FIGURE_BACKGROUND_COLOR,
    edgecolor=FIGURE_BACKGROUND_COLOR,
    labelcolor=TEXT_COLOR,
)

# Save the lift comparison figure.
lift_figure.savefig(
    fname=VALIDATION_DIRECTORY / "lift_validation.png",
    dpi=300,
    bbox_inches="tight",
)

# Delete the extraneous pointers.
del airplanes
del stackSimForces_W
del step

# Calculate the lift mean absolute error (MAE). The experimental and simulated lift
# comparison here is valid because, due to the interpolation steps, the experimental and
# simulated lifts time histories are discretized so that they are with respect to the
# same timescale.
lift_absolute_errors = np.abs(final_flap_sim_lifts - exp_lifts)
lift_mean_absolute_error = np.mean(lift_absolute_errors)

sim_lift_rms = math.sqrt(np.mean(final_flap_sim_lifts**2))
exp_lift_rms = math.sqrt(np.mean(exp_lifts**2))
lift_rmsape = 100 * abs((sim_lift_rms - exp_lift_rms) / exp_lift_rms)

# Print and log the RMS lift results.
lift_rmsape_message = (
    "Lift RMS Absolute Percent Error: " + str(np.round(lift_rmsape, 2)) + "%"
)
sim_lift_rms_message = "Simulated Lift RMS: " + str(np.round(sim_lift_rms, 4)) + " N"
exp_lift_rms_message = "Experimental Lift RMS: " + str(np.round(exp_lift_rms, 4)) + " N"
print("\n" + lift_rmsape_message)
print(sim_lift_rms_message)
print(exp_lift_rms_message)
validation_logger.info(lift_rmsape_message)
validation_logger.info(sim_lift_rms_message)
validation_logger.info(exp_lift_rms_message)

# Print and log the MAE.
lift_mean_absolute_error_message = (
    "Mean Absolute Error on Lift: " + str(np.round(lift_mean_absolute_error, 4)) + "N"
)
print("\n" + lift_mean_absolute_error_message)
validation_logger.info(lift_mean_absolute_error_message)

ps.output.draw(
    solver=converged_solver,
    show_wake_vortices=True,
    scalar_type="lift",
    save=True,
    path=VALIDATION_DIRECTORY / "draw.webp",
)
