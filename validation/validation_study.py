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

# This wing planform has a rounded tip so the outermost wing cross section needs to
# be inset some amount. This value is in meters.
tip_inset = 0.005

# Import the extracted coordinates from the paper’s diagram of the planform. The
# resulting array is of the form [spanwise coordinate, chordwise coordinate],
# and is ordered from the leading edge root, to the tip, to the trailing edge root.
# The origin is the trailing edge root point. The positive spanwise axis extends from
# root to tip and the positive chordwise axis from trailing edge to leading edge. The
# coordinates are in millimeters.
planform_coords = np.genfromtxt("extracted_planform_coordinates.csv", delimiter=",")

# Convert the coordinates to meters.
planform_coords = planform_coords / 1000

# Set the origin to the leading edge root point.
planform_coords = planform_coords - np.array([0, chord])

# Switch the sign of the chordwise coordinates.
planform_coords = planform_coords * np.array([1, -1])

# Swap the axes to the form [chordwise coordinate, spanwise coordinate]. The
# coordinates are now in the geometry frame projected on the XY plane.
planform_coords[:, [0, 1]] = planform_coords[:, [1, 0]]

# Find the index of the point where the planform x-coordinate equals the half span.
tip_index = np.where(planform_coords[:, 1] == half_span)[0][0]

# Using the tip index, split the coordinates into two arrays of leading and trailing
# edge coordinates.
leading_coords = planform_coords[:tip_index, :]
trailing_coords = np.flip(planform_coords[tip_index:, :], axis=0)

# Set the number of flap cycles to run the simulation for. The converged result is 3
# flaps.
num_flaps = 3

# Set the number of chordwise panels. The converged result is 5 panels.
num_chordwise_panels = 5

# Set the number of sections to map on each wing half. There will be this number +1
# wing cross sections per wing half. The converged result is 18 spanwise sections.
num_spanwise_sections = 18

# Set the chordwise spacing scheme for the panels. This is set to uniform,
# as is standard for UVLM simulations.
chordwise_spacing = "uniform"

# Calculate the spanwise difference between the wing cross sections.
spanwise_step = (half_span - tip_inset) / num_spanwise_sections

# Define four arrays to hold the coordinates of the front and back points of each
# section’s left and right wing cross sections.
front_left_vertices = np.zeros((num_spanwise_sections, 2))
front_right_vertices = np.zeros((num_spanwise_sections, 2))
back_left_vertices = np.zeros((num_spanwise_sections, 2))
back_right_vertices = np.zeros((num_spanwise_sections, 2))

# Iterate through the locations of the future sections to populate the wing cross
# section coordinates.
for spanwise_loc in range(num_spanwise_sections):
    # Find the y-coordinates of the vertices.
    front_left_vertices[spanwise_loc, 1] = spanwise_loc * spanwise_step
    back_left_vertices[spanwise_loc, 1] = spanwise_loc * spanwise_step
    front_right_vertices[spanwise_loc, 1] = (spanwise_loc + 1) * spanwise_step
    back_right_vertices[spanwise_loc, 1] = (spanwise_loc + 1) * spanwise_step

    # Interpolate between the leading edge coordinates to find the x-coordinate of
    # the front left vertex.
    front_left_vertices[spanwise_loc, 0] = np.interp(
        spanwise_loc * spanwise_step,
        leading_coords[:, 1],
        leading_coords[:, 0],
    )

    # Interpolate between the trailing edge coordinates to find the x-coordinate of
    # the back left vertex.
    back_left_vertices[spanwise_loc, 0] = np.interp(
        spanwise_loc * spanwise_step,
        trailing_coords[:, 1],
        trailing_coords[:, 0],
    )

    # Interpolate between the leading edge coordinates to find the x-coordinate of
    # the front right vertex.
    front_right_vertices[spanwise_loc, 0] = np.interp(
        (spanwise_loc + 1) * spanwise_step,
        leading_coords[:, 1],
        leading_coords[:, 0],
    )

    # Interpolate between the trailing edge coordinates to find the x-coordinate of
    # the back right vertex.
    back_right_vertices[spanwise_loc, 0] = np.interp(
        (spanwise_loc + 1) * spanwise_step,
        trailing_coords[:, 1],
        trailing_coords[:, 0],
    )

# Define an empty list to hold the wing cross sections.
validation_airplane_wing_cross_sections = []

# Iterate through the wing cross section vertex arrays to create the wing cross
# section objects.
for i in range(num_spanwise_sections):

    # Get the left wing cross section’s vertices at this position.
    this_front_left_vertex = front_left_vertices[i, :]
    this_back_left_vertex = back_left_vertices[i, :]

    # Get this wing cross section’s leading and trailing edge x-coordinates.
    this_x_le = this_front_left_vertex[0]
    this_x_te = this_back_left_vertex[0]

    # Get this wing cross section’s leading edge y-coordinate.
    this_y_le = this_front_left_vertex[1]

    # Calculate this wing cross section’s chord.
    this_chord = this_x_te - this_x_le

    # Define this wing cross section object.
    this_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        x_le=this_x_le,
        y_le=this_y_le,
        chord=this_chord,
        airfoil=ps.geometry.airfoil.Airfoil(
            name="naca0000",
        ),
        num_spanwise_panels=1,
    )

    # Append this wing cross section to the list of wing cross sections.
    validation_airplane_wing_cross_sections.append(this_wing_cross_section)

    # Check if this the last section.
    if i == num_spanwise_sections - 1:
        # If so, get the right wing cross section vertices at this position.
        this_front_right_vertex = front_right_vertices[i, :]
        this_back_right_vertex = back_right_vertices[i, :]

        # Get this wing cross section’s leading and trailing edge x-coordinates.
        this_x_le = this_front_right_vertex[0]
        this_x_te = this_back_right_vertex[0]

        # Get this wing cross section’s leading edge y-coordinate.
        this_y_le = this_front_right_vertex[1]

        # Calculate this wing cross section’s chord.
        this_chord = this_x_te - this_x_le

        # Define this wing cross section object.
        this_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            x_le=this_x_le,
            y_le=this_y_le,
            chord=this_chord,
            airfoil=ps.geometry.airfoil.Airfoil(
                name="naca0000",
            ),
            num_spanwise_panels=1,
        )

        # Append this wing cross section to the list of wing cross sections.
        validation_airplane_wing_cross_sections.append(this_wing_cross_section)

# Define the validation airplane object.
validation_airplane = ps.geometry.airplane.Airplane(
    name="Validation Airplane",
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=validation_airplane_wing_cross_sections,
            symmetric=True,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
        ),
    ],
)

# Delete the extraneous pointer.
del validation_airplane_wing_cross_sections

# Initialize an empty list to hold each wing cross section movement object.
validation_wing_cross_section_movements = []

# Define the first wing cross section movement, which is stationary.
first_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=validation_airplane.wings[0].wing_cross_sections[0],
    )
)

# Append the first wing cross section movement object to the list.
validation_wing_cross_section_movements.append(first_wing_cross_section_movement)

# Delete the extraneous pointer.
del first_wing_cross_section_movement


def validation_geometry_sweep_function(time):
    """This function takes in the time during a flap cycle and returns the flap angle
    in degrees. It uses the flapping frequency defined in the encompassing script,
    and is based on a fourth-order Fourier series. The coefficients were calculated
    by the authors of Yeo et al., 2011.

    :param time: float or 1D array of floats
        This is a single time or an array of time values at which to calculate the
        flap angle. The units are seconds.
    :return flap_angle: float or 1D array of floats
        This is a single flap angle or an array of flap angle values at the inputted
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
    f = validation_flapping_frequency

    # Calculate and return the flap angle(s).
    flap_angle = (
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
    return flap_angle


def normalized_validation_geometry_sweep_function_rad(time):
    """This function takes in the time during a flap cycle and returns the flap angle
    in radians. It uses a normalized flapping frequency of 1 Hertz, and is based on a
    fourth-order Fourier series. The coefficients were calculated by the authors of
    Yeo et al., 2011.

    :param time: float or 1D array of floats
        This is a single time or an array of time values at which to calculate the
        flap angle. The units are seconds.
    :return flap_angle: float or 1D array of floats
        This is a single flap angle or an array of flap angle values at the inputted
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
    flap_angle = -(
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
    return flap_angle


# Iterate through each of the wing cross sections.
for j in range(1, num_spanwise_sections + 1):
    # Define the wing cross section movement for this wing cross section. The
    # amplitude and period are both set to one because the true amplitude and period
    # are already accounted for in the custom sweep function.
    this_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=validation_airplane.wings[0].wing_cross_sections[j],
            sweeping_amplitude=1,
            sweeping_period=1,
            sweeping_spacing="custom",
            custom_sweep_function=validation_geometry_sweep_function,
        )
    )

    # Append this wing cross section movement to the list of wing cross section
    # movements.
    validation_wing_cross_section_movements.append(this_wing_cross_section_movement)

# Define the wing movement object that contains the wing cross section movements.
validation_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=validation_airplane.wings[0],
    wing_cross_section_movements=validation_wing_cross_section_movements,
)

# Delete the extraneous pointer.
del validation_wing_cross_section_movements

# Define the airplane movement that contains the wing movement.
validation_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=validation_airplane,
    wing_movements=[
        validation_main_wing_movement,
    ],
)

# Delete the extraneous pointers.
del validation_airplane
del validation_main_wing_movement

# Define an operating point corresponding to the conditions of the validation study.
validation_operating_point = ps.operating_point.OperatingPoint(
    vCg__E=validation_velocity, alpha=validation_alpha
)

# Define an operating point movement that contains the operating point.
validation_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=validation_operating_point
    )
)

# Get the geometry-to-wind transformation matrix, which will later be used to convert
# the experimental pressure-based-lift measurements to wind axes.
T_pas_G_Cg_to_W_Cg = validation_operating_point.T_pas_G_Cg_to_W_Cg

# Delete the extraneous pointer.
del validation_operating_point

# Calculate the period of this case’s flapping motion. The units are in seconds.
validation_flapping_period = 1 / validation_flapping_frequency

# Calculate the time step (in seconds) so that the area of the wake ring vortices
# roughly equal the area of the bound ring vortices.
validation_delta_time = (
    validation_airplane_movement.base_airplane.c_ref
    / num_chordwise_panels
    / validation_velocity
)

# Calculate the number of steps required for the wing to have flapped the prescribed
# number of times.
validation_num_steps = math.ceil(
    num_flaps / validation_flapping_frequency / validation_delta_time
)

# Define the overall movement.
validation_movement = ps.movement.Movement(
    airplane_movements=[validation_airplane_movement],
    operating_point_movement=validation_operating_point_movement,
    num_steps=validation_num_steps,
    delta_time=validation_delta_time,
)

# Delete the extraneous pointers.
del validation_airplane_movement
del validation_operating_point_movement

# Define the validation problem.
validation_problem = ps.problems.UnsteadyProblem(
    movement=validation_movement,
)

# Delete the extraneous pointer.
del validation_movement

# Define the validation solver.
validation_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=validation_problem,
    )
)

# Delete the extraneous pointer.
del validation_problem

# Define the position of the coordinates of interest and the area of their
# rectangles. These values were extracted by digitizing the figures in Yeo et al.,
# 2011.
blue_trailing_point_coords = [0.060, 0.036]
blue_trailing_area = 0.072 * 0.024
blue_middle_point_coords = [0.036, 0.036]
blue_middle_area = 0.072 * 0.024
blue_leading_point_coords = [0.012, 0.036]
blue_leading_area = 0.072 * 0.024
orange_trailing_point_coords = [0.05532, 0.107]
orange_trailing_area = 0.07 * 0.02112
orange_middle_point_coords = [0.0342, 0.107]
orange_middle_area = 0.07 * 0.02112
orange_leading_point_coords = [0.01308, 0.107]
orange_leading_area = 0.07 * 0.02112
green_trailing_point_coords = [0.04569, 0.162825]
green_trailing_area = 0.04165 * 0.015
green_middle_point_coords = [0.03069, 0.176]
green_middle_area = 0.06565 * 0.015
green_leading_point_coords = [0.01569, 0.1775]
green_leading_area = 0.071 * 0.015

# Run the validation solver. This validation study was run using a prescribed wake.
validation_solver.run(prescribed_wake=True)

# Create a variable to hold the time in seconds at each of the simulation’s time steps.
times = np.linspace(
    0,
    validation_num_steps * validation_delta_time,
    validation_num_steps,
    endpoint=False,
)

# Discretize the time period of the final flap analyzed into 100 steps. Store this to
# an array.
final_flap_times = np.linspace(
    validation_flapping_period * (num_flaps - 1),
    validation_flapping_period * num_flaps,
    100,
    endpoint=False,
)

# Discretize the normalized flap cycle times into 100 steps. Store this to an array.
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

# Convert each experimental panel’s normal force time history to a time history of
# the force in the z axis of the geometry frame. This is done by finding the vertical
# component given the wing’s sweep angle at each time step.
exp_blue_trailing_z_forces_geometry_axes = exp_blue_trailing_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_blue_middle_z_forces_geometry_axes = exp_blue_middle_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_blue_leading_z_forces_geometry_axes = exp_blue_leading_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_orange_trailing_z_forces_geometry_axes = exp_orange_trailing_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_orange_middle_z_forces_geometry_axes = exp_orange_middle_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_orange_leading_z_forces_geometry_axes = exp_orange_leading_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_green_trailing_z_forces_geometry_axes = exp_green_trailing_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_green_middle_z_forces_geometry_axes = exp_green_middle_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)
exp_green_leading_z_forces_geometry_axes = exp_green_leading_normal_forces * np.cos(
    normalized_validation_geometry_sweep_function_rad(normalized_times)
)

# Calculate the net experimental force in the z direction of the geometry axes. This
# is the sum of all the z-direction-forces on each of the experimental panels
# multiplied by two ( because the experimental panels only cover one of the symmetric
# wing halves).
exp_net_z_forces_geometry_axes = 2 * (
    exp_blue_trailing_z_forces_geometry_axes
    + exp_blue_middle_z_forces_geometry_axes
    + exp_blue_leading_z_forces_geometry_axes
    + exp_orange_trailing_z_forces_geometry_axes
    + exp_orange_middle_z_forces_geometry_axes
    + exp_orange_leading_z_forces_geometry_axes
    + exp_green_trailing_z_forces_geometry_axes
    + exp_green_middle_z_forces_geometry_axes
    + exp_green_leading_z_forces_geometry_axes
)

# Initialize an array to hold the wind frame z-direction forces.
exp_net_z_forces_wind_axes = np.zeros(exp_net_z_forces_geometry_axes.size)

# To eventually compare the experimental and simulated lifts, use the rotation matrix
# to shift the experimental force z-direction force to the wind frame.
for force_id, force in enumerate(exp_net_z_forces_geometry_axes):
    exp_net_force_geometry_axes = np.array([0, 0, force])
    exp_net_force_wind_axes = ps._transformations.apply_T_to_vectors(
        T_pas_G_Cg_to_W_Cg, exp_net_force_geometry_axes, has_point=False
    )
    exp_net_z_force_wind_axes = exp_net_force_wind_axes[2]
    exp_net_z_forces_wind_axes[force_id] = exp_net_z_force_wind_axes

# Get the experimental lift forces. Lift is defined as the z-direction wind-frame
# force multiplied by negative one.
exp_lifts = -1 * exp_net_z_forces_wind_axes

# Get this solver’s problem’s airplanes.
airplanes = []
for steady_problem in validation_solver.steady_problems:
    airplanes.append(steady_problem.airplanes[0])

# Initialize matrices to hold the forces at each time step.
sim_forces_wind_axes = np.zeros((3, validation_num_steps))

# Iterate through the time steps and add the results to their respective matrices.
for step in range(validation_num_steps):
    # Get the airplane at this time step.
    airplane = airplanes[step]
    # Add the total forces on the airplane at this time step to the list
    # of simulated forces.
    sim_forces_wind_axes[:, step] = airplane.total_near_field_force_wind_axes

# Initialize the figure and axes of the experimental versus simulated lift plot.
lift_figure, lift_axes = plt.subplots(figsize=(5, 4))

# Get the simulated lift forces. Lift is defined as the z-direction wind-frame force
# multiplied by negative one.
sim_lifts = -1 * sim_forces_wind_axes[2, :]

# Interpolate the simulated lift forces to find them with respect to the normalized
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

# Plot the simulated lift forces. The x-axis is set to the normalized times,
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

# Plot the experimental lift forces.
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
del sim_forces_wind_axes
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
