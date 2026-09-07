"""Demonstrates running Ptera Software's FreeFlightUnsteadyRingVortexLatticeMethodSolver
on a static-geometry configuration.

The simulation releases a statically stable glider into an unpowered glide, coupling the
unsteady aerodynamics to MuJoCo's rigid body dynamics so the airplane flies a free six-
degree-of-freedom trajectory through the scene.

The script will likely take several minutes to run, and will log simulation progress and
results in a log file.
"""

import logging

# First, import the software's main package. Note that if you wished to import this
# software into another package, you would first install it by running "pip install
# pterasoftware" in your terminal.
import pterasoftware as ps

# Configure logging to write info level messages to a file. To display log messages on
# the console alongside progress bars instead, omit the handler argument.
ps.set_up_logging(level="Info", handler=logging.FileHandler("example_solver.log"))

# Create an Airplane with our custom geometry. I am going to declare every parameter for
# Airplane, even though most of them have usable default values. This is for educational
# purposes, but keep in mind that it makes the code much longer than it needs to be. For
# details about each parameter, read the detailed class docstring. The same caveats
# apply to the other classes, methods, and functions I call in this script.
#
# This is a simple glider: a cambered main wing, a symmetric horizontal stabilizer set
# at a negative incidence relative to the main wing, and a single vertical stabilizer.
# The planform, center of gravity, and inertia (defined further below) were tuned for
# static pitch and yaw stability and verified in XFLR5. The negative horizontal
# stabilizer incidence provides the restoring pitch moment that keeps the glider on a
# bounded, damped trajectory rather than diverging once it is released into free flight.
example_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=10,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                        outline_A_Lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                        outline_A_Lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.0, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=6,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_Lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_Lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="Horizontal Stabilizer",
            Ler_Gs_Cgs=(5.0, 0.0, 0.5),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=6,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type=None,
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_Lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 2.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type=None,
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_Lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="Vertical Stabilizer",
            Ler_Gs_Cgs=(5.0, 0.0, 1.0),
            angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=4,
            chordwise_spacing="uniform",
        ),
    ],
    name="Example Airplane",
    Cg_GP1_CgP1=(0.0, 0.0, 0.0),
    weight=420.0,
    s_ref=None,
    c_ref=None,
    b_ref=None,
)

# Now define each Wing's WingMovements (and the WingCrossSectionMovements they contain).
# The main Wing was defined with symmetric=True, mirror_only=False, and a symmetry plane
# coincident with its axes' xz plane, giving it type 4 symmetry (see the Wing class
# documentation for more details on symmetry types), so it stays a single Wing and is
# not split into a separate reflected Wing. The Airplane therefore has exactly three
# Wings, and we need one WingMovement per Wing.
#
# This is a free flight of a rigid glider, so there is no prescribed flapping or
# deformation. Every prescribed-motion amplitude is left at its zero default. The only
# motion in the simulation comes from the rigid body dynamics during the free flight
# phase. We build the movements with a loop, since with all-default (static) motion they
# are identical in form across the three Wings.
wing_movements = []
for wing in example_airplane.wings:
    wing_cross_section_movements = [
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section,
        )
        for wing_cross_section in wing.wing_cross_sections
    ]
    wing_movements.append(
        ps.movements.wing_movement.WingMovement(
            base_wing=wing,
            wing_cross_section_movements=wing_cross_section_movements,
        )
    )

# Now define the example Airplane's AirplaneMovement. As with the Wings, the
# airplane-level prescribed motion is left at its zero defaults.
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=wing_movements,
)

# Delete the extraneous pointer to the WingMovements, as these are now contained within
# the AirplaneMovement. This is optional, but it can make debugging easier.
del wing_movements

# Define a new OperatingPoint describing the trimmed gliding condition. The speed and
# angle of attack are the trimmed glide found for this airframe. The initial body
# orientation (angles_E_to_BP1_izyx) pitches the airplane nose up by the angle of attack
# with zero sideslip, which places the trim velocity along the horizontal Earth x axis
# at the start of free flight. The glider flies an unpowered glide: externalFX_W must be
# zero in free flight (the dynamics never apply it), and thrust would instead be modeled
# with FreeFlightUnsteadyProblem's external_loads_fn. Standard gravity is set explicitly
# via g_E (the default is no gravitational field), while the zero initial body rates
# (omegas_BP1__E) are left at their default.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225,
    vCg__E=12.9,
    alpha=3.3,
    beta=0.0,
    angles_E_to_BP1_izyx=(0.0, 3.3, 0.0),
    externalFX_W=0.0,
    nu=15.06e-6,
    g_E=(0.0, 0.0, 9.80665),
)

# Define the OperatingPoint's FreeFlightOperatingPointMovement. It holds only the
# initial OperatingPoint. The solver populates its operating_points list with the body
# state from the dynamics integration at each time step.
operating_point_movement = (
    ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
        base_operating_point=example_operating_point,
    )
)

# Delete the extraneous pointer.
del example_operating_point

# Define the FreeFlightMovement. This contains the AirplaneMovement and the
# FreeFlightOperatingPointMovement. The glider first holds its trimmed condition for
# prescribed_num_steps time steps so the wake can develop, then the solver releases the
# rigid body dynamics for the remaining free_num_steps time steps.
movement = ps.movements.free_flight_movement.FreeFlightMovement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    delta_time=0.01292,
    prescribed_num_steps=5,
    free_num_steps=195,
)

# Delete the extraneous pointers.
del airplane_movement
del operating_point_movement

# Define the FreeFlightUnsteadyProblem. The inertia matrix is the one tuned alongside
# the planform geometry for static stability and verified in XFLR5. It is expressed in
# the first Airplane's body axes relative to the first Airplane's center of gravity,
# which is at the geometry origin. The off-diagonal terms are the body-axes products of
# inertia. The mass equals the Airplane's weight (420.0 N) divided by the gravitational
# acceleration magnitude (9.80665 m/s^2), which keeps the weight, mass, and gravity
# consistent. The solver applies the gravitational force as mass * g_E. No external
# loads are applied (external_loads_fn=None), so the glider flies an unpowered glide
# driven only by its aerodynamics, gravity, and inertia.
example_problem = ps.problems.FreeFlightUnsteadyProblem(
    movement=movement,
    mass=420.0 / 9.80665,
    I_BP1_CgP1=(
        (155.614, 0.0, -45.658),
        (0.0, 398.513, 0.0),
        (-45.658, 0.0, 508.699),
    ),
    external_loads_fn=None,
)

# Delete the extraneous pointer.
del movement

# Define a new solver. We'll create a FreeFlightUnsteadyRingVortexLatticeMethodSolver,
# which requires a FreeFlightUnsteadyProblem.
example_solver = ps.free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver(
    free_flight_unsteady_problem=example_problem,
)

# Delete the extraneous pointer.
del example_problem

# Run the solver.
example_solver.run(
    prescribed_wake=False,
    show_progress=True,
)

# Save the solved solver to a .psz file. This allows us to load the results later
# without re-running the simulation.
ps.save("example_solver.psz", example_solver)

# Load the saved solver. The loaded object is identical to the original and can be
# passed to any output function.
loaded_solver = ps.load("example_solver.psz")

# The load function is annotated as returning object because a saved file can hold any
# Ptera Software object. This assert narrows the type for type checkers and guards
# against loading the wrong file.
assert isinstance(
    loaded_solver,
    ps.free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
)

# Log the loaded solver's loads. For a free flight solver, this also logs the first
# Airplane's initial and final six-degree-of-freedom state: its position, velocity,
# orientation, angular velocity, and aerodynamic angles.
ps.output.log_results(solver=loaded_solver)

# Call the draw function on the loaded solver. For a free flight solver, the geometry is
# drawn in Earth axes at the final time step's true flight pose, so the airplane appears
# where and how it ended up after final time step. Press any key to close the plotter
# after it draws the output.
ps.output.draw(
    solver=loaded_solver,
    scalar_type="lift",
    show_streamlines=True,
    show_wake_vortices=False,
    save=True,
)

# Call the animate function on the loaded solver. For a free flight solver, this flies
# the airplane through the scene along its computed trajectory, producing an animated
# WebP saved in the same directory as this script. Press any key, after orienting the
# view, to begin the animation.
ps.output.animate(
    unsteady_solver=loaded_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
)

# Call the plotting function on the solver. This produces graphs of the loads with
# respect to time. For a free flight solver, it also plots the first Airplane's
# six-degree-of-freedom state history.
#
# Setting save_csv exports the same series the figures plot as CSV files, which is what
# you want in order to compare quantities that no single figure holds together, or to
# take the results into another tool. It writes two files here, because the loads and
# the state history do not share a time base: the loads begin at the solver's first
# results step, while the state history begins at the first time step. Note that
# save_csv is independent of save, so the data can be exported without rendering any
# images.
ps.output.plot_results_versus_time(
    unsteady_solver=loaded_solver,
    show=True,
    save=True,
    save_csv=True,
)
