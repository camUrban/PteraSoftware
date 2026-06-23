"""This script is an example of how to run Ptera Software's
FreeFlightUnsteadyRingVortexLatticeMethodSolver with a flapping wing. The main wing
flaps symmetrically while the unsteady aerodynamics are coupled to MuJoCo's rigid body
dynamics, so the airplane flies a free six-degree-of-freedom trajectory through the
scene under the loads produced by its own flapping motion."""

# First, import the software's main package. Note that if you wished to import this
# software into another package, you would first install it by running "pip install
# pterasoftware" in your terminal.
import pterasoftware as ps

# Configure logging to display info level messages. This is important for seeing the
# output from the log_results function.
ps.set_up_logging(level="Info")

# Create an Airplane with our custom geometry. I am going to declare every parameter
# for Airplane, even though most of them have usable default values. This is for
# educational purposes, but keep in mind that it makes the code much longer than it
# needs to be. For details about each parameter, read the detailed class docstring.
# The same caveats apply to the other classes, methods, and functions I call in this
# script.
#
# This is a flapping-wing demonstration airframe: a cambered main wing that flaps
# symmetrically (defined further below, a 15 degree amplitude at 1 Hz) plus a tail.
# Unlike the glider example, this airframe is not tuned for static stability or trim. It
# exists to exercise the strongly coupled free-flight solver under the large, oscillatory
# loads that flapping produces, so the trajectory is a driven flapping flight rather than
# a settled equilibrium.
example_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=8,
                    chord=1.75,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="cosine",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.5,
                    Lp_Wcsp_Lpp=(0.75, 6.0, 1.0),
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.5, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
        ),
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
                    Lp_Wcsp_Lpp=(0.5, 2.0, 1.0),
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
            num_chordwise_panels=6,
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

# Now define the main wing's root and tip WingCrossSections' WingCrossSectionMovements.
main_wing_root_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
    phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
)
main_wing_tip_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
    phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
)

# The main Wing was defined to have symmetric=True, mirror_only=False, and with a
# symmetry plane offset non coincident with the Wing's axes yz plane. Therefore,
# that Wing had type 5 symmetry (see the Wing class documentation for more details on
# symmetry types). Therefore, it was actually split into two Wings, the with the
# second Wing being a reflected version of the first. Therefore, we need to define a
# WingMovement for this reflected Wing. To start, we'll first define the reflected
# main wing's root and tip WingCrossSections' WingCrossSectionMovements.
reflected_main_wing_root_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[0],
    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
    phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
)
reflected_main_wing_tip_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[1],
    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
    ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
    phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
)

# Now define the v tail's root and tip WingCrossSections' WingCrossSectionMovements.
v_tail_root_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
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
v_tail_tip_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
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

# Now define the main wing's WingMovement, the reflected main wing's WingMovement and
# the v tail's WingMovement.
main_wing_movement = ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_section_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)
reflected_main_wing_movement = (
    ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
        base_wing=example_airplane.wings[1],
        wing_cross_section_movements=[
            reflected_main_wing_root_wing_cross_section_movement,
            reflected_main_wing_tip_wing_cross_section_movement,
        ],
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )
)
v_tail_movement = ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
    base_wing=example_airplane.wings[2],
    wing_cross_section_movements=[
        v_tail_root_wing_cross_section_movement,
        v_tail_tip_wing_cross_section_movement,
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

# Delete the extraneous pointers to the WingCrossSectionMovements, as these are now
# contained within the WingMovements. This is optional, but it can make debugging
# easier.
del main_wing_root_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement
del reflected_main_wing_root_wing_cross_section_movement
del reflected_main_wing_tip_wing_cross_section_movement
del v_tail_root_wing_cross_section_movement
del v_tail_tip_wing_cross_section_movement

# Now define the example airplane's AirplaneMovement.
airplane_movement = (
    ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement(
        base_airplane=example_airplane,
        wing_movements=[
            main_wing_movement,
            reflected_main_wing_movement,
            v_tail_movement,
        ],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )
)

# Delete the extraneous pointers to the WingMovements.
del main_wing_movement
del reflected_main_wing_movement
del v_tail_movement

# Define a new OperatingPoint describing the initial flight condition. The initial body
# orientation (angles_E_to_BP1_izyx) pitches the airplane nose up by the angle of attack
# with zero sideslip, which places the initial velocity along the horizontal Earth x axis
# at the start of free flight. No external thrust is applied (externalFX_W=0.0); in free
# flight the dynamics never apply externalFX_W, and propulsion here comes from the
# flapping motion rather than an external force. Standard gravity is set explicitly via
# g_E (the default is no gravitational field), while the zero initial body rates
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

# Define the operating point's FreeFlightOperatingPointMovement. It holds only the
# initial OperatingPoint; the solver populates its operating_points list with the body
# state from the dynamics integration at each time step.
operating_point_movement = (
    ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
        base_operating_point=example_operating_point,
    )
)

# Delete the extraneous pointer.
del example_operating_point

# Define the FreeFlightMovement. This contains the FreeFlightAirplaneMovement and the
# FreeFlightOperatingPointMovement. The airplane first holds its initial flight condition
# for prescribed_num_steps time steps so the wake can develop, then the solver releases
# the rigid body dynamics for the remaining free_num_steps time steps.
movement = ps.movements.free_flight_movement.FreeFlightMovement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    delta_time=0.01292,
    prescribed_num_steps=10,
    free_num_steps=290,
)

# Delete the extraneous pointers.
del airplane_movement
del operating_point_movement

# Define the FreeFlightUnsteadyProblem. The inertia matrix holds representative mass
# properties reused from the glider example; this airframe is a coupling demonstration and
# is not separately tuned or verified for static stability. It is expressed in the first
# Airplane's body axes relative to the first Airplane's center of gravity, which is at the
# geometry origin. The off-diagonal terms are the body-axes products of inertia. The mass
# equals the Airplane's weight (420.0 N) divided by the gravitational acceleration
# magnitude (9.80665 m/s^2), which keeps the weight, mass, and gravity consistent; the
# solver applies the gravitational force as mass * g_E. No external loads are applied
# (external_loads_fn=None), so the airplane is driven only by its flapping aerodynamics,
# gravity, and inertia.
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

# Save the solved solver to a compressed JSON file. This allows us to load the results
# later without re-running the simulation. Use ".json.gz" for gzip compression, which is
# recommended over plain JSONs for all but the smallest, unmeshed geometry objects.
ps.save("example_free_flight_solver_flapping.json.gz", example_solver)

# Load the saved solver. The loaded object is identical to the original and can be
# passed to any output function.
loaded_solver = ps.load("example_free_flight_solver_flapping.json.gz")

# Log the loaded solver's loads. For a free flight solver, this also logs the first
# Airplane's initial and final six-degree-of-freedom state: its position, velocity,
# orientation, angular velocity, and aerodynamic angles.
ps.output.log_results(solver=loaded_solver)

# Call the draw function on the loaded solver. For a free flight solver, the geometry is
# drawn in Earth axes at the final time step's true flight pose, so the airplane appears
# where and how it ended up after the final time step. Press any key to close the plotter
# after it draws the output.
ps.output.draw(
    solver=loaded_solver,
    scalar_type="lift",
    show_streamlines=True,
    show_wake_vortices=False,
    save=True,
)

# Call the animate function on the loaded solver. For a free flight solver, this flies
# the airplane through the scene along its computed trajectory, producing a GIF saved in
# the same directory as this script. Press any key, after orienting the view, to begin
# the animation.
ps.output.animate(
    unsteady_solver=loaded_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
)

# Call the plotting function on the solver. This produces graphs of the loads with
# respect to time. For a free flight solver, it also plots the first Airplane's
# six-degree-of-freedom state history.
ps.output.plot_results_versus_time(
    unsteady_solver=loaded_solver,
    show=True,
    save=True,
)
