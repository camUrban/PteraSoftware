"""This script is an example of how to run Ptera Software's
AeroelasticUnsteadyRingVortexLatticeMethodSolver with a flapping airplane whose main
wing deforms under its own aerodynamic loads."""

# First, import the software's main package. Note that if you wished to import this
# software into another package, you would first install it by running "pip install
# pterasoftware" in your terminal.
import pterasoftware as ps

# Create an Airplane with our custom geometry. I am going to declare every parameter
# for Airplane, even though most of them have usable default values. This is for
# educational purposes, but keep in mind that it makes the code much longer than it
# needs to be. For details about each parameter, read the detailed class docstring.
# The same caveats apply to the other classes, methods, and functions I call in this
# script.

# Wing cross section initialization
# offsets for the spacing
num_spanwise_panels = 2
Lp_Wcsp_Lpp_Offsets = (0.1, 0.5, 0.0)
cross_section_chords = [1.75, 1.75, 1.75, 1.75, 1.65, 1.55, 1.4, 1.2, 1.0]
wing_cross_sections = []

# Initialization loop for our wing cross sections. Here we are defining automatically
# wing cross sections with a variable set of chords. All of the wing cross sections for
# deformation simulation are defined to have num_spanwise_panels=1 (except the wing tip which
# is always None). This is because we deform each strip of wing cross section independently by
# modeling them as torsional springs, and that model only really works if those strips are thin.
# Note that if you want to go thinner for the same base definition, you can increase the number
# of spanwise panels and ensure that in Wing you set the explode_into_strips parameter to True,
# which will ensure that the wing is split back up into single strips for deformation.
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
            spanwise_spacing="uniform" if i < len(cross_section_chords) - 1 else None,
            airfoil=ps.geometry.airfoil.Airfoil(
                name="naca2412",
                outline_A_lp=None,
                resample=True,
                n_points_per_side=400,
            ),
        )
    )

# Primary wing definition. Note that the explode_into_strips parameter is set to True,
# which means that the wing will be split into strips for deformation, and each
# strip will be modeled as a torsional spring.
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

# Actually generate the airplane. The V-tail is added as a second lifting surface but
# is not split into deformation strips, because it follows prescribed rigid motion
# rather than deforming. The solver deforms each wing from its own aerodynamic loads,
# but only wings backed by an AeroelasticWingMovement; the V-tail uses a standard
# WingMovement and stays rigid. Leaving the tail rigid also keeps its movement
# definition simpler, since it then needs no per-strip cross sections.
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

# The main Wing was defined to have symmetric=True, mirror_only=False, and with a
# symmetry plane offset non-coincident with the Wing's axes yz-plane. Therefore,
# that Wing had type 5 symmetry (see the Wing class documentation for more details on
# symmetry types). Therefore, it was actually split into two Wings, the with the
# second Wing being a reflected version of the first. Therefore, we need to define a
# WingMovement for this reflected Wing. To start, we'll first define the reflected
# main wing's root and tip WingCrossSections' WingCrossSectionMovements.

# definitions for wing cross section movement parameters
dephase_x = 0.0
period_x = 0.0
amplitude_x = 0.0

dephase_y = 0.0
period_y = 0.0
amplitude_y = 0.0

dephase_z = 0.0
period_z = 0.0
amplitude_z = 0.0

# A list of wing cross section movements for the main wing
main_wcs_movements_list = []

# A list of wing cross section movements for the reflected wing
reflected_wcs_movements_list = []

# A loop for defining the movement for the main wing and its reflected counterpart's wing
# cross sections. Each wing cross section has its own AeroelasticWingCrossSectionMovement
# which allows the solver to apply deformation angles at each time step based on the
# aerodynamic loads.
for i in range(len(example_airplane.wings[0].wing_cross_sections)):
    if i == 0:
        wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[i],
        )
        main_wcs_movements_list.append(wcs_movement)
        reflected_wcs_movements_list.append(wcs_movement)

    else:
        wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[i],
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(amplitude_x, amplitude_y, amplitude_z),
            periodAngles_Wcsp_to_Wcs_ixyz=(period_x, period_y, period_z),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(dephase_x, dephase_y, dephase_z),
        )
        main_wcs_movements_list.append(wcs_movement)
        reflected_wcs_movements_list.append(wcs_movement)


# Now define the v-tail's root and tip WingCrossSections' WingCrossSectionMovements.
# The V-tail is not an aeroelastic surface, so we use standard WingCrossSectionMovements
# and a standard WingMovement. This keeps the V-tail mesh consistent across all time
# steps without applying any structural deformation.
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
# Reflected V-tail WingCrossSectionMovements reuse the same static motion as the
# original V-tail. Both halves are symmetric and neither deforms.

# This dephase parameter is used to make the wing start in a flat position
dephase = 169.0

# Now define the main wing's AeroelasticWingMovement, the reflected main wing's
# AeroelasticWingMovement, and the v-tail's AeroelasticWingMovement.
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

# Delete the extraneous pointers to the WingCrossSectionMovements, as these are now
# contained within the WingMovements. This is optional, but it can make debugging
# easier.
del v_tail_root_wcs_movement
del v_tail_tip_wcs_movement

# Now define the example airplane's AeroelasticAirplaneMovement.
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

del main_wing_movement
del reflected_main_wing_movement
del v_tail_wing_movement

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define the operating point's AeroelasticOperatingPointMovement.
example_operating_point_movement = (
    ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
        base_operating_point=example_operating_point,
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
    )
)

# Delete the extraneous pointer.
del example_operating_point

# Define the AeroelasticMovement. This contains the AeroelasticAirplaneMovement and
# the AeroelasticOperatingPointMovement. The delta_time and num_steps must be specified
# explicitly. With a flapping period of 1.0s, 3 cycles at dt=0.03s gives 100 steps.
example_movement = ps.movements.aeroelastic_movement.AeroelasticMovement(
    airplane_movements=[example_airplane_movement],
    operating_point_movement=example_operating_point_movement,
    delta_time=0.03,
    num_steps=100,
)

# Define the AeroelasticUnsteadyProblem.
# The deformation parameters are set here.
# The wing_density, spring_constant and damping_constant are the primary parameters
# you should expect to change. The rest are more for considering numerical issues
# with our integrator and debugging. Plotting the flap cycle can give good data as well.
example_problem = ps.problems.AeroelasticUnsteadyProblem(
    movement=example_movement,
    wing_density=0.012,
    spring_constant=10.0,
    damping_constant=1.0,
    aero_scaling=1.0,
    step_discards=5,
    moment_scaling_factor=1.0,
    plot_flap_cycle=False,
)

# Define a new solver. We'll create an AeroelasticUnsteadyRingVortexLatticeMethodSolver,
# which requires an AeroelasticUnsteadyProblem.
example_solver = ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver(
    aeroelastic_unsteady_problem=example_problem,
)

# Delete the extraneous pointer.
del example_problem

# Run the solver.
example_solver.run(
    prescribed_wake=False,
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
