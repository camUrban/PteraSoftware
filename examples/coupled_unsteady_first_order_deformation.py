"""This is script is an example of how to run Ptera Software's
UnsteadyRingVortexLatticeMethodSolver with a custom Airplane with a non-static
Movement."""

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
# of spanwise panels and ensure that in Wing you set the single_step_wing parameter to True,
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
            spanwise_spacing="cosine" if i < len(cross_section_chords) - 1 else None,
            airfoil=ps.geometry.airfoil.Airfoil(
                name="naca2412",
                outline_A_lp=None,
                resample=True,
                n_points_per_side=400,
            ),
        )
    )

# Primary wing definition. Note that the single_step_wing parameter is set to True,
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
    single_step_wing=True,
    num_chordwise_panels=6,
    chordwise_spacing="uniform",
)

# Actually generating the airplane. A tail is added to the airplane, but it is not
# split into strips for deformation as currently only the first wing is considered
# for deformation in the codebase. Fututre versions of this feature could allow for
# the deformation of multiple wings. For now, it is convenient to not split the tail
# into single strip wing cross sections as it reduces the number of movement variables
# that need to be defined.
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
            single_step_wing=False,
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

# defintions for wing cross section movement parameters
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
main_movements_list = []
main_single_step_movements_list = []

# A list of movements for the reflected wing
reflected_movements_list = []
reflected_single_step_movements_list = []

# A loop for defining the movement for the main wing and its reflected counterpart's wing
# cross sections. Here, we are defining single step wing cross movement, a movement class
# that functions differently from the standard Movement class, by giving the next
# position of wing cross section from the previous instead of attempting to precompute
# the entire movement beforehand as that is impossible in scenarios where the deformation
# is dependent on the aerodynamic loads.
for i in range(len(example_airplane.wings[0].wing_cross_sections)):
    if i == 0:
        single_step_movement = ps.movements.single_step.single_step_wing_cross_section_movement.SingleStepWingCrossSectionMovement(
            base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[i],
        )
        main_single_step_movements_list.append(single_step_movement)
        reflected_single_step_movements_list.append(single_step_movement)

    else:
        single_step_movement = ps.movements.single_step.single_step_wing_cross_section_movement.SingleStepWingCrossSectionMovement(
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
        main_single_step_movements_list.append(single_step_movement)
        reflected_single_step_movements_list.append(single_step_movement)


# Now define the v-tail's root and tip WingCrossSections' WingCrossSectionMovements.
single_step_v_tail_root_wing_cross_section_movement = ps.movements.single_step.single_step_wing_cross_section_movement.SingleStepWingCrossSectionMovement(
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
single_step_v_tail_tip_wing_cross_section_movement = ps.movements.single_step.single_step_wing_cross_section_movement.SingleStepWingCrossSectionMovement(
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

# This dephase parameter is used to make the wing start in a flat position
dephase = 169.0

# Now define the main wing's SingleStepWingMovement, the reflected main wing's SingleStepWingMovement and
# the v-tail's SingleStepWingMovement.
single_step_main_wing_movement = (
    ps.movements.single_step.single_step_wing_movement.SingleStepWingMovement(
        base_wing=example_airplane.wings[0],
        single_step_wing_cross_section_movements=main_single_step_movements_list,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),  # (0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),  # (0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(dephase, 0.0, 0.0),
    )
)

single_step_reflected_main_wing_movement = (
    ps.movements.single_step.single_step_wing_movement.SingleStepWingMovement(
        base_wing=example_airplane.wings[1],
        single_step_wing_cross_section_movements=reflected_single_step_movements_list,
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

single_step_v_tail_movement = (
    ps.movements.single_step.single_step_wing_movement.SingleStepWingMovement(
        base_wing=example_airplane.wings[2],
        single_step_wing_cross_section_movements=[
            single_step_v_tail_root_wing_cross_section_movement,
            single_step_v_tail_tip_wing_cross_section_movement,
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
)

# Delete the extraneous pointers to the WingCrossSectionMovements, as these are now
# contained within the WingMovements. This is optional, but it can make debugging
# easier.
del single_step_v_tail_root_wing_cross_section_movement
del single_step_v_tail_tip_wing_cross_section_movement

# Now define the example airplane's SingleStepAirplaneMovement.
single_step_airplane_movement = (
    ps.movements.single_step.single_step_airplane_movement.SingleStepAirplaneMovement(
        base_airplane=example_airplane,
        single_step_wing_movements=[
            single_step_main_wing_movement,
            single_step_reflected_main_wing_movement,
            single_step_v_tail_movement,
        ],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )
)

# Delete the extraneous pointers to the WingMovements.
del single_step_main_wing_movement
del single_step_reflected_main_wing_movement
del single_step_v_tail_movement

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define the operating point's OperatingPointMovement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point, periodVCg__E=0.0, spacingVCg__E="sine"
)

single_step_operating_point_movement = ps.movements.single_step.single_step_operating_point_movement.SingleStepOperatingPointMovement(
    base_operating_point=example_operating_point,
    ampVCg__E=0.0,
    periodVCg__E=0.0,
    spacingVCg__E="sine",
)

# Delete the extraneous pointer.
del example_operating_point

# Define the SingleStepMovement. This contains the SingleStepAirplaneMovement and the
# SingleStepOperatingPointMovement.

single_step_movement = ps.movements.single_step.single_step_movement.SingleStepMovement(
    single_step_airplane_movements=[single_step_airplane_movement],
    single_step_operating_point_movement=single_step_operating_point_movement,
    delta_time=0.03,
    num_cycles=3,
)

# Delete the extraneous pointers.
del operating_point_movement

# Define the UnsteadyProblem.
# The deformation parameters are set here
# The wing_density, spring_constant and damping_constant are the primary parameters
# you should expect to change. The rest are more for considering numerical issues
# with our integrator and debugging. Plotting the flap cycle can give good data as well.
example_problem = ps.problems.AeroelasticUnsteadyProblem(
    single_step_movement=single_step_movement,
    wing_density=0.012,
    spring_constant=5.0,
    damping_constant=1.0,
    aero_scaling=1.0,
    moment_scaling_factor=1.0,
    damping_eps=1e-3,
    plot_flap_cycle=False,
    custom_spacing_second_derivative=None,
    only_final_results=False,
)

# Define a new solver. The available solver classes are
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver,
# and UnsteadyRingVortexLatticeMethodSolver. We'll create an
# UnsteadyRingVortexLatticeMethodSolver, which requires a UnsteadyProblem.
example_solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
    coupled_unsteady_problem=example_problem,
)

# Delete the extraneous pointer.
del example_problem

# Run the solver.
example_solver.run(
    prescribed_wake=True,
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
