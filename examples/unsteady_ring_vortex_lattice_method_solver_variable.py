"""This is script is an example of how to run Ptera Software's
UnsteadyRingVortexLatticeMethodSolver with a custom Airplane with a non static
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
    weight=0.0,
    s_ref=None,
    c_ref=None,
    b_ref=None,
)

# Now define the main wing's root and tip WingCrossSections' WingCrossSectionMovements.
main_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
)
main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
)

# The main Wing was defined to have symmetric=True, mirror_only=False, and with a
# symmetry plane offset non coincident with the Wing's axes yz plane. Therefore,
# that Wing had type 5 symmetry (see the Wing class documentation for more details on
# symmetry types). Therefore, it was actually split into two Wings, the with the
# second Wing being a reflected version of the first. Therefore, we need to define a
# WingMovement for this reflected Wing. To start, we'll first define the reflected
# main wing's root and tip WingCrossSections' WingCrossSectionMovements.
reflected_main_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
)
reflected_main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
)

# Now define the v tail's root and tip WingCrossSections' WingCrossSectionMovements.
v_tail_root_wing_cross_section_movement = (
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
v_tail_tip_wing_cross_section_movement = (
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

# Now define the main wing's WingMovement, the reflected main wing's WingMovement and
# the v tail's WingMovement.
main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_section_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),  # (0.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),  # (0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)
reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[1],
    wing_cross_section_movements=[
        reflected_main_wing_root_wing_cross_section_movement,
        reflected_main_wing_tip_wing_cross_section_movement,
    ],
    ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
    periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
    spacingLer_Gs_Cgs=("sine", "sine", "sine"),
    phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
    ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),  # (0.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),  # (0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)
v_tail_movement = ps.movements.wing_movement.WingMovement(
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
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[main_wing_movement, reflected_main_wing_movement, v_tail_movement],
    ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
    periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
    spacingCg_GP1_CgP1=("sine", "sine", "sine"),
    phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
)

# Delete the extraneous pointers to the WingMovements.
del main_wing_movement
del reflected_main_wing_movement
del v_tail_movement

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=1.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
)

# Define the operating point's OperatingPointMovement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point, periodVCg__E=0.0, spacingVCg__E="sine"
)

# Delete the extraneous pointer.
del example_operating_point

# Define the Movement. This contains the AirplaneMovement and the
# OperatingPointMovement.
movement = ps.movements.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    delta_time=None,
    num_cycles=3,
    num_chords=None,
    num_steps=None,
)

# Delete the extraneous pointers.
del airplane_movement
del operating_point_movement

# Define the UnsteadyProblem.
example_problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

# Define a new solver. The available solver classes are
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver,
# and UnsteadyRingVortexLatticeMethodSolver. We'll create an
# UnsteadyRingVortexLatticeMethodSolver, which requires a UnsteadyProblem.
example_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)

# Delete the extraneous pointer.
del example_problem

# Run the solver.
example_solver.run(
    logging_level="Warning",
    prescribed_wake=True,
)

# Call the animate function on the solver. This produces a GIF of the wake being
# shed. The GIF is saved in the same directory as this script. Press "q",
# after orienting the view, to begin the animation.
ps.output.animate(
    unsteady_solver=example_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=False,
)
