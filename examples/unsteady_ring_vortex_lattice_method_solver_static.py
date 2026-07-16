"""Demonstrates Ptera Software's UnsteadyRingVortexLatticeMethodSolver with a custom
airplane with static geometry.

The script will log simulation results in a log file.
"""

import logging
from pathlib import Path

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
                        outline_A_Lp=None,
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
                        outline_A_Lp=None,
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
                        outline_A_Lp=None,
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

# Now define the main Wing's root and tip WingCrossSections' WingCrossSectionMovements.
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

# Now define the V-tail's root and tip WingCrossSections' WingCrossSectionMovements.
v_tail_root_wing_cross_section_movement = (
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
v_tail_tip_wing_cross_section_movement = (
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

# Now define the main Wing's WingMovement and the V-tail's WingMovement.
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
    ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
)
v_tail_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[1],
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
del v_tail_root_wing_cross_section_movement
del v_tail_tip_wing_cross_section_movement

# Now define the example Airplane's AirplaneMovement.
airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[main_wing_movement, v_tail_movement],
    ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
    periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
    spacingCg_GP1_CgP1=("sine", "sine", "sine"),
    phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
)

# Delete the extraneous pointers to the WingMovements.
del main_wing_movement
del v_tail_movement

# Define a new OperatingPoint.
example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=1.0, beta=0.0, nu=15.06e-6
)

# Define the OperatingPoint's OperatingPointMovement.
operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point,
    ampVCg__E=0.0,
    periodVCg__E=0.0,
    spacingVCg__E="sine",
    phaseVCg__E=0.0,
)

# Delete the extraneous pointer.
del example_operating_point

# Define the Movement. This contains the AirplaneMovement and the
# OperatingPointMovement.
movement = ps.movements.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    delta_time=None,
    num_cycles=None,
    num_chords=10,
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
# SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver, and
# UnsteadyRingVortexLatticeMethodSolver. We'll create an
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
    prescribed_wake=True,
    show_progress=True,
)

# Save the solved solver to a compressed JSON file. This allows us to load the results
# later without re-running the simulation. Use ".json.gz" for gzip compression, which is
# recommended over plain JSONs for all but the smallest, unmeshed geometry objects.
ps.save("example_solver.json.gz", example_solver)

# Load the saved solver. The loaded object is identical to the original and can be
# passed to any output function.
loaded_solver = ps.load("example_solver.json.gz")

# The load function is annotated as returning object because a saved file can hold any
# Ptera Software object. This assert narrows the type for type checkers and guards
# against loading the wrong file.
assert isinstance(
    loaded_solver,
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
)

ps.output.log_results(solver=loaded_solver)

# Choose where the output files will go. Each output function writes to the working
# directory by default, which is fine for a script that produces a file or two, but this
# one produces six. Sorting them into their own directories keeps them together and
# keeps them from crowding whatever else is in the working directory. The output
# functions do not create a directory that does not exist, so that a mistyped
# destination is reported rather than silently scattering files, which means the script
# creates them itself.
Path("renders").mkdir(exist_ok=True)
Path("plots").mkdir(exist_ok=True)

# Call the draw function on the loaded solver. The path parameter names the file to save
# the drawing to. Press any key to close the plotter after it draws the output.
ps.output.draw(
    solver=loaded_solver,
    scalar_type="lift",
    show_streamlines=True,
    show_wake_vortices=False,
    save=True,
    path="renders/draw.webp",
)

# Call the animate function on the loaded solver. This produces an animated WebP of the
# wake being shed, saved to the path given below. Press any key, after orienting the
# view, to begin the animation.
ps.output.animate(
    unsteady_solver=loaded_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
    path="renders/animate.webp",
)

# Call the plotting function on the solver. This produces graphs of the loads with
# respect to time. It writes several files whose names it derives from each Airplane's
# name and the quantity plotted, so it takes the directory to write them into rather
# than a single path.
ps.output.plot_results_versus_time(
    unsteady_solver=loaded_solver,
    show=True,
    save=True,
    directory="plots",
)
