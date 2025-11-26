"""This is a debugging script for running an unsteady free flight simulation on a
flat plate with one Panel."""

import logging

import numpy as np
import matplotlib.pyplot as plt

import pterasoftware as ps

# Configure logging to display INFO messages.
# Some libraries configure logging before we can, so we need to directly set the
# root logger level and update any existing handlers.
logging.root.setLevel(logging.INFO)
for handler in logging.root.handlers:
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("[%(name)s] %(levelname)s: %(message)s"))

script_logger = logging.getLogger("script")
script_logger.setLevel(logging.INFO)

trim_vCg__E = 10.0
trim_alpha = 5.0
trim_beta = 0.0
trim_externalFX_W = 0.0

delta_time = 0.1
converged_num_steps = 5

free_num_steps = 40
this_g = -9.80665
this_airplane_name = "Simple Glider"
this_initial_key_frame_name = "Start"
this_weight = 100

this_mass = this_weight / abs(this_g)
I_BP1_CgP1 = np.array(
    [[100.0, 0.0, 0.0], [0.0, 400.0, 0.0], [0.0, 0.0, 400.0]], dtype=float
)

flat_plate_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=1,
                    control_surface_symmetry_type=None,
                    spanwise_spacing="uniform",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type=None,
                    spanwise_spacing=None,
                ),
            ],
            chordwise_spacing="uniform",
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=1,
        ),
    ],
    weight=this_weight,
)

flat_plate_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=flat_plate_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=flat_plate_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=flat_plate_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=flat_plate_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

del flat_plate_airplane

flat_plate_coupled_operating_point = ps.operating_point.CoupledOperatingPoint(
    vCg__E=trim_vCg__E,
    alpha=trim_alpha,
    beta=trim_beta,
    externalFX_W=trim_externalFX_W,
)

flat_plate_coupled_movement = ps.movements.movement.CoupledMovement(
    airplane_movement=flat_plate_airplane_movement,
    initial_coupled_operating_point=flat_plate_coupled_operating_point,
    delta_time=delta_time,
    prescribed_num_steps=converged_num_steps,
    free_num_steps=free_num_steps,
)

del flat_plate_airplane_movement
del flat_plate_coupled_operating_point

flat_plate_coupled_unsteady_problem = ps.problems.CoupledUnsteadyProblem(
    coupled_movement=flat_plate_coupled_movement, I_BP1_CgP1=I_BP1_CgP1
)

del flat_plate_coupled_movement

flat_plate_coupled_solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
    coupled_unsteady_problem=flat_plate_coupled_unsteady_problem,
)

# Run the coupled simulation.
script_logger.info("Starting free flight simulation.")
script_logger.info(f"    Time steps: {flat_plate_coupled_unsteady_problem.num_steps}")
script_logger.info(
    f"    Delta time: {flat_plate_coupled_unsteady_problem.delta_time} s"
)
script_logger.info(
    f"    Total duration: "
    f"{flat_plate_coupled_solver.num_steps * flat_plate_coupled_solver.delta_time} s"
)

del flat_plate_coupled_unsteady_problem

script_logger.info("Free flight simulation completed successfully.")

flat_plate_coupled_solver.run(
    logging_level="Info",
    prescribed_wake=True,
)

times = np.linspace(
    0,
    (free_num_steps + converged_num_steps - 1) * delta_time,
    (free_num_steps + converged_num_steps),
)
plt.plot(times, flat_plate_coupled_solver.stackPosition_E_E[:, 0])
plt.title("CgP1_X")
plt.show()

plt.plot(times, flat_plate_coupled_solver.stackPosition_E_E[:, 2])
plt.title("CgP1_Z")
plt.show()

ps.output.animate_free_flight(
    coupled_solver=flat_plate_coupled_solver,
    scalar_type="lift",
    show_wake_vortices=True,
)

ps.output.draw(
    solver=flat_plate_coupled_solver, scalar_type="lift", show_wake_vortices=True
)
