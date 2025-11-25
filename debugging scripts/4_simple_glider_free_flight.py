"""This script is for running an unsteady free flight simulation on the simple glider
after finding its converged and trimmed parameters.

This example creates a simple rectangular Wing and couples Ptera Software's
aerodynamic solver with MuJoCo's rigid body dynamics to simulate free flight. This is
the simplest possible test case for the coupled solver."""

import logging
import math

import numpy as np
import matplotlib.pyplot as plt

import pterasoftware as ps

# noinspection PyProtectedMember
import pterasoftware._transformations as _transformations


# REFACTOR: This function (and its use) is duplicated in flat_plate_free_flight.py.
#  It should somehow be incorporated into the MuJoCoModel class or otherwise
#  standardized and then tested with unit tests.
def R_to_quat_wxyz(R):
    # R is 3x3 list/ndarray
    m00, m01, m02 = R[0]
    m10, m11, m12 = R[1]
    m20, m21, m22 = R[2]
    t = m00 + m11 + m22
    if t > 0.0:
        s = math.sqrt(t + 1.0) * 2.0
        qw = 0.25 * s
        qx = (m21 - m12) / s
        qy = (m02 - m20) / s
        qz = (m10 - m01) / s
    elif (m00 > m11) and (m00 > m22):
        s = math.sqrt(1.0 + m00 - m11 - m22) * 2.0
        qw = (m21 - m12) / s
        qx = 0.25 * s
        qy = (m01 + m10) / s
        qz = (m02 + m20) / s
    elif m11 > m22:
        s = math.sqrt(1.0 + m11 - m00 - m22) * 2.0
        qw = (m02 - m20) / s
        qx = (m01 + m10) / s
        qy = 0.25 * s
        qz = (m12 + m21) / s
    else:
        s = math.sqrt(1.0 + m22 - m00 - m11) * 2.0
        qw = (m10 - m01) / s
        qx = (m02 + m20) / s
        qy = (m12 + m21) / s
        qz = 0.25 * s
    # normalize to be safe
    n = math.sqrt(qw * qw + qx * qx + qy * qy + qz * qz)
    return qw / n, qx / n, qy / n, qz / n


# Configure logging to display INFO messages.
# Some libraries configure logging before we can, so we need to directly set the
# root logger level and update any existing handlers.
logging.root.setLevel(logging.INFO)
for handler in logging.root.handlers:
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("[%(name)s] %(levelname)s: %(message)s"))

script_logger = logging.getLogger("script")
script_logger.setLevel(logging.INFO)

converged_prescribed_wake = True
converged_num_chords = 13
converged_num_chordwise_panels = 6
wing_1_converged_num_spanwise_panels = 30
wing_2_converged_num_spanwise_panels = 6
wing_3_converged_num_spanwise_panels = 12

trim_vCg__E = 12.9
trim_alpha = 3.3
trim_beta = 0.0
trim_externalFX_W = 6.5

delta_time = 0.012920
converged_num_steps = 78

free_num_steps = 20
this_g = -9.80665
this_airplane_name = "Simple Glider"
this_initial_key_frame_name = "Start"
this_weight = 420

this_mass = this_weight / abs(this_g)

R_act_E_to_B_0 = _transformations.generate_rot_T(
    angles=(0.0, trim_alpha, trim_beta), passive=False, intrinsic=True, order="zyx"
)[:3, :3]
quat0w, quat0x, quat0y, quat0z = R_to_quat_wxyz(R_act_E_to_B_0)

mujoco_xml = f"""
<mujoco model="{this_airplane_name}">
  <option timestep="{delta_time:.6f}" integrator="RK4" gravity="0.0 0.0 0.0"/>

  <worldbody>
    <light pos="0.0 0.0 10.0"/>

    <body name="{this_airplane_name}" pos="0.0 0.0 0.0" >
      <freejoint/>
      <inertial pos="0.0 0.0 0.0" mass="{this_mass:.6f}" fullinertia="155.614 398.513 508.699 0 -45.658 0"/>
    </body>
  </worldbody>
  
  <keyframe>
    <key name="{this_initial_key_frame_name}" qpos="0.0 0.0 0.0 {quat0w:.6f} {quat0x:.6f} {quat0y:.6f} {quat0z:.6f}" qvel="{trim_vCg__E:.6f} 0.0 0.0 0.0 0.0 0.0" />
  </keyframe>
</mujoco>
"""

# Create the MuJoCo model wrapper.
simple_glider_mujoco_model = ps.mujoco_model.MuJoCoModel(
    xml=mujoco_xml,
    body_name=this_airplane_name,
    delta_time=delta_time,
    initial_key_frame_name=this_initial_key_frame_name,
)

simple_glider_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=wing_1_converged_num_spanwise_panels,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            chordwise_spacing="uniform",
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=wing_2_converged_num_spanwise_panels,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 0.5),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=wing_3_converged_num_spanwise_panels,
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 2.0, 0.0),
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 1.0),
            angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
    ],
    weight=this_weight,
)

simple_glider_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=simple_glider_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[2],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

del simple_glider_airplane

simple_glider_coupled_operating_point = ps.operating_point.CoupledOperatingPoint(
    vCg__E=trim_vCg__E,
    alpha=trim_alpha,
    beta=trim_beta,
    externalFX_W=trim_externalFX_W,
)

simple_glider_coupled_movement = ps.movements.movement.CoupledMovement(
    airplane_movement=simple_glider_airplane_movement,
    initial_coupled_operating_point=simple_glider_coupled_operating_point,
    delta_time=delta_time,
    prescribed_num_steps=converged_num_steps,
    free_num_steps=free_num_steps,
)

del simple_glider_airplane_movement
del simple_glider_coupled_operating_point

simple_glider_coupled_unsteady_problem = ps.problems.CoupledUnsteadyProblem(
    coupled_movement=simple_glider_coupled_movement,
)

del simple_glider_coupled_movement

simple_glider_coupled_solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
    coupled_unsteady_problem=simple_glider_coupled_unsteady_problem,
    mujoco_model=simple_glider_mujoco_model,
)

# Run the coupled simulation.
script_logger.info("Starting free flight simulation.")
script_logger.info(
    f"    Time steps: {simple_glider_coupled_unsteady_problem.num_steps}"
)
script_logger.info(
    f"    Delta time: {simple_glider_coupled_unsteady_problem.delta_time} s"
)
script_logger.info(
    f"    Total duration: "
    f"{simple_glider_coupled_solver.num_steps * simple_glider_coupled_solver.delta_time} s"
)

del simple_glider_coupled_unsteady_problem

script_logger.info("Free flight simulation completed successfully.")

simple_glider_coupled_solver.run(
    logging_level="Warning",
    prescribed_wake=converged_prescribed_wake,
)

times = np.linspace(
    0,
    (free_num_steps + converged_num_steps - 1) * delta_time,
    (free_num_steps + converged_num_steps),
)
plt.plot(times, simple_glider_coupled_solver.stackPosition_E_E[:, 0])
plt.title("CgP1_X")
plt.show()

plt.plot(times, simple_glider_coupled_solver.stackPosition_E_E[:, 2])
plt.title("CgP1_Z")
plt.show()

ps.output.animate_free_flight(
    coupled_solver=simple_glider_coupled_solver,
    scalar_type="lift",
    show_wake_vortices=True,
)

ps.output.draw(
    solver=simple_glider_coupled_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=True,
)
