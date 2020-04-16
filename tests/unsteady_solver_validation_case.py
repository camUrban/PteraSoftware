# ToDo: Properly document this module.
"""

"""
import aerosandbox_legacy_v0 as asl
import numpy as np

import aviansoftwareminimumviableproduct as asmvp

unsteady_solver_validation_airplane = asl.Airplane(
    name="Unsteady Solver Validation Airplane",
    xyz_ref=[0.25, 0, 0],
    wings=[
        asl.Wing(
            name="Wing",
            xyz_le=[0, 0, 0],
            symmetric=True,
            xsecs=[
                asl.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1,
                    airfoil=asl.Airfoil(name="naca2412")
                ),
                asl.WingXSec(
                    xyz_le=[0, 5, 0],
                    chord=1,
                    airfoil=asl.Airfoil(name="naca2412")
                )
            ]
        )
    ]
)

unsteady_solver_validation_operating_point = asl.OperatingPoint(
    velocity=10
)

unsteady_solver_validation_movement = asmvp.movement.Movement(
    movement_period=0.5,
    sweeping_amplitude=np.pi / 4
)

unsteady_solver_validation_problem = asmvp.unsteady_vortex_lattice_method.UnsteadyVortexLatticeMethod(
    unsteady_solver_validation_airplane,
    operating_point=unsteady_solver_validation_operating_point,
    movement=unsteady_solver_validation_movement,
    simulation_duration=1,
    simulation_time_step=0.01
)

unsteady_solver_validation_problem.run()
