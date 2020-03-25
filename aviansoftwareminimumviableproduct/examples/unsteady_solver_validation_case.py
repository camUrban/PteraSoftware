from aerosandbox_legacy_v0 import *
from aviansoftwareminimumviableproduct import *
import numpy


unsteady_solver_validation_airplane = Airplane(
    name="Unsteady Solver Validation Airplane",
    xyz_ref=[0.25, 0, 0],
    wings=[
        Wing(
            name="Wing",
            xyz_le=[0, 0, 0],
            symmetric=True,
            xsecs=[
                WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1,
                    airfoil=Airfoil(name="naca2412")
                ),
                WingXSec(
                    xyz_le=[0, 5, 0],
                    chord=1,
                    airfoil=Airfoil(name="naca2412")
                )
            ]
        ),
        Wing(
            name="Horizontal Stabilizer",
            xyz_le=[5, 0, 0.5],
            symmetric=True,
            xsecs=[
                WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1,
                    twist=-5,
                    airfoil=Airfoil(name="naca0012")
                ),
                WingXSec(
                    xyz_le=[0, 1, 0],
                    chord=1,
                    twist=-5,
                    airfoil=Airfoil(name="naca0012")
                )
            ]

        ),
        Wing(
            name="Vertical Stabilizer",
            xyz_le=[5, 0, 0.5],
            symmetric=False,
            xsecs=[
                WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1,
                    airfoil=Airfoil(name="naca0012")
                ),
                WingXSec(
                    xyz_le=[0, 0, 1.5],
                    chord=1,
                    airfoil=Airfoil(name="naca0012")
                )
            ]
        )
    ]
)

unsteady_solver_validation_operating_point = OperatingPoint(
    velocity=10
)

unsteady_solver_validation_movement = Movement(
    movement_period=1,
    sweeping_amplitude=numpy.pi / 2
)

unsteady_solver_validation_problem = UnsteadyAerodynamicsProblem(
    airplane=unsteady_solver_validation_airplane,
    operating_point=unsteady_solver_validation_operating_point,
    movement=unsteady_solver_validation_movement,
    simulation_duration=10,
    simulation_time_step=0.1
)
