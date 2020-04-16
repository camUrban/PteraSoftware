# ToDo: Properly document this module.
"""

"""
import aerosandbox_legacy_v0 as asl

import aviansoftwareminimumviableproduct as asmvp

steady_solver_validation_airplane = asl.Airplane(
    name="Steady Solver Testing Airplane",
    xyz_ref=[0, 0, 0],
    wings=[
        asl.Wing(
            name="Wing",
            xyz_le=[0, 0, 0],
            symmetric=False,
            xsecs=[
                asl.WingXSec(
                    xyz_le=[0, 0, 0],
                    twist=5,
                    chord=1,
                    airfoil=asl.Airfoil(name="naca0010")
                ),
                asl.WingXSec(
                    xyz_le=[-1, 1, 0],
                    twist=10,
                    chord=1,
                    airfoil=asl.Airfoil(name="naca0010")
                )
            ]
        )
    ]
)

steady_solver_validation_problem = asmvp.steady_vortex_lattice_method.SteadyVortexLatticeMethod(
    airplane=steady_solver_validation_airplane,
    operating_point=asl.OperatingPoint(
        velocity=10
    )
)

steady_solver_validation_problem.run()
