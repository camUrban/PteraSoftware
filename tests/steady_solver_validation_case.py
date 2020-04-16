# ToDo: Properly document this module.
"""

"""
import aerosandbox_legacy_v0 as asl

import aviansoftwareminimumviableproduct as asmvp

steady_solver_validation_airplane = asl.Airplane(
    name="Steady Solver Validation Airplane",
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
        ),
        asl.Wing(
            name="Horizontal Stabilizer",
            xyz_le=[5, 0, 0.5],
            symmetric=True,
            xsecs=[
                asl.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1,
                    twist=-5,
                    airfoil=asl.Airfoil(name="naca0012")
                ),
                asl.WingXSec(
                    xyz_le=[0, 1, 0],
                    chord=1,
                    twist=-5,
                    airfoil=asl.Airfoil(name="naca0012")
                )
            ]

        ),
        asl.Wing(
            name="Vertical Stabilizer",
            xyz_le=[5, 0, 0.5],
            symmetric=False,
            xsecs=[
                asl.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=1,
                    airfoil=asl.Airfoil(name="naca0012")
                ),
                asl.WingXSec(
                    xyz_le=[0, 0, 1.5],
                    chord=1,
                    airfoil=asl.Airfoil(name="naca0012")
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
