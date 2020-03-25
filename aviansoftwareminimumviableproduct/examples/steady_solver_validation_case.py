from aerosandbox_legacy_v0 import *

steady_solver_validation_airplane = Airplane(
    name="Steady Solver Validation Airplane",
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

steady_solver_validation_problem = vlm3(
    airplane=steady_solver_validation_airplane,
    op_point=OperatingPoint(
        velocity=10
    )
)

steady_solver_validation_problem.run()
steady_solver_validation_problem.draw()
