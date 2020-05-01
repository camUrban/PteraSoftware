# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp

steady_solver_validation_airplane = asmvp.geometry.Airplane(
    name="Steady Solver Testing Airplane",
    x_ref=0,
    y_ref=0,
    z_ref=0,
    wings=[
        asmvp.geometry.Wing(
            name="Wing",
            x_le=0,
            y_le=0,
            z_le=0,
            symmetric=False,
            cross_sections=[
                asmvp.geometry.WingCrossSection(
                    x_le=0,
                    y_le=0,
                    z_le=0,
                    twist=5,
                    chord=1,
                    airfoil=asmvp.geometry.Airfoil(name="naca0010")
                ),
                asmvp.geometry.WingCrossSection(
                    x_le=-1,
                    y_le=1,
                    z_le=0,
                    twist=10,
                    chord=1,
                    airfoil=asmvp.geometry.Airfoil(name="naca0010")
                )
            ]
        )
    ]
)
