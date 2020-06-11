
# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly document this method.
def make_steady_validation_airplane():
    """

    :return:
    """

    return asmvp.geometry.Airplane(

        # Name the current_airplane.
        name="Steady Solver Testing Airplane",

        # Define a list of the current_airplane's wings.
        wings=[

            # Initialize the wing object.
            asmvp.geometry.Wing(

                # Name the wing.
                name="Wing",

                # This will be a symmetrical wing.
                symmetric=True,

                # Define a list of the wing's cross sections.
                wing_cross_sections=[

                    # Initialize the root cross section object.
                    asmvp.geometry.WingCrossSection(
                        # Define the cross section's twist chord.
                        chord=1.0,

                        # Initialize this cross section's airfoil object.
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    ),

                    # Initialize the tip cross section object.
                    asmvp.geometry.WingCrossSection(

                        # Define the cross section's leading edge placement.
                        x_le=1.0,
                        y_le=5.0,
                        z_le=0.0,

                        # Define the cross section's twist and chord.
                        twist=5.0,
                        chord=0.75,

                        # Initialize this cross section's airfoil object.
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    )
                ]
            )
        ]
    )
