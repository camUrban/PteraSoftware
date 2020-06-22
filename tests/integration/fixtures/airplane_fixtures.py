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
                    ),
                ],
            )
        ],
    )


# ToDo: Properly document this method.
def make_multiple_wing_steady_validation_airplane():
    """

    :return:
    """

    return asmvp.geometry.Airplane(
        # Name the current_airplane.
        name="Multiple Wing Steady Solver Testing Airplane",
        # Define a list of the current_airplane's wings.
        wings=[
            # Initialize the wing object.
            asmvp.geometry.Wing(
                # Name the wing.
                name="Main Wing",
                # This will be a symmetrical wing.
                symmetric=True,
                # Define a list of the wing's cross sections.
                wing_cross_sections=[
                    # Initialize the root cross section object.
                    asmvp.geometry.WingCrossSection(
                        # Define the cross section's twist chord.
                        chord=1.0,
                        # Initialize this cross section's airfoil object.
                        airfoil=asmvp.geometry.Airfoil(name="naca23012"),
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
                        airfoil=asmvp.geometry.Airfoil(name="naca23012"),
                    ),
                ],
            ),
            # Initialize the wing object for the airplane's horizontal stabilizer.
            asmvp.geometry.Wing(
                # Name the wing.
                name="Horizontal Stabilizer",
                # This will be a symmetrical wing.
                symmetric=True,
                # Set this wing to be 5 meters back from the main wing.
                x_le=5.0,
                # Define a list of the wing's cross sections.
                wing_cross_sections=[
                    # Initialize the root cross section object.
                    asmvp.geometry.WingCrossSection(
                        # Define the cross section's twist chord.
                        chord=1.0,
                        # Initialize this cross section's airfoil object.
                        airfoil=asmvp.geometry.Airfoil(name="naca0010"),
                        twist=-5.0,
                    ),
                    # Initialize the tip cross section object.
                    asmvp.geometry.WingCrossSection(
                        # Define the cross section's leading edge placement.
                        x_le=1.0,
                        y_le=1.0,
                        z_le=0.0,
                        # Define the cross section's twist and chord.
                        twist=-5.0,
                        chord=0.75,
                        # Initialize this cross section's airfoil object.
                        airfoil=asmvp.geometry.Airfoil(name="naca0010"),
                    ),
                ],
            ),
        ],
    )


# ToDo: Properly document this method.
def make_asymmetric_unsteady_validation_airplane():
    """

    :return:
    """

    return asmvp.geometry.Airplane(
        name="Unsteady Solver Testing Airplane",
        y_ref=5.0,
        # Define a list of the current_airplane's wings.
        wings=[
            # Initialize the wing object.
            asmvp.geometry.Wing(
                # Name the wing.
                name="Wing",
                symmetric=False,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                # Define a list of the wing's cross sections.
                wing_cross_sections=[
                    # Initialize the root cross section object.
                    asmvp.geometry.WingCrossSection(
                        chord=1.0,
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                        num_spanwise_panels=16,
                    ),
                    # Initialize the tip cross section object.
                    asmvp.geometry.WingCrossSection(
                        y_le=10.0,
                        chord=1.0,
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                        num_spanwise_panels=16,
                    ),
                ],
            )
        ],
    )


# ToDo: Properly document this method.
def make_symmetric_unsteady_validation_airplane():
    """

    :return:
    """

    return asmvp.geometry.Airplane(
        name="Unsteady Solver Testing Airplane",
        # Define a list of the current_airplane's wings.
        wings=[
            # Initialize the wing object.
            asmvp.geometry.Wing(
                # Name the wing.
                name="Wing",
                symmetric=True,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                # Define a list of the wing's cross sections.
                wing_cross_sections=[
                    # Initialize the root cross section object.
                    asmvp.geometry.WingCrossSection(
                        chord=1.0,
                        airfoil=asmvp.geometry.Airfoil(name="naca0012"),
                        num_spanwise_panels=8,
                    ),
                    # Initialize the tip cross section object.
                    asmvp.geometry.WingCrossSection(
                        y_le=5.0,
                        chord=1.0,
                        airfoil=asmvp.geometry.Airfoil(name="naca0012"),
                        num_spanwise_panels=8,
                    ),
                ],
            )
        ],
    )
