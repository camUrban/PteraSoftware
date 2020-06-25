""" This module creates airplane objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_steady_validation_airplane: This function creates an airplane object to be used as a fixture for testing steady
                                     solvers.
    make_multiple_wing_steady_validation_airplane: This function creates a multi-wing airplane object to be used as a
                                                   fixture for testing steady solvers.
    make_asymmetric_unsteady_validation_airplane: This function creates an asymmetric airplane object to be used as a
                                                  fixture for testing unsteady solvers.
    make_symmetric_unsteady_validation_airplane: This function creates a symmetric airplane object to be used as a
                                                 fixture for testing unsteady solvers.

"""

import aviansoftwareminimumviableproduct as asmvp


def make_steady_validation_airplane():
    """ This function creates an airplane object to be used as a fixture for testing steady solvers.

    :return steady_validation_airplane: Airplane
        This is the airplane fixture.
    """

    # Create and return the airplane object.
    steady_validation_airplane = asmvp.geometry.Airplane(
        wings=[
            asmvp.geometry.Wing(
                symmetric=True,
                wing_cross_sections=[
                    asmvp.geometry.WingCrossSection(
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    ),
                    asmvp.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=5.0,
                        twist=5.0,
                        chord=0.75,
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                    ),
                ],
            )
        ],
    )
    return steady_validation_airplane


def make_multiple_wing_steady_validation_airplane():
    """ This function creates a multi-wing airplane object to be used as a fixture for testing steady solvers.

    :return multiple_wing_steady_validation_airplane: Airplane
        This is the airplane fixture.
    """

    # Create and return the airplane object.
    multiple_wing_steady_validation_airplane = asmvp.geometry.Airplane(
        wings=[
            asmvp.geometry.Wing(
                symmetric=True,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    asmvp.geometry.WingCrossSection(
                        airfoil=asmvp.geometry.Airfoil(name="naca23012"),
                        spanwise_spacing="uniform",
                    ),
                    asmvp.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=5.0,
                        chord=0.75,
                        airfoil=asmvp.geometry.Airfoil(name="naca23012"),
                        spanwise_spacing="uniform",
                    ),
                ],
            ),
            asmvp.geometry.Wing(
                symmetric=True,
                x_le=5.0,
                wing_cross_sections=[
                    asmvp.geometry.WingCrossSection(
                        airfoil=asmvp.geometry.Airfoil(name="naca0010"), twist=-5.0,
                    ),
                    asmvp.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=1.0,
                        twist=-5.0,
                        chord=0.75,
                        airfoil=asmvp.geometry.Airfoil(name="naca0010"),
                    ),
                ],
            ),
        ],
    )
    return multiple_wing_steady_validation_airplane


def make_asymmetric_unsteady_validation_airplane():
    """ This function creates an asymmetric airplane object to be used as a fixture for testing unsteady solvers.

    :return asymmetric_unsteady_validation_airplane: Airplane
        This is the airplane fixture.
    """

    # Create and return the airplane object.
    asymmetric_unsteady_validation_airplane = asmvp.geometry.Airplane(
        y_ref=5.0,
        wings=[
            asmvp.geometry.Wing(
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    asmvp.geometry.WingCrossSection(
                        airfoil=asmvp.geometry.Airfoil(name="naca2412"),
                        num_spanwise_panels=16,
                    ),
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
    return asymmetric_unsteady_validation_airplane


def make_symmetric_unsteady_validation_airplane():
    """ This function creates a symmetric airplane object to be used as a fixture for testing unsteady solvers.

    :return symmetric_unsteady_validation_airplane: Airplane
        This is the airplane fixture.
    """

    # Create and return the airplane object.
    symmetric_unsteady_validation_airplane = asmvp.geometry.Airplane(
        wings=[
            asmvp.geometry.Wing(
                symmetric=True,
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    asmvp.geometry.WingCrossSection(
                        airfoil=asmvp.geometry.Airfoil(name="naca0012"),
                    ),
                    asmvp.geometry.WingCrossSection(
                        y_le=5.0, airfoil=asmvp.geometry.Airfoil(name="naca0012"),
                    ),
                ],
            )
        ],
    )
    return symmetric_unsteady_validation_airplane
