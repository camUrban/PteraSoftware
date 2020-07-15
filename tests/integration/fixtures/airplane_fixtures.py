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
    make_symmetric_multiple_wing_unsteady_validation_airplane: This function creates a multi-wing, symmetric airplane
                                                               object to be used as a fixture for testing unsteady
                                                               solvers.

"""

import pterasoftware as ps


def make_steady_validation_airplane():
    """ This function creates an airplane object to be used as a fixture for testing steady solvers.

    :return steady_validation_airplane: Airplane
        This is the airplane fixture.
    """

    # Create and return the airplane object.
    steady_validation_airplane = ps.geometry.Airplane(
        wings=[
            ps.geometry.Wing(
                symmetric=True,
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=5.0,
                        twist=5.0,
                        chord=0.75,
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
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
    multiple_wing_steady_validation_airplane = ps.geometry.Airplane(
        wings=[
            ps.geometry.Wing(
                symmetric=True,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca23012"),
                        spanwise_spacing="uniform",
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=5.0,
                        chord=0.75,
                        airfoil=ps.geometry.Airfoil(name="naca23012"),
                        spanwise_spacing="uniform",
                    ),
                ],
            ),
            ps.geometry.Wing(
                symmetric=True,
                x_le=5.0,
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca0010"), twist=-5.0,
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=1.0,
                        twist=-5.0,
                        chord=0.75,
                        airfoil=ps.geometry.Airfoil(name="naca0010"),
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
    asymmetric_unsteady_validation_airplane = ps.geometry.Airplane(
        y_ref=5.0,
        wings=[
            ps.geometry.Wing(
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                        num_spanwise_panels=16,
                        spanwise_spacing="cosine",
                        chord=1.0,
                    ),
                    ps.geometry.WingCrossSection(
                        y_le=10.0,
                        chord=1.0,
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                        num_spanwise_panels=16,
                        spanwise_spacing="cosine",
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
    symmetric_unsteady_validation_airplane = ps.geometry.Airplane(
        wings=[
            ps.geometry.Wing(
                symmetric=True,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                        chord=2.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.WingCrossSection(
                        y_le=5.0,
                        chord=2.0,
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                        spanwise_spacing="cosine",
                    ),
                ],
            ),
        ],
    )
    return symmetric_unsteady_validation_airplane


def make_symmetric_multiple_wing_unsteady_validation_airplane():
    """ This function creates a multi-wing, symmetric airplane object to be used as a fixture for testing unsteady
    solvers.

    :return symmetric_multiple_wing_steady_validation_airplane: Airplane
        This is the airplane fixture.
    """

    # Create and return the airplane object.
    symmetric_multiple_wing_steady_validation_airplane = ps.geometry.Airplane(
        wings=[
            ps.geometry.Wing(
                symmetric=True,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=1.0,
                        y_le=5.0,
                        chord=0.75,
                        z_le=0.5,
                        airfoil=ps.geometry.Airfoil(name="naca2412"),
                        spanwise_spacing="cosine",
                    ),
                ],
            ),
            ps.geometry.Wing(
                symmetric=True,
                z_le=1.75,
                x_le=6.25,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca0010"),
                        spanwise_spacing="cosine",
                        twist=-5,
                        x_le=0.25,
                        chord=0.75,
                    ),
                    ps.geometry.WingCrossSection(
                        y_le=1.25,
                        twist=-5,
                        chord=0.5,
                        x_le=0.5,
                        airfoil=ps.geometry.Airfoil(name="naca0010"),
                        spanwise_spacing="cosine",
                    ),
                ],
            ),
            ps.geometry.Wing(
                symmetric=False,
                z_le=0.125,
                x_le=6.25,
                chordwise_spacing="uniform",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        airfoil=ps.geometry.Airfoil(name="naca0010"),
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.WingCrossSection(
                        z_le=1.5,
                        chord=0.75,
                        x_le=0.25,
                        airfoil=ps.geometry.Airfoil(name="naca0010"),
                        spanwise_spacing="cosine",
                    ),
                ],
            ),
        ],
    )

    return symmetric_multiple_wing_steady_validation_airplane
