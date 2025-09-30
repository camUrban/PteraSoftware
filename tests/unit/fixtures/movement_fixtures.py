"""This module contains functions to create movement objects for use in tests.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_sine_spacing_Lp_wing_cross_section_movement_fixture: This method makes
    a fixture that is a WingCrossSectionMovement with sine spacing for
    Lp_Wcsp_Lpp.

    make_uniform_spacing_Lp_wing_cross_section_movement_fixture: This method
    makes a fixture that is a WingCrossSectionMovement with uniform spacing for
    Lp_Wcsp_Lpp.

    make_mixed_spacing_Lp_wing_cross_section_movement_fixture: This method makes
    a fixture that is a WingCrossSectionMovement with mixed spacing for
    Lp_Wcsp_Lpp.

    make_sine_spacing_angles_wing_cross_section_movement_fixture: This method
    makes a fixture that is a WingCrossSectionMovement with sine spacing for
    angles_Wcsp_to_Wcs_izyx.

    make_uniform_spacing_angles_wing_cross_section_movement_fixture: This method
    makes a fixture that is a WingCrossSectionMovement with uniform spacing for
    angles_Wcsp_to_Wcs_izyx.

    make_mixed_spacing_angles_wing_cross_section_movement_fixture: This method
    makes a fixture that is a WingCrossSectionMovement with mixed spacing for
    angles_Wcsp_to_Wcs_izyx.
"""

import pterasoftware as ps
from . import geometry_fixtures


def make_sine_spacing_Lp_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with sine
    spacing for Lp_Wcsp_Lpp.

    :return sine_spacing_Lp_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with sine spacing for Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = (
        geometry_fixtures.make_root_wing_cross_section_fixture()
    )

    # Create the WingCrossSectionMovement with sine spacing.
    sine_spacing_Lp_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
    )

    # Delete the constructing fixture.
    del base_wing_cross_section

    # Return the WingCrossSectionMovement fixture.
    return sine_spacing_Lp_wing_cross_section_movement_fixture


def make_uniform_spacing_Lp_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with uniform
    spacing for Lp_Wcsp_Lpp.

    :return uniform_spacing_Lp_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with uniform spacing for
        Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = (
        geometry_fixtures.make_root_wing_cross_section_fixture()
    )

    # Create the WingCrossSectionMovement with uniform spacing.
    uniform_spacing_Lp_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("uniform", "uniform", "uniform"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
    )

    # Delete the constructing fixture.
    del base_wing_cross_section

    # Return the WingCrossSectionMovement fixture.
    return uniform_spacing_Lp_wing_cross_section_movement_fixture


def make_mixed_spacing_Lp_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with mixed
    spacing for Lp_Wcsp_Lpp.

    :return mixed_spacing_Lp_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with mixed spacing for Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    # Use tip fixture which has Lp_Wcsp_Lpp[1] = 2.0, allowing for amplitude of 1.5.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the WingCrossSectionMovement with mixed spacing.
    mixed_spacing_Lp_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(1.0, 1.5, 0.5),
            periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
            spacingLp_Wcsp_Lpp=("sine", "uniform", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
    )

    # Delete the constructing fixture.
    del base_wing_cross_section

    # Return the WingCrossSectionMovement fixture.
    return mixed_spacing_Lp_wing_cross_section_movement_fixture


def make_sine_spacing_angles_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with sine
    spacing for angles_Wcsp_to_Wcs_izyx.

    :return sine_spacing_angles_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with sine spacing for
        angles_Wcsp_to_Wcs_izyx.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = (
        geometry_fixtures.make_root_wing_cross_section_fixture()
    )

    # Create the WingCrossSectionMovement with sine spacing for angles.
    sine_spacing_angles_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_izyx=(10.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_izyx=(1.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Delete the constructing fixture.
    del base_wing_cross_section

    # Return the WingCrossSectionMovement fixture.
    return sine_spacing_angles_wing_cross_section_movement_fixture


def make_uniform_spacing_angles_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with uniform
    spacing for angles_Wcsp_to_Wcs_izyx.

    :return uniform_spacing_angles_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with uniform spacing for
        angles_Wcsp_to_Wcs_izyx.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = (
        geometry_fixtures.make_root_wing_cross_section_fixture()
    )

    # Create the WingCrossSectionMovement with uniform spacing for angles.
    uniform_spacing_angles_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_izyx=(10.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_izyx=(1.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_izyx=("uniform", "uniform", "uniform"),
            phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Delete the constructing fixture.
    del base_wing_cross_section

    # Return the WingCrossSectionMovement fixture.
    return uniform_spacing_angles_wing_cross_section_movement_fixture


def make_mixed_spacing_angles_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with mixed
    spacing for angles_Wcsp_to_Wcs_izyx.

    :return mixed_spacing_angles_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with mixed spacing for
        angles_Wcsp_to_Wcs_izyx.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = (
        geometry_fixtures.make_root_wing_cross_section_fixture()
    )

    # Create the WingCrossSectionMovement with mixed spacing for angles.
    mixed_spacing_angles_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_izyx=(10.0, 20.0, 5.0),
            periodAngles_Wcsp_to_Wcs_izyx=(1.0, 1.0, 1.0),
            spacingAngles_Wcsp_to_Wcs_izyx=("sine", "uniform", "sine"),
            phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Delete the constructing fixture.
    del base_wing_cross_section

    # Return the WingCrossSectionMovement fixture.
    return mixed_spacing_angles_wing_cross_section_movement_fixture