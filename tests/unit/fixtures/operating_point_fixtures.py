"""This module contains functions to create OperatingPoints for use in tests.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_basic_operating_point_fixture: This method makes a fixture that is an
    OperatingPoint with standard atmospheric conditions for general testing.

    make_zero_alpha_beta_operating_point_fixture: This method makes a fixture
    that is an OperatingPoint with zero alpha and beta for aligned flow testing.

    make_high_alpha_operating_point_fixture: This method makes a fixture that is
    an OperatingPoint with large positive alpha for large angle transformation
    testing.

    make_negative_alpha_operating_point_fixture: This method makes a fixture that
    is an OperatingPoint with negative alpha for negative angle handling testing.

    make_nonzero_beta_operating_point_fixture: This method makes a fixture that is
    an OperatingPoint with non-zero sideslip angle for 3D wind axes testing.

    make_high_speed_operating_point_fixture: This method makes a fixture that is
    an OperatingPoint with high velocity for dynamic pressure scaling testing.

    make_low_density_operating_point_fixture: This method makes a fixture that is
    an OperatingPoint with low density for altitude effects testing.

    make_with_external_force_operating_point_fixture: This method makes a fixture
    that is an OperatingPoint with non-zero external force for trim analysis
    testing.

    make_custom_viscosity_operating_point_fixture: This method makes a fixture
    that is an OperatingPoint with custom kinematic viscosity for vortex core
    growth parameter testing.

    make_boundary_alpha_operating_point_fixture: This method makes a fixture that
    is an OperatingPoint with alpha at boundary values for angle wrapping edge
    case testing.
"""

import pterasoftware as ps


def make_basic_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with standard
    atmospheric conditions for general testing.

    :return basic_operating_point_fixture: OperatingPoint
        This is the OperatingPoint configured for general testing with standard
        sea level conditions, moderate speed, small positive alpha, and zero beta.
    """
    basic_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return basic_operating_point_fixture


def make_zero_alpha_beta_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with zero alpha and
    beta for aligned flow testing.

    :return zero_alpha_beta_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with alpha and beta both set to zero to test
        baseline wind axes alignment.
    """
    zero_alpha_beta_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=0.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return zero_alpha_beta_operating_point_fixture


def make_high_alpha_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with large positive
    alpha for large angle transformation testing.

    :return high_alpha_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with large positive angle of attack to test
        large angle transformations.
    """
    high_alpha_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=45.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return high_alpha_operating_point_fixture


def make_negative_alpha_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with negative alpha
    for negative angle handling testing.

    :return negative_alpha_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with negative angle of attack to test negative
        angle handling.
    """
    negative_alpha_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=-15.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return negative_alpha_operating_point_fixture


def make_nonzero_beta_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with non-zero sideslip
    angle for 3D wind axes testing.

    :return nonzero_beta_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with non-zero sideslip angle to test 3D wind
        axes orientation.
    """
    nonzero_beta_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=10.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return nonzero_beta_operating_point_fixture


def make_high_speed_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with high velocity for
    dynamic pressure scaling testing.

    :return high_speed_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with high velocity to test dynamic pressure
        scaling.
    """
    high_speed_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=100.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return high_speed_operating_point_fixture


def make_low_density_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with low density for
    altitude effects testing.

    :return low_density_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with low air density to test high altitude
        effects on dynamic pressure.
    """
    low_density_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=0.3,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return low_density_operating_point_fixture


def make_with_external_force_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with non-zero external
    force for trim analysis testing.

    :return with_external_force_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with non-zero external force to test trim
        analysis configuration.
    """
    with_external_force_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=50.0,
        nu=15.06e-6,
    )

    return with_external_force_operating_point_fixture


def make_custom_viscosity_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with custom kinematic
    viscosity for vortex core growth parameter testing.

    :return custom_viscosity_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with custom kinematic viscosity value to test
        vortex core growth parameter.
    """
    custom_viscosity_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=20.0e-6,
    )

    return custom_viscosity_operating_point_fixture


def make_boundary_alpha_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with alpha at boundary
    values for angle wrapping edge case testing.

    :return boundary_alpha_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with alpha at the boundary value of 180
        degrees to test angle wrapping edge cases.
    """
    boundary_alpha_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=180.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return boundary_alpha_operating_point_fixture
