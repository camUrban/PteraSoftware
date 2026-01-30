"""This module contains functions to create OperatingPoints for use in tests."""

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


def make_negative_beta_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with negative sideslip
    angle for negative angle handling testing.

    :return negative_beta_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with negative sideslip angle to test negative
        angle handling in wind axes orientation.
    """
    negative_beta_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=-15.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return negative_beta_operating_point_fixture


def make_boundary_beta_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with beta at boundary
    values for angle wrapping edge case testing.

    :return boundary_beta_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with beta at the boundary value of 180
        degrees to test angle wrapping edge cases.
    """
    boundary_beta_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=0.0,
        beta=180.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return boundary_beta_operating_point_fixture


def make_combined_boundary_angles_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with both alpha and
    beta at boundary values for combined boundary edge case testing.

    :return combined_boundary_angles_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with both alpha and beta at boundary
        values to test combined boundary edge cases.
    """
    combined_boundary_angles_operating_point_fixture = (
        ps.operating_point.OperatingPoint(
            rho=1.225,
            vCg__E=10.0,
            alpha=180.0,
            beta=180.0,
            externalFX_W=0.0,
            nu=15.06e-6,
        )
    )

    return combined_boundary_angles_operating_point_fixture


def make_very_low_speed_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with very low but valid
    velocity for testing numerical stability at low speeds.

    :return very_low_speed_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with very low velocity to test numerical
        stability at low dynamic pressures.
    """
    very_low_speed_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=0.01,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return very_low_speed_operating_point_fixture


def make_integer_parameters_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint initialized with integer
    values to verify type conversion to floats.

    :return integer_parameters_operating_point_fixture: OperatingPoint
        This is the OperatingPoint initialized with integer values to test
        internal conversion to floats.
    """
    integer_parameters_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1,
        vCg__E=10,
        alpha=5,
        beta=0,
        externalFX_W=0,
        nu=1,
    )

    return integer_parameters_operating_point_fixture


def make_negative_external_force_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with negative external
    force for drag simulation testing.

    :return negative_external_force_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with negative external force to test drag
        simulation.
    """
    negative_external_force_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=-25.0,
        nu=15.06e-6,
    )

    return negative_external_force_operating_point_fixture


def make_near_boundary_alpha_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with alpha near the
    lower boundary for testing near boundary behavior.

    :return near_boundary_alpha_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with alpha near the lower boundary value
        of negative 180 degrees to test near boundary behavior.
    """
    near_boundary_alpha_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=-179.999,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return near_boundary_alpha_operating_point_fixture


def make_near_boundary_beta_operating_point_fixture():
    """This method makes a fixture that is an OperatingPoint with beta near the
    lower boundary for testing near boundary behavior.

    :return near_boundary_beta_operating_point_fixture: OperatingPoint
        This is the OperatingPoint with beta near the lower boundary value
        of negative 180 degrees to test near boundary behavior.
    """
    near_boundary_beta_operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=0.0,
        beta=-179.999,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    return near_boundary_beta_operating_point_fixture
