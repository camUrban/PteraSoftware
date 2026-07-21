"""This module creates OperatingPoints objects to be used as fixtures."""

import pterasoftware as ps


def make_validation_operating_point() -> ps.operating_point.OperatingPoint:
    """This method makes an OperatingPoint for use in tests.

    :return operating_point_fixture: OperatingPoint This is an OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225, vCg__E=10.0, alpha=5.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
    )
    return operating_point_fixture


def make_surface_effect_operating_point() -> ps.operating_point.OperatingPoint:
    """This function creates an OperatingPoint with an image surface for surface effect
    testing.

    The ground plane is at the Earth origin with a downward normal (+z is down). The
    first Airplane's CG is at (0, 0, -5), which is 5 meters above the ground plane. The
    body attitude pitches the airplane nose up by the 5 degree angle of attack (with
    zero sideslip), which places the flight path level along the horizontal Earth x
    axis, so the airplane flies level in ground effect rather than descending toward the
    surface.

    :return operating_point_fixture: OperatingPoint This is an OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        angles_E_to_BP1_izyx=(0.0, 5.0, 0.0),
        surfaceNormal_E=(0.0, 0.0, 1.0),
        surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        CgP1_E_Eo=(0.0, 0.0, -5.0),
    )
    return operating_point_fixture


def make_surface_effect_free_air_operating_point() -> ps.operating_point.OperatingPoint:
    """This function creates an OperatingPoint without an image surface, matching the
    flow conditions of the surface effect operating point, for use as a free-air
    baseline in surface effect validation tests.

    :return operating_point_fixture: OperatingPoint This is an OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
    )
    return operating_point_fixture


def make_simple_glider_operating_point() -> ps.operating_point.OperatingPoint:
    """This function creates the trimmed OperatingPoint for the simple glider's free
    flight test.

    The speed and angle of attack are the trimmed gliding condition found for the simple
    glider. The initial body orientation pitches the airplane nose up by the angle of
    attack (with zero sideslip), which places the trim velocity along the horizontal
    Earth x axis at the start of free flight. No external thrust is applied, so the
    glider flies an unpowered glide.

    :return operating_point_fixture: OperatingPoint This is the simple glider
        OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=12.9,
        alpha=3.3,
        beta=0.0,
        angles_E_to_BP1_izyx=(0.0, 3.3, 0.0),
        g_E=(0.0, 0.0, 9.80665),
    )
    return operating_point_fixture


def make_flapping_free_flight_operating_point() -> ps.operating_point.OperatingPoint:
    """This function creates the initial OperatingPoint for the flapping-wing free
    flight test.

    The initial body orientation pitches the airplane nose up by the angle of attack
    (with zero sideslip), which places the initial velocity along the horizontal Earth x
    axis at the start of free flight. No external thrust is applied; in free flight the
    propulsion comes from the flapping motion rather than an external force.

    :return operating_point_fixture: OperatingPoint This is the flapping-wing
        OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=12.9,
        alpha=3.3,
        beta=0.0,
        angles_E_to_BP1_izyx=(0.0, 3.3, 0.0),
        g_E=(0.0, 0.0, 9.80665),
    )
    return operating_point_fixture
