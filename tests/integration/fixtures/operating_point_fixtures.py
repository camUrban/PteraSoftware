"""This module creates OperatingPoints objects to be used as fixtures."""

import pterasoftware as ps


def make_validation_operating_point():
    """This method makes an OperatingPoint for use in tests.

    :return operating_point_fixture: OperatingPoint
        This is an OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225, vCg__E=10.0, alpha=5.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
    )
    return operating_point_fixture


def make_surface_effect_operating_point():
    """This function creates an OperatingPoint with an image surface for surface effect
    testing.

    The ground plane is at the Earth origin with a downward normal (+z is down). The
    first Airplane's CG is at (0, 0, -5), which is 5 meters above the
    ground plane.

    :return operating_point_fixture: OperatingPoint
        This is an OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        surfaceNormal_E=(0.0, 0.0, 1.0),
        surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        CgP1_E_Eo=(0.0, 0.0, -5.0),
    )
    return operating_point_fixture


def make_surface_effect_free_air_operating_point():
    """This function creates an OperatingPoint without an image surface, matching the
    flow conditions of the surface effect operating point, for use as a free-air
    baseline in surface effect validation tests.

    :return operating_point_fixture: OperatingPoint
        This is an OperatingPoint fixture.
    """
    operating_point_fixture = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
    )
    return operating_point_fixture
