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
