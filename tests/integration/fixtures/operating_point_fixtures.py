# NOTE: I haven't yet started refactoring this module.
"""This module creates operating point objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_validation_operating_point: This method makes an operating point object for
    use in tests.
"""

import pterasoftware as ps


def make_validation_operating_point():
    """This method makes an operating point object for use in tests.

    :return operating_point_fixture: OperatingPoint
        This is an operating point fixture.
    """

    # Create and return an operating point fixture.
    operating_point_fixture = ps.operating_point.OperatingPoint(rho=1.225, vCg__E=10.0,
                                                                alpha=5.0, beta=0.0,
                                                                externalFX_W=0.0,
                                                                nu=15.06e-6)
    return operating_point_fixture
