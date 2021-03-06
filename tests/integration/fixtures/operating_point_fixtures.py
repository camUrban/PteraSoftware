""" This module creates operating point objects to be used as fixtures.

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
    """ This method makes an operating point object for use in tests.

    :return operating_point_fixture: OperatingPoint
        This is an operating point fixture.
    """

    # Create and return an operating point fixture.
    operating_point_fixture = ps.operating_point.OperatingPoint()
    return operating_point_fixture
