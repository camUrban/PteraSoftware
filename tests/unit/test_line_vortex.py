
# ToDo: Properly document this module.
"""

"""

import unittest
import numpy as np

import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly document this class.
class TestLineVortex(unittest.TestCase):
    """

    """

    # ToDo: Properly document this method.
    def test_class(self):
        """

        :return:
        """

        origin_fixture = np.zeros(3)
        termination_fixture = np.ones(3)
        strength_fixture = 10
        line_vortex_fixture = asmvp.aerodynamics.LineVortex(
            origin=origin_fixture,
            termination=termination_fixture,
            strength=strength_fixture
        )

        self.assertIsInstance(line_vortex_fixture, asmvp.aerodynamics.LineVortex)
        self.assertTrue(np.allclose(line_vortex_fixture.origin, origin_fixture))
        self.assertTrue(np.allclose(line_vortex_fixture.termination, termination_fixture))
        self.assertTrue(np.allclose(line_vortex_fixture.strength, strength_fixture))
        self.assertTrue(np.allclose(line_vortex_fixture.center, np.array([0.5, 0.5, 0.5])))
        self.assertTrue(np.allclose(line_vortex_fixture.vector, np.array([1, 1, 1])))
