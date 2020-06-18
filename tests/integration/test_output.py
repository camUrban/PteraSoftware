
# ToDo: Properly document this module.

import unittest
import matplotlib.pyplot as plt
import unittest.mock as mock
import pyvista as pv

import numpy as np

import aviansoftwareminimumviableproduct as asmvp
import tests.integration


# ToDo: Properly document this class.
class TestUnsteadyRingVortexLatticeMethod(unittest.TestCase):
    """

    """

    # ToDo: Properly document this method.
    def setUp(self):
        """

        :return:
        """

        self.movement = (
            tests.integration.fixtures.movement_fixtures.make_static_validation_movement()
        )

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return:
        """

        del self.movement

    # ToDo: Properly document this method.
    def test_plot_results_versus_time(self):
        """

        :return:
        """
        num_figs_before = plt.gcf().number

        asmvp.output.plot_results_versus_time(movement=self.movement, verbose=False)

        num_figs_after = plt.gcf().number

        self.assertEqual(num_figs_before + 4, num_figs_after)
