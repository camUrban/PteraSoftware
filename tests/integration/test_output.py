
# ToDo: Properly document this module.

import unittest
import matplotlib.pyplot as plt

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
            tests.integration.fixtures.movement_fixtures.make_validation_movement()
        )

        num_airplanes = len(self.movement.airplanes)
        times = np.linspace(0, 1, num_airplanes)

        for airplane_num in range(num_airplanes):
            airplane = self.movement.airplanes[airplane_num]
            time = times[airplane_num]

            airplane.total_near_field_force_wind_axes = np.ones(3) * np.sin(time)
            airplane.total_near_field_force_coefficients_wind_axes = - np.ones(3) * np.sin(time)
            airplane.total_near_field_moment_wind_axes = np.ones(3) * np.cos(time)
            airplane.total_near_field_moment_coefficients_wind_axes = - np.ones(3) * np.cos(time)

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return:
        """

        del self.movement

    # ToDo: Properly document this method.
    def test_method(self):
        """

        :return:
        """
        num_figs_before = plt.gcf().number

        asmvp.output.plot_results_versus_time(movement=self.movement, verbose=False)

        num_figs_after = plt.gcf().number

        self.assertEqual(num_figs_before + 4, num_figs_after)
