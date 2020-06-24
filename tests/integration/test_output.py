# ToDo: Properly document this module.

import unittest

import matplotlib.pyplot as plt

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

        self.unsteady_solver = (
            tests.integration.fixtures.solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return:
        """

        del self.unsteady_solver

    # ToDo: Properly document this method.
    def test_plot_results_versus_time(self):
        """

        :return:
        """
        num_figs_before = plt.gcf().number

        asmvp.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver, testing=True
        )

        num_figs_after = plt.gcf().number

        self.assertEqual(num_figs_before + 4, num_figs_after)

    # ToDo: Properly document this method.
    def test_animate_does_not_throw(self):
        """

        :return:
        """

        asmvp.output.animate(
            unsteady_solver=self.unsteady_solver,
            show_delta_pressures=False,
            show_wake_vortices=False,
        )

    # ToDo: Properly document this method.
    def test_draw_does_not_throw(self):
        """

        :return:
        """

        asmvp.output.draw(
            airplane=self.unsteady_solver.steady_problems[0].airplane,
            show_delta_pressures=False,
            show_wake_vortices=False,
            show_streamlines=False,
        )
