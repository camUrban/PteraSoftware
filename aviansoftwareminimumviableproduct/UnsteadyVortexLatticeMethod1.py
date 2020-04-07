from aviansoftwareminimumviableproduct import *
from aviansoftwareminimumviableproduct import mesh_tools
from aviansoftwareminimumviableproduct import output_tools
import numpy as np


class UnsteadyVortexLatticeMethod1(UnsteadyAerodynamicsProblem):

    def __init__(self, airplane, operating_point, movement, simulation_duration, simulation_time_step):
        super().__init__(airplane, operating_point, movement, simulation_duration, simulation_time_step)
        self.verbose = True
        self.front_left_vertices = None
        self.front_right_vertices = None
        self.back_left_vertices = None
        self.back_right_vertices = None
        self.front_wing_vortex_centers = None
        self.back_wing_vortex_centers = None
        self.left_wing_vortex_centers = None
        self.right_wing_vortex_centers = None
        self.front_wing_vortex_legs = None
        self.right_wing_vortex_legs = None
        self.back_wing_vortex_legs = None
        self.left_wing_vortex_legs = None
        self.areas = None
        self.is_trailing_edge = None
        self.is_leading_edge = None
        self.collocation_points = None
        self.normal_directions = None
        self.n_panels = None

        self.time_steps = None
        self.sweep_angles = None
        self.current_time = 0
        self.current_sweep_angle = 0

        self.initial_front_left_vertices = None
        self.initial_front_right_vertices = None
        self.initial_back_left_vertices = None
        self.initial_back_right_vertices = None

    def run(self):

        self.time_steps = np.arange(0, self.simulation_duration + self.simulation_time_step, self.simulation_time_step)
        self.sweep_angles = np.arcsin(np.sin(self.movement.sweeping_amplitude / 2) * np.sin(
            2 * np.pi * self.time_steps / self.movement.movement_period))

        for i in range(len(self.time_steps)):
            self.current_time = self.time_steps[i]
            self.current_sweep_angle = self.sweep_angles[i]
            if self.current_time == 0:
                mesh_tools.make_panels(self)

                self.initial_front_left_vertices = self.front_left_vertices
                self.initial_front_right_vertices = self.front_right_vertices
                self.initial_back_left_vertices = self.back_left_vertices
                self.initial_back_right_vertices = self.back_right_vertices
                output_tools.draw(self)
            else:
                mesh_tools.move_panels(self)
                output_tools.draw(self)
