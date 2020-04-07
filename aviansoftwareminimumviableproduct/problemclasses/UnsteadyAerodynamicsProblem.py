from aviansoftwareminimumviableproduct.problemclasses import AerodynamicsProblem


class UnsteadyAerodynamicsProblem(AerodynamicsProblem):

    def __init__(self, airplane, operating_point, movement, simulation_duration, simulation_time_step):
        super().__init__(airplane, operating_point)
        self.movement = movement
        self.simulation_duration = simulation_duration
        self.simulation_time_step = simulation_time_step
