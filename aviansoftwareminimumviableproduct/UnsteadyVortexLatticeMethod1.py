from aviansoftwareminimumviableproduct import *


class UnsteadyVortexLatticeMethod1(UnsteadyAerodynamicsProblem):

    def __init__(self, airplane, op_point, movement, simulation_duration, simulation_time_step):
        super().__init__(airplane, op_point, movement, simulation_duration, simulation_time_step)
