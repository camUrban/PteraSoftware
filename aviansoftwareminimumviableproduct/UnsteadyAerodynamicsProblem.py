class UnsteadyAerodynamicsProblem:

    def __init__(self, airplane, operating_point, movement, simulation_duration, simulation_time_step):
        self.airplane = airplane,
        self.op_point = operating_point,
        self.movement = movement
        self.simulation_duration = simulation_duration
        self.simulation_time_step = simulation_time_step
