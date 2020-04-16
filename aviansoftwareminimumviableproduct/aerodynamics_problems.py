# ToDo: Properly document this module.
"""

"""


# ToDo: Properly document this class.
class AerodynamicsProblem:
    """

    """

    # ToDo: Properly document this method.
    def __init__(self, airplane, operating_point):
        """

        :param airplane:
        :param operating_point:
        """
        self.airplane = airplane
        self.operating_point = operating_point


# ToDo: Properly document this class.
class SteadyAerodynamicsProblem(AerodynamicsProblem):
    """

    """

    # ToDo: Properly document this method.
    def __init__(self, airplane, operating_point):
        """

        :param airplane:
        :param operating_point:
        """
        super().__init__(airplane=airplane, operating_point=operating_point)


# ToDo: Properly document this class.
class MeshedSteadyAerodynamicsProblem(SteadyAerodynamicsProblem):
    """

    """

    # ToDo: Properly document this method.
    def __init__(self, meshed_airplane, operating_point):
        """

        :param meshed_airplane:
        :param operating_point:
        """
        super().__init__(airplane=meshed_airplane, operating_point=operating_point)


# ToDo: Properly document this class.
class UnsteadyAerodynamicsProblem(AerodynamicsProblem):
    """

    """

    # ToDo: Properly document this method.
    def __init__(self, airplane, operating_point, movement, simulation_duration, simulation_time_step):
        """

        :param airplane:
        :param operating_point:
        :param movement:
        :param simulation_duration:
        :param simulation_time_step:
        """
        super().__init__(airplane=airplane, operating_point=operating_point)
        self.movement = movement
        self.simulation_duration = simulation_duration
        self.simulation_time_step = simulation_time_step


# ToDo: Properly document this class.
class MeshedUnsteadyAerodynamicsProblem(UnsteadyAerodynamicsProblem):
    """

    """

    # ToDo: Properly document this method.
    def __init__(self, meshed_airplane, operating_point, movement, simulation_duration, simulation_time_step):
        """

        :param meshed_airplane:
        :param operating_point:
        :param movement:
        :param simulation_duration:
        :param simulation_time_step:
        """
        super().__init__(airplane=meshed_airplane, operating_point=operating_point, movement=movement,
                         simulation_duration=simulation_duration, simulation_time_step=simulation_time_step)
