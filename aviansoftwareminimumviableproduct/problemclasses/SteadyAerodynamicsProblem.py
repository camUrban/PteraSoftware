from aviansoftwareminimumviableproduct.problemclasses import AerodynamicsProblem


class SteadyAerodynamicsProblem(AerodynamicsProblem):

    def __init__(self, airplane, operating_point):
        super().__init__(airplane, operating_point)
