
# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp
import tests as tests


# ToDo: Properly document this method.
def make_steady_validation_problem():
    """

    :return:
    """

    steady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_steady_validation_airplane()
    )
    steady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_steady_validation_operating_point()
    )

    steady_validation_problem = asmvp.problems.SteadyProblem(
        airplane=steady_validation_airplane, operating_point=steady_validation_operating_point
    )

    del steady_validation_airplane
    del steady_validation_operating_point

    return steady_validation_problem
