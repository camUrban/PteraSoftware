
# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp
import tests as tests


# ToDo: Properly document this method.
def make_validation_movement():
    """

    :return:
    """
    
    unsteady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )
    
    unsteady_validation_root_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
        base_wing_cross_section=unsteady_validation_airplane.wings[0].wing_cross_sections[0]
    )

    unsteady_validation_tip_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
        base_wing_cross_section=unsteady_validation_airplane.wings[0].wing_cross_sections[1],
    )

    unsteady_validation_wing_movement = asmvp.movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_sections_movements=[unsteady_validation_root_wing_cross_section_movement,
                                       unsteady_validation_tip_wing_cross_section_movement]
    )

    del unsteady_validation_root_wing_cross_section_movement
    del unsteady_validation_tip_wing_cross_section_movement

    unsteady_validation_airplane_movement = asmvp.movement.AirplaneMovement(
        base_airplane=unsteady_validation_airplane,
        wing_movements=[unsteady_validation_wing_movement]
    )

    del unsteady_validation_airplane
    del unsteady_validation_wing_movement

    unsteady_validation_operating_point_movement = asmvp.movement.OperatingPointMovement(
        base_operating_point=unsteady_validation_operating_point
    )
    
    del unsteady_validation_operating_point

    unsteady_validation_movement = asmvp.movement.Movement(
        airplane_movement=unsteady_validation_airplane_movement,
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=15,
        delta_time=1 / 6 / 10
    )

    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement
    
    return unsteady_validation_movement
