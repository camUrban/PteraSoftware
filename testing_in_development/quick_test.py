import tests.integration.fixtures.airplane_fixtures
import aviansoftwareminimumviableproduct as asmvp


airplane = tests.integration.fixtures.airplane_fixtures.make_steady_validation_airplane()

asmvp.output.draw_geometry(airplane)
