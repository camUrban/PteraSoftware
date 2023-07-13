import math

import pterasoftware as ps
from tests.integration.fixtures import problem_fixtures

steady_validation_problem = problem_fixtures.make_steady_validation_problem()

converged_parameters = ps.convergence.analyze_steady_convergence(
    ref_problem=steady_validation_problem,
    solver_type="steady horseshoe vortex lattice method",
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 20),
    convergence_criteria=0.1,
)

[panel_aspect_ratio, num_chordwise_panels] = converged_parameters
section_length = 5
section_standard_mean_chord = 0.875

num_spanwise_panels = round(
    (section_length * num_chordwise_panels)
    / (section_standard_mean_chord * panel_aspect_ratio)
)

num_spanwise_panels = math.ceil(num_spanwise_panels)
print(num_spanwise_panels)
