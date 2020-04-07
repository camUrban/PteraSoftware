from aviansoftwareminimumviableproduct import *
from aviansoftwareminimumviableproduct import mesh_tools
from aviansoftwareminimumviableproduct import output_tools


class SteadyVortexLatticeMethod1(SteadyAerodynamicsProblem):

    def __init__(self, airplane, operating_point):
        super().__init__(airplane, operating_point)
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

        self.initial_front_left_vertices = None
        self.initial_front_right_vertices = None
        self.initial_back_left_vertices = None
        self.initial_back_right_vertices = None

    def run(self):
        mesh_tools.make_panels(self)
        output_tools.draw(self)
