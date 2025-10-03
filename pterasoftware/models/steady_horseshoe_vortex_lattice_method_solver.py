# NOTE: I haven't yet started refactoring this module.
import pterasoftware as ps


class steadyHorseshoeVortexLatticeMethodSolver:
    def __init__(self):
        var = "Variables"

    def runSolver(self):
        example_airplane = ps.geometry.airplane.Airplane(
            name="Example Airplane",
            x_ref=0.0,
            y_ref=0.0,
            z_ref=0.0,
            s_ref=None,
            b_ref=None,
            c_ref=None,
            wings=[
                ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        ps.geometry.wing_cross_section.WingCrossSection(
                            x_le=0.0,
                            y_le=0.0,
                            z_le=0.0,
                            twist=0.0,
                            control_surface_type="symmetric",
                            control_surface_hinge_point=0.75,
                            control_surface_deflection=0.0,
                            num_spanwise_panels=8,
                            spanwise_spacing="cosine",
                            chord=1.75,
                            airfoil=ps.geometry.airfoil.Airfoil(
                                name="naca2412",
                                coordinates=None,
                                repanel=True,
                                n_points_per_side=400,
                            ),
                        ),
                        ps.geometry.wing_cross_section.WingCrossSection(
                            x_le=0.75,
                            y_le=6.0,
                            z_le=1.0,
                            chord=1.5,
                            twist=5.0,
                            airfoil=ps.geometry.airfoil.Airfoil(
                                name="naca2412",
                            ),
                        ),
                    ],
                    name="Main Wing",
                    symmetric=True,
                    num_chordwise_panels=8,
                    chordwise_spacing="cosine",
                ),
                ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        ps.geometry.wing_cross_section.WingCrossSection(
                            chord=1.5,
                            airfoil=ps.geometry.airfoil.Airfoil(
                                name="naca0012",
                            ),
                            twist=-5.0,
                        ),
                        ps.geometry.wing_cross_section.WingCrossSection(
                            x_le=0.5,
                            y_le=2.0,
                            z_le=1.0,
                            chord=1.0,
                            twist=-5.0,
                            airfoil=ps.geometry.airfoil.Airfoil(
                                name="naca0012",
                            ),
                        ),
                    ],
                    name="V-Tail",
                    symmetric=True,
                ),
            ],
        )

        example_operating_point = ps.operating_point.OperatingPoint(
            rho=1.225, vCg__E=10.0, alpha=1.0, beta=0.0
        )

        example_problem = ps.problems.SteadyProblem(
            airplanes=[example_airplane],
            operating_point=example_operating_point,
        )

        del example_airplane
        del example_operating_point

        example_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
            steady_problem=example_problem
        )

        del example_problem

        example_solver.run(
            logging_level="Warning",
        )

        ps.output.print_results()

        ps.output.draw(
            solver=example_solver,
            scalar_type="lift",
            show_streamlines=True,
            show_wake_vortices=False,
            save=False,
        )
