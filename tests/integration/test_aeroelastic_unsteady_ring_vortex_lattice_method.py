"""Integration tests for the AeroelasticUnsteadyRingVortexLatticeMethodSolver.

These tests verify three things:

1. The solver runs to completion and populates the expected output state.
2. A wing with higher density deforms more than a wing with lower density when all
   other parameters are held constant.
3. A wing with a softer torsional spring deforms more than one with a stiffer spring
   when all other parameters are held constant.

The second property follows from the torsional spring-damper model. The prescribed
flapping motion applies an inertial torque tau_inertial = I * d^2(theta_prescribed)/dt^2
that is proportional to the wing's rotational inertia, and therefore to wing_density. A
heavier wing applies a larger inertial torque and so deforms more. Because that torque
scales with mass while the aerodynamic forcing does not, the density comparison is run
with an elevated damping constant that suppresses the aerodynamically driven motion of
the light wing, isolating the mass dependence and making it robustly monotonic.

The third property follows from the same model at low flapping frequency: when the
excitation is slow relative to the torsional natural frequency, the response is
quasi-static and the deformation amplitude is approximately tau / k, so halving the
spring constant roughly doubles the deformation. That comparison uses low damping to
stay in this stiffness-dominated regime rather than the damping-limited one.
"""

import unittest

import numpy as np

import pterasoftware as ps


def _make_aeroelastic_solver(
    wing_density: float,
    spring_constant: float = 20.0,
    damping_constant: float = 1.0,
) -> (
    ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver
):
    """Create a minimal AeroelasticUnsteadyRingVortexLatticeMethodSolver.

    Uses a symmetric two-strip wing with sinusoidal flapping to produce measurable
    aeroelastic deformation. The same geometry and motion are used for every test
    condition; only the structural parameters passed here vary. Different physics
    tests probe different structural regimes, so spring_constant and damping_constant
    are exposed as parameters rather than hard-coded.

    The main wing is declared symmetric so that the reflected Wing is automatically
    generated. AeroelasticUnsteadyProblem.initialize_next_problem requires at least
    two WingMovements (wings[0] and wings[1]) for the deformation update.

    :param wing_density: The wing mass per unit area (kg/m^2).
    :param spring_constant: The torsional spring stiffness (N*m/rad). The default is
        20.0.
    :param damping_constant: The torsional damping constant (N*m*s/rad). The default
        is 1.0.
    :return: A configured AeroelasticUnsteadyRingVortexLatticeMethodSolver ready
        to run.
    """
    airfoil = ps.geometry.airfoil.Airfoil(name="naca2412")

    # Three WingCrossSections: root, one intermediate, and tip. Each non-tip
    # WingCrossSection has num_spanwise_panels=1 so each strip is modeled as a
    # single torsional spring element.
    root_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        num_spanwise_panels=1,
        chord=1.0,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing="uniform",
        airfoil=airfoil,
    )
    mid_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        num_spanwise_panels=1,
        chord=0.8,
        Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing="uniform",
        airfoil=airfoil,
    )
    tip_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        num_spanwise_panels=None,
        chord=0.6,
        Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing=None,
        airfoil=airfoil,
    )

    # Wing root is offset slightly from the symmetry plane (y=0), which produces
    # type-5 (non-coincident) symmetry and causes the Airplane constructor to generate
    # a reflected Wing at airplane.wings[1]. AeroelasticUnsteadyProblem requires at
    # least two WingMovements so that wings[0] and wings[1] both receive deformation.
    main_wing = ps.geometry.wing.Wing(
        wing_cross_sections=[
            root_wing_cross_section,
            mid_wing_cross_section,
            tip_wing_cross_section,
        ],
        name="Main Wing",
        Ler_Gs_Cgs=(0.0, 0.1, 0.0),
        angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        explode_into_strips=False,
        num_chordwise_panels=2,
        chordwise_spacing="uniform",
    )

    airplane = ps.geometry.airplane.Airplane(
        wings=[main_wing],
        name="Test Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )

    # Build WingCrossSectionMovements from wings[0]'s WingCrossSections. Following
    # the example's pattern, the same movement objects are reused for the reflected
    # Wing (wings[1]) since both halves flap symmetrically.
    wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section,
        )
        for wing_cross_section in airplane.wings[0].wing_cross_sections
    ]

    main_wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=wing_cross_section_movements,
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(169.0, 0.0, 0.0),
    )

    reflected_wing_movement = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=airplane.wings[1],
            wing_cross_section_movements=wing_cross_section_movements,
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(169.0, 0.0, 0.0),
        )
    )

    airplane_movement = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=airplane,
            wing_movements=[main_wing_movement, reflected_wing_movement],
        )
    )

    operating_point = ps.operating_point.OperatingPoint(
        rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
    )

    operating_point_movement = ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
        base_operating_point=operating_point,
    )

    movement = ps.movements.aeroelastic_movement.AeroelasticMovement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        delta_time=0.05,
        num_steps=20,
    )

    problem = ps.problems.AeroelasticUnsteadyProblem(
        movement=movement,
        wing_density=wing_density,
        spring_constant=spring_constant,
        damping_constant=damping_constant,
    )

    return ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver(
        aeroelastic_unsteady_problem=problem,
    )


class TestAeroelasticUnsteadySolverCompletion(unittest.TestCase):
    """Verifies that the AeroelasticUnsteadyRingVortexLatticeMethodSolver runs to
    completion and populates the expected output state."""

    def setUp(self) -> None:
        """Create and run the solver.

        :return: None
        """
        self.solver = _make_aeroelastic_solver(wing_density=0.01)
        self.solver.run(
            prescribed_wake=False,
            calculate_streamlines=False,
            show_progress=False,
        )

    def test_solver_completes_and_populates_data(self) -> None:
        """The solver produces an AeroelasticUnsteadyProblem with net_data after a
        successful run.

        net_data is appended once per initialize_next_problem call, which is invoked
        num_steps - 1 times (steps 1 through num_steps - 1).

        :return: None
        """
        self.assertIsInstance(
            self.solver.unsteady_problem, ps.problems.AeroelasticUnsteadyProblem
        )
        problem = self.solver._aeroelastic_unsteady_problem
        # initialize_next_problem is called for steps 1 through num_steps - 1.
        expected_net_data_length = 19  # num_steps - 1 = 20 - 1
        self.assertEqual(len(problem.net_data_per_wing[0]), expected_net_data_length)


class TestAeroelasticUnsteadySolverPhysics(unittest.TestCase):
    """Verifies physically consistent behavior across wing density values.

    A heavier wing accumulates more inertial torque from the prescribed flapping
    motion and therefore deforms more than a lighter wing when all other parameters
    are held constant. The comparison uses an elevated damping constant so the model
    sits in the damping-limited regime: at low density the light wing's motion is
    strongly suppressed by damping, while the heavier wing's larger mass-proportional
    inertial torque overcomes that damping and produces a substantially larger
    torsional deformation. This separates the density signal cleanly from the
    aerodynamic forcing, which is independent of wing mass. At the default low damping
    the density effect is present but small and saturates quickly, so DAMPING_CONSTANT
    is raised here to make the monotonic mass dependence robust.
    """

    LOW_DENSITY = 0.5
    HIGH_DENSITY = 6.0
    DAMPING_CONSTANT = 10.0

    def setUp(self) -> None:
        """Create and run both a low-density and a high-density solver.

        :return: None
        """
        self.low_density_solver = _make_aeroelastic_solver(
            self.LOW_DENSITY, damping_constant=self.DAMPING_CONSTANT
        )
        self.high_density_solver = _make_aeroelastic_solver(
            self.HIGH_DENSITY, damping_constant=self.DAMPING_CONSTANT
        )
        self.low_density_solver.run(
            prescribed_wake=False,
            calculate_streamlines=False,
            show_progress=False,
        )
        self.high_density_solver.run(
            prescribed_wake=False,
            calculate_streamlines=False,
            show_progress=False,
        )

    def test_higher_density_deforms_more(self) -> None:
        """A heavier wing deforms more than a lighter wing under identical conditions.

        net_data entries have shape (num_spanwise_panels + 1, 3). The y-component
        (index 1) holds the torsional angle in radians; the last row corresponds to
        the outermost strip. We compare the peak absolute torsional angle over the
        full simulation so the result is independent of the sign convention at any
        particular step.

        :return: None
        """
        low_density_problem = self.low_density_solver._aeroelastic_unsteady_problem
        high_density_problem = self.high_density_solver._aeroelastic_unsteady_problem

        low_density_outermost_thetas = np.array(
            low_density_problem.net_data_per_wing[0]
        )[:, -1, 1]
        high_density_outermost_thetas = np.array(
            high_density_problem.net_data_per_wing[0]
        )[:, -1, 1]

        max_theta_low_density = float(np.max(np.abs(low_density_outermost_thetas)))
        max_theta_high_density = float(np.max(np.abs(high_density_outermost_thetas)))

        self.assertGreater(max_theta_high_density, max_theta_low_density)


class TestAeroelasticUnsteadySolverSpringStiffness(unittest.TestCase):
    """Verifies physically consistent behavior across torsional spring stiffness.

    At low flapping frequency the torsional response is quasi-static, so the
    deformation amplitude is approximately tau / spring_constant. A softer spring
    therefore yields a larger torsional deformation than a stiffer one when every
    other parameter is held constant. Low damping is used so the response stays in
    this stiffness-dominated regime rather than the damping-limited regime exercised
    by the density test.
    """

    SOFT_SPRING = 5.0
    STIFF_SPRING = 100.0
    DENSITY = 0.05
    DAMPING_CONSTANT = 1.0

    def setUp(self) -> None:
        """Create and run both a soft-spring and a stiff-spring solver.

        :return: None
        """
        self.soft_spring_solver = _make_aeroelastic_solver(
            self.DENSITY,
            spring_constant=self.SOFT_SPRING,
            damping_constant=self.DAMPING_CONSTANT,
        )
        self.stiff_spring_solver = _make_aeroelastic_solver(
            self.DENSITY,
            spring_constant=self.STIFF_SPRING,
            damping_constant=self.DAMPING_CONSTANT,
        )
        self.soft_spring_solver.run(
            prescribed_wake=False,
            calculate_streamlines=False,
            show_progress=False,
        )
        self.stiff_spring_solver.run(
            prescribed_wake=False,
            calculate_streamlines=False,
            show_progress=False,
        )

    def test_softer_spring_deforms_more(self) -> None:
        """A softer torsional spring deforms more than a stiffer one under identical
        conditions.

        net_data entries have shape (num_spanwise_panels + 1, 3). The y-component
        (index 1) holds the torsional angle in radians; the last row corresponds to
        the outermost strip. We compare the peak absolute torsional angle over the
        full simulation so the result is independent of the sign convention at any
        particular step.

        :return: None
        """
        soft_spring_problem = self.soft_spring_solver._aeroelastic_unsteady_problem
        stiff_spring_problem = self.stiff_spring_solver._aeroelastic_unsteady_problem

        soft_spring_outermost_thetas = np.array(
            soft_spring_problem.net_data_per_wing[0]
        )[:, -1, 1]
        stiff_spring_outermost_thetas = np.array(
            stiff_spring_problem.net_data_per_wing[0]
        )[:, -1, 1]

        max_theta_soft_spring = float(np.max(np.abs(soft_spring_outermost_thetas)))
        max_theta_stiff_spring = float(np.max(np.abs(stiff_spring_outermost_thetas)))

        self.assertGreater(max_theta_soft_spring, max_theta_stiff_spring)
