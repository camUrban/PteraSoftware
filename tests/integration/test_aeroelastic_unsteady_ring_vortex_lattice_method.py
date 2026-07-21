"""Integration tests for the AeroelasticUnsteadyRingVortexLatticeMethodSolver.

These tests verify four things:

1. The solver runs to completion and populates the expected output state.
2. A wing with higher density deforms more than a wing with lower density when all
   other parameters are held constant.
3. A wing with a softer torsional spring deforms more than one with a stiffer spring
   when all other parameters are held constant.
4. A geometrically and kinematically symmetric flapping configuration produces a
   mirrored aeroelastic response: the per-strip torsional aerodynamic moments (the y
   components, in the first Airplane's geometry axes, of the moments relative to each
   strip's leading edge point) and the deformation angle y component histories (the
   entries of listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz) of the main and reflected
   wings match.

The second property follows from the torsional spring-damper model. The prescribed
flapping motion applies an inertial moment M_inertial = I * d^2(theta_prescribed)/dt^2
that is proportional to the wing's rotational inertia, and therefore to wing_density. A
heavier wing applies a larger inertial moment and so deforms more. Because that moment
scales with mass while the aerodynamic forcing does not, the density comparison is run
with an elevated damping constant that suppresses the aerodynamically driven motion of
the light wing, isolating the mass dependence and making it robustly monotonic.

The third property follows from the same model at low flapping frequency: when the
excitation is slow relative to the torsional natural frequency, the response is
quasi-static and the deformation amplitude is approximately M_external / k, so halving the
spring constant roughly doubles the deformation. That comparison uses low damping to
stay in this stiffness-dominated regime rather than the damping-limited one.

The fourth property follows from mirror symmetry. The fixture wing is type-5
symmetric, so the Airplane constructor generates the reflected half as a mirror-meshed
Wing whose panel grid runs tip to root spanwise while the structural solve runs root
to tip. Under a symmetric flap, each strip's torsional forcing and torsional
deformation must match its mirror image's once the solve maps between the two
orderings correctly, both in
the strip pairing itself and in the SLEP moment reference points, which must be
mirror-image corners on the two halves. With both mappings correct, the halves match
to machine precision; mismatched chord, width, moment-arm, or reference-point pairings
show up many orders of magnitude above the test's tolerance.
"""

import unittest

import numpy as np

import pterasoftware as ps


def _make_aeroelastic_solver(
    wing_density: float,
    spring_constant_rad: float = 20.0,
    damping_constant_rad: float = 1.0,
) -> (
    ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver
):
    """Create a minimal AeroelasticUnsteadyRingVortexLatticeMethodSolver.

    Uses a symmetric two-strip wing with sinusoidal flapping to produce measurable
    aeroelastic deformation. The same geometry and motion are used for every test
    condition; only the structural parameters passed here vary. Different physics tests
    probe different structural regimes, so spring_constant_rad and damping_constant_rad
    are exposed as parameters rather than hard-coded.

    The main wing is declared symmetric so that the reflected Wing is automatically
    generated. AeroelasticUnsteadyProblem.initialize_next_problem requires at least two
    WingMovements (wings[0] and wings[1]) for the deformation update.

    :param wing_density: The wing mass per unit area (kg/m^2).
    :param spring_constant_rad: The torsional spring stiffness (N*m/rad). The default is
        20.0.
    :param damping_constant_rad: The torsional damping constant (N*m*s/rad). The default
        is 1.0.
    :return: A configured AeroelasticUnsteadyRingVortexLatticeMethodSolver ready to run.
    """
    airfoil = ps.geometry.airfoil.Airfoil(name="naca2412")

    # Three WingCrossSections: root, one intermediate, and tip. Each non-tip
    # WingCrossSection has num_spanwise_panels=1 so each strip is modeled as a single
    # torsional spring element.
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

    # Wing root is offset slightly from the symmetry plane (y=0), which produces type-5
    # (non-coincident) symmetry and causes the Airplane constructor to generate a
    # reflected Wing at airplane.wings[1]. AeroelasticUnsteadyProblem requires at least
    # two WingMovements so that wings[0] and wings[1] both receive deformation.
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

    # Build WingCrossSectionMovements from wings[0]'s WingCrossSections. Following the
    # example's pattern, the same movement objects are reused for the reflected Wing
    # (wings[1]) since both halves flap symmetrically.
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

    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point,
        )
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
        spring_constant_rad=spring_constant_rad,
        damping_constant_rad=damping_constant_rad,
    )

    return ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver(
        aeroelastic_unsteady_problem=problem,
    )


def _per_strip_torsional_aero_moments(
    problem: ps.problems.AeroelasticUnsteadyProblem, wing_idx: int
) -> np.ndarray:
    """Derive each strip's torsional aerodynamic moment from the per-step Panels.

    For every solved time step, re-references each Panel's stored moment (in the first
    Airplane's geometry axes, relative to the first Airplane's CG) to be relative to its
    strip's leading edge point, using the moment transfer relation, then sums the y
    components over each strip's chordwise Panels. Each strip's leading edge point is
    the outboard front point of its first-chordwise-row Panel: the front right panel
    point on a root-to-tip grid, and the front left panel point on a mirror-meshed (tip-
    to-root) grid. The strip axis is returned in root-to-tip order, so mirror-meshed
    grids' columns are reversed.

    :param problem: The solved AeroelasticUnsteadyProblem.
    :param wing_idx: Index of the wing in airplane.wings.
    :return: A (num_steps, num_spanwise_panels) ndarray of torsional moments (N*m).
    """
    per_step_strip_moments = []
    for steady_problem in problem.steady_problems:
        wing = steady_problem.airplanes[0].wings[wing_idx]
        assert wing.panels is not None
        num_spanwise_panels = wing.num_spanwise_panels
        assert num_spanwise_panels is not None
        strip_moments = []
        for span in range(num_spanwise_panels):
            first_row_panel = wing.panels[0][span]
            if wing.mirror_only:
                Slep_GP1_CgP1 = first_row_panel.Flpp_GP1_CgP1
            else:
                Slep_GP1_CgP1 = first_row_panel.Frpp_GP1_CgP1
            strip_moment = 0.0
            for row in wing.panels:
                panel = row[span]
                assert panel.forces_GP1 is not None
                assert panel.moments_GP1_CgP1 is not None
                moment_GP1_Slep = panel.moments_GP1_CgP1 - np.cross(
                    Slep_GP1_CgP1, panel.forces_GP1
                )
                strip_moment += float(moment_GP1_Slep[1])
            strip_moments.append(strip_moment)
        if wing.mirror_only:
            strip_moments.reverse()
        per_step_strip_moments.append(strip_moments)
    return np.array(per_step_strip_moments)


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
        """The solver produces an AeroelasticUnsteadyProblem with a populated
        deformation angle time series (listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz) after
        a successful run.

        The time series is seeded with one initial-state entry at construction.
        initialize_next_problem is invoked on every time step and returns early on the
        final one, so num_steps - 1 of its calls (steps 1 through num_steps - 1) each
        append one entry.

        :return: None
        """
        self.assertIsInstance(
            self.solver.unsteady_problem, ps.problems.AeroelasticUnsteadyProblem
        )
        problem = self.solver._aeroelastic_unsteady_problem
        # One seed entry plus one entry for steps 1 through num_steps - 1.
        expected_time_series_length = 20  # 1 + (num_steps - 1) = 1 + (20 - 1)
        self.assertEqual(
            len(problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]),
            expected_time_series_length,
        )


class TestAeroelasticUnsteadySolverPhysics(unittest.TestCase):
    """Verifies physically consistent behavior across wing density values.

    A heavier wing accumulates more inertial moment from the prescribed flapping motion
    and therefore deforms more than a lighter wing when all other parameters are held
    constant. The comparison uses an elevated damping constant so the model sits in the
    damping-limited regime: at low density the light wing's motion is strongly
    suppressed by damping, while the heavier wing's larger mass-proportional inertial
    moment overcomes that damping and produces a substantially larger torsional
    deformation. This separates the density signal cleanly from the aerodynamic forcing,
    which is independent of wing mass. At the default low damping the density effect is
    present but small and saturates quickly, so DAMPING_CONSTANT is raised here to make
    the monotonic mass dependence robust.
    """

    LOW_DENSITY = 0.5
    HIGH_DENSITY = 6.0
    DAMPING_CONSTANT = 10.0

    def setUp(self) -> None:
        """Create and run both a low-density and a high-density solver.

        :return: None
        """
        self.low_density_solver = _make_aeroelastic_solver(
            self.LOW_DENSITY, damping_constant_rad=self.DAMPING_CONSTANT
        )
        self.high_density_solver = _make_aeroelastic_solver(
            self.HIGH_DENSITY, damping_constant_rad=self.DAMPING_CONSTANT
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

        listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz entries have shape
        (num_spanwise_panels + 1,) and hold the torsional angles in radians; the last
        element corresponds to the outermost strip. We compare the peak absolute
        torsional angle over the full simulation so the result is independent of the
        sign convention at any particular step.

        :return: None
        """
        low_density_problem = self.low_density_solver._aeroelastic_unsteady_problem
        high_density_problem = self.high_density_solver._aeroelastic_unsteady_problem

        low_density_outermost_thetas = np.array(
            low_density_problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]
        )[:, -1]
        high_density_outermost_thetas = np.array(
            high_density_problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]
        )[:, -1]

        max_theta_low_density = float(np.max(np.abs(low_density_outermost_thetas)))
        max_theta_high_density = float(np.max(np.abs(high_density_outermost_thetas)))

        self.assertGreater(max_theta_high_density, max_theta_low_density)


class TestAeroelasticUnsteadySolverSpringStiffness(unittest.TestCase):
    """Verifies physically consistent behavior across torsional spring stiffness.

    At low flapping frequency the torsional response is quasi-static, so the deformation
    amplitude is approximately M_external / spring_constant_rad. A softer spring
    therefore yields a larger torsional deformation than a stiffer one when every other
    parameter is held constant. Low damping is used so the response stays in this
    stiffness-dominated regime rather than the damping-limited regime exercised by the
    density test.
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
            spring_constant_rad=self.SOFT_SPRING,
            damping_constant_rad=self.DAMPING_CONSTANT,
        )
        self.stiff_spring_solver = _make_aeroelastic_solver(
            self.DENSITY,
            spring_constant_rad=self.STIFF_SPRING,
            damping_constant_rad=self.DAMPING_CONSTANT,
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

        listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz entries have shape
        (num_spanwise_panels + 1,) and hold the torsional angles in radians; the last
        element corresponds to the outermost strip. We compare the peak absolute
        torsional angle over the full simulation so the result is independent of the
        sign convention at any particular step.

        :return: None
        """
        soft_spring_problem = self.soft_spring_solver._aeroelastic_unsteady_problem
        stiff_spring_problem = self.stiff_spring_solver._aeroelastic_unsteady_problem

        soft_spring_outermost_thetas = np.array(
            soft_spring_problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]
        )[:, -1]
        stiff_spring_outermost_thetas = np.array(
            stiff_spring_problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]
        )[:, -1]

        max_theta_soft_spring = float(np.max(np.abs(soft_spring_outermost_thetas)))
        max_theta_stiff_spring = float(np.max(np.abs(stiff_spring_outermost_thetas)))

        self.assertGreater(max_theta_soft_spring, max_theta_stiff_spring)


class TestAeroelasticUnsteadySolverMirrorSymmetry(unittest.TestCase):
    """Verifies that a geometrically and kinematically symmetric flapping configuration
    produces a mirrored aeroelastic response between the main Wing and the reflected
    Wing the Airplane constructor generates from it.

    The reflected Wing is mirror-meshed, so its panel grid runs tip to root spanwise
    while the structural solve runs root to tip. A correct mapping between the two
    orderings, in both the strip pairing and the SLEP moment reference points, makes the
    two halves' per-strip torsional moments and deformation angle histories match to
    machine precision, so the tolerance is set six orders of magnitude above that floor
    while remaining six or more below any mispairing's signature (mispaired chords,
    strip widths, moment arms, or reference corners produce relative asymmetries between
    about one percent and order one).
    """

    RELATIVE_TOLERANCE = 1e-9
    DENSITY = 0.5
    SPRING_CONSTANT = 500.0
    DAMPING_CONSTANT = 10.0

    def setUp(self) -> None:
        """Create and run the symmetric flapping solver.

        :return: None
        """
        self.solver = _make_aeroelastic_solver(
            self.DENSITY,
            spring_constant_rad=self.SPRING_CONSTANT,
            damping_constant_rad=self.DAMPING_CONSTANT,
        )
        self.solver.run(
            prescribed_wake=False,
            calculate_streamlines=False,
            show_progress=False,
        )
        self.problem = self.solver._aeroelastic_unsteady_problem

    def test_per_strip_torsional_moments_mirror(self) -> None:
        """The per-strip torsional aerodynamic moments of the two halves match.

        The moments (in the first Airplane's geometry axes) are derived from each solved
        step's retained Panel loads, re-referenced to each strip's leading edge point
        and summed chordwise into per-strip torsional forcing in root-to-tip order for
        both halves, so the arrays must match strip for strip.

        :return: None
        """
        aero_main = _per_strip_torsional_aero_moments(self.problem, 0)
        aero_reflected = _per_strip_torsional_aero_moments(self.problem, 1)

        scale = float(np.max(np.abs(aero_main)))
        self.assertGreater(scale, 0.0)
        max_mismatch = float(np.max(np.abs(aero_main - aero_reflected)))
        self.assertLess(max_mismatch, self.RELATIVE_TOLERANCE * scale)

    def test_deformation_angle_histories_mirror(self) -> None:
        """The two halves' deformation angle y component histories
        (listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz) match at every recorded step.

        listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz entries have shape
        (num_spanwise_panels + 1,) in root-to-tip order for both halves, holding the
        torsional angles in radians, so the time series must match WingCrossSection for
        WingCrossSection.

        :return: None
        """
        theta_main = np.array(
            self.problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]
        )
        theta_reflected = np.array(
            self.problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[1]
        )

        peak = float(np.max(np.abs(theta_main)))
        self.assertGreater(peak, 0.0)
        max_mismatch = float(np.max(np.abs(theta_main - theta_reflected)))
        self.assertLess(max_mismatch, self.RELATIVE_TOLERANCE * peak)
