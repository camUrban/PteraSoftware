"""Contains the SteadyProblem, UnsteadyProblem, AeroelasticUnsteadyProblem, and
FreeFlightUnsteadyProblem classes.

**Contains the following classes:**

SteadyProblem: A class used to contain steady aerodynamics problems.

UnsteadyProblem: A class used to contain unsteady aerodynamics problems.

AeroelasticUnsteadyProblem: A class used to couple unsteady aerodynamics with wing
structural dynamics (torsional spring-mass-damper model) for aeroelastic simulations.

FreeFlightUnsteadyProblem: A class used to contain problems with coupled unsteady
aerodynamics and rigid body dynamics.

**Contains the following functions:**

None
"""

from __future__ import annotations

import copy
from collections.abc import Callable, Sequence
from pathlib import PureWindowsPath
from typing import TYPE_CHECKING, cast

import numpy as np
from scipy.integrate import solve_ivp

from . import (
    _core,
    _fixed_point_relaxation,
    _logging,
    _mujoco_model,
    _parameter_validation,
    _private_access,
    _transformations,
    geometry,
    movements,
)
from . import operating_point as operating_point_mod
from .movements import aeroelastic_movement as aeroelastic_movement_mod
from .movements import aeroelastic_wing_movement as aeroelastic_wing_movement_mod
from .movements import free_flight_movement

if TYPE_CHECKING:
    from ._coupled_unsteady_ring_vortex_lattice_method import (
        CoupledUnsteadyRingVortexLatticeMethodSolver,
    )
    from .aeroelastic_unsteady_ring_vortex_lattice_method import (
        AeroelasticUnsteadyRingVortexLatticeMethodSolver,
    )
    from .free_flight_unsteady_ring_vortex_lattice_method import (
        FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    )

_logger = _logging.get_logger("problems")


class SteadyProblem:
    """A class used to contain steady aerodynamics problems.

    **Contains the following methods:**

    reynolds_numbers: A tuple of Reynolds numbers, one for each Airplane in the
    SteadyProblem.
    """

    __slots__ = (
        "_airplanes",
        "_operating_point",
        "_reynolds_numbers",
    )

    def __init__(
        self,
        airplanes: list[geometry.airplane.Airplane],
        operating_point: operating_point_mod.OperatingPoint,
    ) -> None:
        """The initialization method.

        :param airplanes: The list of the Airplanes for this SteadyProblem.
        :param operating_point: The OperatingPoint for this SteadyProblem.
        :return: None
        """
        # Validate and store immutable attributes.
        if not isinstance(airplanes, list):
            raise TypeError("airplanes must be a list.")
        if len(airplanes) < 1:
            raise ValueError("airplanes must have at least one element.")
        for airplane in airplanes:
            if not isinstance(airplane, geometry.airplane.Airplane):
                raise TypeError("Every element in airplanes must be an Airplane.")
        # Store as tuple to prevent external mutation via .append(), .pop(), etc.
        self._airplanes: tuple[geometry.airplane.Airplane, ...] = tuple(airplanes)

        if not isinstance(operating_point, operating_point_mod.OperatingPoint):
            raise TypeError("operating_point must be an OperatingPoint.")
        self._operating_point = operating_point

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._reynolds_numbers: tuple[float, ...] | None = None

        # Validate that the first Airplane has Cg_GP1_CgP1 set to zeros.
        self._airplanes[0].validate_first_airplane_constraints()

        # Populate GP1_CgP1 coordinates for all Airplanes' Panels. This finds the
        # Panels' positions in the first Airplane's geometry axes, relative to the first
        # Airplane's CG based on their locally defined positions.
        for airplane in self._airplanes:
            # Compute the passive transformation matrix from this Airplane's local
            # geometry axes, relative to its CG, to the first Airplane's geometry axes,
            # relative to the first Airplane's CG.
            T_pas_G_Cg_to_GP1_CgP1 = airplane.T_pas_G_Cg_to_GP1_CgP1

            for wing in airplane.wings:
                assert wing.panels is not None

                for panel in np.ravel(wing.panels):
                    panel.Frpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Frpp_G_Cg, is_position=True
                    )
                    panel.Flpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Flpp_G_Cg, is_position=True
                    )
                    panel.Blpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Blpp_G_Cg, is_position=True
                    )
                    panel.Brpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                        T_pas_G_Cg_to_GP1_CgP1, panel.Brpp_G_Cg, is_position=True
                    )

    # --- Immutable: read only properties ---
    @property
    def airplanes(self) -> tuple[geometry.airplane.Airplane, ...]:
        return self._airplanes

    @property
    def operating_point(self) -> operating_point_mod.OperatingPoint:
        return self._operating_point

    # --- Immutable derived: manual lazy caching ---
    @property
    def reynolds_numbers(self) -> tuple[float, ...]:
        """A tuple of Reynolds numbers, one for each Airplane in the SteadyProblem.

        **Notes:**

        The Reynolds number is calculated as: Re = (V x L) / nu, where V is the
        freestream speed, observed from the Earth frame (vCg__E from OperatingPoint,
        m/s), L is the characteristic length (c_ref from Airplane, m), and nu is the
        kinematic viscosity (nu from OperatingPoint, m^2/s).

        These Reynolds numbers only consider the freestream speed, not any apparent
        velocity due to prescribed motion, so be careful interpreting it for cases where
        this SteadyProblem corresponds to one time step in an UnsteadyProblem.

        :return: A tuple of Reynolds numbers, one for each Airplane.
        """
        if self._reynolds_numbers is None:
            v = self._operating_point.vCg__E
            nu = self._operating_point.nu

            reynolds_list = []
            for airplane in self._airplanes:
                c_ref = airplane.c_ref
                assert c_ref is not None, "Airplane c_ref must be set to calculate Re"
                re = (v * c_ref) / nu
                reynolds_list.append(re)

            # Store as tuple to prevent external mutation.
            self._reynolds_numbers = tuple(reynolds_list)
        return self._reynolds_numbers


class UnsteadyProblem(_core.CoreUnsteadyProblem):
    """A class used to contain unsteady aerodynamics problems.

    **Contains the following methods:**

    only_final_results: Determines whether the solver will only calculate loads for the
    final time step or final cycle.

    num_steps: The number of time steps.

    delta_time: The time step size in seconds.

    first_averaging_step: The first time step included in cycle averaging.

    first_results_step: The first time step for which loads are calculated.

    max_wake_rows: The maximum chordwise wake rows per Wing.

    movement: The Movement that contains this UnsteadyProblem's OperatingPointMovement
    and AirplaneMovements.

    steady_problems: A tuple of SteadyProblems, one for each time step.

    finalInducedDrags_W: The final induced drag force experienced by each Airplane (in
    wind axes).

    finalSideForces_W: The final side force experienced by each Airplane (in wind axes).

    finalLifts_W: The final lift force experienced by each Airplane (in wind axes).

    finalInducedDragCoefficients_W: The final induced drag force coefficient experienced
    by each Airplane (in wind axes).

    finalSideForceCoefficients_W: The final side force coefficient experienced by each
    Airplane (in wind axes).

    finalLiftCoefficients_W: The final lift force coefficient experienced by each
    Airplane (in wind axes).

    finalRollingMoments_W_Cg: The final rolling moment experienced by each Airplane (in
    wind axes, relative to its own CG).

    finalPitchingMoments_W_Cg: The final pitching moment experienced by each Airplane
    (in wind axes, relative to its own CG).

    finalYawingMoments_W_Cg: The final yawing moment experienced by each Airplane (in
    wind axes, relative to its own CG).

    finalRollingMomentCoefficients_W_Cg: The final rolling moment coefficient
    experienced by each Airplane (in wind axes, relative to its own CG).

    finalPitchingMomentCoefficients_W_Cg: The final pitching moment coefficient
    experienced by each Airplane (in wind axes, relative to its own CG).

    finalYawingMomentCoefficients_W_Cg: The final yawing moment coefficient experienced
    by each Airplane (in wind axes, relative to its own CG).

    finalMeanInducedDrags_W: The final cycle averaged induced drag force experienced by
    each Airplane (in wind axes).

    finalMeanSideForces_W: The final cycle averaged side force experienced by each
    Airplane (in wind axes).

    finalMeanLifts_W: The final cycle averaged lift force experienced by each Airplane
    (in wind axes).

    finalMeanInducedDragCoefficients_W: The final cycle averaged induced drag force
    coefficient experienced by each Airplane (in wind axes).

    finalMeanSideForceCoefficients_W: The final cycle averaged side force coefficient
    experienced by each Airplane (in wind axes).

    finalMeanLiftCoefficients_W: The final cycle averaged lift force coefficient
    experienced by each Airplane (in wind axes).

    finalMeanRollingMoments_W_Cg: The final cycle averaged rolling moment experienced by
    each Airplane (in wind axes, relative to its own CG).

    finalMeanPitchingMoments_W_Cg: The final cycle averaged pitching moment experienced
    by each Airplane (in wind axes, relative to its own CG).

    finalMeanYawingMoments_W_Cg: The final cycle averaged yawing moment experienced by
    each Airplane (in wind axes, relative to its own CG).

    finalMeanRollingMomentCoefficients_W_Cg: The final cycle averaged rolling moment
    coefficient experienced by each Airplane (in wind axes, relative to its own CG).

    finalMeanPitchingMomentCoefficients_W_Cg: The final cycle averaged pitching moment
    coefficient experienced by each Airplane (in wind axes, relative to its own CG).

    finalMeanYawingMomentCoefficients_W_Cg: The final cycle averaged yawing moment
    coefficient experienced by each Airplane (in wind axes, relative to its own CG).

    **Notes:**

    The solver populates the mutable load lists during simulation, with one entry per
    Airplane: the final lists for static motion, the final cycle averaged lists for
    cyclic motion, and the final cycle root mean squared lists for variable geometry
    problems. The lists with names ending in Forces_W, ForceCoefficients_W,
    Moments_W_CgP1, and MomentCoefficients_W_CgP1 hold the total force (in wind axes),
    the total moment (in wind axes, relative to the first Airplane's CG), and their
    coefficients. The lists with names ending in Forces_G, ForceCoefficients_G,
    Moments_G_Cg, MomentCoefficients_G_Cg, Moments_W_Cg, and MomentCoefficients_W_Cg
    hold the total force (in geometry axes), the total moment (in geometry axes,
    relative to the CG), the total moment (in wind axes, relative to the CG), and their
    coefficients. In these per-Airplane lists, the G and Cg IDs are used without an
    Airplane index to implicitly mean "in the entry's own Airplane's geometry axes" and
    "relative to the entry's own Airplane's CG". For the first Airplane's entries, the
    Moments_W_Cg and MomentCoefficients_W_Cg lists equal the Moments_W_CgP1 and
    MomentCoefficients_W_CgP1 lists.

    The named load properties derive the named load components and coefficients, as
    defined in AXES_POINTS_AND_FRAMES.md, from the wind axes lists. Each returns one
    float per Airplane, or an empty list until its source list has been populated.
    """

    __slots__ = (
        "_movement",
        "_steady_problems",
    )

    def __init__(
        self,
        movement: movements.movement.Movement,
        only_final_results: bool | np.bool_ = False,
    ) -> None:
        """The initialization method.

        :param movement: The Movement that contains this UnsteadyProblem's
            OperatingPointMovement and AirplaneMovements.
        :param only_final_results: Determines whether the Solver will only calculate
            loads for the final time step (for static Movements) or (for non static
            Movements) for will only calculate loads for the time steps in the final
            complete motion cycle (of the Movement's sub Movement with the longest
            period), which increases simulation speed. Can be a bool or a numpy bool and
            will be converted internally to a bool. The default is False.
        :return: None
        """
        # Validate and store the Movement before calling super().__init__() because the
        # Movement provides the parameters that the core class needs.
        if not isinstance(movement, movements.movement.Movement):
            raise TypeError("movement must be a Movement.")
        self._movement = movement

        # Delegate shared initialization (validation, first_averaging_step computation,
        # load list initialization) to the core class.
        super().__init__(
            only_final_results=only_final_results,
            delta_time=self._movement.delta_time,
            num_steps=self._movement.num_steps,
            max_wake_rows=self._movement.max_wake_rows,
            lcm_period=self._movement.lcm_period,
        )

        # Initialize an empty list to hold the SteadyProblems as they are generated.
        steady_problems_temp: list[SteadyProblem] = []

        # Iterate through the UnsteadyProblem's time steps.
        for step_id in range(self._num_steps):

            # Get the Airplanes and the OperatingPoint associated with this time step.
            these_airplanes = []
            for this_base_airplane in movement.airplanes:
                these_airplanes.append(this_base_airplane[step_id])
            this_operating_point = movement.operating_points[step_id]

            # Initialize the SteadyProblem at this time step.
            this_steady_problem = SteadyProblem(
                airplanes=these_airplanes, operating_point=this_operating_point
            )

            # Append this SteadyProblem to the temporary list.
            steady_problems_temp.append(this_steady_problem)

        # Store as tuple to prevent external mutation via .append(), .pop(), etc.
        self._steady_problems: tuple[SteadyProblem, ...] = tuple(steady_problems_temp)

    # --- Immutable: read only properties ---
    @property
    def movement(self) -> movements.movement.Movement:
        return self._movement

    @property
    def steady_problems(self) -> tuple[SteadyProblem, ...]:
        return self._steady_problems


class _CoupledUnsteadyProblem(_core.CoreUnsteadyProblem):
    """A class for coupled unsteady aerodynamics problems.

    This class extends CoreUnsteadyProblem to manage SteadyProblems for coupled
    simulations where the geometry at each time step depends on the solver's results
    from previous time steps.

    **Contains the following methods:**

    movement: The CoreMovement that defines the motion parameters for this problem.

    steady_problems: A tuple of SteadyProblems, one for each time step that has been
    initialized so far.

    get_steady_problem: Gets the SteadyProblem at a specified time step.

    initialize_next_problem: Initializes the next time step's SteadyProblem. Must be
    overridden by subclasses.
    """

    __slots__ = (
        "_movement",
        "_steady_problems",
    )

    def __init__(
        self,
        movement: _core.CoreMovement,
        initial_airplanes: list[geometry.airplane.Airplane],
        initial_operating_point: operating_point_mod.OperatingPoint,
    ) -> None:
        """The initialization method.

        Initializes the coupled unsteady problem with the first time step's geometry and
        the motion parameters from the provided CoreMovement.

        :param movement: A CoreMovement object that defines the motion parameters
            (delta_time, num_steps, max_wake_rows, lcm_period) for this problem.
        :param initial_airplanes: The list of Airplanes at the first time step.
        :param initial_operating_point: The OperatingPoint at the first time step.
        :return: None
        """
        self._movement = movement

        # Delegate shared initialization (validation, first_averaging_step computation,
        # load list initialization) to the core class. _CoupledUnsteadyProblems require
        # per step results to feed the coupling hook, so only_final_results is always
        # False.
        super().__init__(
            only_final_results=False,
            delta_time=self._movement.delta_time,
            num_steps=self._movement.num_steps,
            max_wake_rows=self._movement.max_wake_rows,
            lcm_period=self._movement.lcm_period,
        )

        # Coupled-specific state: a mutable list of SteadyProblems that grows as the
        # solver advances. Subclass initialize_next_problem overrides append to this
        # list. External code reads through the steady_problems tuple-view property to
        # preserve the read-only contract inherited from UnsteadyProblem. Seed with a
        # SteadyProblem built from the initial geometry so step zero is always ready.
        self._steady_problems: list[SteadyProblem] = [
            SteadyProblem(
                airplanes=initial_airplanes,
                operating_point=initial_operating_point,
            )
        ]

    # --- Immutable: read only properties ---
    @property
    def movement(self) -> _core.CoreMovement:
        return self._movement

    @property
    def steady_problems(self) -> tuple[SteadyProblem, ...]:
        return tuple(self._steady_problems)

    def get_steady_problem(self, step: int) -> SteadyProblem:
        """Get the SteadyProblem at a given time step.

        :param step: The time step index (zero indexed). Must be greater than or equal
            to zero and less than the total number of time steps.
        :return: The SteadyProblem at the specified time step.
        """
        step = _parameter_validation.int_in_range_return_int(
            step, "step", 0, True, len(self._steady_problems), False
        )

        return self._steady_problems[step]

    def initialize_next_problem(
        self,
        solver: CoupledUnsteadyRingVortexLatticeMethodSolver,
        step: int,
    ) -> None:
        """Initialize the next time step's SteadyProblem and perform per step work.

        Subclasses must override this method. It is invoked by the solver on every step,
        so subclasses are responsible for guarding any work that depends on a next step
        existing (such as building the next SteadyProblem) with ``step < self.num_steps
        - 1``. Per step work that should run on every step (such as recording the
        current step's loads) belongs outside that guard.

        :param solver: The CoupledUnsteadyRingVortexLatticeMethodSolver instance
            providing aerodynamic data from the current time step.
        :param step: The current time step index (zero indexed).
        :return: None
        :raises NotImplementedError: Always. Subclasses must override this method.
        """
        raise NotImplementedError("Subclasses must implement initialize_next_problem.")


# These are the permitted top-level keys for FreeFlightUnsteadyProblem's extra_xml
# injection-point dict. Each maps to an XML fragment that MuJoCoModel injects into the
# generated model XML at the matching location.
_EXTRA_XML_INJECTION_POINTS = frozenset(
    {"default", "asset", "visual", "worldbody", "body"}
)

# These are the permitted values for FreeFlightUnsteadyProblem's integrator parameter.
# Each names a MuJoCo numerical integrator that MuJoCoModel sets in the generated model
# XML's option element.
_MUJOCO_INTEGRATORS = frozenset({"Euler", "RK4", "implicit", "implicitfast"})

# These are the strongly coupled free-flight sub-iteration tunables. The relative and
# absolute tolerances form the mixed convergence test on the nondimensionalized state
# residual. The divergence tolerance guards the Aitken relaxation factor against a
# collapsing denominator. The initial relaxation factor under-relaxes the first update
# before the Aitken formula takes over.
_SUBITERATION_RELATIVE_TOLERANCE = 1e-6
_SUBITERATION_ABSOLUTE_TOLERANCE = 1e-10
_SUBITERATION_DIVERGENCE_TOLERANCE = 1e-20
_SUBITERATION_INITIAL_RELAXATION_FACTOR = 0.5


class FreeFlightUnsteadyProblem(_CoupledUnsteadyProblem):
    """A class used to contain problems with coupled unsteady aerodynamics and rigid
    body dynamics.

    **Contains the following methods:**

    only_final_results: Determines whether the solver will only calculate loads for the
    final time step or final cycle.

    num_steps: The number of time steps.

    delta_time: The time step size in seconds.

    first_averaging_step: The first time step included in cycle averaging.

    first_results_step: The first time step for which loads are calculated.

    max_wake_rows: The maximum chordwise wake rows per Wing.

    movement: The FreeFlightMovement that defines the motion parameters for this
    FreeFlightUnsteadyProblem.

    steady_problems: A tuple of SteadyProblems, one for each time step that has been
    initialized so far.

    get_steady_problem: Gets the SteadyProblem at a specified time step.

    initialize_next_problem: Initializes the next time step's SteadyProblem from rigid
    body dynamics.

    mass: The mass of the Airplane in kilograms.

    I_BP1_CgP1: The inertia matrix of the Airplane (in the first Airplane's body axes,
    relative to the first Airplane's CG) in kilogram square meters.

    k_max: The maximum number of strongly coupled sub-iterations per free-flight time
    step.

    external_loads_fn: A callable that computes additional forces and moments to apply
    to the Airplane during the simulation, or None.
    """

    __slots__ = (
        "_I_BP1_CgP1",
        "_mass",
        "_external_loads_fn",
        "_external_loads_validated",
        "_k_max",
        "_mujoco_model",
    )

    def __init__(
        self,
        movement: movements.free_flight_movement.FreeFlightMovement,
        mass: float | int,
        I_BP1_CgP1: np.ndarray | Sequence[Sequence[float | int]],
        external_loads_fn: (
            Callable[
                [
                    operating_point_mod.OperatingPoint,
                    geometry.airplane.Airplane,
                ],
                tuple[np.ndarray, np.ndarray],
            ]
            | None
        ) = None,
        k_max: int = 20,
        integrator: str = "RK4",
        extra_xml: dict[str, str] | None = None,
        mujoco_assets: dict[str, bytes] | None = None,
    ) -> None:
        """The initialization method.

        :param movement: The FreeFlightMovement that defines the prescribed Airplane
            geometry for this FreeFlightUnsteadyProblem. The initial Airplane and
            OperatingPoint are derived from the FreeFlightMovement at the first time
            step. The FreeFlightMovement must contain exactly one AirplaneMovement;
            multi-airplane free flight is not supported in this release.
        :param mass: A number (int or float) representing the mass of the Airplane. It
            must be greater than zero and will be converted internally to a float. The
            units are in kilograms. It must satisfy weight == mass * |g_E| within
            floating point tolerance, where weight is the Airplane's weight and g_E is
            the OperatingPoint's gravitational acceleration, which keeps the Airplane's
            weight, the supplied mass, and the gravitational field mutually consistent.
        :param I_BP1_CgP1: An array-like object of numbers (int or float) with shape
            (3,3) representing the inertia matrix of the Airplane (in the first
            Airplane's body axes, relative to the first Airplane's CG). It must be
            symmetric. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. The units are in kilogram square meters.
        :param external_loads_fn: A callable that computes additional forces and moments
            to apply to the Airplane during the simulation. It takes an OperatingPoint
            and an Airplane and returns a tuple of two (3,) ndarrays of floats: the
            additional force (in wind axes, in Newtons) and the additional moment (in
            wind axes, relative to the first Airplane's CG, in Newton meters). The
            return value is validated on the callable's first invocation; a return that
            is not a pair of (3,) finite numeric vectors raises a descriptive error. The
            physical correctness of the forces and moments themselves is not checked.
            Setting this to None applies no additional loads. The default is None. This
            is the only mechanism for non-aerodynamic loads in free flight: the
            OperatingPoint's externalFX_W is never applied and must be zero.
        :param k_max: An int giving the maximum number of strongly coupled sub-
            iterations per free-flight time step. Each time step drives the aerodynamic
            loads and the rigid body state to mutual consistency with an Aitken-relaxed
            fixed-point sub-iteration; this caps the iterations spent before the step is
            accepted. A step that reaches the cap without converging is accepted anyway,
            with a warning. It must be greater than zero and will be converted
            internally to an int. Raising it trades more work per step for a tighter
            solve when convergence is slow, which can happen for a very light vehicle (a
            large added-mass ratio). The default is 20.
        :param integrator: A str naming the MuJoCo integrator used to advance the rigid
            body dynamics. It must be one of "Euler", "RK4", "implicit", or
            "implicitfast". The choice is baked into the generated MuJoCo model XML, so
            it survives saving and loading. "RK4" is accurate for smooth dynamics but
            handles contacts poorly, so runs with collision geometry (injected via
            extra_xml) should prefer "implicit" or "implicitfast". The default is "RK4".
        :param extra_xml: A dict mapping injection point names to XML fragment strings
            to inject into the MuJoCo model's XML. Supported keys are "default",
            "asset", "visual", "worldbody", and "body". Setting this to None injects no
            extra XML. The default is None. The argument is checked to be a dict (or
            None) whose keys are supported injection points and whose values are
            strings; the XML fragments themselves are not validated, which is left to
            MuJoCo, so this is an advanced-user parameter.
        :param mujoco_assets: A dict mapping virtual filenames to their binary contents
            for the MuJoCo model. Setting this to None provides no extra assets. The
            default is None. The argument is checked to be a dict (or None) mapping
            string filenames to bytes; whether a referenced asset is actually supplied
            is left to MuJoCo, so this is an advanced-user parameter. A
            FreeFlightUnsteadyProblem built with mujoco_assets cannot be saved: save()
            raises, because the saved engine is rebuilt on load from the stored XML
            alone, whose asset references would be unresolvable.
        :return: None
        """
        if not isinstance(movement, free_flight_movement.FreeFlightMovement):
            raise TypeError("movement must be a FreeFlightMovement.")

        if len(movement.airplane_movements) != 1:
            raise ValueError(
                "movement must have exactly one AirplaneMovement. "
                "Multi-airplane free flight is not supported in this release."
            )

        # Derive the initial Airplane and OperatingPoint from the FreeFlightMovement's
        # first time step. SteadyProblem mutates each Panel's GP1_CgP1 attributes once,
        # so the initial Airplane must be a fresh object, which
        # generate_airplane_at_time_step returns.
        initial_airplane_movement = movement.airplane_movements[0]
        initial_airplane = initial_airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=movement.delta_time
        )
        initial_operating_point = movement.operating_point_movement.operating_points[0]

        super().__init__(
            movement=movement,
            initial_airplanes=[initial_airplane],
            initial_operating_point=initial_operating_point,
        )

        I_BP1_CgP1 = _parameter_validation.m_by_n_number_arrayLike_return_float(
            I_BP1_CgP1, "I_BP1_CgP1", 3, 3
        )
        if not np.allclose(I_BP1_CgP1, I_BP1_CgP1.T):
            raise ValueError("I_BP1_CgP1 must be symmetric.")
        self._I_BP1_CgP1 = I_BP1_CgP1
        self._I_BP1_CgP1.flags.writeable = False

        self._mass = _parameter_validation.number_in_range_return_float(
            mass,
            "mass",
            min_val=0.0,
            min_inclusive=False,
        )

        # The Airplane's weight, the supplied mass, and the gravitational field must be
        # mutually consistent: weight == mass * |g_E|. The free-flight dynamics derive
        # the gravitational force as mass * g_E, while the Airplane's weight is used
        # elsewhere (for example by the trim functions), so a disagreement would be a
        # silent modeling error. g_E is constant across the run, so checking the initial
        # OperatingPoint suffices. A zero g_E (the default, no gravitational field)
        # consistently requires a zero weight.
        expected_weight = self._mass * float(
            np.linalg.norm(initial_operating_point.g_E)
        )
        if not np.isclose(initial_airplane.weight, expected_weight):
            raise ValueError(
                "The Airplane's weight must equal mass * |g_E| within floating point "
                f"tolerance: the Airplane's weight is {initial_airplane.weight} N, but "
                f"mass * |g_E| is {expected_weight} N. Set the Airplane's weight, the "
                "mass, and the OperatingPoint's g_E so they agree (for a zero-gravity "
                "simulation, leave both g_E and the weight at zero)."
            )

        # The free-flight dynamics never apply externalFX_W: Non-aerodynamic loads enter
        # free flight only through external_loads_fn, which is strictly more capable
        # (full force and moment, time varying). A nonzero value would be silently
        # ignored, so raise instead.
        if initial_operating_point.externalFX_W != 0.0:
            raise ValueError(
                "The OperatingPoint's externalFX_W must be zero for free flight. The "
                "free-flight dynamics never apply externalFX_W; use external_loads_fn "
                "to model thrust or other additional loads."
            )

        if external_loads_fn is not None and not callable(external_loads_fn):
            raise TypeError("external_loads_fn must be callable or None.")
        self._external_loads_fn = external_loads_fn

        # Tracks whether the external_loads_fn's return structure has been validated.
        # The return contract is static across time steps, so it is checked once, on the
        # function's first invocation in initialize_next_problem, rather than every
        # step.
        self._external_loads_validated = False

        # This is the cap on the strongly coupled sub-iteration count per free-flight
        # time step. The remaining sub-iteration tunables stay module constants by
        # design (they are universal dimensionless constants), so the cap is the only
        # one exposed.
        self._k_max = _parameter_validation.int_in_range_return_int(
            k_max,
            "k_max",
            min_val=0,
            min_inclusive=False,
        )

        # Validate the integrator name against the supported MuJoCo integrators. The
        # MuJoCoModel sets the validated name in the generated model XML's option
        # element, so the choice survives saving and loading.
        integrator = _parameter_validation.str_return_str(integrator, "integrator")
        if integrator not in _MUJOCO_INTEGRATORS:
            raise ValueError(
                f"integrator '{integrator}' is not a supported MuJoCo integrator; "
                f"expected one of {sorted(_MUJOCO_INTEGRATORS)}."
            )

        # Validate the extra_xml injection-point dict (it must be a dict, or None, whose
        # keys are permitted injection points and whose values are XML fragment strings)
        # and rebuild it from the validated values. Deeper XML correctness is left to
        # MuJoCo's own parser. The integrator, extra_xml, and mujoco_assets arguments
        # are the only raw user input forwarded to the MuJoCoModel, which performs no
        # validation of its own.
        if extra_xml is not None:
            if not isinstance(extra_xml, dict):
                raise TypeError("extra_xml must be a dict or None.")
            validated_extra_xml: dict[str, str] = {}
            for key, value in extra_xml.items():
                if key not in _EXTRA_XML_INJECTION_POINTS:
                    raise ValueError(
                        f"extra_xml key '{key}' is not a permitted injection point; "
                        f"expected one of {sorted(_EXTRA_XML_INJECTION_POINTS)}."
                    )
                validated_extra_xml[key] = _parameter_validation.str_return_str(
                    value, f"extra_xml['{key}']"
                )
            extra_xml = validated_extra_xml

        # Validate the mujoco_assets dict (it must be a dict, or None, mapping str
        # filenames to bytes, where each filename is a bare basename with a nonempty
        # extension). Whether a referenced asset is actually supplied is left to
        # MuJoCo's own parser.
        if mujoco_assets is not None:
            if not isinstance(mujoco_assets, dict):
                raise TypeError("mujoco_assets must be a dict or None.")
            for filename, contents in mujoco_assets.items():
                if not isinstance(filename, str):
                    raise TypeError(
                        "mujoco_assets keys must be str filenames, not "
                        f"{type(filename).__name__}."
                    )

                # Windows path rules recognize both separator styles and drive prefixes,
                # so a single PureWindowsPath check rejects every path shape on all
                # platforms. The basename requirement keeps machine-specific paths out
                # of saved files.
                if PureWindowsPath(filename).name != filename:
                    raise ValueError(
                        f"mujoco_assets key '{filename}' must be a bare filename "
                        "with no path separators or drive prefixes."
                    )

                # MuJoCo selects its asset decoder from the extension, so an
                # extension-less virtual filename fails to compile even with an explicit
                # content_type attribute.
                stem, _, extension = filename.rpartition(".")
                if not stem or not extension:
                    raise ValueError(
                        f"mujoco_assets key '{filename}' must be a filename with a "
                        "nonempty extension."
                    )

                if not isinstance(contents, bytes):
                    raise TypeError(
                        f"mujoco_assets['{filename}'] must be bytes, not "
                        f"{type(contents).__name__}."
                    )

        self._mujoco_model = _mujoco_model.MuJoCoModel(
            name=initial_airplane.name,
            mass=self._mass,
            omegas_BP1__E=initial_operating_point.omegas_BP1__E,
            T_pas_BP1_CgP1_to_E_CgP1=initial_operating_point.T_pas_BP1_CgP1_to_E_CgP1,
            vCg_E__E=-1
            * _transformations.apply_T_to_vectors(
                initial_operating_point.T_pas_GP1_CgP1_to_E_CgP1,
                initial_operating_point.vInf_GP1__E,
                is_position=False,
            ),
            I_BP1_CgP1=self._I_BP1_CgP1,
            delta_time=movement.delta_time,
            integrator=integrator,
            extra_xml=extra_xml,
            mujoco_assets=mujoco_assets,
        )

        # A file reference that mujoco_assets does not cover was resolved from the local
        # filesystem when the model compiled, which ties the simulation to a
        # machine-specific path and would make it unsaveable. Reject it now, before an
        # expensive run, rather than surprising the user at save time. This check has to
        # run after the MuJoCoModel is constructed because the references live in the
        # generated model XML.
        uncovered_file_references = self._mujoco_model.uncovered_file_references()
        if uncovered_file_references:
            raise ValueError(
                "The MuJoCo model XML references files that mujoco_assets does not "
                f"cover: {uncovered_file_references}. Supply each referenced file's "
                "contents through mujoco_assets instead of referencing it by path, "
                "which would tie the simulation to files on this machine and make it "
                "unsaveable."
            )

    # --- Immutable: read only properties ---
    @property
    def I_BP1_CgP1(self) -> np.ndarray:
        return self._I_BP1_CgP1

    @property
    def mass(self) -> float:
        return self._mass

    @property
    def k_max(self) -> int:
        return self._k_max

    @property
    def external_loads_fn(
        self,
    ) -> (
        Callable[
            [
                operating_point_mod.OperatingPoint,
                geometry.airplane.Airplane,
            ],
            tuple[np.ndarray, np.ndarray],
        ]
        | None
    ):
        return self._external_loads_fn

    @property
    def _free_flight_movement(self) -> free_flight_movement.FreeFlightMovement:
        """Type narrowed view of the inherited _movement attribute.

        The parent stores _movement as a CoreMovement (widened to let subclasses pass
        their own variants). __init__ accepts only a FreeFlightMovement, so the cast
        here is safe.

        :return: The _movement narrowed to FreeFlightMovement.
        """
        return cast(free_flight_movement.FreeFlightMovement, self._movement)

    @staticmethod
    def _validate_external_loads_return(result: object) -> None:
        """Validates the structure of a value returned by the external_loads_fn.

        The external_loads_fn must return a (force, moment) pair, where each is a (3,)
        array-like of finite numbers (the force in wind axes in Newtons, the moment in
        wind axes relative to the first Airplane's CG in Newton meters). This checks
        only the shape and finiteness of the return value, not the physical correctness
        of the forces and moments, which cannot be validated in general. It is called
        once, on the external_loads_fn's first invocation, because the return contract
        is static across time steps.

        :param result: The object returned by a call to the external_loads_fn.
        :return: None
        """
        if not isinstance(result, (tuple, list, np.ndarray)):
            raise TypeError(
                "external_loads_fn must return a (force, moment) pair as a tuple, "
                f"list, or ndarray, but got {type(result).__name__}."
            )
        if len(result) != 2:
            raise ValueError(
                "external_loads_fn must return a (force, moment) pair, but its return "
                f"value had {len(result)} items."
            )

        # Each component must be a (3,) vector of finite numbers. Defer the shape,
        # finiteness, and numeric-type checks to the package's standard validator.
        _parameter_validation.threeD_number_vectorLike_return_float(
            result[0], "external_loads_fn's returned force"
        )
        _parameter_validation.threeD_number_vectorLike_return_float(
            result[1], "external_loads_fn's returned moment"
        )

    def _assemble_interval_loads(
        self,
        operating_point: operating_point_mod.OperatingPoint,
        airplane: geometry.airplane.Airplane,
        step: int,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Assembles the Earth-axes interval load to apply over the next time step.

        In the free flight phase, returns the total force and moment (in Earth axes,
        relative to the first Airplane's CG) to apply to the MuJoCo model: the
        aerodynamic loads plus any external loads, transformed from wind axes to Earth
        axes, with the weight added. In the prescribed phase, returns None so the loads
        are withheld and the body coasts at its initial trimmed condition (MuJoCo's
        internal gravity is off and weight is applied through the loads, so withholding
        them leaves zero net force and the body keeps its initial velocity and
        orientation), letting the wake develop at a steady operating condition before
        the rigid body dynamics are released. The external_loads_fn is still invoked
        once, on the first step, so its return structure is validated fail-fast rather
        than only when the free flight phase begins.

        :param operating_point: The OperatingPoint at which the loads are evaluated.
        :param airplane: The first Airplane, carrying the aerodynamic loads the solver
            calculated at this operating point.
        :param step: The current time step index (zero indexed).
        :return: A tuple of two (3,) ndarrays of floats representing the total force and
            moment (in Earth axes, relative to the first Airplane's CG) to apply, in
            Newtons and Newton meters, or None to withhold the loads in the prescribed
            phase.
        """
        in_free_phase = step >= self._free_flight_movement.prescribed_num_steps

        if not in_free_phase:
            # In the prescribed phase the loads are withheld, but the external_loads_fn
            # is still invoked once, on the first step, so its return structure is
            # validated fail-fast rather than only when the free flight phase begins.
            if (
                step == 0
                and self._external_loads_fn is not None
                and not self._external_loads_validated
            ):
                self._validate_external_loads_return(
                    self._external_loads_fn(operating_point, airplane)
                )
                self._external_loads_validated = True
            return None

        aeroForces_W = airplane.forces_W
        assert aeroForces_W is not None
        aeroMoments_W_CgP1 = airplane.moments_W_CgP1
        assert aeroMoments_W_CgP1 is not None

        # Start from the aerodynamic loads and add any external loads from the
        # external_loads_fn.
        if self._external_loads_fn is not None:
            external_result = self._external_loads_fn(operating_point, airplane)
            # Validate the return structure once, on the first invocation.
            if not self._external_loads_validated:
                self._validate_external_loads_return(external_result)
                self._external_loads_validated = True
            externalForces_W, externalMoments_W_CgP1 = external_result
            totalForces_W = aeroForces_W + externalForces_W
            totalMoments_W_CgP1 = aeroMoments_W_CgP1 + externalMoments_W_CgP1
        else:
            totalForces_W = aeroForces_W
            totalMoments_W_CgP1 = aeroMoments_W_CgP1

        # Transform loads from wind axes to Earth axes.
        T_pas_W_CgP1_to_E_CgP1 = operating_point.T_pas_W_CgP1_to_E_CgP1
        totalForces_E = _transformations.apply_T_to_vectors(
            T_pas_W_CgP1_to_E_CgP1, totalForces_W, is_position=False
        )
        totalMoments_E_CgP1 = _transformations.apply_T_to_vectors(
            T_pas_W_CgP1_to_E_CgP1, totalMoments_W_CgP1, is_position=False
        )

        # Add the weight force in Earth axes. The gravitational force is mass * g_E,
        # which is zero when g_E is zero (no gravitational field).
        totalForces_E = totalForces_E + self._mass * operating_point.g_E

        return totalForces_E, totalMoments_E_CgP1

    def _advance_body(
        self, interval_loads_E: tuple[np.ndarray, np.ndarray] | None
    ) -> _mujoco_model.MuJoCoState:
        """Applies the interval load (if any) and steps the dynamics, returning the new
        state.

        Applies the given Earth-axes force and moment to the MuJoCo model, steps the
        rigid body dynamics forward by one time step, and returns the new state. When
        the load is None (the prescribed phase) no load is applied and the body coasts.
        The MuJoCo model is stepped from its current state; restoring it to a frozen
        snapshot, when re-stepping within a strongly coupled sub-iteration, is the
        caller's responsibility.

        :param interval_loads_E: A tuple of two (3,) ndarrays of floats representing the
            force and moment (in Earth axes, relative to the first Airplane's CG) to
            apply over the interval, in Newtons and Newton meters, or None to apply no
            load.
        :return: The MuJoCo model's state after the step, as returned by
            MuJoCoModel.get_state.
        """
        if interval_loads_E is not None:
            forces_E, moments_E_CgP1 = interval_loads_E
            self._mujoco_model.apply_loads(forces_E, moments_E_CgP1)
        self._mujoco_model.step()

        return self._mujoco_model.get_state()

    def _operating_point_from_state(
        self,
        state: _mujoco_model.MuJoCoState,
        reference_operating_point: operating_point_mod.OperatingPoint,
    ) -> operating_point_mod.OperatingPoint:
        """Builds the next OperatingPoint from a MuJoCo state.

        Derives the speed, angle of attack, sideslip angle, and Earth-to-body Euler
        angles from the rigid body state, and carries the environmental quantities
        (fluid density, surface geometry, external force, kinematic viscosity, and
        gravity) across from the reference OperatingPoint unchanged.

        :param state: A MuJoCo state, as returned by MuJoCoModel.get_state.
        :param reference_operating_point: The current OperatingPoint, whose
            environmental quantities are carried across to the new OperatingPoint.
        :return: The OperatingPoint describing the new rigid body state.
        """
        angles_E_to_BP1_izyx = _transformations.R_to_angles_izyx(
            state["R_pas_E_to_BP1"]
        )

        return self._build_next_operating_point(
            position_E_Eo=state["position_E_Eo"],
            angles_E_to_BP1_izyx=angles_E_to_BP1_izyx,
            velocity_E__E=state["velocity_E__E"],
            omegas_BP1__E=state["omegas_BP1__E"],
            reference_operating_point=reference_operating_point,
        )

    def _operating_point_from_vector(
        self,
        x: np.ndarray,
        reference_operating_point: operating_point_mod.OperatingPoint,
    ) -> operating_point_mod.OperatingPoint:
        """Builds the next OperatingPoint from a rigid body state 12-vector.

        The 12-vector is a relaxed trial state of the sub-iteration, expressed in Earth
        axes with angular quantities in radians (see _state_to_vector). Unlike a MuJoCo
        state, a relaxed trial state need not correspond to any state MuJoCo produced,
        so the OperatingPoint is built from the 12-vector directly. The attitude is
        converted from radians to degrees and the angular rate is rotated from Earth
        axes back into the first Airplane's body axes and converted to degrees per
        second, the form the OperatingPoint expects.

        :param x: A (12,) ndarray of floats representing the rigid body state 12-vector,
            as produced by _state_to_vector.
        :param reference_operating_point: The current OperatingPoint, whose
            environmental quantities are carried across to the new OperatingPoint.
        :return: The OperatingPoint describing the trial rigid body state.
        """
        position_E_Eo = x[0:3]
        anglesRad_E_to_BP1_izyx = x[3:6]
        velocity_E__E = x[6:9]
        omegasRad_E__E = x[9:12]

        angles_E_to_BP1_izyx = np.rad2deg(anglesRad_E_to_BP1_izyx)

        # Rotate the angular rate from Earth axes back into the first Airplane's body
        # axes and convert to degrees per second for the OperatingPoint.
        T_pas_E_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
            angles=angles_E_to_BP1_izyx,
            passive=True,
            intrinsic=True,
            order="zyx",
        )
        omegas_BP1__E = np.rad2deg(
            _transformations.apply_T_to_vectors(
                T_pas_E_CgP1_to_BP1_CgP1, omegasRad_E__E, is_position=False
            )
        )

        return self._build_next_operating_point(
            position_E_Eo=position_E_Eo,
            angles_E_to_BP1_izyx=angles_E_to_BP1_izyx,
            velocity_E__E=velocity_E__E,
            omegas_BP1__E=omegas_BP1__E,
            reference_operating_point=reference_operating_point,
        )

    @staticmethod
    def _build_next_operating_point(
        position_E_Eo: np.ndarray,
        angles_E_to_BP1_izyx: np.ndarray,
        velocity_E__E: np.ndarray,
        omegas_BP1__E: np.ndarray,
        reference_operating_point: operating_point_mod.OperatingPoint,
    ) -> operating_point_mod.OperatingPoint:
        """Builds the next OperatingPoint from already-extracted rigid body components.

        Derives the speed, angle of attack, and sideslip angle from the velocity and
        attitude, and carries the environmental quantities (fluid density, surface
        geometry, external force, kinematic viscosity, and gravity) across from the
        reference OperatingPoint unchanged. The attitude and angular rate are supplied
        in the form the OperatingPoint expects: the intrinsic zyx Euler angles from
        Earth axes to the first Airplane's body axes in degrees, and the body angular
        rate in the first Airplane's body axes in degrees per second.

        :param position_E_Eo: A (3,) ndarray of floats representing the first Airplane's
            CG position (in Earth axes, relative to the Earth origin) in meters.
        :param angles_E_to_BP1_izyx: A (3,) ndarray of floats representing the intrinsic
            zyx Euler angles from Earth axes to the first Airplane's body axes, in
            degrees.
        :param velocity_E__E: A (3,) ndarray of floats representing the first Airplane's
            CG velocity (in Earth axes, observed from the Earth frame) in meters per
            second.
        :param omegas_BP1__E: A (3,) ndarray of floats representing the body angular
            rate (in the first Airplane's body axes, observed from the Earth frame) in
            degrees per second.
        :param reference_operating_point: The current OperatingPoint, whose
            environmental quantities are carried across to the new OperatingPoint.
        :return: The OperatingPoint describing the new rigid body state.
        """
        vCg__E = float(np.linalg.norm(velocity_E__E))
        T_pas_E_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
            angles=angles_E_to_BP1_izyx,
            passive=True,
            intrinsic=True,
            order="zyx",
        )
        vInf_E__E = -velocity_E__E
        vInf_BP1__E = _transformations.apply_T_to_vectors(
            T_pas_E_CgP1_to_BP1_CgP1, vInf_E__E, is_position=False
        )
        alpha, beta = _transformations.alpha_and_beta_from_vInf_BP1(vInf_BP1__E, vCg__E)

        return operating_point_mod.OperatingPoint(
            rho=reference_operating_point.rho,
            vCg__E=vCg__E,
            alpha=alpha,
            beta=beta,
            angles_E_to_BP1_izyx=angles_E_to_BP1_izyx,
            CgP1_E_Eo=position_E_Eo,
            surfaceNormal_E=reference_operating_point.surfaceNormal_E,
            surfacePoint_E_Eo=reference_operating_point.surfacePoint_E_Eo,
            externalFX_W=reference_operating_point.externalFX_W,
            nu=reference_operating_point.nu,
            g_E=reference_operating_point.g_E,
            omegas_BP1__E=omegas_BP1__E,
        )

    @staticmethod
    def _state_to_vector(state: _mujoco_model.MuJoCoState) -> np.ndarray:
        """Converts a MuJoCo state to the rigid body state 12-vector.

        The 12-vector concatenates the four blocks of the rigid body state, all
        expressed in Earth axes with angular quantities in radians, so the sub-
        iteration's residuals and updates are plain vector arithmetic: the CG position
        (relative to the Earth origin, in meters), the intrinsic zyx Euler angles from
        Earth axes to the first Airplane's body axes (in radians), the CG velocity
        (observed from the Earth frame, in meters per second), and the body angular rate
        (observed from the Earth frame, in radians per second). The MuJoCo state carries
        the attitude as an orientation matrix and the angular rate in the first
        Airplane's body axes in degrees per second, so this method extracts the Euler
        angles, rotates the rate into Earth axes, and converts both angular quantities
        to radians.

        :param state: A MuJoCo state, as returned by MuJoCoModel.get_state.
        :return: A (12,) ndarray of floats representing the rigid body state 12-vector.
        """
        R_pas_E_to_BP1 = state["R_pas_E_to_BP1"]

        position_E_Eo = state["position_E_Eo"]
        anglesRad_E_to_BP1_izyx = np.deg2rad(
            _transformations.R_to_angles_izyx(R_pas_E_to_BP1)
        )
        velocity_E__E = state["velocity_E__E"]
        # Rotate the angular rate from the first Airplane's body axes into Earth axes
        # (v_E = R_pas_E_to_BP1.T @ v_BP1, a proper rotation that preserves the
        # axial-vector sign) and convert to radians per second.
        omegasRad_E__E = R_pas_E_to_BP1.T @ np.deg2rad(state["omegas_BP1__E"])

        return np.concatenate(
            [position_E_Eo, anglesRad_E_to_BP1_izyx, velocity_E__E, omegasRad_E__E]
        )

    def _build_relaxation_weights(self) -> np.ndarray:
        """Builds the diagonal of the sub-iteration weighting matrix D.

        D nondimensionalizes the rigid body state 12-vector so the convergence test and
        the Aitken inner products combine its four blocks (position, attitude, velocity,
        and angular rate) on a common scale despite their differing units and their
        differing powers of the time step under a load perturbation. It is diagonal and
        block constant, ordered like the state as (position, attitude, velocity, angular
        rate):

            D = diag( (1 / L_ref) I3, (1 / theta_ref) I3,
                      (delta_time / L_ref) I3, (delta_time / theta_ref) I3 )

        The reference length L_ref is the first Airplane's reference chord. The reference
        angle theta_ref is mass * L_ref**2 / I_bar, with I_bar the mean of the body
        inertia's principal moments (one third of its trace), a rotation invariant chosen
        so D stays constant as the body reorients. The delta_time factors on the velocity
        and angular-rate blocks convert a rate into the increment it produces over one
        step, lifting those blocks to match the configuration blocks. All inputs are
        fixed for a run, so D is constant.

        :return: A (12,) ndarray of floats representing the diagonal of the weighting
            matrix D.
        """
        airplane = self._free_flight_movement.airplanes[0][0]
        L_ref = airplane.c_ref
        assert L_ref is not None

        mean_inertia = float(np.trace(self._I_BP1_CgP1)) / 3.0
        theta_ref = self._mass * L_ref**2 / mean_inertia

        delta_time = self.delta_time

        return np.concatenate(
            [
                np.full(3, 1.0 / L_ref),
                np.full(3, 1.0 / theta_ref),
                np.full(3, delta_time / L_ref),
                np.full(3, delta_time / theta_ref),
            ]
        )

    def _commit_next_problem(
        self,
        next_operating_point: operating_point_mod.OperatingPoint,
        step: int,
    ) -> None:
        """Commits the next time step's OperatingPoint and SteadyProblem.

        Appends the new OperatingPoint to the operating point movement's history and
        appends the next SteadyProblem (the prescribed Airplane geometry for the next
        step paired with the new OperatingPoint) to the problem's steady problems. These
        are the once-per-time-step side effects of advancing to the next step.

        :param next_operating_point: The OperatingPoint for the next time step.
        :param step: The current time step index (zero indexed).
        :return: None
        """
        self._free_flight_movement.operating_point_movement.operating_points.append(
            next_operating_point
        )

        next_airplane = self._free_flight_movement.airplanes[0][step + 1]

        next_steady_problem = SteadyProblem(
            airplanes=[next_airplane],
            operating_point=next_operating_point,
        )
        self._steady_problems.append(next_steady_problem)

    def _advance_strongly_coupled(
        self,
        solver: FreeFlightUnsteadyRingVortexLatticeMethodSolver,
        current_airplane: geometry.airplane.Airplane,
        current_operating_point: operating_point_mod.OperatingPoint,
        step: int,
    ) -> None:
        """Advances the body to the next step by the strongly coupled sub-iteration.

        Within one free flight time step, drives the aerodynamic loads and the rigid
        body state to mutual consistency by an Aitken-relaxed fixed-point iteration on
        the rigid body state 12-vector. The snapshot MuJoCo state and the solver's
        frozen step data are saved once. Each iteration evaluates the aerodynamic loads
        at the current trial state, integrates the body from the snapshot under those
        loads, and relaxes the trial state toward the dynamics-propagated state until
        the weighted residual passes the convergence test or the iteration cap is
        reached. The accepted state is the final dynamics-propagated state (a genuine
        MuJoCo trajectory output), and its OperatingPoint is committed as the next
        step's. The solver's current-step state is then restored for the inherited run
        loop's official wake build.

        :param solver: The FreeFlightUnsteadyRingVortexLatticeMethodSolver providing the
            aerodynamic evaluation at each trial state.
        :param current_airplane: The first Airplane at the current step, carrying the
            current step's solved aerodynamic loads, used to seed the iteration.
        :param current_operating_point: The current step's OperatingPoint, the reference
            for the environmental quantities carried into each trial OperatingPoint.
        :param step: The current time step index (zero indexed).
        :return: None
        """
        # Freeze: while the model is still at it, snapshot the MuJoCo state and the
        # current step's state 12-vector, the solver's frozen step data, and the
        # weighting. Build the transient next-step SteadyProblem over a scratch copy of
        # the prescribed next-step Airplane, so the trials evaluate aerodynamics on the
        # copy and the canonical Airplane's set-once Panel coordinates are reserved for
        # the official SteadyProblem committed once the solve accepts. Its
        # OperatingPoint is a placeholder. The trial OperatingPoint is supplied to each
        # trial.
        snapshot = self._mujoco_model.save_state()
        snapshot_x = self._state_to_vector(self._mujoco_model.get_state())
        next_airplane = self._free_flight_movement.airplanes[0][step + 1]
        transient_next_steady_problem = SteadyProblem(
            airplanes=[copy.deepcopy(next_airplane)],
            operating_point=current_operating_point,
        )
        solver.freeze_substep(transient_next_steady_problem)
        weights = self._build_relaxation_weights()

        # Seed the iteration with the loose step, S(x_n, l_n), evaluated from the
        # snapshot.
        snapshot_interval_loads = self._assemble_interval_loads(
            current_operating_point, current_airplane, step
        )
        trial_x = self._state_to_vector(self._advance_body(snapshot_interval_loads))

        relaxation_factor = _SUBITERATION_INITIAL_RELAXATION_FACTOR
        previous_residual: np.ndarray | None = None
        residual_norm = 0.0
        converged = False

        for iteration in range(self.k_max + 1):
            # Aerodynamic loads at the trial state: build the trial OperatingPoint,
            # evaluate the aerodynamics at it, and assemble the Earth-axes interval
            # load.
            trial_operating_point = self._operating_point_from_vector(
                trial_x, current_operating_point
            )
            trial_airplane = solver.evaluate_trial_aero_loads(
                trial_operating_point, step
            )
            interval_loads = self._assemble_interval_loads(
                trial_operating_point, trial_airplane, step
            )

            # Integrate the body from the snapshot under the trial load.
            self._mujoco_model.restore_state(snapshot)
            propagated_x = self._state_to_vector(self._advance_body(interval_loads))

            residual = propagated_x - trial_x
            increment = trial_x - snapshot_x
            residual_norm = _fixed_point_relaxation.weighted_norm(weights, residual)

            # The relaxation factor applied this iteration: the initial factor on the
            # first iteration, the Aitken factor thereafter.
            if iteration > 0:
                assert previous_residual is not None
                relaxation_factor = _fixed_point_relaxation.aitken_relaxation_factor(
                    weights,
                    residual,
                    previous_residual,
                    relaxation_factor,
                    _SUBITERATION_INITIAL_RELAXATION_FACTOR,
                    _SUBITERATION_DIVERGENCE_TOLERANCE,
                )

            _logger.debug(
                _logging.indent()
                + "Free flight step %d, sub-iteration %d: weighted residual norm "
                "%#.3G, relaxation factor %#.3G",
                step,
                iteration,
                residual_norm,
                relaxation_factor,
            )

            if _fixed_point_relaxation.is_converged(
                weights,
                residual,
                increment,
                _SUBITERATION_RELATIVE_TOLERANCE,
                _SUBITERATION_ABSOLUTE_TOLERANCE,
            ):
                converged = True
                break

            if iteration == self.k_max:
                break

            trial_x = trial_x + relaxation_factor * residual
            previous_residual = residual

        if not converged:
            _logger.warning(
                _logging.indent()
                + "Free flight sub-iteration at step %d reached the %d-iteration "
                "cap without converging",
                step,
                self.k_max,
            )
            _logger.warning(
                _logging.indent()
                + "Accepting the capped iterate with weighted residual norm %#.3G",
                residual_norm,
            )

        # Accept: The final dynamics-propagated state, at which the model now sits, is
        # the next step's state. Commit its OperatingPoint, then restore the solver's
        # current-step state for the official wake build in the inherited run loop.
        accepted_operating_point = self._operating_point_from_state(
            self._mujoco_model.get_state(), current_operating_point
        )
        self._commit_next_problem(accepted_operating_point, step)
        solver.restore_substep(step)

    def initialize_next_problem(
        self,
        solver: CoupledUnsteadyRingVortexLatticeMethodSolver,
        step: int,
    ) -> None:
        """Initializes the next time step's SteadyProblem from rigid body dynamics.

        On every step except the last one, advances the MuJoCo dynamics and creates the
        next SteadyProblem with the new OperatingPoint and the prescribed Airplane
        geometry for the next step. During the prescribed phase, the loads are withheld
        so the body coasts at its initial trimmed condition while the wake develops, and
        the body advances once. During the free flight phase (once the step index
        reaches the movement's prescribed_num_steps), the aerodynamic loads and the
        rigid body state are driven to mutual consistency by the strongly coupled sub-
        iteration before the next step is committed.

        :param solver: The CoupledUnsteadyRingVortexLatticeMethodSolver instance
            providing aerodynamic data from the current time step.
        :param step: The current time step index (zero indexed).
        :return: None
        """
        current_airplane = solver.current_airplanes[0]
        current_operating_point = solver.current_operating_point

        if step < self.num_steps - 1:
            if step >= self._free_flight_movement.prescribed_num_steps:
                # Free flight phase: drive the loads and the body state to mutual
                # consistency before committing the next step.
                self._advance_strongly_coupled(
                    cast(
                        "FreeFlightUnsteadyRingVortexLatticeMethodSolver",
                        solver,
                    ),
                    current_airplane,
                    current_operating_point,
                    step,
                )
            else:
                # Prescribed phase: The loads are withheld, so the body coasts at its
                # trimmed condition while the wake develops. Advance once and commit the
                # next step's SteadyProblem from the new state.
                interval_loads_E = self._assemble_interval_loads(
                    current_operating_point, current_airplane, step
                )
                new_state = self._advance_body(interval_loads_E)
                next_operating_point = self._operating_point_from_state(
                    new_state, current_operating_point
                )
                self._commit_next_problem(next_operating_point, step)


# These are tolerances for the torsional spring-damper ODE integration in
# _spring_numerical_ode. They are a few orders of magnitude stricter than scipy's
# defaults because the integration is re-seeded from the previous state at every outer
# time step, so its local errors compound across a simulation, and because a loose
# absolute tolerance would swamp small torsional responses. The absolute tolerance
# bounds both state components, the torsional angle (radians) and its time derivative
# (rad/s). Loosen these only for local debugging (for example, a nearly massless wing
# makes the ODE stiff and the strict tolerances expensive).
_SPRING_ODE_RELATIVE_TOLERANCE = 1e-6
_SPRING_ODE_ABSOLUTE_TOLERANCE_RAD = 1e-9


class AeroelasticUnsteadyProblem(_CoupledUnsteadyProblem):
    """A subclass of _CoupledUnsteadyProblem used to couple aeroelastic wing
    deformations with unsteady aerodynamics.

    This class couples aerodynamic loads with wing structural dynamics (spring-mass-
    damper system) to simulate aeroelastic deformation. Each time step, wing
    deformations are calculated based on the combined effects of each strip's
    aerodynamic, inertial, and spring-damper restoring torsional moments, all taken as y
    components (in the first Airplane's geometry axes) of moments relative to the
    strip's leading edge point.

    **Contains the following methods:**

    only_final_results: Determines whether the solver will only calculate loads for the
    final time step or final cycle.

    num_steps: The number of time steps.

    delta_time: The time step size in seconds.

    first_averaging_step: The first time step included in cycle averaging.

    first_results_step: The first time step for which loads are calculated.

    max_wake_rows: The maximum chordwise wake rows per Wing.

    movement: The AeroelasticMovement that defines the motion parameters for this
    AeroelasticUnsteadyProblem.

    steady_problems: A tuple of SteadyProblems, one for each time step that has been
    initialized so far.

    get_steady_problem: Gets the SteadyProblem at a specified time step.

    initialize_next_problem: Initializes the next time step's SteadyProblem from the
    deformed geometry.

    wing_density: The mass per unit span area of the wing (kg/m^2).

    spring_constant_rad: The torsional spring stiffness for the spring-mass-damper model
    (N*m/rad).

    damping_constant_rad: The torsional damping coefficient (N*m*s/rad).

    step_discards: The number of initial time steps to discard for numerical stability.

    **Notes:**

    The aeroelastic coupling assumes a torsional spring-mass-damper model for each
    spanwise section. Wing motion is prescribed through wing flapping, and each strip's
    aerodynamic moment y components (in the first Airplane's geometry axes, relative to
    the strip's leading edge point) are combined with its inertial and spring restoring
    moments about the same axis and point via ODE integration to produce structural
    deformations.
    """

    __slots__ = (
        "_wing_density",
        "_spring_constant_rad",
        "_damping_constant_rad",
        "_step_discards",
        "listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz",
        "_listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz",
    )

    def __init__(
        self,
        movement: aeroelastic_movement_mod.AeroelasticMovement,
        wing_density: float,
        spring_constant_rad: float,
        damping_constant_rad: float,
        step_discards: int = 5,
    ) -> None:
        """The initialization method.

        Sets up the aeroelastic problem with structural parameters for the torsional
        spring-mass-damper model applied to each wing spanwise section. Initializes the
        per-wing structural state time series.

        See _CoupledUnsteadyProblem's initialization method for descriptions of
        inherited parameters.

        :param movement: An AeroelasticMovement containing the prescribed motion and
            aerodynamic setup for the aeroelastic simulation.
        :param wing_density: The mass per unit span area of the wing (kg/m^2). Used to
            distribute wing mass across Panels for inertial calculations.
        :param spring_constant_rad: The torsional spring stiffness for the spring-mass-
            damper model (N*m/rad). Controls the restoring moment opposing deformation.
        :param damping_constant_rad: The torsional damping coefficient (N*m*s/rad).
            Controls the viscous damping in the spring-mass-damper system.
        :param step_discards: The number of initial time steps to discard for numerical
            stability (there are inconsistent startup effects from the UVLM solver).
            During these time steps, the solver will run but the results will not be
            applied to the deformation of the Wings. The default is 5.
        :return: None
        """
        if not isinstance(movement, aeroelastic_movement_mod.AeroelasticMovement):
            raise TypeError("movement must be an AeroelasticMovement.")

        # Generate the initial Airplane at step 0 with no deformation.
        initial_airplane = movement.generate_airplane_at_time_step(
            airplane_movement_index=0, step=0
        )

        super().__init__(
            movement=movement,
            initial_airplanes=[initial_airplane],
            initial_operating_point=movement.operating_points[0],
        )

        # These are the tunable parameters.
        self._wing_density = wing_density  # This is per unit height, in kg/m^2.
        self._spring_constant_rad = spring_constant_rad
        self._damping_constant_rad = damping_constant_rad

        # These are the permanent parameters.
        self._step_discards = step_discards

        # This is the per-wing time series of the cumulative torsional deformation
        # state, indexed as [wing_idx][entry], where entry 0 is the zero-valued initial
        # state and every generated next step appends one entry, so the last entry is
        # always the state. Each entry is a (num_spanwise_panels + 1,) ndarray whose
        # element n is the y component (radians) of the deformation angle vector that
        # perturbs the corresponding WingCrossSection's angles_Wcsp_to_Wcs_ixyz (angles
        # describing the orientation of the wing cross section axes relative to the wing
        # cross section parent axes using an intrinsic xy'z" sequence). The
        # perturbations' x and z components are structurally zero, so only the y
        # components are stored. The derivative entries are the angle elements' time
        # derivatives (rad/s). They are rates of change of scalar angle components, not
        # angular velocity vector components, and they carry no reference frame ID
        # because differentiating a scalar coordinate involves no rotating basis (see
        # the Angle Component Time Derivatives section of
        # ANGLE_VECTORS_AND_TRANSFORMATIONS.md).
        self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz: list[list[np.ndarray]] = []
        self._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz: list[
            list[np.ndarray]
        ] = []

        # Initialize per-wing state now that we have the initial Airplane's geometry.
        self._initialize_per_wing_state(initial_airplane)

    # --- Immutable: read only properties ---
    @property
    def _aeroelastic_movement(
        self,
    ) -> aeroelastic_movement_mod.AeroelasticMovement:
        # The parent stores the Movement as a CoreMovement in _movement. The constructor
        # guarantees it is an AeroelasticMovement, so the cast here is safe.
        return cast(
            aeroelastic_movement_mod.AeroelasticMovement,
            self._movement,
        )

    @property
    def wing_density(self) -> float:
        return self._wing_density

    @property
    def spring_constant_rad(self) -> float:
        return self._spring_constant_rad

    @property
    def damping_constant_rad(self) -> float:
        return self._damping_constant_rad

    @property
    def step_discards(self) -> int:
        return self._step_discards

    def _initialize_per_wing_state(self, airplane: geometry.airplane.Airplane) -> None:
        """Allocate per-wing time series lists sized to the Airplane's geometry.

        Called once from __init__ after the initial Airplane is generated. Iterates over
        every Wing in the Airplane and appends one entry per Wing to each per-wing time
        series list. listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz and
        _listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz are each seeded with a
        zero-valued (num_spanwise_panels + 1,) initial-state entry, so their last
        entries always hold the current structural state.

        :param airplane: The initial Airplane whose Wings define the geometry.
        :return: None
        """
        for wing in airplane.wings:
            num_spanwise_panels = wing.num_spanwise_panels
            assert num_spanwise_panels is not None
            self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz.append(
                [np.zeros(num_spanwise_panels + 1, dtype=float)]
            )
            self._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz.append(
                [np.zeros(num_spanwise_panels + 1, dtype=float)]
            )

    def _record_null_step_for_wing(self, wing_idx: int) -> None:
        """Append zero-valued entries to the time series lists for a non-aeroelastic
        Wing.

        Called once per time step for Wings backed by a standard WingMovement (no
        deformation). Keeps the time series lists length-consistent with the aeroelastic
        Wings so that indexing by time step is always valid.

        :param wing_idx: Index of the Wing in airplane.wings (and the per-wing lists).
        :return: None
        """
        num_wing_cross_sections = self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[
            wing_idx
        ][-1].shape[0]
        self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx].append(
            np.zeros(num_wing_cross_sections, dtype=float)
        )
        self._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx].append(
            np.zeros(num_wing_cross_sections, dtype=float)
        )

    def _calculate_mass_matrix(self, wing: geometry.wing.Wing) -> np.ndarray:
        """Generate the mass distribution matrix for all of a Wing's Panels.

        Distributes the total spanwise mass (wing_density) across Panel areas to form a
        panel-by-panel mass matrix. Each Panel's mass is proportional to its area times
        the specified wing_density.

        :param wing: The Wing whose Panels define the mass distribution.
        :return: A (num_chordwise_panels, num_spanwise_panels, 3) ndarray of floats
            representing the mass at each Panel. The three components are identical
            (mass scalar replicated for x, y, z axes).
        """
        assert wing.panels is not None
        areas = np.array([[panel.area for panel in row] for row in wing.panels])
        return np.repeat(areas[:, :, None], 3, axis=2) * self.wing_density

    def initialize_next_problem(
        self,
        solver: CoupledUnsteadyRingVortexLatticeMethodSolver,
        step: int,
    ) -> None:
        # The solver invokes this on every step. Aeroelasticity only builds the next
        # step's deformed geometry, so there is nothing to do on the final step.
        if step >= self.num_steps - 1:
            return

        aeroelastic_solver = cast(
            "AeroelasticUnsteadyRingVortexLatticeMethodSolver", solver
        )

        next_step = len(self._steady_problems)

        # _calculate_wing_deformation returns a per-wing list: Each element is either
        # the deformation ndarray for an aeroelastic Wing or None for a standard Wing.
        deformationAnglesRad_Wcsp_to_Wcs_ixyz = self._calculate_wing_deformation(
            aeroelastic_solver, next_step
        )

        # The structural model works in radians. The geometry API expects degrees.
        deformationAngles_Wcsp_to_Wcs_ixyz = [
            np.rad2deg(arr) if arr is not None else None
            for arr in deformationAnglesRad_Wcsp_to_Wcs_ixyz
        ]

        # Generate the deformed Airplane at this step.
        airplane = self._aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0,
            step=next_step,
            deformationAngles_Wcsp_to_Wcs_ixyz=deformationAngles_Wcsp_to_Wcs_ixyz,
        )
        operating_point = self._aeroelastic_movement.operating_points[next_step]

        self._steady_problems.append(
            SteadyProblem(
                airplanes=[airplane],
                operating_point=operating_point,
            )
        )

    def _calculate_wing_deformation(
        self,
        solver: AeroelasticUnsteadyRingVortexLatticeMethodSolver,
        step: int,
    ) -> list[np.ndarray | None]:
        """Compute cumulative wing deformation for the current time step.

        Loops over every Wing in the current Airplane. For Wings backed by an
        AeroelasticWingMovement, orchestrates the aero moment extraction and spring ODE
        integration, records the per-wing structural state, and returns the deformation
        array. For Wings backed by a standard WingMovement, records null entries in the
        time series lists and returns None for that Wing.

        :param solver: The solver instance providing aerodynamic moment data
            (moments_GP1_Slep).
        :param step: The current time step index (0-indexed).
        :return: A list of length len(airplane.wings) where each element is either a
            (num_spanwise_panels + 1, 3) ndarray of cumulative deformation angles for an
            aeroelastic Wing or None for a non-aeroelastic Wing. Each row is a (3,)
            angle vector that perturbs the corresponding WingCrossSection's
            angles_Wcsp_to_Wcs_ixyz (angles describing the orientation of the wing cross
            section axes relative to the wing cross section parent axes using an
            intrinsic xy'z" sequence). The units are in radians.
        """
        curr_problem: SteadyProblem = self._steady_problems[-1]
        airplane = curr_problem.airplanes[0]
        wing_movements = self._aeroelastic_movement.airplane_movements[0].wing_movements

        results: list[np.ndarray | None] = []
        panel_offset = 0

        for wing_idx, wing in enumerate(airplane.wings):
            wing_movement = wing_movements[wing_idx]
            num_chordwise_panels = wing.num_chordwise_panels
            num_spanwise_panels = wing.num_spanwise_panels
            assert num_spanwise_panels is not None
            num_panels = num_chordwise_panels * num_spanwise_panels

            if isinstance(
                wing_movement,
                aeroelastic_wing_movement_mod.AeroelasticWingMovement,
            ):
                mass_matrix = self._calculate_mass_matrix(wing)

                aeroMoments_GP1_Slep = self._extract_aero_moments(
                    solver,
                    num_chordwise_panels,
                    num_spanwise_panels,
                    num_panels,
                    panel_offset,
                )

                # A mirror-meshed Wing's Panel grid runs tip to root spanwise, while the
                # structural solve and the WingCrossSectionMovements that consume its
                # output run root to tip, so flip the spanwise axis of the per-Panel
                # arrays to put them in root-to-tip strip order.
                if wing.mirror_only:
                    mass_matrix = mass_matrix[:, ::-1, :]
                    aeroMoments_GP1_Slep = aeroMoments_GP1_Slep[:, ::-1, :]

                (
                    newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz,
                    newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz,
                ) = self._calculate_spring_moments(
                    num_spanwise_panels=num_spanwise_panels,
                    wing=wing,
                    mass_matrix=mass_matrix,
                    aeroMoments_GP1_Slep=aeroMoments_GP1_Slep,
                    step=step,
                    wing_idx=wing_idx,
                    wing_movement=wing_movement,
                )

                self._record_structural_state(
                    step=step,
                    newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz=(
                        newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz
                    ),
                    newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz=(
                        newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz
                    ),
                    wing_idx=wing_idx,
                )

                # Expand the newly recorded deformation angle y components into the
                # (num_spanwise_panels + 1, 3) deformation angle vectors that perturb
                # the WingCrossSections' angles_Wcsp_to_Wcs_ixyz, which the movement
                # chain consumes.
                results.append(
                    self._build_deformation_vector(
                        self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx][-1],
                        num_spanwise_panels,
                    )
                )
            else:
                self._record_null_step_for_wing(wing_idx)
                results.append(None)

            panel_offset += num_panels

        return results

    @staticmethod
    def _extract_aero_moments(
        solver: AeroelasticUnsteadyRingVortexLatticeMethodSolver,
        num_chordwise_panels: int,
        num_spanwise_panels: int,
        num_panels: int,
        panel_offset: int,
    ) -> np.ndarray:
        """Extract the aerodynamic moments (in the first Airplane's geometry axes, each
        relative to its strip's leading edge point) from the solver output.

        Uses the strip leading edge points as the moments' reference points, consistent
        with the assumption of a torsional spring at the leading edge. Each strip's
        leading edge point is the leading edge point of its outboard bounding
        WingCrossSection, the same WingCrossSection that receives the strip's torsional
        deformation.

        :param solver: The solver instance with moments_GP1_Slep data.
        :param num_chordwise_panels: Number of chordwise panel rows.
        :param num_spanwise_panels: Number of spanwise panel rows.
        :param num_panels: Total number of panels (num_chordwise * num_spanwise).
        :param panel_offset: The flat panel index offset into solver.moments_GP1_Slep at
            which this Wing's data begins.
        :return: A (num_chordwise_panels, num_spanwise_panels, 3) ndarray of aerodynamic
            moments (in the first Airplane's geometry axes, each relative to its strip's
            leading edge point).
        """
        aeroMoments_GP1_Slep = np.array(
            solver.moments_GP1_Slep[panel_offset : panel_offset + num_panels]
        ).reshape(num_chordwise_panels, num_spanwise_panels, 3)
        return aeroMoments_GP1_Slep

    @staticmethod
    def _build_deformation_vector(
        deformationAnglesYRad_Wcsp_to_Wcs_ixyz: np.ndarray,
        num_spanwise_panels: int,
    ) -> np.ndarray:
        """Construct the deformation angle vectors from scalar torsional angles.

        Expands the stored torsional angle state (one scalar per WingCrossSection) into
        the full (num_spanwise_panels + 1, 3) array of deformation angle vectors that
        perturb the WingCrossSections' angles_Wcsp_to_Wcs_ixyz, using an intrinsic xy'z"
        sequence, which the movement chain consumes.

        :param deformationAnglesYRad_Wcsp_to_Wcs_ixyz: A (num_spanwise_panels + 1,)
            ndarray holding, for each WingCrossSection, the y component (radians) of the
            cumulative deformation angle vector that perturbs its
            angles_Wcsp_to_Wcs_ixyz, using an intrinsic xy'z" sequence.
        :param num_spanwise_panels: Number of spanwise panel rows.
        :return: A (num_spanwise_panels + 1, 3) ndarray of the deformation angle vectors
            (radians) that perturb the WingCrossSections' angles_Wcsp_to_Wcs_ixyz, with
            zero-valued x and z components and the given torsional angles as the y
            components.
        """
        deformationAnglesRad_Wcsp_to_Wcs_ixyz = np.array(
            [
                np.array(
                    [
                        0,
                        deformationAnglesYRad_Wcsp_to_Wcs_ixyz[i + 1],
                        0,
                    ]
                )
                for i in range(num_spanwise_panels)
            ]
        )
        deformationAnglesRad_Wcsp_to_Wcs_ixyz = np.insert(
            deformationAnglesRad_Wcsp_to_Wcs_ixyz, 0, np.array([0, 0, 0]), axis=0
        )
        return deformationAnglesRad_Wcsp_to_Wcs_ixyz

    def _record_structural_state(
        self,
        step: int,
        newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz: np.ndarray,
        newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz: np.ndarray,
        wing_idx: int,
    ) -> None:
        """Append one Wing's new structural state to its time series.

        The recorded state is the y components (radians) of the cumulative deformation
        angle vectors that perturb the Wing's WingCrossSections' angles_Wcsp_to_Wcs_ixyz
        (angles describing the orientation of each wing cross section's axes relative to
        its parent axes using an intrinsic xy'z" sequence), together with those y
        components' time derivatives (rad/s). Both are recorded together, and both are
        held at their previous state during the early unstable time steps.

        :param step: The current time step index.
        :param newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz: A (num_spanwise_panels + 1,)
            ndarray of this time step's new values (radians) of the y components of the
            deformation angle vectors that perturb the WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz, with one element per WingCrossSection, in root-to-
            tip order.
        :param newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz: A
            (num_spanwise_panels + 1,) ndarray of the new time derivatives (rad/s) of
            the y components of the deformation angle vectors that perturb the
            WingCrossSections' angles_Wcsp_to_Wcs_ixyz, with one element per
            WingCrossSection, in root-to-tip order.
        :param wing_idx: Index of the Wing in airplane.wings (and the per-wing lists).
        :return: None
        """
        # Append the cumulative torsional angles and their time derivatives together,
        # holding both at their previous state during the first step_discards steps,
        # whose numerical startup effects cause large aerodynamic forces. The strips'
        # ODEs are re-seeded from the time series' last entries, so recording only one
        # component would seed a mixed state. Each recorded entry is a fresh array, so
        # no time series entry aliases a later one.
        deformationAnglesYRad_Wcsp_to_Wcs_ixyz = (
            self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = (
            self._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        if step > self.step_discards:
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz.append(
                np.array(newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz, dtype=float)
            )
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz.append(
                np.array(
                    newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz, dtype=float
                )
            )
        else:
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz.append(
                deformationAnglesYRad_Wcsp_to_Wcs_ixyz[-1].copy()
            )
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz.append(
                deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[-1].copy()
            )

    def _calculate_spring_moments(
        self,
        num_spanwise_panels: int,
        wing: geometry.wing.Wing,
        mass_matrix: np.ndarray,
        aeroMoments_GP1_Slep: np.ndarray,
        step: int,
        wing_idx: int,
        wing_movement: aeroelastic_wing_movement_mod.AeroelasticWingMovement,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve the torsional spring-damper ODE for each spanwise section.

        Solves the torsional spring-damper ODE independently for each spanwise section,
        accounting for aerodynamic moments, inertial moments, and structural properties.
        The state integrated for each strip is the y component (radians) of the
        cumulative deformation angle vector that perturbs its outboard bounding
        WingCrossSection's angles_Wcsp_to_Wcs_ixyz, using an intrinsic xy'z" sequence,
        together with that y component's time derivative (rad/s). Each strip's two
        rotational inertias (the torsional inertia in its ODE and the flapping-axis
        inertia that scales its prescribed-motion inertial forcing) come from
        _calculate_strip_inertias.

        The spanwise section index runs root to tip, matching the order of the Wing's
        wing_cross_sections, and the mass_matrix and aeroMoments_GP1_Slep arrays must be
        supplied with their spanwise axes in that same root-to-tip order (the caller
        flips them for mirror-meshed Wings, whose panel grids run tip to root).

        :param num_spanwise_panels: Number of spanwise panel rows in the Wing.
        :param wing: The Wing containing geometric and structural definitions.
        :param mass_matrix: A (num_chordwise_panels, num_spanwise_panels, 3) ndarray of
            panel masses, with the spanwise axis in root-to-tip order.
        :param aeroMoments_GP1_Slep: A (num_chordwise_panels, num_spanwise_panels, 3)
            ndarray of aerodynamic moments from the aerodynamic solver (in the first
            Airplane's geometry axes, each relative to its strip's leading edge point),
            with the spanwise axis in root-to-tip order.
        :param step: The current time step index.
        :param wing_idx: Index of the Wing in airplane.wings.
        :param wing_movement: The AeroelasticWingMovement providing the prescribed
            flapping parameters used for inertial moment generation.
        :return: A tuple of two (num_spanwise_panels + 1,) ndarrays: the new values
            (radians) of the y components of the deformation angle vectors that perturb
            the WingCrossSections' angles_Wcsp_to_Wcs_ixyz, and those y components' new
            time derivatives (rad/s), each with one element per WingCrossSection, in
            root-to-tip order.
        """
        newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz = np.zeros(
            num_spanwise_panels + 1, dtype=float
        )
        newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = np.zeros(
            num_spanwise_panels + 1, dtype=float
        )
        # The current structural state is the last entry of each time series list. Both
        # reads happen before this step's entries are appended, so they see the state at
        # the end of the previous time step.
        lastDeformationAnglesYRad_Wcsp_to_Wcs_ixyz = (
            self.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx][-1]
        )
        lastDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = (
            self._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx][-1]
        )
        dt = self.movement.delta_time
        torsional_inertias, flapping_axis_inertias = self._calculate_strip_inertias(
            wing=wing,
            num_spanwise_panels=num_spanwise_panels,
            mass_matrix=mass_matrix,
        )
        for span_panel in range(num_spanwise_panels):
            aeroMomentY_GP1_Slep = np.sum(aeroMoments_GP1_Slep[:, span_panel, 1])
            # Re-seed the strip's ODE from its own state at the end of the previous time
            # step. The state arrays hold one entry per WingCrossSection, and a strip's
            # state lives at its outboard bounding WingCrossSection, index span_panel +
            # 1. The WingCrossSection at index span_panel bounds this strip on its
            # inboard side, and the root WingCrossSection (index 0) is clamped.
            (
                newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[span_panel + 1],
                newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[span_panel + 1],
            ) = self._calculate_torsional_spring_moment(
                dt,
                torsional_inertia=torsional_inertias[span_panel],
                lastDeformationAngleYRad_Wcsp_to_Wcs_ixyz=(
                    lastDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[span_panel + 1]
                ),
                lastDeformationAngleDerivativeYRad_Wcsp_to_Wcs_ixyz=(
                    lastDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[span_panel + 1]
                ),
                aeroMomentY_GP1_Slep=aeroMomentY_GP1_Slep,
                step=step,
                flapping_axis_inertia=flapping_axis_inertias[span_panel],
                wing_movement=wing_movement,
            )

        return (
            newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz,
            newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz,
        )

    @staticmethod
    def _calculate_strip_inertias(
        wing: geometry.wing.Wing,
        num_spanwise_panels: int,
        mass_matrix: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the two rotational inertias of each of a Wing's strips.

        The structural model uses two distinct inertias per strip. The torsional
        inertia, (1 / 2) * M * L^2, is the inertia paired with the strip's torsional
        deformation state: it divides the net moment in the strip's spring-damper ODE.
        The flapping-axis inertia, (1 / 12) * M * (L^2 + W^2) + M * d^2, is the inertia
        of the strip (modeled as a rectangular prism) about the flapping axis, with the
        M * d^2 term applying the parallel axis theorem: it never appears in the ODE
        directly, and instead scales the prescribed flapping acceleration into the ODE's
        inertial forcing moment. In both expressions, M is the strip's mass, L is its
        mean chord, W is its span width, and d is the distance from the flapping axis to
        its centroid, accumulated across strips in half-span-width increments.

        The strip index runs root to tip, matching the order of the Wing's
        wing_cross_sections, and the mass_matrix must be supplied with its spanwise axis
        in that same root-to-tip order (the panel geometry reads flip the index for
        mirror-meshed Wings, whose panel grids run tip to root).

        :param wing: The Wing containing geometric and structural definitions.
        :param num_spanwise_panels: Number of spanwise panel rows in the Wing.
        :param mass_matrix: A (num_chordwise_panels, num_spanwise_panels, 3) ndarray of
            panel masses, with the spanwise axis in root-to-tip order.
        :return: A tuple of two (num_spanwise_panels,) ndarrays of floats: the strips'
            torsional inertias and their flapping-axis inertias (both kg*m^2), in root-
            to-tip order.
        """
        torsional_inertias = np.zeros(num_spanwise_panels, dtype=float)
        flapping_axis_inertias = np.zeros(num_spanwise_panels, dtype=float)
        assert wing.panels is not None
        # This is the distance from the flapping axis to the strip's centroid,
        # accumulated in half-span-width increments.
        d = 0.0
        for span_panel in range(num_spanwise_panels):
            mass = mass_matrix[:, span_panel, :].sum()
            L = (
                wing.wing_cross_sections[span_panel].chord
                + wing.wing_cross_sections[span_panel + 1].chord
            ) / 2
            # The span_panel index runs root to tip, matching the wing_cross_sections
            # list, but a mirror-meshed Wing's Panel grid runs tip to root spanwise, so
            # map the index when reading Panel geometry.
            if wing.mirror_only:
                panel_span_index = num_spanwise_panels - 1 - span_panel
            else:
                panel_span_index = span_panel
            W: float = float(
                np.linalg.norm(wing.panels[0][panel_span_index].frontLeg_G)
            )
            d += W / 2
            torsional_inertias[span_panel] = 1 / 2 * mass * (L**2)
            flapping_axis_inertias[span_panel] = 1 / 12 * mass * (
                L**2 + W**2
            ) + mass * (d**2)
            d += W / 2
        return torsional_inertias, flapping_axis_inertias

    def _calculate_torsional_spring_moment(
        self,
        dt: float,
        torsional_inertia: float,
        lastDeformationAngleYRad_Wcsp_to_Wcs_ixyz: float,
        lastDeformationAngleDerivativeYRad_Wcsp_to_Wcs_ixyz: float,
        aeroMomentY_GP1_Slep: float,
        step: int,
        flapping_axis_inertia: float,
        wing_movement: aeroelastic_wing_movement_mod.AeroelasticWingMovement,
        num_steps: int = 2,
    ) -> tuple[float, float]:
        """Solve the torsional spring-damper ODE for a single strip.

        Integrates the forced torsional damped harmonic oscillator equation:
        I * d^2(theta)/dt^2 = M_aero + M_inertial - k * theta - c * d(theta)/dt, where
        I is the strip's torsional inertia and k and c are the problem's
        spring_constant_rad and damping_constant_rad.

        The state integrated is the y component (radians) of the cumulative deformation
        angle vector that perturbs the strip's outboard bounding WingCrossSection's
        angles_Wcsp_to_Wcs_ixyz, using an intrinsic xy'z" sequence, together with that
        y component's time derivative (rad/s). Returns their new values at the end of
        the time step.

        :param dt: The time step duration (seconds).
        :param torsional_inertia: The strip's torsional rotational inertia (kg*m^2),
            the I dividing the net moment in the ODE above.
        :param lastDeformationAngleYRad_Wcsp_to_Wcs_ixyz: The value (radians), at the
            end of the previous time step, of the y component of the deformation angle
            vector that perturbs the strip's outboard bounding WingCrossSection's
            angles_Wcsp_to_Wcs_ixyz. It seeds the integration.
        :param lastDeformationAngleDerivativeYRad_Wcsp_to_Wcs_ixyz: The time derivative
            (rad/s), at the end of the previous time step, of the y component of the
            deformation angle vector that perturbs the strip's outboard bounding
            WingCrossSection's angles_Wcsp_to_Wcs_ixyz. It seeds the integration.
        :param aeroMomentY_GP1_Slep: The strip's torsional aerodynamic forcing (N*m):
            the sum, over its chordwise Panels, of the y components of their
            aerodynamic moments (in the first Airplane's geometry axes, relative to the
            strip's leading edge point).
        :param step: The current time step index (used for inertial moment evaluation).
        :param flapping_axis_inertia: The strip's rotational inertia about the flapping
            axis (kg*m^2), including the parallel axis theorem term. It never divides
            the net moment in the ODE above; it only scales the prescribed flapping
            acceleration into the inertial forcing moment M_inertial.
        :param wing_movement: The AeroelasticWingMovement whose prescribed flapping
            parameters are used.
        :param num_steps: Number of time sub-steps for numerical integration. The
            default is 2.
        :return: A tuple of two floats: the new value (radians) of the y component of
            the deformation angle vector that perturbs the strip's outboard bounding
            WingCrossSection's angles_Wcsp_to_Wcs_ixyz, and that y component's new time
            derivative (rad/s), both at the end of the time step.
        """
        t = np.linspace(dt * (step - 1), dt * step, num_steps)

        # Perform the forced numerical integration of the spring-damper ODE.
        return self._spring_numerical_ode(
            t,
            self.spring_constant_rad,
            self.damping_constant_rad,
            torsional_inertia,
            lastDeformationAngleYRad_Wcsp_to_Wcs_ixyz,
            lastDeformationAngleDerivativeYRad_Wcsp_to_Wcs_ixyz,
            aeroMomentY_GP1_Slep,
            self._generate_inertial_moment_function(
                flapping_axis_inertia, wing_movement=wing_movement
            ),
        )

    @staticmethod
    def _generate_inertial_moment_function(
        flapping_axis_inertia: float,
        wing_movement: aeroelastic_wing_movement_mod.AeroelasticWingMovement,
    ) -> Callable[[float], float]:
        """Generate the prescribed wing motion inertial moment function.

        Extracts the prescribed flapping motion from the wing_movement definition and
        creates a callable inertial moment function M_inertial = I *
        d^2(theta_prescribed)/dt^2, where I is the strip's flapping-axis inertia.
        Supports sinusoidal and custom spacing functions.

        :param flapping_axis_inertia: The strip's rotational inertia about the
            flapping axis (kg*m^2).
        :param wing_movement: The AeroelasticWingMovement whose prescribed flapping
            parameters are used.
        :return: A callable function that accepts time and returns the inertial moment
            (N*m) due to the prescribed wing motion acceleration. For a zero flapping
            amplitude (no prescribed flapping), the returned function is identically
            zero regardless of the spacing. For sinusoidal spacing: M_inertial = -I * b^2 *
            sin(b * t + h) * A, where I = flapping_axis_inertia, b = 2 * pi / period,
            h = phase, A = amplitude. For custom spacing, uses the
            AeroelasticWingMovement's spacingAnglesSecondDerivative_Gs_to_Wn_ixyz, which
            its constructor guarantees is present whenever the spacing is a custom
            callable.
        """
        amp = wing_movement.ampAngles_Gs_to_Wn_ixyz[0]

        # A wing with no prescribed flapping applies no inertial moment. Return the zero
        # function before computing the flapping frequency, whose 2 * pi / period
        # expression is undefined for the zero period that a zero amplitude requires.
        if amp == 0.0:
            return lambda time: 0.0

        b_rad = 2 * np.pi / wing_movement.periodAngles_Gs_to_Wn_ixyz[0]
        h_rad = np.deg2rad(wing_movement.phaseAngles_Gs_to_Wn_ixyz[0])
        spacing = wing_movement.spacingAngles_Gs_to_Wn_ixyz[0]
        if spacing == "sine":
            # Since amp is in degrees (ampAngles_Gs_to_Wn_ixyz), convert to radians so
            # the inertial moment (N*m) is consistent with the SI spring-damper ODE.
            moment_func = (
                lambda time: -1
                * (b_rad**2)
                * np.sin(b_rad * time + h_rad)
                * np.deg2rad(amp)
                * flapping_axis_inertia
            )
        elif spacing == "uniform":
            raise ValueError(
                "Sawtooth function (uniform spacing) is not differentiable, "
                "cannot be used for inertial moment function."
            )
        elif callable(spacing):
            # The AeroelasticWingMovement's constructor guarantees a matching second
            # derivative whenever the spacing component is a custom callable.
            deriv = wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz[0]
            assert deriv is not None
            # deriv is the second derivative before amplitude scaling, and amp is in
            # degrees, so convert to radians to keep the moment in N*m.
            moment_func = (
                lambda time: np.deg2rad(amp) * deriv(time) * flapping_axis_inertia
            )

        return moment_func

    @staticmethod
    def _spring_numerical_ode(
        t: np.ndarray,
        spring_constant_rad: float,
        damping_constant_rad: float,
        torsional_inertia: float,
        theta0_rad: float,
        theta_derivative0_rad: float,
        aero_moment: float,
        inertial_moment_func: Callable[[float], float],
    ) -> tuple[float, float]:
        """Numerically integrate the torsional spring-damper ODE.

        Solves the second-order forced ODE: I * d^2(theta)/dt^2 = M_aero +
        M_inertial(t) - k * theta - c * d(theta)/dt

        using scipy.integrate.solve_ivp with tolerances a few orders of magnitude
        stricter than scipy's defaults. This integration is re-seeded from the previous
        state at every outer time step, so its local errors compound across a
        simulation rather than being controlled within a single integration. The
        stricter tolerances guard against that accumulation, and against the absolute
        tolerance floor swamping small torsional responses, at a cost that is small
        relative to the aerodynamic solve.

        :param t: A (N,) ndarray of time points for integration evaluation.
        :param spring_constant_rad: The spring constant (N*m/rad), the k in the ODE
            above.
        :param damping_constant_rad: The damping constant (N*m*s/rad), the c in the ODE
            above.
        :param torsional_inertia: The strip's torsional rotational inertia (kg*m^2),
            the I in the ODE above.
        :param theta0_rad: Initial angular displacement (radians).
        :param theta_derivative0_rad: The initial angular displacement's time
            derivative (rad/s).
        :param aero_moment: Constant aerodynamic moment acting on the section (N*m).
        :param inertial_moment_func: A callable function of time that returns the
            inertial moment from prescribed motion acceleration (N*m).
        :return: A tuple of two floats: the final angle (radians) and its time
            derivative (rad/s) at the last time point in t.
        """

        def external_moment(time: float) -> float:
            """Total external moment M_external (aerodynamic + inertial from prescribed
            motion)."""
            return float(aero_moment + inertial_moment_func(time))

        def ode(time: float, y: np.ndarray) -> np.ndarray:
            """ODE system: the state is (theta, d(theta)/dt), and d^2(theta)/dt^2 =
            (M_external - c * d(theta)/dt - k * theta) / I."""
            theta_rad, theta_derivative_rad = y
            return np.array(
                [
                    theta_derivative_rad,
                    (
                        external_moment(time)
                        - damping_constant_rad * theta_derivative_rad
                        - spring_constant_rad * theta_rad
                    )
                    / torsional_inertia,
                ]
            )

        sol = solve_ivp(
            ode,
            (t[0], t[-1]),
            np.array([theta0_rad, theta_derivative0_rad], dtype=float),
            t_eval=t,
            rtol=_SPRING_ODE_RELATIVE_TOLERANCE,
            atol=_SPRING_ODE_ABSOLUTE_TOLERANCE_RAD,
        )

        final_theta_rad = float(sol.y[0][-1])
        final_theta_derivative_rad = float(sol.y[1][-1])

        return final_theta_rad, final_theta_derivative_rad


def _get_mujoco_model(
    problem: FreeFlightUnsteadyProblem,
) -> _mujoco_model.MuJoCoModel:
    """Returns a FreeFlightUnsteadyProblem's MuJoCoModel.

    Defined here so the read of the FreeFlightUnsteadyProblem's private slot stays
    inside the module that owns it. Registering it with _private_access at import time
    lets the rendering layer reach the MuJoCoModel without any cross-module private
    access.

    :param problem: The FreeFlightUnsteadyProblem whose MuJoCoModel will be returned.
    :return: The FreeFlightUnsteadyProblem's MuJoCoModel.
    """
    return problem._mujoco_model


# Register the getter at import time so the rendering layer can look up a
# FreeFlightUnsteadyProblem's MuJoCoModel through _private_access.
_private_access.register_mujoco_model_getter(_get_mujoco_model)
