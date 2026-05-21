"""Contains the Movement class.

**Contains the following classes:**

Movement: A class used to contain an UnsteadyProblem's movement.

**Contains the following functions:**

None
"""

from __future__ import annotations

import copy
import math
from typing import cast

import numpy as np
import scipy.optimize as sp_opt

from .. import _core, _logging, _parameter_validation, geometry
from .. import operating_point as operating_point_mod
from .. import problems
from . import airplane_movement as airplane_movement_mod
from . import operating_point_movement as operating_point_movement_mod

movement_logger = _logging.get_logger("movements.movement")


class Movement(_core.CoreMovement):
    """A class used to contain an UnsteadyProblem's movement.

    **Contains the following methods:**

    lcm_period: The least common multiple of all motion periods, ensuring all motions
    complete an integer number of cycles when cycle averaging forces and moments.

    max_period: The longest period of motion of Movement's sub movement objects, the
    motion(s) of its sub sub movement object(s), and the motions of its sub sub sub
    movement objects.

    min_period: The shortest non zero period of motion of Movement's sub movement
    objects, the motion(s) of its sub sub movement object(s), and the motions of its sub
    sub sub movement objects.

    static: Flags if Movement's sub movement objects, its sub sub movement object(s),
    and its sub sub sub movement objects all represent no motion.
    """

    __slots__ = (
        "_num_cycles",
        "_num_chords",
        "_max_wake_chords",
        "_max_wake_cycles",
        "_airplanes",
        "_operating_points",
    )

    def __init__(
        self,
        airplane_movements: list[airplane_movement_mod.AirplaneMovement],
        operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
        delta_time: float | int | str | None = None,
        num_cycles: int | None = None,
        num_chords: int | None = None,
        num_steps: int | None = None,
        max_wake_rows: int | None = None,
        max_wake_chords: int | None = None,
        max_wake_cycles: int | None = None,
    ) -> None:
        """The initialization method.

        This method checks that all Wings maintain their symmetry type across all time
        steps. See the WingMovement class documentation for more details on this
        requirement. See the Wing class documentation for more information on symmetry
        types.

        :param airplane_movements: A list of the AirplaneMovements associated with each
            of the UnsteadyProblem's Airplanes.
        :param operating_point_movement: An OperatingPointMovement characterizing any
            changes to the UnsteadyProblem's operating conditions.
        :param delta_time: The time between each time step. Accepts the following: None
            (default): Movement analytically estimates the delta_time that produces wake
            ring vortices with roughly the same chord length as the bound trailing edge
            ring vortices, accounting for both freestream and geometry motion
            velocities. This provides good results across all Strouhal numbers.
            "optimize": Movement first runs the analytical estimation, then uses that
            result as an initial guess for an iterative optimization that minimizes the
            area mismatch between wake ring vortices and their parent bound trailing
            edge ring vortices. This is slower but may produce slightly more accurate
            results. Positive number (int or float): Use the specified value directly.
            All values are converted internally to floats. The units are in seconds.
        :param num_cycles: The number of cycles of the maximum period motion used to
            calculate a num_steps parameter initialized as None if Movement isn't
            static. If num_steps is not None or if Movement is static, this must be
            None. If num_steps is initialized as None and the Movement isn't static,
            num_cycles must be a positive int. In that case, I recommend setting
            num_cycles to 3. The default is None.
        :param num_chords: The number of chord lengths used to calculate a num_steps
            parameter initialized as None if Movement is static. If num_steps is not
            None or if Movement isn't static, this must be None. If num_steps is
            initialized as None and Movement is static, num_chords must be a positive
            int. In that case, I recommend setting num_chords to 10. For cases with
            multiple Airplanes, the num_chords will reference the largest reference
            chord length. The default is None.
        :param num_steps: The number of time steps of the unsteady simulation. If
            initialized as None, and Movement isn't static, Movement will calculate a
            value for num_steps such that the simulation will cover some number of
            cycles of the maximum period of all the motion described in Movement's sub
            movement objects, sub sub movement object(s), and sub sub sub movement
            objects. If num_steps is initialized as None, and Movement is static,
            Movement will calculate a value for num_steps such that the simulation will
            result in a wake extending back by some number of reference chord lengths.
        :param max_wake_rows: The maximum number of chordwise wake ring vortex rows per
            Wing. Works for both static and non static Movements. Must be a positive int
            if set. At most one of max_wake_rows, max_wake_chords, and max_wake_cycles
            can be non None. The default is None (no truncation).
        :param max_wake_chords: The maximum wake length in chord lengths, converted to
            max_wake_rows internally. Only valid for static Movements (mirrors
            num_chords). Must be a positive int if set. The default is None.
        :param max_wake_cycles: The maximum wake length in motion cycles, converted to
            max_wake_rows internally. Only valid for non static Movements (mirrors
            num_cycles). Must be a positive int if set. The default is None.
        :return: None
        """
        # --- Validate types early ---
        # CoreMovement.__init__() validates these, but Movement needs to use
        # the parameters before calling super().__init__() (for delta_time
        # estimation and period computation). Validate here so that bad types
        # produce clear TypeErrors rather than AttributeErrors.
        if not isinstance(airplane_movements, list):
            raise TypeError("airplane_movements must be a list.")
        if len(airplane_movements) < 1:
            raise ValueError("airplane_movements must have at least one element.")
        for airplane_movement in airplane_movements:
            if not isinstance(
                airplane_movement, airplane_movement_mod.AirplaneMovement
            ):
                raise TypeError(
                    "Every element in airplane_movements must be an "
                    "AirplaneMovement."
                )
        if not isinstance(
            operating_point_movement,
            operating_point_movement_mod.OperatingPointMovement,
        ):
            raise TypeError(
                "operating_point_movement must be an OperatingPointMovement."
            )

        # --- Resolve delta_time ---

        # Track whether iterative optimization should run after analytical.
        _should_iteratively_optimize_delta_time: bool = False

        if isinstance(delta_time, str):
            if delta_time == "optimize":
                _should_iteratively_optimize_delta_time = True
                # Fall through to calculate analytical estimate first.
                delta_time = None
            else:
                raise ValueError('delta_time string must be "optimize".')

        if delta_time is not None:
            delta_time = _parameter_validation.number_in_range_return_float(
                delta_time, "delta_time", min_val=0.0, min_inclusive=False
            )
        else:
            # Calculate a fast initial delta_time estimate based on freestream
            # velocity. This is used as a fallback for static movements and as
            # a starting point for the analytical optimization.
            delta_times = []
            for airplane_movement in airplane_movements:
                # TODO: Consider making this also average across each Airplane's Wings.
                c_ref = airplane_movement.base_airplane.c_ref
                assert c_ref is not None
                delta_times.append(
                    c_ref
                    / airplane_movement.base_airplane.wings[0].num_chordwise_panels
                    / operating_point_movement.base_operating_point.vCg__E
                )
            fast_estimate = sum(delta_times) / len(delta_times)

            # Run analytical optimization to get a better delta_time that accounts
            # for both freestream and geometry motion velocities.
            delta_time = _analytically_optimize_delta_time(
                airplane_movements=list(airplane_movements),
                operating_point_movement=operating_point_movement,
                initial_delta_time=fast_estimate,
            )

        # Run iterative optimization if requested, using the analytical result
        # as the initial guess.
        if _should_iteratively_optimize_delta_time:
            delta_time = _optimize_delta_time(
                airplane_movements=list(airplane_movements),
                operating_point_movement=operating_point_movement,
                initial_delta_time=delta_time,
            )

        # --- Compute period properties locally ---
        # These are needed for num_steps and max_wake calculations before
        # super().__init__() is called. The same values will be lazily cached
        # by CoreMovement when accessed via properties.
        _airplane_movement_max_periods = []
        for airplane_movement in airplane_movements:
            _airplane_movement_max_periods.append(airplane_movement.max_period)
        _max_period = max(
            max(_airplane_movement_max_periods),
            operating_point_movement.max_period,
        )
        _static = _max_period == 0

        _lcm_period: float = 0.0
        if not _static:
            _all_periods: list[float] = []
            for airplane_movement in airplane_movements:
                _all_periods.extend(airplane_movement.all_periods)
            _all_periods.append(operating_point_movement.max_period)
            _lcm_period = _core.lcm_multiple(_all_periods)

        # --- Resolve num_steps ---
        if num_steps is None:
            if _static:
                if num_cycles is not None:
                    raise ValueError(
                        "If num_steps is None and the Movement is static, num_cycles "
                        "must be left as None."
                    )
            else:
                if num_cycles is None:
                    raise ValueError(
                        "If num_steps is None and the Movement isn't static, "
                        "num_cycles must be set."
                    )
        if num_cycles is not None:
            num_cycles = _parameter_validation.int_in_range_return_int(
                num_cycles,
                "num_cycles",
                min_val=1,
                min_inclusive=True,
            )

        if num_steps is None:
            if _static:
                if num_chords is None:
                    raise ValueError(
                        "If num_steps is None and the Movement is static, num_chords "
                        "must be set."
                    )
            else:
                if num_chords is not None:
                    raise ValueError(
                        "If num_steps is None and the Movement isn't static, "
                        "num_chords must be left as None."
                    )
        if num_chords is not None:
            num_chords = _parameter_validation.int_in_range_return_int(
                num_chords,
                "num_chords",
                min_val=1,
                min_inclusive=True,
            )

        if num_cycles is not None or num_chords is not None:
            if num_steps is not None:
                raise ValueError(
                    "If either num_cycles or num_chords is not None, num_steps must "
                    "be None."
                )
        if num_steps is not None:
            num_steps = _parameter_validation.int_in_range_return_int(
                num_steps,
                "num_steps",
                min_val=1,
                min_inclusive=True,
            )
        else:
            if _static:
                # Find the value of the largest reference chord length of all the
                # base Airplanes.
                c_refs = []
                for airplane_movement in airplane_movements:
                    c_ref = airplane_movement.base_airplane.c_ref
                    assert c_ref is not None
                    c_refs.append(c_ref)
                max_c_ref = max(c_refs)

                # Set the number of time steps such that the wake extends back by
                # some number of reference chord lengths.
                assert num_chords is not None
                wake_length = num_chords * max_c_ref
                distance_per_time_step = (
                    delta_time * operating_point_movement.base_operating_point.vCg__E
                )
                num_steps = math.ceil(wake_length / distance_per_time_step)
            else:
                # Set the number of time steps such that the simulation runs for
                # some number of cycles of all motions. Use the LCM of all periods
                # to ensure each motion completes an integer number of cycles.
                assert num_cycles is not None
                num_steps = math.ceil(num_cycles * _lcm_period / delta_time)

        # --- Resolve max_wake_rows ---
        # Validate max_wake_* parameters. At most one can be non None.
        _num_max_wake_set = sum(
            x is not None for x in (max_wake_rows, max_wake_chords, max_wake_cycles)
        )
        if _num_max_wake_set > 1:
            raise ValueError(
                "At most one of max_wake_rows, max_wake_chords, and max_wake_cycles "
                "can be non None."
            )

        if max_wake_chords is not None:
            if not _static:
                raise ValueError("max_wake_chords is only valid for static Movements.")
            max_wake_chords = _parameter_validation.int_in_range_return_int(
                max_wake_chords,
                "max_wake_chords",
                min_val=1,
                min_inclusive=True,
            )

        if max_wake_cycles is not None:
            if _static:
                raise ValueError(
                    "max_wake_cycles is only valid for non static Movements."
                )
            max_wake_cycles = _parameter_validation.int_in_range_return_int(
                max_wake_cycles,
                "max_wake_cycles",
                min_val=1,
                min_inclusive=True,
            )

        # Convert max_wake_chords to max_wake_rows using the same formula as
        # num_chords to num_steps.
        if max_wake_chords is not None:
            c_refs = []
            for airplane_movement in airplane_movements:
                c_ref = airplane_movement.base_airplane.c_ref
                assert c_ref is not None
                c_refs.append(c_ref)
            max_c_ref = max(c_refs)

            distance_per_time_step = (
                delta_time * operating_point_movement.base_operating_point.vCg__E
            )
            max_wake_rows = math.ceil(
                max_wake_chords * max_c_ref / distance_per_time_step
            )

        # Convert max_wake_cycles to max_wake_rows using the same formula as
        # num_cycles to num_steps.
        if max_wake_cycles is not None:
            max_wake_rows = math.ceil(max_wake_cycles * _lcm_period / delta_time)

        # --- Initialize CoreMovement ---
        super().__init__(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=delta_time,
            num_steps=num_steps,
            max_wake_rows=max_wake_rows,
        )

        # Pre-populate the lazy caches with values already computed above so
        # that accessing the inherited properties does not redundantly
        # recompute them. _max_period, _lcm_period, and _static are set here
        # because they are always needed during __init__ (via the static
        # check). _min_period remains lazy.
        self._max_period = _max_period
        self._lcm_period = _lcm_period
        self._static = _static

        # --- Store Movement only attributes ---
        self._num_cycles = num_cycles
        self._num_chords = num_chords
        self._max_wake_chords = max_wake_chords
        self._max_wake_cycles = max_wake_cycles

        # --- Batch generate ---
        # Generate a list of lists of Airplanes that are the steps through each
        # AirplaneMovement. The first index identifies the AirplaneMovement, and the
        # second index identifies the time step.
        airplanes_temp: list[list[geometry.airplane.Airplane]] = []
        for airplane_movement in self.airplane_movements:
            airplanes_temp.append(
                airplane_movement.generate_airplanes(
                    num_steps=self._num_steps, delta_time=self._delta_time
                )
            )

        # Validate that all Wings maintain their symmetry type across all time steps.
        for airplane_movement_id, airplane_list in enumerate(airplanes_temp):
            # Get the base Airplane (first time step).
            base_airplane = airplane_list[0]

            # Store the symmetry types of the base Wings.
            base_wing_symmetry_types = []
            for wing in base_airplane.wings:
                base_wing_symmetry_types.append(wing.symmetry_type)

            # Validate all subsequent time steps.
            for step_id, airplane in enumerate(airplane_list):
                # Check that Wings maintain their symmetry types.
                for wing_id, wing in enumerate(airplane.wings):
                    base_symmetry_type = base_wing_symmetry_types[wing_id]
                    if wing.symmetry_type != base_symmetry_type:
                        raise ValueError(
                            f"Wing {wing_id} in AirplaneMovement "
                            f"{airplane_movement_id} changed from type "
                            f"{base_symmetry_type} symmetry at time step 0 to type "
                            f"{wing.symmetry_type} symmetry at time step {step_id}. "
                            f"Wings cannot undergo motion that changes their symmetry "
                            f"type. This happens when a symmetric Wing moves such "
                            f"that its symmetry plane is no longer coincident with "
                            f"the wing axes' yz plane or vice versa."
                        )

        # Store as tuple of tuples to prevent external mutation.
        self._airplanes: tuple[tuple[geometry.airplane.Airplane, ...], ...] = tuple(
            tuple(airplane_list) for airplane_list in airplanes_temp
        )

        # Generate a list of OperatingPoints that are the steps through the
        # OperatingPointMovement.
        operating_points_temp = operating_point_movement.generate_operating_points(
            num_steps=self._num_steps, delta_time=self._delta_time
        )
        # Store as tuple to prevent external mutation.
        self._operating_points: tuple[operating_point_mod.OperatingPoint, ...] = tuple(
            operating_points_temp
        )

    # --- Immutable: read only properties ---
    @property
    def operating_point_movement(
        self,
    ) -> operating_point_movement_mod.OperatingPointMovement:
        assert isinstance(
            self._operating_point_movement,
            operating_point_movement_mod.OperatingPointMovement,
        )
        return self._operating_point_movement

    @property
    def airplane_movements(
        self,
    ) -> tuple[airplane_movement_mod.AirplaneMovement, ...]:
        return cast(
            tuple[airplane_movement_mod.AirplaneMovement, ...],
            self._airplane_movements,
        )

    @property
    def num_cycles(self) -> int | None:
        return self._num_cycles

    @property
    def num_chords(self) -> int | None:
        return self._num_chords

    @property
    def max_wake_chords(self) -> int | None:
        return self._max_wake_chords

    @property
    def max_wake_cycles(self) -> int | None:
        return self._max_wake_cycles

    @property
    def airplanes(self) -> tuple[tuple[geometry.airplane.Airplane, ...], ...]:
        return self._airplanes

    @property
    def operating_points(self) -> tuple[operating_point_mod.OperatingPoint, ...]:
        return self._operating_points


# Oversampling factor for the non static cached path. The high resolution
# Movement is built with _NON_STATIC_CACHE_OVERSAMPLE * max_num_steps
# intervals (and therefore one more snapshot) so the maximum and half maximum
# candidates have integer strides and other candidates stay within roughly
# half a high resolution step of the nominal time.
_NON_STATIC_CACHE_OVERSAMPLE: int = 2


def _compute_wake_area_mismatch(
    delta_time: float,
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    num_steps: int | None = None,
) -> float:
    """Computes the average area mismatch between wake and bound ring vortices.

    Builds a temporary Movement and UnsteadyProblem so the panel vertex points are
    populated in the first Airplane's geometry axes (GP1_CgP1), then computes bound
    trailing edge ring vortex areas and wake first row ring vortex areas analytically
    from the panel corners. Skips solver instantiation, the per step
    initialize_step_geometry walk, the bound corner stack population, and the wake
    corner stack population.

    The area mismatch metric measures how well the wake ring vortex sizing matches the
    bound ring vortex sizing. A lower value indicates better matching.

    The number of time steps checked is picked to capture the full range of differences
    in areas for the wake and bound ring vortex child parent pairs. For static cases,
    this is just a single time step. For non static cases, it is enough time steps to
    cover one full maximum length period of motion.

    Bound trailing edge ring vortex back vertices are computed exactly as
    UnsteadyRingVortexLatticeMethodSolver._initialize_panel_vortices_at would: at step
    0, panel.B[lr]pp + 0.25 * vInf * dt; at step k > 0, 0.75 * panel.B[lr]pp(k) + 0.25 *
    panel.B[lr]pp(k - 1) + 0.25 * vInf(k) * dt. For the first wake-row case
    corresponding to step == 1, the wake ring vortex area reduces algebraically to dt *
    |cross(Brrvp(k) - Blrvp(k), vInf(k - 1))|. For later steps (step >= 2), the wake
    first row ring vortex area depends on both current and previous bound trailing edge
    vertices and is computed using the general quadrilateral-diagonal formula. The
    advection velocity uses step k - 1 because the solver populates the wake using the
    previous step's freestream velocity.

    :param delta_time: The delta_time value to test. It must be a positive float. Its
        units are in seconds.
    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param num_steps: The number of time steps to simulate. If None (the default),
        derived from ceil(max_period / delta_time). Callers that already know the
        intended integer step count should pass it explicitly to match the semantics of
        the cached non-static evaluator and avoid floating point rounding in the
        inferred step count.
    :return: The average area mismatch. The absolute percent error between the area of
        shed wake ring vortices and the area of their parent bound ring vortices (at the
        time step where they were shed). Averaged across all time steps and all pairs of
        child and parent ring vortices. A lower value indicates better matching.
    """
    airplane_movements_copy = copy.deepcopy(airplane_movements)
    operating_point_movement_copy = copy.deepcopy(operating_point_movement)

    if num_steps is None:
        max_airplane_movement_period = 0.0
        for airplane_movement in airplane_movements_copy:
            max_airplane_movement_period = max(
                airplane_movement.max_period, max_airplane_movement_period
            )

        max_period = max(
            max_airplane_movement_period, operating_point_movement_copy.max_period
        )

        # Calculate the number of time steps to traverse the max period (or just a
        # single time step if there is no movement).
        num_steps = 1
        if max_period > 0.0:
            num_steps = math.ceil(max_period / delta_time)

    temp_movement = Movement(
        airplane_movements=airplane_movements_copy,
        operating_point_movement=operating_point_movement_copy,
        delta_time=delta_time,
        num_steps=num_steps,
    )

    # UnsteadyProblem populates panel.*_GP1_CgP1 attributes on every Airplane snapshot,
    # which is the frame the solver's bound ring vortex construction reads.
    temp_problem = problems.UnsteadyProblem(movement=temp_movement)

    if num_steps < 2:
        return 0.0

    # Cache vInf in GP1_CgP1 frame per step.
    v_inf_per_step = [
        steady_problem.operating_point.vInf_GP1__E
        for steady_problem in temp_problem.steady_problems
    ]

    total_mismatch = 0.0
    num_comparisons = 0

    for step in range(1, num_steps):
        current_problem = temp_problem.steady_problems[step]
        previous_problem = temp_problem.steady_problems[step - 1]
        v_inf_curr = v_inf_per_step[step]
        v_inf_prev = v_inf_per_step[step - 1]

        for airplane_id, current_airplane in enumerate(current_problem.airplanes):
            previous_airplane = previous_problem.airplanes[airplane_id]

            for wing_id, current_wing in enumerate(current_airplane.wings):
                previous_wing = previous_airplane.wings[wing_id]

                current_panels = current_wing.panels
                if current_panels is None:
                    continue
                previous_panels = previous_wing.panels
                if previous_panels is None:
                    continue

                num_chordwise_panels = previous_wing.num_chordwise_panels
                trailing_edge_chordwise_index = num_chordwise_panels - 1
                num_spanwise_panels = current_panels.shape[1]

                # Stack TE panel corners across spanwise positions for vectorized math.
                # Shape: (num_spanwise_panels, 3) per array.
                Flpp_curr = np.stack(
                    [
                        current_panels[trailing_edge_chordwise_index, s].Flpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Frpp_curr = np.stack(
                    [
                        current_panels[trailing_edge_chordwise_index, s].Frpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Blpp_curr = np.stack(
                    [
                        current_panels[trailing_edge_chordwise_index, s].Blpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Brpp_curr = np.stack(
                    [
                        current_panels[trailing_edge_chordwise_index, s].Brpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Flpp_prev = np.stack(
                    [
                        previous_panels[trailing_edge_chordwise_index, s].Flpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Frpp_prev = np.stack(
                    [
                        previous_panels[trailing_edge_chordwise_index, s].Frpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Blpp_prev = np.stack(
                    [
                        previous_panels[trailing_edge_chordwise_index, s].Blpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )
                Brpp_prev = np.stack(
                    [
                        previous_panels[trailing_edge_chordwise_index, s].Brpp_GP1_CgP1
                        for s in range(num_spanwise_panels)
                    ]
                )

                # Bound RingVortex front vertices at step - 1 (panel quarter chord).
                Flrvp_prev = 0.75 * Flpp_prev + 0.25 * Blpp_prev
                Frrvp_prev = 0.75 * Frpp_prev + 0.25 * Brpp_prev

                # Bound RingVortex back vertices at step - 1 (TE special case).
                if step - 1 == 0:
                    Blrvp_prev = Blpp_prev + 0.25 * v_inf_prev * delta_time
                    Brrvp_prev = Brpp_prev + 0.25 * v_inf_prev * delta_time
                else:
                    previous_previous_problem = temp_problem.steady_problems[step - 2]
                    previous_previous_panels = (
                        previous_previous_problem.airplanes[airplane_id]
                        .wings[wing_id]
                        .panels
                    )
                    assert previous_previous_panels is not None
                    Blpp_prev_prev = np.stack(
                        [
                            previous_previous_panels[
                                trailing_edge_chordwise_index, s
                            ].Blpp_GP1_CgP1
                            for s in range(num_spanwise_panels)
                        ]
                    )
                    Brpp_prev_prev = np.stack(
                        [
                            previous_previous_panels[
                                trailing_edge_chordwise_index, s
                            ].Brpp_GP1_CgP1
                            for s in range(num_spanwise_panels)
                        ]
                    )
                    Blrvp_prev = (
                        0.75 * Blpp_prev
                        + 0.25 * Blpp_prev_prev
                        + 0.25 * v_inf_prev * delta_time
                    )
                    Brrvp_prev = (
                        0.75 * Brpp_prev
                        + 0.25 * Brpp_prev_prev
                        + 0.25 * v_inf_prev * delta_time
                    )

                # Bound area = 0.5 * |cross(Frrvp - Blrvp, Flrvp - Brrvp)|.
                bound_diag1 = Frrvp_prev - Blrvp_prev
                bound_diag2 = Flrvp_prev - Brrvp_prev
                bound_areas = 0.5 * np.linalg.norm(
                    np.cross(bound_diag1, bound_diag2), axis=1
                )

                # Bound TE RingVortex back vertices at the current step. These define
                # the front of the wake first row RingVortex at the current step.
                Blrvp_curr = (
                    0.75 * Blpp_curr + 0.25 * Blpp_prev + 0.25 * v_inf_curr * delta_time
                )
                Brrvp_curr = (
                    0.75 * Brpp_curr + 0.25 * Brpp_prev + 0.25 * v_inf_curr * delta_time
                )

                # Wake first row area depends on how the solver's wake populator
                # constructs the back vertices. At step 1 (first wake population, run
                # at solver step 0), there is no preexisting wake grid: the populator
                # sets row 0 to the bound TE back vertices at step 1 and row 1 to
                # row 0 + vInf(0) * dt. The area collapses to
                # dt * |cross(Brrvp_curr - Blrvp_curr, vInf_prev)|. At step k >= 2,
                # the populator copies the previous step's wake grid forward, advects
                # it by vInf(k - 1) * dt, then prepends a new front row from bound TE
                # back vertices at step k. The first wake row's back vertices are
                # therefore Blrvp_prev + vInf_prev * dt and Brrvp_prev +
                # vInf_prev * dt, where Blrvp_prev and Brrvp_prev are the previous
                # step's bound TE back vertices already computed above.
                if step == 1:
                    span_vec = Brrvp_curr - Blrvp_curr
                    wake_areas = delta_time * np.linalg.norm(
                        np.cross(span_vec, v_inf_prev), axis=1
                    )
                else:
                    Flwrvp = Blrvp_curr
                    Frwrvp = Brrvp_curr
                    Blwrvp = Blrvp_prev + v_inf_prev * delta_time
                    Brwrvp = Brrvp_prev + v_inf_prev * delta_time
                    wake_diag1 = Frwrvp - Blwrvp
                    wake_diag2 = Flwrvp - Brwrvp
                    wake_areas = 0.5 * np.linalg.norm(
                        np.cross(wake_diag1, wake_diag2), axis=1
                    )

                epsilon = 1e-12
                valid = np.abs(bound_areas) > epsilon
                if np.any(valid):
                    total_mismatch += float(
                        np.sum(
                            np.abs(wake_areas[valid] - bound_areas[valid])
                            / bound_areas[valid]
                        )
                    )
                    num_comparisons += int(np.sum(valid))

    if num_comparisons == 0:
        return 0.0

    return total_mismatch / num_comparisons


def _compute_wake_area_mismatches_cached_non_static(
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    lcm_period: float,
    num_steps_candidates: list[int],
) -> dict[int, float]:
    """Computes wake area mismatch for many integer num_steps candidates from one shared
    high resolution Movement.

    Builds a single Movement with M = _NON_STATIC_CACHE_OVERSAMPLE *
    max(num_steps_candidates) snapshots covering one LCM period of motion. For each
    candidate N, the snapshots needed at delta_time = lcm_period / N are looked up via
    linear interpolation between adjacent high resolution samples, then the same closed
    form bound and wake area formulas as _compute_wake_area_mismatch are applied. The
    Movement and UnsteadyProblem cost is paid once per optimizer call rather than once
    per candidate, which collapses the dominant per evaluation cost into a single up
    front build.

    The linear interpolation introduces a per panel position error bounded by O(omega^2
    * dt^2 * amplitude) where dt is the high resolution sample spacing in time. For
    typical flapping kinematics this is small enough that the optimizer's selected
    num_steps matches the result that would be obtained from a fresh Movement and
    UnsteadyProblem build per candidate.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param lcm_period: The LCM of all motion periods. It must be a positive float. Its
        units are in seconds.
    :param num_steps_candidates: The integer num_steps values to score. Must be non
        empty and contain only positive ints.
    :return: A dict mapping each candidate num_steps value to its computed mismatch.
    """
    max_candidate = max(num_steps_candidates)
    high_res_num_intervals = _NON_STATIC_CACHE_OVERSAMPLE * max_candidate
    high_res_dt = lcm_period / high_res_num_intervals
    # The Movement covers high_res_num_intervals + 1 snapshots so that one full
    # LCM period worth of intervals can be sampled.
    high_res_num_steps = high_res_num_intervals + 1

    airplane_movements_copy = copy.deepcopy(airplane_movements)
    operating_point_movement_copy = copy.deepcopy(operating_point_movement)

    temp_movement = Movement(
        airplane_movements=airplane_movements_copy,
        operating_point_movement=operating_point_movement_copy,
        delta_time=high_res_dt,
        num_steps=high_res_num_steps,
    )
    temp_problem = problems.UnsteadyProblem(movement=temp_movement)

    first_problem = temp_problem.steady_problems[0]

    # Per (airplane, wing) cache of TE panel corners. Each stored array has
    # shape (high_res_num_steps, num_spanwise_panels, 3) so that strided
    # advanced indexing produces (N + 1, num_spanwise_panels, 3) views per
    # candidate.
    cache_per_wing: list[dict] = []
    for airplane_id, airplane in enumerate(first_problem.airplanes):
        for wing_id, wing in enumerate(airplane.wings):
            if wing.panels is None:
                continue
            num_chordwise = wing.num_chordwise_panels
            te_idx = num_chordwise - 1
            num_spanwise = wing.panels.shape[1]

            Flpp = np.empty((high_res_num_steps, num_spanwise, 3), dtype=float)
            Frpp = np.empty_like(Flpp)
            Blpp = np.empty_like(Flpp)
            Brpp = np.empty_like(Flpp)

            for step in range(high_res_num_steps):
                step_panels = (
                    temp_problem.steady_problems[step]
                    .airplanes[airplane_id]
                    .wings[wing_id]
                    .panels
                )
                assert step_panels is not None
                for s in range(num_spanwise):
                    panel = step_panels[te_idx, s]
                    Flpp[step, s] = panel.Flpp_GP1_CgP1
                    Frpp[step, s] = panel.Frpp_GP1_CgP1
                    Blpp[step, s] = panel.Blpp_GP1_CgP1
                    Brpp[step, s] = panel.Brpp_GP1_CgP1

            cache_per_wing.append(
                {
                    "Flpp": Flpp,
                    "Frpp": Frpp,
                    "Blpp": Blpp,
                    "Brpp": Brpp,
                }
            )

    v_inf_high_res = np.array(
        [
            steady_problem.operating_point.vInf_GP1__E
            for steady_problem in temp_problem.steady_problems
        ]
    )

    results: dict[int, float] = {}
    for num_steps in num_steps_candidates:
        results[num_steps] = _evaluate_cached_wake_area_mismatch(
            cache_per_wing=cache_per_wing,
            v_inf_high_res=v_inf_high_res,
            lcm_period=lcm_period,
            high_res_num_intervals=high_res_num_intervals,
            num_steps=num_steps,
        )
    return results


def _evaluate_cached_wake_area_mismatch(
    cache_per_wing: list[dict],
    v_inf_high_res: np.ndarray,
    lcm_period: float,
    high_res_num_intervals: int,
    num_steps: int,
) -> float:
    """Evaluates wake area mismatch at one candidate num_steps using a shared cache.

    Maps each candidate step k in [0, num_steps - 1] to a fractional high resolution
    index k * high_res_num_intervals / num_steps and linearly interpolates the cached
    panel corners and freestream between adjacent high resolution samples. Linear
    interpolation has O(dt^2) error in panel position for smooth motion, which is small
    enough at the oversampling factors used here that the integer num_steps picked by
    the optimizer matches the result obtained from a fresh per candidate Movement and
    UnsteadyProblem build. When the fractional index is an integer (e.g. for candidates
    that divide high_res_num_intervals exactly) the interpolation reduces to a direct
    lookup. The resulting per panel positions feed the same closed form bound and wake
    area formulas as _compute_wake_area_mismatch.

    :param cache_per_wing: A list of per (airplane, wing) caches produced by
        _compute_wake_area_mismatches_cached_non_static. Each entry is a dict with keys
        "Flpp", "Frpp", "Blpp", "Brpp" mapping to (high_res_num_intervals + 1,
        num_spanwise_panels, 3) ndarrays of trailing edge panel corner positions in the
        first Airplane's geometry axes (GP1_CgP1).
    :param v_inf_high_res: A (high_res_num_intervals + 1, 3) ndarray of freestream
        velocity samples in the first Airplane's geometry axes (GP1_CgP1), one per high
        resolution snapshot.
    :param lcm_period: The LCM of all motion periods. It must be a positive float. Its
        units are in seconds.
    :param high_res_num_intervals: The number of high resolution intervals covering one
        LCM period. The cached arrays have high_res_num_intervals + 1 samples along
        their leading axis.
    :param num_steps: The candidate integer number of steps per LCM period to score.
        Must be a positive int. Returns 0.0 immediately for num_steps < 2 since at least
        one step pair is needed for a comparison.
    :return: The average area mismatch at this candidate. The absolute percent error
        between the area of shed wake RingVortices and the area of their parent bound
        RingVortices (at the time step where they were shed). Averaged across all time
        steps and all pairs of child and parent RingVortices. A lower value indicates
        better matching.
    """
    if num_steps < 2:
        return 0.0

    delta_time = lcm_period / num_steps

    # Linear interpolation mapping from candidate steps to high resolution
    # samples. The floor index gives the lower bracketing sample; the fractional
    # part is the interpolation weight on the next sample. Only num_steps
    # candidate samples are generated since the comparison loop reads indices
    # [0, num_steps - 1]; the t = lcm_period sample would be dead.
    fractional_indices = np.arange(num_steps) * high_res_num_intervals / num_steps
    floor_indices = np.floor(fractional_indices).astype(int)
    weights_next = fractional_indices - floor_indices
    next_indices = floor_indices + 1
    weights_floor = 1.0 - weights_next
    # Reshape weights for broadcasting against (num_samples, num_spanwise, 3)
    # panel arrays.
    weights_floor_b = weights_floor[:, None, None]
    weights_next_b = weights_next[:, None, None]
    # Reshape weights for broadcasting against (num_samples, 3) vInf arrays.
    weights_floor_v = weights_floor[:, None]
    weights_next_v = weights_next[:, None]

    total_mismatch = 0.0
    num_comparisons = 0

    for cache in cache_per_wing:
        # Sample and linearly interpolate the cached panel attribute arrays at
        # the candidate step times. Each result has shape (num_steps,
        # num_spanwise_panels, 3); v_inf has shape (num_steps, 3).
        Flpp = (
            weights_floor_b * cache["Flpp"][floor_indices]
            + weights_next_b * cache["Flpp"][next_indices]
        )
        Frpp = (
            weights_floor_b * cache["Frpp"][floor_indices]
            + weights_next_b * cache["Frpp"][next_indices]
        )
        Blpp = (
            weights_floor_b * cache["Blpp"][floor_indices]
            + weights_next_b * cache["Blpp"][next_indices]
        )
        Brpp = (
            weights_floor_b * cache["Brpp"][floor_indices]
            + weights_next_b * cache["Brpp"][next_indices]
        )
        v_inf = (
            weights_floor_v * v_inf_high_res[floor_indices]
            + weights_next_v * v_inf_high_res[next_indices]
        )

        # Comparison axis: index i in [0, num_steps - 2] maps to
        # step = i + 1, prev_step = i. Slice the panel arrays accordingly so
        # bound and wake area calculations can run as one vectorized batch
        # rather than a Python step loop.
        Flpp_prev_arr = Flpp[: num_steps - 1]
        Frpp_prev_arr = Frpp[: num_steps - 1]
        Blpp_prev_arr = Blpp[: num_steps - 1]
        Brpp_prev_arr = Brpp[: num_steps - 1]
        Blpp_curr_arr = Blpp[1:num_steps]
        Brpp_curr_arr = Brpp[1:num_steps]
        v_inf_prev_arr = v_inf[: num_steps - 1, None, :]
        v_inf_curr_arr = v_inf[1:num_steps, None, :]

        # Bound RingVortex front vertices at the previous step (panel quarter
        # chord, no derivative term).
        Flrvp_prev = 0.75 * Flpp_prev_arr + 0.25 * Blpp_prev_arr
        Frrvp_prev = 0.75 * Frpp_prev_arr + 0.25 * Brpp_prev_arr

        # Bound RingVortex back vertices at the previous step. For comparison
        # index 0 (prev_step == 0) there is no panel at step - 2; substituting
        # the step - 1 panel makes the derivative formula collapse to the
        # solver's no-derivative form 0.75*Blpp + 0.25*Blpp + 0.25*vInf*dt =
        # Blpp + 0.25*vInf*dt, so we get the correct value for both i == 0 and
        # i >= 1 from a single vectorized expression.
        Blpp_prev_prev_arr = np.concatenate(
            [Blpp_prev_arr[:1], Blpp_prev_arr[:-1]], axis=0
        )
        Brpp_prev_prev_arr = np.concatenate(
            [Brpp_prev_arr[:1], Brpp_prev_arr[:-1]], axis=0
        )
        Blrvp_prev = (
            0.75 * Blpp_prev_arr
            + 0.25 * Blpp_prev_prev_arr
            + 0.25 * v_inf_prev_arr * delta_time
        )
        Brrvp_prev = (
            0.75 * Brpp_prev_arr
            + 0.25 * Brpp_prev_prev_arr
            + 0.25 * v_inf_prev_arr * delta_time
        )

        # Bound areas across all comparison steps. Shape: (num_steps - 1,
        # num_spanwise_panels).
        bound_diag1 = Frrvp_prev - Blrvp_prev
        bound_diag2 = Flrvp_prev - Brrvp_prev
        bound_areas = 0.5 * np.linalg.norm(np.cross(bound_diag1, bound_diag2), axis=-1)

        # Bound RingVortex back vertices at the current step. These are the
        # front vertices of the wake first row at the current step.
        Blrvp_curr = (
            0.75 * Blpp_curr_arr
            + 0.25 * Blpp_prev_arr
            + 0.25 * v_inf_curr_arr * delta_time
        )
        Brrvp_curr = (
            0.75 * Brpp_curr_arr
            + 0.25 * Brpp_prev_arr
            + 0.25 * v_inf_curr_arr * delta_time
        )

        # Wake first row back vertices. At comparison index 0 (step == 1, first
        # wake population) the solver sets row 1 = row 0 + vInf * dt, so the
        # wake's back row is the bound TE back at the *current* step advected.
        # At i >= 1 the wake's back row is the bound TE back at the *previous*
        # step advected. Substituting Blrvp_curr[0] for index 0 of the previous
        # step array makes a single vectorized expression cover both cases.
        Blrvp_prev_for_wake = np.concatenate([Blrvp_curr[:1], Blrvp_prev[1:]], axis=0)
        Brrvp_prev_for_wake = np.concatenate([Brrvp_curr[:1], Brrvp_prev[1:]], axis=0)
        Flwrvp = Blrvp_curr
        Frwrvp = Brrvp_curr
        Blwrvp = Blrvp_prev_for_wake + v_inf_prev_arr * delta_time
        Brwrvp = Brrvp_prev_for_wake + v_inf_prev_arr * delta_time
        wake_diag1 = Frwrvp - Blwrvp
        wake_diag2 = Flwrvp - Brwrvp
        wake_areas = 0.5 * np.linalg.norm(np.cross(wake_diag1, wake_diag2), axis=-1)

        epsilon = 1e-12
        valid = np.abs(bound_areas) > epsilon
        if np.any(valid):
            total_mismatch += float(
                np.sum(
                    np.abs(wake_areas[valid] - bound_areas[valid]) / bound_areas[valid]
                )
            )
            num_comparisons += int(np.sum(valid))

    if num_comparisons == 0:
        return 0.0
    return total_mismatch / num_comparisons


def _optimize_delta_time(
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    initial_delta_time: float,
    mismatch_cutoff: float = 0.01,
) -> float:
    """Finds an optimal delta_time by minimizing wake area mismatch.

    Optimizes delta_time to minimize the area mismatch between wake ring vortices and
    their parent bound trailing edge ring vortices. This produces better results at high
    Strouhal numbers where motion induced velocity is significant.

    For non static Movements, the optimization uses a brute force search over integer
    num_steps_per_lcm_cycle values (from 0.5x to 2x the initial estimate). This ensures
    the resulting delta_time divides the LCM period evenly.

    For static Movements, the optimization uses scipy.optimize.minimize_scalar with
    early termination if the mismatch falls below the specified cutoff value.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate, typically the result of
        _analytically_optimize_delta_time. It must be a positive float. Its units are in
        seconds.
    :param mismatch_cutoff: A positive float for the optimization's convergence
        threshold. Only used for static Movements. When the average area mismatch (which
        is an absolute percent error) falls below this value, the search terminates
        early. The default is 0.01.
    :return: The optimized delta_time value. Its units are in seconds.
    """
    movement_logger.info("Starting delta_time optimization.")

    # Collect all non zero periods to determine if static.
    all_periods: list[float] = []
    for airplane_movement in airplane_movements:
        all_periods.extend(airplane_movement.all_periods)
    op_period = operating_point_movement.max_period
    if op_period != 0.0:
        all_periods.append(op_period)
    non_zero_periods = [p for p in all_periods if p != 0.0]

    if not non_zero_periods:
        # Static case: use scipy continuous optimization.
        return _optimize_delta_time_static(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
            mismatch_cutoff=mismatch_cutoff,
        )
    else:
        # Non static case: brute force search over integer num_steps_per_lcm_cycle.
        lcm_period = _core.lcm_multiple(non_zero_periods)
        return _optimize_delta_time_non_static(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
            lcm_period=lcm_period,
        )


def _optimize_delta_time_static(
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    initial_delta_time: float,
    mismatch_cutoff: float,
) -> float:
    """Optimizes delta_time for static Movements using scipy.

    Uses scipy.optimize.minimize_scalar to find the delta_time that minimizes wake area
    mismatch. The search terminates early if the mismatch falls below the specified
    cutoff value.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate. It must be a positive float. Its
        units are in seconds.
    :param mismatch_cutoff: A positive float for early termination. When the average
        area mismatch falls below this value, the search terminates early.
    :return: The optimized delta_time value. Its units are in seconds.
    """
    lower_bound = initial_delta_time / 2.0
    upper_bound = initial_delta_time * 2.0

    # Check initial estimate first before running optimizer.
    initial_mismatch = _compute_wake_area_mismatch(
        initial_delta_time, airplane_movements, operating_point_movement
    )

    dt_str = str(round(initial_delta_time, 6))
    mismatch_str = str(round(initial_mismatch, 6))

    state_msg = "\tState: delta_time=" + dt_str
    obj_msg = "\t\tMismatch: " + mismatch_str

    movement_logger.info(state_msg)
    movement_logger.info(obj_msg)

    if initial_mismatch < mismatch_cutoff:
        movement_logger.info("Acceptable value reached.")
        movement_logger.info("Optimization complete.")
        return initial_delta_time

    best_delta_time = initial_delta_time
    best_mismatch = initial_mismatch

    def objective(dt: float) -> float:
        nonlocal best_delta_time, best_mismatch
        mismatch = _compute_wake_area_mismatch(
            dt, airplane_movements, operating_point_movement
        )

        this_dt_str = str(round(dt, 6))
        this_mismatch_str = str(round(mismatch, 6))

        this_state_msg = "\tState: delta_time=" + this_dt_str
        this_obj_msg = "\t\tMismatch: " + this_mismatch_str

        movement_logger.info(this_state_msg)
        movement_logger.info(this_obj_msg)

        if mismatch < best_mismatch:
            best_mismatch = mismatch
            best_delta_time = dt

        if mismatch < mismatch_cutoff:
            raise StopIteration

        return mismatch

    try:
        result = sp_opt.minimize_scalar(
            objective,
            bounds=(lower_bound, upper_bound),
            method="bounded",
            options={"xatol": 0.001},
        )

        if not result.success:
            raise RuntimeError("delta_time optimization failed to converge.")

        optimized_delta_time = float(result.x)
    except StopIteration:
        optimized_delta_time = best_delta_time
        movement_logger.info("Acceptable value reached.")

    # Warn if the optimized value is at one of the bounds.
    bound_tolerance = 1e-6
    if abs(optimized_delta_time - lower_bound) < bound_tolerance:
        movement_logger.warning(
            "Optimized delta_time is at the lower bound. A better value may exist "
            "below the search range."
        )
    elif abs(optimized_delta_time - upper_bound) < bound_tolerance:
        movement_logger.warning(
            "Optimized delta_time is at the upper bound. A better value may exist "
            "above the search range."
        )

    movement_logger.info("Optimization complete.")

    return optimized_delta_time


def _optimize_delta_time_non_static(
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    initial_delta_time: float,
    lcm_period: float,
) -> float:
    """Optimizes delta_time for non static Movements using brute force search.

    Searches over integer num_steps_per_lcm_cycle values from 0.5x to 2x the initial
    estimate. This ensures the resulting delta_time divides the LCM period evenly, which
    is important for clean cycle boundaries when averaging.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate. It must be a positive float. Its
        units are in seconds.
    :param lcm_period: The LCM of all motion periods. It must be a positive float. Its
        units are in seconds.
    :return: The optimized delta_time value. Its units are in seconds.
    """
    # Convert initial_delta_time to num_steps_per_lcm_cycle.
    initial_num_steps = lcm_period / initial_delta_time

    # Search from 0.5x to 2x the initial estimate.
    min_num_steps = max(1, int(initial_num_steps / 2))
    max_num_steps = int(initial_num_steps * 2) + 1

    movement_logger.info(
        "\tSearching num_steps_per_lcm_cycle from "
        + str(min_num_steps)
        + " to "
        + str(max_num_steps)
    )

    candidates = list(range(min_num_steps, max_num_steps + 1))

    # Build one high resolution Movement and score every integer candidate against it
    # by linearly interpolating the cached panel corners. This collapses the per
    # candidate Movement and UnsteadyProblem cost into a single up front build for
    # the entire bracket.
    cached_mismatches = _compute_wake_area_mismatches_cached_non_static(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        lcm_period=lcm_period,
        num_steps_candidates=candidates,
    )
    best_num_steps = min_num_steps
    best_mismatch = float("inf")
    for num_steps in candidates:
        mismatch = cached_mismatches[num_steps]
        delta_time = lcm_period / num_steps
        movement_logger.info(
            "\tnum_steps="
            + str(num_steps)
            + ", delta_time="
            + str(round(delta_time, 6))
            + ", mismatch="
            + str(round(mismatch, 6))
        )
        if mismatch < best_mismatch:
            best_mismatch = mismatch
            best_num_steps = num_steps

    optimized_delta_time = lcm_period / best_num_steps

    movement_logger.info(
        "Optimization complete. Best: num_steps_per_lcm_cycle="
        + str(best_num_steps)
        + ", delta_time="
        + str(round(optimized_delta_time, 6))
        + ", mismatch="
        + str(round(best_mismatch, 6))
    )

    # Warn if the optimized value is at one of the bounds.
    if best_num_steps == min_num_steps:
        movement_logger.warning(
            "Optimized num_steps_per_lcm_cycle is at the lower bound ("
            + str(min_num_steps)
            + "). A better value may exist below the search range."
        )
    elif best_num_steps == max_num_steps:
        movement_logger.warning(
            "Optimized num_steps_per_lcm_cycle is at the upper bound ("
            + str(max_num_steps)
            + "). A better value may exist above the search range."
        )

    return optimized_delta_time


def _analytically_optimize_delta_time(
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    initial_delta_time: float,
) -> float:
    """Analytically estimates the optimal delta_time from wake displacement.

    Estimates the delta_time that produces wake ring vortices with roughly the same
    chord length as the bound trailing edge ring vortices, accounting for both
    freestream and geometry motion velocities. This is faster than _optimize_delta_time
    but may be slightly less accurate.

    The algorithm works by: (1) computing a very small preliminary delta_time as the
    minimum motion period divided by 100 (capped by a maximum of 1000 total time steps
    for one cycle of the motion's LCM period), (2) creating a temporary Movement and
    UnsteadyProblem to generate Airplane geometry at each time step (with all Panel
    coordinates transformed into the first Airplane's geometry axes), (3) measuring the
    average wake displacement per time step for each Wing's trailing edge Panels
    (combining freestream velocity and geometry motion, both in the first Airplane's
    geometry axes), and (4) choosing a delta_time such that each wake ring vortex has
    approximately the same chord as its parent bound ring vortex (averaged across the
    one LCM period).

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate from the chord-based seed. It must
        be a positive float. Its units are in seconds. Used as a fallback for static
        Movements or degenerate cases.
    :return: The analytically optimized delta_time value. Its units are in seconds.
    """
    movement_logger.info("Starting analytical delta_time optimization.")

    # Collect all non zero periods.
    all_periods: list[float] = []
    for airplane_movement in airplane_movements:
        all_periods.extend(airplane_movement.all_periods)
    op_period = operating_point_movement.max_period
    if op_period != 0.0:
        all_periods.append(op_period)
    non_zero_periods = [p for p in all_periods if p != 0.0]

    # If there is no motion, fall back to the initial estimate.
    if not non_zero_periods:
        movement_logger.info(
            "All motion is static. Returning initial delta_time estimate."
        )
        return initial_delta_time

    min_period = min(non_zero_periods)

    lcm_period = _core.lcm_multiple(non_zero_periods)

    # Step 1: Compute a preliminary delta_time that divides the LCM period into roughly
    # min_period / 100 sized steps. Cap at 1000 steps to prevent excessive computation
    # for cases with very large LCM to min period ratios.
    preliminary_delta_time = min_period / 100.0
    preliminary_num_steps = round(lcm_period / preliminary_delta_time)
    max_preliminary_steps = 1000
    if preliminary_num_steps > max_preliminary_steps:
        movement_logger.warning(
            "Preliminary num_steps ("
            + str(preliminary_num_steps)
            + ") exceeds cap of "
            + str(max_preliminary_steps)
            + ". Capping to prevent excessive computation."
        )
        preliminary_num_steps = max_preliminary_steps
    if preliminary_num_steps < 1:
        preliminary_num_steps = 1
    # Adjust to ensure an integer number of steps fits the LCM period.
    preliminary_delta_time = lcm_period / preliminary_num_steps

    # Step 2: Deep copy the movement objects and create a temporary Movement and
    # UnsteadyProblem. The UnsteadyProblem creates SteadyProblems that transform all
    # Panel coordinates into the first Airplane's geometry axes (GP1_CgP1), ensuring
    # consistent coordinate frames for multi airplane formation flight cases.
    airplane_movements_copy = copy.deepcopy(airplane_movements)
    operating_point_movement_copy = copy.deepcopy(operating_point_movement)

    temp_movement = Movement(
        airplane_movements=airplane_movements_copy,
        operating_point_movement=operating_point_movement_copy,
        delta_time=preliminary_delta_time,
        num_steps=preliminary_num_steps + 1,
    )

    temp_problem = problems.UnsteadyProblem(movement=temp_movement)

    # Step 3: For each Airplane and Wing, measure the average wake displacement of
    # trailing edge Panels across all time steps. All coordinates are in the first
    # Airplane's geometry axes (GP1), relative to the first Airplane's CG (CgP1).
    wing_num_steps_values: list[float] = []
    wing_num_spanwise_panels_values: list[int] = []

    num_airplanes = len(temp_problem.steady_problems[0].airplanes)
    for airplane_id in range(num_airplanes):
        # Use the first generated Airplane to get Wing properties, since the base
        # Airplane may have different panel counts (e.g., due to symmetry processing
        # that occurs during generate_airplanes).
        first_airplane = temp_problem.steady_problems[0].airplanes[airplane_id]
        for wing_id, wing in enumerate(first_airplane.wings):
            num_chordwise = wing.num_chordwise_panels
            _panels = wing.panels
            assert _panels is not None
            num_spanwise = _panels.shape[1]

            # Compute the mean chordwise width of trailing edge Panels. This is
            # the target chord length for wake ring vortices. We use the trailing
            # edge Panel chord directly (rather than standard_mean_chord /
            # num_chordwise) because non uniform chordwise spacing can cause
            # trailing edge Panels to have different chords than average.
            total_te_panel_chord = 0.0
            for spanwise_id in range(num_spanwise):
                te_panel = _panels[num_chordwise - 1, spanwise_id]
                _leftLeg_GP1 = te_panel.leftLeg_GP1
                _rightLeg_GP1 = te_panel.rightLeg_GP1
                assert _leftLeg_GP1 is not None
                assert _rightLeg_GP1 is not None
                # Average the left and right leg magnitudes for this Panel.
                panel_chord = (
                    float(np.linalg.norm(_leftLeg_GP1))
                    + float(np.linalg.norm(_rightLeg_GP1))
                ) / 2.0
                total_te_panel_chord += panel_chord
            mean_te_panel_chord = total_te_panel_chord / num_spanwise

            # Accumulate displacement distance across all trailing edge Panels
            # across all consecutive time step pairs.
            total_distance = 0.0
            num_measurements = 0

            for spanwise_id in range(num_spanwise):
                spanwise_distance = 0.0

                for step in range(preliminary_num_steps):
                    # Get this and next time step's trailing edge Panel.
                    this_airplane = temp_problem.steady_problems[step].airplanes[
                        airplane_id
                    ]
                    next_airplane = temp_problem.steady_problems[step + 1].airplanes[
                        airplane_id
                    ]

                    _this_panels = this_airplane.wings[wing_id].panels
                    _next_panels = next_airplane.wings[wing_id].panels
                    assert _this_panels is not None
                    assert _next_panels is not None

                    this_panel = _this_panels[num_chordwise - 1, spanwise_id]
                    next_panel = _next_panels[num_chordwise - 1, spanwise_id]

                    # Compute the trailing edge center at each time step (in the first
                    # Airplane's geometry axes, relative to the first Airplane's CG).
                    _thisBlpp_GP1_CgP1 = this_panel.Blpp_GP1_CgP1
                    _thisBrpp_GP1_CgP1 = this_panel.Brpp_GP1_CgP1
                    _nextBlpp_GP1_CgP1 = next_panel.Blpp_GP1_CgP1
                    _nextBrpp_GP1_CgP1 = next_panel.Brpp_GP1_CgP1
                    assert _thisBlpp_GP1_CgP1 is not None
                    assert _thisBrpp_GP1_CgP1 is not None
                    assert _nextBlpp_GP1_CgP1 is not None
                    assert _nextBrpp_GP1_CgP1 is not None

                    thisCenter_GP1_CgP1 = (
                        _thisBlpp_GP1_CgP1 + _thisBrpp_GP1_CgP1
                    ) / 2.0
                    nextCenter_GP1_CgP1 = (
                        _nextBlpp_GP1_CgP1 + _nextBrpp_GP1_CgP1
                    ) / 2.0

                    # Geometry displacement of the trailing edge (in the first
                    # Airplane's geometry axes).
                    geometryDisplacement_GP1 = nextCenter_GP1_CgP1 - thisCenter_GP1_CgP1

                    # Flow displacement during this time step (in the first Airplane's
                    # geometry axes).
                    vInf_GP1__E = temp_movement.operating_points[step].vInf_GP1__E
                    flowDisplacement_GP1 = vInf_GP1__E * preliminary_delta_time

                    # Wake displacement is the combination of freestream flow and the
                    # opposite of geometry motion (the wake stays where the trailing
                    # edge was, then convects with the flow). Both vectors are in the
                    # first Airplane's geometry axes.
                    wakeDisplacement_GP1 = (
                        flowDisplacement_GP1 - geometryDisplacement_GP1
                    )

                    spanwise_distance += float(np.linalg.norm(wakeDisplacement_GP1))

                total_distance += spanwise_distance
                num_measurements += 1

            # Average distance over one LCM period per spanwise Panel.
            average_distance = total_distance / num_measurements
            # The desired number of steps per LCM period is such that each wake
            # ring vortex chord equals the trailing edge bound Panel chord.
            wing_num_steps_per_lcm = average_distance / mean_te_panel_chord
            wing_num_steps_values.append(wing_num_steps_per_lcm)
            wing_num_spanwise_panels_values.append(num_spanwise)

    # Step 4: Compute the weighted average of num_steps across all Wings.
    if not wing_num_steps_values:
        movement_logger.info(
            "No valid wake displacement data. Returning initial delta_time estimate."
        )
        return initial_delta_time

    total_weight = sum(wing_num_spanwise_panels_values)
    weighted_num_steps = (
        sum(
            n * w
            for n, w in zip(wing_num_steps_values, wing_num_spanwise_panels_values)
        )
        / total_weight
    )

    if weighted_num_steps <= 0.0:
        movement_logger.info(
            "Computed num_steps is non positive. Returning initial delta_time estimate."
        )
        return initial_delta_time

    # Round to an integer number of steps that fits the LCM period.
    final_num_steps = round(weighted_num_steps)
    if final_num_steps < 1:
        final_num_steps = 1
    optimized_delta_time = lcm_period / final_num_steps

    dt_str = str(round(optimized_delta_time, 6))
    movement_logger.info(
        "\tResult: delta_time="
        + dt_str
        + " ("
        + str(final_num_steps)
        + " steps per LCM period)"
    )

    # Warn if the result implies fewer than 20 time steps per minimum period of motion.
    # This indicates the trailing edge Panels are large relative to the motion, so
    # matching wake ring vortex area to bound ring vortex area requires a time step that
    # is large compared to the minimum period. This results in poor temporal resolution
    # of the fastest motion. Increasing the number of chordwise Panels or switching to
    # cosine chordwise spacing (which concentrates Panels at the leading and trailing
    # edges) will reduce the trailing edge Panels' chord lengths and allow finer time
    # stepping.
    steps_per_min_period = min_period / optimized_delta_time
    if steps_per_min_period < 20:
        movement_logger.warning(
            "Analytical optimization result implies only "
            + str(round(steps_per_min_period, 1))
            + " time steps per minimum period of motion. This may cause poor temporal "
            + "resolution. Consider increasing the number of chordwise Panels or "
            + "switching to cosine chordwise spacing to reduce the trailing edge "
            + "Panels' chord lengths."
        )

    movement_logger.info("Analytical optimization complete.")

    return optimized_delta_time
