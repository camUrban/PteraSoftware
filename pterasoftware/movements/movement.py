"""Contains the Movement class.

**Contains the following classes:**

Movement: A class used to contain an UnsteadyProblem's movement.

**Contains the following functions:**

None
"""

from __future__ import annotations

import copy
import math

import numpy as np
import scipy.optimize as sp_opt

from .. import _logging, _parameter_validation, _vortices, geometry
from .. import operating_point as operating_point_mod
from .. import problems
from . import airplane_movement as airplane_movement_mod
from . import operating_point_movement as operating_point_movement_mod

movement_logger = _logging.get_logger("movements.movement")


def _lcm(a: float, b: float) -> float:
    """Calculates the least common multiple of two numbers.

    :param a: First number (period in seconds)
    :param b: Second number (period in seconds)
    :return: LCM of a and b. Returns 0.0 if either input is 0.0.
    """
    if a == 0.0 or b == 0.0:
        return 0.0
    # Convert to integers (periods are typically whole multiples of delta_time)
    # Use sufficiently large multiplier to preserve precision
    multiplier = 1000000
    a_int = int(round(a * multiplier))
    b_int = int(round(b * multiplier))
    lcm_int = abs(a_int * b_int) // math.gcd(a_int, b_int)
    return lcm_int / multiplier


def _lcm_multiple(periods: list[float]) -> float:
    """Calculates the least common multiple of multiple periods.

    :param periods: List of periods in seconds
    :return: LCM of all periods. Returns 0.0 if all periods are 0.0.
    """
    if not periods or all(p == 0.0 for p in periods):
        return 0.0
    # Filter out zero periods and calculate LCM
    non_zero_periods = [p for p in periods if p != 0.0]
    if not non_zero_periods:
        return 0.0
    result = non_zero_periods[0]
    for period in non_zero_periods[1:]:
        result = _lcm(result, period)
    return result


class Movement:
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

    def __init__(
        self,
        airplane_movements: list[airplane_movement_mod.AirplaneMovement],
        operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
        delta_time: float | int | str | None = None,
        num_cycles: int | None = None,
        num_chords: int | None = None,
        num_steps: int | None = None,
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
            RingVortices with roughly the same chord length as the bound trailing edge
            RingVortices, accounting for both freestream and geometry motion velocities.
            This provides good results across all Strouhal numbers. "optimize": Movement
            first runs the analytical estimation, then uses that result as an initial
            guess for an iterative optimization that minimizes the area mismatch between
            wake RingVortices and their parent bound trailing edge RingVortices. This is
            slower but may produce slightly more accurate results. Positive number (int
            or float): Use the specified value directly. All values are converted
            internally to floats. The units are in seconds.
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
        :return: None
        """
        # Validate and store immutable attributes.
        if not isinstance(airplane_movements, list):
            raise TypeError("airplane_movements must be a list.")
        if len(airplane_movements) < 1:
            raise ValueError("airplane_movements must have at least one element.")
        for airplane_movement in airplane_movements:
            if not isinstance(
                airplane_movement, airplane_movement_mod.AirplaneMovement
            ):
                raise TypeError(
                    "Every element in airplane_movements must be an AirplaneMovement."
                )
        # Store as tuple to prevent external mutation.
        self._airplane_movements: tuple[airplane_movement_mod.AirplaneMovement, ...] = (
            tuple(airplane_movements)
        )

        if not isinstance(
            operating_point_movement,
            operating_point_movement_mod.OperatingPointMovement,
        ):
            raise TypeError(
                "operating_point_movement must be an OperatingPointMovement."
            )
        self._operating_point_movement = operating_point_movement

        # Initialize the caches for the properties derived from the immutable
        # attributes. These are initialized early because static is accessed below
        # during __init__ to determine num_steps calculation.
        self._lcm_period: float | None = None
        self._max_period: float | None = None
        self._min_period: float | None = None
        self._static: bool | None = None

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
            for airplane_movement in self._airplane_movements:
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
                airplane_movements=list(self._airplane_movements),
                operating_point_movement=self._operating_point_movement,
                initial_delta_time=fast_estimate,
            )

        # Run iterative optimization if requested, using the analytical result
        # as the initial guess.
        if _should_iteratively_optimize_delta_time:
            delta_time = _optimize_delta_time(
                airplane_movements=list(self._airplane_movements),
                operating_point_movement=self._operating_point_movement,
                initial_delta_time=delta_time,
            )
        self._delta_time: float = delta_time

        _static = self.static

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
        self._num_cycles = num_cycles

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
        self._num_chords = num_chords

        if self._num_cycles is not None or self._num_chords is not None:
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
                for airplane_movement in self._airplane_movements:
                    c_ref = airplane_movement.base_airplane.c_ref
                    assert c_ref is not None
                    c_refs.append(c_ref)
                max_c_ref = max(c_refs)

                # Set the number of time steps such that the wake extends back by
                # some number of reference chord lengths.
                assert self._num_chords is not None
                wake_length = self._num_chords * max_c_ref
                distance_per_time_step = (
                    delta_time
                    * self._operating_point_movement.base_operating_point.vCg__E
                )
                num_steps = math.ceil(wake_length / distance_per_time_step)
            else:
                # Set the number of time steps such that the simulation runs for some
                # number of cycles of all motions. Use the LCM of all periods to ensure
                # each motion completes an integer number of cycles.
                assert self._num_cycles is not None
                num_steps = math.ceil(
                    self._num_cycles * self.lcm_period / self._delta_time
                )
        self._num_steps: int = num_steps

        # Generate a list of lists of Airplanes that are the steps through each
        # AirplaneMovement. The first index identifies the AirplaneMovement, and the
        # second index identifies the time step.
        airplanes_temp: list[list[geometry.airplane.Airplane]] = []
        for airplane_movement in self._airplane_movements:
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

        # Generate a lists of OperatingPoints that are the steps through the
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
    def airplane_movements(self) -> tuple[airplane_movement_mod.AirplaneMovement, ...]:
        return self._airplane_movements

    @property
    def operating_point_movement(
        self,
    ) -> operating_point_movement_mod.OperatingPointMovement:
        return self._operating_point_movement

    @property
    def delta_time(self) -> float:
        return self._delta_time

    @property
    def num_cycles(self) -> int | None:
        return self._num_cycles

    @property
    def num_chords(self) -> int | None:
        return self._num_chords

    @property
    def num_steps(self) -> int:
        return self._num_steps

    @property
    def airplanes(self) -> tuple[tuple[geometry.airplane.Airplane, ...], ...]:
        return self._airplanes

    @property
    def operating_points(self) -> tuple[operating_point_mod.OperatingPoint, ...]:
        return self._operating_points

    # --- Immutable derived: manual lazy caching ---
    @property
    def lcm_period(self) -> float:
        """The least common multiple of all motion periods, ensuring all motions
        complete an integer number of cycles when cycle averaging forces and moments.

        Using the LCM ensures that when cycle averaging forces and moments, we capture a
        complete cycle of all motions, not just the longest one. For example, if one
        motion has a period of 2.0 s and another has a period of 3.0 s, the LCM is 6.0,
        which contains exactly 3 cycles of the first motion and 2 cycles of the second.

        :return: The LCM period in seconds. If all the motion is static, this will be
            0.0.
        """
        if self._lcm_period is None:
            # Collect all periods from AirplaneMovements
            all_periods: list[float] = []
            for airplane_movement in self._airplane_movements:
                all_periods.extend(airplane_movement.all_periods)

            # Add the OperatingPointMovement period
            all_periods.append(self._operating_point_movement.max_period)

            self._lcm_period = _lcm_multiple(all_periods)
        return self._lcm_period

    @property
    def max_period(self) -> float:
        """The longest period of motion of Movement's sub movement objects, the
        motion(s) of its sub sub movement object(s), and the motions of its sub sub sub
        movement objects.

        Note: For cycle averaging calculations, lcm_period should be used instead of
        max_period to ensure all motions complete an integer number of cycles.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        if self._max_period is None:
            # Iterate through the AirplaneMovements and find the one with the largest
            # max period.
            airplane_movement_max_periods = []
            for airplane_movement in self._airplane_movements:
                airplane_movement_max_periods.append(airplane_movement.max_period)
            max_airplane_period = max(airplane_movement_max_periods)

            # The global max period is the maximum of the max AirplaneMovement period
            # and the OperatingPointMovement max period.
            self._max_period = max(
                max_airplane_period,
                self._operating_point_movement.max_period,
            )
        return self._max_period

    @property
    def min_period(self) -> float:
        """The shortest non zero period of motion of Movement's sub movement objects,
        the motion(s) of its sub sub movement object(s), and the motions of its sub sub
        sub movement objects.

        :return: The shortest non zero period in seconds. If all the motion is static,
            this will be 0.0.
        """
        if self._min_period is None:
            # Collect all periods from AirplaneMovements.
            all_periods: list[float] = []
            for airplane_movement in self._airplane_movements:
                all_periods.extend(airplane_movement.all_periods)

            # Add the OperatingPointMovement period.
            op_period = self._operating_point_movement.max_period
            if op_period != 0.0:
                all_periods.append(op_period)

            # Filter out zero periods and find the minimum.
            non_zero_periods = [p for p in all_periods if p != 0.0]
            if not non_zero_periods:
                self._min_period = 0.0
            else:
                self._min_period = min(non_zero_periods)
        return self._min_period

    @property
    def static(self) -> bool:
        """Flags if the Movement's sub movement objects, its sub sub movement object(s),
        and its sub sub sub movement objects all represent no motion.

        :return: True if Movement's sub movement objects, its sub sub movement
            object(s), and its sub sub sub movement objects all represent no motion.
            False otherwise.
        """
        if self._static is None:
            self._static = self.max_period == 0
        return self._static


def _compute_wake_area_mismatch(
    delta_time: float,
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
) -> float:
    """Computes the average area mismatch between wake and bound RingVortices.

    Creates a temporary Problem and solver, steps through some number of time steps
    (geometry only, no aerodynamic solve), and computes the average area mismatch at
    each step.

    The area mismatch metric measures how well the wake RingVortex sizing matches the
    bound RingVortex sizing. A lower value indicates better matching.

    The number of time steps checked is picked to capture the full range of differences
    in areas for the wake and bound RingVortex child parent pairs. For static cases,
    this is just a single time step. For non static cases, it is enough time steps to
    cover one full maximum length period of motion.

    :param delta_time: The delta_time value to test. It must be a positive float. Its
        units are in seconds.
    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :return: The average area mismatch. The absolute percent error between the area of
        shed wake RingVortices and the area of their parent bound RingVortices (at time
        step where they were shed). Averaged across all time steps and all pairs of
        child and parent RingVortices. A lower value indicates better matching.
    """
    from .. import problems, unsteady_ring_vortex_lattice_method

    # Deep copy the movement objects to avoid mutating originals during optimization.
    airplane_movements_copy = copy.deepcopy(airplane_movements)
    operating_point_movement_copy = copy.deepcopy(operating_point_movement)

    max_airplane_movement_period = 0.0
    for airplane_movement in airplane_movements_copy:
        max_airplane_movement_period = max(
            airplane_movement.max_period, max_airplane_movement_period
        )

    max_period = max(
        max_airplane_movement_period, operating_point_movement_copy.max_period
    )

    # Calculate the number of steps to traverse the max period (or just a single step if
    # there is no movement).
    num_steps = 1
    if max_period > 0.0:
        num_steps = math.ceil(max_period / delta_time)

    # Create a temporary Movement with the trial delta_time.
    temp_movement = Movement(
        airplane_movements=airplane_movements_copy,
        operating_point_movement=operating_point_movement_copy,
        delta_time=delta_time,
        num_steps=num_steps,
    )

    # Create an UnsteadyProblem and solver.
    temp_problem = problems.UnsteadyProblem(movement=temp_movement)
    temp_solver = (
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            temp_problem
        )
    )

    # Accumulate area mismatch across all steps > 0.
    total_mismatch = 0.0
    num_comparisons = 0

    # Step through the simulation using geometry only initialization.
    for step in range(num_steps):
        temp_solver.initialize_step_geometry(step)

        # At step > 0, compare wake first row RingVortex areas (current step)
        # to bound trailing edge RingVortex areas (previous step).
        if step > 0:
            # Get the current Airplanes (at step) for wake RingVortices.
            current_airplanes = temp_solver.steady_problems[step].airplanes
            # Get the previous Airplanes (at step - 1) for bound RingVortices.
            previous_airplanes = temp_solver.steady_problems[step - 1].airplanes

            for airplane_id, airplane in enumerate(current_airplanes):
                previous_airplane = previous_airplanes[airplane_id]

                for wing_id, wing in enumerate(airplane.wings):
                    previous_wing = previous_airplane.wings[wing_id]

                    # Get the wake RingVortices (first row, chordwise index 0).
                    wake_ring_vortices = wing.wake_ring_vortices

                    assert wake_ring_vortices is not None

                    # First row of wake is at chordwise index 0.
                    num_spanwise = wake_ring_vortices.shape[1]

                    # Get the trailing edge bound RingVortices from previous step.
                    previous_panels = previous_wing.panels
                    if previous_panels is None:
                        continue

                    num_chordwise_panels = previous_wing.num_chordwise_panels
                    trailing_edge_chordwise_index = num_chordwise_panels - 1

                    for spanwise_id in range(num_spanwise):
                        # Get wake RingVortex area (first row, current step).
                        wake_rv: _vortices.ring_vortex.RingVortex = wake_ring_vortices[
                            0, spanwise_id
                        ]
                        wake_area = wake_rv.area

                        # Get bound trailing edge RingVortex area (previous step).
                        trailing_edge_panel = previous_panels[
                            trailing_edge_chordwise_index, spanwise_id
                        ]
                        _bound_rv = trailing_edge_panel.ring_vortex

                        assert _bound_rv is not None
                        bound_rv: _vortices.ring_vortex.RingVortex = _bound_rv

                        bound_area = bound_rv.area

                        # Accumulate the absolute percent area difference.
                        epsilon = 1e-12
                        if abs(bound_area) > epsilon:
                            total_mismatch += abs(wake_area - bound_area) / bound_area
                            num_comparisons += 1

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

    Optimizes delta_time to minimize the area mismatch between wake RingVortices and
    their parent bound trailing edge RingVortices. This produces better results at high
    Strouhal numbers where motion induced velocity is significant.

    For non static Movements, the optimization uses a brute force search over integer
    num_steps_per_lcm_cycle values (from 0.5x to 2x the initial estimate). This ensures
    the resulting delta_time divides the LCM period evenly.

    For static Movements, the optimization uses scipy.optimize.minimize_scalar with
    early termination if the mismatch falls below the specified cutoff value.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate from the fast calculation. It must
        be a positive float. Its units are in seconds.
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
        lcm_period = _lcm_multiple(non_zero_periods)
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

    best_num_steps = min_num_steps
    best_mismatch = float("inf")

    for num_steps in range(min_num_steps, max_num_steps + 1):
        delta_time = lcm_period / num_steps
        mismatch = _compute_wake_area_mismatch(
            delta_time, airplane_movements, operating_point_movement
        )

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

    Estimates the delta_time that produces wake RingVortices with roughly the same chord
    length as the bound trailing edge RingVortices, accounting for both freestream and
    geometry motion velocities. This is faster than _optimize_delta_time but may be
    slightly less accurate.

    The algorithm works by: (1) computing a very small preliminary delta_time as the
    minimum motion period divided by 100 (capped by a maximum of 1000 total time steps
    for one cycle of the motion's LCM period), (2) creating a temporary Movement and
    UnsteadyProblem to generate Airplane geometry at each time step (with all Panel
    coordinates transformed into the first Airplane's geometry axes), (3) measuring the
    average wake displacement per time step for each Wing's trailing edge Panels
    (combining freestream velocity and geometry motion, both in the first Airplane's
    geometry axes), and (4) choosing a delta_time such that each wake RingVortex has
    approximately the same chord as its parent bound RingVortex (averaged across the one
    LCM period).

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate from the fast calculation. It must
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

    lcm_period = _lcm_multiple(non_zero_periods)

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
            # the target chord length for wake RingVortices. We use the trailing
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
            # RingVortex chord equals the trailing edge bound Panel chord.
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
    # matching wake RingVortex area to bound RingVortex area requires a time step that
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
