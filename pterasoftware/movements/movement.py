"""Contains the Movement class.

**Contains the following classes:**

Movement: A class used to contain an UnsteadyProblem's movement.

**Contains the following functions:**

None
"""

from __future__ import annotations

import copy
import math

import scipy.optimize as sp_opt

from . import airplane_movement as airplane_movement_mod
from . import operating_point_movement as operating_point_movement_mod

from .. import _aerodynamics
from .. import _logging
from .. import _parameter_validation

movement_logger = _logging.get_logger("movements.movement")


class Movement:
    """A class used to contain an UnsteadyProblem's movement.

    **Contains the following methods:**

    lcm_period: The least common multiple of all motion periods, ensuring all motions
    complete an integer number of cycles when cycle averaging forces and moments.

    max_period: The longest period of motion of Movement's sub movement objects, the
    motion(s) of its sub sub movement object(s), and the motions of its sub sub sub
    movement objects.

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
            (default): Movement calculates a fast estimate based on freestream velocity
            alone. This works well when forward velocity dominates, but may give poor
            results at high Strouhal numbers where motion velocity is significant. The
            estimate is based on the first base Airplane's reference chord length, its
            first Wing's number of chordwise panels, and its base OperatingPoint's
            velocity. "optimize": Movement runs an iterative optimization to find the
            delta_time that minimizes the area mismatch between wake RingVortices and
            their parent bound trailing edge RingVortices. This is slower but produces
            better results at high Strouhal numbers. Positive number (int or float): Use
            the specified value directly. All values are converted internally to floats.
            The units are in seconds.
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
        self.airplane_movements = airplane_movements

        if not isinstance(
            operating_point_movement,
            operating_point_movement_mod.OperatingPointMovement,
        ):
            raise TypeError(
                "operating_point_movement must be an OperatingPointMovement."
            )
        self.operating_point_movement = operating_point_movement

        # Track whether optimization should run after the initial setup.
        _should_optimize_delta_time: bool = False

        if isinstance(delta_time, str):
            if delta_time != "optimize":
                raise ValueError('delta_time string must be "optimize".')
            _should_optimize_delta_time = True
            # Fall through to calculate initial estimate for optimization.
            delta_time = None

        if delta_time is not None:
            delta_time = _parameter_validation.number_in_range_return_float(
                delta_time, "delta_time", min_val=0.0, min_inclusive=False
            )
        else:
            # Calculate initial delta_time estimate based on freestream velocity.
            # This works well when forward velocity dominates, but may give poor
            # results at high Strouhal numbers where motion velocity is significant.
            delta_times = []
            for airplane_movement in self.airplane_movements:
                # TODO: Consider making this also average across each Airplane's Wings.
                # For a given Airplane, the ideal time step length is that which
                # sheds RingVortices off the first Wing that have roughly the same
                # chord length as the RingVortices on the first Wing. This is based
                # on the base Airplane's reference chord length, its first Wing's
                # number of chordwise panels, and its base OperatingPoint's velocity.
                c_ref = airplane_movement.base_airplane.c_ref
                assert c_ref is not None
                delta_times.append(
                    c_ref
                    / airplane_movement.base_airplane.wings[0].num_chordwise_panels
                    / operating_point_movement.base_operating_point.vCg__E
                )

            # Set the delta_time to be the average of the Airplanes' ideal delta times.
            delta_time = sum(delta_times) / len(delta_times)

        # Run delta_time optimization if requested.
        if _should_optimize_delta_time:
            delta_time = _optimize_delta_time(
                airplane_movements=self.airplane_movements,
                operating_point_movement=self.operating_point_movement,
                initial_delta_time=delta_time,
            )
        self.delta_time: float = delta_time

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
        self.num_cycles = num_cycles

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
        self.num_chords = num_chords

        if self.num_cycles is not None or self.num_chords is not None:
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
                for airplane_movement in self.airplane_movements:
                    c_ref = airplane_movement.base_airplane.c_ref
                    assert c_ref is not None
                    c_refs.append(c_ref)
                max_c_ref = max(c_refs)

                # Set the number of time steps such that the wake extends back by
                # some number of reference chord lengths.
                assert self.num_chords is not None
                wake_length = self.num_chords * max_c_ref
                distance_per_time_step = (
                    delta_time
                    * self.operating_point_movement.base_operating_point.vCg__E
                )
                num_steps = math.ceil(wake_length / distance_per_time_step)
            else:
                # Set the number of time steps such that the simulation runs for some
                # number of cycles of all motions. Use the LCM of all periods to ensure
                # each motion completes an integer number of cycles.
                assert self.num_cycles is not None
                num_steps = math.ceil(
                    self.num_cycles * self.lcm_period / self.delta_time
                )
        self.num_steps: int = num_steps

        # Generate a list of lists of Airplanes that are the steps through each
        # AirplaneMovement. The first index identifies the AirplaneMovement, and the
        # second index identifies the time step.
        self.airplanes = []
        for airplane_movement in self.airplane_movements:
            self.airplanes.append(
                airplane_movement.generate_airplanes(
                    num_steps=self.num_steps, delta_time=self.delta_time
                )
            )

        # Validate that all Wings maintain their symmetry type across all time steps.
        for airplane_movement_id, airplane_list in enumerate(self.airplanes):
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

        # Generate a lists of OperatingPoints that are the steps through the
        # OperatingPointMovement.
        self.operating_points = operating_point_movement.generate_operating_points(
            num_steps=self.num_steps, delta_time=self.delta_time
        )

    @staticmethod
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

    @staticmethod
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
            result = Movement._lcm(result, period)
        return result

    @property
    def lcm_period(self) -> float:
        """The least common multiple of all motion periods, ensuring all motions
        complete an integer number of cycles when cycle averaging forces and moments.

        Using the LCM ensures that when cycle-averaging forces and moments, we capture a
        complete cycle of all motions, not just the longest one. For example, if one
        motion has a period of 2.0 s and another has a period of 3.0 s, the LCM is 6.0,
        which contains exactly 3 cycles of the first motion and 2 cycles of the second.

        :return: The LCM period in seconds. If all the motion is static, this will be
            0.0.
        """
        # Collect all periods from AirplaneMovements
        all_periods = []
        for airplane_movement in self.airplane_movements:
            all_periods.extend(airplane_movement.all_periods)

        # Add the OperatingPointMovement period
        all_periods.append(self.operating_point_movement.max_period)

        return self._lcm_multiple(all_periods)

    @property
    def max_period(self) -> float:
        """The longest period of motion of Movement's sub movement objects, the
        motion(s) of its sub sub movement object(s), and the motions of its sub sub sub
        movement objects.

        Note: For cycle-averaging calculations, lcm_period should be used instead of
        max_period to ensure all motions complete an integer number of cycles.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        # Iterate through the AirplaneMovements and find the one with the largest max
        # period.
        airplane_movement_max_periods = []
        for airplane_movement in self.airplane_movements:
            airplane_movement_max_periods.append(airplane_movement.max_period)
        max_airplane_period = max(airplane_movement_max_periods)

        # The global max period is the maximum of the max AirplaneMovement period and
        # the OperatingPointMovement max period.
        return max(
            max_airplane_period,
            self.operating_point_movement.max_period,
        )

    @property
    def static(self) -> bool:
        """Flags if the Movement's sub movement objects, its sub sub movement object(s),
        and its sub sub sub movement objects all represent no motion.

        :return: True if Movement's sub movement objects, its sub sub movement
            object(s), and its sub sub sub movement objects all represent no motion.
            False otherwise.
        """
        return self.max_period == 0


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
    from .. import problems
    from .. import unsteady_ring_vortex_lattice_method

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
                        wake_rv: _aerodynamics.RingVortex = wake_ring_vortices[
                            0, spanwise_id
                        ]
                        wake_area = wake_rv.area

                        # Get bound trailing edge RingVortex area (previous step).
                        trailing_edge_panel = previous_panels[
                            trailing_edge_chordwise_index, spanwise_id
                        ]
                        _bound_rv = trailing_edge_panel.ring_vortex

                        assert _bound_rv is not None
                        bound_rv: _aerodynamics.RingVortex = _bound_rv

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
    """Finds an optimal delta_time using sp_opt.minimize_scalar.

    Optimizes delta_time to minimize the area mismatch between wake RingVortices and
    their parent bound trailing edge RingVortices. This produces better results at high
    Strouhal numbers where motion induced velocity is significant.

    The search terminates early if the mismatch falls below the specified cutoff value.
    Otherwise, it will return the locally minimum with an absolute convergence tolerance
    of 0.001.

    The optimization search is bounded within one order of magnitude, centered at the
    specified starting value.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The initial estimate from the fast calculation. It must
        be a positive float. Its units are in seconds.
    :param mismatch_cutoff: A positive float for the optimization's convergence
        threshold. When the average area mismatch (which is an absolute percent error)
        falls below this value, the search terminates early. The default is 0.01.
    :return: The optimized delta_time value. Its units are in seconds.
    :raises RuntimeError: If optimization fails to converge.
    """
    lower_bound = initial_delta_time / math.sqrt(10)
    upper_bound = initial_delta_time * math.sqrt(10)

    movement_logger.info("Starting delta_time optimization.")

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
