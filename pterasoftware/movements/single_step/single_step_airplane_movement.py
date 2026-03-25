"""Contains the SingleStepAirplaneMovement class.

**Contains the following classes:**

SingleStepAirplaneMovement: A single step variant of AirplaneMovement that generates
one Airplane per time step instead of all at once. Uses composition to wrap an
AirplaneMovement.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np


from ..airplane_movement import AirplaneMovement
from ... import geometry

from ..._parameter_validation import (
    int_in_range_return_int,
    number_in_range_return_float,
)

from .._functions import (
    oscillating_sinspaces,
    oscillating_linspaces,
    oscillating_customspaces,
)


class SingleStepAirplaneMovement:
    """A single step variant of AirplaneMovement for coupled simulations.

    This class wraps an AirplaneMovement via composition and generates one Airplane per
    time step (via generate_next_airplane) rather than generating all Airplanes at once.
    The composed AirplaneMovement is accessible via corresponding_airplane_movement.

    **Contains the following methods:**

    all_periods: All unique non zero periods from this
    SingleStepAirplaneMovement, its SingleStepWingMovements, and their
    SingleStepWingCrossSectionMovements.

    generate_next_airplane: Creates the Airplane at a single time step.

    max_period: The longest period of this SingleStepAirplaneMovement's own motion, the
    motion(s) of its sub movement object(s), and the motions of its sub sub movement
    objects.
    """

    __slots__ = (
        "wing_movements",
        "ampCg_GP1_CgP1",
        "periodCg_GP1_CgP1",
        "spacingCg_GP1_CgP1",
        "phaseCg_GP1_CgP1",
        "listCg_GP1_CgP1",
        "corresponding_airplane_movement",
    )

    def __init__(
        self,
        single_step_wing_movements,
        base_airplane,
        ampCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1: (
            np.ndarray | Sequence[str | Callable[[np.ndarray], np.ndarray]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
    ) -> None:
        """The initialization method.

        :param single_step_wing_movements: A list of the SingleStepWingMovements
            associated with each of the base Airplane's Wings. It must have the same
            length as the base Airplane's list of Wings.
        :param base_airplane: The base Airplane from which the Airplane at each time
            step will be created.
        :param ampCg_GP1_CgP1: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            SingleStepAirplaneMovement's changes in its Airplanes' Cg_GP1_CgP1
            parameters. Can be a tuple, list, or ndarray. Values are converted to
            floats internally. Each amplitude must be low enough that it doesn't drive
            its base value out of the range of valid values. Otherwise, this
            SingleStepAirplaneMovement will try to create Airplanes with invalid
            parameter values. Because the first Airplane's Cg_GP1_CgP1 parameter must
            be all zeros, this means that the first Airplane's ampCg_GP1_CgP1 parameter
            must also be all zeros. The units are in meters. The default is (0.0, 0.0,
            0.0).
        :param periodCg_GP1_CgP1: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            SingleStepAirplaneMovement's changes in its Airplanes' Cg_GP1_CgP1
            parameters. Can be a tuple, list, or ndarray. Values are converted to
            floats internally. Each element must be 0.0 if the corresponding element in
            ampCg_GP1_CgP1 is 0.0 and non zero if not. The units are in seconds. The
            default is (0.0, 0.0, 0.0).
        :param spacingCg_GP1_CgP1: An array-like object of strs or callables with shape
            (3,) representing the spacing of the SingleStepAirplaneMovement's changes
            in its Airplanes' Cg_GP1_CgP1 parameters. Can be a tuple, list, or ndarray.
            Each element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2*pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a
            ndarray as input and return a ndarray of the same shape. Custom functions
            are scaled by ampCg_GP1_CgP1, shifted horizontally and vertically by
            phaseCg_GP1_CgP1 and the base value, and have a period set by
            periodCg_GP1_CgP1. The default is ("sine", "sine", "sine").
        :param phaseCg_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's Airplane's Cg_GP1_CgP1 parameter relative to the base Airplane's
            Cg_GP1_CgP1 parameter. Can be a tuple, list, or ndarray. Elements must lie
            in the range (-180.0, 180.0]. Each element must be 0.0 if the corresponding
            element in ampCg_GP1_CgP1 is 0.0 and non zero if not. Values are converted
            to floats internally. The units are in degrees. The default is (0.0, 0.0,
            0.0).
        :return: None
        """
        self.wing_movements = single_step_wing_movements

        # Create the corresponding AirplaneMovement, which validates all oscillation
        # parameters and is also needed by coupled unsteady problems.
        corresponding_wing_movements = [
            wm.corresponding_wing_movement for wm in self.wing_movements
        ]
        self.corresponding_airplane_movement = AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=corresponding_wing_movements,
            ampCg_GP1_CgP1=ampCg_GP1_CgP1,
            periodCg_GP1_CgP1=periodCg_GP1_CgP1,
            spacingCg_GP1_CgP1=spacingCg_GP1_CgP1,
            phaseCg_GP1_CgP1=phaseCg_GP1_CgP1,
        )

        # Copy validated attributes from the corresponding AirplaneMovement.
        self.ampCg_GP1_CgP1 = self.corresponding_airplane_movement.ampCg_GP1_CgP1
        self.periodCg_GP1_CgP1 = self.corresponding_airplane_movement.periodCg_GP1_CgP1
        self.spacingCg_GP1_CgP1 = self.corresponding_airplane_movement.spacingCg_GP1_CgP1
        self.phaseCg_GP1_CgP1 = self.corresponding_airplane_movement.phaseCg_GP1_CgP1

        self.listCg_GP1_CgP1 = None

    @property
    def all_periods(self) -> list[float]:
        """All unique non zero periods from this SingleStepAirplaneMovement, its
        SingleStepWingMovements, and their SingleStepWingCrossSectionMovements.

        :return: A list of all unique non zero periods in seconds. If all motion is
            static, this will be an empty list.
        """
        periods = []

        # Collect all periods from WingMovement(s).
        for wing_movement in self.wing_movements:
            periods.extend(wing_movement.all_periods)

        # Collect all periods from AirplaneMovement's own motion.
        for period in self.periodCg_GP1_CgP1:
            if period > 0.0:
                periods.append(float(period))
        return periods

    def generate_next_airplane(
        self, base_airplane, delta_time: float | int, num_steps: int, step: int, deformation_matrices,
    ) -> geometry.airplane.Airplane:
        """Creates the Airplane at a single time step.

        :param base_airplane: The base Airplane from which the new Airplane will be
            created.
        :param delta_time: The time between each time step. It must be a positive number
            (float or int), and will be converted internally to a float. The units are
            in seconds.
        :param num_steps: The total number of time steps in this movement. It must be a
            positive int.
        :param step: The index of the current time step.
        :param deformation_matrices: Deformation matrices to apply to the Wings, or
            None.
        :return: The Airplane at the specified time step.
        """
        num_steps = int_in_range_return_int(
            num_steps,
            "num_steps",
            min_val=1,
            min_inclusive=True,
        )
        delta_time = number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )

        # Generate oscillating values for each dimension of Cg_GP1_CgP1.
        if self.listCg_GP1_CgP1 is None:
            self._initialize_oscillating_dimensions(delta_time, num_steps, base_airplane)

        wings = []

        # Iterate through the WingMovements.
        for wing_movement_id, wing_movement in enumerate(self.wing_movements):

            wings.append(
                wing_movement.generate_next_wing(
                    base_wing=base_airplane.wings[wing_movement_id],
                    delta_time=delta_time,
                    num_steps=num_steps,
                    step=step,
                    deformation_matrices=deformation_matrices,
                )
            )

        # Get the non changing Airplane attributes.
        this_name = base_airplane.name
        this_weight = base_airplane.weight

        thisCg_GP1_CgP1 = self.listCg_GP1_CgP1[:, step]

        # Make a new Airplane for this time step.
        this_airplane = geometry.airplane.Airplane(
            wings=wings,
            name=this_name,
            Cg_GP1_CgP1=thisCg_GP1_CgP1,
            weight=this_weight,
        )

        return this_airplane

    def _initialize_oscillating_dimensions(self, delta_time, num_steps, base_airplane):
        """Pre computes the oscillating Cg_GP1_CgP1 values for all time steps.

        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps.
        :param base_airplane: The base Airplane providing the base Cg_GP1_CgP1 values.
        :return: None
        """
        self.listCg_GP1_CgP1 = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingCg_GP1_CgP1[dim]
            if spacing == "sine":
                self.listCg_GP1_CgP1[dim, :] = oscillating_sinspaces(
                    amps=self.ampCg_GP1_CgP1[dim],
                    periods=self.periodCg_GP1_CgP1[dim],
                    phases=self.phaseCg_GP1_CgP1[dim],
                    bases=base_airplane.Cg_GP1_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listCg_GP1_CgP1[dim, :] = oscillating_linspaces(
                    amps=self.ampCg_GP1_CgP1[dim],
                    periods=self.periodCg_GP1_CgP1[dim],
                    phases=self.phaseCg_GP1_CgP1[dim],
                    bases=base_airplane.Cg_GP1_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listCg_GP1_CgP1[dim, :] = oscillating_customspaces(
                    amps=self.ampCg_GP1_CgP1[dim],
                    periods=self.periodCg_GP1_CgP1[dim],
                    phases=self.phaseCg_GP1_CgP1[dim],
                    bases=base_airplane.Cg_GP1_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

    @property
    def max_period(self) -> float:
        """The longest period of this SingleStepAirplaneMovement's own motion, the
        motion(s) of its sub movement object(s), and the motions of its sub sub movement
        objects.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        wing_movement_max_periods = []
        for wing_movement in self.wing_movements:
            wing_movement_max_periods.append(wing_movement.max_period)
        max_wing_movement_period = max(wing_movement_max_periods)

        return float(
            max(
                max_wing_movement_period,
                np.max(self.periodCg_GP1_CgP1),
            )
        )
