"""Contains the AirplaneMovement class.

Contains the following classes:
    AirplaneMovement: A class used to contain an Airplane's movement.

Contains the following functions:
    None
"""

from __future__ import annotations

from collections.abc import Sequence, Callable

import numpy as np

from . import _functions
from . import wing_movement as wing_movement_mod

from .. import _parameter_validation
from .. import geometry


class AirplaneMovement:
    """A class used to contain an Airplane's movement.

    Contains the following methods:
        generate_airplanes: Creates the Airplane at each time step, and returns them
        in a list.

        max_period: Defines a property for the longest period of AirplaneMovement's
        own motion, the motion(s) of its sub-movement object(s), and the motions of
        its sub-sub-movement objects.
    """

    def __init__(
        self,
        base_airplane: geometry.airplane.Airplane,
        wing_movements: list[wing_movement_mod.WingMovement],
        ampCg_E_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodCg_E_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingCg_E_CgP1: np.ndarray | Sequence[str] = ("sine", "sine", "sine"),
        phaseCg_E_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        ampAngles_E_to_B_izyx: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodAngles_E_to_B_izyx: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingAngles_E_to_B_izyx: np.ndarray | Sequence[str | Callable] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseAngles_E_to_B_izyx: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
    ) -> None:
        """The initialization method.

        :param base_airplane: The base Airplane from which the Airplane at each time
            step will be created.
        :param wing_movements: A list of the WingMovements associated with each of the
            base Airplane's Wings. It must have the same length as the base Airplane's
            list of Wings.
        :param ampCg_E_CgP1: An array-like object of 3 non-negative numbers (floats or
            ints) representing the amplitudes of the AirplaneMovement's changes in its
            Airplanes' Cg_E_CgP1 parameters. Can be a tuple, list, or ndarray. Values
            are converted to floats internally. Each amplitude must be low enough that
            it doesn't drive its base value out of the range of valid values. Otherwise,
            this AirplaneMovement will try to create Airplanes with invalid parameters
            values. Because the first Airplane's Cg_E_CgP1 parameter must be all zeros,
            this means that the first Airplane's ampCg_E_CgP1 parameter must also be all
            zeros. The units are in meters. The default is (0.0, 0.0, 0.0).
        :param periodCg_E_CgP1: An array-like object of 3 non-negative numbers (floats
            or ints) representing the periods of the AirplaneMovement's changes in its
            Airplanes' Cg_E_CgP1 parameters. Can be a tuple, list, or ndarray. Values
            are converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampCg_E_CgP1 is 0.0 and non-zero if not. The units
            are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingCg_E_CgP1: An array-like object of 3 strs or callables
            representing the spacing of the AirplaneMovement's change in its Airplanes'
            Cg_E_CgP1 parameters. Can be a tuple, list, or ndarray. Each element can be
            the str "sine", the str "uniform", or a callable custom spacing function.
            Custom spacing functions are for advanced users and must start at 0.0,
            return to 0.0 after one period of 2*pi radians, have amplitude of 1.0, be
            periodic, return finite values only, and accept a ndarray as input and
            return a ndarray of the same shape. Custom functions are scaled by
            ampCg_E_CgP1, shifted horizontally and vertically by phaseCg_E_CgP1 the base
            value, and have a period set by periodCg_E_CgP1. The default is ("sine",
            "sine", "sine").
        :param phaseCg_E_CgP1: An array-like object of 3 numbers (floats or ints)
            representing the phase offsets of the elements in the first time step's
            Airplane's Cg_E_CgP1 parameter relative to the base Airplane's Cg_E_CgP1
            parameter. Can be a tuple, list, or ndarray. Elements must lie in the range
            (-180.0, 180.0]. Each element must be 0.0 if the corresponding element in
            ampCg_E_CgP1 is 0.0 and non-zero if not. Values are converted to floats
            internally. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param ampAngles_E_to_B_izyx: An array-like object of 3 non-negative numbers
            (floats or ints) representing the amplitudes of the AirplaneMovement's
            changes in its Airplanes' angles_E_to_B_izyx parameters. Can be a tuple,
            list, or ndarray. Its elements must lie in the range [0.0, 360.0). Values
            are converted to floats internally. Each amplitude must be low enough that
            it doesn't drive its base value out of the range of valid values. Otherwise,
            this AirplaneMovement will try to create Airplanes with invalid parameters
            values. The units are in degrees. The default value is (0.0, 0.0, 0.0).
        :param periodAngles_E_to_B_izyx: An array-like  object of 3 non-negative
            numbers (floats or ints) representing the periods of the AirplaneMovement's
            changes in its Airplanes' angles_E_to_B_izyx parameters. Can be a tuple,
            list, or ndarray. Values are converted to floats internally. Each element
            must be 0.0 if the corresponding element in ampAngles_E_to_B_izyx is 0.0 and
            non-zero if not. The units are in seconds. The default value is (0.0, 0.0,
            0.0).
        :param spacingAngles_E_to_B_izyx: An array-like object of 3 strs or callables
            representing the spacing of the AirplaneMovement's change in its
            Airplanes' angles_E_to_B_izyx parameters. Can be a tuple, list, or ndarray.
            Each element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2*pi radians, have amplitude
            of 1.0, be periodic, return finite values only, and accept a ndarray as
            input and return a ndarray of the same shape. The custom function is scaled
            by ampAngles_E_to_B_izyx, shifted horizontally and vertically by
            phaseAngles_E_to_B_izyx and the base value, and have a period set by
            periodAngles_E_to_B_izyx. The default value is ("sine", "sine", "sine").
        :param phaseAngles_E_to_B_izyx: An array-like object of 3 numbers (floats or
            ints) representing the phase offsets of the elements in the first time
            step's Airplane's angles_E_to_B_izyx parameter relative to the base
            Airplane's angles_E_to_B_izyx parameter. Can be a tuple, list, or ndarray.
            Elements must lie in the range (-180.0, 180.0]. Each element must be 0.0 if
            the corresponding element in ampAngles_E_to_B_izyx is 0.0 and non-zero if
            not. Values are converted to floats internally. The units are in degrees.
            The default value is (0.0, 0.0, 0.0).
        :return: None
        """
        if not isinstance(base_airplane, geometry.airplane.Airplane):
            raise TypeError("base_airplane must be an Airplane.")
        self.base_airplane = base_airplane

        if not isinstance(wing_movements, list):
            raise TypeError("wing_movements must be a list.")
        if len(wing_movements) != len(self.base_airplane.wings):
            raise ValueError(
                "wing_movements must have the same length as base_airplane.wings."
            )
        for wing_movement in wing_movements:
            if not isinstance(wing_movement, wing_movement_mod.WingMovement):
                raise TypeError(
                    "Every element in wing_movements must be a WingMovement."
                )
        self.wing_movements = wing_movements

        ampCg_E_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            ampCg_E_CgP1, "ampCg_E_CgP1"
        )
        if not np.all(ampCg_E_CgP1 >= 0.0):
            raise ValueError("All elements in ampCg_E_CgP1 must be non-negative.")
        self.ampCg_E_CgP1 = ampCg_E_CgP1

        periodCg_E_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            periodCg_E_CgP1, "periodCg_E_CgP1"
        )
        if not np.all(periodCg_E_CgP1 >= 0.0):
            raise ValueError("All elements in periodCg_E_CgP1 must be non-negative.")
        for period_index, period in enumerate(periodCg_E_CgP1):
            amp = self.ampCg_E_CgP1[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampCg_E_CgP1 is 0.0, the corresponding element "
                    "in periodCg_E_CgP1 must be also be 0.0."
                )
        self.periodCg_E_CgP1 = periodCg_E_CgP1

        spacingCg_E_CgP1 = _parameter_validation.threeD_spacing_vectorLike_return_tuple(
            spacingCg_E_CgP1, "spacingCg_E_CgP1"
        )
        self.spacingCg_E_CgP1 = spacingCg_E_CgP1

        phaseCg_E_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            phaseCg_E_CgP1, "phaseCg_E_CgP1"
        )
        if not (np.all(phaseCg_E_CgP1 > -180.0) and np.all(phaseCg_E_CgP1 <= 180.0)):
            raise ValueError(
                "All elements in phaseCg_E_CgP1 must be in the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseCg_E_CgP1):
            amp = self.ampCg_E_CgP1[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampCg_E_CgP1 is 0.0, the corresponding element "
                    "in phaseCg_E_CgP1 must be also be 0.0."
                )
        self.phaseCg_E_CgP1 = phaseCg_E_CgP1

        ampAngles_E_to_B_izyx = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                ampAngles_E_to_B_izyx, "ampAngles_E_to_B_izyx"
            )
        )
        if not (
            np.all(ampAngles_E_to_B_izyx >= 0.0)
            and np.all(ampAngles_E_to_B_izyx < 360.0)
        ):
            raise ValueError(
                "All elements in ampAngles_E_to_B_izyx must be in the range [0.0, "
                "360.0)."
            )
        self.ampAngles_E_to_B_izyx = ampAngles_E_to_B_izyx

        periodAngles_E_to_B_izyx = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_E_to_B_izyx, "periodAngles_E_to_B_izyx"
            )
        )
        if not np.all(periodAngles_E_to_B_izyx >= 0.0):
            raise ValueError(
                "All elements in periodAngles_E_to_B_izyx must be non-negative."
            )
        for period_index, period in enumerate(periodAngles_E_to_B_izyx):
            amp = self.ampAngles_E_to_B_izyx[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_E_to_B_izyx is 0.0, the corresponding "
                    "element in periodAngles_E_to_B_izyx must be also be 0.0."
                )
        self.periodAngles_E_to_B_izyx = periodAngles_E_to_B_izyx

        spacingAngles_E_to_B_izyx = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingAngles_E_to_B_izyx,
                "spacingAngles_E_to_B_izyx",
            )
        )
        self.spacingAngles_E_to_B_izyx = spacingAngles_E_to_B_izyx

        phaseAngles_E_to_B_izyx = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                phaseAngles_E_to_B_izyx, "phaseAngles_E_to_B_izyx"
            )
        )
        if not (
            np.all(phaseAngles_E_to_B_izyx > -180.0)
            and np.all(phaseAngles_E_to_B_izyx <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_E_to_B_izyx must be in the range ("
                "-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_E_to_B_izyx):
            amp = self.ampAngles_E_to_B_izyx[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_E_to_B_izyx is 0.0, the corresponding "
                    "element in phaseAngles_E_to_B_izyx must be also be 0.0."
                )
        self.phaseAngles_E_to_B_izyx = phaseAngles_E_to_B_izyx

    def generate_airplanes(
        self, num_steps: int, delta_time: float | int
    ) -> list[geometry.airplane.Airplane]:
        """Creates the Airplane at each time step, and returns them in a list.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (float or int), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of Airplanes associated with this AirplaneMovement.
        """
        num_steps = _parameter_validation.positive_int_return_int(
            num_steps, "num_steps"
        )
        delta_time = _parameter_validation.positive_number_return_float(
            delta_time, "delta_time"
        )

        # Generate oscillating values for each dimension of Cg_E_CgP1.
        listCg_E_CgP1 = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingCg_E_CgP1[dim]
            if spacing == "sine":
                listCg_E_CgP1[dim, :] = _functions.oscillating_sinspaces(
                    amps=self.ampCg_E_CgP1[dim],
                    periods=self.periodCg_E_CgP1[dim],
                    phases=self.phaseCg_E_CgP1[dim],
                    bases=self.base_airplane.Cg_E_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                listCg_E_CgP1[dim, :] = _functions.oscillating_linspaces(
                    amps=self.ampCg_E_CgP1[dim],
                    periods=self.periodCg_E_CgP1[dim],
                    phases=self.phaseCg_E_CgP1[dim],
                    bases=self.base_airplane.Cg_E_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                listCg_E_CgP1[dim, :] = _functions.oscillating_customspaces(
                    amps=self.ampCg_E_CgP1[dim],
                    periods=self.periodCg_E_CgP1[dim],
                    phases=self.phaseCg_E_CgP1[dim],
                    bases=self.base_airplane.Cg_E_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

        # Generate oscillating values for each dimension of angles_E_to_B_izyx.
        listAngles_E_to_B_izyx = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingAngles_E_to_B_izyx[dim]
            if spacing == "sine":
                listAngles_E_to_B_izyx[dim, :] = _functions.oscillating_sinspaces(
                    amps=self.ampAngles_E_to_B_izyx[dim],
                    periods=self.periodAngles_E_to_B_izyx[dim],
                    phases=self.phaseAngles_E_to_B_izyx[dim],
                    bases=self.base_airplane.angles_E_to_B_izyx[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                listAngles_E_to_B_izyx[dim, :] = _functions.oscillating_linspaces(
                    amps=self.ampAngles_E_to_B_izyx[dim],
                    periods=self.periodAngles_E_to_B_izyx[dim],
                    phases=self.phaseAngles_E_to_B_izyx[dim],
                    bases=self.base_airplane.angles_E_to_B_izyx[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                listAngles_E_to_B_izyx[dim, :] = _functions.oscillating_customspaces(
                    amps=self.ampAngles_E_to_B_izyx[dim],
                    periods=self.periodAngles_E_to_B_izyx[dim],
                    phases=self.phaseAngles_E_to_B_izyx[dim],
                    bases=self.base_airplane.angles_E_to_B_izyx[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

        # Create an empty 2D ndarray that will hold each of the Airplane's Wing's vector
        # of Wings representing its changing state at each time step. The first index
        # denotes a particular base Wing, and the second index denotes the time step.
        wings = np.empty((len(self.wing_movements), num_steps), dtype=object)

        # Iterate through the WingMovements.
        for wing_movement_id, wing_movement in enumerate(self.wing_movements):

            # Generate this Wing's vector of Wings representing its changing state at
            # each time step.
            this_wings_list_of_wings = np.array(
                wing_movement.generate_wings(num_steps=num_steps, delta_time=delta_time)
            )

            # Add this vector the Airplane's 2D ndarray of Wings' Wings.
            wings[wing_movement_id, :] = this_wings_list_of_wings

        # Create an empty list to hold each time step's Airplane.
        airplanes = []

        # Get the non-changing Airplane attributes.
        this_name = self.base_airplane.name
        this_weight = self.base_airplane.weight

        # Iterate through the time steps.
        for step in range(num_steps):
            thisCg_E_CgP1 = listCg_E_CgP1[:, step]
            theseAngles_E_to_B_izyx = listAngles_E_to_B_izyx[:, step]
            these_wings = list(wings[:, step])

            # Make a new Airplane for this time step.
            this_airplane = geometry.airplane.Airplane(
                wings=these_wings,
                name=this_name,
                Cg_E_CgP1=thisCg_E_CgP1,
                angles_E_to_B_izyx=theseAngles_E_to_B_izyx,
                weight=this_weight,
            )

            # Add this new Airplane to the list of Airplanes.
            airplanes.append(this_airplane)

        return airplanes

    @property
    def max_period(self) -> float:
        """Defines a property for the longest period of AirplaneMovement's own motion,
        the motion(s) of its sub-movement object(s), and the motions of its
        sub-sub-movement objects.

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
                np.max(self.periodCg_E_CgP1),
                np.max(self.periodAngles_E_to_B_izyx),
            )
        )
