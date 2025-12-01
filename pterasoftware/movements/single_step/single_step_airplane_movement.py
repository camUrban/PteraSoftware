import numpy as np

from ..._parameter_validation import (
    threeD_number_vectorLike_return_float,
    threeD_spacing_vectorLike_return_tuple,
    positive_number_return_float,
    positive_int_return_int,
)

from .._functions import (
    oscillating_sinspaces, 
    oscillating_linspaces, 
    oscillating_customspaces
)

from ... import geometry


class SingleStepAirplaneMovement:
    def __init__(
        self,
        single_step_wing_movements,
        ampCg_E_CgP1=(0.0, 0.0, 0.0),
        periodCg_E_CgP1=(0.0, 0.0, 0.0),
        spacingCg_E_CgP1=("sine", "sine", "sine"),
        phaseCg_E_CgP1=(0.0, 0.0, 0.0),
        ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
        phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
    ):

        """
        :param wing_movements: list of SingleStepWingMovement

            This is a list of the WingMovement associated with each of the base
            Airplane's Wings. It must have the same length as the base Airplane's
            list of Wings.

        :param ampCg_E_CgP1: array-like of 3 numbers, optional

            The amplitudes of the AirplaneMovement's changes in its Airplanes'
            Cg_E_CgP1 parameters. Can be a tuple, list, or numpy array of
            non-negative numbers (int or float). Also, each amplitude must be low
            enough that it doesn't drive its base value out of the range of valid
            values. Otherwise, this AirplaneMovement will try to create Airplanes
            with invalid parameters values. Because the first Airplane's Cg_E_CgP1
            parameter must be all zeros, this means that the first Airplane's
            ampCg_E_CgP1 parameter must also be all zeros. Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). The units are in
            meters.

        :param periodCg_E_CgP1: array-like of 3 numbers, optional

            The periods of the AirplaneMovement's changes in its Airplanes' Cg_E_CgP1
            parameters. Can be a tuple, list, or numpy array of non-negative numbers
            (int or float). Values are converted to floats internally. The default
            value is (0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding
            element in ampCg_E_CgP1 is 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingCg_E_CgP1: array-like of 3 strs or callables, optional

            The value determines the spacing of the AirplaneMovement's change in its
            Airplanes' Cg_E_CgP1 parameters. Can be a tuple, list, or numpy array.
            Each element can be the string "sine", the string "uniform",
            or a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0, return to 0 after one period of 2*pi
            radians, have amplitude of 1, be periodic, return finite values only,
            and accept a ndarray as input and return a ndarray of the same shape. The
            custom function is scaled by ampCg_E_CgP1, shifted horizontally by
            phaseCg_E_CgP1, and vertically by the base value, with the period
            controlled by periodCg_E_CgP1. The default value is ("sine", "sine",
            "sine").

        :param phaseCg_E_CgP1: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Airplane's
            Cg_E_CgP1 parameter relative to the base Airplane's Cg_E_CgP1 parameter.
            Can be a tuple, list, or numpy array of non-negative numbers (int or
            float) in the range (-180.0, 180.0]. Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). Each element must be
            0.0 if the corresponding element in ampCg_E_CgP1 is 0.0 and non-zero if
            not. The units are in degrees.

        :param ampAngles_E_to_B_izyx: array-like of 3 numbers, optional

            The amplitudes of the AirplaneMovement's changes in its Airplanes'
            angles_E_to_B_izyx parameters. Can be a tuple, list, or numpy array of
            numbers (int or float) in the range [0.0, 360.0). Also, each amplitude
            must be low enough that it doesn't drive its base value out of the range
            of valid values. Otherwise, this AirplaneMovement will try to create
            Airplanes with invalid parameters values. Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). The units are in degrees.

        :param periodAngles_E_to_B_izyx: array-like of 3 numbers, optional

            The periods of the AirplaneMovement's changes in its Airplanes'
            angles_E_to_B_izyx parameters. Can be a tuple, list, or numpy array of
            non-negative numbers (int or float). Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). Each element must be
            0.0 if the corresponding element in ampAngles_E_to_B_izyx is 0.0 and
            non-zero if not. The units are in seconds.

        :param spacingAngles_E_to_B_izyx: array-like of 3 strs or callables, optional

            The value determines the spacing of the AirplaneMovement's change in its
            Airplanes' angles_E_to_B_izyx parameters. Can be a tuple, list, or numpy
            array. Each element can be the string "sine", the string "uniform",
            or a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0, return to 0 after one period of 2*pi
            radians, have amplitude of 1, be periodic, return finite values only,
            and accept a ndarray as input and return a ndarray of the same shape. The
            custom function is scaled by ampAngles_E_to_B_izyx, shifted horizontally
            by phaseAngles_E_to_B_izyx, and vertically by the base value, with the
            period controlled by periodAngles_E_to_B_izyx. The default value is (
            "sine", "sine", "sine").

        :param phaseAngles_E_to_B_izyx: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Airplane's
            angles_E_to_B_izyx parameter relative to the base Airplane's
            angles_E_to_B_izyx parameter. Can be a tuple, list, or numpy array of
            numbers (int or float) in the range (-180.0, 180.0]. Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampAngles_E_to_B_izyx is 0.0
            and non-zero if not. The units are in degrees.
        """

        self.wing_movements = single_step_wing_movements

        ampCg_E_CgP1 = threeD_number_vectorLike_return_float(
            ampCg_E_CgP1, "ampCg_E_CgP1"
        )
        if not np.all(ampCg_E_CgP1 >= 0.0):
            raise ValueError("All elements in ampCg_E_CgP1 must be non-negative.")
        self.ampCg_E_CgP1 = ampCg_E_CgP1

        periodCg_E_CgP1 = threeD_number_vectorLike_return_float(
            periodCg_E_CgP1, "periodCg_E_CgP1"
        )
        if not np.all(periodCg_E_CgP1 >= 0.0):
            raise ValueError("All elements in periodCg_E_CgP1 must be non-negative.")
        for period_index, period in enumerate(periodCg_E_CgP1):
            amp = self.ampCg_E_CgP1[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampCg_E_CgP1 is 0.0, the corresponding element in periodCg_E_CgP1 must be also be 0.0."
                )
        self.periodCg_E_CgP1 = periodCg_E_CgP1

        spacingCg_E_CgP1 = threeD_spacing_vectorLike_return_tuple(
            spacingCg_E_CgP1, "spacingCg_E_CgP1"
        )
        self.spacingCg_E_CgP1 = spacingCg_E_CgP1

        phaseCg_E_CgP1 = threeD_number_vectorLike_return_float(
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
                    "If an element in ampCg_E_CgP1 is 0.0, the corresponding element in phaseCg_E_CgP1 must be also be 0.0."
                )
        self.phaseCg_E_CgP1 = phaseCg_E_CgP1

        ampAngles_E_to_B_izyx = (
            threeD_number_vectorLike_return_float(
                ampAngles_E_to_B_izyx, "ampAngles_E_to_B_izyx"
            )
        )
        if not (
            np.all(ampAngles_E_to_B_izyx >= 0.0)
            and np.all(ampAngles_E_to_B_izyx < 360.0)
        ):
            raise ValueError(
                "All elements in ampAngles_E_to_B_izyx must be in the range [0.0, 360.0)."
            )
        self.ampAngles_E_to_B_izyx = ampAngles_E_to_B_izyx

        periodAngles_E_to_B_izyx = (
            threeD_number_vectorLike_return_float(
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
                    "If an element in ampAngles_E_to_B_izyx is 0.0, the corresponding element in periodAngles_E_to_B_izyx must be also be 0.0."
                )
        self.periodAngles_E_to_B_izyx = periodAngles_E_to_B_izyx

        spacingAngles_E_to_B_izyx = (
            threeD_spacing_vectorLike_return_tuple(
                spacingAngles_E_to_B_izyx,
                "spacingAngles_E_to_B_izyx",
            )
        )
        self.spacingAngles_E_to_B_izyx = spacingAngles_E_to_B_izyx

        phaseAngles_E_to_B_izyx = (
            threeD_number_vectorLike_return_float(
                phaseAngles_E_to_B_izyx, "phaseAngles_E_to_B_izyx"
            )
        )
        if not (
            np.all(phaseAngles_E_to_B_izyx > -180.0)
            and np.all(phaseAngles_E_to_B_izyx <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_E_to_B_izyx must be in the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_E_to_B_izyx):
            amp = self.ampAngles_E_to_B_izyx[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_E_to_B_izyx is 0.0, the corresponding element in phaseAngles_E_to_B_izyx must be also be 0.0."
                )
        self.phaseAngles_E_to_B_izyx = phaseAngles_E_to_B_izyx

        self.listCg_E_CgP1 = None
        self.listAngles_E_to_B_izyx = None

    def generate_next_airplane(
        self, base_airplane, delta_time, num_steps, step
    ):
        """Creates the Airplane at the next timestep

                :param delta_time: number

                    This is the time between each time step. It must be a positive number (
                    int or float), and will be converted internally to a float. The units are
                    in seconds.
        angles_E_to_B_izyx
                :return: Airplanes

                    This is the Airplanes associated with this AirplaneMovement and deformation.
        """
        num_steps = positive_int_return_int(
            num_steps, "num_steps"
        )
        delta_time = positive_number_return_float(
            delta_time, "delta_time"
        )
        # Generate oscillating values for each dimension of Cg_E_CgP1.
        if self.listCg_E_CgP1 is None:
            self._initialize_oscilating_dimensions(delta_time, num_steps, base_airplane)

        # Generate oscillating values for each dimension of angles_E_to_B_izyx.
        if self.listAngles_E_to_B_izyx is None:
            self._initialize_oscilating_angles(delta_time, num_steps, base_airplane)

        wings = []

        # Iterate through the WingMovements.
        for wing_movement_id, wing_movement in enumerate(self.wing_movements):

            # Add this vector the Airplane's 2D ndarray of Wings' Wings.
            wings.append(
                wing_movement.generate_next_wing(
                    base_wing=base_airplane.wings[wing_movement_id],
                    delta_time=delta_time,
                    num_steps=num_steps,
                    step=step,
                )
            )

        # Get the non-changing Airplane attributes.
        this_name = base_airplane.name
        this_weight = base_airplane.weight

        # the 1 is for not the base step, but 1 step deep
        thisCg_E_CgP1 = self.listCg_E_CgP1[:, step]
        theseAngles_E_to_B_izyx = self.listAngles_E_to_B_izyx[:, step]

        # Make a new Airplane for this time step.
        this_airplane = geometry.airplane.Airplane(
            wings=wings,
            name=this_name,
            Cg_E_CgP1=thisCg_E_CgP1,
            angles_E_to_B_izyx=theseAngles_E_to_B_izyx,
            weight=this_weight,
        )

        return this_airplane

    def _initialize_oscilating_dimensions(self, delta_time, num_steps, base_airplane):
        """Initializes the oscillating dimensions for Cg_E_CgP1 and angles_E_to_B_izyx.
        :param delta_time: number

            This is the time between each time step. It must be a positive number (
            int or float), and will be converted internally to a float. The units are
            in seconds.
        :param num_steps: int

            This is the number of time steps in this movement. It must be a positive
            int.
        """
        self.listCg_E_CgP1 = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingCg_E_CgP1[dim]
            if spacing == "sine":
                self.listCg_E_CgP1[dim, :] = oscillating_sinspaces(
                    amps=self.ampCg_E_CgP1[dim],
                    periods=self.periodCg_E_CgP1[dim],
                    phases=self.phaseCg_E_CgP1[dim],
                    bases=base_airplane.Cg_E_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listCg_E_CgP1[dim, :] = oscillating_linspaces(
                    amps=self.ampCg_E_CgP1[dim],
                    periods=self.periodCg_E_CgP1[dim],
                    phases=self.phaseCg_E_CgP1[dim],
                    bases=base_airplane.Cg_E_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listCg_E_CgP1[dim, :] = oscillating_customspaces(
                    amps=self.ampCg_E_CgP1[dim],
                    periods=self.periodCg_E_CgP1[dim],
                    phases=self.phaseCg_E_CgP1[dim],
                    bases=base_airplane.Cg_E_CgP1[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

    def _initialize_oscilating_angles(self, delta_time, num_steps, base_airplane):
        """Initializes the oscillating angles for angles_E_to_B_izyx.
        :param delta_time: number

            This is the time between each time step. It must be a positive number (
            int or float), and will be converted internally to a float. The units are
            in seconds.
        :param num_steps: int 
            This is the number of time steps in this movement. It must be a positive
            int.
        :param base_airplane: Airplane

            This is the base Airplane from which the AirplaneMovement will generate
            its Airplanes.
        """
        self.listAngles_E_to_B_izyx = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingAngles_E_to_B_izyx[dim]
            if spacing == "sine":
                self.listAngles_E_to_B_izyx[dim, :] = oscillating_sinspaces(
                    amps=self.ampAngles_E_to_B_izyx[dim],
                    periods=self.periodAngles_E_to_B_izyx[dim],
                    phases=self.phaseAngles_E_to_B_izyx[dim],
                    bases=base_airplane.angles_E_to_B_izyx[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listAngles_E_to_B_izyx[dim, :] = oscillating_linspaces(
                    amps=self.ampAngles_E_to_B_izyx[dim],
                    periods=self.periodAngles_E_to_B_izyx[dim],
                    phases=self.phaseAngles_E_to_B_izyx[dim],
                    bases=base_airplane.angles_E_to_B_izyx[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listAngles_E_to_B_izyx[dim, :] = oscillating_customspaces(
                    amps=self.ampAngles_E_to_B_izyx[dim],
                    periods=self.periodAngles_E_to_B_izyx[dim],
                    phases=self.phaseAngles_E_to_B_izyx[dim],
                    bases=base_airplane.angles_E_to_B_izyx[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")
