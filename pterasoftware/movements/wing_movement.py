"""This module contains the WingMovement class.

This module contains the following classes:
    WingMovement: This is a class used to contain the Wing movements.

This module contains the following functions:
    None
"""

import numpy as np

from . import _functions
from .wing_cross_section_movement import WingCrossSectionMovement

from .. import geometry
from .. import _parameter_validation


class WingMovement:
    """This is a class used to contain the Wing movements.

    Note: Wings cannot undergo motion that causes them to switch symmetry types. A
    transition between types could change the number of Wings and the panel
    structure, which is incompatible with the unsteady solver. This happens when a
    WingMovement defines motion that causes its base Wing's wing axes' yz-plane and
    its symmetry plane to transition from coincident to non-coincident, or vice
    versa. This is checked by this WingMovement's parent AirplaneMovement's parent
    Movement.

    This class contains the following public methods:

        generate_wings: Creates the Wing at each time step, and returns them in a list.

        max_period: Defines a property for the longest period of WingMovement's own
        motion and that of its sub-movement objects.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        base_wing,
        wing_cross_section_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    ):
        """This is the initialization method.

        :param base_wing: Wing

            This is the base Wing, from which the Wing at each time step will be
            created.

        :param wing_cross_section_movements: list of WingCrossSectionMovements

            This is a list of the WingCrossSectionMovements associated with each of
            the base Wing's WingCrossSections. It must have the same length as the
            base Wing's list of WingCrossSections.

        :param ampLer_Gs_Cgs: array-like of 3 numbers, optional

            The amplitudes of the WingMovement's changes in its Wings' Ler_Gs_Cgs
            parameters. Can be a tuple, list, or numpy array of non-negative numbers
            (int or float). Also, each amplitude must be low enough that it doesn't
            drive its base value out of the range of valid values. Otherwise,
            this WingMovement will try to create Wings with invalid parameters
            values. Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). The units are in meters.

        :param periodLer_Gs_Cgs: array-like of 3 numbers, optional

            The periods of the WingMovement's changes in its Wings' Ler_Gs_Cgs
            parameters. Can be a tuple, list, or numpy array of non-negative numbers
            (int or float). Values are converted to floats internally. The default
            value is (0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding
            element in ampLer_Gs_Cgs is 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingLer_Gs_Cgs: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingMovement's change in its
            Wings' Ler_Gs_Cgs parameters. Can be a tuple, list, or numpy array. Each
            element can be the string "sine", the string "uniform", or a callable
            custom spacing function. Custom spacing functions are for advanced users
            and must start at 0, return to 0 after one period of 2*pi radians,
            have zero mean, have amplitude of 1, be periodic, return finite values
            only, and accept a ndarray as input and return a ndarray of the same
            shape. The custom function is scaled by ampLer_Gs_Cgs, shifted by
            phaseLer_Gs_Cgs, and centered around the base value, with the period
            controlled by periodLer_Gs_Cgs. The default value is ("sine", "sine",
            "sine").

        :param phaseLer_Gs_Cgs: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Wing's
            Ler_Gs_Cgs parameter relative to the base Wing's Ler_Gs_Cgs parameter.
            Can be a tuple, list, or numpy array of non-negative numbers (int or
            float) in the range (-180.0, 180.0]. Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). Each element must be
            0.0 if the corresponding element in ampLer_Gs_Cgs is 0.0 and non-zero if
            not. The units are in degrees.

        :param ampAngles_Gs_to_Wn_ixyz: array-like of 3 numbers, optional

            The amplitudes of the WingMovement's changes in its Wings'
            angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list, or numpy array of
            numbers (int or float) in the range [0.0, 180.0]. Also, each amplitude
            must be low enough that it doesn't drive its base value out of the range
            of valid values. Otherwise, this WingMovement will try to create Wings
            with invalid parameters values. Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). The units are in degrees.

        :param periodAngles_Gs_to_Wn_ixyz: array-like of 3 numbers, optional

            The periods of the WingMovement's changes in its Wings'
            angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list, or numpy array of
            non-negative numbers (int or float). Values are converted to floats
            internally. The default value is (0.0, 0.0, 0.0). Each element must be
            0.0 if the corresponding element in ampAngles_Gs_to_Wn_ixyz is 0.0 and
            non-zero if not. The units are in seconds.

        :param spacingAngles_Gs_to_Wn_ixyz: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingMovement's change in its
            Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list, or numpy
            array. Each element can be the string "sine", the string "uniform",
            or a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0, return to 0 after one period of 2*pi
            radians, have zero mean, have amplitude of 1, be periodic, return finite
            values only, and accept a ndarray as input and return a ndarray of the
            same shape. The custom function is scaled by ampAngles_Gs_to_Wn_ixyz,
            shifted by phaseAngles_Gs_to_Wn_ixyz, and centered around the base value,
            with the period controlled by periodAngles_Gs_to_Wn_ixyz. The default
            value is ("sine", "sine", "sine").

        :param phaseAngles_Gs_to_Wn_ixyz: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's Wing's
            angles_Gs_to_Wn_ixyz parameter relative to the base Wing's
            angles_Gs_to_Wn_ixyz parameter. Can be a tuple, list, or numpy array of
            numbers (int or float) in the range (-180.0, 180.0]. Values are converted
            to floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampAngles_Gs_to_Wn_ixyz is
            0.0 and non-zero if not. The units are in degrees.
        """
        if not isinstance(base_wing, geometry.wing.Wing):
            raise TypeError("base_wing must be a Wing.")
        self.base_wing = base_wing

        if not isinstance(wing_cross_section_movements, list):
            raise TypeError("wing_cross_section_movements must be a list.")
        if len(wing_cross_section_movements) != len(self.base_wing.wing_cross_sections):
            raise ValueError(
                "wing_cross_section_movements must have the same length as base_wing.wing_cross_sections."
            )
        for wing_cross_section_movement in wing_cross_section_movements:
            if not isinstance(wing_cross_section_movement, WingCrossSectionMovement):
                raise TypeError(
                    "Every element in wing_cross_section_movements must be a WingCrossSectionMovement."
                )
        self.wing_cross_section_movements = wing_cross_section_movements

        ampLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            ampLer_Gs_Cgs, "ampLer_Gs_Cgs"
        )
        if not np.all(ampLer_Gs_Cgs >= 0.0):
            raise ValueError("All elements in ampLer_Gs_Cgs must be non-negative.")
        self.ampLer_Gs_Cgs = ampLer_Gs_Cgs

        periodLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            periodLer_Gs_Cgs, "periodLer_Gs_Cgs"
        )
        if not np.all(periodLer_Gs_Cgs >= 0.0):
            raise ValueError("All elements in periodLer_Gs_Cgs must be non-negative.")
        for period_index, period in enumerate(periodLer_Gs_Cgs):
            amp = self.ampLer_Gs_Cgs[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampLer_Gs_Cgs is 0.0, the corresponding element in periodLer_Gs_Cgs must be also be 0.0."
                )
        self.periodLer_Gs_Cgs = periodLer_Gs_Cgs

        spacingLer_Gs_Cgs = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingLer_Gs_Cgs,
                "spacingLer_Gs_Cgs",
            )
        )
        self.spacingLer_Gs_Cgs = spacingLer_Gs_Cgs

        phaseLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            phaseLer_Gs_Cgs, "phaseLer_Gs_Cgs"
        )
        if not (np.all(phaseLer_Gs_Cgs > -180.0) and np.all(phaseLer_Gs_Cgs <= 180.0)):
            raise ValueError(
                "All elements in phaseLer_Gs_Cgs must be in the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseLer_Gs_Cgs):
            amp = self.ampLer_Gs_Cgs[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampLer_Gs_Cgs is 0.0, the corresponding element in phaseLer_Gs_Cgs must be also be 0.0."
                )
        self.phaseLer_Gs_Cgs = phaseLer_Gs_Cgs

        ampAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                ampAngles_Gs_to_Wn_ixyz, "ampAngles_Gs_to_Wn_ixyz"
            )
        )
        if not (
            np.all(ampAngles_Gs_to_Wn_ixyz >= 0.0)
            and np.all(ampAngles_Gs_to_Wn_ixyz <= 180.0)
        ):
            raise ValueError(
                "All elements in ampAngles_Gs_to_Wn_ixyz must be in the range [0.0, 180.0]."
            )
        self.ampAngles_Gs_to_Wn_ixyz = ampAngles_Gs_to_Wn_ixyz

        periodAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_Gs_to_Wn_ixyz, "periodAngles_Gs_to_Wn_ixyz"
            )
        )
        if not np.all(periodAngles_Gs_to_Wn_ixyz >= 0.0):
            raise ValueError(
                "All elements in periodAngles_Gs_to_Wn_ixyz must be non-negative."
            )
        for period_index, period in enumerate(periodAngles_Gs_to_Wn_ixyz):
            amp = self.ampAngles_Gs_to_Wn_ixyz[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_Gs_to_Wn_ixyz is 0.0, the corresponding element in periodAngles_Gs_to_Wn_ixyz must be also be 0.0."
                )
        self.periodAngles_Gs_to_Wn_ixyz = periodAngles_Gs_to_Wn_ixyz

        spacingAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingAngles_Gs_to_Wn_ixyz,
                "spacingAngles_Gs_to_Wn_ixyz",
            )
        )
        self.spacingAngles_Gs_to_Wn_ixyz = spacingAngles_Gs_to_Wn_ixyz

        phaseAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                phaseAngles_Gs_to_Wn_ixyz, "phaseAngles_Gs_to_Wn_ixyz"
            )
        )
        if not (
            np.all(phaseAngles_Gs_to_Wn_ixyz > -180.0)
            and np.all(phaseAngles_Gs_to_Wn_ixyz <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_Gs_to_Wn_ixyz must be in the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_Gs_to_Wn_ixyz):
            amp = self.ampAngles_Gs_to_Wn_ixyz[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_Gs_to_Wn_ixyz is 0.0, the corresponding element in phaseAngles_Gs_to_Wn_ixyz must be also be 0.0."
                )
        self.phaseAngles_Gs_to_Wn_ixyz = phaseAngles_Gs_to_Wn_ixyz

    def generate_wings(self, num_steps, delta_time):
        """Creates the Wing at each time step, and returns them in a list.

        :param num_steps: int

            This is the number of time steps in this movement. It must be a positive
            int.

        :param delta_time: number

            This is the time between each time step. It must be a positive number (
            int or float), and will be converted internally to a float. The units are
            in seconds.

        :return: list of Wings

            This is the list of Wings associated with this WingMovement.
        """
        num_steps = _parameter_validation.positive_int_return_int(
            num_steps, "num_steps"
        )
        delta_time = _parameter_validation.positive_number_return_float(
            delta_time, "delta_time"
        )

        # Generate oscillating values for each dimension of Ler_Gs_Cgs.
        listLer_Gs_Cgs = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingLer_Gs_Cgs[dim]
            if spacing == "sine":
                listLer_Gs_Cgs[dim, :] = _functions.oscillating_sinspaces(
                    amps=self.ampLer_Gs_Cgs[dim],
                    periods=self.periodLer_Gs_Cgs[dim],
                    phases=self.phaseLer_Gs_Cgs[dim],
                    bases=self.base_wing.Ler_Gs_Cgs[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                listLer_Gs_Cgs[dim, :] = _functions.oscillating_linspaces(
                    amps=self.ampLer_Gs_Cgs[dim],
                    periods=self.periodLer_Gs_Cgs[dim],
                    phases=self.phaseLer_Gs_Cgs[dim],
                    bases=self.base_wing.Ler_Gs_Cgs[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                listLer_Gs_Cgs[dim, :] = _functions.oscillating_customspaces(
                    amps=self.ampLer_Gs_Cgs[dim],
                    periods=self.periodLer_Gs_Cgs[dim],
                    phases=self.phaseLer_Gs_Cgs[dim],
                    bases=self.base_wing.Ler_Gs_Cgs[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

        # Generate oscillating values for each dimension of angles_Gs_to_Wn_ixyz.
        listAngles_Gs_to_Wn_ixyz = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingAngles_Gs_to_Wn_ixyz[dim]
            if spacing == "sine":
                listAngles_Gs_to_Wn_ixyz[dim, :] = _functions.oscillating_sinspaces(
                    amps=self.ampAngles_Gs_to_Wn_ixyz[dim],
                    periods=self.periodAngles_Gs_to_Wn_ixyz[dim],
                    phases=self.phaseAngles_Gs_to_Wn_ixyz[dim],
                    bases=self.base_wing.angles_Gs_to_Wn_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                listAngles_Gs_to_Wn_ixyz[dim, :] = _functions.oscillating_linspaces(
                    amps=self.ampAngles_Gs_to_Wn_ixyz[dim],
                    periods=self.periodAngles_Gs_to_Wn_ixyz[dim],
                    phases=self.phaseAngles_Gs_to_Wn_ixyz[dim],
                    bases=self.base_wing.angles_Gs_to_Wn_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                listAngles_Gs_to_Wn_ixyz[dim, :] = _functions.oscillating_customspaces(
                    amps=self.ampAngles_Gs_to_Wn_ixyz[dim],
                    periods=self.periodAngles_Gs_to_Wn_ixyz[dim],
                    phases=self.phaseAngles_Gs_to_Wn_ixyz[dim],
                    bases=self.base_wing.angles_Gs_to_Wn_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

        # Create an empty 2D ndarray that will hold each of the Wings's
        # WingCrossSection's vector of WingCrossSections representing its changing
        # state at each time step. The first index denotes a particular base
        # WingCrossSection, and the second index denotes the time step.
        wing_cross_sections = np.empty(
            (len(self.wing_cross_section_movements), num_steps), dtype=object
        )

        # Iterate through the WingCrossSectionMovements.
        for (
            wing_cross_section_movement_id,
            wing_cross_section_movement,
        ) in enumerate(self.wing_cross_section_movements):

            # Generate this WingCrossSection's vector of WingCrossSections
            # representing its changing state at each time step.
            this_wing_cross_sections_list_of_wing_cross_sections = np.array(
                wing_cross_section_movement.generate_wing_cross_sections(
                    num_steps=num_steps, delta_time=delta_time
                )
            )

            # Add this vector the Wing's 2D ndarray of WingCrossSections'
            # WingCrossSections.
            wing_cross_sections[wing_cross_section_movement_id, :] = (
                this_wing_cross_sections_list_of_wing_cross_sections
            )

        # Create an empty list to hold each time step's Wing.
        wings = []

        # Get the non-changing Wing attributes.
        this_name = self.base_wing.name
        this_symmetric = self.base_wing.symmetric
        this_mirror_only = self.base_wing.mirror_only
        this_symmetryNormal_G = self.base_wing.symmetryNormal_G
        this_symmetryPoint_G_Cg = self.base_wing.symmetryPoint_G_Cg
        this_num_chordwise_panels = self.base_wing.num_chordwise_panels
        this_chordwise_spacing = self.base_wing.chordwise_spacing

        # Iterate through the time steps.
        for step in range(num_steps):
            thisLer_Gs_Cgs = listLer_Gs_Cgs[:, step]
            theseAngles_Gs_to_Wn_ixyz = listAngles_Gs_to_Wn_ixyz[:, step]
            these_wing_cross_sections = list(wing_cross_sections[:, step])

            # Make a new Wing for this time step.
            this_wing = geometry.wing.Wing(
                wing_cross_sections=these_wing_cross_sections,
                name=this_name,
                Ler_Gs_Cgs=thisLer_Gs_Cgs,
                angles_Gs_to_Wn_ixyz=theseAngles_Gs_to_Wn_ixyz,
                symmetric=this_symmetric,
                mirror_only=this_mirror_only,
                symmetryNormal_G=this_symmetryNormal_G,
                symmetryPoint_G_Cg=this_symmetryPoint_G_Cg,
                num_chordwise_panels=this_num_chordwise_panels,
                chordwise_spacing=this_chordwise_spacing,
            )

            # Add this new Wing to the list of Wings.
            wings.append(this_wing)

        return wings

    @property
    def max_period(self):
        """Defines a property for the longest period of WingMovement's own motion and
        that of its sub-movement objects.

        :return: float

            The longest period in seconds. If the all the motion is static, this will
            be 0.0.
        """
        wing_cross_section_movement_max_periods = []
        for wing_cross_section_movement in self.wing_cross_section_movements:
            wing_cross_section_movement_max_periods.append(
                wing_cross_section_movement.max_period
            )
        max_wing_cross_section_movement_period = max(
            wing_cross_section_movement_max_periods
        )

        return float(
            max(
                max_wing_cross_section_movement_period,
                np.max(self.periodLer_Gs_Cgs),
                np.max(self.periodAngles_Gs_to_Wn_ixyz),
            )
        )
