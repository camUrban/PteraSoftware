"""This module contains the WingCrossSectionMovement class.

This module contains the following classes:
    WingCrossSectionMovement: This is a class used to contain the WingCrossSection
    movements.

This module contains the following functions:
    None
"""

import numpy as np

from . import _functions

from .. import geometry
from .. import _parameter_validation


class WingCrossSectionMovement:
    """This is a class used to contain the WingCrossSection movements.

    This class contains the following public methods:

        generate_wing_cross_sections: Creates the WingCrossSection at each time step,
        and returns them in a list.

        max_period: Defines a property for the longest period of
        WingCrossSectionMovement's own motion.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    ):
        """This is the initialization method.

        :param base_wing_cross_section: WingCrossSection

            This is the base WingCrossSection, from which the WingCrossSection at
            each time step will be created.

        :param ampLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The amplitudes of the WingCrossSectionMovement's changes in its
            WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Also, each amplitude must
            be low enough that it doesn't drive its base value out of the range of
            valid values. Otherwise, this WingCrossSectionMovement will try to create
            WingCrossSections with invalid parameters values. Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). The units are in
            meters.

        :param periodLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The periods of the WingCrossSectionMovement's changes in its
            WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampLp_Wcsp_Lpp is 0.0 and
            non-zero if not. The units are in seconds.

        :param spacingLp_Wcsp_Lpp: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingCrossSectionMovement's change
            in its WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list,
            or numpy array. Each element can be the string "sine", the string
            "uniform", or a callable custom spacing function. Custom spacing functions
            are for advanced users and must start at 0, return to 0 after one period
            of 2*pi radians, have zero mean, have amplitude of 1, be periodic, return
            finite values only, and accept a ndarray as input and return a ndarray of
            the same shape. The custom function is scaled by ampLp_Wcsp_Lpp, shifted
            by phaseLp_Wcsp_Lpp, and centered around the base value, with the period
            controlled by periodLp_Wcsp_Lpp. The default value is ("sine", "sine",
            "sine").

        :param phaseLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's
            WingCrossSection's Lp_Wcsp_Lpp parameter relative to the base
            WingCrossSection's Lp_Wcsp_Lpp parameter. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float) in the range (-180.0,
            180.0]. Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding element in
            ampLp_Wcsp_Lpp is 0.0 and non-zero if not. The units are in degrees.

        :param ampAngles_Wcsp_to_Wcs_ixyz: array-like of 3 numbers, optional

            The amplitudes of the WingCrossSectionMovement's changes in its
            WingCrossSections' angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple,
            list, or numpy array of numbers (int or float) in the range [0.0,
            180.0]. Also, each amplitude must be low enough that it doesn't drive its
            base value out of the range of valid values. Otherwise,
            this WingCrossSectionMovement will try to create WingCrossSections with
            invalid parameters values. Values are converted to floats internally. The
            default value is (0.0, 0.0, 0.0). The units are in degrees.

        :param periodAngles_Wcsp_to_Wcs_ixyz: array-like of 3 numbers, optional

            The periods of the WingCrossSectionMovement's changes in its
            WingCrossSections' angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple,
            list, or numpy array of non-negative numbers (int or float). Values are
            converted to floats internally. The default value is (0.0, 0.0,
            0.0). Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingAngles_Wcsp_to_Wcs_ixyz: array-like of 3 strs or callables, optional

            The value determines the spacing of the WingCrossSectionMovement's change
            in its WingCrossSections' angles_Wcsp_to_Wcs_ixyz parameters. Can be a
            tuple, list, or numpy array. Each element can be the string "sine", the
            string "uniform", or a callable custom spacing function. Custom spacing
            functions are for advanced users and must start at 0, return to 0 after
            one period of 2*pi radians, have zero mean, have amplitude of 1, be
            periodic, return finite values only, and accept a ndarray as input and
            return a ndarray of the same shape. The custom function is scaled by
            ampAngles_Wcsp_to_Wcs_ixyz, shifted by phaseAngles_Wcsp_to_Wcs_ixyz, and
            centered around the base value, with the period controlled by
            periodAngles_Wcsp_to_Wcs_ixyz. The default value is ("sine", "sine",
            "sine").

        :param phaseAngles_Wcsp_to_Wcs_ixyz: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's
            WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter relative to the base
            WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter. Can be a tuple,
            list, or numpy array of numbers (int or float) in the range (-180.0,
            180.0]. Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non-zero if not. The units are in
            degrees.
        """
        if not isinstance(
            base_wing_cross_section, geometry.wing_cross_section.WingCrossSection
        ):
            raise TypeError("base_wing_cross_section must be a WingCrossSection.")
        self.base_wing_cross_section = base_wing_cross_section

        ampLp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            ampLp_Wcsp_Lpp, "ampLp_Wcsp_Lpp"
        )
        if not np.all(ampLp_Wcsp_Lpp >= 0.0):
            raise ValueError("All elements in ampLp_Wcsp_Lpp must be non-negative.")
        self.ampLp_Wcsp_Lpp = ampLp_Wcsp_Lpp

        periodLp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            periodLp_Wcsp_Lpp, "periodLp_Wcsp_Lpp"
        )
        if not np.all(periodLp_Wcsp_Lpp >= 0.0):
            raise ValueError("All elements in periodLp_Wcsp_Lpp must be non-negative.")
        for period_index, period in enumerate(periodLp_Wcsp_Lpp):
            amp = self.ampLp_Wcsp_Lpp[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampLp_Wcsp_Lpp is 0.0, the corresponding element in periodLp_Wcsp_Lpp must be also be 0.0."
                )
        self.periodLp_Wcsp_Lpp = periodLp_Wcsp_Lpp

        spacingLp_Wcsp_Lpp = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingLp_Wcsp_Lpp, "spacingLp_Wcsp_Lpp"
            )
        )
        self.spacingLp_Wcsp_Lpp = spacingLp_Wcsp_Lpp

        phaseLp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            phaseLp_Wcsp_Lpp, "phaseLp_Wcsp_Lpp"
        )
        if not (
            np.all(phaseLp_Wcsp_Lpp > -180.0) and np.all(phaseLp_Wcsp_Lpp <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseLp_Wcsp_Lpp must be in the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseLp_Wcsp_Lpp):
            amp = self.ampLp_Wcsp_Lpp[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampLp_Wcsp_Lpp is 0.0, the corresponding element in phaseLp_Wcsp_Lpp must be also be 0.0."
                )
        self.phaseLp_Wcsp_Lpp = phaseLp_Wcsp_Lpp

        ampAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                ampAngles_Wcsp_to_Wcs_ixyz, "ampAngles_Wcsp_to_Wcs_ixyz"
            )
        )
        if not (
            np.all(ampAngles_Wcsp_to_Wcs_ixyz >= 0.0)
            and np.all(ampAngles_Wcsp_to_Wcs_ixyz <= 180.0)
        ):
            raise ValueError(
                "All elements in ampAngles_Wcsp_to_Wcs_ixyz must be in the range [0.0, 180.0]."
            )
        self.ampAngles_Wcsp_to_Wcs_ixyz = ampAngles_Wcsp_to_Wcs_ixyz

        periodAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_Wcsp_to_Wcs_ixyz, "periodAngles_Wcsp_to_Wcs_ixyz"
            )
        )
        if not np.all(periodAngles_Wcsp_to_Wcs_ixyz >= 0.0):
            raise ValueError(
                "All elements in periodAngles_Wcsp_to_Wcs_ixyz must be non-negative."
            )
        for period_index, period in enumerate(periodAngles_Wcsp_to_Wcs_ixyz):
            amp = self.ampAngles_Wcsp_to_Wcs_ixyz[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0, the corresponding element in periodAngles_Wcsp_to_Wcs_ixyz must be also be 0.0."
                )
        self.periodAngles_Wcsp_to_Wcs_ixyz = periodAngles_Wcsp_to_Wcs_ixyz

        spacingAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingAngles_Wcsp_to_Wcs_ixyz,
                "spacingAngles_Wcsp_to_Wcs_ixyz",
            )
        )
        self.spacingAngles_Wcsp_to_Wcs_ixyz = spacingAngles_Wcsp_to_Wcs_ixyz

        phaseAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                phaseAngles_Wcsp_to_Wcs_ixyz, "phaseAngles_Wcsp_to_Wcs_ixyz"
            )
        )
        if not (
            np.all(phaseAngles_Wcsp_to_Wcs_ixyz > -180.0)
            and np.all(phaseAngles_Wcsp_to_Wcs_ixyz <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_Wcsp_to_Wcs_ixyz must be in the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_Wcsp_to_Wcs_ixyz):
            amp = self.ampAngles_Wcsp_to_Wcs_ixyz[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0, the corresponding element in phaseAngles_Wcsp_to_Wcs_ixyz must be also be 0.0."
                )
        self.phaseAngles_Wcsp_to_Wcs_ixyz = phaseAngles_Wcsp_to_Wcs_ixyz

    def generate_wing_cross_sections(
        self,
        num_steps,
        delta_time,
    ):
        """Creates the WingCrossSection at each time step, and returns them in a list.

        :param num_steps: int

            This is the number of time steps in this movement. It must be a positive
            int.

        :param delta_time: number

            This is the time between each time step. It must be a positive number (
            int or float), and will be converted internally to a float. The units are
            in seconds.

        :return: list of WingCrossSections

            This is the list of WingCrossSections associated with this
            WingCrossSectionMovement.
        """
        num_steps = _parameter_validation.positive_int_return_int(
            num_steps, "num_steps"
        )
        delta_time = _parameter_validation.positive_number_return_float(
            delta_time, "delta_time"
        )

        # Generate oscillating values for each dimension of Lp_Wcsp_Lpp.
        listLp_Wcsp_Lpp = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingLp_Wcsp_Lpp[dim]
            if spacing == "sine":
                listLp_Wcsp_Lpp[dim, :] = _functions.oscillating_sinspaces(
                    amps=self.ampLp_Wcsp_Lpp[dim],
                    periods=self.periodLp_Wcsp_Lpp[dim],
                    phases=self.phaseLp_Wcsp_Lpp[dim],
                    bases=self.base_wing_cross_section.Lp_Wcsp_Lpp[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                listLp_Wcsp_Lpp[dim, :] = _functions.oscillating_linspaces(
                    amps=self.ampLp_Wcsp_Lpp[dim],
                    periods=self.periodLp_Wcsp_Lpp[dim],
                    phases=self.phaseLp_Wcsp_Lpp[dim],
                    bases=self.base_wing_cross_section.Lp_Wcsp_Lpp[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                listLp_Wcsp_Lpp[dim, :] = _functions.oscillating_customspaces(
                    amps=self.ampLp_Wcsp_Lpp[dim],
                    periods=self.periodLp_Wcsp_Lpp[dim],
                    phases=self.phaseLp_Wcsp_Lpp[dim],
                    bases=self.base_wing_cross_section.Lp_Wcsp_Lpp[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

        # Generate oscillating values for each dimension of angles_Wcsp_to_Wcs_ixyz.
        listAngles_Wcsp_to_Wcs_ixyz = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingAngles_Wcsp_to_Wcs_ixyz[dim]
            if spacing == "sine":
                listAngles_Wcsp_to_Wcs_ixyz[dim, :] = _functions.oscillating_sinspaces(
                    amps=self.ampAngles_Wcsp_to_Wcs_ixyz[dim],
                    periods=self.periodAngles_Wcsp_to_Wcs_ixyz[dim],
                    phases=self.phaseAngles_Wcsp_to_Wcs_ixyz[dim],
                    bases=self.base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                listAngles_Wcsp_to_Wcs_ixyz[dim, :] = _functions.oscillating_linspaces(
                    amps=self.ampAngles_Wcsp_to_Wcs_ixyz[dim],
                    periods=self.periodAngles_Wcsp_to_Wcs_ixyz[dim],
                    phases=self.phaseAngles_Wcsp_to_Wcs_ixyz[dim],
                    bases=self.base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                listAngles_Wcsp_to_Wcs_ixyz[dim, :] = (
                    _functions.oscillating_customspaces(
                        amps=self.ampAngles_Wcsp_to_Wcs_ixyz[dim],
                        periods=self.periodAngles_Wcsp_to_Wcs_ixyz[dim],
                        phases=self.phaseAngles_Wcsp_to_Wcs_ixyz[dim],
                        bases=self.base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim],
                        num_steps=num_steps,
                        delta_time=delta_time,
                        custom_function=spacing,
                    )
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

        # Create an empty list to hold each time step's WingCrossSection.
        wing_cross_sections = []

        # Get the non-changing WingCrossSectionAttributes.
        this_airfoil = self.base_wing_cross_section.airfoil
        this_num_spanwise_panels = self.base_wing_cross_section.num_spanwise_panels
        this_chord = self.base_wing_cross_section.chord
        this_control_surface_symmetry_type = (
            self.base_wing_cross_section.control_surface_symmetry_type
        )
        this_control_surface_hinge_point = (
            self.base_wing_cross_section.control_surface_hinge_point
        )
        this_control_surface_deflection = (
            self.base_wing_cross_section.control_surface_deflection
        )
        this_spanwise_spacing = self.base_wing_cross_section.spanwise_spacing

        # Iterate through the time steps.
        for step in range(num_steps):
            thisLp_Wcsp_Lpp = listLp_Wcsp_Lpp[:, step]
            theseAngles_Wcsp_to_Wcs_ixyz = listAngles_Wcsp_to_Wcs_ixyz[:, step]

            # Make a new WingCrossSection for this time step.
            this_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                airfoil=this_airfoil,
                num_spanwise_panels=this_num_spanwise_panels,
                chord=this_chord,
                Lp_Wcsp_Lpp=thisLp_Wcsp_Lpp,
                angles_Wcsp_to_Wcs_ixyz=theseAngles_Wcsp_to_Wcs_ixyz,
                control_surface_symmetry_type=this_control_surface_symmetry_type,
                control_surface_hinge_point=this_control_surface_hinge_point,
                control_surface_deflection=this_control_surface_deflection,
                spanwise_spacing=this_spanwise_spacing,
            )

            # Add this new WingCrossSection to the list of WingCrossSections.
            wing_cross_sections.append(this_wing_cross_section)

        return wing_cross_sections

    @property
    def max_period(self):
        """Defines a property for the longest period of WingCrossSectionMovement's
        own motion.

        :return: float

            The longest period in seconds. If the all the motion is static, this will
            be 0.0.
        """
        return float(
            max(
                np.max(self.periodLp_Wcsp_Lpp),
                np.max(self.periodAngles_Wcsp_to_Wcs_ixyz),
            )
        )
