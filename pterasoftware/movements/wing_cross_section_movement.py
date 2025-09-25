"""This module contains the WingCrossSectionMovement class.

This module contains the following classes:
    WingCrossSectionMovement: This is a class used to contain the WingCrossSection
    movements.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from . import functions

from .. import geometry
from .. import parameter_validation


# TODO: Add unit tests for this class.
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
        ampAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
    ):
        """This is the initialization method.

        :param base_wing_cross_section: WingCrossSection

            This is the base WingCrossSection, from which the WingCrossSection at
            each time step will be created.

        :param ampLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The amplitudes of the WingCrossSectionMovement's changes in its
            WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). The units are in
            meters.

        :param periodLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The periods of the WingCrossSectionMovement's changes in its
            WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float). Values are converted to
            floats internally. The default value is (0.0, 0.0, 0.0). Each element
            must be 0.0 if the corresponding element in ampLp_Wcsp_Lpp is 0.0 and
            non-zero if not. The units are in seconds.

        :param spacingLp_Wcsp_Lpp: array-like of 3 strs, optional

            The value determines the spacing of the WingCrossSectionMovement's change
            in its WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple, list,
            or numpy array of strings. Each element can be "sine" or "uniform". The
            default value is ("sine", "sine", "sine").

        :param phaseLp_Wcsp_Lpp: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's
            WingCrossSection's Lp_Wcsp_Lpp parameter relative to the base
            WingCrossSection's Lp_Wcsp_Lpp parameter. Can be a tuple, list, or numpy
            array of non-negative numbers (int or float) in the range [0.0, 360.0).
            Values are converted to floats internally. The default value is (0.0,
            0.0, 0.0). Each element must be 0.0 if the corresponding element in
            ampLp_Wcsp_Lpp is 0.0 and non-zero if not. The units are in degrees.

        :param ampAngles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            The amplitudes of the WingCrossSectionMovement's changes in its
            WingCrossSections' angles_Wcsp_to_Wcs_izyx parameters. Can be a tuple,
            list, or numpy array of numbers (int or float) in the range [0.0,
            180.0). Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). The units are in degrees.

        :param periodAngles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            The periods of the WingCrossSectionMovement's changes in its
            WingCrossSections' angles_Wcsp_to_Wcs_izyx parameters. Can be a tuple,
            list, or numpy array of non-negative numbers (int or float). Values are
            converted to floats internally. The default value is (0.0, 0.0,
            0.0). Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_izyx is 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingAngles_Wcsp_to_Wcs_izyx: array-like of 3 strs, optional

            The value determines the spacing of the WingCrossSectionMovement's change
            in its WingCrossSections' angles_Wcsp_to_Wcs_izyx parameters. Can be a
            tuple, list, or numpy array of strings. Each element can be "sine" or
            "uniform". The default value is ("sine", "sine", "sine").

        :param phaseAngles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            The phase offsets of the elements in the first time step's
            WingCrossSection's angles_Wcsp_to_Wcs_izyx parameter relative to the base
            WingCrossSection's angles_Wcsp_to_Wcs_izyx parameter. Can be a tuple,
            list, or numpy array of numbers (int or float) in the range [0.0,
            360.0). Values are converted to floats internally. The default value is (
            0.0, 0.0, 0.0). Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_izyx is 0.0 and non-zero if not. The units are in
            degrees.
        """
        if not isinstance(
            base_wing_cross_section, geometry.wing_cross_section.WingCrossSection
        ):
            raise TypeError("base_wing_cross_section must be a WingCrossSection.")
        self.base_wing_cross_section = base_wing_cross_section

        ampLp_Wcsp_Lpp = parameter_validation.threeD_number_vectorLike_return_float(
            ampLp_Wcsp_Lpp, "ampLp_Wcsp_Lpp"
        )
        if not np.all(ampLp_Wcsp_Lpp >= 0.0):
            raise ValueError("All elements in ampLp_Wcsp_Lpp must be non-negative.")
        self.ampLp_Wcsp_Lpp = ampLp_Wcsp_Lpp

        periodLp_Wcsp_Lpp = parameter_validation.threeD_number_vectorLike_return_float(
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

        spacingLp_Wcsp_Lpp = parameter_validation.list_return_list(
            spacingLp_Wcsp_Lpp, "spacingLp_Wcsp_Lpp"
        )
        if not np.all(elem in ["sine", "uniform"] for elem in spacingLp_Wcsp_Lpp):
            raise ValueError(
                'All elements in spacingLp_Wcsp_Lpp must be "sine" or "uniform".'
            )
        self.spacingLp_Wcsp_Lpp = spacingLp_Wcsp_Lpp

        phaseLp_Wcsp_Lpp = parameter_validation.threeD_number_vectorLike_return_float(
            phaseLp_Wcsp_Lpp, "phaseLp_Wcsp_Lpp"
        )
        if not (np.all(phaseLp_Wcsp_Lpp >= 0.0) and np.all(phaseLp_Wcsp_Lpp < 360.0)):
            raise ValueError(
                "All elements in phaseLp_Wcsp_Lpp must be in the range [0.0, 360.0)."
            )
        for phase_index, phase in enumerate(phaseLp_Wcsp_Lpp):
            amp = self.ampLp_Wcsp_Lpp[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampLp_Wcsp_Lpp is 0.0, the corresponding element in phaseLp_Wcsp_Lpp must be also be 0.0."
                )
        self.phaseLp_Wcsp_Lpp = phaseLp_Wcsp_Lpp

        ampAngles_Wcsp_to_Wcs_izyx = (
            parameter_validation.threeD_number_vectorLike_return_float(
                ampAngles_Wcsp_to_Wcs_izyx, "ampAngles_Wcsp_to_Wcs_izyx"
            )
        )
        if not (
            np.all(ampAngles_Wcsp_to_Wcs_izyx >= 0.0)
            and np.all(ampAngles_Wcsp_to_Wcs_izyx < 180.0)
        ):
            raise ValueError(
                "All elements in ampAngles_Wcsp_to_Wcs_izyx must be in the range [0.0, 180.0)."
            )
        self.ampAngles_Wcsp_to_Wcs_izyx = ampAngles_Wcsp_to_Wcs_izyx

        periodAngles_Wcsp_to_Wcs_izyx = (
            parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_Wcsp_to_Wcs_izyx, "periodAngles_Wcsp_to_Wcs_izyx"
            )
        )
        if not np.all(periodAngles_Wcsp_to_Wcs_izyx >= 0.0):
            raise ValueError(
                "All elements in periodAngles_Wcsp_to_Wcs_izyx must be non-negative."
            )
        for period_index, period in enumerate(periodAngles_Wcsp_to_Wcs_izyx):
            amp = self.ampAngles_Wcsp_to_Wcs_izyx[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_Wcsp_to_Wcs_izyx is 0.0, the corresponding element in periodAngles_Wcsp_to_Wcs_izyx must be also be 0.0."
                )
        self.periodAngles_Wcsp_to_Wcs_izyx = periodAngles_Wcsp_to_Wcs_izyx

        spacingAngles_Wcsp_to_Wcs_izyx = parameter_validation.list_return_list(
            spacingAngles_Wcsp_to_Wcs_izyx, "spacingAngles_Wcsp_to_Wcs_izyx"
        )
        if not np.all(
            elem in ["sine", "uniform"] for elem in spacingAngles_Wcsp_to_Wcs_izyx
        ):
            raise ValueError(
                'All elements in spacingAngles_Wcsp_to_Wcs_izyx must be "sine" or "uniform".'
            )
        self.spacingAngles_Wcsp_to_Wcs_izyx = spacingAngles_Wcsp_to_Wcs_izyx

        phaseAngles_Wcsp_to_Wcs_izyx = (
            parameter_validation.threeD_number_vectorLike_return_float(
                phaseAngles_Wcsp_to_Wcs_izyx, "phaseAngles_Wcsp_to_Wcs_izyx"
            )
        )
        if not (
            np.all(phaseAngles_Wcsp_to_Wcs_izyx >= 0.0)
            and np.all(phaseAngles_Wcsp_to_Wcs_izyx < 360.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_Wcsp_to_Wcs_izyx must be in the range [0.0, 360.0)."
            )
        for phase_index, phase in enumerate(phaseAngles_Wcsp_to_Wcs_izyx):
            amp = self.ampAngles_Wcsp_to_Wcs_izyx[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_Wcsp_to_Wcs_izyx is 0.0, the corresponding element in phaseAngles_Wcsp_to_Wcs_izyx must be also be 0.0."
                )
        self.phaseAngles_Wcsp_to_Wcs_izyx = phaseAngles_Wcsp_to_Wcs_izyx

    # TODO: Add unit tests for this method.
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
        num_steps = parameter_validation.positive_int_return_int(num_steps, "num_steps")
        delta_time = parameter_validation.positive_number_return_float(
            delta_time, "delta_time"
        )

        # FIXME: I'm pretty sure the oscillating_* functions are going to break if I
        #  use them like this. They expect single values for the inputs but I'm
        #  passing vectors. I need to modify them to accept and return vectors.
        if self.spacingLp_Wcsp_Lpp == "sine":
            # FIXME: Add a phase offset parameter to oscillating_sinspace
            listLp_Wcsp_Lpp = functions.oscillating_sinspace(
                amplitude=self.ampLp_Wcsp_Lpp,
                period=self.periodLp_Wcsp_Lpp,
                base_value=self.base_wing_cross_section.Lp_Wcsp_Lpp,
                phase_offset=self.phaseLp_Wcsp_Lpp,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:
            # FIXME: Add a phase offset parameter to oscillating_linspace
            listLp_Wcsp_Lpp = functions.oscillating_linspace(
                amplitude=self.ampLp_Wcsp_Lpp,
                period=self.periodLp_Wcsp_Lpp,
                base_value=self.base_wing_cross_section.Lp_Wcsp_Lpp,
                phase_offset=self.phaseLp_Wcsp_Lpp,
                num_steps=num_steps,
                delta_time=delta_time,
            )

        if self.spacingAngles_Wcsp_to_Wcs_izyx == "sine":
            listAngles_Wcsp_to_Wcs_izyx = functions.oscillating_sinspace(
                amplitude=self.ampAngles_Wcsp_to_Wcs_izyx,
                period=self.periodAngles_Wcsp_to_Wcs_izyx,
                base_value=self.base_wing_cross_section.angles_Wcsp_to_Wcs_izyx,
                phase_offset=self.phaseAngles_Wcsp_to_Wcs_izyx,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:
            listAngles_Wcsp_to_Wcs_izyx = functions.oscillating_linspace(
                amplitude=self.ampAngles_Wcsp_to_Wcs_izyx,
                period=self.periodAngles_Wcsp_to_Wcs_izyx,
                base_value=self.base_wing_cross_section.angles_Wcsp_to_Wcs_izyx,
                phase_offset=self.phaseAngles_Wcsp_to_Wcs_izyx,
                num_steps=num_steps,
                delta_time=delta_time,
            )

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
            thisLp_Wcsp_Lpp = listLp_Wcsp_Lpp[step]
            theseAngles_Wcsp_to_Wcs_izyx = listAngles_Wcsp_to_Wcs_izyx[step]

            # Make a new WingCrossSection for this time step.
            this_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                airfoil=this_airfoil,
                num_spanwise_panels=this_num_spanwise_panels,
                chord=this_chord,
                Lp_Wcsp_Lpp=thisLp_Wcsp_Lpp,
                angles_Wcsp_to_Wcs_izyx=theseAngles_Wcsp_to_Wcs_izyx,
                control_surface_symmetry_type=this_control_surface_symmetry_type,
                control_surface_hinge_point=this_control_surface_hinge_point,
                control_surface_deflection=this_control_surface_deflection,
                spanwise_spacing=this_spanwise_spacing,
            )

            # Add this new WingCrossSection to the list of WingCrossSections.
            wing_cross_sections.append(this_wing_cross_section)

        return wing_cross_sections

    # TODO: Add unit tests for this method.
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
                np.max(self.periodAngles_Wcsp_to_Wcs_izyx),
            )
        )
