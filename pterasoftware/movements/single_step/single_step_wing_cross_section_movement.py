"""Contains the SingleStepWingCrossSectionMovement class.

**Contains the following classes:**

SingleStepWingCrossSectionMovement: A single step variant of
WingCrossSectionMovement that generates one WingCrossSection per time step instead of
all at once. Uses composition to wrap a WingCrossSectionMovement.

**Contains the following functions:**

None
"""

import numpy as np

from ..wing_cross_section_movement import WingCrossSectionMovement

from ..._parameter_validation import (
    int_in_range_return_int,
    number_in_range_return_float
)

from .._functions import (
    oscillating_sinspaces,
    oscillating_linspaces,
    oscillating_customspaces,
)

from ... import geometry


class SingleStepWingCrossSectionMovement:
    """A single step variant of WingCrossSectionMovement for coupled simulations.

    This class wraps a WingCrossSectionMovement via composition and generates one
    WingCrossSection per time step (via generate_next_wing_cross_sections) rather than
    generating all WingCrossSections at once. The composed WingCrossSectionMovement is
    accessible via corresponding_wing_cross_section_movement.

    **Contains the following methods:**

    generate_next_wing_cross_sections: Creates the WingCrossSection at a single time
    step.

    max_period: The longest period of this SingleStepWingCrossSectionMovement's own
    motion.
    """

    __slots__ = (
        "ampLp_Wcsp_Lpp",
        "periodLp_Wcsp_Lpp",
        "spacingLp_Wcsp_Lpp",
        "phaseLp_Wcsp_Lpp",
        "ampAngles_Wcsp_to_Wcs_ixyz",
        "periodAngles_Wcsp_to_Wcs_ixyz",
        "spacingAngles_Wcsp_to_Wcs_ixyz",
        "phaseAngles_Wcsp_to_Wcs_ixyz",
        "listLp_Wcsp_Lpp",
        "listAngles_Wcsp_to_Wcs_ixyz",
        "corresponding_wing_cross_section_movement",
    )

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
        """The initialization method.

        :param base_wing_cross_section: The base WingCrossSection from which the
            WingCrossSection at each time step will be created.
        :param ampLp_Wcsp_Lpp: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            SingleStepWingCrossSectionMovement's changes in its WingCrossSections'
            Lp_Wcsp_Lpp parameters. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Each amplitude must be low enough that it
            doesn't drive its base value out of the range of valid values. Otherwise,
            this SingleStepWingCrossSectionMovement will try to create WingCrossSections
            with invalid parameter values. The units are in meters. The default is
            (0.0, 0.0, 0.0).
        :param periodLp_Wcsp_Lpp: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            SingleStepWingCrossSectionMovement's changes in its WingCrossSections'
            Lp_Wcsp_Lpp parameters. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampLp_Wcsp_Lpp is 0.0 and non zero if not. The
            units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingLp_Wcsp_Lpp: An array-like object of strs or callables with shape
            (3,) representing the spacing of the SingleStepWingCrossSectionMovement's
            changes in its WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple,
            list, or ndarray. Each element can be the str "sine", the str "uniform", or
            a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0.0, return to 0.0 after one period of
            2*pi radians, have amplitude of 1.0, be periodic, return finite values
            only, and accept a ndarray as input and return a ndarray of the same shape.
            Custom functions are scaled by ampLp_Wcsp_Lpp, shifted horizontally and
            vertically by phaseLp_Wcsp_Lpp and the base value, and have a period set by
            periodLp_Wcsp_Lpp. The default is ("sine", "sine", "sine").
        :param phaseLp_Wcsp_Lpp: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's WingCrossSection's Lp_Wcsp_Lpp parameter relative to the base
            WingCrossSection's Lp_Wcsp_Lpp parameter. Can be a tuple, list, or ndarray.
            Elements must lie in the range (-180.0, 180.0]. Each element must be 0.0 if
            the corresponding element in ampLp_Wcsp_Lpp is 0.0 and non zero if not.
            Values are converted to floats internally. The units are in degrees. The
            default is (0.0, 0.0, 0.0).
        :param ampAngles_Wcsp_to_Wcs_ixyz: An array-like object of non negative numbers
            (int or float) with shape (3,) representing the amplitudes of the
            SingleStepWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray.
            Values are converted to floats internally. Each amplitude must be low enough
            that it doesn't drive its base value out of the range of valid values.
            Otherwise, this SingleStepWingCrossSectionMovement will try to create
            WingCrossSections with invalid parameter values. The units are in degrees.
            The default is (0.0, 0.0, 0.0).
        :param periodAngles_Wcsp_to_Wcs_ixyz: An array-like object of non negative
            numbers (int or float) with shape (3,) representing the periods of the
            SingleStepWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray.
            Values are converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non zero if
            not. The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Wcsp_to_Wcs_ixyz: An array-like object of strs or
            callables with shape (3,) representing the spacing of the
            SingleStepWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Each
            element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2*pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a
            ndarray as input and return a ndarray of the same shape. Custom functions
            are scaled by ampAngles_Wcsp_to_Wcs_ixyz, shifted horizontally and
            vertically by phaseAngles_Wcsp_to_Wcs_ixyz and the base value, and have a
            period set by periodAngles_Wcsp_to_Wcs_ixyz. The default is ("sine",
            "sine", "sine").
        :param phaseAngles_Wcsp_to_Wcs_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the phase offsets of the elements in
            the first time step's WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter
            relative to the base WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter.
            Can be a tuple, list, or ndarray. Elements must lie in the range (-180.0,
            180.0]. Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non zero if not. Values are converted
            to floats internally. The units are in degrees. The default is (0.0, 0.0,
            0.0).
        :return: None
        """

        # Warn about potential deformation issues with multiple spanwise panels.
        if base_wing_cross_section.num_spanwise_panels is not None and base_wing_cross_section.num_spanwise_panels > 1:
            print("base_wing_cross_section must have num_spanwise_panels equal to None or 1 to do deformation. " + \
                  "This wing cross section has " + str(base_wing_cross_section.num_spanwise_panels) + " spanwise panels. Please be sure this is intended. " + \
                    "Applications that make sense for this are tails and non-primary wings.")

        # Create the corresponding WingCrossSectionMovement, which validates all
        # oscillation parameters and is also needed by coupled unsteady problems.
        self.corresponding_wing_cross_section_movement = WingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=ampLp_Wcsp_Lpp,
            periodLp_Wcsp_Lpp=periodLp_Wcsp_Lpp,
            spacingLp_Wcsp_Lpp=spacingLp_Wcsp_Lpp,
            phaseLp_Wcsp_Lpp=phaseLp_Wcsp_Lpp,
            ampAngles_Wcsp_to_Wcs_ixyz=ampAngles_Wcsp_to_Wcs_ixyz,
            periodAngles_Wcsp_to_Wcs_ixyz=periodAngles_Wcsp_to_Wcs_ixyz,
            spacingAngles_Wcsp_to_Wcs_ixyz=spacingAngles_Wcsp_to_Wcs_ixyz,
            phaseAngles_Wcsp_to_Wcs_ixyz=phaseAngles_Wcsp_to_Wcs_ixyz,
        )

        # Copy validated attributes from the corresponding WingCrossSectionMovement.
        self.ampLp_Wcsp_Lpp = self.corresponding_wing_cross_section_movement.ampLp_Wcsp_Lpp
        self.periodLp_Wcsp_Lpp = self.corresponding_wing_cross_section_movement.periodLp_Wcsp_Lpp
        self.spacingLp_Wcsp_Lpp = self.corresponding_wing_cross_section_movement.spacingLp_Wcsp_Lpp
        self.phaseLp_Wcsp_Lpp = self.corresponding_wing_cross_section_movement.phaseLp_Wcsp_Lpp
        self.ampAngles_Wcsp_to_Wcs_ixyz = self.corresponding_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz
        self.periodAngles_Wcsp_to_Wcs_ixyz = self.corresponding_wing_cross_section_movement.periodAngles_Wcsp_to_Wcs_ixyz
        self.spacingAngles_Wcsp_to_Wcs_ixyz = self.corresponding_wing_cross_section_movement.spacingAngles_Wcsp_to_Wcs_ixyz
        self.phaseAngles_Wcsp_to_Wcs_ixyz = self.corresponding_wing_cross_section_movement.phaseAngles_Wcsp_to_Wcs_ixyz

        self.listLp_Wcsp_Lpp = None
        self.listAngles_Wcsp_to_Wcs_ixyz = None

    def generate_next_wing_cross_sections(
        self,
        base_wing_cross_section,
        delta_time,
        num_steps,
        step,
        deformation_matrix,
    ):
        """Creates the WingCrossSection at a single time step.

        :param base_wing_cross_section: The base WingCrossSection from which the new
            WingCrossSection will be created.
        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps in this movement.
        :param step: The index of the current time step.
        :param deformation_matrix: Deformation matrix to apply to the
            WingCrossSection's angles, or None.
        :return: A list containing the WingCrossSection at the specified time step.
        """
        num_steps = int_in_range_return_int(
            num_steps, "num_steps", min_val=1, min_inclusive=True
        )
        delta_time = number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )

        # Generate oscillating values for each dimension of Lp_Wcsp_Lpp.
        if self.listLp_Wcsp_Lpp is None:
            self._initialize_oscillating_dimensions(
                delta_time,
                num_steps,
                base_wing_cross_section,
            )

        # Generate oscillating values for each dimension of angles_Wcsp_to_Wcs_ixyz.
        if self.listAngles_Wcsp_to_Wcs_ixyz is None:
            self._initialize_oscillating_angles(
                delta_time,
                num_steps,
                base_wing_cross_section,
            )

        # Create an empty list to hold each time step's WingCrossSection.
        wing_cross_sections = []

        # Get the non-changing WingCrossSectionAttributes.
        this_airfoil = base_wing_cross_section.airfoil
        this_num_spanwise_panels = base_wing_cross_section.num_spanwise_panels
        this_chord = base_wing_cross_section.chord
        this_control_surface_symmetry_type = (
            base_wing_cross_section.control_surface_symmetry_type
        )
        this_control_surface_hinge_point = (
            base_wing_cross_section.control_surface_hinge_point
        )
        this_control_surface_deflection = (
            base_wing_cross_section.control_surface_deflection
        )
        this_spanwise_spacing = base_wing_cross_section.spanwise_spacing

        thisLp_Wcsp_Lpp = self.listLp_Wcsp_Lpp[:, step] 
        theseAngles_Wcsp_to_Wcs_ixyz = self.listAngles_Wcsp_to_Wcs_ixyz[
            :, step
        ] + deformation_matrix

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

    def _initialize_oscillating_dimensions(
        self,
        delta_time,
        num_steps,
        base_wing_cross_section,
    ):
        """Pre computes the oscillating Lp_Wcsp_Lpp values for all time steps.

        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps.
        :param base_wing_cross_section: The base WingCrossSection providing the base
            Lp_Wcsp_Lpp values.
        :return: None
        """

        # Generate oscillating values for each dimension of Lp_Wcsp_Lpp.
        self.listLp_Wcsp_Lpp = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingLp_Wcsp_Lpp[dim]
            if spacing == "sine":
                self.listLp_Wcsp_Lpp[dim, :] = oscillating_sinspaces(
                    amps=self.ampLp_Wcsp_Lpp[dim],
                    periods=self.periodLp_Wcsp_Lpp[dim],
                    phases=self.phaseLp_Wcsp_Lpp[dim],
                    bases=base_wing_cross_section.Lp_Wcsp_Lpp[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listLp_Wcsp_Lpp[dim, :] = oscillating_linspaces(
                    amps=self.ampLp_Wcsp_Lpp[dim],
                    periods=self.periodLp_Wcsp_Lpp[dim],
                    phases=self.phaseLp_Wcsp_Lpp[dim],
                    bases=base_wing_cross_section.Lp_Wcsp_Lpp[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listLp_Wcsp_Lpp[dim, :] = oscillating_customspaces(
                    amps=self.ampLp_Wcsp_Lpp[dim],
                    periods=self.periodLp_Wcsp_Lpp[dim],
                    phases=self.phaseLp_Wcsp_Lpp[dim],
                    bases=base_wing_cross_section.Lp_Wcsp_Lpp[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

    def _initialize_oscillating_angles(
        self,
        delta_time,
        num_steps,
        base_wing_cross_section,
    ):
        """Pre computes the oscillating angles_Wcsp_to_Wcs_ixyz values for all time steps.

        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps.
        :param base_wing_cross_section: The base WingCrossSection providing the base
            angles_Wcsp_to_Wcs_ixyz values.
        :return: None
        """
        self.listAngles_Wcsp_to_Wcs_ixyz = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingAngles_Wcsp_to_Wcs_ixyz[dim]
            if spacing == "sine":
                self.listAngles_Wcsp_to_Wcs_ixyz[dim, :] = oscillating_sinspaces(
                    amps=self.ampAngles_Wcsp_to_Wcs_ixyz[dim],
                    periods=self.periodAngles_Wcsp_to_Wcs_ixyz[dim],
                    phases=self.phaseAngles_Wcsp_to_Wcs_ixyz[dim],
                    bases=base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listAngles_Wcsp_to_Wcs_ixyz[dim, :] = oscillating_linspaces(
                    amps=self.ampAngles_Wcsp_to_Wcs_ixyz[dim],
                    periods=self.periodAngles_Wcsp_to_Wcs_ixyz[dim],
                    phases=self.phaseAngles_Wcsp_to_Wcs_ixyz[dim],
                    bases=base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listAngles_Wcsp_to_Wcs_ixyz[dim, :] = oscillating_customspaces(
                    amps=self.ampAngles_Wcsp_to_Wcs_ixyz[dim],
                    periods=self.periodAngles_Wcsp_to_Wcs_ixyz[dim],
                    phases=self.phaseAngles_Wcsp_to_Wcs_ixyz[dim],
                    bases=base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")
