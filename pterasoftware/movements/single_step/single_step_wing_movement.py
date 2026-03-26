"""Contains the SingleStepWingMovement class.

**Contains the following classes:**

SingleStepWingMovement: A single step variant of WingMovement that generates one Wing
per time step instead of all at once. Uses composition to wrap a WingMovement.

**Contains the following functions:**

None
"""

import numpy as np

from ... import geometry
from ..._parameter_validation import (
    int_in_range_return_int,
    number_in_range_return_float,
)
from .._functions import (
    oscillating_customspaces,
    oscillating_linspaces,
    oscillating_sinspaces,
)
from ..wing_movement import WingMovement


class SingleStepWingMovement:
    """A single step variant of WingMovement for coupled simulations.

    This class wraps a WingMovement via composition and generates one Wing per time step
    (via generate_next_wing) rather than generating all Wings at once. The composed
    WingMovement is accessible via corresponding_wing_movement.

    Wings cannot undergo motion that causes them to switch symmetry types. See the
    WingMovement class documentation for details.

    **Contains the following methods:**

    generate_next_wing: Creates the Wing at a single time step.

    max_period: The longest period of this SingleStepWingMovement's own motion and that
    of its sub movement objects.
    """

    __slots__ = (
        "wing_cross_section_movements",
        "ampLer_Gs_Cgs",
        "periodLer_Gs_Cgs",
        "spacingLer_Gs_Cgs",
        "phaseLer_Gs_Cgs",
        "ampAngles_Gs_to_Wn_ixyz",
        "periodAngles_Gs_to_Wn_ixyz",
        "spacingAngles_Gs_to_Wn_ixyz",
        "phaseAngles_Gs_to_Wn_ixyz",
        "listLer_Gs_Cgs",
        "listAngles_Gs_to_Wn_ixyz",
        "corresponding_wing_movement",
    )

    def __init__(
        self,
        base_wing,
        single_step_wing_cross_section_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    ):
        """The initialization method.

        :param base_wing: The base Wing from which the Wing at each time step will be
            created.
        :param single_step_wing_cross_section_movements: A list of
            SingleStepWingCrossSectionMovements associated with each of the base Wing's
            WingCrossSections. It must have the same length as the base Wing's list of
            WingCrossSections.
        :param ampLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            SingleStepWingMovement's changes in its Wings' Ler_Gs_Cgs parameters. Can be
            a tuple, list, or ndarray. Values are converted to floats internally. Each
            amplitude must be low enough that it doesn't drive its base value out of the
            range of valid values. Otherwise, this SingleStepWingMovement will try to
            create Wings with invalid parameter values. The units are in meters. The
            default is (0.0, 0.0, 0.0).
        :param periodLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            SingleStepWingMovement's changes in its Wings' Ler_Gs_Cgs parameters. Can be
            a tuple, list, or ndarray. Values are converted to floats internally. Each
            element must be 0.0 if the corresponding element in ampLer_Gs_Cgs is 0.0 and
            non zero if not. The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingLer_Gs_Cgs: An array-like object of strs or callables with shape
            (3,) representing the spacing of the SingleStepWingMovement's change in its
            Wings' Ler_Gs_Cgs parameters. Can be a tuple, list, or ndarray. Each element
            can be the str "sine", the str "uniform", or a callable custom spacing
            function. Custom spacing functions are for advanced users and must start at
            0.0, return to 0.0 after one period of 2*pi radians, have amplitude of 1.0,
            be periodic, return finite values only, and accept a ndarray as input and
            return a ndarray of the same shape. The custom function is scaled by
            ampLer_Gs_Cgs, shifted horizontally and vertically by phaseLer_Gs_Cgs and
            the base value, and have a period set by periodLer_Gs_Cgs. The default is
            ("sine", "sine", "sine").
        :param phaseLer_Gs_Cgs: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's Wing's Ler_Gs_Cgs parameter relative to the base Wing's Ler_Gs_Cgs
            parameter. Can be a tuple, list, or ndarray. Values must lie in the range
            (-180.0, 180.0] and will be converted to floats internally. Each element
            must be 0.0 if the corresponding element in ampLer_Gs_Cgs is 0.0 and non
            zero if not. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param ampAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or float)
            with shape (3,) representing the amplitudes of the SingleStepWingMovement's
            changes in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Values must lie in the range [0.0, 180.0] and will be converted
            to floats internally. Each amplitude must be low enough that it doesn't
            drive its base value out of the range of valid values. Otherwise, this
            SingleStepWingMovement will try to create Wings with invalid parameter
            values. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param periodAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the periods of the
            SingleStepWingMovement's changes in its Wings' angles_Gs_to_Wn_ixyz
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampAngles_Gs_to_Wn_ixyz is 0.0 and non zero if not. The units are in
            seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Gs_to_Wn_ixyz: An array-like object of strs or callables
            with shape (3,) representing the spacing of the SingleStepWingMovement's
            change in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Each element can be the str "sine", the str "uniform", or a
            callable custom spacing function. Custom spacing functions are for advanced
            users and must start at 0.0, return to 0.0 after one period of 2*pi radians,
            have amplitude of 1.0, be periodic, return finite values only, and accept a
            ndarray as input and return a ndarray of the same shape. The custom function
            is scaled by ampAngles_Gs_to_Wn_ixyz, shifted horizontally and vertically by
            phaseAngles_Gs_to_Wn_ixyz and the base value, with the period set by
            periodAngles_Gs_to_Wn_ixyz. The default is ("sine", "sine", "sine").
        :param phaseAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or float)
            with shape (3,) representing the phase offsets of the elements in the first
            time step's Wing's angles_Gs_to_Wn_ixyz parameter relative to the base
            Wing's angles_Gs_to_Wn_ixyz parameter. Can be a tuple, list, or ndarray.
            Values must lie in the range (-180.0, 180.0] and will be converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampAngles_Gs_to_Wn_ixyz is 0.0 and non zero if not. The units are in
            degrees. The default is (0.0, 0.0, 0.0).
        :return: None
        """
        self.wing_cross_section_movements = single_step_wing_cross_section_movements

        # Create the corresponding WingMovement, which validates all oscillation
        # parameters and is also needed by coupled unsteady problems.
        corresponding_wing_cross_section_movements = [
            wing_cross_section_movement.corresponding_wing_cross_section_movement
            for wing_cross_section_movement in single_step_wing_cross_section_movements
        ]
        self.corresponding_wing_movement = WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=corresponding_wing_cross_section_movements,
            ampLer_Gs_Cgs=ampLer_Gs_Cgs,
            periodLer_Gs_Cgs=periodLer_Gs_Cgs,
            spacingLer_Gs_Cgs=spacingLer_Gs_Cgs,
            phaseLer_Gs_Cgs=phaseLer_Gs_Cgs,
            ampAngles_Gs_to_Wn_ixyz=ampAngles_Gs_to_Wn_ixyz,
            periodAngles_Gs_to_Wn_ixyz=periodAngles_Gs_to_Wn_ixyz,
            spacingAngles_Gs_to_Wn_ixyz=spacingAngles_Gs_to_Wn_ixyz,
            phaseAngles_Gs_to_Wn_ixyz=phaseAngles_Gs_to_Wn_ixyz,
        )

        # Copy validated attributes from the corresponding WingMovement.
        self.ampLer_Gs_Cgs = self.corresponding_wing_movement.ampLer_Gs_Cgs
        self.periodLer_Gs_Cgs = self.corresponding_wing_movement.periodLer_Gs_Cgs
        self.spacingLer_Gs_Cgs = self.corresponding_wing_movement.spacingLer_Gs_Cgs
        self.phaseLer_Gs_Cgs = self.corresponding_wing_movement.phaseLer_Gs_Cgs
        self.ampAngles_Gs_to_Wn_ixyz = (
            self.corresponding_wing_movement.ampAngles_Gs_to_Wn_ixyz
        )
        self.periodAngles_Gs_to_Wn_ixyz = (
            self.corresponding_wing_movement.periodAngles_Gs_to_Wn_ixyz
        )
        self.spacingAngles_Gs_to_Wn_ixyz = (
            self.corresponding_wing_movement.spacingAngles_Gs_to_Wn_ixyz
        )
        self.phaseAngles_Gs_to_Wn_ixyz = (
            self.corresponding_wing_movement.phaseAngles_Gs_to_Wn_ixyz
        )

        self.listLer_Gs_Cgs = None
        self.listAngles_Gs_to_Wn_ixyz = None

    def generate_next_wing(
        self, base_wing, delta_time, num_steps, step, deformation_matrices
    ):
        """Creates the Wing at a single time step.

        :param base_wing: The base Wing from which the new Wing will be created.
        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps in this movement.
        :param step: The index of the current time step.
        :param deformation_matrices: Deformation matrices to apply to the
            WingCrossSections, or None.
        :return: The Wing at the specified time step.
        """
        num_steps = int_in_range_return_int(
            num_steps, "num_steps", min_val=1, min_inclusive=True
        )
        delta_time = number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )
        # Account for null deformation_matrices input.
        if deformation_matrices is None:
            deformation_matrices = np.zeros(len(self.wing_cross_section_movements))

        # Generate oscillating values for each dimension of Ler_Gs_Cgs.
        if self.listLer_Gs_Cgs is None:
            self._initialize_oscillating_dimensions(delta_time, num_steps, base_wing)

        # Generate oscillating values for each dimension of angles_Gs_to_Wn_ixyz.
        if self.listAngles_Gs_to_Wn_ixyz is None:
            self._initialize_oscillating_angles(delta_time, num_steps, base_wing)

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
                wing_cross_section_movement.generate_next_wing_cross_sections(
                    base_wing_cross_section=base_wing.wing_cross_sections[
                        wing_cross_section_movement_id
                    ],
                    delta_time=delta_time,
                    num_steps=num_steps,
                    step=step,
                    deformation_matrix=deformation_matrices[
                        wing_cross_section_movement_id
                    ],
                )
            )

            # Add this vector the Wing's 2D ndarray of WingCrossSections'
            # WingCrossSections.
            wing_cross_sections[wing_cross_section_movement_id, :] = (
                this_wing_cross_sections_list_of_wing_cross_sections
            )

        # Get the non-changing Wing attributes.
        this_name = base_wing.name
        this_symmetric = base_wing.symmetric
        this_mirror_only = base_wing.mirror_only
        this_symmetryNormal_G = base_wing.symmetryNormal_G
        this_symmetryPoint_G_Cg = base_wing.symmetryPoint_G_Cg
        this_num_chordwise_panels = base_wing.num_chordwise_panels
        this_chordwise_spacing = base_wing.chordwise_spacing

        thisLer_Gs_Cgs = self.listLer_Gs_Cgs[:, step]
        theseAngles_Gs_to_Wn_ixyz = self.listAngles_Gs_to_Wn_ixyz[:, step]
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

        return this_wing

    def _initialize_oscillating_dimensions(self, delta_time, num_steps, base_wing):
        """Pre computes the oscillating Ler_Gs_Cgs values for all time steps.

        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps.
        :param base_wing: The base Wing providing the base Ler_Gs_Cgs values.
        :return: None
        """
        self.listLer_Gs_Cgs = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingLer_Gs_Cgs[dim]
            if spacing == "sine":
                self.listLer_Gs_Cgs[dim, :] = oscillating_sinspaces(
                    amps=self.ampLer_Gs_Cgs[dim],
                    periods=self.periodLer_Gs_Cgs[dim],
                    phases=self.phaseLer_Gs_Cgs[dim],
                    bases=base_wing.Ler_Gs_Cgs[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listLer_Gs_Cgs[dim, :] = oscillating_linspaces(
                    amps=self.ampLer_Gs_Cgs[dim],
                    periods=self.periodLer_Gs_Cgs[dim],
                    phases=self.phaseLer_Gs_Cgs[dim],
                    bases=base_wing.Ler_Gs_Cgs[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listLer_Gs_Cgs[dim, :] = oscillating_customspaces(
                    amps=self.ampLer_Gs_Cgs[dim],
                    periods=self.periodLer_Gs_Cgs[dim],
                    phases=self.phaseLer_Gs_Cgs[dim],
                    bases=base_wing.Ler_Gs_Cgs[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")

    def _initialize_oscillating_angles(self, delta_time, num_steps, base_wing):
        """Pre computes the oscillating angles_Gs_to_Wn_ixyz values for all time steps.

        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps.
        :param base_wing: The base Wing providing the base angles_Gs_to_Wn_ixyz values.
        :return: None
        """
        self.listAngles_Gs_to_Wn_ixyz = np.zeros((3, num_steps), dtype=float)
        for dim in range(3):
            spacing = self.spacingAngles_Gs_to_Wn_ixyz[dim]
            if spacing == "sine":
                self.listAngles_Gs_to_Wn_ixyz[dim, :] = oscillating_sinspaces(
                    amps=self.ampAngles_Gs_to_Wn_ixyz[dim],
                    periods=self.periodAngles_Gs_to_Wn_ixyz[dim],
                    phases=self.phaseAngles_Gs_to_Wn_ixyz[dim],
                    bases=base_wing.angles_Gs_to_Wn_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif spacing == "uniform":
                self.listAngles_Gs_to_Wn_ixyz[dim, :] = oscillating_linspaces(
                    amps=self.ampAngles_Gs_to_Wn_ixyz[dim],
                    periods=self.periodAngles_Gs_to_Wn_ixyz[dim],
                    phases=self.phaseAngles_Gs_to_Wn_ixyz[dim],
                    bases=base_wing.angles_Gs_to_Wn_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                )
            elif callable(spacing):
                self.listAngles_Gs_to_Wn_ixyz[dim, :] = oscillating_customspaces(
                    amps=self.ampAngles_Gs_to_Wn_ixyz[dim],
                    periods=self.periodAngles_Gs_to_Wn_ixyz[dim],
                    phases=self.phaseAngles_Gs_to_Wn_ixyz[dim],
                    bases=base_wing.angles_Gs_to_Wn_ixyz[dim],
                    num_steps=num_steps,
                    delta_time=delta_time,
                    custom_function=spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {spacing}")
