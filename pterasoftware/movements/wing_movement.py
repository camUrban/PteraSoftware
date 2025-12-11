"""Contains the WingMovement class.

**Contains the following classes:**

WingMovement: A class used to contain a Wing's movement.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence, Callable

import numpy as np

from . import _functions
from . import wing_cross_section_movement as wing_cross_section_movement_mod

from .. import geometry
from .. import _parameter_validation


class WingMovement:
    """A class used to contain a Wing's movement.

    **Contains the following methods:**

    all_periods: All unique non zero periods from this WingMovement and its
    WingCrossSectionMovements.

    generate_wings: Creates the Wing at each time step, and returns them in a list.

    max_period: The longest period of WingMovement's own motion and that of its sub
    movement objects.

    **Notes:**

    Wings cannot undergo motion that causes them to switch symmetry types. A transition
    between types could change the number of Wings and the Panel structure, which is
    incompatible with the unsteady solver. This happens when a WingMovement defines
    motion that causes its base Wing's wing axes' yz plane and its symmetry plane to
    transition from coincident to non coincident, or vice versa. This is checked by this
    WingMovement's parent AirplaneMovement's parent Movement.
    """

    def __init__(
        self,
        base_wing: geometry.wing.Wing,
        wing_cross_section_movements: list[
            wing_cross_section_movement_mod.WingCrossSectionMovement
        ],
        ampLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs: (
            np.ndarray | Sequence[str | Callable[[np.ndarray], np.ndarray]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        spacingAngles_Gs_to_Wn_ixyz: (
            np.ndarray | Sequence[str | Callable[[np.ndarray], np.ndarray]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
    ) -> None:
        """The initialization method.

        :param base_wing: The base Wing from which the Wing at each time step will be
            created.
        :param wing_cross_section_movements: A list of WingCrossSectionMovements
            associated with each of the base Wing's WingCrossSections. It must have the
            same length as the base Wing's list of WingCrossSections.
        :param ampLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the WingMovement's
            changes in its Wings' Ler_Gs_Cgs parameters. Can be a tuple, list, or
            ndarray. Values are converted to floats internally. Each amplitude must be
            low enough that it doesn't drive its base value out of the range of valid
            values. Otherwise, this WingMovement will try to create Wings with invalid
            parameters values. The units are in meters. The default is (0.0, 0.0, 0.0).
        :param periodLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the WingMovement's
            changes in its Wings' Ler_Gs_Cgs parameters. Can be a tuple, list, or
            ndarray. Values are converted to floats internally. Each element must be 0.0
            if the corresponding element in ampLer_Gs_Cgs is 0.0 and non zero if not.
            The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingLer_Gs_Cgs: An array-like object of strs or callables with shape
            (3,) representing the spacing of the WingMovement's change in its Wings'
            Ler_Gs_Cgs parameters. Can be a tuple, list, or ndarray. Each element can be
            the string "sine", the string "uniform", or a callable custom spacing
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
            with shape (3,) representing the amplitudes of the WingMovement's changes in
            its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list, or
            ndarray. Values must lie in the range [0.0, 180.0] and will be converted to
            floats internally. Each amplitude must be low enough that it doesn't drive
            its base value out of the range of valid values. Otherwise, this
            WingMovement will try to create Wings with invalid parameters values. The
            units are in degrees. The default is (0.0, 0.0, 0.0).
        :param periodAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the periods of the WingMovement's
            changes in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Values are converted to floats internally. Each element must be
            0.0 if the corresponding element in ampAngles_Gs_to_Wn_ixyz is 0.0 and non
            zero if not. The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Gs_to_Wn_ixyz: An array-like object of strs or callables
            with shape (3,) representing the spacing of the WingMovement's change in its
            Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list, or ndarray.
            Each element can be the string "sine", the string "uniform", or a callable
            custom spacing function. Custom spacing functions are for advanced users and
            must start at 0.0, return to 0.0 after one period of 2*pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a
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
        """
        if not isinstance(base_wing, geometry.wing.Wing):
            raise TypeError("base_wing must be a Wing.")
        self.base_wing = base_wing

        if not isinstance(wing_cross_section_movements, list):
            raise TypeError("wing_cross_section_movements must be a list.")
        if len(wing_cross_section_movements) != len(self.base_wing.wing_cross_sections):
            raise ValueError(
                "wing_cross_section_movements must have the same length as "
                "base_wing.wing_cross_sections."
            )
        for wing_cross_section_movement in wing_cross_section_movements:
            if not isinstance(
                wing_cross_section_movement,
                wing_cross_section_movement_mod.WingCrossSectionMovement,
            ):
                raise TypeError(
                    "Every element in wing_cross_section_movements must be a "
                    "WingCrossSectionMovement."
                )
        self.wing_cross_section_movements = wing_cross_section_movements

        ampLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            ampLer_Gs_Cgs, "ampLer_Gs_Cgs"
        )
        if not np.all(ampLer_Gs_Cgs >= 0.0):
            raise ValueError("All elements in ampLer_Gs_Cgs must be non negative.")
        self.ampLer_Gs_Cgs = ampLer_Gs_Cgs

        periodLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            periodLer_Gs_Cgs, "periodLer_Gs_Cgs"
        )
        if not np.all(periodLer_Gs_Cgs >= 0.0):
            raise ValueError("All elements in periodLer_Gs_Cgs must be non negative.")
        for period_index, period in enumerate(periodLer_Gs_Cgs):
            amp = self.ampLer_Gs_Cgs[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampLer_Gs_Cgs is 0.0, the corresponding element "
                    "in periodLer_Gs_Cgs must be also be 0.0."
                )
        self.periodLer_Gs_Cgs = periodLer_Gs_Cgs

        spacingLer_Gs_Cgs = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingLer_Gs_Cgs, "spacingLer_Gs_Cgs"
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
                    "If an element in ampLer_Gs_Cgs is 0.0, the corresponding element "
                    "in phaseLer_Gs_Cgs must be also be 0.0."
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
                "All elements in ampAngles_Gs_to_Wn_ixyz must be in the range [0.0, "
                "180.0]."
            )
        self.ampAngles_Gs_to_Wn_ixyz = ampAngles_Gs_to_Wn_ixyz

        periodAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_Gs_to_Wn_ixyz, "periodAngles_Gs_to_Wn_ixyz"
            )
        )
        if not np.all(periodAngles_Gs_to_Wn_ixyz >= 0.0):
            raise ValueError(
                "All elements in periodAngles_Gs_to_Wn_ixyz must be non negative."
            )
        for period_index, period in enumerate(periodAngles_Gs_to_Wn_ixyz):
            amp = self.ampAngles_Gs_to_Wn_ixyz[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_Gs_to_Wn_ixyz is 0.0, "
                    "the corresponding element in periodAngles_Gs_to_Wn_ixyz must be "
                    "also be 0.0."
                )
        self.periodAngles_Gs_to_Wn_ixyz = periodAngles_Gs_to_Wn_ixyz

        spacingAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingAngles_Gs_to_Wn_ixyz, "spacingAngles_Gs_to_Wn_ixyz"
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
                "All elements in phaseAngles_Gs_to_Wn_ixyz must be in the range ("
                "-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_Gs_to_Wn_ixyz):
            amp = self.ampAngles_Gs_to_Wn_ixyz[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_Gs_to_Wn_ixyz is 0.0, "
                    "the corresponding element in phaseAngles_Gs_to_Wn_ixyz must be "
                    "also be 0.0."
                )
        self.phaseAngles_Gs_to_Wn_ixyz = phaseAngles_Gs_to_Wn_ixyz

    @property
    def all_periods(self) -> list[float]:
        """All unique non zero periods from this WingMovement and its
        WingCrossSectionMovements.

        :return: A list of all unique non zero periods in seconds. If all motion is
            static, this will be an empty list.
        """
        periods = []

        # Collect all periods from WingCrossSectionMovements
        for wing_cross_section_movement in self.wing_cross_section_movements:
            periods.extend(wing_cross_section_movement.all_periods)

        # Collect all periods from WingMovement's own motion
        for period in self.periodLer_Gs_Cgs:
            if period > 0.0:
                periods.append(float(period))
        for period in self.periodAngles_Gs_to_Wn_ixyz:
            if period > 0.0:
                periods.append(float(period))
        return periods

    def generate_wings(
        self, num_steps: int, delta_time: float | int
    ) -> list[geometry.wing.Wing]:
        """Creates the Wing at each time step, and returns them in a list.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (int or float), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of Wings associated with this WingMovement.
        """
        num_steps = _parameter_validation.int_in_range_return_int(
            num_steps,
            "num_steps",
            min_val=1,
            min_inclusive=True,
        )
        delta_time = _parameter_validation.number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
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

        # Get the non changing Wing attributes.
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
    def max_period(self) -> float:
        """The longest period of WingMovement's own motion and that of its sub movement
        objects.

        :return: The longest period in seconds. If all the motion is static, this will
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
