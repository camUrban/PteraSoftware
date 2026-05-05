"""Contains the core classes for the movement and problem hierarchies."""

from __future__ import annotations

import copy
import math
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING

import numpy as np

from . import _oscillation, _parameter_validation, _transformations, geometry
from . import operating_point as operating_point_mod

if TYPE_CHECKING:
    from . import problems


def lcm(a: float, b: float) -> float:
    """Calculates the least common multiple of two numbers.

    :param a: First number (period in seconds).
    :param b: Second number (period in seconds).
    :return: LCM of a and b. Returns 0.0 if either input is 0.0.
    """
    if a == 0.0 or b == 0.0:
        return 0.0
    # Convert to integers (periods are typically whole multiples of delta_time).
    # Use sufficiently large multiplier to preserve precision.
    multiplier = 1000000
    a_int = int(round(a * multiplier))
    b_int = int(round(b * multiplier))
    lcm_int = abs(a_int * b_int) // math.gcd(a_int, b_int)
    return lcm_int / multiplier


def lcm_multiple(periods: list[float]) -> float:
    """Calculates the least common multiple of multiple periods.

    :param periods: A list of periods in seconds.
    :return: LCM of all periods. Returns 0.0 if all periods are 0.0.
    """
    if not periods or all(p == 0.0 for p in periods):
        return 0.0
    non_zero_periods = [p for p in periods if p != 0.0]
    if not non_zero_periods:
        return 0.0
    result = non_zero_periods[0]
    for period in non_zero_periods[1:]:
        result = lcm(result, period)
    return result


class CoreOperatingPointMovement:
    """A core class used to contain the shared foundation of OperatingPointMovement and
    its feature variant siblings.

    See OperatingPointMovement for full documentation of the shared interface.

    CoreOperatingPointMovement holds the base OperatingPoint and oscillation parameters,
    and provides generate_operating_point_at_time_step() for creating OperatingPoints
    one step at a time.
    """

    __slots__ = (
        "_base_operating_point",
        "_ampVCg__E",
        "_periodVCg__E",
        "_spacingVCg__E",
        "_phaseVCg__E",
        "_max_period",
    )

    def __init__(
        self,
        base_operating_point: operating_point_mod.OperatingPoint,
        ampVCg__E: float | int = 0.0,
        periodVCg__E: float | int = 0.0,
        spacingVCg__E: str | Callable[[float], float] = "sine",
        phaseVCg__E: float | int = 0.0,
    ) -> None:
        """The initialization method.

        See OperatingPointMovement's initialization method for full parameter
        descriptions.

        :param base_operating_point: The base OperatingPoint.
        :param ampVCg__E: The amplitude of vCg__E oscillation in meters per second.
        :param periodVCg__E: The period of vCg__E oscillation in seconds.
        :param spacingVCg__E: The spacing type: "sine", "uniform", or a callable.
        :param phaseVCg__E: The phase offset in degrees.
        :return: None
        """
        # Validate and store immutable attributes.
        if not isinstance(base_operating_point, operating_point_mod.OperatingPoint):
            raise TypeError("base_operating_point must be an OperatingPoint")
        self._base_operating_point = base_operating_point

        self._ampVCg__E = _parameter_validation.number_in_range_return_float(
            ampVCg__E, "ampVCg__E", min_val=0.0, min_inclusive=True
        )

        periodVCg__E = _parameter_validation.number_in_range_return_float(
            periodVCg__E, "periodVCg__E", min_val=0.0, min_inclusive=True
        )
        if self._ampVCg__E == 0 and periodVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then periodVCg__E must also be 0.0.")
        self._periodVCg__E = periodVCg__E

        if isinstance(spacingVCg__E, str):
            if spacingVCg__E not in ["sine", "uniform"]:
                raise ValueError(
                    f"spacingVCg__E must be 'sine', 'uniform', or a callable, "
                    f"got string '{spacingVCg__E}'."
                )
        elif not callable(spacingVCg__E):
            raise TypeError(
                f"spacingVCg__E must be 'sine', 'uniform', or a callable, got "
                f"{type(spacingVCg__E).__name__}."
            )
        self._spacingVCg__E = spacingVCg__E

        phaseVCg__E = _parameter_validation.number_in_range_return_float(
            phaseVCg__E, "phaseVCg__E", -180.0, False, 180.0, True
        )
        if self._ampVCg__E == 0 and phaseVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then phaseVCg__E must also be 0.0.")
        self._phaseVCg__E = phaseVCg__E

        # Initialize the cache for the property derived from the immutable
        # attributes.
        self._max_period: float | None = None

    # --- Immutable: read only properties ---
    @property
    def base_operating_point(self) -> operating_point_mod.OperatingPoint:
        return self._base_operating_point

    @property
    def ampVCg__E(self) -> float:
        return self._ampVCg__E

    @property
    def periodVCg__E(self) -> float:
        return self._periodVCg__E

    @property
    def spacingVCg__E(self) -> str | Callable[[float], float]:
        return self._spacingVCg__E

    @property
    def phaseVCg__E(self) -> float:
        return self._phaseVCg__E

    # --- Immutable derived: manual lazy caching ---
    @property
    def max_period(self) -> float:
        """The longest period of this OperatingPoint's motion.

        :return: The longest period in seconds. If the motion is static, this will be
            0.0.
        """
        if self._max_period is None:
            self._max_period = self._periodVCg__E
        return self._max_period

    # --- Other methods ---
    def generate_operating_point_at_time_step(
        self, step: int, delta_time: float | int
    ) -> operating_point_mod.OperatingPoint:
        """Creates the OperatingPoint at a single time step.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :return: The OperatingPoint at this time step.
        """
        time = step * delta_time

        # Evaluate the oscillating function for VCg__E.
        if self._spacingVCg__E == "sine":
            this_vCg__E = _oscillation.oscillating_sin_at_time(
                amp=self._ampVCg__E,
                period=self._periodVCg__E,
                phase=self._phaseVCg__E,
                base=self._base_operating_point.vCg__E,
                time=time,
            )
        elif self._spacingVCg__E == "uniform":
            this_vCg__E = _oscillation.oscillating_lin_at_time(
                amp=self._ampVCg__E,
                period=self._periodVCg__E,
                phase=self._phaseVCg__E,
                base=self._base_operating_point.vCg__E,
                time=time,
            )
        elif callable(self._spacingVCg__E):
            this_vCg__E = _oscillation.oscillating_custom_at_time(
                amp=self._ampVCg__E,
                period=self._periodVCg__E,
                phase=self._phaseVCg__E,
                base=self._base_operating_point.vCg__E,
                time=time,
                custom_function=self._spacingVCg__E,
            )
        else:
            raise ValueError(f"Invalid spacing value: {self._spacingVCg__E}")

        return operating_point_mod.OperatingPoint(
            rho=self._base_operating_point.rho,
            vCg__E=this_vCg__E,
            alpha=self._base_operating_point.alpha,
            beta=self._base_operating_point.beta,
            externalFX_W=self._base_operating_point.externalFX_W,
            nu=self._base_operating_point.nu,
            angles_E_to_BP1_izyx=self._base_operating_point.angles_E_to_BP1_izyx,
            CgP1_E_Eo=self._base_operating_point.CgP1_E_Eo,
            surfaceNormal_E=self._base_operating_point.surfaceNormal_E,
            surfacePoint_E_Eo=self._base_operating_point.surfacePoint_E_Eo,
        )

    def generate_operating_points(
        self, num_steps: int, delta_time: float | int
    ) -> list[operating_point_mod.OperatingPoint]:
        """Creates the OperatingPoint at each time step, and returns them in a list.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (int or float), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of OperatingPoints associated with this
            CoreOperatingPointMovement.
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

        operating_points = []
        for step in range(num_steps):
            operating_points.append(
                self.generate_operating_point_at_time_step(step, delta_time)
            )

        return operating_points


class CoreWingCrossSectionMovement:
    """A core class used to contain the shared foundation of WingCrossSectionMovement
    and its feature variant siblings.

    See WingCrossSectionMovement for full documentation of the shared interface.

    CoreWingCrossSectionMovement holds the base WingCrossSection and oscillation
    parameters, and provides generate_wing_cross_section_at_time_step() for creating
    WingCrossSections one step at a time.
    """

    __slots__ = (
        "_base_wing_cross_section",
        "_ampLp_Wcsp_Lpp",
        "_periodLp_Wcsp_Lpp",
        "_spacingLp_Wcsp_Lpp",
        "_phaseLp_Wcsp_Lpp",
        "_ampAngles_Wcsp_to_Wcs_ixyz",
        "_periodAngles_Wcsp_to_Wcs_ixyz",
        "_spacingAngles_Wcsp_to_Wcs_ixyz",
        "_phaseAngles_Wcsp_to_Wcs_ixyz",
        "_all_periods",
        "_max_period",
    )

    def __init__(
        self,
        base_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
        ampLp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseLp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        periodAngles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        spacingAngles_Wcsp_to_Wcs_ixyz: (
            np.ndarray | Sequence[str | Callable[[float], float]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseAngles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
    ) -> None:
        """The initialization method.

        See WingCrossSectionMovement's initialization method for full parameter
        descriptions.

        :param base_wing_cross_section: The base WingCrossSection.
        :param ampLp_Wcsp_Lpp: The amplitudes of Lp_Wcsp_Lpp oscillation in meters.
        :param periodLp_Wcsp_Lpp: The periods of Lp_Wcsp_Lpp oscillation in seconds.
        :param spacingLp_Wcsp_Lpp: The spacing types for Lp_Wcsp_Lpp oscillation.
        :param phaseLp_Wcsp_Lpp: The phase offsets of Lp_Wcsp_Lpp oscillation in
            degrees.
        :param ampAngles_Wcsp_to_Wcs_ixyz: The amplitudes of angles_Wcsp_to_Wcs_ixyz
            oscillation in degrees.
        :param periodAngles_Wcsp_to_Wcs_ixyz: The periods of angles_Wcsp_to_Wcs_ixyz
            oscillation in seconds.
        :param spacingAngles_Wcsp_to_Wcs_ixyz: The spacing types for
            angles_Wcsp_to_Wcs_ixyz oscillation.
        :param phaseAngles_Wcsp_to_Wcs_ixyz: The phase offsets of
            angles_Wcsp_to_Wcs_ixyz oscillation in degrees.
        :return: None
        """
        # Validate and store immutable attributes. Set those that are numpy
        # arrays to be read only.
        if not isinstance(
            base_wing_cross_section,
            geometry.wing_cross_section.WingCrossSection,
        ):
            raise TypeError("base_wing_cross_section must be a WingCrossSection.")
        self._base_wing_cross_section = base_wing_cross_section

        ampLp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            ampLp_Wcsp_Lpp, "ampLp_Wcsp_Lpp"
        )
        if not np.all(ampLp_Wcsp_Lpp >= 0.0):
            raise ValueError("All elements in ampLp_Wcsp_Lpp must be non negative.")
        self._ampLp_Wcsp_Lpp = ampLp_Wcsp_Lpp
        self._ampLp_Wcsp_Lpp.flags.writeable = False

        periodLp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            periodLp_Wcsp_Lpp, "periodLp_Wcsp_Lpp"
        )
        if not np.all(periodLp_Wcsp_Lpp >= 0.0):
            raise ValueError("All elements in periodLp_Wcsp_Lpp must be non negative.")
        for period_index, period in enumerate(periodLp_Wcsp_Lpp):
            amp = self._ampLp_Wcsp_Lpp[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampLp_Wcsp_Lpp is 0.0, the "
                    "corresponding element in periodLp_Wcsp_Lpp must be "
                    "also be 0.0."
                )
        self._periodLp_Wcsp_Lpp = periodLp_Wcsp_Lpp
        self._periodLp_Wcsp_Lpp.flags.writeable = False

        # Store as tuple to prevent external mutation.
        self._spacingLp_Wcsp_Lpp = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingLp_Wcsp_Lpp, "spacingLp_Wcsp_Lpp"
            )
        )

        phaseLp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            phaseLp_Wcsp_Lpp, "phaseLp_Wcsp_Lpp"
        )
        if not (
            np.all(phaseLp_Wcsp_Lpp > -180.0) and np.all(phaseLp_Wcsp_Lpp <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseLp_Wcsp_Lpp must be in the range "
                "(-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseLp_Wcsp_Lpp):
            amp = self._ampLp_Wcsp_Lpp[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampLp_Wcsp_Lpp is 0.0, the "
                    "corresponding element in phaseLp_Wcsp_Lpp must be "
                    "also be 0.0."
                )
        self._phaseLp_Wcsp_Lpp = phaseLp_Wcsp_Lpp
        self._phaseLp_Wcsp_Lpp.flags.writeable = False

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
                "All elements in ampAngles_Wcsp_to_Wcs_ixyz must be in "
                "the range [0.0, 180.0]."
            )
        self._ampAngles_Wcsp_to_Wcs_ixyz = ampAngles_Wcsp_to_Wcs_ixyz
        self._ampAngles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        periodAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_Wcsp_to_Wcs_ixyz,
                "periodAngles_Wcsp_to_Wcs_ixyz",
            )
        )
        if not np.all(periodAngles_Wcsp_to_Wcs_ixyz >= 0.0):
            raise ValueError(
                "All elements in periodAngles_Wcsp_to_Wcs_ixyz must be " "non negative."
            )
        for period_index, period in enumerate(periodAngles_Wcsp_to_Wcs_ixyz):
            amp = self._ampAngles_Wcsp_to_Wcs_ixyz[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0, "
                    "the corresponding element in "
                    "periodAngles_Wcsp_to_Wcs_ixyz must be also be 0.0."
                )
        self._periodAngles_Wcsp_to_Wcs_ixyz = periodAngles_Wcsp_to_Wcs_ixyz
        self._periodAngles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        # Store as tuple to prevent external mutation.
        self._spacingAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingAngles_Wcsp_to_Wcs_ixyz,
                "spacingAngles_Wcsp_to_Wcs_ixyz",
            )
        )

        phaseAngles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                phaseAngles_Wcsp_to_Wcs_ixyz,
                "phaseAngles_Wcsp_to_Wcs_ixyz",
            )
        )
        if not (
            np.all(phaseAngles_Wcsp_to_Wcs_ixyz > -180.0)
            and np.all(phaseAngles_Wcsp_to_Wcs_ixyz <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_Wcsp_to_Wcs_ixyz must be in "
                "the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_Wcsp_to_Wcs_ixyz):
            amp = self._ampAngles_Wcsp_to_Wcs_ixyz[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0, "
                    "the corresponding element in "
                    "phaseAngles_Wcsp_to_Wcs_ixyz must be also be 0.0."
                )
        self._phaseAngles_Wcsp_to_Wcs_ixyz = phaseAngles_Wcsp_to_Wcs_ixyz
        self._phaseAngles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._all_periods: tuple[float, ...] | None = None
        self._max_period: float | None = None

    # --- Deep copy method ---
    def __deepcopy__(self, memo: dict) -> CoreWingCrossSectionMovement:
        """Creates a deep copy of this CoreWingCrossSectionMovement.

        See WingCrossSectionMovement for full documentation.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new instance with copied attributes.
        """
        # Create a new instance without calling __init__ to avoid redundant
        # validation. Use type(self) so subclasses get the correct type.
        new_movement = object.__new__(type(self))

        # Store this instance in memo to handle potential circular references.
        memo[id(self)] = new_movement

        # Deep copy the base WingCrossSection to ensure independence
        # (immutable).
        new_movement._base_wing_cross_section = copy.deepcopy(
            self._base_wing_cross_section, memo
        )

        # Copy numpy arrays and make them read only.
        new_movement._ampLp_Wcsp_Lpp = self._ampLp_Wcsp_Lpp.copy()
        new_movement._ampLp_Wcsp_Lpp.flags.writeable = False

        new_movement._periodLp_Wcsp_Lpp = self._periodLp_Wcsp_Lpp.copy()
        new_movement._periodLp_Wcsp_Lpp.flags.writeable = False

        new_movement._phaseLp_Wcsp_Lpp = self._phaseLp_Wcsp_Lpp.copy()
        new_movement._phaseLp_Wcsp_Lpp.flags.writeable = False

        new_movement._ampAngles_Wcsp_to_Wcs_ixyz = (
            self._ampAngles_Wcsp_to_Wcs_ixyz.copy()
        )
        new_movement._ampAngles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        new_movement._periodAngles_Wcsp_to_Wcs_ixyz = (
            self._periodAngles_Wcsp_to_Wcs_ixyz.copy()
        )
        new_movement._periodAngles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        new_movement._phaseAngles_Wcsp_to_Wcs_ixyz = (
            self._phaseAngles_Wcsp_to_Wcs_ixyz.copy()
        )
        new_movement._phaseAngles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        # Copy tuples directly (they are immutable).
        new_movement._spacingLp_Wcsp_Lpp = self._spacingLp_Wcsp_Lpp
        new_movement._spacingAngles_Wcsp_to_Wcs_ixyz = (
            self._spacingAngles_Wcsp_to_Wcs_ixyz
        )

        # Initialize cache variables to None (caches will be recomputed on
        # access).
        new_movement._all_periods = None
        new_movement._max_period = None

        return new_movement

    # --- Immutable: read only properties ---
    @property
    def base_wing_cross_section(
        self,
    ) -> geometry.wing_cross_section.WingCrossSection:
        return self._base_wing_cross_section

    @property
    def ampLp_Wcsp_Lpp(self) -> np.ndarray:
        return self._ampLp_Wcsp_Lpp

    @property
    def periodLp_Wcsp_Lpp(self) -> np.ndarray:
        return self._periodLp_Wcsp_Lpp

    @property
    def spacingLp_Wcsp_Lpp(
        self,
    ) -> tuple[str | Callable[[float], float], ...]:
        return self._spacingLp_Wcsp_Lpp

    @property
    def phaseLp_Wcsp_Lpp(self) -> np.ndarray:
        return self._phaseLp_Wcsp_Lpp

    @property
    def ampAngles_Wcsp_to_Wcs_ixyz(self) -> np.ndarray:
        return self._ampAngles_Wcsp_to_Wcs_ixyz

    @property
    def periodAngles_Wcsp_to_Wcs_ixyz(self) -> np.ndarray:
        return self._periodAngles_Wcsp_to_Wcs_ixyz

    @property
    def spacingAngles_Wcsp_to_Wcs_ixyz(
        self,
    ) -> tuple[str | Callable[[float], float], ...]:
        return self._spacingAngles_Wcsp_to_Wcs_ixyz

    @property
    def phaseAngles_Wcsp_to_Wcs_ixyz(self) -> np.ndarray:
        return self._phaseAngles_Wcsp_to_Wcs_ixyz

    # --- Immutable derived: manual lazy caching ---
    @property
    def all_periods(self) -> tuple[float, ...]:
        """All unique non zero periods of this WingCrossSection's motion.

        :return: A tuple of all unique non zero periods in seconds. If the motion is
            static, this will be an empty tuple.
        """
        if self._all_periods is None:
            periods = []

            # Collect all periods from positional motion.
            for period in self._periodLp_Wcsp_Lpp:
                if period > 0.0:
                    periods.append(float(period))

            # Collect all periods from angular motion.
            for period in self._periodAngles_Wcsp_to_Wcs_ixyz:
                if period > 0.0:
                    periods.append(float(period))

            self._all_periods = tuple(periods)
        return self._all_periods

    @property
    def max_period(self) -> float:
        """The longest period of this WingCrossSection's motion.

        :return: The longest period in seconds. If the motion is static, this will be
            0.0.
        """
        if self._max_period is None:
            self._max_period = float(
                max(
                    np.max(self._periodLp_Wcsp_Lpp),
                    np.max(self._periodAngles_Wcsp_to_Wcs_ixyz),
                )
            )
        return self._max_period

    # --- Other methods ---
    def generate_wing_cross_section_at_time_step(
        self, step: int, delta_time: float | int
    ) -> geometry.wing_cross_section.WingCrossSection:
        """Creates the WingCrossSection at a single time step.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :return: The WingCrossSection at this time step.
        """
        time = step * delta_time

        # Evaluate the oscillating value for each dimension of Lp_Wcsp_Lpp.
        thisLp_Wcsp_Lpp = np.zeros(3, dtype=float)
        for dim in range(3):
            this_spacing = self._spacingLp_Wcsp_Lpp[dim]
            this_amp = self._ampLp_Wcsp_Lpp[dim]
            this_period = self._periodLp_Wcsp_Lpp[dim]
            this_phase = self._phaseLp_Wcsp_Lpp[dim]
            this_base = self._base_wing_cross_section.Lp_Wcsp_Lpp[dim]

            if this_spacing == "sine":
                thisLp_Wcsp_Lpp[dim] = _oscillation.oscillating_sin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif this_spacing == "uniform":
                thisLp_Wcsp_Lpp[dim] = _oscillation.oscillating_lin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif callable(this_spacing):
                thisLp_Wcsp_Lpp[dim] = _oscillation.oscillating_custom_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                    custom_function=this_spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {this_spacing}")

        # Evaluate the oscillating value for each dimension of
        # angles_Wcsp_to_Wcs_ixyz.
        theseAngles_Wcsp_to_Wcs_ixyz = np.zeros(3, dtype=float)
        for dim in range(3):
            this_spacing = self._spacingAngles_Wcsp_to_Wcs_ixyz[dim]
            this_amp = self._ampAngles_Wcsp_to_Wcs_ixyz[dim]
            this_period = self._periodAngles_Wcsp_to_Wcs_ixyz[dim]
            this_phase = self._phaseAngles_Wcsp_to_Wcs_ixyz[dim]
            this_base = self._base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[dim]

            if this_spacing == "sine":
                theseAngles_Wcsp_to_Wcs_ixyz[dim] = (
                    _oscillation.oscillating_sin_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                    )
                )
            elif this_spacing == "uniform":
                theseAngles_Wcsp_to_Wcs_ixyz[dim] = (
                    _oscillation.oscillating_lin_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                    )
                )
            elif callable(this_spacing):
                theseAngles_Wcsp_to_Wcs_ixyz[dim] = (
                    _oscillation.oscillating_custom_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                        custom_function=this_spacing,
                    )
                )
            else:
                raise ValueError(f"Invalid spacing value: {this_spacing}")

        return geometry.wing_cross_section.WingCrossSection(
            airfoil=self._base_wing_cross_section.airfoil,
            num_spanwise_panels=self._base_wing_cross_section.num_spanwise_panels,
            chord=self._base_wing_cross_section.chord,
            Lp_Wcsp_Lpp=thisLp_Wcsp_Lpp,
            angles_Wcsp_to_Wcs_ixyz=theseAngles_Wcsp_to_Wcs_ixyz,
            control_surface_symmetry_type=(
                self._base_wing_cross_section.control_surface_symmetry_type
            ),
            control_surface_hinge_point=(
                self._base_wing_cross_section.control_surface_hinge_point
            ),
            control_surface_deflection=(
                self._base_wing_cross_section.control_surface_deflection
            ),
            spanwise_spacing=self._base_wing_cross_section.spanwise_spacing,
        )

    def generate_wing_cross_sections(
        self,
        num_steps: int,
        delta_time: float | int,
    ) -> list[geometry.wing_cross_section.WingCrossSection]:
        """Creates the WingCrossSection at each time step, and returns them in a list.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (int or float), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of WingCrossSections associated with this
            CoreWingCrossSectionMovement.
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

        wing_cross_sections = []
        for step in range(num_steps):
            wing_cross_sections.append(
                self.generate_wing_cross_section_at_time_step(step, delta_time)
            )

        return wing_cross_sections


class CoreWingMovement:
    """A core class used to contain the shared foundation of WingMovement and its
    feature variant siblings.

    See WingMovement for full documentation of the shared interface.

    CoreWingMovement holds the base Wing, WingCrossSectionMovements, and oscillation
    parameters, and provides generate_wing_at_time_step() for creating Wings one step at
    a time.
    """

    __slots__ = (
        "_base_wing",
        "_wing_cross_section_movements",
        "_ampLer_Gs_Cgs",
        "_periodLer_Gs_Cgs",
        "_spacingLer_Gs_Cgs",
        "_phaseLer_Gs_Cgs",
        "_ampAngles_Gs_to_Wn_ixyz",
        "_periodAngles_Gs_to_Wn_ixyz",
        "_spacingAngles_Gs_to_Wn_ixyz",
        "_phaseAngles_Gs_to_Wn_ixyz",
        "_rotationPointOffset_Gs_Ler",
        "_all_periods",
        "_max_period",
    )

    def __init__(
        self,
        base_wing: geometry.wing.Wing,
        wing_cross_section_movements: Sequence[CoreWingCrossSectionMovement],
        ampLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        periodAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        spacingAngles_Gs_to_Wn_ixyz: (
            np.ndarray | Sequence[str | Callable[[float], float]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        rotationPointOffset_Gs_Ler: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
    ) -> None:
        """The initialization method.

        See WingMovement's initialization method for full parameter descriptions.

        :param base_wing: The base Wing.
        :param wing_cross_section_movements: The CoreWingCrossSectionMovements for each
            WingCrossSection.
        :param ampLer_Gs_Cgs: The amplitudes of Ler_Gs_Cgs oscillation in meters.
        :param periodLer_Gs_Cgs: The periods of Ler_Gs_Cgs oscillation in seconds.
        :param spacingLer_Gs_Cgs: The spacing types for Ler_Gs_Cgs oscillation.
        :param phaseLer_Gs_Cgs: The phase offsets of Ler_Gs_Cgs oscillation in degrees.
        :param ampAngles_Gs_to_Wn_ixyz: The amplitudes of angles_Gs_to_Wn_ixyz
            oscillation in degrees.
        :param periodAngles_Gs_to_Wn_ixyz: The periods of angles_Gs_to_Wn_ixyz
            oscillation in seconds.
        :param spacingAngles_Gs_to_Wn_ixyz: The spacing types for angles_Gs_to_Wn_ixyz
            oscillation.
        :param phaseAngles_Gs_to_Wn_ixyz: The phase offsets of angles_Gs_to_Wn_ixyz
            oscillation in degrees.
        :param rotationPointOffset_Gs_Ler: The position of the rotation point for
            angular motion in meters.
        :return: None
        """
        # Validate and store immutable attributes. Set those that are numpy
        # arrays to be read only.
        if not isinstance(base_wing, geometry.wing.Wing):
            raise TypeError("base_wing must be a Wing.")
        self._base_wing = base_wing

        if not isinstance(wing_cross_section_movements, list):
            raise TypeError("wing_cross_section_movements must be a list.")
        if len(wing_cross_section_movements) != len(
            self._base_wing.wing_cross_sections
        ):
            raise ValueError(
                "wing_cross_section_movements must have the same length "
                "as base_wing.wing_cross_sections."
            )
        for wing_cross_section_movement in wing_cross_section_movements:
            if not isinstance(
                wing_cross_section_movement,
                CoreWingCrossSectionMovement,
            ):
                raise TypeError(
                    "Every element in wing_cross_section_movements must "
                    "be a CoreWingCrossSectionMovement."
                )
        # Store as tuple to prevent external mutation.
        self._wing_cross_section_movements: tuple[CoreWingCrossSectionMovement, ...] = (
            tuple(wing_cross_section_movements)
        )

        ampLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            ampLer_Gs_Cgs, "ampLer_Gs_Cgs"
        )
        if not np.all(ampLer_Gs_Cgs >= 0.0):
            raise ValueError("All elements in ampLer_Gs_Cgs must be non negative.")
        self._ampLer_Gs_Cgs = ampLer_Gs_Cgs
        self._ampLer_Gs_Cgs.flags.writeable = False

        periodLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            periodLer_Gs_Cgs, "periodLer_Gs_Cgs"
        )
        if not np.all(periodLer_Gs_Cgs >= 0.0):
            raise ValueError("All elements in periodLer_Gs_Cgs must be non negative.")
        for period_index, period in enumerate(periodLer_Gs_Cgs):
            amp = self._ampLer_Gs_Cgs[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampLer_Gs_Cgs is 0.0, the "
                    "corresponding element in periodLer_Gs_Cgs must be "
                    "also be 0.0."
                )
        self._periodLer_Gs_Cgs = periodLer_Gs_Cgs
        self._periodLer_Gs_Cgs.flags.writeable = False

        # Store as tuple to prevent external mutation.
        self._spacingLer_Gs_Cgs = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingLer_Gs_Cgs, "spacingLer_Gs_Cgs"
            )
        )

        phaseLer_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            phaseLer_Gs_Cgs, "phaseLer_Gs_Cgs"
        )
        if not (np.all(phaseLer_Gs_Cgs > -180.0) and np.all(phaseLer_Gs_Cgs <= 180.0)):
            raise ValueError(
                "All elements in phaseLer_Gs_Cgs must be in the range "
                "(-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseLer_Gs_Cgs):
            amp = self._ampLer_Gs_Cgs[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampLer_Gs_Cgs is 0.0, the "
                    "corresponding element in phaseLer_Gs_Cgs must be "
                    "also be 0.0."
                )
        self._phaseLer_Gs_Cgs = phaseLer_Gs_Cgs
        self._phaseLer_Gs_Cgs.flags.writeable = False

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
                "All elements in ampAngles_Gs_to_Wn_ixyz must be in "
                "the range [0.0, 180.0]."
            )
        self._ampAngles_Gs_to_Wn_ixyz = ampAngles_Gs_to_Wn_ixyz
        self._ampAngles_Gs_to_Wn_ixyz.flags.writeable = False

        periodAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                periodAngles_Gs_to_Wn_ixyz,
                "periodAngles_Gs_to_Wn_ixyz",
            )
        )
        if not np.all(periodAngles_Gs_to_Wn_ixyz >= 0.0):
            raise ValueError(
                "All elements in periodAngles_Gs_to_Wn_ixyz must be " "non negative."
            )
        for period_index, period in enumerate(periodAngles_Gs_to_Wn_ixyz):
            amp = self._ampAngles_Gs_to_Wn_ixyz[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampAngles_Gs_to_Wn_ixyz is 0.0, "
                    "the corresponding element in "
                    "periodAngles_Gs_to_Wn_ixyz must be also be 0.0."
                )
        self._periodAngles_Gs_to_Wn_ixyz = periodAngles_Gs_to_Wn_ixyz
        self._periodAngles_Gs_to_Wn_ixyz.flags.writeable = False

        # Store as tuple to prevent external mutation.
        self._spacingAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingAngles_Gs_to_Wn_ixyz,
                "spacingAngles_Gs_to_Wn_ixyz",
            )
        )

        phaseAngles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                phaseAngles_Gs_to_Wn_ixyz,
                "phaseAngles_Gs_to_Wn_ixyz",
            )
        )
        if not (
            np.all(phaseAngles_Gs_to_Wn_ixyz > -180.0)
            and np.all(phaseAngles_Gs_to_Wn_ixyz <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseAngles_Gs_to_Wn_ixyz must be in "
                "the range (-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseAngles_Gs_to_Wn_ixyz):
            amp = self._ampAngles_Gs_to_Wn_ixyz[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampAngles_Gs_to_Wn_ixyz is 0.0, "
                    "the corresponding element in "
                    "phaseAngles_Gs_to_Wn_ixyz must be also be 0.0."
                )
        self._phaseAngles_Gs_to_Wn_ixyz = phaseAngles_Gs_to_Wn_ixyz
        self._phaseAngles_Gs_to_Wn_ixyz.flags.writeable = False

        rotationPointOffset_Gs_Ler = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                rotationPointOffset_Gs_Ler,
                "rotationPointOffset_Gs_Ler",
            )
        )
        self._rotationPointOffset_Gs_Ler = rotationPointOffset_Gs_Ler
        self._rotationPointOffset_Gs_Ler.flags.writeable = False

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._all_periods: tuple[float, ...] | None = None
        self._max_period: float | None = None

    # --- Deep copy method ---
    def __deepcopy__(self, memo: dict) -> CoreWingMovement:
        """Creates a deep copy of this CoreWingMovement.

        See WingMovement for full documentation.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new instance with copied attributes.
        """
        # Create a new instance without calling __init__ to avoid redundant
        # validation. Use type(self) so subclasses get the correct type.
        new_movement = object.__new__(type(self))

        # Store this instance in memo to handle potential circular references.
        memo[id(self)] = new_movement

        # Deep copy the base Wing to ensure independence (immutable).
        new_movement._base_wing = copy.deepcopy(self._base_wing, memo)

        # Deep copy the WingCrossSectionMovements and store as tuple.
        new_movement._wing_cross_section_movements = tuple(
            copy.deepcopy(wing_cross_section_movement, memo)
            for wing_cross_section_movement in self._wing_cross_section_movements
        )

        # Copy numpy arrays and make them read only.
        new_movement._ampLer_Gs_Cgs = self._ampLer_Gs_Cgs.copy()
        new_movement._ampLer_Gs_Cgs.flags.writeable = False

        new_movement._periodLer_Gs_Cgs = self._periodLer_Gs_Cgs.copy()
        new_movement._periodLer_Gs_Cgs.flags.writeable = False

        new_movement._phaseLer_Gs_Cgs = self._phaseLer_Gs_Cgs.copy()
        new_movement._phaseLer_Gs_Cgs.flags.writeable = False

        new_movement._ampAngles_Gs_to_Wn_ixyz = self._ampAngles_Gs_to_Wn_ixyz.copy()
        new_movement._ampAngles_Gs_to_Wn_ixyz.flags.writeable = False

        new_movement._periodAngles_Gs_to_Wn_ixyz = (
            self._periodAngles_Gs_to_Wn_ixyz.copy()
        )
        new_movement._periodAngles_Gs_to_Wn_ixyz.flags.writeable = False

        new_movement._phaseAngles_Gs_to_Wn_ixyz = self._phaseAngles_Gs_to_Wn_ixyz.copy()
        new_movement._phaseAngles_Gs_to_Wn_ixyz.flags.writeable = False

        new_movement._rotationPointOffset_Gs_Ler = (
            self._rotationPointOffset_Gs_Ler.copy()
        )
        new_movement._rotationPointOffset_Gs_Ler.flags.writeable = False

        # Copy tuples directly (they are immutable).
        new_movement._spacingLer_Gs_Cgs = self._spacingLer_Gs_Cgs
        new_movement._spacingAngles_Gs_to_Wn_ixyz = self._spacingAngles_Gs_to_Wn_ixyz

        # Initialize cache variables to None (caches will be recomputed on
        # access).
        new_movement._all_periods = None
        new_movement._max_period = None

        return new_movement

    # --- Immutable: read only properties ---
    @property
    def base_wing(self) -> geometry.wing.Wing:
        return self._base_wing

    @property
    def wing_cross_section_movements(
        self,
    ) -> tuple[CoreWingCrossSectionMovement, ...]:
        return self._wing_cross_section_movements

    @property
    def ampLer_Gs_Cgs(self) -> np.ndarray:
        return self._ampLer_Gs_Cgs

    @property
    def periodLer_Gs_Cgs(self) -> np.ndarray:
        return self._periodLer_Gs_Cgs

    @property
    def spacingLer_Gs_Cgs(
        self,
    ) -> tuple[str | Callable[[float], float], ...]:
        return self._spacingLer_Gs_Cgs

    @property
    def phaseLer_Gs_Cgs(self) -> np.ndarray:
        return self._phaseLer_Gs_Cgs

    @property
    def ampAngles_Gs_to_Wn_ixyz(self) -> np.ndarray:
        return self._ampAngles_Gs_to_Wn_ixyz

    @property
    def periodAngles_Gs_to_Wn_ixyz(self) -> np.ndarray:
        return self._periodAngles_Gs_to_Wn_ixyz

    @property
    def spacingAngles_Gs_to_Wn_ixyz(
        self,
    ) -> tuple[str | Callable[[float], float], ...]:
        return self._spacingAngles_Gs_to_Wn_ixyz

    @property
    def phaseAngles_Gs_to_Wn_ixyz(self) -> np.ndarray:
        return self._phaseAngles_Gs_to_Wn_ixyz

    @property
    def rotationPointOffset_Gs_Ler(self) -> np.ndarray:
        return self._rotationPointOffset_Gs_Ler

    # --- Immutable derived: manual lazy caching ---
    @property
    def all_periods(self) -> tuple[float, ...]:
        """All unique non zero periods of motion from this Wing and each
        WingCrossSection's movement class.

        :return: A tuple of all unique non zero periods in seconds. If all motion is
            static, this will be an empty tuple.
        """
        if self._all_periods is None:
            periods: list[float] = []

            # Collect all periods from WingCrossSectionMovements.
            for wing_cross_section_movement in self._wing_cross_section_movements:
                periods.extend(wing_cross_section_movement.all_periods)

            # Collect all periods from CoreWingMovement's own motion.
            for period in self._periodLer_Gs_Cgs:
                if period > 0.0:
                    periods.append(float(period))
            for period in self._periodAngles_Gs_to_Wn_ixyz:
                if period > 0.0:
                    periods.append(float(period))

            self._all_periods = tuple(periods)
        return self._all_periods

    @property
    def max_period(self) -> float:
        """The longest period of motion across this Wing and each WingCrossSection's
        movement class.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        if self._max_period is None:
            wing_cross_section_movement_max_periods = []
            for wing_cross_section_movement in self._wing_cross_section_movements:
                wing_cross_section_movement_max_periods.append(
                    wing_cross_section_movement.max_period
                )
            max_wing_cross_section_movement_period = max(
                wing_cross_section_movement_max_periods
            )

            self._max_period = float(
                max(
                    max_wing_cross_section_movement_period,
                    np.max(self._periodLer_Gs_Cgs),
                    np.max(self._periodAngles_Gs_to_Wn_ixyz),
                )
            )
        return self._max_period

    # --- Other methods ---
    def generate_wing_at_time_step(
        self, step: int, delta_time: float | int
    ) -> geometry.wing.Wing:
        """Creates the Wing at a single time step.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :return: The Wing at this time step.
        """
        time = step * delta_time

        # Evaluate the oscillating value for each dimension of Ler_Gs_Cgs.
        thisLer_Gs_Cgs = np.zeros(3, dtype=float)
        for dim in range(3):
            this_spacing = self._spacingLer_Gs_Cgs[dim]
            this_amp = self._ampLer_Gs_Cgs[dim]
            this_period = self._periodLer_Gs_Cgs[dim]
            this_phase = self._phaseLer_Gs_Cgs[dim]
            this_base = self._base_wing.Ler_Gs_Cgs[dim]

            if this_spacing == "sine":
                thisLer_Gs_Cgs[dim] = _oscillation.oscillating_sin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif this_spacing == "uniform":
                thisLer_Gs_Cgs[dim] = _oscillation.oscillating_lin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif callable(this_spacing):
                thisLer_Gs_Cgs[dim] = _oscillation.oscillating_custom_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                    custom_function=this_spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {this_spacing}")

        # Evaluate the oscillating value for each dimension of
        # angles_Gs_to_Wn_ixyz.
        theseAngles_Gs_to_Wn_ixyz = np.zeros(3, dtype=float)
        for dim in range(3):
            this_spacing = self._spacingAngles_Gs_to_Wn_ixyz[dim]
            this_amp = self._ampAngles_Gs_to_Wn_ixyz[dim]
            this_period = self._periodAngles_Gs_to_Wn_ixyz[dim]
            this_phase = self._phaseAngles_Gs_to_Wn_ixyz[dim]
            this_base = self._base_wing.angles_Gs_to_Wn_ixyz[dim]

            if this_spacing == "sine":
                theseAngles_Gs_to_Wn_ixyz[dim] = _oscillation.oscillating_sin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif this_spacing == "uniform":
                theseAngles_Gs_to_Wn_ixyz[dim] = _oscillation.oscillating_lin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif callable(this_spacing):
                theseAngles_Gs_to_Wn_ixyz[dim] = (
                    _oscillation.oscillating_custom_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                        custom_function=this_spacing,
                    )
                )
            else:
                raise ValueError(f"Invalid spacing value: {this_spacing}")

        # Generate the WingCrossSections for this time step.
        these_wing_cross_sections = []
        for wing_cross_section_movement in self._wing_cross_section_movements:
            these_wing_cross_sections.append(
                wing_cross_section_movement.generate_wing_cross_section_at_time_step(
                    step, delta_time
                )
            )

        # If there is a non zero rotation point offset, adjust the position
        # to account for rotation about the offset point instead of the
        # leading edge root.
        if not np.allclose(self._rotationPointOffset_Gs_Ler, np.zeros(3, dtype=float)):
            # Get the active rotation matrix for this step's angles.
            rot_T_act = _transformations.generate_rot_T(
                theseAngles_Gs_to_Wn_ixyz,
                passive=False,
                intrinsic=True,
                order="xyz",
            )
            rot_R_act = rot_T_act[:3, :3]

            # Compute the position adjustment due to the offset rotation
            # point.
            offsetRotationPointAdjustment_Gs = (
                _transformations.compute_offset_rotation_adjustment(
                    rotation_matrix=rot_R_act,
                    offset=self._rotationPointOffset_Gs_Ler,
                )
            )

            # Apply the position adjustment to the leading edge root.
            thisLer_Gs_Cgs = thisLer_Gs_Cgs + offsetRotationPointAdjustment_Gs

        return geometry.wing.Wing(
            wing_cross_sections=these_wing_cross_sections,
            name=self._base_wing.name,
            Ler_Gs_Cgs=thisLer_Gs_Cgs,
            angles_Gs_to_Wn_ixyz=theseAngles_Gs_to_Wn_ixyz,
            symmetric=self._base_wing.symmetric,
            mirror_only=self._base_wing.mirror_only,
            symmetryNormal_G=self._base_wing.symmetryNormal_G,
            symmetryPoint_G_Cg=self._base_wing.symmetryPoint_G_Cg,
            num_chordwise_panels=self._base_wing.num_chordwise_panels,
            chordwise_spacing=self._base_wing.chordwise_spacing,
        )

    def generate_wings(
        self, num_steps: int, delta_time: float | int
    ) -> list[geometry.wing.Wing]:
        """Creates the Wing at each time step, and returns them in a list.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (int or float), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of Wings associated with this CoreWingMovement.
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

        wings = []
        for step in range(num_steps):
            wings.append(self.generate_wing_at_time_step(step, delta_time))

        return wings


class CoreAirplaneMovement:
    """A core class used to contain the shared foundation of AirplaneMovement and its
    feature variant siblings.

    See AirplaneMovement for full documentation of the shared interface.

    CoreAirplaneMovement holds the base Airplane, WingMovements, and oscillation
    parameters, and provides generate_airplane_at_time_step() for creating Airplanes one
    step at a time.
    """

    __slots__ = (
        "_base_airplane",
        "_wing_movements",
        "_ampCg_GP1_CgP1",
        "_periodCg_GP1_CgP1",
        "_spacingCg_GP1_CgP1",
        "_phaseCg_GP1_CgP1",
        "_all_periods",
        "_max_period",
    )

    def __init__(
        self,
        base_airplane: geometry.airplane.Airplane,
        wing_movements: Sequence[CoreWingMovement],
        ampCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        spacingCg_GP1_CgP1: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
    ) -> None:
        """The initialization method.

        See AirplaneMovement's initialization method for full parameter descriptions.

        :param base_airplane: The base Airplane.
        :param wing_movements: The CoreWingMovements for each Wing.
        :param ampCg_GP1_CgP1: The amplitudes of Cg_GP1_CgP1 oscillation in meters.
        :param periodCg_GP1_CgP1: The periods of Cg_GP1_CgP1 oscillation in seconds.
        :param spacingCg_GP1_CgP1: The spacing types for Cg_GP1_CgP1 oscillation.
        :param phaseCg_GP1_CgP1: The phase offsets of Cg_GP1_CgP1 oscillation in
            degrees.
        :return: None
        """
        # Validate and store immutable attributes. Set those that are numpy
        # arrays to be read only.
        if not isinstance(base_airplane, geometry.airplane.Airplane):
            raise TypeError("base_airplane must be an Airplane.")
        self._base_airplane = base_airplane

        if not isinstance(wing_movements, list):
            raise TypeError("wing_movements must be a list.")
        if len(wing_movements) != len(self._base_airplane.wings):
            raise ValueError(
                "wing_movements must have the same length as " "base_airplane.wings."
            )
        for wing_movement in wing_movements:
            if not isinstance(wing_movement, CoreWingMovement):
                raise TypeError(
                    "Every element in wing_movements must be a " "CoreWingMovement."
                )
        # Store as tuple to prevent external mutation.
        self._wing_movements: tuple[CoreWingMovement, ...] = tuple(wing_movements)

        ampCg_GP1_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            ampCg_GP1_CgP1, "ampCg_GP1_CgP1"
        )
        if not np.all(ampCg_GP1_CgP1 >= 0.0):
            raise ValueError("All elements in ampCg_GP1_CgP1 must be non negative.")
        self._ampCg_GP1_CgP1 = ampCg_GP1_CgP1
        self._ampCg_GP1_CgP1.flags.writeable = False

        periodCg_GP1_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            periodCg_GP1_CgP1, "periodCg_GP1_CgP1"
        )
        if not np.all(periodCg_GP1_CgP1 >= 0.0):
            raise ValueError("All elements in periodCg_GP1_CgP1 must be non negative.")
        for period_index, period in enumerate(periodCg_GP1_CgP1):
            amp = self._ampCg_GP1_CgP1[period_index]
            if amp == 0 and period != 0:
                raise ValueError(
                    "If an element in ampCg_GP1_CgP1 is 0.0, the "
                    "corresponding element in periodCg_GP1_CgP1 must "
                    "be also be 0.0."
                )
        self._periodCg_GP1_CgP1 = periodCg_GP1_CgP1
        self._periodCg_GP1_CgP1.flags.writeable = False

        # Store as tuple to prevent external mutation.
        self._spacingCg_GP1_CgP1 = (
            _parameter_validation.threeD_spacing_vectorLike_return_tuple(
                spacingCg_GP1_CgP1, "spacingCg_GP1_CgP1"
            )
        )

        phaseCg_GP1_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            phaseCg_GP1_CgP1, "phaseCg_GP1_CgP1"
        )
        if not (
            np.all(phaseCg_GP1_CgP1 > -180.0) and np.all(phaseCg_GP1_CgP1 <= 180.0)
        ):
            raise ValueError(
                "All elements in phaseCg_GP1_CgP1 must be in the range "
                "(-180.0, 180.0]."
            )
        for phase_index, phase in enumerate(phaseCg_GP1_CgP1):
            amp = self._ampCg_GP1_CgP1[phase_index]
            if amp == 0 and phase != 0:
                raise ValueError(
                    "If an element in ampCg_GP1_CgP1 is 0.0, the "
                    "corresponding element in phaseCg_GP1_CgP1 must "
                    "be also be 0.0."
                )
        self._phaseCg_GP1_CgP1 = phaseCg_GP1_CgP1
        self._phaseCg_GP1_CgP1.flags.writeable = False

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._all_periods: tuple[float, ...] | None = None
        self._max_period: float | None = None

    # --- Deep copy method ---
    def __deepcopy__(self, memo: dict) -> CoreAirplaneMovement:
        """Creates a deep copy of this CoreAirplaneMovement.

        See AirplaneMovement for full documentation.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new instance with copied attributes.
        """
        # Create a new instance without calling __init__ to avoid redundant
        # validation. Use type(self) so subclasses get the correct type.
        new_movement = object.__new__(type(self))

        # Store this instance in memo to handle potential circular references.
        memo[id(self)] = new_movement

        # Deep copy the base Airplane to ensure independence (immutable).
        new_movement._base_airplane = copy.deepcopy(self._base_airplane, memo)

        # Deep copy WingMovements and store as tuple.
        new_movement._wing_movements = tuple(
            copy.deepcopy(wing_movement, memo) for wing_movement in self._wing_movements
        )

        # Copy numpy arrays and make them read only.
        new_movement._ampCg_GP1_CgP1 = self._ampCg_GP1_CgP1.copy()
        new_movement._ampCg_GP1_CgP1.flags.writeable = False

        new_movement._periodCg_GP1_CgP1 = self._periodCg_GP1_CgP1.copy()
        new_movement._periodCg_GP1_CgP1.flags.writeable = False

        new_movement._phaseCg_GP1_CgP1 = self._phaseCg_GP1_CgP1.copy()
        new_movement._phaseCg_GP1_CgP1.flags.writeable = False

        # Copy tuple directly (it is immutable).
        new_movement._spacingCg_GP1_CgP1 = self._spacingCg_GP1_CgP1

        # Initialize cache variables to None (caches will be recomputed on
        # access).
        new_movement._all_periods = None
        new_movement._max_period = None

        return new_movement

    # --- Immutable: read only properties ---
    @property
    def base_airplane(self) -> geometry.airplane.Airplane:
        return self._base_airplane

    @property
    def wing_movements(self) -> tuple[CoreWingMovement, ...]:
        return self._wing_movements

    @property
    def ampCg_GP1_CgP1(self) -> np.ndarray:
        return self._ampCg_GP1_CgP1

    @property
    def periodCg_GP1_CgP1(self) -> np.ndarray:
        return self._periodCg_GP1_CgP1

    @property
    def spacingCg_GP1_CgP1(
        self,
    ) -> tuple[str | Callable[[float], float], ...]:
        return self._spacingCg_GP1_CgP1

    @property
    def phaseCg_GP1_CgP1(self) -> np.ndarray:
        return self._phaseCg_GP1_CgP1

    # --- Immutable derived: manual lazy caching ---
    @property
    def all_periods(self) -> tuple[float, ...]:
        """All unique non zero periods of motion from this Airplane, each Wing's
        movement class, and each WingCrossSection's movement class.

        :return: A tuple of all unique non zero periods in seconds. If all motion is
            static, this will be an empty tuple.
        """
        if self._all_periods is None:
            periods: list[float] = []

            # Collect all periods from WingMovement(s).
            for wing_movement in self._wing_movements:
                periods.extend(wing_movement.all_periods)

            # Collect all periods from CoreAirplaneMovement's own motion.
            for period in self._periodCg_GP1_CgP1:
                if period > 0.0:
                    periods.append(float(period))

            self._all_periods = tuple(periods)
        return self._all_periods

    @property
    def max_period(self) -> float:
        """The longest period of motion across this Airplane, each Wing's movement
        class, and each WingCrossSection's movement class.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        if self._max_period is None:
            wing_movement_max_periods = []
            for wing_movement in self._wing_movements:
                wing_movement_max_periods.append(wing_movement.max_period)
            max_wing_movement_period = max(wing_movement_max_periods)

            self._max_period = float(
                max(
                    max_wing_movement_period,
                    np.max(self._periodCg_GP1_CgP1),
                )
            )
        return self._max_period

    # --- Other methods ---
    def generate_airplane_at_time_step(
        self, step: int, delta_time: float | int
    ) -> geometry.airplane.Airplane:
        """Creates the Airplane at a single time step.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :return: The Airplane at this time step.
        """
        time = step * delta_time

        # Evaluate the oscillating value for each dimension of Cg_GP1_CgP1.
        thisCg_GP1_CgP1 = np.zeros(3, dtype=float)
        for dim in range(3):
            this_spacing = self._spacingCg_GP1_CgP1[dim]
            this_amp = self._ampCg_GP1_CgP1[dim]
            this_period = self._periodCg_GP1_CgP1[dim]
            this_phase = self._phaseCg_GP1_CgP1[dim]
            this_base = self._base_airplane.Cg_GP1_CgP1[dim]

            if this_spacing == "sine":
                thisCg_GP1_CgP1[dim] = _oscillation.oscillating_sin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif this_spacing == "uniform":
                thisCg_GP1_CgP1[dim] = _oscillation.oscillating_lin_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                )
            elif callable(this_spacing):
                thisCg_GP1_CgP1[dim] = _oscillation.oscillating_custom_at_time(
                    amp=this_amp,
                    period=this_period,
                    phase=this_phase,
                    base=this_base,
                    time=time,
                    custom_function=this_spacing,
                )
            else:
                raise ValueError(f"Invalid spacing value: {this_spacing}")

        # Generate the Wings for this time step.
        these_wings = []
        for wing_movement in self._wing_movements:
            these_wings.append(
                wing_movement.generate_wing_at_time_step(step, delta_time)
            )

        return geometry.airplane.Airplane(
            wings=these_wings,
            name=self._base_airplane.name,
            Cg_GP1_CgP1=thisCg_GP1_CgP1,
            weight=self._base_airplane.weight,
        )

    def generate_airplanes(
        self, num_steps: int, delta_time: float | int
    ) -> list[geometry.airplane.Airplane]:
        """Creates the Airplane at each time step, and returns them in a list.

        For static geometry (no periodic motion), this method optimizes performance by
        creating the first Airplane with full mesh generation, then using deepcopy for
        subsequent time steps. This avoids redundant mesh generation when the geometry
        is identical across all steps.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (float or int), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of Airplanes associated with this CoreAirplaneMovement.
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

        # Check if geometry is static (no periodic motion).
        is_static_geometry = self.max_period == 0.0

        if is_static_geometry:
            # Optimization for static geometry: create first Airplane with full mesh
            # generation, then deepcopy for subsequent steps.
            return self._generate_airplanes_static(num_steps, delta_time)
        else:
            # For variable geometry, use the standard approach.
            return self._generate_airplanes_variable(num_steps, delta_time)

    def _generate_airplanes_static(
        self, num_steps: int, delta_time: float
    ) -> list[geometry.airplane.Airplane]:
        """Generates Airplanes for static geometry using deepcopy optimization.

        Creates the first Airplane with full mesh generation, then uses deepcopy for
        subsequent time steps to avoid redundant mesh generation.

        :param num_steps: The number of time steps.
        :param delta_time: The time between each time step in seconds.
        :return: The list of Airplanes.
        """
        # Create the first Airplane (triggers full mesh generation).
        first_airplane = self.generate_airplane_at_time_step(0, delta_time)

        # Create list with first Airplane.
        airplanes = [first_airplane]

        # Create copies for remaining steps with different Cg_GP1_CgP1 positions.
        for step in range(1, num_steps):
            time = step * delta_time

            # Evaluate Cg_GP1_CgP1 at this step.
            thisCg_GP1_CgP1 = np.zeros(3, dtype=float)
            for dim in range(3):
                this_spacing = self._spacingCg_GP1_CgP1[dim]
                this_amp = self._ampCg_GP1_CgP1[dim]
                this_period = self._periodCg_GP1_CgP1[dim]
                this_phase = self._phaseCg_GP1_CgP1[dim]
                this_base = self._base_airplane.Cg_GP1_CgP1[dim]

                if this_spacing == "sine":
                    thisCg_GP1_CgP1[dim] = _oscillation.oscillating_sin_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                    )
                elif this_spacing == "uniform":
                    thisCg_GP1_CgP1[dim] = _oscillation.oscillating_lin_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                    )
                elif callable(this_spacing):
                    thisCg_GP1_CgP1[dim] = _oscillation.oscillating_custom_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                        custom_function=this_spacing,
                    )
                else:
                    raise ValueError(f"Invalid spacing value: {this_spacing}")

            copied_airplane = first_airplane.deep_copy_with_Cg_GP1_CgP1(thisCg_GP1_CgP1)
            airplanes.append(copied_airplane)

        return airplanes

    def _generate_airplanes_variable(
        self, num_steps: int, delta_time: float
    ) -> list[geometry.airplane.Airplane]:
        """Generates Airplanes for variable (periodic) geometry.

        Uses a conservative optimization approach that validates periodicity before
        applying deepcopy optimization. If validation fails, falls back to standard
        generation.

        :param num_steps: The number of time steps.
        :param delta_time: The time between each time step in seconds.
        :return: The list of Airplanes.
        """
        # Step 1: Calculate geometry LCM period.
        geometry_lcm = self._geometry_lcm_period()

        # Step 2: Pre-validation checks.
        # Check if period aligns cleanly with delta_time.
        steps_per_period_float = geometry_lcm / delta_time
        steps_per_period = int(round(steps_per_period_float))

        alignment_error = abs(steps_per_period_float - steps_per_period)
        if alignment_error > 1e-6:
            # Period doesn't align with time steps. Fall back to standard generation.
            return self._generate_airplanes_variable_standard(num_steps, delta_time)

        # Check if there's meaningful benefit.
        if num_steps <= steps_per_period:
            # No repetition occurs. No benefit from optimization.
            return self._generate_airplanes_variable_standard(num_steps, delta_time)

        # Step 3: Generate first period + one validation step.
        validation_num_steps = steps_per_period + 1

        # Create an empty 2D ndarray to hold Wings for validation steps.
        wings = np.empty(
            (len(self._wing_movements), validation_num_steps), dtype=object
        )

        # Iterate through the CoreWingMovements.
        for wing_movement_id, wing_movement in enumerate(self._wing_movements):
            this_wings_list_of_wings = np.array(
                wing_movement.generate_wings(
                    num_steps=validation_num_steps, delta_time=delta_time
                )
            )
            wings[wing_movement_id, :] = this_wings_list_of_wings

        # Step 4: Validate periodicity.
        # Compare step 0 geometry to step steps_per_period.
        if not self._geometry_matches(
            wings_step_a=wings[:, 0],
            wings_step_b=wings[:, steps_per_period],
            tolerance=1e-9,
        ):
            # Geometry doesn't repeat as expected. Fall back to standard generation.
            return self._generate_airplanes_variable_standard(num_steps, delta_time)

        # Step 5: Create Airplanes for first period.
        first_period_airplanes = []
        for step in range(steps_per_period):
            this_airplane = self.generate_airplane_at_time_step(step, delta_time)
            first_period_airplanes.append(this_airplane)

        # Step 6: Create copies for remaining steps with different Cg_GP1_CgP1.
        airplanes = list(first_period_airplanes)
        for step in range(steps_per_period, num_steps):
            source_step = step % steps_per_period

            # Evaluate Cg_GP1_CgP1 at this step.
            time = step * delta_time
            thisCg_GP1_CgP1 = np.zeros(3, dtype=float)
            for dim in range(3):
                this_spacing = self._spacingCg_GP1_CgP1[dim]
                this_amp = self._ampCg_GP1_CgP1[dim]
                this_period = self._periodCg_GP1_CgP1[dim]
                this_phase = self._phaseCg_GP1_CgP1[dim]
                this_base = self._base_airplane.Cg_GP1_CgP1[dim]

                if this_spacing == "sine":
                    thisCg_GP1_CgP1[dim] = _oscillation.oscillating_sin_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                    )
                elif this_spacing == "uniform":
                    thisCg_GP1_CgP1[dim] = _oscillation.oscillating_lin_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                    )
                elif callable(this_spacing):
                    thisCg_GP1_CgP1[dim] = _oscillation.oscillating_custom_at_time(
                        amp=this_amp,
                        period=this_period,
                        phase=this_phase,
                        base=this_base,
                        time=time,
                        custom_function=this_spacing,
                    )
                else:
                    raise ValueError(f"Invalid spacing value: {this_spacing}")

            copied_airplane = first_period_airplanes[
                source_step
            ].deep_copy_with_Cg_GP1_CgP1(thisCg_GP1_CgP1)
            airplanes.append(copied_airplane)

        return airplanes

    def _generate_airplanes_variable_standard(
        self, num_steps: int, delta_time: float
    ) -> list[geometry.airplane.Airplane]:
        """Generates Airplanes for variable geometry using standard approach.

        Creates new Airplanes for each time step without optimization.

        :param num_steps: The number of time steps.
        :param delta_time: The time between each time step in seconds.
        :return: The list of Airplanes.
        """
        airplanes = []
        for step in range(num_steps):
            airplanes.append(self.generate_airplane_at_time_step(step, delta_time))

        return airplanes

    def _geometry_lcm_period(self) -> float:
        """Calculates the LCM of all geometry related periods.

        Excludes OperatingPoint periods since those don't affect geometry. Uses the
        all_periods property which already collects only geometry related periods.

        :return: The LCM of all geometry related periods in seconds. Returns 0.0 if all
            geometry is static.
        """
        return lcm_multiple(list(self.all_periods))

    @staticmethod
    def _geometry_matches(
        wings_step_a: np.ndarray,
        wings_step_b: np.ndarray,
        tolerance: float = 1e-9,
    ) -> bool:
        """Compares two sets of Wings to verify their geometry matches within tolerance.

        Checks Wing position (Ler_Gs_Cgs), Wing angles (angles_Gs_to_Wn_ixyz), and Panel
        corner positions (Frpp_G_Cg, Flpp_G_Cg, Blpp_G_Cg, Brpp_G_Cg).

        :param wings_step_a: A (num_wings,) ndarray of Wings from the first time step.
        :param wings_step_b: A (num_wings,) ndarray of Wings from the second time step.
        :param tolerance: The tolerance for floating point comparison. The default is
            1e-9.
        :return: True if all geometry attributes match within tolerance, False
            otherwise.
        """
        if len(wings_step_a) != len(wings_step_b):
            return False

        for wing_a, wing_b in zip(wings_step_a, wings_step_b):
            # Check Wing position.
            if not np.allclose(
                wing_a.Ler_Gs_Cgs, wing_b.Ler_Gs_Cgs, atol=tolerance, rtol=0.0
            ):
                return False

            # Check Wing angles.
            if not np.allclose(
                wing_a.angles_Gs_to_Wn_ixyz,
                wing_b.angles_Gs_to_Wn_ixyz,
                atol=tolerance,
                rtol=0.0,
            ):
                return False

            # Check Panel corner positions if Wings are meshed.
            if wing_a.panels is not None and wing_b.panels is not None:
                if wing_a.panels.shape != wing_b.panels.shape:
                    return False

                for i in range(wing_a.panels.shape[0]):
                    for j in range(wing_a.panels.shape[1]):
                        panel_a = wing_a.panels[i, j]
                        panel_b = wing_b.panels[i, j]

                        # Check all four corner positions.
                        if not np.allclose(
                            panel_a.Frpp_G_Cg,
                            panel_b.Frpp_G_Cg,
                            atol=tolerance,
                            rtol=0.0,
                        ):
                            return False
                        if not np.allclose(
                            panel_a.Flpp_G_Cg,
                            panel_b.Flpp_G_Cg,
                            atol=tolerance,
                            rtol=0.0,
                        ):
                            return False
                        if not np.allclose(
                            panel_a.Blpp_G_Cg,
                            panel_b.Blpp_G_Cg,
                            atol=tolerance,
                            rtol=0.0,
                        ):
                            return False
                        if not np.allclose(
                            panel_a.Brpp_G_Cg,
                            panel_b.Brpp_G_Cg,
                            atol=tolerance,
                            rtol=0.0,
                        ):
                            return False

        return True


class CoreMovement:
    """A core class used to contain the shared foundation of Movement and its feature
    variant siblings.

    See Movement for full documentation of the shared interface.

    CoreMovement holds the fundamental parameters and shared derived properties that all
    Movement variants need. Unlike Movement, CoreMovement requires delta_time and
    num_steps to be provided directly and does not perform automatic estimation or batch
    generation.
    """

    __slots__ = (
        "_airplane_movements",
        "_operating_point_movement",
        "_delta_time",
        "_num_steps",
        "_max_wake_rows",
        "_lcm_period",
        "_max_period",
        "_min_period",
        "_static",
    )

    def __init__(
        self,
        airplane_movements: Sequence[CoreAirplaneMovement],
        operating_point_movement: CoreOperatingPointMovement,
        delta_time: float | int,
        num_steps: int,
        max_wake_rows: int | None = None,
    ) -> None:
        """The initialization method.

        See Movement's initialization method for full parameter descriptions.

        :param airplane_movements: The CoreAirplaneMovements for each Airplane.
        :param operating_point_movement: The CoreOperatingPointMovement for operating
            conditions.
        :param delta_time: The time step size in seconds. Must be positive.
        :param num_steps: The number of time steps. Must be a positive int.
        :param max_wake_rows: The maximum chordwise wake rows per Wing. The default is
            None (no truncation).
        :return: None
        """
        # Validate and store the CoreAirplaneMovements.
        if not isinstance(airplane_movements, list):
            raise TypeError("airplane_movements must be a list.")
        if len(airplane_movements) < 1:
            raise ValueError("airplane_movements must have at least one element.")
        for airplane_movement in airplane_movements:
            if not isinstance(airplane_movement, CoreAirplaneMovement):
                raise TypeError(
                    "Every element in airplane_movements must be a "
                    "CoreAirplaneMovement."
                )
        # Store as tuple to prevent external mutation.
        self._airplane_movements: tuple[CoreAirplaneMovement, ...] = tuple(
            airplane_movements
        )

        # Validate and store the CoreOperatingPointMovement.
        if not isinstance(
            operating_point_movement,
            CoreOperatingPointMovement,
        ):
            raise TypeError(
                "operating_point_movement must be a " "CoreOperatingPointMovement."
            )
        self._operating_point_movement = operating_point_movement

        # Initialize the caches for the properties derived from the immutable
        # attributes. These are initialized early because static is accessed
        # below during __init__ to validate max_wake_rows.
        self._lcm_period: float | None = None
        self._max_period: float | None = None
        self._min_period: float | None = None
        self._static: bool | None = None

        # Validate and store delta_time.
        delta_time = _parameter_validation.number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )
        self._delta_time: float = delta_time

        # Validate and store num_steps.
        num_steps = _parameter_validation.int_in_range_return_int(
            num_steps,
            "num_steps",
            min_val=1,
            min_inclusive=True,
        )
        self._num_steps: int = num_steps

        # Validate and store max_wake_rows.
        if max_wake_rows is not None:
            max_wake_rows = _parameter_validation.int_in_range_return_int(
                max_wake_rows,
                "max_wake_rows",
                min_val=1,
                min_inclusive=True,
            )
        self._max_wake_rows = max_wake_rows

    # --- Immutable: read only properties ---
    @property
    def airplane_movements(
        self,
    ) -> tuple[CoreAirplaneMovement, ...]:
        return self._airplane_movements

    @property
    def operating_point_movement(
        self,
    ) -> CoreOperatingPointMovement:
        return self._operating_point_movement

    @property
    def delta_time(self) -> float:
        return self._delta_time

    @property
    def num_steps(self) -> int:
        return self._num_steps

    @property
    def max_wake_rows(self) -> int | None:
        return self._max_wake_rows

    # --- Immutable derived: manual lazy caching ---
    @property
    def lcm_period(self) -> float:
        """The least common multiple of all motion periods, ensuring all motions
        complete an integer number of cycles when cycle averaging forces and moments.

        Using the LCM ensures that when cycle averaging forces and moments, we capture a
        complete cycle of all motions, not just the longest one. For example, if one
        motion has a period of 2.0 s and another has a period of 3.0 s, the LCM is 6.0,
        which contains exactly 3 cycles of the first motion and 2 cycles of the second.

        :return: The LCM period in seconds. If all the motion is static, this will be
            0.0.
        """
        if self._lcm_period is None:
            # Collect all periods from AirplaneMovements.
            all_periods: list[float] = []
            for airplane_movement in self._airplane_movements:
                all_periods.extend(airplane_movement.all_periods)

            # Add the OperatingPointMovement period.
            all_periods.append(self._operating_point_movement.max_period)

            self._lcm_period = lcm_multiple(all_periods)
        return self._lcm_period

    @property
    def max_period(self) -> float:
        """The longest period of motion across each Airplane's movement class, each
        Wing's movement class, each WingCrossSection's movement class, and the
        OperatingPoint's movement class.

        For cycle averaging calculations, lcm_period should be used instead of
        max_period to ensure all motions complete an integer number of cycles.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        if self._max_period is None:
            # Iterate through the AirplaneMovements and find the one with the
            # largest max period.
            airplane_movement_max_periods = []
            for airplane_movement in self._airplane_movements:
                airplane_movement_max_periods.append(airplane_movement.max_period)
            max_airplane_period = max(airplane_movement_max_periods)

            # The global max period is the maximum of the max
            # AirplaneMovement period and the OperatingPointMovement max
            # period.
            self._max_period = max(
                max_airplane_period,
                self._operating_point_movement.max_period,
            )
        return self._max_period

    @property
    def min_period(self) -> float:
        """The shortest non zero period of motion across each Airplane's movement class,
        each Wing's movement class, each WingCrossSection's movement class, and the
        OperatingPoint's movement class.

        :return: The shortest non zero period in seconds. If all the motion is static,
            this will be 0.0.
        """
        if self._min_period is None:
            # Collect all periods from AirplaneMovements.
            all_periods: list[float] = []
            for airplane_movement in self._airplane_movements:
                all_periods.extend(airplane_movement.all_periods)

            # Add the OperatingPointMovement period.
            op_period = self._operating_point_movement.max_period
            if op_period != 0.0:
                all_periods.append(op_period)

            # Filter out zero periods and find the minimum.
            non_zero_periods = [p for p in all_periods if p != 0.0]
            if not non_zero_periods:
                self._min_period = 0.0
            else:
                self._min_period = min(non_zero_periods)
        return self._min_period

    @property
    def static(self) -> bool:
        """Flags if all motion across each Airplane's movement class, each Wing's
        movement class, each WingCrossSection's movement class, and the OperatingPoint's
        movement class is static.

        :return: True if all motion is static. False otherwise.
        """
        if self._static is None:
            self._static = self.max_period == 0
        return self._static


class CoreUnsteadyProblem:
    """A core class used to contain the shared foundation of UnsteadyProblem and its
    feature variant siblings.

    See UnsteadyProblem for full documentation of the shared interface.

    CoreUnsteadyProblem holds the time stepping parameters, wake truncation setting,
    result storage mode, and the mutable load result lists that the solver populates.
    Unlike UnsteadyProblem, it does not take a Movement or pre create SteadyProblems.
    Feature variants (FreeFlightUnsteadyProblem, AeroelasticUnsteadyProblem) extend this
    class and provide SteadyProblems dynamically at each time step.
    """

    __slots__ = (
        "_only_final_results",
        "_num_steps",
        "_delta_time",
        "_max_wake_rows",
        "_first_averaging_step",
        "_first_results_step",
        "finalForces_W",
        "finalForceCoefficients_W",
        "finalMoments_W_CgP1",
        "finalMomentCoefficients_W_CgP1",
        "finalMeanForces_W",
        "finalMeanForceCoefficients_W",
        "finalMeanMoments_W_CgP1",
        "finalMeanMomentCoefficients_W_CgP1",
        "finalRmsForces_W",
        "finalRmsForceCoefficients_W",
        "finalRmsMoments_W_CgP1",
        "finalRmsMomentCoefficients_W_CgP1",
    )

    def __init__(
        self,
        only_final_results: bool | np.bool_,
        delta_time: float | int,
        num_steps: int,
        max_wake_rows: int | None,
        lcm_period: float | int,
    ) -> None:
        """The initialization method.

        See UnsteadyProblem's initialization method for full parameter descriptions.

        :param only_final_results: Determines whether the solver will only calculate
            loads for the final time step or final cycle.
        :param delta_time: The time step size in seconds. Must be positive.
        :param num_steps: The number of time steps. Must be a positive int.
        :param max_wake_rows: The maximum chordwise wake rows per Wing. None means no
            truncation.
        :param lcm_period: The least common multiple of all motion periods in seconds.
            Used to compute the first averaging step. Must be non negative. A value of
            0.0 indicates static motion.
        :return: None
        """
        # Validate and store immutable attributes.
        self._only_final_results: bool = _parameter_validation.boolLike_return_bool(
            only_final_results, "only_final_results"
        )

        self._delta_time: float = _parameter_validation.number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )

        self._num_steps: int = _parameter_validation.int_in_range_return_int(
            num_steps, "num_steps", min_val=1, min_inclusive=True
        )

        if max_wake_rows is not None:
            max_wake_rows = _parameter_validation.int_in_range_return_int(
                max_wake_rows,
                "max_wake_rows",
                min_val=1,
                min_inclusive=True,
            )
        self._max_wake_rows: int | None = max_wake_rows

        lcm_period = _parameter_validation.number_in_range_return_float(
            lcm_period, "lcm_period", min_val=0.0, min_inclusive=True
        )

        # For CoreUnsteadyProblems with a static CoreMovement, we are typically
        # interested in the final time step's forces and moments, which, assuming
        # convergence, will be the most accurate. For CoreUnsteadyProblems with cyclic
        # movement (e.g., flapping wings), we are typically interested in the forces
        # and moments averaged over the last cycle simulated. Use the LCM of all motion
        # periods to ensure we average over a complete cycle of all motions.
        self._first_averaging_step: int
        if lcm_period == 0:
            self._first_averaging_step = self._num_steps - 1
        else:
            self._first_averaging_step = max(
                0,
                math.floor(self._num_steps - (lcm_period / self._delta_time)),
            )

        # If we only want to calculate forces and moments for the final cycle (for
        # cyclic motion) or for the final time step (for static motion), set the first
        # step to calculate results to the first averaging step. Otherwise, set it to
        # zero, which is the first time step.
        self._first_results_step: int
        if self._only_final_results:
            self._first_results_step = self._first_averaging_step
        else:
            self._first_results_step = 0

        # Initialize empty lists to hold the final loads and load coefficients each
        # Airplane experiences. These will only be populated if this
        # CoreUnsteadyProblem's motion is static. These are mutable and populated by
        # the solver.
        self.finalForces_W: list[np.ndarray] = []
        self.finalForceCoefficients_W: list[np.ndarray] = []
        self.finalMoments_W_CgP1: list[np.ndarray] = []
        self.finalMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize empty lists to hold the final cycle averaged loads and load
        # coefficients each Airplane experiences. These will only be populated if this
        # CoreUnsteadyProblem's motion is cyclic. These are mutable and populated by
        # the solver.
        self.finalMeanForces_W: list[np.ndarray] = []
        self.finalMeanForceCoefficients_W: list[np.ndarray] = []
        self.finalMeanMoments_W_CgP1: list[np.ndarray] = []
        self.finalMeanMomentCoefficients_W_CgP1: list[np.ndarray] = []

        # Initialize empty lists to hold the final cycle root mean squared loads and
        # load coefficients each Airplane experiences. These will only be populated for
        # variable geometry problems. These are mutable and populated by the solver.
        self.finalRmsForces_W: list[np.ndarray] = []
        self.finalRmsForceCoefficients_W: list[np.ndarray] = []
        self.finalRmsMoments_W_CgP1: list[np.ndarray] = []
        self.finalRmsMomentCoefficients_W_CgP1: list[np.ndarray] = []

    # --- Immutable: read only properties ---
    @property
    def only_final_results(self) -> bool:
        return self._only_final_results

    @property
    def num_steps(self) -> int:
        return self._num_steps

    @property
    def delta_time(self) -> float:
        return self._delta_time

    @property
    def first_averaging_step(self) -> int:
        return self._first_averaging_step

    @property
    def first_results_step(self) -> int:
        return self._first_results_step

    @property
    def max_wake_rows(self) -> int | None:
        return self._max_wake_rows

    @property
    def movement(self) -> CoreMovement:
        # This stub lets the UnsteadyRingVortexLatticeMethodSolver access movement and
        # steady_problems on any CoreUnsteadyProblem without knowing the concrete
        # subclass.
        raise NotImplementedError(
            "Subclasses of CoreUnsteadyProblem must override the movement property."
        )

    @property
    def steady_problems(self) -> tuple[problems.SteadyProblem, ...]:
        # This stub lets the UnsteadyRingVortexLatticeMethodSolver access movement and
        # steady_problems on any CoreUnsteadyProblem without knowing the concrete
        # subclass.
        raise NotImplementedError(
            "Subclasses of CoreUnsteadyProblem must override the steady_problems "
            "property."
        )
