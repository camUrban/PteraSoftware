"""Contains the AeroelasticWingCrossSectionMovement class.

**Contains the following classes:**

AeroelasticWingCrossSectionMovement: A class used to contain a WingCrossSection's
movement in an aeroelastic simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from .. import _core, _oscillation, geometry


class AeroelasticWingCrossSectionMovement(_core.CoreWingCrossSectionMovement):
    """A class used to contain a WingCrossSection's movement in an aeroelastic
    simulation.

    In aeroelastic simulations, wing cross section geometry is prescribed via
    oscillation parameters (the same oscillation based generation as
    WingCrossSectionMovement), but the solver adds structural deformation angles at each
    time step. This class overrides generate_wing_cross_section_at_time_step to accept
    an optional deformation that is added to the prescribed angles_Wcsp_to_Wcs_ixyz.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this AeroelasticWingCrossSectionMovement.

    all_periods: All unique non zero periods from this
    AeroelasticWingCrossSectionMovement.

    max_period: AeroelasticWingCrossSectionMovement's longest period of motion.

    generate_wing_cross_section_at_time_step: Creates the WingCrossSection at a single
    time step, optionally applying structural deformation.

    generate_wing_cross_sections: Creates the WingCrossSection at each time step, and
    returns them in a list.
    """

    __slots__ = ()

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

        :param base_wing_cross_section: The base WingCrossSection from which the
            WingCrossSection at each time step will be created.
        :param ampLp_Wcsp_Lpp: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            AeroelasticWingCrossSectionMovement's changes in its WingCrossSections'
            Lp_Wcsp_Lpp parameters. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Each amplitude must be low enough that it
            doesn't drive its base value out of the range of valid values. Otherwise,
            this AeroelasticWingCrossSectionMovement will try to create
            WingCrossSections with invalid parameter values. The units are in meters.
            The default is (0.0, 0.0, 0.0).
        :param periodLp_Wcsp_Lpp: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            AeroelasticWingCrossSectionMovement's changes in its WingCrossSections'
            Lp_Wcsp_Lpp parameters. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampLp_Wcsp_Lpp is 0.0 and non zero if not. The
            units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingLp_Wcsp_Lpp: An array-like object of strs or callables with shape
            (3,) representing the spacing of the AeroelasticWingCrossSectionMovement's
            changes in its WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple,
            list, or ndarray. Each element can be the str "sine", the str "uniform", or
            a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0.0, return to 0.0 after one period of 2.0
            * pi radians, have amplitude of 1.0, be periodic, return finite values only,
            and accept a float as input and return a float. Custom functions are scaled
            by ampLp_Wcsp_Lpp, shifted horizontally and vertically by phaseLp_Wcsp_Lpp
            and the base value, and have a period set by periodLp_Wcsp_Lpp. The default
            is ("sine", "sine", "sine").
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
            AeroelasticWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Values
            are converted to floats internally. Each amplitude must be low enough that
            it doesn't drive its base value out of the range of valid values. Otherwise,
            this AeroelasticWingCrossSectionMovement will try to create
            WingCrossSections with invalid parameter values. The units are in degrees.
            The default is (0.0, 0.0, 0.0).
        :param periodAngles_Wcsp_to_Wcs_ixyz: An array-like object of non negative
            numbers (int or float) with shape (3,) representing the periods of the
            AeroelasticWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Values
            are converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non zero if
            not. The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Wcsp_to_Wcs_ixyz: An array-like object of strs or callables
            with shape (3,) representing the spacing of the
            AeroelasticWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Each
            element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2.0 * pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a float
            as input and return a float. Custom functions are scaled by
            ampAngles_Wcsp_to_Wcs_ixyz, shifted horizontally and vertically by
            phaseAngles_Wcsp_to_Wcs_ixyz and the base value, and have a period set by
            periodAngles_Wcsp_to_Wcs_ixyz. The default is ("sine", "sine", "sine").
        :param phaseAngles_Wcsp_to_Wcs_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the phase offsets of the elements in the
            first time step's WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter
            relative to the base WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter.
            Can be a tuple, list, or ndarray. Elements must lie in the range (-180.0,
            180.0]. Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non zero if not. Values are converted
            to floats internally. The units are in degrees. The default is (0.0, 0.0,
            0.0).
        :return: None
        """
        super().__init__(
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

    def generate_wing_cross_section_at_time_step(
        self,
        step: int,
        delta_time: float | int,
        deformation_angles_ixyz: np.ndarray | None = None,
    ) -> geometry.wing_cross_section.WingCrossSection:
        """Creates the WingCrossSection at a single time step, optionally applying
        structural deformation.

        Computes the prescribed WingCrossSection using the inherited oscillation logic,
        then adds the deformation angles to the prescribed angles_Wcsp_to_Wcs_ixyz. This
        is the integration point where the aeroelastic solver's structural response
        modifies the wing cross section geometry.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :param deformation_angles_ixyz: A (3,) ndarray of floats representing the
            structural deformation angles to add to the prescribed
            angles_Wcsp_to_Wcs_ixyz, using an intrinsic xy'z" sequence. The units are in
            degrees. When None, no deformation is applied and the result is identical to
            the prescribed WingCrossSection. The default is None.
        :return: The WingCrossSection at this time step, with structural deformation
            applied if provided.
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

        # Apply structural deformation if provided. The deformation angles
        # are added to the prescribed angles, representing the structural
        # response computed by the aeroelastic solver.
        if deformation_angles_ixyz is not None:
            theseAngles_Wcsp_to_Wcs_ixyz = (
                theseAngles_Wcsp_to_Wcs_ixyz + deformation_angles_ixyz
            )

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
