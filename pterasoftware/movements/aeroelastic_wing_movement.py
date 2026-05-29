"""Contains the AeroelasticWingMovement class.

**Contains the following classes:**

AeroelasticWingMovement: A class used to contain a Wing's movement in an aeroelastic
simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import cast

import numpy as np

from .. import _core, _oscillation, _transformations, geometry
from . import (
    aeroelastic_wing_cross_section_movement as aeroelastic_wing_cross_section_movement_mod,
)


class AeroelasticWingMovement(_core.CoreWingMovement):
    """A class used to contain a Wing's movement in an aeroelastic simulation.

    In aeroelastic simulations, wing geometry is prescribed via oscillation parameters
    (the same oscillation based generation as WingMovement), but the solver adds
    structural deformation at each time step. This class overrides
    generate_wing_at_time_step to accept per WingCrossSection deformation angles that
    are threaded down to its AeroelasticWingCrossSectionMovement children.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this AeroelasticWingMovement.

    all_periods: All unique non zero periods from this AeroelasticWingMovement and its
    AeroelasticWingCrossSectionMovements.

    max_period: The longest period of AeroelasticWingMovement's own motion and that of
    its sub movement objects.

    generate_wing_at_time_step: Creates the Wing at a single time step, optionally
    applying structural deformation to each WingCrossSection.

    generate_wings: Creates the Wing at each time step, and returns them in a list.

    **Notes:**

    Wings cannot undergo motion that causes them to switch symmetry types. A transition
    between types could change the number of Wings and the Panel structure, which is
    incompatible with the unsteady solver. This happens when an AeroelasticWingMovement
    defines motion that causes its base Wing's wing axes' yz plane and its symmetry
    plane to transition from coincident to non coincident, or vice versa. This is
    checked by this AeroelasticWingMovement's parent AeroelasticAirplaneMovement's
    parent AeroelasticMovement.
    """

    __slots__ = ("_spacingAnglesSecondDerivative_Gs_to_Wn_ixyz",)

    def __init__(
        self,
        base_wing: geometry.wing.Wing,
        wing_cross_section_movements: list[
            aeroelastic_wing_cross_section_movement_mod.AeroelasticWingCrossSectionMovement
        ],
        ampLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs: np.ndarray | Sequence[str | Callable[[float], float]] = (
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
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz: (
            Sequence[Callable[[float], float] | None] | None
        ) = None,
    ) -> None:
        """The initialization method.

        :param base_wing: The base Wing from which the Wing at each time step will be
            created.
        :param wing_cross_section_movements: A list of
            AeroelasticWingCrossSectionMovements associated with each of the base Wing's
            WingCrossSections. It must have the same length as the base Wing's list of
            WingCrossSections.
        :param ampLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            AeroelasticWingMovement's changes in its Wings' Ler_Gs_Cgs parameters. Can
            be a tuple, list, or ndarray. Values are converted to floats internally.
            Each amplitude must be low enough that it doesn't drive its base value out
            of the range of valid values. Otherwise, this AeroelasticWingMovement will
            try to create Wings with invalid parameters values. The units are in meters.
            The default is (0.0, 0.0, 0.0).
        :param periodLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            AeroelasticWingMovement's changes in its Wings' Ler_Gs_Cgs parameters. Can
            be a tuple, list, or ndarray. Values are converted to floats internally.
            Each element must be 0.0 if the corresponding element in ampLer_Gs_Cgs is
            0.0 and non zero if not. The units are in seconds. The default is (0.0, 0.0,
            0.0).
        :param spacingLer_Gs_Cgs: An array-like object of strs or callables with shape
            (3,) representing the spacing of the AeroelasticWingMovement's change in its
            Wings' Ler_Gs_Cgs parameters. Can be a tuple, list, or ndarray. Each element
            can be the string "sine", the string "uniform", or a callable custom spacing
            function. Custom spacing functions are for advanced users and must start at
            0.0, return to 0.0 after one period of 2.0 * pi radians, have amplitude of
            1.0, be periodic, return finite values only, and accept a float as input and
            return a float. The custom function is scaled by ampLer_Gs_Cgs, shifted
            horizontally and vertically by phaseLer_Gs_Cgs and the base value, and have
            a period set by periodLer_Gs_Cgs. The default is ("sine", "sine", "sine").
        :param phaseLer_Gs_Cgs: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's Wing's Ler_Gs_Cgs parameter relative to the base Wing's Ler_Gs_Cgs
            parameter. Can be a tuple, list, or ndarray. Values must lie in the range
            (-180.0, 180.0] and will be converted to floats internally. Each element
            must be 0.0 if the corresponding element in ampLer_Gs_Cgs is 0.0 and non
            zero if not. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param ampAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or float)
            with shape (3,) representing the amplitudes of the AeroelasticWingMovement's
            changes in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Values must lie in the range [0.0, 180.0] and will be converted
            to floats internally. Each amplitude must be low enough that it doesn't
            drive its base value out of the range of valid values. Otherwise, this
            AeroelasticWingMovement will try to create Wings with invalid parameters
            values. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param periodAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the periods of the
            AeroelasticWingMovement's changes in its Wings' angles_Gs_to_Wn_ixyz
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampAngles_Gs_to_Wn_ixyz is 0.0 and non zero if not. The units are in
            seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Gs_to_Wn_ixyz: An array-like object of strs or callables
            with shape (3,) representing the spacing of the AeroelasticWingMovement's
            change in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Each element can be the string "sine", the string "uniform", or
            a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0.0, return to 0.0 after one period of 2.0
            * pi radians, have amplitude of 1.0, be periodic, return finite values only,
            and accept a float as input and return a float. The custom function is
            scaled by ampAngles_Gs_to_Wn_ixyz, shifted horizontally and vertically by
            phaseAngles_Gs_to_Wn_ixyz and the base value, with the period set by
            periodAngles_Gs_to_Wn_ixyz. A component set to a custom callable must be
            paired with a matching spacingAnglesSecondDerivative_Gs_to_Wn_ixyz
            component, and a "sine" or "uniform" component must not be; see that
            parameter for the full pairing rule. The default is ("sine", "sine",
            "sine").
        :param phaseAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or float)
            with shape (3,) representing the phase offsets of the elements in the first
            time step's Wing's angles_Gs_to_Wn_ixyz parameter relative to the base
            Wing's angles_Gs_to_Wn_ixyz parameter. Can be a tuple, list, or ndarray.
            Values must lie in the range (-180.0, 180.0] and will be converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampAngles_Gs_to_Wn_ixyz is 0.0 and non zero if not. The units are in
            degrees. The default is (0.0, 0.0, 0.0).
        :param rotationPointOffset_Gs_Ler: An array-like object of 3 numbers (int or
            float) representing the position of the rotation point for the Wing's
            angular motion (in geometry axes after accounting for symmetry, relative to
            the leading edge root point). Can be a tuple, list, or ndarray. Values are
            converted to floats internally. This offset defines where the Wing rotates
            about when angles_Gs_to_Wn_ixyz oscillates. When set to (0, 0, 0), rotation
            occurs about the leading edge root point (default behavior). The units are
            in meters. The default is (0.0, 0.0, 0.0).
        :param spacingAnglesSecondDerivative_Gs_to_Wn_ixyz: An optional sequence with
            shape (3,) holding the analytical second time derivative of each
            spacingAngles_Gs_to_Wn_ixyz component, used by the aeroelastic solver to
            compute the inertial torque from the prescribed flapping acceleration. Each
            element is either a callable that accepts a time (in seconds) and returns
            the second derivative (in radians per second squared, before amplitude
            scaling), or None when the corresponding spacing component does not have
            one. Under the current model only the x (flap) component is consulted. Each
            component must agree with its matching spacingAngles_Gs_to_Wn_ixyz
            component: a custom (callable) spacing must have a non-None derivative here,
            and a "sine" or "uniform" spacing must have None here (their derivatives are
            handled analytically or rejected as non-differentiable when the torque is
            generated, so a supplied derivative would be ignored). Either mismatch
            raises a ValueError. When None, every component is None, which is valid only
            when no spacingAngles_Gs_to_Wn_ixyz component is a custom callable. The
            default is None.
        :return: None
        """
        # Validate that every element is an AeroelasticWingCrossSectionMovement,
        # not just a CoreWingCrossSectionMovement. CoreWingMovement.__init__()
        # validates at the Core level, but AeroelasticWingMovement enforces the
        # stricter type.
        for wing_cross_section_movement in wing_cross_section_movements:
            if not isinstance(
                wing_cross_section_movement,
                aeroelastic_wing_cross_section_movement_mod.AeroelasticWingCrossSectionMovement,
            ):
                raise TypeError(
                    "Every element in wing_cross_section_movements must "
                    "be an AeroelasticWingCrossSectionMovement."
                )

        super().__init__(
            base_wing=base_wing,
            wing_cross_section_movements=wing_cross_section_movements,
            ampLer_Gs_Cgs=ampLer_Gs_Cgs,
            periodLer_Gs_Cgs=periodLer_Gs_Cgs,
            spacingLer_Gs_Cgs=spacingLer_Gs_Cgs,
            phaseLer_Gs_Cgs=phaseLer_Gs_Cgs,
            ampAngles_Gs_to_Wn_ixyz=ampAngles_Gs_to_Wn_ixyz,
            periodAngles_Gs_to_Wn_ixyz=periodAngles_Gs_to_Wn_ixyz,
            spacingAngles_Gs_to_Wn_ixyz=spacingAngles_Gs_to_Wn_ixyz,
            phaseAngles_Gs_to_Wn_ixyz=phaseAngles_Gs_to_Wn_ixyz,
            rotationPointOffset_Gs_Ler=rotationPointOffset_Gs_Ler,
        )

        # Validate the second-derivative companion to spacingAngles_Gs_to_Wn_ixyz.
        if spacingAnglesSecondDerivative_Gs_to_Wn_ixyz is None:
            derivatives: tuple[Callable[[float], float] | None, ...] = (
                None,
                None,
                None,
            )
        else:
            if (
                not isinstance(
                    spacingAnglesSecondDerivative_Gs_to_Wn_ixyz, (list, tuple)
                )
                or len(spacingAnglesSecondDerivative_Gs_to_Wn_ixyz) != 3
            ):
                raise ValueError(
                    "spacingAnglesSecondDerivative_Gs_to_Wn_ixyz must be None or a "
                    "3-element sequence."
                )
            for i, derivative in enumerate(spacingAnglesSecondDerivative_Gs_to_Wn_ixyz):
                if derivative is not None and not callable(derivative):
                    raise TypeError(
                        "Each element of spacingAnglesSecondDerivative_Gs_to_Wn_ixyz "
                        "must be a callable or None, got "
                        f"{type(derivative).__name__} at index {i}."
                    )
            derivatives = tuple(spacingAnglesSecondDerivative_Gs_to_Wn_ixyz)

        # The second derivative is meaningful only for a custom (callable) spacing
        # component, so the two must agree per component. A callable spacing has no
        # analytical derivative the solver can take, so it must be paired with one. A
        # named ("sine" or "uniform") spacing already has its derivative handled
        # analytically or rejected as non-differentiable when the torque is generated,
        # so a supplied derivative would be silently ignored. Reject either mismatch
        # here rather than when the torque is generated.
        for i, spacing in enumerate(self.spacingAngles_Gs_to_Wn_ixyz):
            spacing_is_callable = callable(spacing)
            if spacing_is_callable and derivatives[i] is None:
                raise ValueError(
                    "A custom (callable) spacingAngles_Gs_to_Wn_ixyz requires a "
                    "matching spacingAnglesSecondDerivative_Gs_to_Wn_ixyz, but element "
                    f"{i} has a callable spacing with no derivative."
                )
            if not spacing_is_callable and derivatives[i] is not None:
                raise ValueError(
                    "A spacingAnglesSecondDerivative_Gs_to_Wn_ixyz may only be given "
                    "for a custom (callable) spacingAngles_Gs_to_Wn_ixyz, but element "
                    f"{i} has a '{spacing}' spacing with a derivative."
                )

        self._spacingAnglesSecondDerivative_Gs_to_Wn_ixyz = derivatives

    def __deepcopy__(self, memo: dict) -> AeroelasticWingMovement:
        """Creates a deep copy of this AeroelasticWingMovement.

        Extends the parent deep copy to also copy this class's second-derivative
        companion to the angular spacing.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new instance with copied attributes.
        """
        new_movement = cast(AeroelasticWingMovement, super().__deepcopy__(memo))

        # Copy the tuple directly (immutable; holds callables or None per component).
        new_movement._spacingAnglesSecondDerivative_Gs_to_Wn_ixyz = (
            self._spacingAnglesSecondDerivative_Gs_to_Wn_ixyz
        )

        return new_movement

    @property
    def spacingAnglesSecondDerivative_Gs_to_Wn_ixyz(
        self,
    ) -> tuple[Callable[[float], float] | None, ...]:
        return self._spacingAnglesSecondDerivative_Gs_to_Wn_ixyz

    def generate_wing_at_time_step(
        self,
        step: int,
        delta_time: float | int,
        deformation_angles_ixyz: np.ndarray | None = None,
    ) -> geometry.wing.Wing:
        """Creates the Wing at a single time step, optionally applying structural
        deformation to each WingCrossSection.

        Computes the prescribed Wing using the inherited oscillation logic, then threads
        per WingCrossSection deformation angles down to each
        AeroelasticWingCrossSectionMovement child.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :param deformation_angles_ixyz: A (N, 3) ndarray of floats where N is the number
            of WingCrossSections in this Wing. Each row is a (3,) deformation angle
            vector using an intrinsic xy'z" sequence that is added to the corresponding
            WingCrossSection's prescribed angles_Wcsp_to_Wcs_ixyz. The units are in
            degrees. When None, no deformation is applied. The default is None.
        :return: The Wing at this time step, with structural deformation applied to each
            WingCrossSection if provided.
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

        # Generate the WingCrossSections for this time step, threading
        # deformation to each AeroelasticWingCrossSectionMovement child.
        these_wing_cross_sections = []
        for i, wing_cross_section_movement in enumerate(
            self._wing_cross_section_movements
        ):
            assert isinstance(
                wing_cross_section_movement,
                aeroelastic_wing_cross_section_movement_mod.AeroelasticWingCrossSectionMovement,
            )

            # Extract this WingCrossSection's deformation row, or None.
            this_deformation = None
            if deformation_angles_ixyz is not None:
                this_deformation = deformation_angles_ixyz[i]

            these_wing_cross_sections.append(
                wing_cross_section_movement.generate_wing_cross_section_at_time_step(
                    step,
                    delta_time,
                    deformation_angles_ixyz=this_deformation,
                )
            )

        # If there is a non zero rotation point offset, adjust the position
        # to account for rotation about the offset point instead of the
        # leading edge root.
        if not np.allclose(self._rotationPointOffset_Gs_Ler, np.zeros(3, dtype=float)):
            rot_T_act = _transformations.generate_rot_T(
                theseAngles_Gs_to_Wn_ixyz,
                passive=False,
                intrinsic=True,
                order="xyz",
            )
            rot_R_act = rot_T_act[:3, :3]

            offsetRotationPointAdjustment_Gs = (
                np.eye(3, dtype=float) - rot_R_act
            ) @ self._rotationPointOffset_Gs_Ler

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
