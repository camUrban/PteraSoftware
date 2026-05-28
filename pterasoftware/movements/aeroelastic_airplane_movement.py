"""Contains the AeroelasticAirplaneMovement class.

**Contains the following classes:**

AeroelasticAirplaneMovement: A class used to contain an Airplane's movement in an
aeroelastic simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import cast

import numpy as np

from .. import _core, _oscillation, geometry
from . import aeroelastic_wing_movement as aeroelastic_wing_movement_mod
from . import wing_movement as wing_movement_mod


class AeroelasticAirplaneMovement(_core.CoreAirplaneMovement):
    """A class used to contain an Airplane's movement in an aeroelastic simulation.

    In aeroelastic simulations, airplane geometry is prescribed via oscillation
    parameters (the same oscillation based generation as AirplaneMovement), but the
    solver adds structural deformation at each time step. This class overrides
    generate_airplane_at_time_step to accept per Wing deformation that is threaded down
    to its AeroelasticWingMovement children.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this AeroelasticAirplaneMovement.

    all_periods: All unique non zero periods from this AeroelasticAirplaneMovement, its
    AeroelasticWingMovement(s), and their AeroelasticWingCrossSectionMovements.

    max_period: The longest period of AeroelasticAirplaneMovement's own motion, the
    motion(s) of its sub movement object(s), and the motions of its sub sub movement
    objects.

    generate_airplane_at_time_step: Creates the Airplane at a single time step,
    optionally applying structural deformation to each Wing.

    generate_airplanes: Creates the Airplane at each time step, and returns them in a
    list.
    """

    __slots__ = ()

    def __init__(
        self,
        base_airplane: geometry.airplane.Airplane,
        wing_movements: list[
            aeroelastic_wing_movement_mod.AeroelasticWingMovement
            | wing_movement_mod.WingMovement
        ],
        ampCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
    ) -> None:
        """The initialization method.

        :param base_airplane: The base Airplane from which the Airplane at each time
            step will be created.
        :param wing_movements: A list of the wing movements associated with each of the
            base Airplane's Wings. Each element must be either an
            AeroelasticWingMovement (which will receive structural deformation at each
            time step) or a WingMovement (which will be advanced without deformation).
            The list must have the same length as the base Airplane's list of Wings.
        :param ampCg_GP1_CgP1: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            AeroelasticAirplaneMovement's changes in its Airplanes' Cg_GP1_CgP1
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each amplitude must be low enough that it doesn't drive its base
            value out of the range of valid values. Otherwise, this
            AeroelasticAirplaneMovement will try to create Airplanes with invalid
            parameter values. Because the first Airplane's Cg_GP1_CgP1 parameter must be
            all zeros, this means that the first Airplane's ampCg_GP1_CgP1 parameter
            must also be all zeros. The units are in meters. The default is (0.0, 0.0,
            0.0).
        :param periodCg_GP1_CgP1: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            AeroelasticAirplaneMovement's changes in its Airplanes' Cg_GP1_CgP1
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampCg_GP1_CgP1 is 0.0 and non zero if not. The units are in seconds. The
            default is (0.0, 0.0, 0.0).
        :param spacingCg_GP1_CgP1: An array-like object of strs or callables with shape
            (3,) representing the spacing of the AeroelasticAirplaneMovement's changes
            in its Airplanes' Cg_GP1_CgP1 parameters. Can be a tuple, list, or ndarray.
            Each element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2.0 * pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a float
            as input and return a float. Custom functions are scaled by ampCg_GP1_CgP1,
            shifted horizontally and vertically by phaseCg_GP1_CgP1 and the base value,
            and have a period set by periodCg_GP1_CgP1. The default is ("sine", "sine",
            "sine").
        :param phaseCg_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's Airplane's Cg_GP1_CgP1 parameter relative to the base Airplane's
            Cg_GP1_CgP1 parameter. Can be a tuple, list, or ndarray. Elements must lie
            in the range (-180.0, 180.0]. Each element must be 0.0 if the corresponding
            element in ampCg_GP1_CgP1 is 0.0 and non zero if not. Values are converted
            to floats internally. The units are in degrees. The default is (0.0, 0.0,
            0.0).
        :return: None
        """
        # Validate that every element is an AeroelasticWingMovement or a
        # WingMovement. CoreAirplaneMovement.__init__() validates at the Core level,
        # but AeroelasticAirplaneMovement enforces this stricter requirement.
        for wing_movement in wing_movements:
            if not isinstance(
                wing_movement,
                (
                    aeroelastic_wing_movement_mod.AeroelasticWingMovement,
                    wing_movement_mod.WingMovement,
                ),
            ):
                raise TypeError(
                    "Every element in wing_movements must be an "
                    "AeroelasticWingMovement or a WingMovement."
                )

        super().__init__(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=ampCg_GP1_CgP1,
            periodCg_GP1_CgP1=periodCg_GP1_CgP1,
            spacingCg_GP1_CgP1=spacingCg_GP1_CgP1,
            phaseCg_GP1_CgP1=phaseCg_GP1_CgP1,
        )

    # --- Immutable: read only properties ---
    @property
    def wing_movements(
        self,
    ) -> tuple[
        aeroelastic_wing_movement_mod.AeroelasticWingMovement
        | wing_movement_mod.WingMovement,
        ...,
    ]:
        return cast(
            tuple[
                aeroelastic_wing_movement_mod.AeroelasticWingMovement
                | wing_movement_mod.WingMovement,
                ...,
            ],
            self._wing_movements,
        )

    def generate_airplane_at_time_step(
        self,
        step: int,
        delta_time: float | int,
        wing_deformation_angles_ixyz: list[np.ndarray | None] | None = None,
    ) -> geometry.airplane.Airplane:
        """Creates the Airplane at a single time step, optionally applying structural
        deformation to each Wing.

        Computes the prescribed Airplane using the inherited oscillation logic, then
        threads per Wing deformation down to each AeroelasticWingMovement child.

        :param step: The time step index. Must be a non negative int.
        :param delta_time: The time between each time step in seconds. Must be a
            positive number (int or float).
        :param wing_deformation_angles_ixyz: A list of (N_wcs, 3) ndarrays of floats,
            one per Wing, where N_wcs is the number of WingCrossSections in that Wing.
            Each row is a (3,) deformation angle vector using an intrinsic xy'z"
            sequence. The units are in degrees. When None, no deformation is applied.
            The default is None.
        :return: The Airplane at this time step, with structural deformation applied to
            each Wing if provided.
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

        # Generate the Wings for this time step, threading deformation to
        # each AeroelasticWingMovement child. WingMovement children are advanced
        # without deformation.
        these_wings = []
        for i, wing_movement in enumerate(self.wing_movements):
            # Extract this Wing's deformation, or None.
            this_deformation = None
            if wing_deformation_angles_ixyz is not None:
                this_deformation = wing_deformation_angles_ixyz[i]

            if isinstance(
                wing_movement,
                aeroelastic_wing_movement_mod.AeroelasticWingMovement,
            ):
                these_wings.append(
                    wing_movement.generate_wing_at_time_step(
                        step,
                        delta_time,
                        deformation_angles_ixyz=this_deformation,
                    )
                )
            else:
                these_wings.append(
                    wing_movement.generate_wing_at_time_step(step, delta_time)
                )

        return geometry.airplane.Airplane(
            wings=these_wings,
            name=self._base_airplane.name,
            Cg_GP1_CgP1=thisCg_GP1_CgP1,
            weight=self._base_airplane.weight,
        )
