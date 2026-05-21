"""Contains the AeroelasticMovement class.

**Contains the following classes:**

AeroelasticMovement: A class used to contain an AeroelasticUnsteadyProblem's movement.

**Contains the following functions:**

None
"""

from __future__ import annotations

from typing import cast

import numpy as np

from .. import _core, geometry
from .. import operating_point as operating_point_mod
from . import aeroelastic_airplane_movement as aeroelastic_airplane_movement_mod
from . import (
    aeroelastic_operating_point_movement as aeroelastic_operating_point_movement_mod,
)


class AeroelasticMovement(_core.CoreMovement):
    """A class used to contain an AeroelasticUnsteadyProblem's movement.

    In aeroelastic simulations, wing geometry is prescribed via oscillation parameters
    (flapping, CG oscillation, etc.) but the solver adds structural deformation at each
    time step based on aerodynamic loads. OperatingPoints are prescribed via the same
    oscillation parameters as OperatingPointMovement.

    AeroelasticMovement pre generates all OperatingPoints upfront (since they are
    prescribed) but does not pre generate Airplanes, because the deformed wing geometry
    at each time step depends on the solver's structural response calculation.

    **Contains the following methods:**

    lcm_period: The least common multiple of all motion periods, ensuring all motions
    complete an integer number of cycles when cycle averaging forces and moments.

    max_period: The longest period of motion of AeroelasticMovement's sub movement
    objects, the motion(s) of its sub sub movement object(s), and the motions of its sub
    sub sub movement objects.

    min_period: The shortest non zero period of motion of AeroelasticMovement's sub
    movement objects, the motion(s) of its sub sub movement object(s), and the motions
    of its sub sub sub movement objects.

    static: Flags if AeroelasticMovement's sub movement objects, its sub sub movement
    object(s), and its sub sub sub movement objects all represent no motion.

    generate_airplane_at_time_step: Creates the Airplane at a single time step, applying
    deformation from the solver's structural response.
    """

    __slots__ = ("_operating_points",)

    def __init__(
        self,
        airplane_movements: list[
            aeroelastic_airplane_movement_mod.AeroelasticAirplaneMovement
        ],
        operating_point_movement: aeroelastic_operating_point_movement_mod.AeroelasticOperatingPointMovement,
        delta_time: float | int,
        num_steps: int,
        max_wake_rows: int | None = None,
    ) -> None:
        """The initialization method.

        This method checks that all Wings maintain their symmetry type across all time
        steps (using the undeformed prescribed geometry). See the WingMovement class
        documentation for more details on this requirement. See the Wing class
        documentation for more information on symmetry types.

        :param airplane_movements: A list of the AeroelasticAirplaneMovements associated
            with each of the AeroelasticUnsteadyProblem's Airplanes.
        :param operating_point_movement: An AeroelasticOperatingPointMovement holding
            the oscillation parameters for prescribing OperatingPoints at each time
            step.
        :param delta_time: The time, in seconds, between each time step. It must be a
            positive number (int or float). It will be converted internally to a float.
        :param num_steps: The number of time steps to simulate. It must be a positive
            int.
        :param max_wake_rows: The maximum number of chordwise wake ring vortex rows per
            Wing. Must be a positive int if set. The default is None (no truncation).
        :return: None
        """
        # Validate that every element is an AeroelasticAirplaneMovement, not
        # just a CoreAirplaneMovement. CoreMovement.__init__() validates at
        # the Core level, but AeroelasticMovement enforces the stricter type.
        for airplane_movement in airplane_movements:
            if not isinstance(
                airplane_movement,
                aeroelastic_airplane_movement_mod.AeroelasticAirplaneMovement,
            ):
                raise TypeError(
                    "Every element in airplane_movements must be an "
                    "AeroelasticAirplaneMovement."
                )

        # Validate that operating_point_movement is an
        # AeroelasticOperatingPointMovement.
        if not isinstance(
            operating_point_movement,
            aeroelastic_operating_point_movement_mod.AeroelasticOperatingPointMovement,
        ):
            raise TypeError(
                "operating_point_movement must be an "
                "AeroelasticOperatingPointMovement."
            )

        # --- Initialize CoreMovement ---
        super().__init__(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=delta_time,
            num_steps=num_steps,
            max_wake_rows=max_wake_rows,
        )

        # --- Batch generate OperatingPoints ---
        # OperatingPoints are prescribed in aeroelastic simulations, so
        # generate them all upfront.
        operating_points_list = operating_point_movement.generate_operating_points(
            num_steps=self._num_steps, delta_time=self._delta_time
        )
        self._operating_points: tuple[operating_point_mod.OperatingPoint, ...] = tuple(
            operating_points_list
        )

    # --- Immutable: read only properties ---
    @property
    def operating_point_movement(
        self,
    ) -> aeroelastic_operating_point_movement_mod.AeroelasticOperatingPointMovement:
        assert isinstance(
            self._operating_point_movement,
            aeroelastic_operating_point_movement_mod.AeroelasticOperatingPointMovement,
        )
        return self._operating_point_movement

    @property
    def airplane_movements(
        self,
    ) -> tuple[aeroelastic_airplane_movement_mod.AeroelasticAirplaneMovement, ...]:
        return cast(
            tuple[aeroelastic_airplane_movement_mod.AeroelasticAirplaneMovement, ...],
            self._airplane_movements,
        )

    @property
    def operating_points(self) -> tuple[operating_point_mod.OperatingPoint, ...]:
        return self._operating_points

    def generate_airplane_at_time_step(
        self,
        airplane_movement_index: int,
        step: int,
        wing_deformation_angles_ixyz: list[np.ndarray | None] | None = None,
    ) -> geometry.airplane.Airplane:
        """Creates the Airplane at a single time step for a given
        AeroelasticAirplaneMovement, applying deformation from the solver's structural
        response.

        This is the method the aeroelastic solver calls at each time step to get the
        deformed Airplane geometry.

        :param airplane_movement_index: The index of the AeroelasticAirplaneMovement in
            this AeroelasticMovement's airplane_movements tuple.
        :param step: The time step index. Must be a non negative int.
        :param wing_deformation_angles_ixyz: A list of (N_wcs, 3) ndarrays of floats,
            one per Wing, where N_wcs is the number of WingCrossSections in that Wing.
            Each row is a (3,) deformation angle vector using an intrinsic xy'z"
            sequence. The units are in degrees. When None, no deformation is applied.
            The default is None.
        :return: The Airplane at this time step, with structural deformation applied if
            provided.
        """
        return self.airplane_movements[
            airplane_movement_index
        ].generate_airplane_at_time_step(
            step,
            self._delta_time,
            wing_deformation_angles_ixyz=wing_deformation_angles_ixyz,
        )
