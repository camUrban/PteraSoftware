"""Contains the registration pattern that grants cross-module access to private
attributes, currently a FreeFlightUnsteadyProblem's MuJoCoModel for the rendering
layer."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import _mujoco_model, problems

# The getter that problems.py registers at import time. It is module private, so no
# other module reads or writes it directly.
_mujoco_model_getter: (
    Callable[[problems.FreeFlightUnsteadyProblem], _mujoco_model.MuJoCoModel] | None
) = None


def register_mujoco_model_getter(
    getter: Callable[[problems.FreeFlightUnsteadyProblem], _mujoco_model.MuJoCoModel],
) -> None:
    """Registers the getter that maps a FreeFlightUnsteadyProblem to its MuJoCoModel.

    problems.py calls this once at import time with a getter defined there, which keeps
    the read of the FreeFlightUnsteadyProblem's private slot inside the module that owns
    it while letting the rendering layer reach the MuJoCoModel through get_mujoco_model.

    :param getter: A callable that takes a FreeFlightUnsteadyProblem and returns its
        MuJoCoModel.
    :return: None
    """
    global _mujoco_model_getter
    _mujoco_model_getter = getter


def get_mujoco_model(
    problem: problems.FreeFlightUnsteadyProblem,
) -> _mujoco_model.MuJoCoModel:
    """Returns a FreeFlightUnsteadyProblem's MuJoCoModel.

    :param problem: The FreeFlightUnsteadyProblem whose MuJoCoModel will be returned.
    :return: The FreeFlightUnsteadyProblem's MuJoCoModel.
    """
    # Constructing the FreeFlightUnsteadyProblem passed in requires importing
    # problems.py, and that import registers the getter, so an unregistered getter is a
    # programming error rather than a reachable state.
    assert _mujoco_model_getter is not None
    return _mujoco_model_getter(problem)
