"""Contains the function behind the deprecated modules at the old public paths."""

from __future__ import annotations

import importlib
import warnings
from collections.abc import Callable
from typing import Any


def make_deprecated_module(
    old_module_name: str, new_module_name: str, names: tuple[str, ...]
) -> tuple[Callable[[str], Any], Callable[[], list[str]], list[str]]:
    """Returns the module __getattr__, the module __dir__, and the __all__ list for a
    deprecated module at an old public module path.

    Each deprecated module binds the three to its own dunder names in one statement. The
    __getattr__ reads a forwarded name from its new module and warns that the old path
    is deprecated, so the warning is emitted once per access of a public name, and
    accessing the deprecated module itself never warns. The name is read from its new
    module on every access rather than cached on the deprecated module, so every access
    warns. The __dir__ and the __all__ list hold the forwarded names, so dir() and a
    star import through the deprecated module see exactly those names.

    :param old_module_name: The deprecated module's dotted name.
    :param new_module_name: The dotted name of the module that now defines the names.
    :param names: The public names the deprecated module forwards.
    :return: A tuple of the module __getattr__, the module __dir__, and the __all__
        list, in that order.
    """

    def __getattr__(name: str) -> Any:
        if name not in names:
            raise AttributeError(
                f'module "{old_module_name}" has no attribute "{name}"'
            )
        warnings.warn(
            f"{old_module_name}.{name} is deprecated and will be removed in v6.0.0. Use "
            f"pterasoftware.{name} instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(importlib.import_module(new_module_name), name)

    def __dir__() -> list[str]:
        return list(names)

    return __getattr__, __dir__, list(names)
