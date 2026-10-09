"""Contains the function behind the deprecated modules at the old public paths."""

from __future__ import annotations

import importlib
import warnings
from typing import Any


def get_deprecated_attribute(
    old_module_name: str,
    new_module_name: str,
    names: tuple[str, ...],
    name: str,
) -> Any:
    """Returns a public name through one of the deprecated modules at the old public
    module paths, and warns that the old path is deprecated.

    Each deprecated module's __getattr__ calls this function, so the warning is emitted
    once per access of a public name, and accessing the deprecated module itself never
    warns. The name is read from its new module on every access rather than cached on
    the deprecated module, so every access warns.

    :param old_module_name: The deprecated module's dotted name.
    :param new_module_name: The dotted name of the module that now defines the names.
    :param names: The public names the deprecated module forwards.
    :param name: The name being accessed.
    :return: The object the name refers to in its new module.
    """
    if name not in names:
        raise AttributeError(f'module "{old_module_name}" has no attribute "{name}"')
    warnings.warn(
        f"{old_module_name}.{name} is deprecated and will be removed in v6.0.0. Use "
        f"pterasoftware.{name} instead.",
        DeprecationWarning,
        stacklevel=3,
    )
    return getattr(importlib.import_module(new_module_name), name)
