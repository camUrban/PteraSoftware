"""Contains the function behind the deprecated modules at the old public paths and the
decorator behind the deprecated parameter kinds."""

from __future__ import annotations

import functools
import importlib
import inspect
import warnings
from collections.abc import Callable
from typing import Any, ParamSpec, TypeVar

P = ParamSpec("P")
R = TypeVar("R")


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


def keyword_only(
    positional_only: tuple[str, ...] = (),
) -> Callable[[Callable[P, R]], Callable[P, R]]:
    """Returns a decorator that deprecates passing a callable's parameters with any kind
    other than the kind they will have in v6.0.0.

    In v6.0.0, every parameter of the decorated callable becomes keyword-only, except
    the parameters named in positional_only, which become positional-only. Until then,
    the decorated callable still accepts every parameter both ways, and it emits a
    DeprecationWarning when a parameter outside positional_only is passed positionally,
    or when a parameter in positional_only is passed by keyword. A leading parameter
    named self or cls is left alone, so the decorator applies to methods and, placed
    under classmethod, to class methods. The decorated callable reports its v6.0.0
    signature to inspect.signature, so the API reference and the public-API test see the
    v6.0.0 parameter kinds now.

    :param positional_only: The names of the parameters that become positional-only in
        v6.0.0. They must be the decorated callable's leading parameters, in order, so
        they cannot be given for a method. The default is an empty tuple, meaning every
        parameter becomes keyword-only.
    :return: The decorator.
    """

    def decorator(func: Callable[P, R]) -> Callable[P, R]:
        signature = inspect.signature(func)

        # These are the names that can be bound positionally, in order, including any
        # leading self or cls.
        positional_names = [
            name
            for name, parameter in signature.parameters.items()
            if parameter.kind
            in (
                inspect.Parameter.POSITIONAL_ONLY,
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
            )
        ]
        bound_name = (
            positional_names[0]
            if positional_names and positional_names[0] in ("self", "cls")
            else None
        )

        for name in positional_only:
            if name not in signature.parameters:
                raise ValueError(
                    f'{func.__qualname__} has no parameter named "{name}".'
                )
        for parameter in signature.parameters.values():
            if parameter.kind in (
                inspect.Parameter.VAR_POSITIONAL,
                inspect.Parameter.VAR_KEYWORD,
            ):
                raise TypeError(
                    f"{func.__qualname__} has a variadic parameter, so its parameter "
                    "kinds cannot be deprecated."
                )

        # Build the v6.0.0 signature. Signature raises a ValueError if a name in
        # positional_only is not one of the leading parameters, since a positional-only
        # parameter cannot follow a keyword-only one.
        new_parameters = []
        for name, parameter in signature.parameters.items():
            if name == bound_name:
                new_parameters.append(parameter)
            elif name in positional_only:
                new_parameters.append(
                    parameter.replace(kind=inspect.Parameter.POSITIONAL_ONLY)
                )
            else:
                new_parameters.append(
                    parameter.replace(kind=inspect.Parameter.KEYWORD_ONLY)
                )
        new_signature = signature.replace(parameters=new_parameters)

        # Name the callable as its users do, so a constructor is named by its class.
        display_name = func.__qualname__.removesuffix(".__init__")

        @functools.wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            positional = [
                name
                for name in positional_names[: len(args)]
                if name != bound_name and name not in positional_only
            ]
            if positional:
                joined = ", ".join(positional)
                warnings.warn(
                    f"Passing {joined} positionally to {display_name} is deprecated. "
                    f"Starting in v6.0.0, {display_name} accepts {joined} only by "
                    "keyword.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            by_keyword = [name for name in kwargs if name in positional_only]
            if by_keyword:
                joined = ", ".join(by_keyword)
                warnings.warn(
                    f"Passing {joined} by keyword to {display_name} is deprecated. "
                    f"Starting in v6.0.0, {display_name} accepts {joined} only "
                    "positionally.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            return func(*args, **kwargs)

        setattr(wrapper, "__signature__", new_signature)
        return wrapper

    return decorator
