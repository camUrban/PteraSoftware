"""Contains the function behind the deprecated modules at the old public paths, the
decorator behind the deprecated parameter kinds, and the source scan that finds the
package's deprecated members and parameters."""

from __future__ import annotations

import ast
import functools
import importlib
import inspect
import warnings
from collections.abc import Callable, Iterator
from pathlib import Path
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


def is_deprecation_warning(node: ast.AST) -> bool:
    """Returns whether a node is a statement calling warnings.warn with
    DeprecationWarning.

    :param node: The node to check.
    :return: True if the node is such a statement and False otherwise.
    """
    if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
        return False
    call = node.value
    if not (
        isinstance(call.func, ast.Attribute)
        and call.func.attr == "warn"
        and isinstance(call.func.value, ast.Name)
        and call.func.value.id == "warnings"
    ):
        return False
    category: ast.expr | None = call.args[1] if len(call.args) > 1 else None
    for keyword in call.keywords:
        if keyword.arg == "category":
            category = keyword.value
    return isinstance(category, ast.Name) and category.id == "DeprecationWarning"


def iter_functions(package_dir: Path) -> Iterator[tuple[str, ast.FunctionDef]]:
    """Yields each module level function and method in a package, read from its source.

    Each function is yielded with its fully qualified name, which is the module's dotted
    name followed by the function's qualified name, as in the function's own __module__
    and __qualname__. A constructor is named by its class instead, since the API
    reference renders its parameters on the class signature.

    :param package_dir: The package's directory.
    :return: An iterator of tuples, each holding a function's fully qualified name and
        its ast.FunctionDef node.
    """
    for path in sorted(package_dir.rglob("*.py")):
        parts = path.relative_to(package_dir).with_suffix("").parts
        if parts[-1] == "__init__":
            parts = parts[:-1]
        module = ".".join((package_dir.name, *parts))
        for node in ast.parse(path.read_text()).body:
            if isinstance(node, ast.FunctionDef):
                yield f"{module}.{node.name}", node
            elif isinstance(node, ast.ClassDef):
                class_id = f"{module}.{node.name}"
                for member in node.body:
                    if isinstance(member, ast.FunctionDef):
                        if member.name == "__init__":
                            yield class_id, member
                        else:
                            yield f"{class_id}.{member.name}", member


def get_parameter_names(function: ast.FunctionDef) -> set[str]:
    """Returns the names of a function's parameters, other than its * and ** parameters.

    :param function: The function's ast.FunctionDef node.
    :return: The set of the parameter names.
    """
    return {
        argument.arg
        for argument in (
            *function.args.posonlyargs,
            *function.args.args,
            *function.args.kwonlyargs,
        )
    }


def find_deprecations(package_dir: Path) -> tuple[set[str], dict[str, list[str]]]:
    """Finds the deprecated members and parameters in a package from its source.

    A deprecation is whatever emits a DeprecationWarning, so this reads the package's
    source rather than any docstring marker. A function, method, or property is
    deprecated when its body calls warnings.warn with DeprecationWarning as a top level
    statement. A parameter is deprecated when such a call sits inside a top level if
    statement whose test names that parameter. The API reference and the public-API test
    both use this scan, so they agree on what is deprecated.

    :param package_dir: The package's directory.
    :return: A tuple of the fully qualified names of the deprecated members, and a dict
        mapping the fully qualified name of each function or method with deprecated
        parameters to those parameters' names, in that order. A constructor's parameters
        are keyed by its class, since the API reference renders them on the class
        signature.
    """
    members: set[str] = set()
    parameters: dict[str, list[str]] = {}
    for function_id, function in iter_functions(package_dir):
        if any(is_deprecation_warning(statement) for statement in function.body):
            members.add(function_id)
        parameter_names = get_parameter_names(function)
        for statement in function.body:
            if not isinstance(statement, ast.If):
                continue
            if not any(is_deprecation_warning(node) for node in ast.walk(statement)):
                continue
            for node in ast.walk(statement.test):
                if isinstance(node, ast.Name) and node.id in parameter_names:
                    parameters.setdefault(function_id, []).append(node.id)
    return members, parameters
