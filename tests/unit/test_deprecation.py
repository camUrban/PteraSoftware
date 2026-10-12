"""Tests for the deprecation module and the deprecated modules it backs."""

import importlib
import inspect
import pkgutil
import types
import unittest
import warnings
from typing import Any

import pterasoftware as ps
from pterasoftware import _deprecation, _trim

# The deprecated packages at the old public paths. Each holds deprecated modules.
SHIM_PACKAGE_NAMES = ("pterasoftware.geometry", "pterasoftware.movements")

# The public names that never had an old public module path, so no deprecated module
# forwards them.
NAMES_WITHOUT_OLD_PATHS = {"load", "save", "set_up_logging"}


def import_shim_modules() -> list[types.ModuleType]:
    """Imports every deprecated module at the old public paths.

    :return: A list of the deprecated modules, excluding the two deprecated packages.
    """
    shim_modules = []
    for module_name in ps._LAZY_MODULES.values():
        if module_name in SHIM_PACKAGE_NAMES:
            package = importlib.import_module(module_name)
            for info in pkgutil.iter_modules(package.__path__):
                shim_modules.append(
                    importlib.import_module(f"{module_name}.{info.name}")
                )
        else:
            shim_modules.append(importlib.import_module(module_name))
    return shim_modules


@_deprecation.keyword_only(positional_only=("path", "obj"))
def sample_function(
    path: str, obj: object, quality: int = 75, show: bool = True
) -> tuple[str, object, int, bool]:
    """Returns its arguments, so a test can check how a call bound them."""
    return path, obj, quality, show


class SampleClass:
    """A class whose constructor, method, and class method have deprecated parameter
    kinds."""

    @_deprecation.keyword_only()
    def __init__(self, span: float, chord: float = 1.0) -> None:
        """The initialization method."""
        self.span = span
        self.chord = chord

    @_deprecation.keyword_only()
    def scaled_span(self, factor: float) -> float:
        """Returns the span multiplied by a factor."""
        return self.span * factor

    @classmethod
    @_deprecation.keyword_only()
    def from_area(cls, area: float, chord: float = 1.0) -> "SampleClass":
        """Returns a SampleClass with the given area and chord."""
        return cls(span=area / chord, chord=chord)


class TestMakeDeprecatedModule(unittest.TestCase):
    """Tests for the make_deprecated_module function."""

    def setUp(self) -> None:
        """Make the module __getattr__, module __dir__, and __all__ list for a
        deprecated module that forwards one name."""
        self.module_getattr, self.module_dir, self.module_all = (
            _deprecation.make_deprecated_module(
                "pterasoftware.trim", "pterasoftware._trim", ("analyze_steady_trim",)
            )
        )

    def test_getattr_returns_the_object_from_the_new_module(self) -> None:
        """The module __getattr__ should return the named object from the new module."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = self.module_getattr("analyze_steady_trim")
        self.assertIs(result, _trim.analyze_steady_trim)

    def test_getattr_warns_naming_the_old_path_and_the_flat_replacement(self) -> None:
        """The module __getattr__ should emit a DeprecationWarning that names the old
        path and the package top level replacement."""
        with self.assertWarns(DeprecationWarning) as context:
            self.module_getattr("analyze_steady_trim")
        message = str(context.warning)
        self.assertIn("pterasoftware.trim.analyze_steady_trim", message)
        self.assertIn("v6.0.0", message)
        self.assertIn("pterasoftware.analyze_steady_trim", message)

    def test_getattr_raises_for_a_name_not_forwarded(self) -> None:
        """The module __getattr__ should raise AttributeError, without warning, for a
        name the deprecated module does not forward."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with self.assertRaises(AttributeError) as context:
                self.module_getattr("SEED")
        self.assertEqual(caught, [])
        self.assertIn("pterasoftware.trim", str(context.exception))
        self.assertIn("SEED", str(context.exception))

    def test_dir_and_all_list_the_forwarded_names(self) -> None:
        """The module __dir__ and the __all__ list should both hold exactly the
        forwarded names."""
        self.assertEqual(self.module_dir(), ["analyze_steady_trim"])
        self.assertEqual(self.module_all, ["analyze_steady_trim"])

    def test_attributes_the_warning_to_the_accessing_line(self) -> None:
        """Accessing a name through a deprecated module should attribute the warning to
        the line that accessed it, not to the deprecated module or this function."""
        trim_module = importlib.import_module("pterasoftware.trim")
        with self.assertWarns(DeprecationWarning) as context:
            trim_module.analyze_steady_trim
        self.assertEqual(context.filename, __file__)


class TestKeywordOnly(unittest.TestCase):
    """Tests for the keyword_only decorator."""

    def test_signature_reports_the_v6_parameter_kinds(self) -> None:
        """The decorated function's signature should mark the positional_only parameters
        positional-only and every other parameter keyword-only, while keeping its name
        and docstring."""
        parameters = inspect.signature(sample_function).parameters
        self.assertEqual(
            [parameter.kind for parameter in parameters.values()],
            [
                inspect.Parameter.POSITIONAL_ONLY,
                inspect.Parameter.POSITIONAL_ONLY,
                inspect.Parameter.KEYWORD_ONLY,
                inspect.Parameter.KEYWORD_ONLY,
            ],
        )
        self.assertEqual(list(parameters), ["path", "obj", "quality", "show"])
        self.assertEqual(sample_function.__name__, "sample_function")
        self.assertEqual(
            sample_function.__doc__,
            "Returns its arguments, so a test can check how a call bound them.",
        )

    def test_signatures_of_methods_leave_self_and_cls_alone(self) -> None:
        """The signatures of a decorated constructor, method, and class method should
        mark every parameter keyword-only except self and cls."""
        self.assertEqual(
            str(inspect.signature(SampleClass)),
            "(*, span: float, chord: float = 1.0) -> None",
        )
        self.assertEqual(
            str(inspect.signature(SampleClass.scaled_span)),
            "(self, *, factor: float) -> float",
        )
        self.assertEqual(
            str(inspect.signature(SampleClass.from_area)),
            "(*, area: float, chord: float = 1.0) -> 'SampleClass'",
        )

    def test_calls_with_the_v6_parameter_kinds_do_not_warn(self) -> None:
        """Calls that pass every parameter with its v6.0.0 kind should not warn and
        should bind the arguments unchanged."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            function_result = sample_function("a", 1, quality=50)
            instance = SampleClass(span=2.0)
            method_result = instance.scaled_span(factor=3.0)
            class_method_result = SampleClass.from_area(area=6.0, chord=2.0)
        self.assertEqual(caught, [])
        self.assertEqual(function_result, ("a", 1, 50, True))
        self.assertEqual(method_result, 6.0)
        self.assertEqual(class_method_result.span, 3.0)

    def test_passing_parameters_positionally_warns_once_and_binds(self) -> None:
        """Passing parameters outside positional_only positionally should emit one
        DeprecationWarning naming all of them, and still bind them."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = sample_function("a", 1, 50, False)
        self.assertEqual(len(caught), 1)
        self.assertIs(caught[0].category, DeprecationWarning)
        message = str(caught[0].message)
        self.assertIn("quality, show positionally to sample_function", message)
        self.assertIn("v6.0.0", message)
        self.assertEqual(result, ("a", 1, 50, False))

    def test_passing_positional_only_parameters_by_keyword_warns_and_binds(
        self,
    ) -> None:
        """Passing positional_only parameters by keyword should emit one
        DeprecationWarning naming all of them, and still bind them."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = sample_function(path="a", obj=1)
        self.assertEqual(len(caught), 1)
        self.assertIs(caught[0].category, DeprecationWarning)
        message = str(caught[0].message)
        self.assertIn("path, obj by keyword to sample_function", message)
        self.assertIn("v6.0.0", message)
        self.assertEqual(result, ("a", 1, 75, True))

    def test_constructor_warning_names_the_class(self) -> None:
        """A decorated constructor's warning should name its class rather than its
        __init__ method."""
        with self.assertWarns(DeprecationWarning) as context:
            instance = SampleClass(2.0)
        message = str(context.warning)
        self.assertIn("span positionally to SampleClass is deprecated", message)
        self.assertNotIn("__init__", message)
        self.assertEqual(instance.span, 2.0)

    def test_method_and_class_method_warnings_skip_self_and_cls(self) -> None:
        """Calling a decorated method or class method positionally should warn about its
        own parameters and never about self or cls."""
        instance = SampleClass(span=2.0)
        with self.assertWarns(DeprecationWarning) as method_context:
            method_result = instance.scaled_span(3.0)
        with self.assertWarns(DeprecationWarning) as class_method_context:
            class_method_result = SampleClass.from_area(6.0)
        method_message = str(method_context.warning)
        class_method_message = str(class_method_context.warning)
        self.assertIn("Passing factor positionally", method_message)
        self.assertNotIn("self", method_message)
        self.assertIn("Passing area positionally", class_method_message)
        self.assertNotIn("cls", class_method_message)
        self.assertEqual(method_result, 6.0)
        self.assertEqual(class_method_result.span, 6.0)

    def test_attributes_the_warning_to_the_calling_line(self) -> None:
        """A deprecated call should attribute the warning to the line that made the
        call, not to the decorator."""
        with self.assertWarns(DeprecationWarning) as function_context:
            sample_function(path="a", obj=1)
        with self.assertWarns(DeprecationWarning) as constructor_context:
            SampleClass(2.0)
        self.assertEqual(function_context.filename, __file__)
        self.assertEqual(constructor_context.filename, __file__)

    def test_rejects_a_name_that_is_not_a_parameter(self) -> None:
        """Decorating a callable should raise ValueError when positional_only names a
        parameter the callable does not have."""

        # noinspection unused-parameter
        def target(path: str) -> None:
            pass

        with self.assertRaises(ValueError):
            _deprecation.keyword_only(positional_only=("missing",))(target)

    def test_rejects_a_name_that_is_not_a_leading_parameter(self) -> None:
        """Decorating a callable should raise ValueError when positional_only names a
        parameter that follows one that becomes keyword-only."""

        # noinspection unused-parameter
        def target(path: str, show: bool = True) -> None:
            pass

        with self.assertRaises(ValueError):
            _deprecation.keyword_only(positional_only=("show",))(target)

    def test_rejects_a_callable_with_a_variadic_parameter(self) -> None:
        """Decorating a callable with a variadic parameter should raise TypeError."""

        # noinspection unused-parameter
        def target(path: str, *args: int) -> None:
            pass

        with self.assertRaises(TypeError):
            _deprecation.keyword_only()(target)


class TestDeprecatedModules(unittest.TestCase):
    """Tests for the deprecated modules at the old public module paths."""

    shim_modules: list[types.ModuleType]

    @classmethod
    def setUpClass(cls) -> None:
        """Import every deprecated module."""
        cls.shim_modules = import_shim_modules()

    def test_importing_a_deprecated_module_does_not_warn(self) -> None:
        """Importing a deprecated module should not warn, since only accessing a public
        name through it does."""
        for module in self.shim_modules:
            with self.subTest(module=module.__name__):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    importlib.reload(module)
                self.assertEqual(caught, [])

    def test_every_forwarded_name_warns_once_and_is_the_public_object(self) -> None:
        """Accessing each name through its deprecated module should emit exactly one
        DeprecationWarning and return the same object as the package top level."""
        for module in self.shim_modules:
            for name in module.__all__:
                with self.subTest(module=module.__name__, name=name):
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        result = getattr(module, name)
                    self.assertEqual(len(caught), 1)
                    self.assertIs(caught[0].category, DeprecationWarning)
                    self.assertIs(result, getattr(ps, name))

    def test_dir_lists_exactly_the_forwarded_names(self) -> None:
        """Calling dir() on a deprecated module should list exactly the names it
        forwards."""
        for module in self.shim_modules:
            with self.subTest(module=module.__name__):
                self.assertEqual(dir(module), sorted(module.__all__))

    def test_star_import_binds_exactly_the_forwarded_names(self) -> None:
        """A star import through a deprecated module should bind exactly the names it
        forwards, each to the public object, and warn once per name."""
        for module in self.shim_modules:
            with self.subTest(module=module.__name__):
                namespace: dict[str, Any] = {}
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    exec(f"from {module.__name__} import *", namespace)
                # The exec call adds __builtins__ to the namespace, and warnings.warn
                # records each warning it attributes to the star import's frame in
                # __warningregistry__ there, so neither dunder name came from the
                # import.
                bound = {
                    name: value
                    for name, value in namespace.items()
                    if not name.startswith("__")
                }
                self.assertEqual(set(bound), set(module.__all__))
                for name, result in bound.items():
                    self.assertIs(result, getattr(ps, name))
                self.assertEqual(len(caught), len(module.__all__))
                for warning in caught:
                    self.assertIs(warning.category, DeprecationWarning)

    def test_forwarded_names_cover_every_public_name_with_an_old_path(self) -> None:
        """Together, the deprecated modules should forward every public name that had an
        old public module path, each exactly once."""
        forwarded = [name for module in self.shim_modules for name in module.__all__]
        self.assertEqual(len(forwarded), len(set(forwarded)))
        self.assertEqual(set(forwarded), set(ps.__all__) - NAMES_WITHOUT_OLD_PATHS)


class TestDeprecatedPackages(unittest.TestCase):
    """Tests for the deprecated packages at the old public paths."""

    def test_packages_expose_their_deprecated_modules(self) -> None:
        """Each deprecated package should expose its deprecated modules as attributes,
        as the old packages did."""
        for package_name in SHIM_PACKAGE_NAMES:
            package = importlib.import_module(package_name)
            for info in pkgutil.iter_modules(package.__path__):
                with self.subTest(module=f"{package_name}.{info.name}"):
                    self.assertIs(
                        getattr(package, info.name),
                        importlib.import_module(f"{package_name}.{info.name}"),
                    )
