"""Tests for the deprecation module and the deprecated modules it backs."""

import importlib
import pkgutil
import types
import unittest
import warnings

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


class TestGetDeprecatedAttribute(unittest.TestCase):
    """Tests for the get_deprecated_attribute function."""

    def test_returns_the_object_from_the_new_module(self) -> None:
        """get_deprecated_attribute should return the named object from the new
        module."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = _deprecation.get_deprecated_attribute(
                "pterasoftware.trim",
                "pterasoftware._trim",
                ("analyze_steady_trim",),
                "analyze_steady_trim",
            )
        self.assertIs(result, _trim.analyze_steady_trim)

    def test_warns_naming_the_old_path_and_the_flat_replacement(self) -> None:
        """get_deprecated_attribute should emit a DeprecationWarning that names the old
        path and the package top level replacement."""
        with self.assertWarns(DeprecationWarning) as context:
            _deprecation.get_deprecated_attribute(
                "pterasoftware.trim",
                "pterasoftware._trim",
                ("analyze_steady_trim",),
                "analyze_steady_trim",
            )
        message = str(context.warning)
        self.assertIn("pterasoftware.trim.analyze_steady_trim", message)
        self.assertIn("v6.0.0", message)
        self.assertIn("pterasoftware.analyze_steady_trim", message)

    def test_raises_for_a_name_not_forwarded(self) -> None:
        """get_deprecated_attribute should raise AttributeError, without warning, for a
        name the deprecated module does not forward."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with self.assertRaises(AttributeError) as context:
                _deprecation.get_deprecated_attribute(
                    "pterasoftware.trim",
                    "pterasoftware._trim",
                    ("analyze_steady_trim",),
                    "SEED",
                )
        self.assertEqual(caught, [])
        self.assertIn("pterasoftware.trim", str(context.exception))
        self.assertIn("SEED", str(context.exception))

    def test_attributes_the_warning_to_the_accessing_line(self) -> None:
        """Accessing a name through a deprecated module should attribute the warning to
        the line that accessed it, not to the deprecated module or this function."""
        trim_module = importlib.import_module("pterasoftware.trim")
        with self.assertWarns(DeprecationWarning) as context:
            trim_module.analyze_steady_trim
        self.assertEqual(context.filename, __file__)


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
            for name in module.NAMES:
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
                self.assertEqual(dir(module), sorted(module.NAMES))

    def test_forwarded_names_cover_every_public_name_with_an_old_path(self) -> None:
        """Together, the deprecated modules should forward every public name that had an
        old public module path, each exactly once."""
        forwarded = [name for module in self.shim_modules for name in module.NAMES]
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
