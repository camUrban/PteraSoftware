"""This module contains tests for the pterasoftware package __init__.py."""

import ast
import importlib
import unittest
from pathlib import Path

import pterasoftware as ps


class TestLazyModuleImports(unittest.TestCase):
    """Tests for lazy module imports via __getattr__."""

    def test_lazy_module_import_aeroelastic_unsteady_ring(self) -> None:
        """Accessing aeroelastic_unsteady_ring_vortex_lattice_method module should
        trigger lazy import.

        :return: None
        """
        module = ps.aeroelastic_unsteady_ring_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__,
            "pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method",
        )

    def test_lazy_module_import_convergence(self) -> None:
        """Accessing convergence module should trigger lazy import.

        :return: None
        """
        module = ps.convergence
        self.assertIsNotNone(module)
        self.assertEqual(module.__name__, "pterasoftware.convergence")

    def test_lazy_module_import_output(self) -> None:
        """Accessing output module should trigger lazy import.

        :return: None
        """
        module = ps.output
        self.assertIsNotNone(module)
        self.assertEqual(module.__name__, "pterasoftware.output")

    def test_lazy_module_import_steady_horseshoe(self) -> None:
        """Accessing steady_horseshoe_vortex_lattice_method module should trigger lazy
        import.

        :return: None
        """
        module = ps.steady_horseshoe_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__, "pterasoftware.steady_horseshoe_vortex_lattice_method"
        )

    def test_lazy_module_import_steady_ring(self) -> None:
        """Accessing steady_ring_vortex_lattice_method module should trigger lazy
        import.

        :return: None
        """
        module = ps.steady_ring_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__, "pterasoftware.steady_ring_vortex_lattice_method"
        )

    def test_lazy_module_import_trim(self) -> None:
        """Accessing trim module should trigger lazy import.

        :return: None
        """
        module = ps.trim
        self.assertIsNotNone(module)
        self.assertEqual(module.__name__, "pterasoftware.trim")

    def test_lazy_module_import_unsteady_ring(self) -> None:
        """Accessing unsteady_ring_vortex_lattice_method module should trigger lazy
        import.

        :return: None
        """
        module = ps.unsteady_ring_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__, "pterasoftware.unsteady_ring_vortex_lattice_method"
        )

    def test_lazy_module_caching(self) -> None:
        """Lazy modules should be cached in globals after first access.

        :return: None
        """
        # First access triggers the import
        module1 = ps.convergence

        # Second access should return the cached version
        module2 = ps.convergence

        self.assertIs(module1, module2)


class TestLazyCallableImports(unittest.TestCase):
    """Tests for lazy callable imports via __getattr__."""

    def test_lazy_callable_import_set_up_logging(self) -> None:
        """Accessing set_up_logging should trigger lazy callable import.

        :return: None
        """
        func = ps.set_up_logging
        self.assertIsNotNone(func)
        self.assertTrue(callable(func))
        self.assertEqual(func.__name__, "set_up_logging")

    def test_lazy_callable_import_save(self) -> None:
        """Accessing save should trigger lazy callable import.

        :return: None
        """
        func = ps.save
        self.assertIsNotNone(func)
        self.assertTrue(callable(func))
        self.assertEqual(func.__name__, "save")

    def test_lazy_callable_import_load(self) -> None:
        """Accessing load should trigger lazy callable import.

        :return: None
        """
        func = ps.load
        self.assertIsNotNone(func)
        self.assertTrue(callable(func))
        self.assertEqual(func.__name__, "load")

    def test_lazy_callable_caching(self) -> None:
        """Lazy callables should be cached in globals after first access.

        :return: None
        """
        # First access triggers the import
        func1 = ps.set_up_logging

        # Second access should return the cached version
        func2 = ps.set_up_logging

        self.assertIs(func1, func2)

    def test_lazy_callable_is_correct_function(self) -> None:
        """The lazy imported set_up_logging should be the actual function from _logging.

        :return: None
        """
        logging_module = importlib.import_module("pterasoftware._logging")

        lazy_func = ps.set_up_logging
        direct_func = logging_module.set_up_logging

        self.assertIs(lazy_func, direct_func)

    def test_lazy_callable_save_is_correct_function(self) -> None:
        """The lazy imported save should be the actual function from _serialization.

        :return: None
        """
        serialization_module = importlib.import_module("pterasoftware._serialization")

        self.assertIs(ps.save, serialization_module.save)

    def test_lazy_callable_load_is_correct_function(self) -> None:
        """The lazy imported load should be the actual function from _serialization.

        :return: None
        """
        serialization_module = importlib.import_module("pterasoftware._serialization")

        self.assertIs(ps.load, serialization_module.load)


class TestEagerImports(unittest.TestCase):
    """Tests for eagerly imported modules."""

    def test_geometry_is_eagerly_imported(self) -> None:
        """The geometry subpackage should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.geometry)
        self.assertEqual(ps.geometry.__name__, "pterasoftware.geometry")

    def test_movements_is_eagerly_imported(self) -> None:
        """The movements subpackage should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.movements)
        self.assertEqual(ps.movements.__name__, "pterasoftware.movements")

    def test_operating_point_is_eagerly_imported(self) -> None:
        """The operating_point module should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.operating_point)
        self.assertEqual(ps.operating_point.__name__, "pterasoftware.operating_point")

    def test_problems_is_eagerly_imported(self) -> None:
        """The problems module should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.problems)
        self.assertEqual(ps.problems.__name__, "pterasoftware.problems")


class TestDirFunction(unittest.TestCase):
    """Tests for the __dir__ function."""

    def test_dir_includes_lazy_modules(self) -> None:
        """The dir() function should include lazy module names.

        :return: None
        """
        package_dir = dir(ps)

        # Check lazy modules are listed
        self.assertIn("convergence", package_dir)
        self.assertIn("output", package_dir)
        self.assertIn("steady_horseshoe_vortex_lattice_method", package_dir)
        self.assertIn("steady_ring_vortex_lattice_method", package_dir)
        self.assertIn("trim", package_dir)
        self.assertIn("unsteady_ring_vortex_lattice_method", package_dir)

    def test_dir_includes_lazy_callables(self) -> None:
        """The dir() function should include lazy callable names.

        :return: None
        """
        package_dir = dir(ps)

        self.assertIn("load", package_dir)
        self.assertIn("save", package_dir)
        self.assertIn("set_up_logging", package_dir)

    def test_dir_includes_eager_imports(self) -> None:
        """The dir() function should include eagerly imported module names.

        :return: None
        """
        package_dir = dir(ps)

        self.assertIn("geometry", package_dir)
        self.assertIn("movements", package_dir)
        self.assertIn("operating_point", package_dir)
        self.assertIn("problems", package_dir)


class TestInvalidAttributeAccess(unittest.TestCase):
    """Tests for accessing invalid attributes."""

    def test_invalid_attribute_raises_attribute_error(self) -> None:
        """Accessing a nonexistent attribute should raise AttributeError.

        :return: None
        """
        with self.assertRaises(AttributeError) as context:
            _ = ps.nonexistent_module

        self.assertIn("nonexistent_module", str(context.exception))
        self.assertIn("has no attribute", str(context.exception))


class TestTypeCheckingImportSync(unittest.TestCase):
    """Tests that the package's TYPE_CHECKING imports stay in sync with its lazy import
    tables.

    Type checkers never execute __getattr__, so they resolve each lazily loaded name
    through the static imports in the package __init__.py's TYPE_CHECKING block instead.
    A lazy name missing from that block is silently typed as Any, which reverts every
    use of it through the package namespace to being unchecked. These tests parse the
    __init__.py source and fail whenever the block and the lazy tables drift apart, in
    either direction.
    """

    type_checking_modules: set[str]
    type_checking_callables: dict[str, tuple[str, str]]

    @classmethod
    def setUpClass(cls) -> None:
        """Parse the TYPE_CHECKING block's imports from the package __init__.py."""
        init_path = Path(ps.__file__)
        tree = ast.parse(init_path.read_text())

        cls.type_checking_modules = set()
        cls.type_checking_callables = {}
        for node in tree.body:
            if not isinstance(node, ast.If):
                continue
            if not (
                isinstance(node.test, ast.Name) and node.test.id == "TYPE_CHECKING"
            ):
                continue
            for statement in node.body:
                if not isinstance(statement, ast.ImportFrom):
                    continue
                if statement.module == "pterasoftware":
                    for alias in statement.names:
                        cls.type_checking_modules.add(alias.name)
                else:
                    for alias in statement.names:
                        cls.type_checking_callables[alias.name] = (
                            str(statement.module),
                            alias.name,
                        )

    def test_every_lazy_module_has_a_type_checking_import(self) -> None:
        """Test that the TYPE_CHECKING block imports exactly the lazy modules.

        A lazy module missing from the block is typed as Any by type checkers, so every
        use of it through the package namespace goes unchecked. An extra import in the
        block advertises a name that __getattr__ cannot deliver.
        """
        self.assertEqual(
            self.type_checking_modules,
            set(ps._LAZY_MODULES),
            msg=(
                "The package __init__.py's TYPE_CHECKING block and its "
                "_LAZY_MODULES table no longer import the same module names. Add "
                "any new lazy module to both, so that type checkers can resolve it."
            ),
        )

    def test_every_lazy_callable_has_a_type_checking_import(self) -> None:
        """Test that the TYPE_CHECKING block imports exactly the lazy callables, each
        from the same module that __getattr__ loads it from."""
        self.assertEqual(
            self.type_checking_callables,
            dict(ps._LAZY_CALLABLES),
            msg=(
                "The package __init__.py's TYPE_CHECKING block and its "
                "_LAZY_CALLABLES table no longer import the same callables from "
                "the same modules. Add any new lazy callable to both, so that "
                "type checkers can resolve it."
            ),
        )
