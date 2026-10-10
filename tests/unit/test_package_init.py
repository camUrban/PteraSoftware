"""This module contains tests for the pterasoftware package __init__.py."""

import ast
import concurrent.futures
import importlib
import os
import subprocess
import sys
import unittest
import warnings
from pathlib import Path

import pterasoftware as ps


class TestAll(unittest.TestCase):
    """Tests for the package's __all__."""

    def test_all_has_no_duplicates(self) -> None:
        """Test that no name appears in __all__ more than once.

        :return: None
        """
        self.assertEqual(len(ps.__all__), len(set(ps.__all__)))

    def test_all_matches_the_lazy_callables(self) -> None:
        """Test that __all__ lists exactly the names in the lazy callable table.

        :return: None
        """
        self.assertEqual(set(ps.__all__), set(ps._LAZY_CALLABLES))


class TestLazyCallableImports(unittest.TestCase):
    """Tests for lazy imports of the public names via __getattr__."""

    def test_every_public_name_resolves_to_its_internal_definition(self) -> None:
        """Test that every name in __all__ resolves to the object its internal module
        defines.

        :return: None
        """
        for name in ps.__all__:
            with self.subTest(name=name):
                module_path, attr_name = ps._LAZY_CALLABLES[name]
                module = importlib.import_module(module_path)
                self.assertIs(getattr(ps, name), getattr(module, attr_name))

    def test_every_public_name_keeps_its_own_name(self) -> None:
        """Test that every name in __all__ refers to an object defined under that same
        name.

        :return: None
        """
        for name in ps.__all__:
            with self.subTest(name=name):
                self.assertEqual(getattr(ps, name).__name__, name)

    def test_lazy_callable_caching(self) -> None:
        """Lazy callables should be cached in globals after first access.

        :return: None
        """
        # First access triggers the import
        func1 = ps.set_up_logging

        # Second access should return the cached version
        func2 = ps.set_up_logging

        self.assertIs(func1, func2)


class TestFirstAccess(unittest.TestCase):
    """Tests that every public name loads when it is the first one accessed.

    The package imports nothing eagerly, so the first public name a program accesses
    decides which internal module starts loading, and a circular import between internal
    modules can fail for one starting point while succeeding for every other. Within one
    test process the internal modules are already loaded, which would hide such a
    failure, so each access runs in its own fresh interpreter.
    """

    def test_every_public_name_loads_as_the_first_access(self) -> None:
        """Test that accessing each name in __all__ succeeds in a fresh interpreter that
        has accessed nothing else.

        :return: None
        """

        def access_first(name: str) -> subprocess.CompletedProcess[str]:
            """Accesses one public name in a fresh interpreter.

            :param name: The public name to access.
            :return: The completed interpreter process.
            """
            return subprocess.run(
                [sys.executable, "-c", f"import pterasoftware as ps; ps.{name}"],
                capture_output=True,
                text=True,
            )

        # Each interpreter spends about a second importing the package's dependencies,
        # so running them concurrently keeps this test from dominating the suite's run
        # time. Each one also holds a few hundred megabytes once the solver stack is
        # imported, so the pool is capped rather than sized to the core count, which
        # would start up to 32 of them at once on a large workstation.
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(8, os.cpu_count() or 1)
        ) as executor:
            results = dict(zip(ps.__all__, executor.map(access_first, ps.__all__)))

        for name, result in results.items():
            with self.subTest(name=name):
                self.assertEqual(result.returncode, 0, msg=result.stderr)


class TestDeprecatedModuleAttributes(unittest.TestCase):
    """Tests for the deprecated module attributes at the old public module paths."""

    def test_every_deprecated_module_attribute_resolves_to_its_module(self) -> None:
        """Test that every old module attribute resolves to the deprecated module at
        that path.

        :return: None
        """
        for name, module_path in ps._LAZY_MODULES.items():
            with self.subTest(name=name):
                self.assertEqual(getattr(ps, name).__name__, module_path)

    def test_accessing_a_deprecated_module_attribute_does_not_warn(self) -> None:
        """Test that accessing an old module attribute emits no warning, since only
        accessing a public name through it does.

        :return: None
        """
        for name in ps._LAZY_MODULES:
            with self.subTest(name=name):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    getattr(ps, name)
                self.assertEqual(caught, [])


class TestDirFunction(unittest.TestCase):
    """Tests for the __dir__ function."""

    def test_dir_lists_exactly_the_public_names(self) -> None:
        """The dir() function should list exactly the names in __all__.

        :return: None
        """
        self.assertEqual(dir(ps), sorted(ps.__all__))


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
    """Tests that the package's TYPE_CHECKING imports stay in sync with its lazy
    callable table.

    Type checkers never execute __getattr__, so they resolve each lazily loaded name
    through the static imports in the package __init__.py's TYPE_CHECKING block instead.
    A lazy name missing from that block is silently typed as Any, which reverts every
    use of it through the package namespace to being unchecked. This test parses the
    __init__.py source and fails whenever the block and the table drift apart, in either
    direction.
    """

    type_checking_callables: dict[str, tuple[str, str]]

    @classmethod
    def setUpClass(cls) -> None:
        """Parse the TYPE_CHECKING block's imports from the package __init__.py."""
        init_path = Path(ps.__file__)
        tree = ast.parse(init_path.read_text())

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
                for alias in statement.names:
                    cls.type_checking_callables[alias.name] = (
                        str(statement.module),
                        alias.name,
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
