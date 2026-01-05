"""This module contains tests for the pterasoftware package __init__.py."""

import unittest

import pterasoftware as ps


class TestLazyModuleImports(unittest.TestCase):
    """Tests for lazy module imports via __getattr__."""

    def test_lazy_module_import_convergence(self):
        """Accessing convergence module should trigger lazy import.

        :return: None
        """
        module = ps.convergence
        self.assertIsNotNone(module)
        self.assertEqual(module.__name__, "pterasoftware.convergence")

    def test_lazy_module_import_output(self):
        """Accessing output module should trigger lazy import.

        :return: None
        """
        module = ps.output
        self.assertIsNotNone(module)
        self.assertEqual(module.__name__, "pterasoftware.output")

    def test_lazy_module_import_steady_horseshoe(self):
        """Accessing steady_horseshoe_vortex_lattice_method module should trigger lazy
        import.

        :return: None
        """
        module = ps.steady_horseshoe_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__, "pterasoftware.steady_horseshoe_vortex_lattice_method"
        )

    def test_lazy_module_import_steady_ring(self):
        """Accessing steady_ring_vortex_lattice_method module should trigger lazy
        import.

        :return: None
        """
        module = ps.steady_ring_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__, "pterasoftware.steady_ring_vortex_lattice_method"
        )

    def test_lazy_module_import_trim(self):
        """Accessing trim module should trigger lazy import.

        :return: None
        """
        module = ps.trim
        self.assertIsNotNone(module)
        self.assertEqual(module.__name__, "pterasoftware.trim")

    def test_lazy_module_import_unsteady_ring(self):
        """Accessing unsteady_ring_vortex_lattice_method module should trigger lazy
        import.

        :return: None
        """
        module = ps.unsteady_ring_vortex_lattice_method
        self.assertIsNotNone(module)
        self.assertEqual(
            module.__name__, "pterasoftware.unsteady_ring_vortex_lattice_method"
        )

    def test_lazy_module_caching(self):
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

    def test_lazy_callable_import_set_up_logging(self):
        """Accessing set_up_logging should trigger lazy callable import.

        :return: None
        """
        func = ps.set_up_logging
        self.assertIsNotNone(func)
        self.assertTrue(callable(func))
        self.assertEqual(func.__name__, "set_up_logging")

    def test_lazy_callable_caching(self):
        """Lazy callables should be cached in globals after first access.

        :return: None
        """
        # First access triggers the import
        func1 = ps.set_up_logging

        # Second access should return the cached version
        func2 = ps.set_up_logging

        self.assertIs(func1, func2)

    def test_lazy_callable_is_correct_function(self):
        """The lazy imported set_up_logging should be the actual function from
        _logging.

        :return: None
        """
        # noinspection PyProtectedMember
        from pterasoftware import _logging

        lazy_func = ps.set_up_logging
        direct_func = _logging.set_up_logging

        self.assertIs(lazy_func, direct_func)


class TestEagerImports(unittest.TestCase):
    """Tests for eagerly imported modules."""

    def test_geometry_is_eagerly_imported(self):
        """The geometry subpackage should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.geometry)
        self.assertEqual(ps.geometry.__name__, "pterasoftware.geometry")

    def test_movements_is_eagerly_imported(self):
        """The movements subpackage should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.movements)
        self.assertEqual(ps.movements.__name__, "pterasoftware.movements")

    def test_operating_point_is_eagerly_imported(self):
        """The operating_point module should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.operating_point)
        self.assertEqual(ps.operating_point.__name__, "pterasoftware.operating_point")

    def test_problems_is_eagerly_imported(self):
        """The problems module should be eagerly imported.

        :return: None
        """
        self.assertIsNotNone(ps.problems)
        self.assertEqual(ps.problems.__name__, "pterasoftware.problems")


class TestDirFunction(unittest.TestCase):
    """Tests for the __dir__ function."""

    def test_dir_includes_lazy_modules(self):
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

    def test_dir_includes_lazy_callables(self):
        """The dir() function should include lazy callable names.

        :return: None
        """
        package_dir = dir(ps)

        self.assertIn("set_up_logging", package_dir)

    def test_dir_includes_eager_imports(self):
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

    def test_invalid_attribute_raises_attribute_error(self):
        """Accessing a nonexistent attribute should raise AttributeError.

        :return: None
        """
        with self.assertRaises(AttributeError) as context:
            _ = ps.nonexistent_module

        self.assertIn("nonexistent_module", str(context.exception))
        self.assertIn("has no attribute", str(context.exception))


if __name__ == "__main__":
    unittest.main()
