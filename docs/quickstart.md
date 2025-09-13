# Quickstart

Below is a minimal steady-solver example after installing the package.

```python
import pterasoftware as ps

airplane = ps.geometry.Airplane(
    wings=[
        ps.geometry.Wing(
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    airfoil=ps.geometry.Airfoil(name="naca2412"),
                ),
                ps.geometry.WingCrossSection(
                    y_le=5.0,
                    airfoil=ps.geometry.Airfoil(name="naca2412"),
                ),
            ],
        ),
    ],
)

operating_point = ps.operating_point.OperatingPoint()

problem = ps.problems.SteadyProblem(
    airplanes=[airplane], operating_point=operating_point
)

solver = (
    ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=problem
    )
)
result = solver.run()
print(result)
```

Next: explore the [Examples](examples.md) and the
[API Reference](api.md).

