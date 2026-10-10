# API Reference

This reference covers every public class and function in Ptera Software. Each one is available at the package top level, for example as `ps.Airplane` after `import pterasoftware as ps`. Pages are generated automatically from docstrings.

## Geometry

```{toctree}
:maxdepth: 1

Airfoil <api/Airfoil>
WingCrossSection <api/WingCrossSection>
Wing <api/Wing>
Airplane <api/Airplane>
```

## Operating Point

```{toctree}
:maxdepth: 1

OperatingPoint <api/OperatingPoint>
```

## Movements

```{toctree}
:maxdepth: 1

WingCrossSectionMovement <api/WingCrossSectionMovement>
WingMovement <api/WingMovement>
AirplaneMovement <api/AirplaneMovement>
OperatingPointMovement <api/OperatingPointMovement>
Movement <api/Movement>
AeroelasticWingCrossSectionMovement <api/AeroelasticWingCrossSectionMovement>
AeroelasticWingMovement <api/AeroelasticWingMovement>
AeroelasticAirplaneMovement <api/AeroelasticAirplaneMovement>
AeroelasticMovement <api/AeroelasticMovement>
FreeFlightOperatingPointMovement <api/FreeFlightOperatingPointMovement>
FreeFlightMovement <api/FreeFlightMovement>
```

## Problems

```{toctree}
:maxdepth: 1

SteadyProblem <api/SteadyProblem>
UnsteadyProblem <api/UnsteadyProblem>
AeroelasticUnsteadyProblem <api/AeroelasticUnsteadyProblem>
FreeFlightUnsteadyProblem <api/FreeFlightUnsteadyProblem>
```

## Solvers

```{toctree}
:maxdepth: 1

SteadyHorseshoeVortexLatticeMethodSolver <api/SteadyHorseshoeVortexLatticeMethodSolver>
SteadyRingVortexLatticeMethodSolver <api/SteadyRingVortexLatticeMethodSolver>
UnsteadyRingVortexLatticeMethodSolver <api/UnsteadyRingVortexLatticeMethodSolver>
AeroelasticUnsteadyRingVortexLatticeMethodSolver <api/AeroelasticUnsteadyRingVortexLatticeMethodSolver>
FreeFlightUnsteadyRingVortexLatticeMethodSolver <api/FreeFlightUnsteadyRingVortexLatticeMethodSolver>
```

## Convergence and Trim

```{toctree}
:maxdepth: 1

analyze_steady_convergence() <api/analyze_steady_convergence>
analyze_unsteady_convergence() <api/analyze_unsteady_convergence>
analyze_steady_trim() <api/analyze_steady_trim>
analyze_unsteady_trim() <api/analyze_unsteady_trim>
```

## Output

```{toctree}
:maxdepth: 1

draw() <api/draw>
animate() <api/animate>
plot_results_versus_time() <api/plot_results_versus_time>
log_results() <api/log_results>
```

## Saving, Loading, and Logging

```{toctree}
:maxdepth: 1

save() <api/save>
load() <api/load>
set_up_logging() <api/set_up_logging>
```
