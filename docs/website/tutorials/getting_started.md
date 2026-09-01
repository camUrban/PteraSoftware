
# Getting Started

This page introduces the getting started notebook, which walks through a first
steady horseshoe vortex lattice method (VLM) simulation from scratch.

```{literalinclude} ../../../tutorials/getting_started.ipynb
:language: ipynb
```

## Running the Notebook

The notebook is designed to run top to bottom in any Jupyter environment with
Ptera Software installed:

```shell
pip install pterasoftware
jupyter notebook tutorials/getting_started.ipynb
```

It covers the full object model in small steps:

1. `Airfoil`
2. `WingCrossSection`
3. `Wing`
4. `Airplane`
5. `OperatingPoint`
6. `SteadyProblem`
7. `SteadyHorseshoeVortexLatticeMethodSolver`
8. Results, visualization, logging, and serialization

The same object model is shared by every solver in Ptera Software, so once you
are comfortable with this tutorial you can move on to the examples pages and
try the ring VLM, unsteady, aeroelastic, and free flight solvers.
