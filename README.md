![Logo](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/website/Logo.png)

***

[![DOI](https://zenodo.org/badge/249337717.svg)](https://doi.org/10.5281/zenodo.19229119)
![license](https://img.shields.io/github/license/camUrban/PteraSoftware?color=blue)
![build](https://github.com/camUrban/PteraSoftware/actions/workflows/tests.yml/badge.svg?branch=main)
![coverage](https://img.shields.io/codecov/c/gh/camUrban/PteraSoftware)
![python](https://img.shields.io/pypi/pyversions/pterasoftware)
![types](https://img.shields.io/pypi/types/pterasoftware)
![code style](https://img.shields.io/badge/code%20style-black-black)
![source rank](https://img.shields.io/librariesio/sourcerank/pypi/PteraSoftware?color=blue&label=source%20rank)

***
<!-- docs-include-start -->
![Flapping Wings in Ground Effect](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/hero_graphics/hero_animated.webp)

This is Ptera Software: a fast, easy-to-use, and open-source package for analyzing flapping-wing flight.

## Quick Start

Install the package with pip (requires Python 3.11, 3.12, or 3.13):

```shell
pip install pterasoftware
```

Then run a simulation. The following snippet defines a simple rectangular wing, solves for its aerodynamics using the steady horseshoe vortex lattice method (VLM), and visualizes the results:

```python
import pterasoftware as ps

airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=8,
                    control_surface_symmetry_type="asymmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0, 5, 0),
                    control_surface_symmetry_type="asymmetric",
                ),
            ],
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
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

solver.run()

ps.output.draw(solver=solver, scalar_type="lift", show_streamlines=True)
```

## Features

1. Various Aerodynamic Simulation Methods
    * Steady simulations can be run with a standard horseshoe vortex-lattice method (VLM) or a ring VLM.
    * Unsteady simulations use a ring unsteady VLM (UVLM) solver.
    * Unsteady simulations support both fixed and free wakes.
    * Unsteady simulations implement vortex aging to reduce numerical instabilities.
    * All three solvers support surface effects (e.g., ground effect) via the method of images.
2. Customizable Aircraft Geometry
    * Aircraft can be defined as a collection of one or more wings of any dimensions and positions.
    * Wings can be defined as a collection of two or more wing cross sections of any dimensions and positions.
    * Wing cross sections can be specified to match the mean camber line of an airfoil.
    * The package comes with a massive database of airfoils to choose from, courtesy of the [UIUC Airfoil Coordinates Database](https://m-selig.ae.illinois.edu/ads/coord_database.html).
    * Wings are automatically discretized into panels with customizable sizes and spacings.
3. Customizable Aircraft Motion
    * The relative motion of wings and wing cross sections can be defined using any time-dependent functions of sweep, pitch, and heave angles.
4. Customizable Operating Points
    * Parameters such as the free-stream velocity, density, angle of attack, angle of sideslip, etc. can be changed by the user.
5. High-Speed Simulations
    * Using Just-In-Time compilation, Ptera Software can solve many unsteady flapping-wing simulations in less than a minute!
    * Steady simulations take only seconds!
6. Simulations of Formation Flight
    * Since v2.0.0, Ptera Software has supported simulations with more than one airplane.
    * This feature can be used to analyze the aerodynamics of flapping-wing formation flight!
7. Save and Load Simulation Results
    * Save solved simulations to JSON files and load them back without re-running.
    * Uses JSON serialization instead of pickle, avoiding arbitrary code execution vulnerabilities.
    * Supports gzip compression for reduced file sizes.
    * Loaded objects are fully compatible with all output and visualization functions.
8. Features for Flapping-Wing Vehicle Design
    * Ptera Software is focused on developing features to facilitate designing flapping-wing vehicles.
    * For example, use the functions in the trim module to automatically search for a trim operating point for steady and unsteady simulations of aircraft.
9. A Basic GUI
    * This is still in its alpha stage, but we will be adding more functionality soon.

## Installation

### As a Package

If you haven't already, install Ptera Software from PyPI (see [Quick Start](#quick-start) above):

```shell
pip install pterasoftware
```

Your IDE should automatically provide docstring hints for the available classes and functions. For more detailed documentation, visit the [Ptera Software documentation site](https://docs.pterasoftware.com/).

### From Source

If you want to browse the example scripts or dig into the source code, you will need a local copy of the repository. Follow the environment setup instructions in the [Contributing Guidelines](CONTRIBUTING.md#contributing-code) to clone the repository, create a virtual environment, and install dependencies.

Once set up, the `examples/` directory contains scripts that demonstrate the full range of Ptera Software's features and solvers. These scripts are also available on the [documentation site](https://docs.pterasoftware.com/en/latest/examples.html).

## Example Output

This package currently supports three different solvers, a steady horseshoe VLM, a steady ring VLM, and an unsteady ring VLM (UVLM). Here are examples of the output you can expect to receive from each of them.

### Steady Horseshoe VLM

![Example Steady Horseshoe VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples_expected_output/steady_horseshoe_vortex_lattice_method_solver/Draw.webp)

### Steady Ring VLM

![Example Steady Ring VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples_expected_output/steady_ring_vortex_lattice_method_solver/Draw.webp)

### Unsteady Ring VLM

![Example Unsteady Ring VLM Animation Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples_expected_output/unsteady_ring_vortex_lattice_method_solver_static/Animate.webp)

## Validation

Since the release of version 1.0.0, Ptera Software is now validated against experimental flapping-wing data! See the `validation/` directory to run the test case and read a report on the software's accuracy.

## Documentation

For detailed API documentation and guides, visit the [Ptera Software documentation site](https://docs.pterasoftware.com/).

## How to Contribute

The primary goal of this project is to increase the open-source community's understanding and appreciation for unsteady aerodynamics in general and flapping-wing flight in particular. This will only happen through your participation. Feel free to request features, report bugs or security issues, and provide suggestions.

Before contributing, make sure to read through the [Contributing Guidelines](CONTRIBUTING.md) for how to best help out.

## Contributors

* Cameron Urban ([camUrban](https://github.com/camUrban))
* Zach Tait ([Zach10a](https://github.com/Zach10a))
* Jonah Jaffe ([JonahJ27](https://github.com/JonahJ27))
* Venkata Akhil Mettu ([AKHIL-149](https://github.com/AKHIL-149))
* Savitha N ([Savitha-Akhilu](https://github.com/Savitha-Akhilu))
* Pedro Bornia ([BorniaPedro](https://github.com/BorniaPedro))
* Mohamed Abdulghany ([MohamedMG7](https://github.com/MohamedMG7))
* Hang Haotian ([haotianh9](https://github.com/haotianh9))

### Supporters

* Peter Sharpe
* Suhas Kodali
* Ramesh Agarwal
* E. Farrell Helbling
* Raphael Zufferey
* Joseph Katz
* Allen Plotkin
* Austin Stover

## Background

Ptera Software grew out of a desire to make flapping-wing aerodynamics accessible without expensive commercial CFD tools or hard-to-use open-source alternatives. Initially built on [AeroSandbox](https://github.com/peterdsharpe/AeroSandbox) with the support of Peter Sharpe and Suhas Kodali, it has developed into an actively-maintained UVLM package that is well documented, tested, and validated. We hope that with your help, we will increase the open-source community's interest and understanding of biological flight.
<!-- docs-include-end -->
