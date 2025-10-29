![Logo](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/Logo.png)

***

![build](https://github.com/camUrban/PteraSoftware/actions/workflows/tests.yml/badge.svg?branch=main)
![coverage](https://img.shields.io/codecov/c/gh/camUrban/PteraSoftware)
![code quality](https://img.shields.io/codefactor/grade/github/camUrban/PteraSoftware)
![source rank](https://img.shields.io/librariesio/sourcerank/pypi/PteraSoftware?color=blue&label=source%20rank)
![license](https://img.shields.io/github/license/camUrban/PteraSoftware?color=blue)
![code style](https://img.shields.io/badge/code%20style-black-black)

***

![Example Unsteady Formation Flight](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20formation/Animate.webp)

This is Ptera Software: a fast, easy-to-use, and open-source package for analyzing
flapping-wing flight.

## Motivation

In late 2018, I became curious about biological flight. To sate this curiosity, I
wanted to computationally simulate some flapping-wing fliers. I quickly realized I had
two options:

1. Spend thousands of dollars on a closed-source CFD program, which would take hours to
   solve a simple case.
2. Try to learn someone else's open-source, unsteady solver written in a language I
   didn't know, or using a framework that is overly complicated for my use case.

Neither of these seemed like the right choice.

Thankfully, my friend, Peter Sharpe, had just released his own open-source aerodynamics
solver: AeroSandbox. With his support, I have used AeroSandbox as a jumping-off point
to develop a solver package capable of unsteady simulations.

Through the combined efforts of Peter Sharpe, Suhas Kodali, and me, Ptera Software was
born. It is an easy-to-use, open-source, and actively-maintained UVLM package capable
of analyzing flapping-wing flight. Moreover, it's written in Python, is well
documented, tested, and validated.

Beginning with version 3.0.0, Ptera Software also includes a GUI developed by Zach Tait.
Although it is still rudimentary, we hope that it will help make this tool accessible to
even more users.

With your help, I hope we will increase the open-source community's interest and
understanding of biological flight.

## Features

1. Various Aerodynamic Simulation Methods
    * Steady simulations can be run with a standard horseshoe vortex-lattice method
      (VLM) or a ring VLM.
    * Unsteady simulations use a ring unsteady VLM (UVLM) solver.
    * Unsteady simulations support both fixed and free wakes.
    * Unsteady simulations implement vortex aging to reduce numerical instabilities.
2. Customizable Aircraft Geometry
    * Aircraft can be defined as a collection of one or more wings of any dimensions and
      positions.
    * Wings can be defined as a collection of two or more wing cross sections of any
      dimensions and positions.
    * Wing cross sections can be specified to match the mean camber line of an airfoil.
    * The package comes with a massive database of airfoil to chose from.
    * Wings are automatically discretized into panels with customizable sizes and
      spacings.
3. Customizable Aircraft Motion
    * The relative motion of wings and wing cross sections can be defined using any
      time-dependent functions of sweep, pitch, and heave angles.
4. Customizable Operating Points
    * Parameters such as the free-stream velocity, density, angle of attack, angle of
      sideslip, etc. can be changed by the user.
5. High-Speed Simulations
    * Using Just-In-Time compilation, Ptera Software can solve many unsteady
      flapping-wing simulations in less than a minute!
    * Steady simulations take only seconds!
6. Simulations of Formation Flight
    * Since v2.0.0, Ptera Software has supported simulations with more than one
      airplane.
    * This feature can be used to analyze the aerodynamics of flapping-wing formation
      flight!
7. Features for Flapping-Wing Vehicle Design
    * Ptera Software is focused on developing features to facilitate designing
      flapping-wing vehicles.
    * For example, use the functions in the trim module to automatically search for a
      trim operating point for steady and unsteady simulations of aircraft.
8. A Basic GUI
    * This is still in its beta stage, but we will be adding more functionality over the
      next several releases.

## Installation and Use

First things first, you will need a copy of Python 3.13, which you can download from the
official Python website. At this time, I do not recommend using a version from the 
Anaconda distribution as it could introduce compatibility issues with PyPI.

There are two ways to use Ptera Software. The first is by downloading GitHub release,
which will provide you your own copy of the source code, in which you can get a feel
for how it works (this can also be accomplished by forking the main branch). The second
is by importing the Ptera Software package using PyPI, which will allow you to call
Ptera Software's functions in your own scripts. If you are new to this tool, I
recommend first downloading a release, as this will give you access to the "examples"
directory.

Next, make sure you have an IDE in which you can run Ptera Software. I recommend using
the Community Edition of PyCharm, which is free, powerful, and well documented. If
you've never set up a Python project before, follow
[this guide](https://www.jetbrains.com/help/pycharm/quick-start-guide.html) to set up a
new project in PyCharm. If you'll be downloading a release, follow that tutorial's
"Open an existing project guide." Otherwise, follow the "Create a new project guide."

### Downloading A Release

To download a release, navigate to
[the releases page](https://github.com/camUrban/PteraSoftware/releases) and download
the latest zipped directory. Extract the contents, and set up a python project as
described in the PyCharm tutorial.

Then, open a command prompt window in your project's directory and enter:

```shell
pip install -r requirements.txt
```

via the command prompt in your fork's directory. You may also want to run:

```shell
pip install -r requirements_dev.txt
```

if you plan on making significant changes to the software.

Finally, open the "examples" folder, which contains several heavily commented scripts
that demonstrate different features and simulations. Read through each example, and
then run them to admire their pretty output!

### Importing As A Package

If you wish to use this package as a dependency in your own project, simply run:

```shell
pip install pterasoftware
```

via the command prompt in your project's directory. Then, in a script that you'd like
to use features from Ptera Software, add:

```python
import pterasoftware as ps
```

If you haven't previously downloaded Ptera Software's source code, you can also learn
about the available functions by reading their docstrings, which should be fetched
automatically by many IDEs. Otherwise, you can return to the GitHub and read through
the docstrings there.

I am hoping to implement a web-based documentation guide soon! If you'd like to
contribute to this, feel free to open a feature request issue and start a conversation!

### Using the GUI (Beta)

If you downloaded the source code from GitHub, you can try using the experimental GUI 
to run the example simulations:

```shell
python gui/main.py
```

**Note:** The GUI is currently in beta and is not included in the PyPI package. It is
only available when you download the source code from GitHub. The GUI provides a
graphical interface for running the 10 example simulations but does not yet support
custom aircraft creation. We welcome contributions if you'd like to extend it!

### What If I'm Having Trouble Getting Set Up?

Not to worry! I've made [a video](https://www.youtube.com/watch?v=oX8u2ZflJM4) that
walks through getting Ptera Software up and running. It includes every step, from 
downloading Python for the first time to setting up your IDE to running the software. 
Please note that the video demonstrates installation with Python 3.8, but for Ptera 
Software version 3.2.0 and later, you should use Python 3.13. If you still run into 
problems, feel free to open an issue for guidance.

## Example Code

The following code snippet is all that is needed (after running pip install
pterasoftware) to run the steady horseshoe solver on an airplane with custom geometry.

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

## Example Output

This package currently supports three different solvers, a steady horseshoe vortex
lattice method (VLM), a steady ring VLM, and an unsteady ring VLM (UVLM). Here are
examples of the output you can expect to receive from each of them.

### Steady Horseshoe VLM

![Example Steady Horseshoe VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples%20expected%20output/steady%20horseshoe%20vortex%20lattice%20method%20solver/Draw.webp)

### Steady Ring VLM

![Example Steady Ring VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples%20expected%20output/steady%20ring%20vortex%20lattice%20method%20solver/Draw.webp)

### Unsteady Ring VLM

![Example Unsteady Ring VLM Animation Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20static/Animate.webp)

![Example Unsteady Ring VLM Force Coefficient Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20static/Example%20Airplane%20Force%20Coefficients.png)

![Example Unsteady Ring VLM Moment Coefficient Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/main/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20static/Example%20Airplane%20Moment%20Coefficients.png)

## Requirements

Here are the requirements necessary to run Ptera Software:

* matplotlib >= 3.10.6, < 4.0.0
* numpy >= 2.2.6, < 2.3.0
* pyvista >= 0.46.3, < 1.0.0
* scipy >= 1.16.1, < 2.0.0
* numba >= 0.61.2, < 1.0.0
* cmocean >= 4.0.3, < 5.0.0
* tqdm >= 4.67.1, < 5.0.0
* webp >= 0.4.0, < 1.0.0
* PySide6 >= 6.9.2, < 7.0.0

Additionally, these packages are useful for continued development of the software:

* codecov >= 2.1.13, < 3.0.0
* black >= 25.1.0, < 26.0.0
* codespell >= 2.4.1, < 3.0.0
* pre-commit >= 4.3.0, < 5.0.0
* build >= 1.3.0, < 2.0.0
* twine >= 6.1.0, < 7.0.0
* PyInstaller >= 6.15.0, < 7.0.0
* setuptools >= 80.9.0, < 81.0.0
* wheel >= 0.45.1, < 0.46.0

## Validation

Since the release of version 1.0.0, Ptera Software is now validated against
experimental flapping-wing data! See the "validation" directory to run the test case
and read a report on the software's accuracy.

## How to Contribute

As I said before, the primary goal of this project is to increase the open-source
community's understanding and appreciation for unsteady aerodynamics in general and
flapping-wing flight in particular. This will only happen through your participation.
Feel free to request features, report bugs or security issues, and provide suggestions.
No comment is too big or small!

Before contributing, make sure to read through the
[CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to best help out.

Here is a list of changes I would like to make in the coming releases. If you want to
contribute and don't know where to start, this is for you!

### Testing

* We should make sure that all the integration tests compare output against expected
  results. This means getting rid of all the "test_method_does_not_throw" tests.
* We should maintain the repository's testing coverage to be at least 80%.

### Style and Documentation

* Maintain the repository's A CodeFactor Rating.
* We should fill in any of the "Properly document this..." TODO statements.
* We should ensure that all files be at least 30% comment lines.
* We should continue to ensure that all source code is formatted using Black.

### Features

* We should implement a leading-edge model to account for flow separation. See
  "Modified Unsteady Vortex-Lattice Method to Study Flapping Wings in Hover Flight." by
  Bruno Roccia, Sergio Preidikman, Julio Massa, and Dean Mook for details.
* We should try to implement aeroelastic effects in Ptera Software's solvers.
* Flapping wing controls is both fascinating and complicated. We should try to create a
  workflow in Ptera Software for controls systems identification for flapping-wing
  vehicles.

## Credits

Here is a list of all the people and packages that helped me created Ptera Software in
no particular order. Specific citations can be found in the source code's docstrings
where applicable.

* Suhas Kodali
* Peter Sharpe
* Zach Tait
* Jonah Jaffe
* Ramesh Agarwal
* E. Farrell Helbling
* Raphael Zufferey
* Joseph Katz
* Allen Plotkin
* Austin Stover
* AeroSandbox
* Black
* Codecov
* NumPy
* SciPy
* PyVista
* MatPlotLib
* Numba
* Pre-Commit
* SetupTools
* GitIgnore
* Shields.io
* PyPI
* Wheel
* Twine
* SemVer
* GitFlow
* GitHub Flow
* Cmocean
* Tqdm
* WebP
* Build

## Notes

To the best of my ability, I am following SemVer conventions in naming my releases. I
am also using the 
[GitHub Flow](https://docs.github.com/en/get-started/using-github/github-flow) method 
of branching for this project's development, with a version bump and deployment to 
GitHub and PyPI about once per month, plus on-demand releases for critical bug 
fixes.
