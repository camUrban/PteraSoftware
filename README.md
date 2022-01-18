# Ptera Software

![Ptera Software Logo](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/PteraSoftwareLogo.jpg)

![build](https://img.shields.io/travis/camUrban/PteraSoftware/master)
![coverage](https://img.shields.io/codecov/c/gh/camUrban/PteraSoftware/master)
![code quality](https://img.shields.io/codefactor/grade/github/camUrban/PteraSoftware/master)
![source rank](https://img.shields.io/librariesio/sourcerank/pypi/PteraSoftware?color=blue&label=source%20rank)
![license](https://img.shields.io/github/license/camUrban/PteraSoftware?color=blue)
![code style](https://img.shields.io/badge/code%20style-black-black)

This is Ptera Software: a fast, easy-to-use, and open-source package for analyzing
flapping-wing flight.

## Motivation

In late 2018, I became curious about biological flight. To sate this curiosity, I wanted
 to computationally simulate some flapping-wing fliers. I quickly realized I had two
 options:

1. Spend thousands of dollars on a closed-source CFD program, which would take hours to
solve a simple case.
2. Try to learn someone else's open-source, unsteady solver written in a language I
didn't know, and using a framework that is overly complicated for my use case.

Neither of these seemed like the right choice.

Thankfully, my friend, Peter Sharpe, had just released his own open-source aerodynamics
solver: AeroSandbox. With his blessing, I have used AeroSandbox as a jumping-off point
to develop a solver package capable of unsteady simulations.

Through the combined efforts of Peter Sharpe, Suhas Kodali, and me, Ptera Software 
was born. It is an easy-to-use, open-source, and actively-maintained UVLM package 
capable of analyzing flapping-wing flight. Moreover, it's written in Python, is well 
documented, tested, and validated.

With your help, I hope we will increase the open-source community's interest and
understanding of biological flight.

## How to Install

First things first, you will need a copy of Python 3.7 or 3.8. Python 3.9 is not yet
supported due to a dependency issue in VTK. Download Python 3.7 or 3.8 from the official
 Python website. At this time, I do not recommend using a version from the Anaconda
 distribution as it could introduce compatibility issues with PyPI.

There are a few ways to install Ptera Software. If you wish to use this package as a
dependency in your own projects, simply run:

```pip install pterasoftware```

via the command prompt in your project's directory.

If you just want to play around with the software, feel free to fork this repository and
 open the source code in the IDE of your choice. You will then need to run:

```pip install -r requirements.txt```

via the command prompt in your fork's directory.

### Requirements

Here are the requirements necessary to run Ptera Software:

* matplotlib >= 3.5.0, < 4.0.0
* numpy >= 1.21.0, < 1.22.0
* pyvista >= 0.33.0, < 1.0.0
* scipy >= 1.7.0, < 2.0.0
* numba >= 0.55.0, < 1.0.0
* cmocean >= 2.0.0, < 3.0.0
* tqdm >= 4.62.0, < 5.0.0

### What if I am Having Trouble Getting the Package Up And Running?

Not to worry! I am working on a video that walks through getting Ptera Software up and
running. It will include every step, from downloading Python for the first time to
setting up your IDE to running the software. Feel free to reach out for guidance. You
can reach me at camerongurban@gmail.com.

## How to Use

By reading this file, you are already off to a good start! After installing Ptera
Software in the way that best suits your use case, the next step would be to open the
"examples" directory and read through each heavily commented script. Each one will give
you insight into the software's interface. After you finish reading, try running the
scripts and admiring their pretty output!

## Example Code

The following code snippet is all that is needed (after running pip install
pterasoftware) to run the steady horseshoe solver on a custom airplane object.

```
import pterasoftware as ps

example_airplane = ps.geometry.Airplane(
    wings=[
        ps.geometry.Wing(
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    airfoil=ps.geometry.Airfoil(name="naca2412",),
                ),
                ps.geometry.WingCrossSection(
                    y_le=5.0, airfoil=ps.geometry.Airfoil(name="naca2412",),
                ),
            ],
        ),
    ],
)

example_operating_point = ps.operating_point.OperatingPoint()

example_problem = ps.problems.SteadyProblem(
    airplane=example_airplane, operating_point=example_operating_point,
)

example_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
    steady_problem=example_problem
)

example_solver.run()

ps.output.draw(
    solver=example_solver, show_delta_pressures=True, show_streamlines=True,
)
```

## Example Output

This package currently supports three different solvers, a steady horseshoe vortex
lattice method (VLM), a steady ring VLM, and an unsteady ring VLM (UVLM). Here are
examples of the output you can expect to receive from each of them.

### Steady Horseshoe VLM

![Example Steady Horseshoe VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/steady%20horseshoe%20vortex%20lattice%20method%20solver%20example%20expected%20output/Draw%20Output.jpg)

### Steady Ring VLM

![Example Steady Ring VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/steady%20ring%20vortex%20lattice%20method%20solver%20example%20expected%20output/Draw%20Output.jpg)

### Unsteady Ring VLM

![Example Unsteady Ring VLM Animation Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20example%20expected%20output/Animate%20Output.gif)

![Example Unsteady Ring VLM Force Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20example%20expected%20output/Plot%20Output%201.png)

![Example Unsteady Ring VLM Moment Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20example%20expected%20output/Plot%20Output%203.png)

## Validation

Since the release of version 1.0.0, Ptera Software is now validated against experimental
flapping-wing data! See the "validation" directory to run the test case and read a
report on the software's accuracy.

## How to Contribute

As I said before, the primary goal of this project is to increase the open-source
community's understanding and appreciation for unsteady aerodynamics in general and
flapping-wing flight in particular. This will only happen through your participation.
Feel free to request features, report bugs and security issues, and provide suggestions.
 No comment is too big or small!

Here is a list of changes I would like to make in the coming releases. If you want to
contribute and don't know where to start, this is for you!

### Testing

* We should make sure that all the integration tests compare output against expected
results. This means getting rid of all the "test_method_does_not_throw" tests.
* We should maintain the repository's testing coverage to be at least 80%.

### Style and Documentation

* Maintain the repository's A CodeFactor Rating.
* We should fill in any of the "Properly document this..." TODO statements.
* We should ensure that all files have between 30% and 70% comment lines.
* We should continue to ensure that all source code is formatted using Black.

### Features

* We should create a setup tutorial video and add it to the documentation. This should
be geared toward a user who doesn't have Python, an IDE, or Ptera Software installed on
their computer yet.
* We should implement a leading-edge model to account for flow separation. See 
  "Modified Unsteady Vortex-Lattice Method to Study Flapping Wings in Hover Flight." 
  by Bruno Roccia, Sergio Preidikman, Julio Massa, and Dean Mook for details.
* We should create a command-line interface or GUI.
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
* Ramesh Agarwal
* Joseph Katz
* Allen Plotkin
* Austin Stover
* AeroSandbox
* Black
* Codecov
* Travis CI
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
* Cmocean
* Tqdm

## Notes

To the best of my ability, I am following SemVer conventions in naming my releases. 
I am also using the GitFlow method of branching for this project's development. This 
means that nightly builds will be available on the develop branch. The latest stable 
releases can be found on the master branch.