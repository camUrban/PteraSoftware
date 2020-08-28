# Ptera Software

![Ptera Software Logo](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/PteraSoftwareLogo.jpg)

![Build Status](https://img.shields.io/travis/camUrban/PteraSoftware)
![Percent Coverage](https://img.shields.io/codecov/c/gh/camUrban/PteraSoftware)
![Code Quality Grade](https://img.shields.io/codefactor/grade/github/camUrban/PteraSoftware)
![MIT License](https://img.shields.io/github/license/camUrban/PteraSoftware?color=blue)
![Black Code Style](https://img.shields.io/badge/code%20style-black-black)

This is Ptera Software: a fast, easy-to-use, and open source package for analyzing flapping-wing flight.

## Motivation

Around a year ago, I became curious about biological flight. To sate this curiosity, I wanted to computationally
simulate some flapping-wing fliers. I quickly realized I had two options:

1. Spend thousands of dollars on a closed-source CFD program, which would take hours to solve a simple case.
2. Try to learn someone else's open source, unsteady solver written in a language I didn't know, and using a framework
   that is overly complicated for my use case.

Neither of these seemed like a good choice.

Thankfully, my friend, Peter Sharpe had just released his own open-source aerodynamics solver: AeroSandbox. With his
blessing, I have used AeroSandbox as a jumping off point to develop a solver package capable of unsteady simulations.

Through the combined efforts of Peter Sharpe, Suhas Kodali, and myself, Ptera
Software was born. It is the only fast, easy-to-use, and open-source package I know of that is capable of analyzing
flapping wing flight. Moreover, it's written in Python, is well documented, and is well tested.

With your help, I hope we will improve the open-source communities interest and understanding of biological flight.

## How to Install

### Standard Method

First things first, you will need a copy of Python 3.7.6 or 3.7.7. Download it from the official Python website.

There are a few ways to install Ptera Software. If you wish to use this package as a dependency in your own projects,
simply run:

```pip install PteraSoftware```

via command prompt in your project's directory.

If you just want to play around with the software, feel free to fork this repository and open the source code in the IDE
of your choice. You will then need to run:

```pip install -r REQUIREMENTS.txt```

via command prompt in your fork's directory.

### Anaconda Method

Download a copy of Python 3.7.6 or 3.7.7 via the Anaconda distribution.

To be continued...

### Requirements

Here are the requirements necessary to run Ptera Software:

* matplotlib >= 3.2.2, < 4.0.0
* numpy >= 1.18.5, < 1.19.0
* pyvista >= 0.25.3, < 1.0.0
* scipy >= 1.5, < 2.0

### What if I am Having Trouble Getting the Package Up And Running?

Not to worry! I working on releasing a video that walks through getting Ptera Software up and running. It will include
every step, from downloading Python for the first time, to setting up your IDE to run it. Feel free to reach out for
guidance. You can reach me at camerongurban@gmail.com.

## How to Use

By reading this file, you are already off to a good start! After installing Ptera Software in the way that best suits
your use case, the next step would be to open the examples directory, and read through each heavily commented script.
Each one will give you insight into the software's interface. After you finish reading, try running them and admiring
the pretty output!

## Example Code

The following code snippet is all that is needed (after running pip install pterasoftware), run the steady horseshoe
solver on a custom airplane object.

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

This package currently supports three different solvers, a steady horseshoe vortex lattice method (VLM), a steady ring
VLM, and an unsteady ring VLM. Here are examples of the output you can expect to receive from each of them.

### Steady Horseshoe VLM

![Example Steady Horseshoe VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/steady%20horseshoe%20vortex%20lattice%20method%20solver%20example%20expected%20output/Draw%20Output.jpg)

### Steady Ring VLM

![Example Steady Ring VLM Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/steady%20ring%20vortex%20lattice%20method%20solver%20example%20expected%20output/Draw%20Output.jpg)

### Unsteady Ring VLM

![Example Unsteady Ring VLM Animation Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20example%20expected%20output/Animate%20Output.gif)

![Example Unsteady Ring VLM Force Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20example%20expected%20output/Plot%20Output%201.png)

![Example Unsteady Ring VLM Moment Output](https://raw.githubusercontent.com/camUrban/PteraSoftware/master/docs/examples%20expected%20output/unsteady%20ring%20vortex%20lattice%20method%20solver%20variable%20example%20expected%20output/Plot%20Output%203.png)

## How to Contribute

As I said before, the primary goal of this project is to increase the open source communities understanding and
appreciation for unsteady aerodynamics in general, and flapping-wing flight in particular. This will only happen through
your participation. Feel free to request features, report bugs and security issues, and to provide suggestions. No
comment is to big or small!

Note: This is my first attempt at creating a Python project this size, and at managing a public repository. I am bound
to make mistakes but I am determined to listen to your feedback and to use it to improve.

## Credits

Here is a list, in no particular order, of all the people and packages that helped me created Ptera Software. Specific
citations can be found in the source code's docstrings where applicable.


* Suhas Kodali
* Peter Sharpe
* Joseph Katz
* Allen Plotkin
* AeroSandbox
* Black
* Coverage
* Travis CI
* NumPy
* SciPy
* PyVista
* MatPlotLib
* Pre-Commit
* SetupTools
* GitIgnore
* Shields.io
* PyPI
* Wheel
* Twine
* SemVer
* GitFlow

## Notes

To the best of my ability, I am following SemVer conventions in naming my releases. I am also using the GitFlow method
of branching for this project's development. This means that nightly builds will be available on the develop branch, and
the latest stable releases can be found on the master branch.
