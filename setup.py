import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PteraSoftware",
    version="0.1.1",
    author="Cameron Urban",
    author_email="camerongurban@gmail.com",
    description="This is an open-source, unsteady aerodynamics solver for analyzing flapping-wing flight.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/camurban/pterasoftware",
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        "Natural Language :: English",
        'Programming Language :: Python :: 3.7',
    ],
    keywords='aerospace computational-biology airplane cfd computational-fluid-dynamics aerodynamics aeronautics aerospace-engineering unmanned-aerial-system aircraft-design unmanned-aerial-vehicle ornithopter ornithology vortex-lattice-method unsteady-flows vlm potential-flow',
    python_requires='>= 3.7.7, < 3.8',
    install_requires=[
        'matplotlib >= 3.2.2, < 4'
        'numpy >= 1.18.5, < 1.19.0',
        'pyvista >= 0.25.3, < 1',
        'scipy >= 1.5, < 2',
        'codecov >= 2.1.8, < 3',
        'numexpr >= 2.7.1, < 3',
        'pre-commit >= 2.6, < 3',
        'black >= 19.10b0, < 20',
    ],
)
