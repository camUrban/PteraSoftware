# Removed Airfoils

This file documents airfoils that have been removed from the database and the reasons for their removal. The original UIUC Airfoil Coordinates Database (https://m-selig.ae.illinois.edu/ads/coord_database.html, version updated 11/17/2024) contained some files that are incompatible with Ptera Software's parsing and interpolation requirements.

## Removal Categories

### NACA 4-Series

Ptera Software generates NACA 4-series airfoil coordinates programmatically using the standard equations, so storing these in the database is redundant. The programmatic generation also ensures consistent point density and spacing.

- naca0006.dat
- naca0008.dat
- naca0010.dat
- naca0015.dat
- naca0018.dat
- naca0021.dat
- naca0024.dat
- naca1408.dat
- naca1410.dat
- naca1412.dat
- naca2408.dat
- naca2410.dat
- naca2411.dat
- naca2412.dat
- naca2415.dat
- naca2418.dat
- naca2421.dat
- naca2424.dat
- naca4412.dat
- naca4415.dat
- naca4418.dat
- naca4421.dat
- naca4424.dat
- naca6412.dat

### Multi-Element Airfoils

These files contain multiple aerodynamic elements (slat, main element, flap) separated by comment lines. Ptera Software expects single-element airfoils.

- 30p-30n.dat

### Non-Standard Shapes

These files represent shapes that are not traditional airfoils (e.g., cowlings, nacelles) and have only a single surface rather than upper and lower surfaces.

- naca1.dat