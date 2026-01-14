# Removed Airfoils

This file documents airfoils that have been removed from the database and the reasons for their removal. The original UIUC Airfoil Coordinates Database (https://m-selig.ae.illinois.edu/ads/coord_database.html, version updated 11/17/2024) contained some files that are incompatible with Ptera Software's parsing and interpolation requirements.

## Removal Categories

### NACA 4-Series

Ptera Software generates NACA 4-series airfoil coordinates programmatically using the standard equations, so storing these in the database is redundant. The programmatic generation also ensures consistent point density and spacing.

- n0009sm.dat
- n0012.dat
- n2414.dat
- n2415.dat
- n6409.dat
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

These files contain multiple aerodynamic elements (slat, main element, flap), or are themselves only elements of an airfoil. Ptera Software expects single-element airfoils.

- 30p-30n.dat
- R1145MSM.dat
- R1145MSF.dat
- ua79sfm.dat
- ua79sff.dat

### Non-Standard Shapes

These files represent shapes that are not traditional airfoils (e.g., cowlings, nacelles) and have only a single surface rather than upper and lower surfaces.

- naca1.dat

### Inward Trailing Edges

These files have trailing edges where the outline wraps around, causing x-coordinates to not be strictly monotonic from leading edge to trailing edge. This breaks the interpolation that assumes x always increases along each surface.

- cap21c.dat

### Deflected Control Surfaces

These airfoils have deflected control surfaces (flaps, ailerons). Ptera Software requires the given outline to represent the airfoil without any control surface deflection.

- s1221-4deg-flap.dat

### "Extended" Airfoils

These airfoils are explicitly labeled as "extended" versions of other airfoils, and therefore it isn't clear what value to use for their normalized chord length.

- e664ex.dat

### Self-Intersecting Airfoils

These airfoils have upper and lower surfaces that intersect or cross at some interior points, indicating issues with the coordinates.

- e340.dat
- e378.dat
- fx38153.dat
- fx62k131.dat
- fx63147.dat
- fx72150b.dat
- fx72ls160.dat
- mh150.dat