# Axes, Points, and Frames

In simulating flapping-wing dynamics and aerodynamics, Ptera Software uses several vector-valued quantities. To avoid ambiguity, these vectors often require additional details about their interpretation, such as axis systems, reference points, and reference frames.

This document lays out Ptera Software's notation for defining vectors using these constructs. The notation and terminology I use is an extended version of that introduced in "Flight Vehicle Aerodynamics" by Mark Drela.

## Axis Systems vs. Reference Points vs. Reference Frames

An axis system, also called "axes," contains information about three directions. In Ptera Software, all axes are Cartesian, meaning their three basis directions are linear, not angular.

Reference points, also called "points," contain information about the location of a particular point in space.

Lastly, a reference frame, also called a "frame," contains information about the location of an "observer," and their motion relative to what is observed.

Consider the arbitrary vector **r**, which exists in 3D space. For now, let's say that **r** is a force vector. In order to express **r** using components (x, y, z), we must, at a minimum, pick an axis system. If instead **r** is a position vector, we need both axes and a reference point to serve as an origin, so we must pick both before writing down **r**'s three components. The same is true if **r** is a moment, but now the reference point no longer serves as an origin, but instead the point about which the moment acts. Lastly, if **r** is some time derivative of position, such as a velocity or acceleration vector, then we no longer need a reference point, but we do require both an axis system and a reference frame.

Due to the nested structure of Ptera Software's geometry objects, in practice, many vector-valued quantities like positions and moments, use reference points and axes that are defined locally within a given object. An example of this next structure for an unsteady vortex lattice method simulation is shown below.

<img src="ObjectHierarchy.jpg" alt="Object Hierarchy" width="400"/>

For information on how axes are defined relative to one another, and how vectors can be transformed within an axis system, read through [Angle Vectors and Transformations](https://raw.githubusercontent.com/camUrban/PteraSoftware/feature/improved_geometry_definitions/docs/ANGLE_VECTORS_AND_TRANSFORMATIONS.md)

# Specifying Axes, Points, and Frames

Given the varied requirements for vector-valued quantities, it is important that we are very specific when assigning them variable names or referencing them in text. Also, due to the hierarchical structure of Ptera Software's objects, additional specificity may be required depending on the context. For example, if we use a force vector within the Wing class that references wing axes, we still need to specify that this vector is given in wind axes, but we don't (and can't) specify which of the parent Airplane's Wing's axes we mean. In contrast, if we declare a variable inside the Wing class that references wing cross section axes, we must specify which of the Wing's WingCrossSections's axes we are referring to.

## Patterns

There are three useful combinations of axes, points, and frames. For variables that fall into each of these three cases, we denote them by appending information to their variable names using **IDs**. When referencing the variables in comments and docstrings, we add this additional information parenthetically using **names**:

1. Axes without a point and without a frame  
   \[variable name\]\_\[axes ID\]  
   "\[variable name\] (in \[axes name\])"
2. Axes without a point and with a frame  
   \[variable name\]\_\[axes ID\]\_\_\[frame ID\]  
   "\[variable name\] (in \[axes name\], observed from the \[frame name\])"  
3. Axes with a point and without a frame  
   \[variable name\]\_\[axes ID\]\_\[point ID\]  
   "\[variable name\] (in \[axes name\], relative to the \[point name\])"

The correct name and ID for a particular axis system, point, or frame depends on the level of context. However, in all cases IDs consist of a series of abbreviations, moving in scope from most specific to least specific. By contrast, names move from least specific to most specific. Also, in contrast with IDs, the exact syntax for names is slightly flexible to allow for the description to sound correct in plain English.

The standard abbreviations and names are given below for reference. See the section for a particular axis system, point, or frame for examples of the correct IDs and names in various contexts.

## ID Abbreviations and Names

* E: Earth
* B: body
* P: airplane
* W: wind
* Pr: problem
* G: geometry
* Wn: wing
* Wcs…: wing cross section
  * …i: inner
  * …o: outer
* Wcsp: wing cross section parent
* A…: airfoil
  * …i: inner
  * …o: outer
* I: simulation starting point
* Cgi: starting point
* Cg: CG point
* Ler: leading edge root point
* Lp: leading point
* Lpp: leading point parent
* …pp…: panel point
  * Fr…: front right
  * Fo…: forward outer
  * Fl…: front left
  * Fi…: forward inner
  * Bl…: back left
  * Bi…: backward inner
  * Br…: back right
  * Bo…: backward outer
  * C…: collocation
  * …r\[m\]c\[n\]: (m, n)
* …bhvp…: bound horseshoe vortex point
  * Fr…: front right
  * Fl…: front left
  * Bl…: back left
  * Br…: back right
  * …r\[m\]c\[n\]: (m, n)
* …brvp…: bound ring vortex point
  * Fr…: front right
  * Fl…: front left
  * Bl…: back left
  * Br…: back right
  * …r\[m\]c\[n\]: (m, n)
* …wrvp…: wake ring vortex point
  * Fr…: front right
  * Fl…: front left
  * Bl…: back left
  * Br…: back right
  * …r\[m\]c\[n\]: (m, n)
* …whvp…: wake horseshoe vortex point
  * Fr…: front right
  * Fl…: front left
  * Bl…: back left
  * Br…: back right
  * …\[n\]: n
* …lvp…: line vortex point
  * S…: start
  * E…: end
  * C…: center
  * …f: front leg
  * …l: left leg
  * …b: back leg
  * …r: right leg

# Axis Systems

## 1. The Earth axis system

* Basis directions  
  1. +x: North  
  2. +y: East  
  3. +z: Down  
* Right-handed  
* Ownership: None  
* References  
  * Text: …in Earth axes…  
  * Variables: …\_E…

## 2. Body axes

* Basis directions  
  1. +x: Towards the front of the Airplane  
  2. +y: Towards the right of the Airplane  
  3. +z: Towards the bottom of the Airplane  
* Right-handed
* Ownership: Airplane  
* Local reference examples  
  * Text: …in body axes…  
  * Variables: …\_B…  
* Non-local reference examples  
  * Text: …in the first Airplane's body axes…  
  * Variables: …\_BP1…

## 3. Wind axes

* Caveat: We assume a still airmass, so the freestream velocity observed from the body frame is solely due to the Airplane's velocity observed from the Earth frame.  
* Basis directions  
  1. +x: In line with (parallel, not anti-parallel, to) the freestream velocity vector observed from the body frame  
  2. +y: In the direction perpendicular to first and third components\*  
  3. +z: In the direction perpendicular to first and second components\*  
* \*There are infinite options for the second and third components that satisfy the perpendicularity requirement. Therefore, we define them using a thought experiment: Imagine three unit vectors pointing along the geometry axes' basis directions, as defined in Earth axes. There is exactly one pair of angles, which we'll call α and β, that we can use to perform a y-z series of **extrinsic** rotations to construct wind axes from geometry axes that will exactly align x-axis with the freestream velocity vector observed from the body frame. This series of rotations also constructs the wind axes +y and +z basis directions.  
The two angles α and β are referred to as the angle of attack and the angle of sideslip. Wind axes are commonly defined using these angles. This is because they are intuitively understood by many aerodynamicists: in the simplest scenarios, a positive α corresponds to the Airplane's nose pointing above its direction of travel, and a positive β to its nose pointing to the left of its direction of travel. However, this can seem a bit cyclical, and it obscures some subtlety in their definition: defining α and β based on the convention described previously, allows us to define lift as the aerodynamic force's component in the wind axes' z-axis, thereby making lift independent of sideslip.  
* Right-handed
* Ownership: SteadyProblem or UnsteadyProblem  
* Local reference examples  
  * Text: …in wind axes…  
  * Variables: …\_W…  
* Non-local reference examples  
  * Text: …in the first Problem's wind axes…  
  * Variables: …\_WPr1…

## 4. Geometry axes

* Basis directions  
  1. +x: Towards the back of the Airplane (aft)  
  2. +y: Towards the right of the Airplane  
  3. +z: Towards the top of Airplane  
* Right-handed  
* Ownership: Airplane  
* Local reference examples  
  * Text: …in geometry axes…  
  * Variables: …\_G…  
* Non-local reference examples  
  * Text: …in the first Airplane's geometry axes…  
  * Variables: …\_GP1…

## 5. Wing axes

* Basis directions  
  1. +x: Towards the back of the Wing at its root  
  2. +y: From a Wing's root towards its second WingCrossSection's plane  
  3. +z: Towards the top surface of the Wing  
* Right-handed for non-symmetric and symmetric-continuous Wings. Left-handed for mirror-only Wings.  
* Ownership: Wing  
* Local reference examples  
  * Text: …in wing axes…  
  * Variables: …\_Wn…  
* Airplane-local reference examples  
  * Text: …in the first Wing's axes…  
  * Variables: …\_Wn1…  
* Non-local reference examples  
  * Text: …in the first Airplane's second Wing's axes…  
  * Variables: …\_Wn2P1…

## 6. Wing cross section axes

* Basis directions  
  1. +x: Towards the trailing edge in the WingCrossSection's plane  
  2. +y: Normal to the WingCrossSection's plane in the direction of the next WingCrossSection  
  3. +z: Towards the top surface of the Wing  
* Right-handed for WingCrossSections of non-symmetric and symmetric-continuous Wings. Left-handed for WingCrossSections of mirror-only Wings.  
* Ownership: WingCrossSection  
* Local reference examples  
  * Text: …in wing cross section axes…  
  * Variables: …\_Wcs…  
* Wing-local reference examples  
  * Text: …in the first WingCrossSection's axes…  
  * Variables: …\_Wcs1…  
* Airplane-local reference examples  
  * Text: …in the second Wing's third WingCrossSection's axes…  
  * Variables: …\_Wcs3Wn2…  
* Non-local reference examples  
  * Text: …in the first Airplane's second Wing's first WingCrossSection's axes…  
  * Variables: …\_Wcs1Wn2P1…

## 7. Wing cross section parent axes

* Basis directions identical to Wing axes for a Wing's first WingCrossSection, and identical to the previous WingCrossSection's axes for subsequent ones.  
* Right-handed for WingCrossSections of non-symmetric and symmetric-continuous Wings. Left-handed for WingCrossSections of mirror-only Wings.  
* Ownership: WingCrossSection  
* Local reference examples  
  * Text: …in wing cross section parent axes…  
  * Variables: …\_Wcsp…  
* Wing-local reference examples  
  * Text: …in the second WingCrossSection's parent axes…  
  * Variables: …\_Wcsp1…  
* Airplane-local reference examples  
  * Text: …in the second Wing's third WingCrossSection's parent axes…  
  * Variables: …\_Wcsp3Wn2…  
* Non-local reference examples  
  * Text: …in the first Airplane's second Wing's first WingCrossSection's parent axes…  
  * Variables: …\_Wcsp1Wn2P1…

## 8. Airfoil axes

* Basis directions  
  1. +x: Chordwise towards the Airfoil's trailing point  
  2. +y: Normal to the chord towards the Airfoil's upper line  
* Two-dimensional  
* Ownership: Airfoil  
* Local reference examples  
  * Text: …in airfoil axes…  
  * Variables: …\_A…  
* Wing-local reference examples  
  * Text: …in the second WingCrossSection's Airfoil's axes…  
  * Variables: …\_AWcs2…  
* Airplane-local reference examples  
  * Text: …in the second Wing's third WingCrossSection's Airfoil's axes…  
  * Variables: …\_AWcs3Wn2…  
* Non-local reference examples  
  * Text: …in the first Airplane's second Wing's first WingCrossSection's Airfoil's axes…  
  * Variables: …\_AWcs1Wn2P1…

# Reference Points

## 1. Simulation starting point

* Position of the first Airplane's CG at the start of the simulation  
* Ownership: None  
* Reference examples  
  * Text: …relative to the simulation starting point…  
  * Variables: …\_I

## 2. Starting point

* Position of the Airplane's CG at the start of the simulation  
* Ownership: None  
* Local reference examples  
  * Text: …relative to the starting point…  
  * Variables: …\_Cgi  
* Non-local reference examples  
  * Text: …relative to the first Airplane's starting point…  
  * Variables: …\_CgiP1

## 3. CG point

* Position of the Airplane's CG  
* Ownership: Airplane  
* Local reference examples  
  * Text: …relative to the CG point…  
  * Variables: …\_Cg  
* Non-local reference examples  
  * Text: …relative to the second Airplane's CG point…  
  * Variables: …\_CgP2

## 4. Leading edge root point

* Root point of the Wing's leading edge  
* Ownership: Wing  
* Local reference examples  
  * Text: …relative to the leading edge root point…  
  * Variables: …\_Ler  
* Airplane-local reference examples  
  * Text: …relative to the first Wing's leading edge root point…  
  * Variables: …\_Ler1  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's leading edge root point…  
  * Variables: …\_Ler2P1

## 5. Leading point

* The leading point of the WingCrossSection  
* Ownership: WingCrossSection  
* Local reference examples  
  * Text: …relative to the leading point…  
  * Variables: …\_Lp  
* Wing-local reference examples  
  * Text: …relative to the first WingCrossSection's leading point…  
  * Variables: …\_Lp1  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's first WingCrossSection's leading point…  
  * Variables: …\_Lp1Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's first WingCrossSection's leading point…  
  * Variables: …\_Lp1Wn2P1

## 6. Leading point parent

* For a Wing's first WingCrossSection, this is the Wing's leading edge root point. For subsequent WingCrossSections, this is the previous WingCrossSection's leading point.  
* Ownership: WingCrossSection  
* Local reference examples  
  * Text: …relative to the leading point parent)  
  * Variables: …\_Lpp  
* Wing-local reference examples  
  * Text: …relative to the first WingCrossSection's leading point parent)  
  * Variables: …\_Lpp1  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's first WingCrossSection's leading point parent)  
  * Variables: …\_Lpp1Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's first WingCrossSection's leading parent point…  
  * Variables: …\_Lpp1Wn2P1

## 7. Panel points

* The front right, front left, back left, back right, and collocation points of a Panel  
* Ownership: Panel  
* Local reference examples  
  * Text: …relative to the panel front right point…  
  * Variables: …\_Frpp  
* Wing-local reference examples  
  * Text: …relative to the (3, 2\) Panel's front right point…  
  * Variables: …\_Frppr3c2  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's (3, 2\) Panel's front right point…  
  * Variables: …\_Frppr3c2Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's (3, 2\) Panel's front right point…  
  * Variables: …\_Frppr3c2Wn2P1

## 8. Bound horseshoe vortex points

* Only relevant in steady horseshoe vortex lattice method simulations  
* The front right, front left, back left, and back right points of a bound HorseshoeVortex  
* Ownership: HorseshoeVortex  
* Local reference examples  
  * Text: …relative to the bound horseshoe vortex front right point…  
  * Variables: …\_Frbhvp  
* Wing-local reference examples  
  * Text: …relative to the (3, 2\) Panel's bound HorseshoeVortex's front right point…  
  * Variables: …\_Frbhvpr3c2  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's (3, 2\) Panel's bound HorseshoeVortex's front right point…  
  * Variables: …\_Frbhvpr3c2Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's (3, 2\) Panel's bound HorseshoeVortex's front right point…  
  * Variables: …\_Frbhvpr3c2Wn2P1

## 9. Bound ring vortex points

* The front right, front left, back left, and back right points of a bound RingVortex  
* Ownership: RingVortex  
* Local reference examples  
  * Text: …relative to the bound ring vortex front right point…  
  * Variables: …\_Frbrvp  
* Wing-local reference examples  
  * Text: …relative to the (3, 2\) Panel's bound RingVortex's front right point…  
  * Variables: …\_Frbrvpr3c2  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's (3, 2\) Panel's bound RingVortex's front right point…  
  * Variables: …\_Frbrvpr3c2Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's (3, 2\) Panel's bound RingVortex's front right point…  
  * Variables: …\_Frbrvpr3c2Wn2P1

## 10. Wake horseshoe vortex points

* Only relevant in steady horseshoe and steady ring vortex lattice method simulations  
* The front right, front left, back left, and back right points of a wake HorseshoeVortex  
* Ownership: HorseshoeVortex  
* Local reference examples  
  * Text: …relative to the wake horseshoe vortex front right point…  
  * Variables: …\_Frwhvp  
* Wing-local reference examples  
  * Text: …relative to the third wake HorseshoeVortex's front right point…  
  * Variables: …\_Frwhvp3  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's third wake HorseshoeVortex's front right point…  
  * Variables: …\_Frwhvp3Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's third wake HorseshoeVortex's front right point…  
  * Variables: …\_Frwhvp3Wn2P1

## 11. Wake ring vortex points

* Only relevant in unsteady ring vortex lattice method simulations  
* The front right, front left, back left, and back right points of a wake RingVortex  
* Ownership: RingVortex  
* Local reference examples  
  * Text: …relative to the wake ring vortex front right point…  
  * Variables: …\_Frwrvp  
* Wing-local reference examples  
  * Text: …relative to the (3, 2\) wake RingVortex's front right point…  
  * Variables: …\_Frwrvpr3c2  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's (3, 2\) wake RingVortex's front right point…  
  * Variables: …\_Frwrvpr3c2Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's (3, 2\) wake RingVortex's front right point…  
  * Variables: …\_Frwrvpr3c2Wn2P1

## 12. Line vortex points

* The start, end, and center points of a LineVortex  
* Ownership: LineVortex  
* Local reference examples  
  * Text: …relative to the LineVortex start point…  
  * Variables: …\_Slvp  
* Parent-vortex-local reference examples  
  * Text: …relative to the bound HorseshoeVortex's front LineVortex's center point…  
  * Variables: …\_Clvpf  
* Wing-local reference examples  
  * Text: …relative to the (3, 2\) bound HorseshoeVortex's front LineVortex's center point…  
  * Variables: …\_ClvpfBhvr3c2  
* Airplane-local reference examples  
  * Text: …relative to the second Wing's (3, 2\) bound HorseshoeVortex's front LineVortex's center point…  
  * Variables: …\_ClvpfBhvr3c2Wn2  
* Non-local reference examples  
  * Text: …relative to the first Airplane's second Wing's (3, 2\) bound HorseshoeVortex's front LineVortex's center point…  
  * Variables: …\_ClvpfBhvr3c2Wn2P1

# Reference Frames

## 1. Earth reference frame

* Inertial  
* Attached rigidly to the Earth  
* Ownership: None  
* References  
  * Text: …observed from the Earth frame…  
  * Variables: …\_\_E

## 2. Body reference frame

* Non-inertial  
* Attached rigidly to the Airplane's body  
* Ownership: Airplane  
* Local reference examples  
  * Text: …observed from the body frame…  
  * Variables …\_\_B  
* Non-local reference examples  
  * Text: …observed from the second Airplane's body frame…  
  * Variables …\_\_BP2

## 3. Wing reference frame

* Non-inertial  
* Attached rigidly to the root of a Wing's leading edge  
* Ownership: Wing  
* Local reference examples  
  * Text: …observed from the wing frame…  
  * Variables …\_\_Wn  
* Airplane-local reference examples  
  * Text: …observed from the second Wing's frame…  
  * Variables …\_\_Wn2  
* Non-local reference examples  
  * Text: …observed from the fourth Airplane's second Wing's frame…  
  * Variables …\_\_Wn2P4

## 4. Wing cross section reference frame

* Non-inertial  
* Attached rigidly to the leading point of a WingCrossSection  
* Ownership: WingCrossSection  
* Local reference examples  
  * Text: …observed from the wing cross section frame…  
  * Variables …\_\_Wcs  
* Wing-local reference examples  
  * Text: …observed from the third WingCrossSection's frame…  
  * Variables …\_\_Wcs3  
* Airplane-local reference examples  
  * Text: …observed from the second Wing's third WingCrossSection's frame…  
  * Variables …\_\_Wcs3Wn2  
* Non-local reference examples  
  * Text: …observed from the fourth Airplane's second Wing's third WingCrossSection's frame…  
  * Variables …\_\_Wcs3Wn2P4

## 5. Wing cross section parent reference frame

* Non-inertial  
* Attached rigidly to the leading parent point of a WingCrossSection  
* Ownership: WingCrossSection  
* Local reference examples  
  * Text: …observed from the wing cross section parent frame…  
  * Variables …\_\_Wcsp  
* Wing-local reference examples  
  * Text: …observed from the third WingCrossSection's parent frame…  
  * Variables …\_\_Wcsp3  
* Airplane-local reference examples  
  * Text: …observed from the second Wing's third WingCrossSection's parent frame…  
  * Variables …\_\_Wcsp3Wn2  
* Non-local reference examples  
  * Text: …observed from the fourth Airplane's second Wing's third WingCrossSection's parent frame…  
  * Variables …\_\_Wcsp3Wn2P4
