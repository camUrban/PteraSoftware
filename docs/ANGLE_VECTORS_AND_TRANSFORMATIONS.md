# Angle Vectors and Transformations

As discussed in [Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md), Ptera Software defines vector-valued quantities in a multitude of different axis systems and relative to different reference points. Therefore, we must be able to find a vector's equivalent representations in different axis systems (i.e. perform passive transformations). Also, as flapping-wing flight inherently involves lots of rotational motion, we must be able to rotate vectors within their current axis systems (i.e. perform active transformations). We accomplish both of these tasks using angle vectors, rotation matrices, and transformation matrices.

## Angle Vectors

Angle vectors contain three scalar angles. These angles can either represent the orientation of one axis system with respect to another (passive angle vectors), or they can be the angles that we would like to use to rotate a given vector within one axis system (active angle vectors). Angle vectors always take the form (angleX, angleY, angleZ).

### Sequence IDs and Names

For both active and passive angle vectors, Ptera Software only uses Tait-Bryan rotation sequences. Below are the acceptable sequence IDs and names, which we use for variable naming and for describing angle vectors in text.

* ixyz: intrinsic xy'z"
* ixzy: intrinsic x-z'-y"
* iyxz: intrinsic y-x'-z"
* iyzx: intrinsic y-z'-x"
* izxy: intrinsic z-x'-y"
* izyx: intrinsic zy'x"
* exyz: extrinsic xyz
* exzy: extrinsic x-z-y
* eyxz: extrinsic y-x-z
* eyzx: extrinsic y-z-x
* ezxy: extrinsic z-x-y
* ezyx: extrinsic z-y-x

### Passive Angle Vectors

Unlike force or position vectors which have components in a single axis system, passive angle vectors inherently relate two different axis systems and therefore require special notation to specify both the source and target axes, as well as the rotation sequence convention.

#### Specifying Passive Angle Vectors

Given the dual-axes nature of angle vectors, we denote them by appending information about both axis systems and the rotation sequence:

##### Variable Name Pattern

\[variable name, default to "angles"\]\_\[source axes ID\]\_to\_\[target axes ID\]\_\[sequence ID\]

##### Text Reference Pattern

"\[variable name, default to "angles"\] describing the orientation of \[target axes name\] relative to \[source axes name\] using an \[sequence name\] sequence"

#### Passive Angle Vector Examples

##### Local reference examples

* Variables: angles\_E\_to\_B\_izyx
* Text: ...angles describing the orientation of the body axes relative to the Earth axes using an intrinsic zy'x" sequence...

##### Wing-local reference examples

* Variables: angles\_Wcs1\_to\_Wn\_izyx
* Text: ...angles describing the orientation of the wing axes from the first WingCrossSection's axes using an intrinsic zy'x" sequence...

##### Airplane-local reference examples

* Variables: angles\_Wn2\_to\_G\_izyx
* Text: ...angles describing the orientation of the geometry axes from the second Wing's axes using an intrinsic zy'x" sequence...

##### Non-local reference examples

* Variables: angles\_BP1\_to\_E\_exyz
* Text: ...angles describing the orientation of the Earth axes from the first Airplane's body axes using an extrinsic xyz sequence...

### Active Angle Vectors

Active angle vectors give instructions for rotating a vector within its current axis system. Therefore, they don't require information about the particular axes, only the type or rotation and the sequence.

#### Specifying Active Angle Vectors

##### Variable Name Pattern

\[variable name, default to "angles"\]\_act\_\[sequence ID\]

##### Text Reference Pattern

"\[variable name, default to "angles"\] for rotation using an \[sequence name\] sequence"

#### Active Angle Vector Examples

* Variables: angles\_act\_izyx
* Text: ...angles for rotation using an intrinsic zy'x" sequence...

### Angle Vector Components

Appending a component letter (X, Y, or Z) to the end of the camelCase portion of an angle vector's variable name denotes a variable that holds only that scalar component (or a stack of that component), while keeping the rest of the name so the component stays tied to the angle vector it came from. Match the name's plurality to what the variable holds: a single component value singularizes the name, while a stack of component values keeps it plural. For example, angleY\_E\_to\_B\_izyx would hold only the y component of angles\_E\_to\_B\_izyx, and deformationAnglesYRad\_Wcsp\_to\_Wcs\_ixyz holds one y component per WingCrossSection. This is the same component letter convention that [Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md) defines for all vector-valued quantities, and the stack, grid, and list collection prefixes defined there apply to angle vectors in the same way (e.g., listDeformationAnglesYRad\_Wcsp\_to\_Wcs\_ixyz holds per-Wing lists of stored angle components).

### Angle Component Time Derivatives

Some quantities are time derivatives of an angle vector's individual components, such as the aeroelastic structural model's torsional angle rates and the movement classes' second-derivative spacing functions (e.g., spacingAnglesSecondDerivative\_Gs\_to\_Wn\_ixyz). These are rates of change of scalar coordinates.

An angle vector's componentwise time derivative is not, in general, an angular velocity vector. For Tait-Bryan sequences, the angular velocity vector relates to the component rates through a sequence-dependent kinematic matrix, and the two coincide only in special cases (for example, a rotation about a single fixed axis). Reserve omega-style names (e.g., omegas\_BP1\_\_E) for true angular velocity vectors, which take an axes ID and a frame ID like any other vector observed over time. Never name a stack of angle component derivatives as omegas, and never treat one as an angular velocity vector.

To name an angle component time derivative, start from the name of the angle vector being differentiated, add a derivative marker to the end of the camelCase portion (Derivative for a first derivative, SecondDerivative for a second), add a component letter (X, Y, or Z) when only some components are stored, and keep the source angle vector's full source-to-target and sequence suffix so the derivative is tied to the angle vector it differentiates. When the value is in radians-based units (rad/s, rad/s^2), mark it with Rad as the last element of the camelCase portion. For example, deformationAnglesDerivativeYRad\_Wcsp\_to\_Wcs\_ixyz holds the time derivatives (rad/s) of the y components of deformationAngles\_Wcsp\_to\_Wcs\_ixyz, angles describing the orientation of the wing cross section axes relative to the wing cross section parent axes using an intrinsic xy'z" sequence.

Angle component time derivatives take no reference frame ID. A frame ID records the observer when differentiating a vector, because the result depends on the rotation of the basis the vector is expressed in. A scalar angle component has no basis, so once the component is defined, its time derivative is unambiguous and no observer needs to be named.

### Implementation Notes

1. Angle wrapping: All angles should be wrapped to the range (-180, 180\] for consistency
2. Singularities: Different sequences experience gimbal lock at particular points
3. Units: For all quantities whose units include an angular component, that angular component is in degrees unless explicitly noted otherwise
4. Intrinsic vs Extrinsic: Remember that intrinsic and extrinsic rotations are equivalent with the order reversed (e.g. zy'x" is the same as xyz)
5. Radians naming: For vector-valued quantities and their components, noting radians-based units (rad, rad/s, rad/s^2, etc.) in a docstring or comment is not enough, so the variable's name itself must declare them, with Rad as the last element of its camelCase portion (e.g., omegasRad\_BP1\_\_E), while non-vector-valued variables instead take a \_rad suffix on their snake case names, as described in [Code Style](CODE_STYLE.md)

## Rotation and Transformation Matrices

This section formalizes how Ptera Software represents and composes rotations and more general transformations between the many axis systems and reference points defined elsewhere in the guide. It builds on the notation for axes, points, and frames as well as the notation for passive and active angle vectors.

Like angle vectors, rotation and transformation matrices can either represent the position and orientation of one axis system relative to another (passive rotation and transformation matrices), or be used to transform a vector within its current axes (active rotation and transformation matrices).

By convention, we treat vectors as column vectors and left-multiply by matrices for both active and passive transformations.

### Homogeneous Coordinates

We can use 3x3 rotation matrices to transform our vector-valued quantities as is. However, when working with 4x4 transformation matrices, we must convert vector's three components to four homogeneous coordinates before applying the transformation.

For vector with components (x, y, z) that is given relative to a reference point, its homogeneous form is (x, y, z, 1). If it is a vector independent of any reference point, its homogeneous form is (x, y, z, 0). After applying the transformation, we can convert either type back to non-homogeneous components by dropping the last coordinate.

### Passive Matrices

#### Overview of Passive Matrices

* Ptera Software uses two passive matrix types:
    * R\_pas\_...: 2x2 or 3x3 rotation matrices that relate the orientation of one axis system relative to another.
    * T\_pas\_...: 4x4 transformation matrices in homogeneous coordinates that maps components from a source axis system and reference point to a target axis system and reference point. It applies the orientation change implied by the two axes (rotation or reflection) and, when the quantity is tied to a reference point (e.g., positions), also applies the translation between the points. It only changes how the same physical quantity is expressed (axes and/or point) and never introduces scaling or shear.

#### Passive Matrix Name Patterns

##### 3x3 Rotation Matrices

* R\_pas\_\[source axes ID\]\_to\_\[target axes ID\]
* ...rotation matrix R, which maps from \[source axes name\] to \[target axes name\]...

##### 4x4 General Transformation Matrices

* T\_pas\_\[source axes ID\]\_\[source point ID\]\_to\_\[target axes ID\]\_\[target point ID\]
* ...transformation matrix T, which maps in homogeneous coordinates from \[source axes name\] relative to \[source point ID\] to \[target axes name\] relative to \[target point ID\]...

##### Examples:

* R\_pas\_W\_to\_B: ...rotation matrix R, which maps from wind axes to body axes...
* T\_pas\_Wn\_Ler\_to\_G\_I: ...which maps in homogeneous coordinates from wing axes relative to the leading edge root point to geometry axes relative to the simulation starting point...

See the section on angle vectors for examples that can be adapted to form text references and variable names for matrices in non-local contexts.

### Active Matrices

#### Overview of Active Matrices

* Ptera Software uses two active matrix types:
    * [variable name]\_R\_act: 2x2 or 3x3 rotation matrices that are used to rotate a vector in its current axis system.
    * [variable name]\_T\_act: 4x4 homogeneous transformation that operates within a single axis system. It applies a rigid orientation change (rotation or reflection) and, for quantities tied to a reference point, may also apply a translation. It never changes which axes a vector is expressed in and never introduces scaling or shear. The translation has no effect on non-position vectors (e.g., forces and moments).

#### Active Matrix Name Patterns

##### 3x3 Rotation Matrices

* [variable name]\_R\_act
* ...[variable name], a matrix for active rotations...

##### 4x4 General Transformation Matrices

* [variable name]\_T\_act
* ...[variable name], a matrix for active transformations in homogeneous coordinates...

##### Examples:

* mirror\_T\_act: ...mirror, a matrix for active transformations in homogeneous coordinates...
* translate\_wing\_T\_act: ...translate_wing, a matrix for active transformations in homogeneous coordinates...

See the section on angle vectors for examples that can be adapted to form text references and variable names for matrices in non-local contexts.
