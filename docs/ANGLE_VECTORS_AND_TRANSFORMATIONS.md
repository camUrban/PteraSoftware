# Angle Vectors and Transformations

As discussed in [Axes, Points, and Frames](https://raw.githubusercontent.com/camUrban/PteraSoftware/feature/improved_geometry_definitions/docs/AXES_POINTS_AND_FRAMES.md), Ptera Software defines vector-valued quantities in a multitude of different axis systems and relative to different reference points. Therefore, we must be able to find a vector's equivalent representations in different axis systems (i.e. perform passive transformations). Also, as flapping-wing flight inherently involves lots of rotational motion, we must be able to rotate vectors within their current axis systems (i.e. perform active transformations). We accomplish both of these tasks using angle vectors, rotation matrices, and transformation matrices.

# Angle Vectors

Angle vectors contain three scalar angles. These angles can either represent the orientation of one axis system with respect to another (passive angle vectors), or they can be the angles that we would like to use to rotate a given vector within one axis system (active angle vectors). Angle vectors always take the form (angleX, angleY, angleZ).

## Sequence IDs and Names

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

## Passive Angle Vectors
Unlike force or position vectors which have components in a single axis system, passive angle vectors inherently relate two different axis systems and therefore require special notation to specify both the source and target axes, as well as the rotation sequence convention.

### Specifying Passive Angle Vectors

Given the dual-axes nature of angle vectors, we denote them by appending information about both axis systems and the rotation sequence:

#### Variable Name Pattern

\[variable name, default to "angles"\]\_\[source axes ID\]\_to\_\[target axes ID\]\_\[sequence ID\]

#### Text Reference Pattern

"\[variable name, default to "angles"\] describing the orientation of \[target axes name\] relative to \[source axes name\] using an \[sequence name\] sequence"

### Passive Angle Vector Examples

#### Local reference examples

* Variables: angles\_E\_to\_B\_izyx  
* Text: …angles describing the orientation of the body axes relative to the Earth axes using an intrinsic zy'x" sequence…

#### Wing-local reference examples

* Variables: angles\_Wcs1\_to\_Wn\_izyx  
* Text: …angles describing the orientation of the wing axes from the first WingCrossSection's axes using an intrinsic zy'x" sequence…

#### Airplane-local reference examples

* Variables: angles\_Wn2\_to\_G\_izyx  
* Text: …angles describing the orientation of the geometry axes from the second Wing's axes using an intrinsic zy'x" sequence…

#### Non-local reference examples

* Variables: angles\_BP1\_to\_E\_exyz  
* Text: …angles describing the orientation of the Earth axes from the first Airplane's body axes using an extrinsic xyz sequence…

## Active Angle Vectors
Active angle vectors give instructions for rotating a vector within its current axis system. Therefore, they don't require information about the particular axes, only the type or rotation and the sequence.

### Specifying Active Angle Vectors

#### Variable Name Pattern

\[variable name, default to "angles"\]\_act\_\[sequence ID\]

#### Text Reference Pattern

"\[variable name, default to "angles"\] for rotation using an \[sequence name\] sequence"

### Active Angle Vector Examples

* Variables: angles\_act\_izyx  
* Text: …angles for rotation using an intrinsic zy'x" sequence…

## Implementation Notes

1. Angle wrapping: All angles should be wrapped to the range (-180, 180\] for consistency  
2. Singularities: Different sequences experience gimbal lock at particular points  
3. Units: All angles are in degrees unless explicitly noted otherwise  
4. Intrinsic vs Extrinsic: Remember that intrinsic and extrinsic rotations are equivalent with the order reversed (e.g. zy'x" is the same as xyz)

# Rotation and Transformation Matrices

This section formalizes how Ptera Software represents and composes rotations and more general transformations between the many axis systems and reference points defined elsewhere in the guide. It builds on the notation for axes, points, and frames as well as the notation for passive and active angle vectors.

Like angle vectors, rotation and transformation matrices can either represent the position and orientation of one axis system relative to another (passive rotation and transformation matrices), or be used to transform a vector within its current axes (active rotation and transformation matrices).

By convention, we treat vectors as column vectors and left-multiply by matrices for both active and passive transformations.

## Homogeneous Coordinates
We can use 3x3 rotation matrices to transform our vector-valued quantities as is. However, when working with 4x4 transformation matrices, we must convert vector's three components to four homogeneous coordinates before applying the transformation.

For vector with components (x, y, z) that is given relative to a reference point, its homogeneous form is (x, y, z, 1). If it is a vector independent of any reference point, its homogeneous form is (x, y, z, 0). After applying the transformation, we can convert either type back to non-homogeneous components by dropping the last coordinate.

## Passive Matrices

### Overview of Passive Matrices 

* Ptera Software uses two passive matrix types:  
  * R\_pas\_…: 3x3 rotation matrices that relate the orientation of one axis system relative to another.  
  * T\_pas\_…: 4x4 transformation matrices in homogeneous coordinates that maps components from a source axis system and reference point to a target axis system and reference point. It applies the orientation change implied by the two axes (rotation or reflection) and, when the quantity is tied to a reference point (e.g., positions), also applies the translation between the points. It only changes how the same physical quantity is expressed (axes and/or point) and never introduces scaling or shear.

### Passive Matrix Name Patterns

#### 3x3 Rotation Matrices

* R\_pas\_\[source axes ID\]\_to\_\[target axes ID\]  
* …rotation matrix R, which maps from \[source axes name\] to \[target axes name\]…

#### 4x4 General Transformation Matrices

* T\_pas\_\[source axes ID\]\_\[source point ID\]\_to\_\[target axes ID\]\_\[target point ID\]  
* …transformation matrix T, which maps in homogeneous coordinates from \[source axes name\] relative to \[source point ID\] to \[target axes name\] relative to \[target point ID\]…

#### Examples:

* R\_pas\_W\_to\_B: …rotation matrix R, which maps from wind axes to body axes…  
* T\_pas\_Wn\_Ler\_to\_G\_I: …which maps in homogeneous coordinates from wing axes relative to the leading edge root point to geometry axes relative to the simulation starting point…

See the section on angle vectors for examples that can be adapted to form text references and variable names for matrices in non-local contexts.

## Active Matrices

### Overview of Active Matrices 

* Ptera Software uses two active matrix types:  
  * [variable name]\_R\_act: 3x3 rotation matrices that are used to rotate a vector in its current axis system.  
  * [variable name]\_T\_act: 4x4 homogeneous transformation that operates within a single axis system. It applies a rigid orientation change (rotation or reflection) and, for quantities tied to a reference point, may also apply a translation. It never changes which axes a vector is expressed in and never introduces scaling or shear. The translation has no effect on free vectors (e.g., forces).

### Active Matrix Name Patterns

#### 3x3 Rotation Matrices

* [variable name]\_R\_act  
* …[variable name], a matrix for active rotations…

#### 4x4 General Transformation Matrices

* [variable name]\_T\_act  
* …[variable name], a matrix for active transformations in homogeneous coordinates…

#### Examples:

* mirror\_T\_act: …mirror, a matrix for active transformations in homogeneous coordinates…  
* translate\_wing\_T\_act: …translate_wing, a matrix for active transformations in homogeneous coordinates…

See the section on angle vectors for examples that can be adapted to form text references and variable names for matrices in non-local contexts.
