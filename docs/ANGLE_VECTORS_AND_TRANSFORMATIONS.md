# Angle Vectors and Transformations

As discussed in [Axes, Points, and Frames](https://raw.githubusercontent.com/camUrban/PteraSoftware/feature/improved_geometry_definitions/docs/AXES_COORDINATES_AND_FRAMES.md), Ptera Software defines vector-valued quantities in a multitude of different axis systems and relative to different reference points. Therefore, we must be able to find a vectors equivalent representations in different axis systems (i.e. perform passive transformations). Also, as flapping-wing flight inherently involves lots of rotational motion, we must be able to rotate vectors within their current axis systems (i.e. perform active transformations). We accomplish both of these tasks using angle vectors, rotation matrices, and transformation matrices.

# Angle Vectors

Angle vectors contain three scalar angles. These angles can either represent the orientation of one axis system with respect to another (passive angle vectors), or they can be the angles that would like use to rotate a given vector within one axis system (active angle vectors).

## Sequence IDs and Names

For both active and passive angle vectors, Ptera Software only uses Tait-Bryan rotation sequences. Below are the acceptable sequence IDs and names, which we use for variable naming and for describing angle vectors in text.

* i123: intrinsic 1-2’-3”  
* i132: intrinsic 1-3’-2”  
* i213: intrinsic 2-1’-3”  
* i231: intrinsic 2-3’-1”  
* i312: intrinsic 3-1’-2”  
* i321: intrinsic 3-2’-1”  
* e123: extrinsic 1-2-3  
* e132: extrinsic 1-3-2  
* e213: extrinsic 2-1-3  
* e231: extrinsic 2-3-1  
* e312: extrinsic 3-1-2  
* e321: extrinsic 3-2-1

## Passive Angle Vectors
Unlike force or position vectors which have components in a single axis system, passive angle vectors inherently relate two different axis systems and therefore require special notation to specify both the source and target axes, as well as the rotation sequence convention.

### Specifying Passive Angle Vectors

Given the dual-axes nature of angle vectors, we denote them by appending information about both axis systems and the rotation sequence:

#### Variable Name Pattern

\[variable name, default to “angles”\]\_\[source axes ID\]\_to\_\[target axes ID\]\_\[sequence ID\]

#### Text Reference Pattern

"\[variable name, default to “angles”\] describing the orientation of \[target axes name\] relative to \[source axes name\] using an \[sequence name\] sequence"

### Passive Angle Vector Examples

#### Local reference examples

* Variables: angles\_E\_to\_B\_i321  
* Text: …angles describing the orientation of the Earth axes relative to the body axes using an intrinsic 3-2’-1” sequence…

#### Wing-local reference examples

* Variables: angles\_Wcs1\_to\_Wn\_i321  
* Text: …angles describing the orientation of the wing axes from the first wing cross section's axes using an intrinsic 3-2’-1” sequence…

#### Airplane-local reference examples

* Variables: angles\_Wn2\_to\_G\_i321  
* Text: …angles describing the orientation of the geometry axes from the second wing's axes using an intrinsic 3-2’-1” sequence…

#### Non-local reference examples

* Variables: angles\_BP1\_to\_E\_e123  
* Text: …angles describing the orientation of the Earth axes from the first airplane's body axes using an extrinsic 1-2-3 sequence…

## Active Angle Vectors
Active angle vectors give instructions for rotating a vector within its current axis system. Therefore, they don't require information about the particular axes, only the type or rotation and the sequence.

### Specifying Active Angle Vectors

#### Variable Name Pattern

\[variable name, default to “angles”\]\_act\_\[sequence ID\]

#### Text Reference Pattern

"\[variable name, default to “angles”\] for rotation using an \[sequence name\] sequence"

### Active Angle Vector Examples

* Variables: angles\_act\_i321  
* Text: …angles for rotation using an intrinsic 3-2’-1” sequence…

## Implementation Notes

1. Angle wrapping: All angles should be wrapped to the range (-180, 180\] for consistency  
2. Singularities: Different sequences experience gimbal lock at particular points  
3. Units: All angles are in degrees unless explicitly noted otherwise  
4. Intrinsic vs Extrinsic: Remember that intrinsic and extrinsic rotations are equivalent with the order reversed (e.g. 3-2’-1” is the same as 1-2-3)

# Rotation and Transformation Matrices

This section formalizes how Ptera Software represents and composes rotations and more general transformations between the many axis systems and reference points defined elsewhere in the guide. It builds on the notation for axes, points, and frames as well as the notation for passive and active angle vectors.

Like angle vectors, rotation and transformation matrices can either represent the position and orientation of one axis system relative to another (passive rotation and transformation matrices), or they can be a used to transform a vector within its current axes.

By convention, we treat vectors as column vectors and left-multiply by matrices for both active and passive transformations.

## Passive Matrices

### Overview of Passive Matrices 

* Ptera Software uses two passive matrix types:  
  * R_pas: 3×3 rotation matrices that relate the orientation of one axis system relative to another.  
  * T_pas: 4×4 transformation matrices in homogeneous coordinates that can encode the relative rotation, plus additional transformations like translations and reflections, between two axes.

### Passive Matrix Name Patterns

#### 3x3 Rotation Matrices

* R\_pas\_\[source axes ID\]\_to\_\[target axes ID\]  
* …rotation matrix R, which maps from \[source axes name\] to \[target axes name\]…

#### 4x4 General Transformation Matrices

* T\_pas\_\[source axes ID\]\_\[source point ID\]\_to\_\[target axes ID\]\_\[target point ID\]  
* …transformation matrix T, which maps in homogeneous coordinates from \[source axes name\] relative to \[source point ID\] to \[target axes name\] relative to \[target point ID\]…

#### Examples:

* R\_pas\_W\_to\_B: …rotation matrix R, which maps from wind axes to body axes…  
* T\_pas\_Wn\_Ler\_to\_G\_I: …which maps in homogeneous coordinates from wing axes relative to the leading edge root point to geometry axes relative to the simulation’s starting point…

See the section on angle vectors for examples that can be adapted to form text references and variable names for matrices in non-local contexts.

## Active Matrices

### Overview of Active Matrices 

* Ptera Software uses two active matrix types:  
  * R_act: 3×3 rotation matrices that are used to rotate a vector in its current axis system.  
  * T_act: 4×4 transformation matrices in homogeneous coordinates that can rotate a vector and also perform additional transformations like translations and reflections.

### Active Matrix Name Patterns

#### 3x3 Rotation Matrices

* [variable name]\_R\_act  
* …[variable name], a matrix for active rotations…

#### 4x4 General Transformation Matrices

* [variable name]\_T\_act  
* …[variable name], a matrix for active transformations in homogeneous coordinates…

#### Examples:

* mirror_T\_act: …mirror, a matrix for active transformations in homogeneous coordinates…  
* translate_wing\_T\_act: …translate_wing, a matrix for active transformations in homogeneous coordinates…

See the section on angle vectors for examples that can be adapted to form text references and variable names for matrices in non-local contexts.
