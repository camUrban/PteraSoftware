"""This module contains useful class definitions for different types of geometries.

This module contains the following classes:
    Airplane: This is a class used to contain airplanes.

    Wing: This is a class used to contain the wings of an Airplane.

    WingCrossSection: This class is used to contain the cross sections of a Wing.

    Airfoil: This class is used to contain the airfoil of a WingCrossSection.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import importlib.resources

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sp_interp

from . import functions
from . import meshing


class Airplane:
    """This is a class used to contain airplanes.

    Citation:
        Adapted from:         geometry.Airplane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/23/2020

    This class contains the following public methods:
        set_reference_dimensions_from_main_wing: This method sets the reference
        dimensions of the Airplane from measurements obtained from its first Wing.

        process_wing_symmetry: This method processes Wings with symmetric=True
        and converts them to separate Wings when needed (Scenario 5).

        validate_first_airplane_constraints: This method validates that the first
        Airplane in a simulation has Cgi_E_I set to zeros, as required by the
        definition of the simulation's starting point.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.

    The Airplane class serves as the root coordinate system for the aircraft geometry.
    All coordinate systems are ultimately defined relative to the Airplane's local
    geometry axes, which represent the aircraft's primary reference frame. The Airplane
    class is responsible for:

    1. Defining the local geometry axes as the root coordinate system
    2. Managing Wings and their coordinate transformations
    3. Processing symmetric Wings and converting them to separate wings when
       the symmetry plane is not coincident with the Wing's xz-plane (Scenario 5)
    4. Providing reference dimensions for aerodynamic calculations
    5. Managing the moment reference point for force and moment calculations

    Local geometry axes convention:
    - +x: Points aft along fuselage
    - +y: Points to the right (starboard direction)
    - +z: Points upward (completing right-handed coordinate system)
    """

    def __init__(
        self,
        wings,
        name="Untitled Airplane",
        Cgi_E_I=np.array([0.0, 0.0, 0.0]),
        angles_E_to_B_i321=np.array([0.0, 0.0, 0.0]),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    ):
        """This is the initialization method.

        :param wings: list of Wings
            This is a list of the airplane's wings defined as Wings. It must
            contain at least one Wing. Wings with symmetric=True and non-coincident
            symmetry planes will be automatically processed into separate Wings
            during initialization (Scenario 5).
        :param name: str, optional
            A sensible name for your airplane. The default is "Untitled Airplane".
        :param Cgi_E_I: (3,) ndarray of floats, optional

            Position of this Airplane's starting point (in Earth axes, relative to
            the simulation's starting point). For the first Airplane in a simulation,
            this must be np.array([0.0, 0.0, 0.0]) since the simulation's starting
            point is defined as the first Airplane's starting point (the location of
            its CG at t=0). The default is np.array([0.0, 0.0, 0.0]).

        :param angles_E_to_B_i321: (3,) ndarray of floats, optional

            Angles from Earth axes to body axes using an intrinsic 3-2'-1" sequence.
            This defines the orientation of the airplane's body axes relative to
            Earth axes. Note that body axes differ from geometry axes: body axes
            point forward/right/down while geometry axes point aft/right/up. The
            units are degrees. The first angle must lie in the range (-180,
            180] degrees, the second in [-60, 60] degrees, and the third in (-90,
            90) degrees. This is to reduce the chance of edge cases, and to eliminate
            the risk of gimbal lock. The default is np.array( [0.0, 0.0, 0.0]).

        :param weight: float, optional
            This parameter holds the weight of the aircraft in Newtons. This is used
            by the trim functions. It must be greater than or equal to zero. The
            default value is 0.0.
        :param s_ref: float, optional if more than one Wing is in the wings list.
            This is the reference wetted area. If not set, it populates from first
            Wing. The units are square meters.
        :param c_ref: float, optional if more than one Wing is in the wings list.
            This is the reference chord length. If not set, it populates from first
            Wing. The units are meters.
        :param b_ref: float, optional if more than one Wing is in the wings list.
            This is the reference span. If not set, it populates from first
            Wing. The units are meters.
        """
        # Initialize the list of wings or raise an exception if it is empty.
        if len(wings) > 0:
            self.wings = wings
        else:
            raise Exception("An Airplane's list of wings must have at least one entry.")

        # Initialize the name, starting point coordinates, body orientation, and weight.
        self.name = name
        self.Cgi_E_I = np.array(Cgi_E_I, dtype=float)

        # ToDo: Implement validation methods for these two parameters.
        self.anglesi_E_to_B_i321 = np.array(angles_E_to_B_i321, dtype=float)
        self.weight = weight

        # Set the reference dimensions to be the first Wing's reference dimensions.
        self.set_reference_dimensions_from_main_wing()

        # If any of the passed reference dimensions are not None, set that reference
        # dimension to be what was passed.
        if s_ref is not None:
            self.s_ref = s_ref
        if c_ref is not None:
            self.c_ref = c_ref
        if b_ref is not None:
            self.b_ref = b_ref

        # Calculate the total number of panels of all the Wings in this Airplane.
        self.num_panels = 0
        for wing_position, wing in enumerate(self.wings):
            self.num_panels += wing.num_panels

        # Initialize empty class attributes to hold the force, moment,
        # force coefficients, and moment coefficients this Airplane experiences.
        self.total_near_field_force_wind_axes = None
        self.total_near_field_force_coefficients_wind_axes = None
        self.total_near_field_moment_wind_axes = None
        self.total_near_field_moment_coefficients_wind_axes = None

    def set_reference_dimensions_from_main_wing(self):
        """This method sets the reference dimensions of the Airplane from
        measurements obtained from the first Wing.

        :return: None
        """

        main_wing = self.wings[0]

        # Set the object's reference dimension attributes to be the reference
        # dimension attributes of the first Wing. These attributes are calculated via
        # methods in the Wing class.
        self.s_ref = main_wing.projected_area
        self.b_ref = main_wing.span
        self.c_ref = main_wing.mean_aerodynamic_chord

    def validate_first_airplane_constraints(self):
        """This method validates constraints specific to the first Airplane in a simulation.

        The first Airplane in a simulation must have Cgi_E_I set to zeros since the
        simulation starting point is defined as the first Airplane's starting point.

        This method should be called by SteadyProblem or UnsteadyProblem classes.

        :raises Exception: If first Airplane constraints are violated.
        """
        if not np.allclose(self.Cgi_E_I, np.array([0.0, 0.0, 0.0])):
            raise Exception(
                "The first Airplane in a simulation must have Cgi_E_I set to "
                "np.array([0.0, 0.0, 0.0]) since the simulation starting point "
                "is defined as the first Airplane's CG at t=0."
            )

    def process_wing_symmetry(self):
        """This method processes Wings with symmetric=True for Scenario 5 conversion.

        For Wings with symmetric=True and symmetry planes not coincident with the Wing's
        xz-plane, this method:
        1. Modifies the original Wing to become a "Scenario 1 Wing"
        2. Creates a new reflected Wing as a "Scenario 3 Wing"
        3. Adds the reflected Wing to the wings list immediately after the original
        4. Handles asymmetric control surface deflection sign flipping

        :return: None
        """
        # TODO: Implement Scenario 5 Wing processing
        # This will involve:
        # 1. Identify Wings with symmetric=True and non-coincident symmetry planes
        # 2. Modify original Wing parameters (set symmetric=False, etc.)
        # 3. Create reflected Wing with mirror_only=True
        # 4. Handle control surface deflection sign flipping for asymmetric surfaces
        # 5. Insert reflected Wing into wings list
        pass


class Wing:
    """This is a class used to contain the wings of an Airplane.

    Citation:
        Adapted from:         geometry.Wing in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/24/2020

    This class contains the following public methods:
        unit_up_vector: This method sets a property for the Wing's up orientation
        vector, which is defined as the cross product of its unit chordwise and unit
        normal vectors.

        projected_area: This method defines a property for the area of the Wing
        projected onto the plane defined by the projected unit normal vector.

        wetted_area: This method defines a property for the Wing's wetted area.

        span: This method defines a property for the Wing's span.

        standard_mean_chord: This method calculates the standard mean chord of the
        Wing and assigns it to the standard_mean_chord attribute.

        mean_aerodynamic_chord: This method calculates the mean aerodynamic chord of
        the Wing and assigns it to the mean_aerodynamic_chord attribute.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.

    Every Wing has its axis system, known as wing axes. The user sets the
    relationship between these axes and geometry axes with the local_position and
    local_rotations parameters. However, the steps for transforming a vector from
    geometry axes to wing axes, and the interpretation of the wing axes orientation
    and position relative to a Wing's geometry, also depend on the parameters
    symmetric, mirror_only, symmetry_normal_G, and symmetry_point_G_Cg:
        1. symmetric is False
            A. mirror_only is False
                I. The symmetry plane must be undefined (symmetry_normal_G and
                symmetry_point_G_Cg must be None)
                    Scenario 1:
                    - local_position is the final location of leading edge of this
                    Wing's root WingCrossSection, as defined in geometry axes.
                    - local_position is also the final location of the origin of this
                    Wing's axes, as defined in geometry axes.
                    - Translation by local_position followed by rotations by
                    local_rotations fully define this Wing's axes with respect to the
                    geometry axes. The wing axes will also retain the handedness of
                    the geometry axes.
            B. mirror_only is True
                I. The symmetry plane is coincident with this Wing's axes' xz-plane
                    Scenario 2:
                    - local_position is the final location of leading edge of this
                    Wing's root WingCrossSection, as defined in geometry axes.
                    - local_position is also the final location of the origin of this
                    Wing's axes, as defined in geometry axes.
                    - Translation by local_position followed by rotations by
                    local_rotations does not fully define orientation of this Wing's
                    axes with respect to the geometry axes. After translation and
                    rotation, the coordinate system also needs to be reflected across
                    the symmetry plane, which will flip the wing axes' handedness to
                    be opposite that of geometry axes.
                II. The symmetry plane is not coincident with this Wing's axes'
                xz-plane Scenario 3:
                    - local_position is not final location of leading edge of this
                    Wing's root WingCrossSection, as defined in geometry axes.
                    - local_position is not the final location of the origin of this
                    Wing's axes, as defined in geometry axes.
                    - Translation by local_position followed by rotations by
                    local_rotations does not fully define orientation of this Wing's
                    axes with respect to the geometry axes. After translation and
                    rotation, the coordinate system also needs to be reflected across
                    the symmetry plane, which will flip the wing axes' handedness to
                    be opposite that of geometry axes.
        2. symmetric is True
            A. mirror_only must be False
                I. the symmetry plane is coincident with this Wing's axes' xz-plane
                    Scenario 4:
                    - local_position is the final location of leading edge of this
                    Wing's root WingCrossSection, as defined in geometry axes.
                    However, while the root WingCrossSection is the still the first
                    item in the wing_cross_sections list, when meshed, panels will
                    extend from the root in both the +y and -y wing axis directions.
                    The length of the wing_cross_sections list remains unchanged.
                    - local_position is also the final location of the origin of this
                    Wing's axes, as defined in geometry axes.
                    - Translation by local_position followed by rotations by
                    local_rotations fully define this Wing's axes with respect to the
                    geometry axes. The wing axes will also retain the handedness of
                    the geometry axes.
                II. the symmetry plane is not coincident with this Wing's axes' xz-plane
                    Scenario 5:
                    - This Wing's Airplane will set this Wing's symmetric parameter
                    to False, its mirror_only parameter to False,
                    its symmetry_normal_G parameter to None and its
                    symmetry_point_G_Cg parameter to None. These changes turn this
                    Wing into a "Scenario 1 Wing."
                    - The Airplane will also create a new Wing, and add it to its
                    wings list immediately after the current Wing. The new Wing will
                    have the same name as this Wing, but with the prefix "Reflected
                    ". The new Wing also will have all the same parameters as this
                    Wing, except that symmetric will be False and mirror_only will
                    be True, which means that it will be a "Scenario 3 Wing."
                    - Also, if the control_surface_type is "asymmetric" for any of
                    this Wing's WingCrossSections, the reflected Wing's corresponding
                    WingCrossSections will have their control_surface_deflection
                    values multiplied by negative one.
    """

    def __init__(
        self,
        wing_cross_sections,
        name="Untitled Wing",
        prelim_Ler_G_Cg=np.array([0.0, 0.0, 0.0]),
        angles_G_to_prelim_Wn=np.array([0.0, 0.0, 0.0]),
        symmetric=False,
        mirror_only=False,
        symmetry_normal_G=None,
        symmetry_point_G_Cg=None,
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    ):
        """This is the initialization method.

        :param wing_cross_sections: list of WingCrossSections
            This is a list of WingCrossSections that represent the wing's cross
            sections in order from root to tip. It must contain at least two
            WingCrossSections.

        :param name: str, optional
            This is a sensible name for the Wing. The default is "Untitled Wing".
        :param prelim_Ler_G_Cg: (3,) ndarray of floats, optional
            This is the position [x, y, z] of the origin of this Wing's axes (in
            geometry axes, relative to the starting point) before any symmetry or
            mirror has been applied. It may differ from the actual position as
            explained in the class docstring. The units are meters. The default is
            np.array([0.0, 0.0, 0.0]).
        :param angles_G_to_prelim_Wn: (3,) ndarray of floats, optional
            This is the rotation angles [roll, pitch, yaw] in degrees that define the
            orientation of this Wing's axes relative to the geometry axes before any
            symmetry or mirror has been applied. All angles must be in the range (
            -90, 90) degrees. Roll is rotation about the x-axis, pitch is rotation
            about the y-axis, and yaw is rotation about the z-axis. Rotations are
            intrinsic, and proceed in the z-y'-x'' order conventional for Euler
            angles. It may differ from the actual position as explained in the class
            docstring. The units are meters. The default is np.array([0.0, 0.0,
            0.0]). The units are degrees. The default is np.array([0.0, 0.0, 0.0]).
        :param symmetric: bool, optional
            Set this to True if the Wing's geometry should be mirrored across the
            symmetry plane while retaining the non-mirrored side. If mirror_only is
            True, symmetric must be False. If symmetric is true, then neither
            symmetry_normal_G nor symmetry_point_G_Cg can be None. If the symmetry
            plane is coincident with this Wing's wing axes' xz-plane, the mirrored
            and non-mirrored geometry will be meshed as a single wing. If not,
            this Wing's Airplane will automatically create another Wing with the
            mirrored geometry, modify both Wings' parameters, and add the reflected
            Wing to its list of wings immediately following this one. For more
            details on how that process, and how this parameter interacts with
            symmetry_normal_G, symmetry_point_G_Cg, and mirror_only, see the class
            docstring. The default is False.
        :param mirror_only: bool, optional
            Set this to True if the Wing's geometry should be reflected about the
            symmetry plane without retaining the non-reflected geometry. If symmetric
            is True, mirror_only must be False. If mirror_only is true, then neither
            symmetry_normal_G nor symmetry_point_G_Cg can be None. For more details
            on how this parameter interacts with symmetry_normal_G,
            symmetry_point_G_Cg, and symmetric, see the class docstring. The default
            is False.
        :param symmetry_normal_G: (3,) ndarray of floats or None, optional
            The unit normal vector (in geometry axes) that, together with
            symmetry_point_G_Cg, defines the plane used for symmetry or mirroring.
            Note that reversing the normal direction (using the antiparallel vector)
            defines the same plane and produces the same result. This value must be
            None if both symmetric and mirror_only are False, and cannot be None if
            either are True. For more details on how this parameter interacts with
            symmetry_point_G_Cg, symmetric, and mirror_only, see the class docstring.
            The default is None.
        :param symmetry_point_G_Cg: (3,) ndarray of floats or None, optional
            A point (in geometry axes, relative to the starting CG point) that,
            along with symmetry_normal_G, defines the location of the plane about
            which symmetry or mirroring is applied. This value must be None if both
            symmetric and mirror_only are False, and cannot be None if either are
            True. For more details on how this parameter interacts with
            symmetry_normal_G, symmetric, and mirror_only, see the class docstring.
            The units are meters. The default is None.
        :param num_chordwise_panels: int, optional
            This is the number of chordwise panels to be used on this wing,
            which must be set to a positive integer. The default is 8.
        :param chordwise_spacing: str, optional
            This is the type of spacing between the wing's chordwise panels. It can
            be "cosine" or "uniform". Using cosine spacing is highly recommended for
            steady simulations and uniform spacing is highly recommended for unsteady
            simulations. The default is "cosine".
        """
        # Initialize the list of wing cross sections or raise an exception if it
        # contains less than two entries.
        if len(wing_cross_sections) >= 2:
            self.wing_cross_sections = wing_cross_sections
        else:
            raise Exception(
                "A wing's list of wing cross sections must have at least "
                "two entries."
            )

        # Initialize the wing's attributes.
        self.name = name
        self.prelim_Ler_G_Cg = np.array(prelim_Ler_G_Cg, dtype=float)
        self.angles_G_to_prelim_Wn = np.array(angles_G_to_prelim_Wn, dtype=float)
        self.symmetric = symmetric
        self.mirror_only = mirror_only
        self.symmetry_normal_G = symmetry_normal_G
        self.symmetry_point_G_Cg = symmetry_point_G_Cg
        self.num_chordwise_panels = num_chordwise_panels

        # If the value for the chordwise spacing is valid, initialize it. Otherwise,
        # raise an exception.
        if chordwise_spacing in ["cosine", "uniform"]:
            self.chordwise_spacing = chordwise_spacing
        else:
            raise Exception('The chordwise spacing must be "cosine" or "uniform".')

        # Catch invalid preliminary rotation angles.
        if not np.all(
            (-90.0 < self.angles_G_to_prelim_Wn) & (self.angles_G_to_prelim_Wn < 90.0)
        ):
            raise Exception(
                "All local rotation angles must be in the range (-90, 90) degrees."
            )

        # Add mutual exclusivity validation for symmetric and mirror_only.
        if self.symmetric and self.mirror_only:
            raise Exception("symmetric and mirror_only cannot both be True.")

        # Validate symmetry plane parameters based on symmetric/mirror_only flags.
        if self.symmetric or self.mirror_only:
            if self.symmetry_normal_G is None:
                raise Exception(
                    "symmetry_normal_G cannot be None when symmetric or mirror_only is True."
                )
            if self.symmetry_point_G_Cg is None:
                raise Exception(
                    "symmetry_point_G_Cg cannot be None when symmetric or mirror_only is True."
                )
            # Convert to numpy arrays and validate unit vector.
            self.symmetry_normal_G = np.array(self.symmetry_normal_G, dtype=float)
            self.symmetry_point_G_Cg = np.array(self.symmetry_point_G_Cg, dtype=float)
            if not np.isclose(np.linalg.norm(self.symmetry_normal_G), 1.0):
                raise Exception("The symmetry_normal_G must be a unit vector.")
        else:
            if self.symmetry_normal_G is not None:
                raise Exception(
                    "symmetry_normal_G must be None when both symmetric and mirror_only are False."
                )
            if self.symmetry_point_G_Cg is not None:
                raise Exception(
                    "symmetry_point_G_Cg must be None when both symmetric and mirror_only are False."
                )

        # Validate root WingCrossSection constraints.
        self.wing_cross_sections[0].validate_root_constraints()

        # Validate tip WingCrossSection constraints.
        self.wing_cross_sections[-1].validate_tip_constraints()

        # TODO: Determine which reflection "Wing Scenario" applies

        # TODO: Uncomment this once we've determined which "Wing Scenario" applies
        # # Find the number of spanwise panels on the wing by adding each cross
        # # section's number of spanwise panels. Exclude the last cross section's
        # # number of spanwise panels as this is irrelevant. If the wing is symmetric,
        # # multiply the summation by two.
        # self.num_spanwise_panels = 0
        # for wing_cross_section in self.wing_cross_sections[:-1]:
        #     self.num_spanwise_panels += wing_cross_section.num_spanwise_panels
        # if self.symmetric:
        #     self.num_spanwise_panels *= 2
        #
        # # Calculate the number of panels on this wing.
        # self.num_panels = self.num_spanwise_panels * self.num_chordwise_panels

        # TODO: Change these if we are a "Scenario 5" Wing.
        self.Ler_G_Cg = self.prelim_Ler_G_Cg
        self.angles_G_to_Wn = self.angles_G_to_prelim_Wn

        # TODO: Check that meshing works correctly
        # # Initialize the panels attribute. Then mesh the wing, which will populate
        # # this attribute.
        # self.panels = None
        # meshing.mesh_wing(self)

        # Initialize an empty array to hold this wing's wake ring vortices and its
        # wake ring vortex vertices.
        self.wake_ring_vortex_vertices = np.empty((0, self.num_spanwise_panels + 1, 3))
        self.wake_ring_vortices = np.zeros((0, self.num_spanwise_panels), dtype=object)

    @property
    def T_G_Cg_to_Wn_Ler(self):
        """This method makes the transformation matrix T, which maps in homogeneous
        coordinates from geometry axes relative to the CG point to wing axes relative
        to the leading edge root point.

        The transformation applies the 5 scenarios described in the class docstring.

        :return: (4,4) ndarray of floats
            4x4 transformation matrix
        """

        # Step 1: Create T_G_Cg_to_G_prelim_Ler = [I, prelim_Ler_G_Cg; 0, 1]. This
        # matrix maps in homogeneous coordinates from geometry axes relative to the
        # CG point to geometry axes relative to the preliminary (i.e. non-reflected)
        # leading edge root point. This is the translation step.
        T_G_Cg_to_G_prelim_Ler = np.eye(4)
        T_G_Cg_to_G_prelim_Ler[:3, 3] = self.prelim_Ler_G_Cg

        # Step 2: Create R_G_to_prelim_Wn, which maps from which maps from geometry
        # axes to preliminary (i.e. non-reflected) wing axes. It will be formed from
        # using the angles from geometry axes to preliminary wing axes, applied
        # intrinsically in the order z-y'-x".
        roll_rad, pitch_rad, yaw_rad = np.radians(self.angles_G_to_prelim_Wn)

        R_z = np.array(
            [
                [np.cos(yaw_rad), -np.sin(yaw_rad), 0],
                [np.sin(yaw_rad), np.cos(yaw_rad), 0],
                [0, 0, 1],
            ]
        )
        R_y = np.array(
            [
                [np.cos(pitch_rad), 0, np.sin(pitch_rad)],
                [0, 1, 0],
                [-np.sin(pitch_rad), 0, np.cos(pitch_rad)],
            ]
        )
        R_x = np.array(
            [
                [1, 0, 0],
                [0, np.cos(roll_rad), -np.sin(roll_rad)],
                [0, np.sin(roll_rad), np.cos(roll_rad)],
            ]
        )

        # Combined rotation: R_G_to_prelim_Wn = R_x * R_y * R_z (for intrinsic
        # rotations)
        R_G_to_prelim_Wn = R_x @ R_y @ R_z

        # Step 3: Create T_G_to_prelim_Wn = [R_G_to_prelim_Wn, 0; 0, 1]. This matrix
        # maps in homogeneous coordinates from geometry axes to preliminary (i.e.
        # non-reflected) wing axes.
        T_G_to_prelim_Wn = np.eye(4)
        T_G_to_prelim_Wn[:3, :3] = R_G_to_prelim_Wn

        # Step 3: Create reflection matrix T_prelim_Wn_Ler_to_Wn_Ler = [H, c; 0,
        # 1] if needed
        if self.mirror_only:
            # Reflection across arbitrary plane defined by normal n and point P
            # Using H = I - 2*n@n^T and c = 2*n@n^T@P
            n = self.symmetry_normal_G  # Already normalized unit vector
            P = self.symmetry_point_G_Cg

            # Create reflection matrix H = I - 2*n@n^T
            I = np.eye(3)
            n_outer = np.outer(n, n)  # n @ n^T
            H = I - 2 * n_outer

            # Calculate translation component c = 2*n@n^T@P
            c = 2 * n_outer @ P

            # Construct 4x4 reflection matrix T_reflect = [H, c; 0, 1]
            T_reflect = np.eye(4)
            T_reflect[:3, :3] = H
            T_reflect[:3, 3] = c
        else:
            T_reflect = np.eye(4)

        # Step 4: Combine transformations: T_G_Cg_to_Wn_Ler = T_reflect @ T_rotate @ T_translate
        return T_reflect @ prelim_T_G_to_Wn @ prelim_T_G_Cg_to_G_Ler

    @property
    def unit_chordwise_vector(self):
        """This method sets a property for the wing's chordwise (x-axis) orientation
        vector in geometry axes.

        :return: (3,) ndarray of floats
            This is the wing's unit chordwise vector in geometry axes.
        """
        # Wing x-axis in wing coordinates is [1, 0, 0, 0] in homogeneous coordinates
        wing_x_axis = np.array([1.0, 0.0, 0.0, 0.0])
        # Transform to geometry coordinates (only rotation, not translation)
        T = self.geometry_to_wing_axes_transformation_matrix
        geometry_vector = T @ wing_x_axis
        return geometry_vector[:3]  # Return only xyz components

    @property
    def unit_spanwise_vector(self):
        """This method sets a property for the wing's spanwise (y-axis) orientation
        vector in geometry axes.

        :return: (3,) ndarray of floats
            This is the wing's unit spanwise vector in geometry axes.
        """
        # Wing y-axis in wing coordinates is [0, 1, 0, 0] in homogeneous coordinates
        wing_y_axis = np.array([0.0, 1.0, 0.0, 0.0])
        # Transform to geometry coordinates (only rotation, not translation)
        T = self.geometry_to_wing_axes_transformation_matrix
        geometry_vector = T @ wing_y_axis
        return geometry_vector[:3]  # Return only xyz components

    @property
    def unit_up_vector(self):
        """This method sets a property for the wing's up (z-axis) orientation
        vector in geometry axes.

        :return: (3,) ndarray of floats
            This is the wing's unit up vector in geometry axes.
        """
        # Wing z-axis in wing coordinates is [0, 0, 1, 0] in homogeneous coordinates
        wing_z_axis = np.array([0.0, 0.0, 1.0, 0.0])
        # Transform to geometry coordinates (only rotation, not translation)
        T = self.geometry_to_wing_axes_transformation_matrix
        geometry_vector = T @ wing_z_axis
        return geometry_vector[:3]  # Return only xyz components

    @property
    def projected_area(self):
        """This method defines a property for the area of the wing projected onto the
        plane defined by the projected unit normal vector.

        If the wing is symmetric, the area of the mirrored half is included.

        :return projected_area: float
            This attribute is the projected area of the wing. It has units of square
            meters.
        """
        projected_area = 0

        # Iterate through the chordwise and spanwise indices of the panels and add
        # their area to the total projected area.
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                projected_area += self.panels[
                    chordwise_location, spanwise_location
                ].calculate_projected_area(self.unit_up_vector)

        return projected_area

    @property
    def wetted_area(self):
        """This method defines a property for the wing's wetted area.

        If the wing is symmetrical, the area of the mirrored half is included.

        :return wetted_area: float
            This attribute is the wetted area of the wing. It has units of square
            meters.
        """
        wetted_area = 0

        # Iterate through the chordwise and spanwise indices of the panels and add
        # their area to the total wetted area.
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                wetted_area += self.panels[chordwise_location, spanwise_location].area

        return wetted_area

    @property
    def span(self):
        """This method defines a property for the wing's span.

        The span is found by first finding vector connecting the leading edges of the
        root and tip wing cross sections. Then, this vector is projected onto the
        symmetry plane's unit normal vector. The span is defined as the magnitude of
        this projection. If the wing is symmetrical, this method includes the span of
        the mirrored half.

        :return span: float
            This is the wing's span. It has units of meters.
        """
        root_to_tip_leading_edge = (
            self.wing_cross_sections[-1].leading_edge
            - self.wing_cross_sections[0].leading_edge
        )

        projected_leading_edge = (
            np.dot(root_to_tip_leading_edge, self.unit_normal_vector)
            * self.unit_normal_vector
        )

        span = np.linalg.norm(projected_leading_edge)

        # If the wing is symmetric, multiply the span by two.
        if self.symmetric:
            span *= 2

        return span

    @property
    def standard_mean_chord(self):
        """This method calculates the standard mean chord of the wing and assigns it
        to the standard_mean_chord attribute. The standard mean chord is defined as
        the projected area divided by the span. See their respective methods for the
        definitions of span and projected area.

        :return: float
            This is the standard mean chord of the wing. It has units of meters.
        """
        return self.projected_area / self.span

    @property
    def mean_aerodynamic_chord(self):
        """This method calculates the mean aerodynamic chord of the wing and assigns
        it to the mean_aerodynamic_chord attribute.

        :return: float
            This is the mean aerodynamic chord of the wing. It has units of meters.
        """
        # This method is based on the equation for the mean aerodynamic chord of a
        # wing, which can be found here:
        # https://en.wikipedia.org/wiki/Chord_(aeronautics)#Mean_aerodynamic_chord.
        # This equation integrates the squared chord from the wing center to the wing
        # tip. We will perform this integral piecewise for each section of the wing.
        integral = 0

        # Iterate through the wing cross sections to add the contribution of their
        # corresponding wing section to the piecewise integral.
        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections[:-1]
        ):
            next_wing_cross_section = self.wing_cross_sections[
                wing_cross_section_id + 1
            ]

            root_chord = wing_cross_section.chord
            tip_chord = next_wing_cross_section.chord

            # Find this section's span by following the same procedure as for the
            # overall wing span.
            section_leading_edge = (
                next_wing_cross_section.leading_edge - wing_cross_section.leading_edge
            )

            projected_section_leading_edge = (
                np.dot(section_leading_edge, self.unit_normal_vector)
                * self.unit_normal_vector
            )

            section_span = np.linalg.norm(projected_section_leading_edge)

            # Each wing section is, by definition, trapezoidal (at least when
            # projected on to the wing's projection plane). For a trapezoid,
            # the integral from the cited equation can be shown to evaluate to the
            # following.
            integral += (
                section_span
                * (root_chord**2 + root_chord * tip_chord + tip_chord**2)
                / 3
            )

        # Multiply the integral's value by the coefficients from the cited equation.
        if self.symmetric:
            return 2 * integral / self.projected_area
        return integral / self.projected_area


class WingCrossSection:
    """This class is used to contain the cross sections of a Wing.

    Citation:
        Adapted from:         geometry.WingXSec in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/26/2020

    This class contains the following public methods:
        validate_root_constraints: This method validates constraints specific to
        root WingCrossSections.

        validate_tip_constraints: This method validates constraints specific to
        tip WingCrossSections.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.

    The first WingCrossSection in a Wing's wing_cross_section list is known as the
    root WingCrossSection. The last is known as the tip WingCrossSection.

    Every WingCrossSection has its own wing cross section axes. For root
    WingCrossSections, their wing cross section axes are identical in position,
    orientation, and handedness to their Wing's wing axes. For all other
    WingCrossSections, their wing cross section axes are defined relative to the axes
    of the previous WingCrossSection. Locally, the x-axis points from a cross
    section's leading point to its trailing point, the y-axis points spanwise in the
    general direction of the next WingCrossSection, and the z-axis points upwards.

    Things can get a little confusing with respect to WingCrossSections for Wings
    with symmetric or mirror_only set to True. For more details, look in the Wing
    class's docstring. Also remember that WingCrossSections themselves aren't used
    for any simulations, they are merely one of the Wing attributes that help explain
    how to the meshing function how we'd like to generate its Panels.
    """

    def __init__(
        self,
        airfoil,
        num_spanwise_panels,
        chord=1.0,
        local_position=np.array([0.0, 0.0, 0.0]),
        local_rotations=np.array([0.0, 0.0, 0.0]),
        control_surface_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing="cosine",
    ):
        """This is the initialization method.

        :param airfoil: Airfoil
            This is the Airfoil to be used at this WingCrossSection.
        :param num_spanwise_panels: int or None
            This is the number of spanwise panels to be used between this
            WingCrossSection and the next one. For tip WingCrossSections,
            this must be None. For all other WingCrossSections, this must be a
            positive integer.
        :param chord: float, optional
            This is the chord of the wing at this WingCrossSection. The units are
            meters. It must be greater than 0.0. The default value is 1.0.
        :param local_position: (3,) ndarray of floats, optional
            This is the position [x, y, z] of this WingCrossSection's leading edge in
            the previous WingCrossSection's wing cross section axes. This is also the
            position of the origin of this WingCrossSection's wing cross section axes
            in the previous WingCrossSection's wing cross section axes. If this is
            the root WingCrossSection, this is relative to its Wing's axes and
            must be np.array([0.0, 0.0, 0.0]). All components must be non-negative.
            The units are meters. The default is np.array([0.0, 0.0, 0.0]).
        :param local_rotations: (3,) ndarray of floats, optional
            This is the rotation angles [roll, pitch, yaw] in degrees that define the
            orientation of this WingCrossSection's wing cross section axes relative
            to the previous WingCrossSection's wing cross section axes. For the root
            WingCrossSection, the orientation is relative to its Wing's wing
            axes and must be np.array([0.0, 0.0, 0.0]). All angles must be in the
            range (-90, 90) degrees. Roll is rotation about the x-axis, pitch is rotation about
            the y-axis, and yaw is rotation about the z-axis. Rotations are intrinsic, and proceed in
            the z-y'-x'' order conventional for Euler angles. The units are degrees.
            The default is np.array([0.0, 0.0, 0.0]).
        :param control_surface_type: str, optional
            This is type of control surfaces for this WingCrossSection. It can
            be "symmetric" or "asymmetric". An example of symmetric control surfaces
            are flaps. An example of asymmetric control surfaces are ailerons. The
            default value is "symmetric".
        :param control_surface_hinge_point: float, optional
            This is the location of the control surface hinge from the leading edge
            as a fraction of chord. It must be on the range (0.0, 1.0). The default
            value is 0.75.
        :param control_surface_deflection: float, optional
            This is the control deflection in degrees. Deflection downwards is
            positive. Must be in the range (-90, 90) degrees. The default value is
            0.0 degrees.
        :param spanwise_spacing: str, optional
            This can be "cosine" or "uniform". Using cosine spacing is highly
            recommended. The default value is "cosine".
        """

        # Initialize all the user-provided attributes.
        self.chord = chord
        self.local_position = np.array(local_position, dtype=float)
        self.local_rotations = np.array(local_rotations, dtype=float)
        self.airfoil = airfoil
        self.control_surface_type = control_surface_type
        self.control_surface_hinge_point = control_surface_hinge_point
        self.control_surface_deflection = control_surface_deflection
        self.num_spanwise_panels = num_spanwise_panels
        self.spanwise_spacing = spanwise_spacing

        # Catch invalid local rotation angles.
        if not np.all((-90.0 < self.local_rotations) & (self.local_rotations < 90.0)):
            raise Exception(
                "All local rotation angles must be in the range (-90, 90) degrees."
            )

        # Catch negative local position components.
        if np.any(self.local_position < 0.0):
            raise Exception("All local_position components must be non-negative.")

        # Define an attribute for the parent wing's unit chordwise vector, which will
        # be set by this wing cross section's parent wing's initialization method.
        self.wing_unit_chordwise_vector = None

        # Catch bad values of the chord length.
        if self.chord <= 0:
            raise Exception("A wing cross section's chord length must be >0.0 m.")

        # Catch invalid control surface hinge points.
        if not (0.0 < self.control_surface_hinge_point < 1.0):
            raise Exception(
                "A wing cross section's control surface hinge point must be on the "
                "range (0.0, 1.0)."
            )

        # Catch invalid control surface deflections.
        if not (-90.0 < self.control_surface_deflection < 90.0):
            raise Exception(
                "A wing cross section's control surface deflection must be on the "
                "range (-90.0, 90.0) degrees."
            )

        # Catch invalid values of the control surface type.
        if self.control_surface_type not in ["symmetric", "asymmetric"]:
            raise Exception(
                'A wing cross section\'s control surface type must be "symmetric" or '
                '"asymmetric".'
            )

        # Catch invalid values of the spanwise spacing.
        if self.spanwise_spacing not in ["cosine", "uniform"]:
            raise Exception(
                'A wing cross section\'s spanwise spacing must be "cosine" or '
                '"uniform".'
            )

    def validate_root_constraints(self):
        """Validate constraints specific to root WingCrossSections.

        Root WingCrossSections must have local_position and local_rotations
        set to zero arrays.

        :raises Exception: If root constraints are violated.
        """
        if not np.allclose(self.local_position, np.array([0.0, 0.0, 0.0])):
            raise Exception(
                "The root wing cross section's local position coordinates must be "
                "np.array([0.0, 0.0, 0.0])."
            )
        if not np.allclose(self.local_rotations, np.array([0.0, 0.0, 0.0])):
            raise Exception(
                "The root wing cross section's local rotation angles must be "
                "np.array([0.0, 0.0, 0.0])."
            )

    def validate_tip_constraints(self):
        """Validate constraints specific to tip WingCrossSections.

        Tip WingCrossSections must have num_spanwise_panels set to None.

        :raises Exception: If tip constraints are violated.
        """
        if self.num_spanwise_panels is not None:
            raise Exception(
                "The tip wing cross section must have num_spanwise_panels=None."
            )

    # @property
    # def unit_chordwise_vector(self):
    #     """This method defines a property for the wing cross section's unit chordwise
    #     vector.
    #
    #     The unit chordwise vector is defined as the parent wing's unit chordwise
    #     vector, rotated by the wing cross section's twist about the wing cross
    #     section's normal vector.
    #
    #     :return: (3,) array of floats
    #         This is the unit vector for the wing cross section's chordwise direction.
    #         The units are meters.
    #     """
    #     # Find the rotation matrix given the cross section's twist.
    #     twist_rotation_matrix = functions.angle_axis_rotation_matrix(
    #         self.twist * np.pi / 180, self.unit_normal_vector
    #     )
    #
    #     # Use the rotation matrix and the leading edge coordinates to calculate the
    #     # unit chordwise vector.
    #     return twist_rotation_matrix @ self.wing_unit_chordwise_vector

    # @property
    # def unit_up_vector(self):
    #     """This method defines a property for the wing cross section's unit up vector.
    #
    #     :return: (3,) array of floats
    #         This is the unit vector for the wing cross section's chordwise direction.
    #         The units are meters.
    #     """
    #     return np.cross(self.unit_chordwise_vector, self.unit_normal_vector)

    # @property
    # def trailing_edge(self):
    #     """This method defines a property for the coordinates of this wing cross
    #     section's trailing edge.
    #
    #     :return: (3,) array of floats
    #         This is an array of the coordinates of this wing cross section's trailing
    #         edge.
    #     """
    #     chordwise_vector = self.chord * self.unit_chordwise_vector
    #
    #     return self.leading_edge + chordwise_vector


class Airfoil:
    """This class is used to contain the airfoil of a WingCrossSection.

    Citation:
        Adapted from:         geometry.Airfoil in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/27/2020

    This class contains the following public methods:
        populate_coordinates: This method populates a variable with the coordinates
        of the airfoil.

        populate_mcl_coordinates: This method creates a list of the airfoil's mean
        camber line coordinates. It also creates two lists of the vectors needed to
        go from the mcl coordinates to the upper and lower surfaces. It also creates
        list of the thicknesses at the x coordinates along the mean camber line.

        leading_edge_index: This method returns the index of the point along the
        leading edge.

        lower_coordinates: This method returns a matrix of x and y coordinates that
        describe the lower surface of the airfoil.

        upper_coordinates: This method returns a matrix of x and y coordinates that
        describe the upper surface of the airfoil.

        get_downsampled_mcl: This method returns the mean camber line in a
        downsampled form.

        get_camber_at_chord_fraction: This method returns the camber of the airfoil
        at a given fraction of the chord.

        repanel_current_airfoil: This method returns a repaneled version of the
        airfoil with cosine-spaced coordinates on the upper and lower surfaces.

        add_control_surface: This method returns a version of the airfoil with a
        control surface added at a given point.

        draw: This method plots this Airfoil's coordinates using PyPlot.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        name="Untitled Airfoil",
        coordinates=None,
        repanel=True,
        n_points_per_side=400,
    ):
        """This is the initialization method.

        :param name: str, optional
            This is the name of the airfoil. It should correspond to the name in the
            airfoils directory unless you are passing in your own coordinates. The
            default is "Untitled Airfoil".
        :param coordinates: array, optional
            This is an N x 2 array of the airfoil's coordinates, where N is the
            number of coordinates. Treat this as an immutable, don't edit directly
            after initialization. If you wish to load coordinates from the airfoil
            directory, leave this as None. The default is None. Make sure that any
            airfoil coordinates used range in x from 0 to 1.
        :param repanel: bool, optional
            This is the variable that determines whether you would like to repanel
            the airfoil coordinates. This applies to coordinates passed in by the
            user or to the directory coordinates. I highly recommended setting this
            to True. The default is True.
        :param n_points_per_side: int, optional
            This is number of points to use when repaneling the airfoil. It is
            ignored if the repanel is False. The default is 400.
        """

        # Initialize the airfoil name.
        self.name = name

        # Check if the user supplied coordinates.
        if coordinates is not None:
            self.coordinates = coordinates
        else:
            # If not, populate the coordinates from the directory.
            self.populate_coordinates()

        # Initialize other user-supplied attributes.
        self.repanel = repanel
        self.n_points_per_side = n_points_per_side

        # Check that the coordinates have been set.
        assert hasattr(self, "coordinates")

        # If repanel is True, repanel the airfoil.
        if self.repanel:
            self.repanel_current_airfoil(n_points_per_side=self.n_points_per_side)

        # Initialize other attributes that will be set by populate_mcl_coordinates.
        self.mcl_coordinates = None
        self.upper_minus_mcl = None
        self.thickness = None

        # Populate the mean camber line attributes.
        self.populate_mcl_coordinates()

    def populate_coordinates(self):
        """This method populates a variable with the coordinates of the airfoil.

        The airfoil coordinates will either be generated, if the airfoil is a NACA
        4-series airfoil, or loaded from the airfoil database (a folder named
        "airfoils" in this directory, that contains a library of dat files for
        airfoil coordinates). NACA 4-series airfoil generation is an adaptation of:
        https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4
        -digit_NACA_airfoil.

        :return: None
        """

        # Sanitize the name input.
        name = self.name.lower().strip()

        # Check if the airfoil name is a NACA 4-series airfoil. If so, generate it.
        if "naca" in name:
            naca_number = name.split("naca")[1]
            if naca_number.isdigit():
                if len(naca_number) == 4:

                    # Parse the characteristics from the name.
                    max_camber = int(naca_number[0]) * 0.01
                    camber_loc = int(naca_number[1]) * 0.1
                    thickness = int(naca_number[2:]) * 0.01

                    # Set the number of points per side.
                    n_points_per_side = 100

                    # Make uncambered coordinates and generate cosine-spaced points.
                    x_t = functions.cosspace(0, 1, n_points_per_side)
                    y_t = (
                        5
                        * thickness
                        * (
                            +0.2969 * np.power(x_t, 0.5)
                            - 0.1260 * x_t
                            - 0.3516 * np.power(x_t, 2)
                            + 0.2843 * np.power(x_t, 3)
                            - 0.1015 * np.power(x_t, 4)
                        )
                    )

                    # Prevent divide by zero errors for airfoils like the NACA 0012.
                    if camber_loc == 0:
                        camber_loc = 0.5

                    # Get the camber.
                    y_c_piece1 = (
                        max_camber
                        / camber_loc**2
                        * (
                            2 * camber_loc * x_t[x_t <= camber_loc]
                            - x_t[x_t <= camber_loc] ** 2
                        )
                    )
                    y_c_piece2 = (
                        max_camber
                        / (1 - camber_loc) ** 2
                        * (
                            (1 - 2 * camber_loc)
                            + 2 * camber_loc * x_t[x_t > camber_loc]
                            - x_t[x_t > camber_loc] ** 2
                        )
                    )
                    y_c = np.hstack((y_c_piece1, y_c_piece2))

                    # Get camber slope.
                    first_piece_slope = (
                        2
                        * max_camber
                        / camber_loc**2
                        * (camber_loc - x_t[x_t <= camber_loc])
                    )
                    second_piece_slope = (
                        2
                        * max_camber
                        / (1 - camber_loc) ** 2
                        * (camber_loc - x_t[x_t > camber_loc])
                    )
                    slope = np.hstack((first_piece_slope, second_piece_slope))
                    theta = np.arctan(slope)

                    # Combine everything.
                    x_u = x_t - y_t * np.sin(theta)
                    x_l = x_t + y_t * np.sin(theta)
                    y_u = y_c + y_t * np.cos(theta)
                    y_l = y_c - y_t * np.cos(theta)

                    # Flip upper surface so it's back to front.
                    x_u, y_u = np.flipud(x_u), np.flipud(y_u)

                    # Trim 1 point from lower surface so there's no overlap.
                    x_l, y_l = x_l[1:], y_l[1:]

                    # Combine and format the coordinates.
                    x = np.hstack((x_u, x_l))
                    y = np.hstack((y_u, y_l))
                    coordinates = np.column_stack((x, y))

                    # Populate the coordinates attribute and return.
                    self.coordinates = coordinates
                    return

        # Try to read from the airfoil directory.
        try:

            # Import the airfoils package as "airfoils".
            airfoils = importlib.import_module(
                name=".airfoils",
                package="pterasoftware",
            )

            # Read the text from the airfoil file.
            raw_text = importlib.resources.read_text(airfoils, name + ".dat")

            # Trim the text at the return characters.
            trimmed_text = raw_text[raw_text.find("\n") :]

            # Input the coordinates into a 1D array.
            coordinates_1d = np.fromstring(trimmed_text, sep="\n")

            # Check to make sure the number of elements in the array is even.
            assert len(coordinates_1d) % 2 == 0, (
                "File was found in airfoil database, "
                "but it could not be read correctly."
            )

            # Reshape the 1D coordinates array into an N x 2 array, where N is the
            # number of rows.
            coordinates = np.reshape(coordinates_1d, (-1, 2))

            # Populate the coordinates attribute and return.
            self.coordinates = coordinates
            return

        # If the airfoil was not a NACA 4-series and was not found in the database,
        # throw an error.
        except FileNotFoundError:
            raise Exception("Airfoil not in database!")

    def populate_mcl_coordinates(self):
        """This method creates a list of the airfoil's mean camber line coordinates.
        It also creates two lists of the vectors needed to go from the mcl
        coordinates to the upper and lower surfaces. It also creates list of the
        thicknesses at the x coordinates along the mean camber line.

        All vectors are listed from the leading edge to the trailing edge of the
        airfoil.

        :return: None
        """

        # Get the upper and lower coordinates. Flip the upper coordinates so that it
        # is ordered from the leading edge to the trailing edge.
        upper = np.flipud(self.upper_coordinates())
        lower = self.lower_coordinates()

        # Calculate the approximate mean camber line and populate the class attribute.
        mcl_coordinates = (upper + lower) / 2
        self.mcl_coordinates = mcl_coordinates

        # Find the vectors from each mean camber line coordinate to its upper
        # coordinate.
        self.upper_minus_mcl = upper - self.mcl_coordinates

        # Create a list of values that are the thickness of the airfoil at each mean
        # camber line.
        thickness = np.sqrt(np.sum(np.power(self.upper_minus_mcl, 2), axis=1)) * 2

        # Populate the class attribute with the thicknesses at their associated x
        # coordinates.
        self.thickness = np.column_stack((self.mcl_coordinates[:, 0], thickness))

    def leading_edge_index(self):
        """Returns the index of the leading edge point.

        :return leading_edge_index: int
            This is the index of the leading edge point.
        """

        # Find the index of the coordinate pair with the minimum value of the x
        # coordinate. This is the leading edge index.
        leading_edge_index = np.argmin(self.coordinates[:, 0])

        # Return the leading edge index.
        return leading_edge_index

    def lower_coordinates(self):
        """This method returns a matrix of x and y coordinates that describe the
        lower surface of the airfoil.

        The order of the returned matrix is from leading edge to trailing edge. This
        matrix includes the leading edge point so be careful about duplicates if
        using this method in conjunction with self.upper_coordinates.

        :return lower_coordinates: array
            This is an N x 2 array of x and y coordinates that describe the lower
            surface of the airfoil, where N is the number of points.
        """

        # Find the lower coordinates.
        lower_coordinates = self.coordinates[self.leading_edge_index() :, :]

        # Return the lower coordinates.
        return lower_coordinates

    def upper_coordinates(self):
        """This method returns a matrix of x and y coordinates that describe the
        upper surface of the airfoil.

        The order of the returned matrix is from trailing edge to leading edge. This
        matrix includes the leading edge point so be careful about duplicates if
        using this method in conjunction with self.lower_coordinates.

        :return upper_coordinates: array
            This is an N x 2 array of x and y coordinates that describe the upper
            surface of the airfoil, where N is the number of points.
        """

        # Find the upper coordinates.
        upper_coordinates = self.coordinates[: self.leading_edge_index() + 1, :]

        # Return the upper coordinates.
        return upper_coordinates

    def get_downsampled_mcl(self, mcl_fractions):
        """This method returns the mean camber line in a downsampled form.

        :param mcl_fractions: 1D array
            This is a 1D array that lists the points along the mean camber line (
            normalized from 0 to 1) at which to return the mean camber line
            coordinates.
        :return mcl_downsampled: 2D array
            This is a 2D array that contains the coordinates of the downsampled mean
            camber line.
        """

        mcl = self.mcl_coordinates

        # Find the distances between points along the mean camber line, assuming
        # linear interpolation.
        mcl_distances_between_points = np.sqrt(
            np.power(mcl[:-1, 0] - mcl[1:, 0], 2)
            + np.power(mcl[:-1, 1] - mcl[1:, 1], 2)
        )

        # Create a horizontal 1D array that contains the distance along the mean
        # camber line of each point.
        mcl_distances_cumulative = np.hstack(
            (0, np.cumsum(mcl_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0 to 1.
        mcl_distances_cumulative_normalized = (
            mcl_distances_cumulative / mcl_distances_cumulative[-1]
        )

        # Linearly interpolate to find the x coordinates of the mean camber line at
        # the given mean camber line fractions.
        mcl_downsampled_x = np.interp(
            x=mcl_fractions, xp=mcl_distances_cumulative_normalized, fp=mcl[:, 0]
        )

        # Linearly interpolate to find the y coordinates of the mean camber line at
        # the given mean camber line fractions.
        mcl_downsampled_y = np.interp(
            x=mcl_fractions, xp=mcl_distances_cumulative_normalized, fp=mcl[:, 1]
        )

        # Combine the x and y coordinates of the downsampled mean camber line.
        mcl_downsampled = np.column_stack((mcl_downsampled_x, mcl_downsampled_y))

        # Return the coordinates of the downsampled mean camber line.
        return mcl_downsampled

    def get_camber_at_chord_fraction(self, chord_fraction):
        """This method returns the camber of the airfoil at a given fraction of the
        chord.

        :param chord_fraction: float
            This is a float of the fraction along the chord (normalized from 0 to 1)
            at which to return the camber.
        :return camber: float
            This is the camber of the airfoil at the requested fraction along the
            chord.
        """

        # Create a function that interpolates between the x and y coordinates of the
        # mean camber line.
        camber_function = sp_interp.interp1d(
            x=self.mcl_coordinates[:, 0],
            y=self.mcl_coordinates[:, 1],
            copy=False,
            fill_value="extrapolate",
        )

        # Find the value of the camber (the y coordinate) of the airfoil at the
        # requested chord fraction.
        camber = camber_function(chord_fraction)

        # Return the camber of the airfoil at the requested chord fraction.
        return camber

    def repanel_current_airfoil(self, n_points_per_side=100):
        """This method returns a repaneled version of the airfoil with cosine-spaced
        coordinates on the upper and lower surfaces.

        The number of points defining the final airfoil will be (n_points_per_side *
        2 - 1), since the leading edge point is shared by both the upper and lower
        surfaces.

        :param n_points_per_side: int, optional
            This is the number of points on the upper and lower surfaces. The default
            value is 100.
        :return: None
        """

        # Get the upper and lower surface coordinates. These both contain the leading
        # edge point.
        upper_original_coordinates = self.upper_coordinates()
        lower_original_coordinates = self.lower_coordinates()

        # Generate a cosine-spaced list of points from 0 to 1.
        cosine_spaced_x_values = functions.cosspace(0, 1, n_points_per_side)

        # Create interpolated functions for the x and y values of the upper and lower
        # surfaces as a function of the chord fractions
        upper_func = sp_interp.PchipInterpolator(
            x=np.flip(upper_original_coordinates[:, 0]),
            y=np.flip(upper_original_coordinates[:, 1]),
        )
        lower_func = sp_interp.PchipInterpolator(
            x=lower_original_coordinates[:, 0], y=lower_original_coordinates[:, 1]
        )

        # Find the x and y coordinates of the upper and lower surfaces at each of the
        # cosine-spaced x values.
        x_coordinates = np.hstack(
            (np.flip(cosine_spaced_x_values), cosine_spaced_x_values[1:])
        )
        y_coordinates = np.hstack(
            (
                upper_func(np.flip(cosine_spaced_x_values)),
                lower_func(cosine_spaced_x_values[1:]),
            )
        )

        # Stack the coordinates together and return them.
        coordinates = np.column_stack((x_coordinates, y_coordinates))
        self.coordinates = coordinates

    def add_control_surface(self, deflection=0.0, hinge_point=0.75):
        """This method returns a version of the airfoil with a control surface added
        at a given point.

        :param deflection: float, optional
            This is the deflection angle in degrees. Deflection downwards is
            positive. The default value is 0.0.
        :param hinge_point: float, optional
            This is the location of the hinge as a fraction of chord length. The
            default value is 0.75.
        :return flapped_airfoil: Airfoil
            This is the new airfoil with the control surface added.
        """

        # Ensure that the airfoil's deflection is not too high, which increases the
        # risk of self intersection.
        if deflection > 90 or deflection < -90:
            raise Exception("Invalid value for deflection!")

        # Make the rotation matrix for the given angle.
        sin_theta = np.sin(np.radians(-deflection))
        cos_theta = np.cos(np.radians(-deflection))
        rotation_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])

        # Find y coordinate at the hinge point x coordinate and make it a vector.
        hinge_point = np.array(
            (hinge_point, self.get_camber_at_chord_fraction(hinge_point))
        )

        # Split the airfoil into the sections before and after the hinge.
        split_index = np.where(self.mcl_coordinates[:, 0] > hinge_point[0])[0][0]
        mcl_coordinates_before = self.mcl_coordinates[:split_index, :]
        mcl_coordinates_after = self.mcl_coordinates[split_index:, :]
        upper_minus_mcl_before = self.upper_minus_mcl[:split_index, :]
        upper_minus_mcl_after = self.upper_minus_mcl[split_index:, :]

        # Rotate the mean camber line coordinates and upper minus mean camber line
        # vectors.
        new_mcl_coordinates_after = (
            np.transpose(
                rotation_matrix @ np.transpose(mcl_coordinates_after - hinge_point)
            )
            + hinge_point
        )
        new_upper_minus_mcl_after = np.transpose(
            rotation_matrix @ np.transpose(upper_minus_mcl_after)
        )

        # Assemble the new, flapped airfoil.
        new_mcl_coordinates = np.vstack(
            (mcl_coordinates_before, new_mcl_coordinates_after)
        )
        new_upper_minus_mcl = np.vstack(
            (upper_minus_mcl_before, new_upper_minus_mcl_after)
        )
        upper_coordinates = np.flipud(new_mcl_coordinates + new_upper_minus_mcl)
        lower_coordinates = new_mcl_coordinates - new_upper_minus_mcl
        coordinates = np.vstack((upper_coordinates, lower_coordinates[1:, :]))

        # Initialize the new, flapped airfoil and return it.
        flapped_airfoil = Airfoil(
            name=self.name + " flapped", coordinates=coordinates, repanel=False
        )
        return flapped_airfoil

    def draw(self):
        """This method plots this Airfoil's coordinates using PyPlot.

        :return: None
        """
        x = self.coordinates[:, 0]
        y = self.coordinates[:, 1]
        plt.plot(x, y)
        plt.xlim(0, 1)
        plt.ylim(-0.5, 0.5)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.show()
