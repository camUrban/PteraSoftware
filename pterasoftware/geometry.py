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
from . import parameter_validation
from . import transformations


class Airplane:
    """This is a class used to contain airplanes.

    Citation:
        Adapted from:         geometry.Airplane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/23/2020

    This class contains the following public methods:

        validate_first_airplane_constraints: This method validates that the first
        Airplane in a simulation has Cgi_E_I set to zeros, as required by the
        definition of the simulation's starting point.

        process_wing_symmetry: This method processes a Wing to determine what type of
        symmetry it has. If necessary, it then modifies the Wing. If type 5 symmetry
        is detected, it also creates a second reflected Wing. Finally, a list of
        Wings is returned. For types 1-4 symmetry this contains only the one modified
        Wing, but for type 5 symmetry it contains the modified Wing followed by the
        new reflected Wing. Before returning them, this method also calls each Wing's
        generate_mesh method, preparing them for use simulation.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.

    The Airplane class is responsible for:

    1. Defining the local geometry axes as the root coordinate system
    2. Managing Wings and their coordinate transformations
    3. Processing symmetric Wings and converting them to separate wings when
       the symmetry plane is not coincident with the Wing's xz-plane (Scenario 5)
    4. Providing reference dimensions for aerodynamic calculations
    5. Managing the moment reference point for force and moment calculations

    Every Airplane has a geometry axis system, where:
    - +x: Points aft along fuselage
    - +y: Points to the right (starboard direction)
    - +z: Points upward (completing right-handed coordinate system)
    """

    def __init__(
        self,
        wings,
        name="Untitled Airplane",
        Cgi_E_I=(0.0, 0.0, 0.0),
        angles_E_to_B_i321=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    ):
        """This is the initialization method.

        :param wings: list of Wings

            This is a list of the airplane's wings defined as Wings. It must contain
            at least one Wing. Wings with symmetric=True and non-coincident symmetry
            planes will be automatically processed into separate Wings during
            initialization (Scenario 5).

        :param name: str, optional

            A sensible name for your airplane. The default is "Untitled Airplane".

        :param Cgi_E_I: array-like of 3 numbers, optional

            Position [x, y, z] of this Airplane's starting point (in Earth axes,
            relative to the simulation's starting point). Can be a list, tuple,
            or numpy array of numbers (int or float). Values are converted to floats
            internally. For the first Airplane in a simulation, this must be [0.0,
            0.0, 0.0] since the simulation's starting point is defined as the first
            Airplane's starting point (the location of its CG at t=0). The default is
            (0.0, 0.0, 0.0).

        :param angles_E_to_B_i321: array-like of 3 numbers, optional

            Angles [angle1, angle2, angle3] from Earth axes to body axes using an
            intrinsic 3-2'-1" sequence. Can be a tuple, list, or numpy array of
            numbers (int or float). Values are converted to floats internally. This
            defines the orientation of the airplane's body axes relative to Earth axes.
            Note that body axes differ from geometry axes: body axes point forward/right/down
            while geometry axes point aft/right/up. The units are degrees. The first angle
            must lie in the range (-180.0, 180.0] degrees, the second in [-60.0, 60.0]
            degrees, and the third in (-90.0, 90.0) degrees. This is to reduce the chance
            of edge cases, and to eliminate the risk of gimbal lock. The default is
            (0.0, 0.0, 0.0).

        :param weight: number, optional

            This parameter is a number (int or float) that represents the weight of
            the aircraft in Newtons. This is used by the trim functions. It must be
            greater than or equal to zero. The default value is 0.0.

        :param s_ref: number, optional

           This parameter is a number (int or float) that represents the reference
           wetted area. If not set or set to None (the default value), it populates
           from first Wing. If set, it must be greater than zero. The units are
           square meters.

        :param c_ref: float, optional

            This parameter is a number (int or float) that represents the reference
            chord length. If not set or set to None (the default value), it populates
            from first Wing. If set, it must be greater than zero. The units are meters.

        :param b_ref: float, optional

            This parameter is a number (int or float) that represents the reference
            span. If not set or set to None (the default value), it populates from
            first Wing. If set, it must be greater than zero. The units are meters.
        """
        wings = parameter_validation.validate_non_empty_list(wings, "wings")
        processed_wings = []
        for wing in wings:
            if not isinstance(wing, Wing):
                raise TypeError("Every element in wings must be a Wing")
            processed_wings.extend(self.process_wing_symmetry(wing))
        self.wings = processed_wings

        self.name = parameter_validation.validate_string(name, "name")
        self.Cgi_E_I = parameter_validation.validate_3d_vector_float(Cgi_E_I, "Cgi_E_I")

        angles_E_to_B_i321 = parameter_validation.validate_3d_vector_float(
            angles_E_to_B_i321, "angles_E_to_B_i321"
        )
        angles_E_to_B_i321[0] = parameter_validation.validate_scalar_in_range_float(
            angles_E_to_B_i321[0], "angles_E_to_B_i321[0]", -180.0, False, 180.0, True
        )
        angles_E_to_B_i321[1] = parameter_validation.validate_scalar_in_range_float(
            angles_E_to_B_i321[1], "angles_E_to_B_i321[1]", -60.0, True, 60.0, True
        )
        angles_E_to_B_i321[2] = parameter_validation.validate_scalar_in_range_float(
            angles_E_to_B_i321[2], "angles_E_to_B_i321[2]", -90.0, False, 90.0, False
        )
        self.angles_E_to_B_i321 = angles_E_to_B_i321

        self.weight = parameter_validation.validate_non_negative_scalar_float(
            weight, "weight"
        )

        # If any of the passed reference dimensions are None, set them to first Wing's
        # corresponding reference. Otherwise, set them to the passed dimension after
        # checking that it is valid.
        if s_ref is None:
            self.s_ref = self.wings[0].projected_area
        else:
            self.s_ref = parameter_validation.validate_non_negative_scalar_float(
                s_ref, "s_ref"
            )
        if c_ref is None:
            self.c_ref = self.wings[0].mean_aerodynamic_chord
        else:
            self.c_ref = parameter_validation.validate_non_negative_scalar_float(
                c_ref, "c_ref"
            )
        if b_ref is None:
            self.b_ref = self.wings[0].span
        else:
            self.b_ref = parameter_validation.validate_non_negative_scalar_float(
                b_ref, "b_ref"
            )

        # Add up the total number of panels of all the Wings.
        self.num_panels = 0
        for wing in self.wings:
            self.num_panels += wing.num_panels

        # Initialize empty class attributes to hold the force, moment,
        # force coefficients, and moment coefficients this Airplane experiences.
        self.total_near_field_force_W = None
        self.total_near_field_force_coefficients_W = None
        self.total_near_field_moment_W = None
        self.total_near_field_moment_coefficients_W = None

    def validate_first_airplane_constraints(self):
        """This method validates constraints specific to the first Airplane in a simulation.

        The first Airplane in a simulation must have Cgi_E_I set to zeros since the
        simulation starting point is defined as the first Airplane's starting point.

        This method should be called by SteadyProblem or UnsteadyProblem classes.

        :raises Exception: If first Airplane constraints are violated.
        """
        if not np.allclose(self.Cgi_E_I, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The first Airplane in a simulation must have Cgi_E_I set to"
                "np.array([0.0, 0.0, 0.0]) since the simulation starting point "
                "is defined as the first Airplane's CG at t=0."
            )

    @staticmethod
    def process_wing_symmetry(wing):
        """This method processes a Wing to determine what type of symmetry it has. If
        necessary, it then modifies the Wing. If type 5 symmetry is detected,
        it also creates a second reflected Wing. Finally, a list of Wings is
        returned. For types 1-4 symmetry this contains only the one modified Wing,
        but for type 5 symmetry it contains the modified Wing followed by the new
        reflected Wing. Before returning them, this method also calls each Wing's
        generate_mesh method, preparing them for use simulation.

        :return: list of Wings
        """
        # Determine if the symmetry plane is coincident with the preliminary wing
        # axes' xz-plane. This is relatively easy their values are either None ( if
        # there isn't any symmetry) or relative to the preliminary wing axes.
        # Therefore, if it exists, the symmetry plane is coincident to the
        # preliminary wing axes' xz-plane if symmetry_point_prelimWn_prelimLer is all
        # zeros (no translational offset), and symmetry_normal_prelimWn is np.array([
        # 0.0, 1.0, 0.0]). We don't need to check types, values, or normalize because
        # this is done in Wing's init method.
        coincident_symmetry_plane = True
        if (
            wing.symmetry_point_prelimWn_prelimLer is None
            or wing.symmetry_normal_prelimWn is None
        ):
            coincident_symmetry_plane = False
        else:
            if not np.allclose(
                wing.symmetry_point_prelimWn_prelimLer, np.array([0.0, 0.0, 0.0])
            ):
                coincident_symmetry_plane = False
            elif not np.allclose(
                wing.symmetry_normal_prelimWn, np.array([0.0, 1.0, 0.0])
            ):
                coincident_symmetry_plane = False

        # See the Wing class docstring for the interpretation of the different
        # symmetry types.
        if not wing.symmetric:
            if not wing.mirror_only:
                # Type 1 Symmetry:
                # symmetric=False, mirror_only=False
                wing.generate_mesh(symmetry_type=1)
                return [wing]
            else:
                if coincident_symmetry_plane:
                    # Type 2 Symmetry:
                    # symmetric=False, mirror_only=True, coincident_symmetry_plane=True
                    wing.generate_mesh(symmetry_type=2)
                    return [wing]
                else:
                    # Type 3 Symmetry:
                    # symmetric=False, mirror_only=True, coincident_symmetry_plane=False
                    wing.generate_mesh(symmetry_type=3)
                    return [wing]
        else:
            if coincident_symmetry_plane:
                # Type 4 Symmetry:
                # symmetric=True, coincident_symmetry_plane=True
                wing.generate_mesh(symmetry_type=4)
                return [wing]
            else:
                # Type 5 Symmetry:
                # symmetric=True, coincident_symmetry_plane=False
                reflected_wing_cross_sections = []
                for wing_cross_section in wing.wing_cross_sections:
                    airfoil = wing_cross_section.airfoil

                    reflected_airfoil = Airfoil(
                        name=airfoil.name,
                        coordinates=np.copy(airfoil.coordinates),
                        repanel=airfoil.repanel,
                        n_points_per_side=airfoil.n_points_per_side,
                    )

                    if wing_cross_section.control_surface_type == "asymmetric":
                        reflected_control_surface_deflection = (
                            -1 * wing_cross_section.control_surface_deflection
                        )
                    else:
                        reflected_control_surface_deflection = (
                            wing_cross_section.control_surface_deflection
                        )

                    reflected_wing_cross_sections.append(
                        WingCrossSection(
                            airfoil=reflected_airfoil,
                            num_spanwise_panels=wing_cross_section.num_spanwise_panels,
                            chord=wing_cross_section.chord,
                            Lp_Wcsp_Lpp=np.copy(wing_cross_section.Lp_Wcsp_Lpp),
                            angles_Wcsp_to_Wcs_i321=np.copy(
                                wing_cross_section.angles_Wcsp_to_Wcs_i321
                            ),
                            control_surface_type=wing_cross_section.control_surface_type,
                            control_surface_hinge_point=wing_cross_section.control_surface_hinge_point,
                            control_surface_deflection=reflected_control_surface_deflection,
                            spanwise_spacing=wing_cross_section.spanwise_spacing,
                        )
                    )

                reflected_wing = Wing(
                    wing_cross_sections=reflected_wing_cross_sections,
                    name=f"Reflected {wing.name}",
                    prelimLer_G_Cg=np.copy(wing.prelimLer_G_Cg),
                    angles_G_to_prelimWn=np.copy(wing.angles_G_to_prelimWn),
                    symmetric=False,
                    mirror_only=True,
                    symmetry_normal_prelimWn=np.copy(wing.symmetry_normal_prelimWn),
                    symmetry_point_prelimWn_prelimLer=np.copy(
                        wing.symmetry_point_prelimWn_prelimLer
                    ),
                    num_chordwise_panels=wing.num_chordwise_panels,
                    chordwise_spacing=wing.chordwise_spacing,
                )

                wing.symmetric = False
                wing.mirror_only = False
                wing.symmetry_normal_prelimWn = None
                wing.symmetry_point_prelimWn_prelimLer = None

                wing.generate_mesh(symmetry_type=1)
                reflected_wing.generate_mesh(symmetry_type=3)
                return [wing, reflected_wing]


class Wing:
    """This is a class used to contain the wings of an Airplane.

    Citation:
        Adapted from:         geometry.Wing in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/24/2020

    This class contains the following public methods:

        generate_mesh: This method generates this Wing's mesh, which finishes the
        process of preparing the Wing to be used in a simulation. It is called by the
        Wing's parent Airplane, after it's determined its symmetry type.

        T_pas_G_Cg_to_Wn_Ler: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from geometry
        axes relative to the CG point to wing axes relative to the leading edge root
        point. This is set to None if the Wing's symmetry type hasn't been defined yet.

        T_pas_Wn_Ler_to_G_Cg: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from wing axes
        relative to the leading edge root point to geometry axes relative to the CG
        point. This is set to None if the Wing's symmetry type hasn't been defined yet.

        WnX_G: This method sets a property for the wing axes' first basis vector (in
        geometry axes).

        WnY_G: This method sets a property for the wing axes' second basis vector (in
        geometry axes).

        WnZ_G: This method sets a property for the wing axes' third basis vector (in
        geometry axes).

        projected_area: This method sets a property for the area of the Wing
        projected onto the plane defined by the wing axes' xy-plane.

        wetted_area: This method sets a property for the Wing's wetted area.

        span: This method sets a property for the Wing's span.

        standard_mean_chord: This method sets a property for the standard mean chord
        of the Wing.

        mean_aerodynamic_chord: This method sets a property for the mean aerodynamic
        chord of the Wing.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.

    Every Wing has its axis system, known as wing axes. The user sets the
    relationship between these axes and geometry axes with the prelimLer_G_Cg and
    angles_G_to_prelimWn parameters. However, the steps for transforming a vector
    from geometry axes to wing axes, and the interpretation of the wing axes
    orientation and position relative to a Wing's geometry, also depend on the
    parameters symmetric, mirror_only, symmetry_normal_prelimWn,
    and symmetry_point_prelimWn_prelimLer:

        1. symmetric is False

            A. mirror_only is False

                I. The symmetry plane must be undefined (symmetry_normal_prelimWn and
                symmetry_point_prelimWn_prelimLer must be None)

                    Type 1 Symmetry:

                    - prelimLer_G_Cg is the final location of the leading edge of
                    this Wing's root WingCrossSection, as defined in geometry axes.

                    - prelimLer_G_Cg is also the final location of the origin of this
                    Wing's axes, as defined in geometry axes.

                    - Translation by prelimLer_G_Cg followed by rotations by
                    angles_G_to_prelimWn fully define this Wing's axes with respect
                    to the geometry axes. The wing axes will also retain the
                    handedness of the geometry axes.

            B. mirror_only is True

                I. The symmetry plane is coincident with this Wing's axes' xz-plane

                    Type 2 Symmetry:

                    - prelimLer_G_Cg is the final location of the leading edge of
                    this Wing's root WingCrossSection, as defined in geometry axes.

                    - prelimLer_G_Cg is also the final location of the origin of this
                    Wing's axes, as defined in geometry axes.

                    - Translation by prelimLer_G_Cg followed by rotations by
                    angles_G_to_prelimWn does not fully define the orientation of
                    this Wing's axes with respect to the geometry axes. After
                    translation and rotation, the coordinate system also needs to be
                    reflected across the symmetry plane, which will flip the wing
                    axes' handedness to be opposite that of geometry axes.

                II. The symmetry plane is not coincident with this Wing's axes'
                xz-plane

                    Type 3 Symmetry:

                    - prelimLer_G_Cg is not final location of the leading edge of
                    this Wing's root WingCrossSection, as defined in geometry axes.

                    - prelimLer_G_Cg is not the final location of the origin of this
                    Wing's axes, as defined in geometry axes.

                    - Translation by prelimLer_G_Cg followed by rotations by
                    angles_G_to_prelimWn does not fully define orientation of this
                    Wing's axes with respect to the geometry axes. After translation
                    and rotation, the coordinate system also needs to be reflected
                    across the symmetry plane, which will flip the wing axes'
                    handedness to be opposite that of geometry axes.

        2. symmetric is True

            A. mirror_only must be False

                I. the symmetry plane is coincident with this Wing's axes' xz-plane

                    Type 4 Symmetry:

                    - prelimLer_G_Cg is the final location of the leading edge of
                    this Wing's root WingCrossSection, as defined in geometry axes.
                    However, while the root WingCrossSection is the still the first
                    item in the wing_cross_sections list, when meshed, panels will
                    extend from the root in both the +y and -y wing axis directions.
                    The length of the wing_cross_sections list remains unchanged.

                    - prelimLer_G_Cg is also the final location of the origin of this
                    Wing's axes, as defined in geometry axes.

                    - Translation by prelimLer_G_Cg followed by rotations by
                    angles_G_to_prelimWn fully define this Wing's axes with respect
                    to the geometry axes. The wing axes will also retain the
                    handedness of the geometry axes.

                II. the symmetry plane is not coincident with this Wing's axes' xz-plane

                    Type 5 Symmetry:

                    - This Wing's Airplane will set this Wing's symmetric parameter
                    to False, its mirror_only parameter to False,
                    its symmetry_normal_prelimWn parameter to None and its
                    symmetry_point_prelimWn_prelimLer parameter to None. These
                    changes turn this Wing into a "Scenario 1 Wing."

                    - The Airplane will also create a new Wing, and add it to its
                    wings list immediately after this Wing. The new Wing will have
                    the same name as this Wing, but with the prefix "Reflected ". The
                    new Wing also will have all the same parameters as this Wing,
                    except that symmetric will be False and mirror_only will be True,
                    which means that it will be a "Scenario 3 Wing."

                    - Also, if the control_surface_type is "asymmetric" for any of
                    this Wing's WingCrossSections, the reflected Wing's corresponding
                    WingCrossSections will have their control_surface_deflection
                    values multiplied by negative one.
    """

    def __init__(
        self,
        wing_cross_sections,
        name="Untitled Wing",
        prelimLer_G_Cg=(0.0, 0.0, 0.0),
        angles_G_to_prelimWn=(0.0, 0.0, 0.0),
        symmetric=False,
        mirror_only=False,
        symmetry_normal_prelimWn=None,
        symmetry_point_prelimWn_prelimLer=None,
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

        :param prelimLer_G_Cg: array-like of 3 numbers, optional

            This is the position [x, y, z] of the origin of this Wing's axes (in
            geometry axes, relative to the starting point) before any symmetry or
            mirror has been applied. Can be a tuple, list, or numpy array of numbers
            (int or float). Values are converted to floats internally. It may differ
            from the actual position as explained in the class docstring. The units are
            meters. The default is (0.0, 0.0, 0.0).

        :param angles_G_to_prelimWn: array-like of 3 numbers, optional

            This is the rotation angles [roll, pitch, yaw] in degrees that define the
            orientation of this Wing's axes relative to the geometry axes before any
            symmetry or mirror has been applied. Can be a tuple, list, or numpy array
            of numbers (int or float). Values are converted to floats internally. All
            angles must be in the range (-90, 90) degrees. Roll is rotation about the
            x-axis, pitch is rotation about the y-axis, and yaw is rotation about the
            z-axis. Rotations are intrinsic, and proceed in the z-y'-x'' order
            conventional for Euler angles. It may differ from the actual position as
            explained in the class docstring. The units are degrees. The default is (
            0.0, 0.0, 0.0).

        :param symmetric: bool, optional

            Set this to True if the Wing's geometry should be mirrored across the
            symmetry plane while retaining the non-mirrored side. If mirror_only is
            True, symmetric must be False. If symmetric is true, then neither
            symmetry_normal_prelimWn nor symmetry_point_prelimWn_prelimLer can be
            None. If the symmetry plane is coincident with this Wing's wing axes'
            xz-plane, the mirrored and non-mirrored geometry will be meshed as a
            single wing. If not, this Wing's Airplane will automatically create
            another Wing with the mirrored geometry, modify both Wings' parameters,
            and add the reflected Wing to its list of wings immediately following
            this one. For more details on how that process, and how this parameter
            interacts with symmetry_normal_prelimWn,
            symmetry_point_prelimWn_prelimLer, and mirror_only, see the class
            docstring. The default is False.

        :param mirror_only: bool, optional

            Set this to True if the Wing's geometry should be reflected about the
            symmetry plane without retaining the non-reflected geometry. If symmetric
            is True, mirror_only must be False. If mirror_only is true, then neither
            symmetry_normal_prelimWn nor symmetry_point_prelimWn_prelimLer can be
            None. For more details on how this parameter interacts with
            symmetry_normal_prelimWn, symmetry_point_prelimWn_prelimLer,
            and symmetric, see the class docstring. The default is False.

        :param symmetry_normal_prelimWn: array-like of 3 numbers or None, optional

            The unit normal vector (in preliminary wing axes) that, together with
            symmetry_point_prelimWn_prelimLer, defines the plane used for symmetry or
            mirroring. Can be a tuple, list, or numpy array of numbers (int or
            float), or None. Values are converted to floats and normalized
            internally. Note that reversing the normal direction (using the
            antiparallel vector) defines the same plane and produces the same result.
            This value must be None if both symmetric and mirror_only are False,
            and cannot be None if either are True. For more details on how this
            parameter interacts with symmetry_point_prelimWn_prelimLer, symmetric,
            and mirror_only, see the class docstring. The default is None.

        :param symmetry_point_prelimWn_prelimLer: array-like of 3 numbers or None,
        optional

            A point [x, y, z] (in preliminary wing axes, relative to the preliminary
            leading edge root point) that, along with symmetry_normal_prelimWn,
            defines the location of the plane about which symmetry or mirroring is
            applied. Can be a list, tuple, or numpy array of numbers (int or float),
            or None. Values are converted to floats internally. This value must be
            None if both symmetric and mirror_only are False, and cannot be None if
            either are True. For more details on how this parameter interacts with
            symmetry_normal_prelimWn, symmetric, and mirror_only, see the class
            docstring. The units are meters. The default is None.

        :param num_chordwise_panels: int, optional

            This is the number of chordwise panels to be used on this wing,
            which must be set to a positive integer. The default is 8.

        :param chordwise_spacing: str, optional

            This is the type of spacing between the wing's chordwise panels. It can
            be "cosine" or "uniform". Using cosine spacing is highly recommended for
            steady simulations and uniform spacing is highly recommended for unsteady
            simulations. The default is "cosine".

        """
        # Validate wing_cross_sections.
        wing_cross_sections = parameter_validation.validate_non_empty_list(
            wing_cross_sections, "wing_cross_sections"
        )
        num_wing_cross_sections = len(wing_cross_sections)
        if num_wing_cross_sections < 2:
            raise ValueError("wing_cross_sections must contain at least two elements.")
        for wing_cross_section_id, wing_cross_section in enumerate(wing_cross_sections):
            if not isinstance(wing_cross_section, WingCrossSection):
                raise TypeError(
                    "Every element in wing_cross_sections must be a WingCrossSection."
                )
            if wing_cross_section_id == 0:
                # Validate root WingCrossSection constraints.
                wing_cross_section.validate_root_constraints()
            elif wing_cross_section_id == num_wing_cross_sections - 1:
                # Validate tip WingCrossSection constraints.
                wing_cross_section.validate_tip_constraints()
            # Set the validated flag for this WingCrossSection.
            wing_cross_section.validated = True
        self.wing_cross_sections = wing_cross_sections

        # Validate name and prelimLer_G_Cg.
        self.name = parameter_validation.validate_string(name, "name")
        self.prelimLer_G_Cg = parameter_validation.validate_3d_vector_float(
            prelimLer_G_Cg, "prelimLer_G_Cg"
        )

        # Validate angles_G_to_prelimWn.
        angles_G_to_prelimWn = parameter_validation.validate_3d_vector_float(
            angles_G_to_prelimWn, "angles_G_to_prelimWn"
        )
        if not np.all((-90.0 < angles_G_to_prelimWn) & (angles_G_to_prelimWn < 90.0)):
            raise ValueError(
                "All elements of angles_G_to_prelimWn must lie in the range (-90, 90) degrees."
            )
        self.angles_G_to_prelimWn = angles_G_to_prelimWn

        # Validate symmetric and mirror_only.
        symmetric = parameter_validation.validate_boolean(symmetric, "symmetric")
        mirror_only = parameter_validation.validate_boolean(mirror_only, "mirror_only")
        if symmetric and mirror_only:
            raise ValueError("symmetric and mirror_only cannot both be True.")
        self.symmetric = symmetric
        self.mirror_only = mirror_only

        # Validate symmetry_normal_prelimWn and symmetry_point_prelimWn_prelimLer.
        if self.symmetric or self.mirror_only:
            if symmetry_normal_prelimWn is None:
                raise ValueError(
                    "symmetry_normal_prelimWn cannot be None when symmetric or mirror_only is True."
                )
            symmetry_normal_prelimWn = (
                parameter_validation.validate_3d_unit_vector_norm_float(
                    symmetry_normal_prelimWn, "symmetry_normal_prelimWn"
                )
            )
            if symmetry_point_prelimWn_prelimLer is None:
                raise ValueError(
                    "symmetry_point_prelimWn_prelimLer cannot be None when symmetric or mirror_only is True."
                )
            symmetry_point_prelimWn_prelimLer = (
                parameter_validation.validate_3d_vector_float(
                    symmetry_point_prelimWn_prelimLer,
                    "symmetry_point_prelimWn_prelimLer",
                )
            )
        else:
            if symmetry_normal_prelimWn is not None:
                raise ValueError(
                    "symmetry_normal_prelimWn must be None when both symmetric and mirror_only are False."
                )
            if symmetry_point_prelimWn_prelimLer is not None:
                raise ValueError(
                    "symmetry_point_prelimWn_prelimLer must be None when both symmetric and mirror_only are False."
                )
        self.symmetry_normal_prelimWn = symmetry_normal_prelimWn
        self.symmetry_point_prelimWn_prelimLer = symmetry_point_prelimWn_prelimLer

        # Validate num_chordwise_panels and chordwise_spacing.
        self.num_chordwise_panels = parameter_validation.validate_positive_scalar_int(
            num_chordwise_panels, "num_chordwise_panels"
        )
        if chordwise_spacing not in ["cosine", "uniform"]:
            raise ValueError('chordwise_spacing must be "cosine" or "uniform".')
        self.chordwise_spacing = chordwise_spacing

        # These attributes will be initialized or populated once this Wing's parent
        # Airplane calls generate_mesh.
        self.symmetry_type = None
        self.num_spanwise_panels = None
        self.num_panels = None
        self.wake_ring_vortex_vertices = None
        self.wake_ring_vortices = None
        self.panels = None

    def generate_mesh(self, symmetry_type):
        """This method generates this Wing's mesh, which finishes the process of
        preparing the Wing to be used in a simulation. It is called by the Wing's
        parent Airplane, after it's determined its symmetry type.

        :param symmetry_type:
        :return:
        """
        # Validate and apply symmetry_type. 5 isn't a valid symmetry type, because
        # the parent Airplane should have modified a Wing that initially had type 5
        # symmetry to have type 1 symmetry, and then made a new reflected Wing with
        # type 3 symmetry.
        symmetry_type = parameter_validation.validate_scalar_int(
            symmetry_type, "symmetry_type"
        )
        valid_symmetry_types = [1, 2, 3, 4]
        if symmetry_type not in valid_symmetry_types:
            raise ValueError(f"symmetry_type must be one of {valid_symmetry_types}")
        self.symmetry_type = symmetry_type

        # Find the number of spanwise panels on the wing by adding each cross
        # section's number of spanwise panels. Exclude the last cross section's
        # number of spanwise panels as this is irrelevant. If the wing has type 4
        # symmetry multiply the summation by two.
        self.num_spanwise_panels = 0
        for wing_cross_section in self.wing_cross_sections[:-1]:
            self.num_spanwise_panels += wing_cross_section.num_spanwise_panels
        if self.symmetry_type == 4:
            self.num_spanwise_panels *= 2

        # Calculate the number of panels on this wing.
        self.num_panels = self.num_spanwise_panels * self.num_chordwise_panels

        # Initialize an empty array to hold this wing's wake ring vortices and its
        # wake ring vortex vertices.
        self.wake_ring_vortex_vertices = np.empty(
            (0, self.num_spanwise_panels + 1, 3), dtype=float
        )
        self.wake_ring_vortices = np.zeros((0, self.num_spanwise_panels), dtype=object)

        # Generate the wing's mesh, which populates the Panels attribute.
        # ToDo: Uncomment this once meshing has been fixed for the refactored
        #  geometry definitions.
        # meshing.mesh_wing(self)

    @property
    def T_pas_G_Cg_to_Wn_Ler(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from geometry axes relative to the CG point
        to wing axes relative to the leading edge root point. This is set to None if
        the Wing's symmetry type hasn't been defined yet.

        :return: (4,4) ndarray of floats or None
            4x4 transformation matrix or None in cases where the Wing's symmetry type
            hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Step 1: Create T_trans_pas_G_Cg_to_G_prelimLer, which maps in homogenous
        # coordinates from geometry axes relative to the CG point to geometry axes
        # relative to the preliminary leading edge root point. This is the
        # translation step.
        T_trans_pas_G_Cg_to_G_prelimLer = transformations.generate_T_trans(
            self.prelimLer_G_Cg, passive=True
        )

        # Step 2: Create T_rot_pas_G_to_prelimWn, which maps in homogeneous
        # coordinates from geometry axes to preliminary wing axes. This is the
        # rotation step.
        T_rot_pas_G_to_prelimWn = transformations.generate_T_rot(
            transformations.generate_R(
                self.angles_G_to_prelimWn, passive=True, intrinsic=True, order="321"
            )
        )

        # Step 3: Create T_reflect_pas_prelimWn_prelimLer_to_Wn_Ler, which maps from
        # which maps in homogeneous coordinates from preliminary wing axes relative
        # to the preliminary leading edge root point to wing axes relative to the
        # leading edge root point. This is the reflection step.
        if (
            self.symmetry_normal_prelimWn is not None
            and self.symmetry_point_prelimWn_prelimLer is not None
        ):
            T_reflect_pas_prelimWn_prelimLer_to_Wn_Ler = (
                transformations.generate_T_reflect(
                    self.symmetry_point_prelimWn_prelimLer,
                    self.symmetry_normal_prelimWn,
                    passive=True,
                )
            )
        else:
            T_reflect_pas_prelimWn_prelimLer_to_Wn_Ler = np.eye(4, dtype=float)

        return (
            T_reflect_pas_prelimWn_prelimLer_to_Wn_Ler
            @ T_rot_pas_G_to_prelimWn
            @ T_trans_pas_G_Cg_to_G_prelimLer
        )

    @property
    def T_pas_Wn_Ler_to_G_Cg(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from wing axes relative to the leading edge
        root point to geometry axe
        s relative to the CG point. This is set to None if
        the Wing's symmetry type hasn't been defined yet.

        :return: (4,4) ndarray of floats or None
            4x4 transformation matrix or None in cases where the Wing's symmetry type
            hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        return np.linalg.inv(self.T_pas_G_Cg_to_Wn_Ler)

    @property
    def WnX_G(self):
        """This method sets a property for the wing axes' first basis vector (in
        geometry axes).

        :return: (3,) ndarray of floats or None
            This is the wing axes' first basis vector (in geometry axes) or None in
            cases where the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        WnX_Wn = np.array([1.0, 0.0, 0.0])

        WnXHomog_Wn = transformations.generate_homog(WnX_Wn, has_point=False)

        WnXHomog_G = self.T_pas_Wn_Ler_to_G_Cg @ WnXHomog_Wn

        return WnXHomog_G[:3]

    @property
    def WnY_G(self):
        """This method sets a property for the wing axes' second basis vector (in
        geometry axes).

        :return: (3,) ndarray of floats or None
            This is the wing axes' second basis vector (in geometry axes) or None in
            cases where the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        WnY_Wn = np.array([0.0, 1.0, 0.0])

        WnYHomog_Wn = transformations.generate_homog(WnY_Wn, has_point=False)

        WnYHomog_G = self.T_pas_Wn_Ler_to_G_Cg @ WnYHomog_Wn

        return WnYHomog_G[:3]

    @property
    def WnZ_G(self):
        """This method sets a property for the wing axes' third basis vector (in
        geometry axes).

        :return: (3,) ndarray of floats or None
            This is the wing axes' third basis vector (in geometry axes) or None in
            cases where the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        WnZ_Wn = np.array([0.0, 0.0, 1.0])

        WnZHomog_Wn = transformations.generate_homog(WnZ_Wn, has_point=False)

        WnZHomog_G = self.T_pas_Wn_Ler_to_G_Cg @ WnZHomog_Wn

        return WnZHomog_G[:3]

    @property
    def projected_area(self):
        """This method sets a property for the area of the Wing projected onto the
        plane defined by the wing axes' xy-plane.

        If the Wing is symmetric and continuous, the area of the mirrored half is
        included.

        :return projected_area: float or None
            This attribute is the projected area of the wing. It has units of square
            meters. If the Wing hasn't been meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self.panels is None:
            return None

        projected_area = 0

        # Iterate through the chordwise and spanwise indices of the panels and add
        # their area to the total projected area.
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                projected_area += self.panels[
                    chordwise_location, spanwise_location
                ].calculate_projected_area(self.WnZ_G)

        return projected_area

    @property
    def wetted_area(self):
        """This method sets a property for the Wing's wetted area.

        If the Wing is symmetric and continuous, the area of the mirrored half is
        included.

        :return wetted_area: float or None
            This attribute is the wetted area of the wing. It has units of square
            meters. If the Wing hasn't been meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self.panels is None:
            return None

        wetted_area = 0

        # Iterate through the chordwise and spanwise indices of the panels and add
        # their area to the total wetted area.
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                wetted_area += self.panels[chordwise_location, spanwise_location].area

        return wetted_area

    @property
    def span(self):
        """This method sets a property for the Wing's span.

        The span is found by first finding vector connecting the leading edges of the
        root and tip WingCrossSections. Then, this vector is projected onto the wing
        axes' third basis vector. The span is defined as the magnitude of this
        projection.

        If the Wing is symmetric and continuous, this method includes the span of
        the mirrored half.

        :return span: float or None
            This is the Wing's span. It has units of meters. None is returned in
            cases where the Wing's symmetry type hasn't been defined yet
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Calculate the tip WingCrossSection's leading point position by chaining
        # transformations through all WingCrossSections from root to tip.

        # Start at root leading point (wing axes origin)
        current_position_Wn_Ler = np.array([0.0, 0.0, 0.0])

        # For each subsequent WingCrossSection, accumulate its displacement
        for i in range(1, len(self.wing_cross_sections)):
            wing_cross_section = self.wing_cross_sections[i]

            # Get displacement in parent axes (WingCrossSection[i-1]'s axes)
            displacement_parent = wing_cross_section.Lp_Wcsp_Lpp

            # Compute transformation from parent axes to wing axes
            T_parent_to_wing = np.eye(4, dtype=float)
            for j in range(i - 1, 0, -1):  # Go backwards from parent to root
                T_parent_to_wing = (
                    self.wing_cross_sections[j].T_pas_Wcs_Lp_to_Wcsp_Lpp
                    @ T_parent_to_wing
                )

            # Transform displacement from parent axes to wing axes
            displacementHomog_parent = transformations.generate_homog(
                displacement_parent, has_point=False
            )
            displacementHomog_wing = T_parent_to_wing @ displacementHomog_parent

            # Add to accumulated position
            current_position_Wn_Ler += displacementHomog_wing[:3]

        # Now we have tip position in wing axes relative to Ler
        tipLp_Wn_Ler = current_position_Wn_Ler

        # Project the tip position onto the wing axes' y-direction (spanwise direction)
        projected_tipLp_Wn_Ler = np.dot(
            tipLp_Wn_Ler, np.array([0.0, 1.0, 0.0])
        ) * np.array([0.0, 1.0, 0.0])

        span = np.linalg.norm(projected_tipLp_Wn_Ler)

        # If the wing is symmetric and continuous, multiply the span by two.
        if self.symmetry_type == 4:
            span *= 2

        return span

    @property
    def standard_mean_chord(self):
        """This method sets a property for the standard mean chord of the Wing.

        The standard mean chord is defined as the projected area divided by the span.
        See their respective methods for the definitions of span and projected area.

        :return: float or None
            This is the standard mean chord of the Wing. It has units of meters. None
            is returned in cases where the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        return self.projected_area / self.span

    @property
    def mean_aerodynamic_chord(self):
        """This method sets a property for the mean aerodynamic chord of the Wing.

        :return: float or None
            This is the mean aerodynamic chord of the Wing. It has units of meters.
            None is returned in cases where the Wing's symmetry type hasn't been
            defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # This method is based on the equation for the mean aerodynamic chord of a
        # wing, which can be found here: https://en.wikipedia.org/wiki/Chord_(
        # aeronautics)#Mean_aerodynamic_chord. This equation integrates the squared
        # chord from the Wing's center to the Wing's tip. We will perform this
        # integral piecewise for each section of the Wing.
        integral = 0

        # Iterate through the WingCrossSections to add the contribution of their
        # corresponding Wing section to the piecewise integral.
        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections[:-1]
        ):
            next_wing_cross_section = self.wing_cross_sections[
                wing_cross_section_id + 1
            ]

            chord = wing_cross_section.chord
            next_chord = next_wing_cross_section.chord

            # Find this section's span by calculating the positions of both
            # WingCrossSections in wing axes, then finding the distance between them.

            # Calculate current WingCrossSection's position in wing axes
            current_position_Wn_Ler = np.array([0.0, 0.0, 0.0])
            for k in range(1, wing_cross_section_id + 1):
                wcs = self.wing_cross_sections[k]
                displacement_parent = wcs.Lp_Wcsp_Lpp

                # Transform from parent axes to wing axes
                T_parent_to_wing = np.eye(4, dtype=float)
                for j in range(k - 1, 0, -1):
                    T_parent_to_wing = (
                        self.wing_cross_sections[j].T_pas_Wcs_Lp_to_Wcsp_Lpp
                        @ T_parent_to_wing
                    )

                displacementHomog_parent = transformations.generate_homog(
                    displacement_parent, has_point=False
                )
                displacementHomog_wing = T_parent_to_wing @ displacementHomog_parent
                current_position_Wn_Ler += displacementHomog_wing[:3]

            # Calculate next WingCrossSection's position in wing axes
            next_position_Wn_Ler = np.array([0.0, 0.0, 0.0])
            for k in range(1, wing_cross_section_id + 2):
                wcs = self.wing_cross_sections[k]
                displacement_parent = wcs.Lp_Wcsp_Lpp

                # Transform from parent axes to wing axes
                T_parent_to_wing = np.eye(4, dtype=float)
                for j in range(k - 1, 0, -1):
                    T_parent_to_wing = (
                        self.wing_cross_sections[j].T_pas_Wcs_Lp_to_Wcsp_Lpp
                        @ T_parent_to_wing
                    )

                displacementHomog_parent = transformations.generate_homog(
                    displacement_parent, has_point=False
                )
                displacementHomog_wing = T_parent_to_wing @ displacementHomog_parent
                next_position_Wn_Ler += displacementHomog_wing[:3]

            # Find the section vector and project it onto spanwise direction
            section_vector_Wn = next_position_Wn_Ler - current_position_Wn_Ler

            # Project section vector onto spanwise direction (wing axes y-direction)
            projected_section_vector = np.dot(
                section_vector_Wn, np.array([0.0, 1.0, 0.0])
            ) * np.array([0.0, 1.0, 0.0])

            section_span = np.linalg.norm(projected_section_vector)

            # Each Wing section is, by definition, trapezoidal (at least when
            # projected on to the wing's projection plane). For a trapezoid,
            # the integral from the cited equation can be shown to evaluate to the
            # following.
            integral += (
                section_span * (chord**2 + chord * next_chord + next_chord**2) / 3
            )

        # Multiply the integral's value by the coefficients from the cited equation.
        # Double if the wing is symmetric and continuous.
        if self.symmetry_type == 4:
            return 2 * integral / self.projected_area
        return integral / self.projected_area


class WingCrossSection:
    """This class is used to contain the cross sections of a Wing.

    Citation:
        Adapted from:         geometry.WingXSec in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/26/2020

    This class contains the following public methods:

        validate_root_constraints: This method is called by the parent Wing to
        validate constraints specific to root WingCrossSections.

        validate_tip_constraints: This method is called by the parent Wing to
        validate constraints specific to tip WingCrossSections.

        T_pas_Wcsp_Lpp_to_Wcs_Lp: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from parent wing
        cross section axes relative to the parent leading point to wing cross section
        axes relative to the leading point. This is set to None if the
        WingCrossSection hasn't been fully validated yet.

        T_pas_Wcs_Lp_to_Wcsp_Lpp: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from wing cross
        section axes relative to the leading point to parent wing cross section axes
        relative to the parent leading point. This is set to None if the
        WingCrossSection hasn't been fully validated yet.

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
    for any simulations, they are merely one of the Wing attributes that tell the
    meshing function how we'd like to generate the Wing's Panels.
    """

    def __init__(
        self,
        airfoil,
        num_spanwise_panels,
        chord=1.0,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_i321=(0.0, 0.0, 0.0),
        control_surface_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing=None,
    ):
        """This is the initialization method.

        :param airfoil: Airfoil

            This is the Airfoil to be used at this WingCrossSection.

        :param num_spanwise_panels: int or None

            This is the number of spanwise panels to be used between this
            WingCrossSection and the next one. For tip WingCrossSections,
            this must be None. For all other WingCrossSections, this must be a
            positive integer.

        :param chord: number, optional

            This is the chord of the wing at this WingCrossSection. The units are
            meters. It must be greater than 0.0 and a number (int or float). The
            default value is 1.0.

        :param Lp_Wcsp_Lpp: array-like of 3 numbers, optional

            This is the position [x, y, z] in meters of this WingCrossSection's
            leading edge in parent wing cross section axes, relative to the parent
            leading edge point. Can be a tuple, list, or numpy array of numbers (int
            or float). Values are converted to floats internally. If this is the root
            WingCrossSection, the parent wing cross section axes are the wing axes
            and the parent leading point is the Wing's leading edge root point. If
            not, the parent axes and point are those of the previous
            WingCrossSection. If this is the root WingCrossSection, it must be a zero
            vector. The second component must be non-negative. The default is (0.0,
            0.0, 0.0).

        :param angles_Wcsp_to_Wcs_i321: array-like of 3 numbers, optional

            This is the angle vector of rotation angles [roll, pitch, yaw] in degrees
            that define the orientation of this WingCrossSection's axes relative to
            the parent wing cross section axes. Can be a tuple, list, or numpy array
            of numbers (int or float). Values are converted to floats internally. If
            this is a root WingCrossSection, these are the wing axes. If not,
            the parent axes are the previous WingCrossSection's axes. For the root
            WingCrossSection, this must be a zero vector. For other
            WingCrossSections, all angles must be in the range (-90, 90) degrees.
            Roll is rotation about the x-axis, pitch is rotation about the y-axis,
            and yaw is rotation about the z-axis. Rotations are intrinsic,
            and proceed in the z-y'-x'' order conventional for Euler angles. The
            units are degrees. The default is (0.0, 0.0, 0.0).

        :param control_surface_type: str, optional

            This is type of control surfaces for this WingCrossSection. It can
            be "symmetric" or "asymmetric". An example of symmetric control surfaces
            are flaps. An example of asymmetric control surfaces are ailerons. The
            default value is "symmetric". This value only affects a Wing's geometry
            if it has type 4 or 5 symmetry.

        :param control_surface_hinge_point: number, optional

            This is the location of the control surface hinge from the leading edge
            as a fraction of chord. It must be a number (int or float) the range (
            0.0, 1.0). The default value is 0.75.

        :param control_surface_deflection: number, optional

            This is the control deflection in degrees. Deflection downwards is
            positive. It must be a number (int or float) in the range (-90.0,
            90.0) degrees. The default value is 0.0 degrees.

        :param spanwise_spacing: str or None, optional

            For non-tip WingCrossSections, this can be "cosine" or "uniform". Using
            cosine spacing is highly recommended. For tip WingCrossSections it must
            be None.
        """
        # Validate airfoil.
        if not isinstance(airfoil, Airfoil):
            raise TypeError("airfoil must be an Airfoil.")
        self.airfoil = airfoil

        # Perform a preliminary validation for num_spanwise_panels. The parent Wing
        # will later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if num_spanwise_panels is not None:
            num_spanwise_panels = parameter_validation.validate_positive_scalar_int(
                num_spanwise_panels, "Non-None num_spanwise"
            )
        self.num_spanwise_panels = num_spanwise_panels

        # Validate chord.
        self.chord = parameter_validation.validate_positive_scalar_float(chord, "chord")

        # Perform a preliminary validation for Lp_Wcsp_Lpp. The parent Wing will
        # later check that this is a zero vector if this WingCrossSection is a root
        # WingCrossSection.
        Lp_Wcsp_Lpp = parameter_validation.validate_3d_vector_float(
            Lp_Wcsp_Lpp, "Lp_Wcsp_Lpp"
        )
        Lp_Wcsp_Lpp[1] = parameter_validation.validate_non_negative_scalar_float(
            Lp_Wcsp_Lpp[1], "Lp_Wcsp_Lpp[1]"
        )
        self.Lp_Wcsp_Lpp = Lp_Wcsp_Lpp

        # Perform a preliminary validation for angles_Wcsp_to_Wcs_i321. The parent
        # Wing will later check that this is a zero vector if this WingCrossSection
        # is a root WingCrossSection.
        angles_Wcsp_to_Wcs_i321 = parameter_validation.validate_3d_vector_float(
            angles_Wcsp_to_Wcs_i321, "angles_Wcsp_to_Wcs_i321"
        )
        for angle_id, angle in enumerate(angles_Wcsp_to_Wcs_i321):
            angles_Wcsp_to_Wcs_i321[angle_id] = (
                parameter_validation.validate_scalar_in_range_float(
                    angle,
                    f"angles_Wcsp_to_Wcs_i321[{angle_id}]",
                    -90.0,
                    False,
                    90.0,
                    False,
                )
            )
        self.angles_Wcsp_to_Wcs_i321 = angles_Wcsp_to_Wcs_i321

        # Validate control surface type.
        control_surface_type = parameter_validation.validate_string(
            control_surface_type, "control_surface_type"
        )
        valid_control_surface_types = ["symmetric", "asymmetric"]
        if control_surface_type not in valid_control_surface_types:
            raise ValueError(
                f"control_surface_type must be one of {valid_control_surface_types}."
            )
        self.control_surface_type = control_surface_type

        # Validate control_surface_hinge_point and control_surface_deflection.
        self.control_surface_hinge_point = (
            parameter_validation.validate_scalar_in_range_float(
                control_surface_hinge_point,
                "control_surface_hinge_point",
                0.0,
                False,
                1.0,
                False,
            )
        )
        self.control_surface_deflection = (
            parameter_validation.validate_scalar_in_range_float(
                control_surface_deflection,
                "control_surface_deflection",
                -90.0,
                False,
                90.0,
                False,
            )
        )

        # Perform a preliminary validation for spanwise_spacing. The parent Wing will
        # later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if spanwise_spacing is not None:
            spanwise_spacing = parameter_validation.validate_string(
                spanwise_spacing, "spanwise_spacing"
            )
            valid_non_none_spanwise_spacings = ["cosine", "uniform"]
            if spanwise_spacing not in valid_non_none_spanwise_spacings:
                raise ValueError(
                    f"Values for non-None spanwise_spacing must be one of "
                    f"{valid_non_none_spanwise_spacings}."
                )
        self.spanwise_spacing = spanwise_spacing

        # Define a flag for if this WingCrossSection has been fully validated. This
        # will be set by the parent Wing after calling its additional validation
        # methods.
        self.validated = False

        # ToDo: Determine if we still need this attribute. If so, uncomment it and
        #  modify Wing to set it.
        # Define an attribute for the parent Wing's unit chordwise vector, which will
        # be set by this WingCrossSection's parent Wing's initialization method.
        # self.wing_unit_chordwise_vector = None

    def validate_root_constraints(self):
        """This method is called by the parent Wing to validate constraints specific
        to root WingCrossSections.

        Root WingCrossSections must have Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_i321
        set to zero vectors.

        :raises ValueError: If root constraints are violated.
        """
        # These checks are sufficient because the types were already validated by the
        # initialization method.
        if not np.allclose(self.Lp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The root WingCrossSection's Lp_Wcsp_Lpp must be np.array([0.0, 0.0, 0.0])."
            )
        if not np.allclose(self.angles_Wcsp_to_Wcs_i321, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The root WingCrossSection's angles_Wcsp_to_Wcs_i321 must be np.array([0.0, 0.0, 0.0])."
            )

    def validate_tip_constraints(self):
        """This method is called by the parent Wing to validate constraints specific
        to tip WingCrossSections.

        Tip WingCrossSections must have num_spanwise_panels and spanwise_spacing set
        to None.

        :raises ValueError: If tip constraints are violated.
        """
        if self.num_spanwise_panels is not None:
            raise ValueError(
                "The tip WingCrossSection must have num_spanwise_panels=None."
            )
        if self.spanwise_spacing is not None:
            raise ValueError(
                "The tip WingCrossSection must have spanwise_spacing=None."
            )

    @property
    def T_pas_Wcsp_Lpp_to_Wcs_Lp(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from parent wing cross section axes relative
        to the parent leading point to wing cross section axes relative to the
        leading point. This is set to None if the WingCrossSection hasn't been
        fully validated yet.

        :return: (4,4) ndarray of floats or None
            4x4 transformation matrix or None if self.validated=False.
        """
        if not self.validated:
            return None

        # Step 1: Create T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp, which maps in homogenous
        # coordinates from parent wing cross section axes relative to the parent
        # leading point to parent wing cross section axes relative to the leading
        # point. This is the translation step.
        T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp = transformations.generate_T_trans(
            self.Lp_Wcsp_Lpp, passive=True
        )

        # Step 2: Create T_rot_pas_Wcsp_to_Wcs, which maps in homogeneous coordinates
        # from parent wing cross section axes to wing cross section axes This is the
        # rotation step.
        T_rot_pas_Wcsp_to_Wcs = transformations.generate_T_rot(
            transformations.generate_R(
                self.angles_Wcsp_to_Wcs_i321, passive=True, intrinsic=True, order="321"
            )
        )

        return T_rot_pas_Wcsp_to_Wcs @ T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp

    @property
    def T_pas_Wcs_Lp_to_Wcsp_Lpp(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from wing cross section axes relative to the
        leading point to parent wing cross section axes relative to the parent
        leading point. This is set to None if the WingCrossSection hasn't been fully
        validated yet.

        :return: (4,4) ndarray of floats or None
            4x4 transformation matrix or None if self.validated=False.
        """
        if not self.validated:
            return None

        return np.linalg.inv(self.T_pas_Wcsp_Lpp_to_Wcs_Lp)

    # ToDo: I'm going to try and replace the need for these properties by calculating
    #  them directly in the meshing function. In the new set up, they require
    #  information about the parent wing cross section axes and parent leading point,
    #  which we don't have access to.
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


# ToDo: Update this class's nomenclature to comply with the new standards discussed
#  in docs\ANGLE_VECTORS_AND_TRANSFORMATIONS.md and docs\AXES_POINTS_AND_FRAMES.md.
#  It works fine as is, so I'm leaving it alone for now.
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

        :param coordinates: (N,2) array-like of numbers, optional

            This is a Nx2 array of the airfoil's coordinates, where N is the number
            of coordinates. Treat this as an immutable, don't edit directly after
            initialization. If you wish to load coordinates from the airfoil
            directory, leave this as None, which is the default. If not, it must be
            an Nx2 numpy array of numbers (int or float). Make sure that any airfoil
            coordinates used range in x from 0 to 1.

        :param repanel: bool, optional

            This is the variable that determines whether you would like to repanel
            the airfoil coordinates. This applies to coordinates passed in by the
            user or to the directory coordinates. I highly recommended setting this
            to True. The default is True.

        :param n_points_per_side: int, optional

            This is number of points to use when repaneling the airfoil. It must be a
            positive int. It is ignored if repanel=False. The default is 400.
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
