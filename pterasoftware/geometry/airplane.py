"""This module contains the Airplane class.

This module contains the following classes:
    Airplane: This is a class used to contain airplanes.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from .airfoil import Airfoil
from .wing import Wing
from .wing_cross_section import WingCrossSection

from .. import parameter_validation


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

    1. Defining the local body axes and geometry axes
    2. Managing Wings and their coordinate transformations
    3. Processing symmetric Wings and converting them to separate wings when the
    symmetry plane is not coincident with the Wing's axes xz-plane (type 5 symmetry)
    4. Providing reference dimensions for aerodynamic calculations

    Every Airplane has a body axis system, where:
    - +x: Points forward along fuselage
    - +y: Points to the right (starboard direction)
    - +z: Points downward (completing right-handed system)

    Every Airplane has a geometry axis system, where:
    - +x: Points aft along fuselage
    - +y: Points to the right (starboard direction)
    - +z: Points upward (completing right-handed system)
    """

    def __init__(
        self,
        wings,
        name="Untitled Airplane",
        Cgi_E_I=(0.0, 0.0, 0.0),
        angles_E_to_B_izyx=(0.0, 0.0, 0.0),
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
            initialization (type 5 symmetry).

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

        :param angles_E_to_B_izyx: array-like of 3 numbers, optional

            Angles [angle1, angle2, angle3] from Earth axes to body axes using an
            intrinsic 3-2'-1" sequence. Can be a tuple, list, or numpy array of
            numbers (int or float). Values are converted to floats internally. This
            defines the orientation of the airplane's body axes relative to Earth
            axes. Note that body axes differ from geometry axes: body axes point
            forward/right/down while geometry axes point aft/right/up. The units are
            degrees. All angles must lie in the range (-180.0, 180.0] degrees. The
            default is (0.0, 0.0, 0.0).

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
        wings = parameter_validation.non_empty_list_return_list(wings, "wings")
        processed_wings = []
        for wing in wings:
            if not isinstance(wing, Wing):
                raise TypeError("Every element in wings must be a Wing")
            processed_wings.extend(self.process_wing_symmetry(wing))
        self.wings = processed_wings

        self.name = parameter_validation.string_return_string(name, "name")
        self.Cgi_E_I = parameter_validation.threeD_number_vectorLike_return_float(
            Cgi_E_I, "Cgi_E_I"
        )

        angles_E_to_B_izyx = parameter_validation.threeD_number_vectorLike_return_float(
            angles_E_to_B_izyx, "angles_E_to_B_izyx"
        )
        angles_E_to_B_izyx[0] = parameter_validation.number_in_range_return_float(
            angles_E_to_B_izyx[0], "angles_E_to_B_izyx[0]", -180.0, False, 180.0, True
        )
        angles_E_to_B_izyx[1] = parameter_validation.number_in_range_return_float(
            angles_E_to_B_izyx[1], "angles_E_to_B_izyx[1]", -180.0, False, 180.0, True
        )
        angles_E_to_B_izyx[2] = parameter_validation.number_in_range_return_float(
            angles_E_to_B_izyx[2], "angles_E_to_B_izyx[2]", -180.0, False, 180.0, True
        )
        self.angles_E_to_B_izyx = angles_E_to_B_izyx

        self.weight = parameter_validation.non_negative_number_return_float(
            weight, "weight"
        )

        # If any of the passed reference dimensions are None, set them to first Wing's
        # corresponding reference. Otherwise, set them to the passed dimension after
        # checking that it is valid.
        if s_ref is None:
            self.s_ref = self.wings[0].projected_area
        else:
            self.s_ref = parameter_validation.positive_number_return_float(
                s_ref, "s_ref"
            )
        if c_ref is None:
            self.c_ref = self.wings[0].mean_aerodynamic_chord
        else:
            self.c_ref = parameter_validation.positive_number_return_float(
                c_ref, "c_ref"
            )
        if b_ref is None:
            self.b_ref = self.wings[0].span
        else:
            self.b_ref = parameter_validation.positive_number_return_float(
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
        """This method validates constraints specific to the first Airplane in a
        simulation.

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
        # preliminary wing axes' xz-plane if symmetry_point_Wn_Ler is all zeros (no
        # translational offset), and symmetry_normal_Wn is np.array([ 0.0, 1.0,
        # 0.0]). We don't need to check types, values, or normalize because this is
        # done in Wing's init method.
        coincident_symmetry_plane = True
        if wing.symmetry_point_Wn_Ler is None or wing.symmetry_normal_Wn is None:
            coincident_symmetry_plane = False
        else:
            if not np.allclose(wing.symmetry_point_Wn_Ler, np.array([0.0, 0.0, 0.0])):
                coincident_symmetry_plane = False
            elif not np.allclose(wing.symmetry_normal_Wn, np.array([0.0, 1.0, 0.0])):
                coincident_symmetry_plane = False

        # See the Wing class docstring for the interpretation of the different
        # symmetry types.
        if not wing.symmetric:
            if not wing.mirror_only:
                # Type 1 Symmetry:
                # symmetric=False, mirror_only=False
                symmetry_type = 1
            else:
                if coincident_symmetry_plane:
                    # Type 2 Symmetry:
                    # symmetric=False, mirror_only=True, coincident_symmetry_plane=True
                    symmetry_type = 2
                else:
                    # Type 3 Symmetry:
                    # symmetric=False, mirror_only=True, coincident_symmetry_plane=False
                    symmetry_type = 3
        else:
            if coincident_symmetry_plane:
                # Type 4 Symmetry:
                # symmetric=True, coincident_symmetry_plane=True
                symmetry_type = 4
            else:
                # Type 5 Symmetry:
                # symmetric=True, coincident_symmetry_plane=False
                symmetry_type = 5

        # Based on the determined symmetry type, validate the Wing's
        # WingCrossSections' control_surface_symmetry types. From the validation done
        # during each WingCrossSection's initialization method, we already know that
        # control_surface_symmetry type is None or a valid string.
        for wing_cross_section in wing.wing_cross_sections:
            control_surface_symmetry_type = (
                wing_cross_section.control_surface_symmetry_type
            )
            if symmetry_type in [1, 2, 3]:
                if control_surface_symmetry_type is not None:
                    raise ValueError(
                        f"control_surface_symmetry_type must be None for symmetry type "
                        f"{symmetry_type}"
                    )
            else:
                if wing_cross_section.control_surface_symmetry_type is None:
                    raise ValueError(
                        f"control_surface_symmetry_type must be specified for symmetry "
                        f"type {symmetry_type}"
                    )

        # Based on symmetry type, generate the mesh and return the wing(s).
        if symmetry_type in [1, 2, 3, 4]:
            wing.generate_mesh(symmetry_type)
            return [wing]
        else:
            reflected_wing_cross_sections = []
            for wing_cross_section in wing.wing_cross_sections:
                airfoil = wing_cross_section.airfoil

                reflected_airfoil = Airfoil(
                    name=airfoil.name,
                    outline_A_lp=np.copy(airfoil.outline_A_lp),
                    resample=airfoil.resample,
                    n_points_per_side=airfoil.n_points_per_side,
                )

                if wing_cross_section.control_surface_symmetry_type == "asymmetric":
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
                        angles_Wcsp_to_Wcs_izyx=np.copy(
                            wing_cross_section.angles_Wcsp_to_Wcs_izyx
                        ),
                        control_surface_symmetry_type="symmetric",
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
                symmetry_normal_Wn=np.copy(wing.symmetry_normal_Wn),
                symmetry_point_Wn_Ler=np.copy(wing.symmetry_point_Wn_Ler),
                num_chordwise_panels=wing.num_chordwise_panels,
                chordwise_spacing=wing.chordwise_spacing,
            )

            wing.symmetric = False
            wing.mirror_only = False
            wing.symmetry_normal_Wn = None
            wing.symmetry_point_Wn_Ler = None

            wing.generate_mesh(symmetry_type=1)
            reflected_wing.generate_mesh(symmetry_type=3)
            return [wing, reflected_wing]
