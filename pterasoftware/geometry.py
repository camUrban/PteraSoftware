"""This module contains useful class definitions for different types of geometries.

This module contains the following classes:
    Airplane: This is a class used to contain airplanes.

    Wing: This is a class used to contain the wings of an Airplane object.

    WingCrossSection: This class is used to contain the cross sections of a Wing
    object.

    Airfoil: This class is used to contain the airfoil of a WingCrossSection object.

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


# ToDo: Add a check that wings has at least one elements.
class Airplane:
    """This is a class used to contain airplanes.

    Citation:
        Adapted from:         geometry.Airplane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/23/2020

    This class contains the following public methods:
        set_reference_dimensions_from_wing: This method sets the reference dimensions
        of the airplane from measurements obtained from the main wing.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        wings,
        name="Untitled Airplane",
        x_ref=0.0,
        y_ref=0.0,
        z_ref=0.0,
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    ):
        """This is the initialization method.

        :param wings: list of Wing objects
            This is a list of the airplane's wings defined as Wing objects.
        :param name: str, optional
            A sensible name for your airplane. The default is "Untitled Airplane".
        :param x_ref: float, optional
            This is the x coordinate of the moment reference point. It should be the
            x coordinate of the center of gravity. The default is 0.0.
        :param y_ref: float, optional
            This is the y coordinate of the moment reference point. It should be the
            y coordinate of the center of gravity. The default is 0.0.
        :param z_ref: float, optional
            This is the z coordinate of the moment reference point. It should be the
            z coordinate of the center of gravity. The default is 0.0.
        :param weight: float, optional
            This parameter holds the weight of the aircraft in Newtons. This is used
            by the trim functions. The default value is 0.0.
        :param s_ref: float, optional if more than one wing is in the wings list.
            This is the reference wetted area. If not set, it populates from first
            wing object.
        :param c_ref: float, optional if more than one wing is in the wings list.
            This is the reference chord length. If not set, it populates from first
            wing object.
        :param b_ref: float, optional if more than one wing is in the wings list.
            This is the reference calculate_span. If not set, it populates from first
            wing object.
        """

        # Initialize the name and the moment reference point.
        self.name = name
        self.x_ref = x_ref
        self.y_ref = y_ref
        self.z_ref = z_ref
        self.xyz_ref = np.array([self.x_ref, self.y_ref, self.z_ref])

        # Initialize the weight.
        self.weight = weight

        self.wings = wings

        # Set the wing reference dimensions to be the main wing's reference dimensions.
        self.set_reference_dimensions_from_main_wing()

        # If any of the passed reference dimensions are not None, set that reference
        # dimension to be what was passed.
        if s_ref is not None:
            self.s_ref = s_ref
        if c_ref is not None:
            self.c_ref = c_ref
        if b_ref is not None:
            self.b_ref = b_ref

        # Calculate the number of panels in the entire current_airplane.
        self.num_panels = 0
        for wing_position, wing in enumerate(self.wings):
            self.num_panels += wing.num_panels

        # Initialize empty class attributes to hold the force, moment,
        # force coefficients, and moment coefficients this airplane experiences after
        self.total_near_field_force_wind_axes = None
        self.total_near_field_force_coefficients_wind_axes = None
        self.total_near_field_moment_wind_axes = None
        self.total_near_field_moment_coefficients_wind_axes = None

    def set_reference_dimensions_from_main_wing(self):
        """This method sets the reference dimensions of the airplane from
        measurements obtained from the main wing.

        This method assumes the main wing to be the first wing in the wings list
        passed by the user.

        :return: None
        """

        # Define the main wing to be the first wing in the wings list.
        main_wing = self.wings[0]

        # Set the objects reference dimension attributes to be the reference
        # dimension attributes of the main wing. These attributes are calculated via
        # methods in the Wing class.
        self.s_ref = main_wing.projected_area
        self.b_ref = main_wing.span
        self.c_ref = main_wing.mean_aerodynamic_chord


class Wing:
    """This is a class used to contain the wings of an Airplane object.

    If the wing is symmetric about some plane, just define the right half,
    supply "symmetric=True" in the constructor, and include the symmetry plane's
    normal vector in the "symmetry_plane_normal" attribute.

    Citation:
        Adapted from:         geometry.Wing in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/24/2020

    This class contains the following public methods:
        projected_area: This method calculates the projected area of the wing and
        assigns it to the projected_area attribute.

        wetted_area: This method calculates the wetted area of the wing based on the
        areas of its panels and assigns it to the wetted_area attribute.

        span: This method calculates the span of the wing and assigns it to the span
        attribute.

        standard_mean_chord: This method calculates the standard mean chord of the
        wing and assigns it to the standard_mean_chord attribute.

        mean_aerodynamic_chord: This method calculates the mean aerodynamic chord of
        the wing and assigns it to the mean_aerodynamic_chord attribute.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # ToDo: Add a check that wing_cross_sections has at least two elements.
    def __init__(
        self,
        wing_cross_sections,
        name="Untitled Wing",
        x_le=0.0,
        y_le=0.0,
        z_le=0.0,
        symmetry_unit_normal_vector=np.array([0.0, 1.0, 0.0]),
        symmetric=False,
        unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    ):
        """This is the initialization method.

        :param wing_cross_sections: list of WingCrossSection objects
            This is a list of WingCrossSection objects, that represent the wing's
            cross sections.
        :param name: str, optional
            This is a sensible name for the wing. The default is "Untitled Wing".
        :param x_le: float, optional
            This is the x coordinate of the leading edge of the wing, relative to the
            airplane's reference point. The default is 0.0.
        :param y_le: float, optional
            This is the y coordinate of the leading edge of the wing, relative to the
            airplane's reference point. The default is 0.0.
        :param z_le: float, optional
            This is the z coordinate of the leading edge of the wing, relative to the
            airplane's reference point. The default is 0.0.
        :param symmetry_unit_normal_vector: ndarray, optional
            This is an (3,) ndarray of floats that represents the unit normal vector
            of the wing's symmetry plane. It is also the direction vector that the
            wing's span will be assessed relative to. Additionally, this vector
            crossed with the "unit_chordwise_vector" defines the normal vector of the
            plane that the wing's projected area will reference. It must be
            equivalent to this wing's root wing cross section's "unit_normal_vector"
            attribute. The default is np.array([0.0, 1.0, 0.0]), which is the XZ
            plane's unit normal vector.
        :param symmetric: bool, optional
            Set this to true if the wing is across the xz plane. Set it to false if
            not. The default is false.
        :param unit_chordwise_vector: ndarray, optional
            This is an (3,) ndarray of floats that represents the unit vector that
            defines the wing's chordwise direction. This vector crossed with the
            "symmetry_unit_normal_vector" defines the normal vector of the plane that
            the wing's projected area will reference. This vector must be parallel to
            the intersection of the wing's symmetry plane with each of its wing cross
            section's planes. The default is np.array([1.0, 0.0, 0.0]), which is the
            X unit vector.
        :param num_chordwise_panels: int, optional
            This is the number of chordwise panels to be used on this wing. The
            default is 8.
        :param chordwise_spacing: str, optional
            This is the type of spacing between the wing's chordwise panels. It can
            be set to "cosine" or "uniform". Using a cosine spacing is highly
            recommended for steady simulations and a uniform spacing is highly
            recommended for unsteady simulations. The default is "cosine".
        """
        # Initialize the name and the position of the wing's leading edge.
        self.name = name
        self.x_le = x_le
        self.y_le = y_le
        self.z_le = z_le
        self.symmetry_unit_normal_vector = symmetry_unit_normal_vector
        self.leading_edge = np.array([self.x_le, self.y_le, self.z_le])

        # Check that the wing's symmetry plane is equal to its root wing cross
        # section's plane.
        if not np.array_equal(
            symmetry_unit_normal_vector, wing_cross_sections[0].unit_normal_vector
        ):
            raise Exception(
                "The wing's symmetry plane must be the same as its root wing cross "
                "section's plane."
            )
        # Check that the root wing cross section's leading edge isn't offset from the
        # wing's leading edge.
        if np.any(wing_cross_sections[0].leading_edge):
            raise Exception(
                "The root wing cross section's leading edge must not be offset from "
                "the wing's leading edge."
            )

        # Initialize the other attributes.
        self.wing_cross_sections = wing_cross_sections
        self.symmetric = symmetric
        self.unit_chordwise_vector = unit_chordwise_vector
        self.num_chordwise_panels = num_chordwise_panels
        self.chordwise_spacing = chordwise_spacing

        # Catch invalid values of chordwise_spacing.
        if self.chordwise_spacing not in ["cosine", "uniform"]:
            raise Exception("Invalid value of chordwise_spacing!")

        # Find the number of spanwise panels on the wing by adding each cross
        # section's number of spanwise panels. Exclude the last cross section's
        # number of spanwise panels as this is irrelevant. If the wing is symmetric,
        # multiply the summation by two.
        self.num_spanwise_panels = 0
        for wing_cross_section in self.wing_cross_sections[:-1]:
            self.num_spanwise_panels += wing_cross_section.num_spanwise_panels
        if self.symmetric:
            self.num_spanwise_panels *= 2

        # ToDo: Document this section.
        for wing_cross_section in self.wing_cross_sections:
            # Find the vector parallel to the intersection of this wing cross
            # section's plane and the wing's symmetry plane.
            orthogonal_vector = np.cross(
                self.symmetry_unit_normal_vector, wing_cross_section.unit_normal_vector
            )

            if np.any(np.cross(orthogonal_vector, self.unit_chordwise_vector)):
                raise Exception(
                    "Every wing cross section's plane must intersect with the wing's "
                    "symmetry plane along a line that is parallel with the wing's "
                    "chordwise direction."
                )
            if (
                np.dot(self.unit_chordwise_vector, self.symmetry_unit_normal_vector)
                != 0
            ):
                raise Exception(
                    "Every wing cross section's plane must intersect with the wing's "
                    "symmetry plane along a line that is parallel with the wing's "
                    "chordwise direction."
                )

        # Calculate the number of panels on this wing.
        self.num_panels = self.num_spanwise_panels * self.num_chordwise_panels

        for wing_cross_section in wing_cross_sections:
            wing_cross_section.unit_chordwise_vector = self.unit_chordwise_vector

        # Initialize the panels attribute. Then mesh the wing, which will
        # populate this attribute.
        self.panels = None
        meshing.mesh_wing(self)

        # Initialize an empty array to hold this wing's wake ring vortices and its
        # wake ring vortex vertices.
        self.wake_ring_vortex_vertices = np.empty((0, self.num_spanwise_panels + 1, 3))
        self.wake_ring_vortices = np.zeros((0, self.num_spanwise_panels), dtype=object)

        # Define an attribute that is the normal vector of the plane that the
        # projected area will reference.
        self.projected_unit_normal_vector = np.cross(
            unit_chordwise_vector, symmetry_unit_normal_vector
        )

    # ToDo: Update this method's documentation.
    @property
    def projected_area(self):
        """This method calculates the projected area of the wing and assigns it to
        the projected_area attribute.

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
                ].calculate_projected_area(self.projected_unit_normal_vector)

        return projected_area

    # ToDo: Update this method's documentation.
    @property
    def wetted_area(self):
        """This method calculates the wetted area of the wing based on the areas of
        its panels and assigns it to the wetted_area attribute.

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

    # ToDo: Update this method's documentation.
    @property
    def span(self):
        """This method calculates the span of the wing and assigns it to the span
        attribute. The span is found first finding vector connecting the leading
        edges of the root and tip wing cross sections. Then, this vector is projected
        onto the symmetry plane's unit normal vector. The span is defined as the
        magnitude of this projection.

        If the wing is symmetrical, this method includes the span of the mirrored half.

        :return span: float
            This attribute is the wingspan. It has units of meters.
        """
        root_to_tip_leading_edge = (
            self.wing_cross_sections[-1].leading_edge
            - self.wing_cross_sections[0].leading_edge
        )

        projected_leading_edge = (
            np.dot(root_to_tip_leading_edge, self.symmetry_unit_normal_vector)
            * self.symmetry_unit_normal_vector
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
                np.dot(section_leading_edge, self.symmetry_unit_normal_vector)
                * self.symmetry_unit_normal_vector
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
    """This class is used to contain the cross sections of a Wing object.

    Citation:
        Adapted from:         geometry.WingXSec in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/26/2020

    This class contains the following public methods:
        trailing_edge: This method calculates the coordinates of the trailing edge of
        this wing cross section and assigns them to the trailing_edge attribute.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        airfoil,
        x_le=0.0,
        y_le=0.0,
        z_le=0.0,
        chord=1.0,
        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
        twist=0.0,
        control_surface_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        num_spanwise_panels=8,
        spanwise_spacing="cosine",
    ):
        """This is the initialization method.

        :param airfoil: Airfoil
            This is the airfoil to be used at this wing cross section.
        :param x_le: float, optional
            This is the x coordinate of the leading edge of the wing cross section
            relative to the wing's datum. The default value is 0.0.
        :param y_le: float, optional
            This is the y coordinate of the leading edge of the wing cross section
            relative to the wing's leading edge. The default value is 0.0.
        :param z_le: float, optional
            This is the z coordinate of the leading edge of the wing cross section
            relative to the wing's datum. The default value is 0.0.
        :param chord: float, optional
            This is the chord of the wing at this wing cross section. The default
            value is 1.0.
        :param unit_normal_vector: ndarray, optional
            This is an (3,) ndarray of floats that represents the unit normal vector
            of the plane this wing cross section lies on. If this wing cross section
            is a wing's root, this vector must be equal to the wing's
            symmetry_unit_normal_vector attribute. Also, every wing cross section
            must have a plane that intersects its parent wing's symmetry plane at a
            line parallel to the parent wing's "unit_chordwise_vector". The default
            is np.array([ 0.0, 1.0, 0.0]), which is the XZ plane's unit normal vector.
        :param twist: float, optional
            This is the twist of the cross section about the leading edge in degrees.
            The default value is 0.0.
        :param control_surface_type: str, optional
            This is type of control surfaces for this wing cross section. It can be
            "symmetric" or "asymmetric". An example of symmetric control surfaces are
            flaps. An example of asymmetric control surfaces are ailerons. The
            default value is "symmetric".
        :param control_surface_hinge_point: float, optional
            This is the location of the control surface hinge from the leading edge
            as a fraction of chord. The default value is 0.75.
        :param control_surface_deflection: float, optional
            This is the control deflection in degrees. Deflection downwards is
            positive. The default value is 0.0 degrees.
        :param num_spanwise_panels: int, optional
            This is the number of spanwise panels to be used between this wing cross
            section and the next one. The default value is 8.
        :param spanwise_spacing: str, optional
            This can be "cosine" or "uniform". Using cosine spacing is highly
            recommended. The default value is "cosine".
        """

        # Initialize all the class attributes.
        self.x_le = x_le
        self.y_le = y_le
        self.z_le = z_le
        self.chord = chord
        self.unit_normal_vector = unit_normal_vector
        self.twist = twist
        self.airfoil = airfoil
        self.control_surface_type = control_surface_type
        self.control_surface_hinge_point = control_surface_hinge_point
        self.control_surface_deflection = control_surface_deflection
        self.num_spanwise_panels = num_spanwise_panels
        self.spanwise_spacing = spanwise_spacing
        self.leading_edge = np.array([x_le, y_le, z_le])

        # ToDo: Document this
        self.unit_chordwise_vector = None

        # Catch bad values of the chord length.
        if self.chord <= 0:
            raise Exception("Invalid value of chord")

        # Catch invalid values of control_surface_type.
        if self.control_surface_type not in ["symmetric", "asymmetric"]:
            raise Exception("Invalid value of control_surface_type")

        # Catch invalid values of spanwise_spacing.
        if self.spanwise_spacing not in ["cosine", "uniform"]:
            raise Exception("Invalid value of spanwise_spacing!")

    @property
    def trailing_edge(self):
        """This method calculates the coordinates of the trailing edge of this wing
        cross section and assigns them to the trailing_edge attribute.

        :return: array
            This is a 1D array that contains the coordinates of this wing cross
            section's trailing edge.
        """

        # Find the rotation matrix given the cross section's twist.
        twist_rotation_matrix = functions.angle_axis_rotation_matrix(
            self.twist * np.pi / 180, np.array([0, 1, 0])
        )

        chordwise_vector = self.chord * self.unit_chordwise_vector

        # Use the rotation matrix and the leading edge coordinates to calculate the
        # trailing edge coordinates.
        return self.leading_edge + twist_rotation_matrix @ chordwise_vector

    # ToDo: Document this
    @property
    def unit_up_vector(self):
        """

        :return:
        """
        return np.cross(self.unit_chordwise_vector, self.unit_normal_vector)


class Airfoil:
    """This class is used to contain the airfoil of a WingCrossSection object.

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
            self.populate_coordinates()  # populates self.coordinates

        # Check that the coordinates have been set.
        assert hasattr(self, "coordinates")

        # Initialize other attributes.
        self.repanel = repanel
        self.mcl_coordinates = None
        self.upper_minus_mcl = None
        self.thickness = None
        self.n_points_per_side = n_points_per_side

        # If repanel is True, repanel the airfoil.
        if self.repanel:
            self.repanel_current_airfoil(n_points_per_side=self.n_points_per_side)

        # Populate the mean camber line attributes.
        self.populate_mcl_coordinates()

    def populate_coordinates(self):
        """This method populates a variable with the coordinates of the airfoil.

        The airfoil coordinates will either be generated, if the airfoil is a NACA
        4-series airfoil, or loaded from the airfoil database (a folder named
        "airfoils" in this directory, that contains a library of dat files for
        airfoil coordinates). NACA 4-series airfoil generation is an adaptation of:
        https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil.

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
