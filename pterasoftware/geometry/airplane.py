"""This module contains the Airplane class.

This module contains the following classes:
    Airplane: This is a class used to contain airplanes.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import pyvista as pv
import time
import webp

from .airfoil import Airfoil
from .wing import Wing
from .wing_cross_section import WingCrossSection

from .. import parameter_validation
from .. import transformations


class Airplane:
    """This is a class used to contain airplanes.

    Citation:
        Adapted from:         geometry.Airplane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/23/2020

    This class contains the following public methods:

        draw: Draw the 3D geometry of this Airplane.

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

    def draw(self, save=False, testing=False):
        """Draw the 3D geometry of this Airplane.

        This method provides a convenient way to visualize the Airplane's Panels
        without needing to create a solver object first. It shows the Panel's
        surfaces in 3D using PyVista.

        :param save: bool, optional
            Set this variable to True to save the image as a WebP. The default
            value is False.
        :param testing: bool, optional
            Set this variable to True to close the image after 1 second, which is
            useful for running test suites. The default value is False.
        :return: None
        """
        # Define visualization constants
        panel_color = "chartreuse"
        plotter_background_color = "black"
        window_size = [1024, 768]
        quality = 75

        # Initialize the plotter and set it to use parallel projection
        plotter = pv.Plotter(window_size=window_size, lighting=None)
        plotter.enable_parallel_projection()

        # Initialize empty arrays to hold the Panels' vertices and faces
        panel_vertices = np.empty((0, 3), dtype=float)
        panel_faces = np.empty(0, dtype=int)

        # Initialize a variable to keep track of how many Panels' data has been added
        # to the arrays
        panel_num = 0

        # Iterate through this Airplane's Wings
        for wing in self.wings:
            # Unravel the Wing's Panel matrix and iterate through it
            panels = np.ravel(wing.panels)
            for panel in panels:
                # Stack this Panel's vertices and faces
                panel_vertices_to_add = np.vstack(
                    (
                        panel.front_left_vertex,
                        panel.front_right_vertex,
                        panel.back_right_vertex,
                        panel.back_left_vertex,
                    )
                )
                panel_face_to_add = np.array(
                    [
                        4,
                        (panel_num * 4),
                        (panel_num * 4) + 1,
                        (panel_num * 4) + 2,
                        (panel_num * 4) + 3,
                    ]
                )

                # Stack this Panel's vertices and faces with the array of all vertices and faces
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # Update the number of previous Panels
                panel_num += 1

        # Convert the Panel vertices and faces to PolyData.
        panel_surfaces = pv.PolyData(panel_vertices, panel_faces)

        # Add the Panels to the plotter
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
            smooth_shading=False,
        )

        # Set the plotter's background color
        plotter.set_background(color=plotter_background_color)

        if not testing:
            # Show the plotter so the user can adjust the camera position and window
            plotter.show(
                title=f"Airplane: {self.name}",
                cpos=(-1, -1, 1),
                full_screen=False,
                auto_close=False,
            )
        else:
            # Show the plotter for 1 second, then proceed automatically (for testing)
            plotter.show(
                title=f"Airplane: {self.name}",
                cpos=(-1, -1, 1),
                full_screen=False,
                interactive=False,
                auto_close=False,
            )
            time.sleep(1)

        # If the user wants to save the image, take a screenshot and save as WebP
        if save:
            image = webp.Image.fromarray(
                plotter.screenshot(
                    filename=None,
                    transparent_background=True,
                    return_img=True,
                )
            )
            webp.save_image(
                img=image,
                file_path=f"{self.name}_geometry.webp",
                lossless=False,
                quality=quality,
            )

        # Close all the plotters
        pv.close_all()

    # ToDo: Document and debug this method and convert it to use the standard Ptera
    #  Software theme for PyVista.
    def get_plottable_data(self, show=False):
        """

        :param show:
        :return:
        """
        # Validate the input flag.
        show = parameter_validation.boolean_return_boolean(show, "show")

        airfoilOutlines_G_Cg = []
        airfoilMcls_G_Cg = []
        for wing_id, wing in enumerate(self.wings):
            [airfoilOutlines_Wn_ler, airfoilMcls_Wn_ler] = wing.get_plottable_data(
                show=False
            )

            these_airfoilOutlines_G_Cg = []
            these_airfoilMcls_G_Cg = []
            for airfoil_id in range(len(airfoilOutlines_Wn_ler)):
                airfoilOutline_Wn_ler = airfoilOutlines_Wn_ler[airfoil_id]
                airfoilMcl_Wn_ler = airfoilMcls_Wn_ler[airfoil_id]

                airfoilOutline_G_Cg = transformations.apply_T_to_vectors(
                    wing.T_pas_Wn_Ler_to_G_Cg, airfoilOutline_Wn_ler, has_point=True
                )
                airfoilMcl_G_Cg = transformations.apply_T_to_vectors(
                    wing.T_pas_Wn_Ler_to_G_Cg, airfoilMcl_Wn_ler, has_point=True
                )

                these_airfoilOutlines_G_Cg.append(airfoilOutline_G_Cg)
                these_airfoilMcls_G_Cg.append(airfoilMcl_G_Cg)

            airfoilOutlines_G_Cg.append(these_airfoilOutlines_G_Cg)
            airfoilMcls_G_Cg.append(these_airfoilMcls_G_Cg)

        if not show:
            return [airfoilOutlines_G_Cg, airfoilMcls_G_Cg]

        plotter = pv.Plotter()

        AxesGCg = pv.AxesAssembly(
            x_label="GX@Cg",
            y_label="GY@Cg",
            z_label="GZ@Cg",
            # labels=None,
            label_color="black",
            show_labels=True,
            # label_position=(1, 1, 1),
            label_size=15,
            x_color="red",
            y_color="green",
            z_color="blue",
            # position=(0.0, 0.0, 0.0),
            # orientation=(0.0, 0.0, 0.0),
            # origin=(0.0, 0.0, 0.0),
            scale=(0.25, 0.25, 0.25),
            user_matrix=np.eye(4, dtype=float),
            name="G",
            shaft_type="cylinder",
            shaft_radius=0.025,
            shaft_length=(0.8, 0.8, 0.8),
            tip_type="cone",
            tip_radius=0.1,
            tip_length=(0.2, 0.2, 0.2),
            symmetric_bounds=False,
        )

        plotter.add_actor(AxesGCg)

        for wing_id, wing in enumerate(self.wings):
            wing_num = wing_id + 1

            AxesWLerWcs1Lp1_G_Cg = pv.AxesAssembly(
                x_label=f"W{wing_num}X@Ler/Wcs1XLp1",
                y_label=f"W{wing_num}Y@Ler/Wcs1YLp1",
                z_label=f"W{wing_num}Z@Ler/Wcs1ZLp1",
                # labels=None,
                label_color="black",
                show_labels=True,
                # label_position=(1, 1, 1),
                label_size=15,
                x_color="red",
                y_color="green",
                z_color="blue",
                # position=(0.0, 0.0, 0.0),
                # orientation=(0.0, 0.0, 0.0),
                # origin=(0.0, 0.0, 0.0),
                scale=(0.25, 0.25, 0.25),
                user_matrix=np.linalg.inv(wing.T_pas_G_Cg_to_Wn_Ler),
                # user_matrix=wingAxes_T_act,
                name=f"W{wing_num}/Wcs1",
                shaft_type="cylinder",
                shaft_radius=0.025,
                shaft_length=(0.8, 0.8, 0.8),
                tip_type="cone",
                tip_radius=0.1,
                tip_length=(0.2, 0.2, 0.2),
                symmetric_bounds=False,
            )

            plotter.add_actor(AxesWLerWcs1Lp1_G_Cg)

            these_airfoilOutlines_G_Cg = airfoilOutlines_G_Cg[wing_id]
            these_airfoilMcls_G_Cg = airfoilMcls_G_Cg[wing_id]

            for wing_cross_section_id, wing_cross_section in enumerate(
                wing.wing_cross_sections
            ):
                airfoilOutline_G_Cg = these_airfoilOutlines_G_Cg[wing_cross_section_id]
                airfoilMcl_G_Cg = these_airfoilMcls_G_Cg[wing_cross_section_id]

                airfoilOutline_faces = np.hstack(
                    [
                        airfoilOutline_G_Cg.shape[0],
                        np.arange(airfoilOutline_G_Cg.shape[0]),
                    ]
                )
                airfoilOutline_mesh = pv.PolyData(
                    airfoilOutline_G_Cg, faces=airfoilOutline_faces
                )
                plotter.add_mesh(airfoilOutline_mesh)
                plotter.add_lines(airfoilMcl_G_Cg)

                if wing_cross_section_id != 0:
                    wing_cross_section_num = wing_cross_section_id + 1

                    AxesWcsLp_G_Cg = pv.AxesAssembly(
                        x_label=f"Wcs{wing_cross_section_num}Wn{wing_num}X@Lp{wing_cross_section_num}Wn{wing_num}",
                        y_label=f"Wcs{wing_cross_section_num}Wn{wing_num}Y@Lp{wing_cross_section_num}Wn{wing_num}",
                        z_label=f"Wcs{wing_cross_section_num}Wn{wing_num}Z@Lp{wing_cross_section_num}Wn{wing_num}",
                        # labels=None,
                        label_color="black",
                        show_labels=True,
                        # label_position=(1, 1, 1),
                        label_size=15,
                        x_color="red",
                        y_color="green",
                        z_color="blue",
                        # position=(0.0, 0.0, 0.0),
                        # orientation=(0.0, 0.0, 0.0),
                        # origin=(0.0, 0.0, 0.0),
                        scale=(0.25, 0.25, 0.25),
                        user_matrix=np.linalg.inv(
                            wing.children_T_pas_G_Cg_to_Wcs_Lp[wing_cross_section_id]
                        ),
                        name=f"Wcs{wing_cross_section_id}Wn{wing_num}",
                        shaft_type="cylinder",
                        shaft_radius=0.025,
                        shaft_length=(0.8, 0.8, 0.8),
                        tip_type="cone",
                        tip_radius=0.1,
                        tip_length=(0.2, 0.2, 0.2),
                        symmetric_bounds=False,
                    )

                    plotter.add_actor(AxesWcsLp_G_Cg)

            if wing.panels is not None:
                # Initialize empty arrays to hold the Panels' vertices and faces
                panel_vertices = np.empty((0, 3), dtype=float)
                panel_faces = np.empty(0, dtype=int)

                # Initialize a variable to keep track of how many Panels' data has been added
                # to the arrays
                panel_num = 0

                # Unravel the Wing's Panel matrix and iterate through it
                panels = np.ravel(wing.panels)
                for panel in panels:
                    # Stack this Panel's vertices and faces
                    panel_vertices_to_add = np.vstack(
                        (
                            panel.front_left_vertex,
                            panel.front_right_vertex,
                            panel.back_right_vertex,
                            panel.back_left_vertex,
                        )
                    )
                    panel_face_to_add = np.array(
                        [
                            4,
                            (panel_num * 4),
                            (panel_num * 4) + 1,
                            (panel_num * 4) + 2,
                            (panel_num * 4) + 3,
                        ]
                    )

                    # Stack this Panel's vertices and faces with the array of all vertices and faces
                    panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                    panel_faces = np.hstack((panel_faces, panel_face_to_add))

                    # Update the number of previous Panels
                    panel_num += 1

                    # Convert the Panel vertices and faces to PolyData.
                    panel_surfaces = pv.PolyData(panel_vertices, panel_faces)

                    # Add the Panels to the plotter
                    plotter.add_mesh(
                        panel_surfaces,
                        show_edges=True,
                        color="chartreuse",
                        smooth_shading=False,
                    )

        plotter.enable_parallel_projection()

        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )

        return None

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
