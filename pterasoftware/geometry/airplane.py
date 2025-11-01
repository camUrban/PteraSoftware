"""Contains the Airplane class.

Contains the following classes:
    Airplane: A class used to contain airplanes.

Contains the following functions:
    None
"""

from collections.abc import Sequence

import numpy as np
import pyvista as pv
import time
import webp

from . import airfoil as airfoil_mod
from . import wing as wing_mod
from . import wing_cross_section as wing_cross_section_mod
from .. import _parameter_validation
from .. import _transformations


class Airplane:
    """A class used to contain airplanes.

    Citation:
        Adapted from:
            geometry.Airplane in AeroSandbox
        Author:
            Peter Sharpe
        Date of Retrieval:
            04/23/2020

    Contains the following methods:
        draw: Draws the 3D geometry of this Airplane.

        get_plottable_data: Returns plottable data for this Airplane's Airfoils'
        outlines and mean camber lines.

        num_panels: Sets a property for the total number of Panels across all Wings.

        validate_first_airplane_constraints: Validates that the first Airplane in a
        simulation has Cg_E_CgP1 set to zeros.

        compute_T_pas_G_Cg_to_GP1_CgP1: Computes the passive transformation matrix
        from this Airplane's geometry axes, relative to this Airplane's CG to the
        first Airplane's geometry axes, relative to the first Airplane's CG.

        process_wing_symmetry: Processes a Wing to determine what type of symmetry it
        has. If necessary, it then modifies the Wing. If type 5 symmetry is detected,
        it also creates a second reflected Wing. Finally, it returns a list of Wings.

    The Airplane class is responsible for:
        1. Defining the local body axes and geometry axes
        2. Managing Wings and their coordinate transformations
        3. Processing symmetric Wings and converting them to separate wings when the
           symmetry plane is not coincident with the Wing's axes xz-plane
           (type 5 symmetry)
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
        wings: list[wing_mod.Wing],
        name: str = "Untitled Airplane",
        Cg_E_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        angles_E_to_B_izyx: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        weight: float | int = 0.0,
        s_ref: float | int | None = None,
        c_ref: float | int | None = None,
        b_ref: float | int | None = None,
    ) -> None:
        """The initialization method.

        :param wings: A list of the airplane's wings defined as Wings. It must contain
            at least one Wing. Wings with symmetric=True and non-coincident symmetry
            planes will be automatically processed into separate Wings during
            initialization (type 5 symmetry).
        :param name: A sensible name for your airplane. The default is "Untitled
            Airplane".
        :param Cg_E_CgP1: An array-like object of 3 numbers representing the position
            of this Airplane's CG (in Earth axes, relative to the first Airplane's CG).
            Can be a list, tuple, or numpy array of numbers (int or float). Values are
            converted to floats internally. For the first Airplane in a simulation, this
            must be equivalent to (0.0, 0.0, 0.0) by definition. Earth axes follow the
            North-East-Down convention. The units are in meters. The default is (0.0,
            0.0, 0.0).
        :param angles_E_to_B_izyx: An array-like object of 3 numbers representing the
            angles from Earth axes to body axes using an intrinsic z-y'-x" sequence. Can
            be a tuple, list, or numpy array of numbers (int or float). Values are
            converted to floats internally. It defines the orientation of the airplane's
            body axes relative to Earth axes. Note that body axes differ from geometry
            axes: body axes point forward/right/down while geometry axes point
            aft/right/up. The units are degrees. All angles must lie in the range
            (-180.0, 180.0] degrees. The default is (0.0, 0.0, 0.0).
        :param weight: A number (int or float) representing the weight of the aircraft
            in Newtons. This is used by the trim functions. It must be greater than or
            equal to zero. The default is 0.0.
        :param s_ref: A number (int or float) representing the reference wetted area. If
            not set or set to None (the default), it populates from first Wing. If set,
            it must be greater than zero, and will be converted to a float internally.
            The units are square meters.
        :param c_ref: A number (int or float) representing the reference chord length.
            If not set or set to None (the default), it populates from first Wing. If
            set, it must be greater than zero, and will be converted to a float
            internally. The units are meters.
        :param b_ref: A number (int or float) representing the reference span. If not
            set or set to None (the default value), it populates from first Wing. If
            set, it must be greater than zero, and will be converted to a float
            internally. The units are meters.
        """
        wings = _parameter_validation.non_empty_list_return_list(wings, "wings")
        processed_wings = []
        for wing in wings:
            if not isinstance(wing, wing_mod.Wing):
                raise TypeError("Every element in wings must be a Wing")
            processed_wings.extend(self.process_wing_symmetry(wing))
        self.wings = processed_wings

        self.name = _parameter_validation.string_return_string(name, "name")
        self.Cg_E_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            Cg_E_CgP1, "Cg_E_CgP1"
        )

        angles_E_to_B_izyx = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                angles_E_to_B_izyx, "angles_E_to_B_izyx"
            )
        )
        angles_E_to_B_izyx[0] = _parameter_validation.number_in_range_return_float(
            angles_E_to_B_izyx[0], "angles_E_to_B_izyx[0]", -180.0, False, 180.0, True
        )
        angles_E_to_B_izyx[1] = _parameter_validation.number_in_range_return_float(
            angles_E_to_B_izyx[1], "angles_E_to_B_izyx[1]", -180.0, False, 180.0, True
        )
        angles_E_to_B_izyx[2] = _parameter_validation.number_in_range_return_float(
            angles_E_to_B_izyx[2], "angles_E_to_B_izyx[2]", -180.0, False, 180.0, True
        )
        self.angles_E_to_B_izyx = angles_E_to_B_izyx

        self.weight = _parameter_validation.non_negative_number_return_float(
            weight, "weight"
        )

        # If any of the passed reference dimensions are None, set them to first Wing's
        # corresponding reference. Otherwise, set them to the passed dimension after
        # checking that it is valid.
        if s_ref is None:
            self.s_ref = self.wings[0].projected_area
        else:
            self.s_ref = _parameter_validation.positive_number_return_float(
                s_ref, "s_ref"
            )
        if c_ref is None:
            self.c_ref = self.wings[0].mean_aerodynamic_chord
        else:
            self.c_ref = _parameter_validation.positive_number_return_float(
                c_ref, "c_ref"
            )
        if b_ref is None:
            self.b_ref = self.wings[0].span
        else:
            self.b_ref = _parameter_validation.positive_number_return_float(
                b_ref, "b_ref"
            )

        # Initialize empty class attributes to hold the force, moment,
        # force coefficients, and moment coefficients this Airplane experiences.
        self.forces_W = None
        self.forceCoefficients_W = None
        self.moments_W_CgP1 = None
        self.momentCoefficients_W_CgP1 = None

    def draw(
        self, save: bool | np.bool_ = False, testing: bool | np.bool_ = False
    ) -> None:
        """Draws the 3D geometry of this Airplane.

        This method provides a convenient way to visualize the Airplane's Panels without
        needing to create a solver object first. It shows the Panel's surfaces in 3D
        using PyVista.

        :param save: Set to True to save the image as a WebP. Can be a bool or a numpy
            bool and will be converted internally to bool. The default value is False.
        :param testing: Set to True to close the image after 1 second, which is useful
            for running test suites. Can be a bool or a numpy bool and will be converted
            internally to bool. The default value is False.
        :return: None
        """
        save = _parameter_validation.boolLike_return_bool(save, "save")
        testing = _parameter_validation.boolLike_return_bool(testing, "testing")

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
                        panel.Flpp_G_Cg,
                        panel.Frpp_G_Cg,
                        panel.Brpp_G_Cg,
                        panel.Blpp_G_Cg,
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

    # TEST: Consider adding unit tests for this method.
    # DOCUMENT: After testing it, document this method.
    def get_plottable_data(
        self, show: bool | np.bool_ = False
    ) -> list[list[np.ndarray]] | None:
        """Returns plottable data for this Airplane's Airfoils' outlines and mean camber
        lines.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: If show is True, returns None. If show is False, returns a list of two
            lists, each containing one ndarray for every one of this Airplane's
            Airfoils. These ndarrays represent points on each Airfoil's outline and mean
            camber lines, respectively. The points are in geometry axes, relative to the
            CG. The units are in meters.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

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

                airfoilOutline_G_Cg = _transformations.apply_T_to_vectors(
                    wing.T_pas_Wn_Ler_to_G_Cg, airfoilOutline_Wn_ler, has_point=True
                )
                airfoilMcl_G_Cg = _transformations.apply_T_to_vectors(
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

                # Initialize a variable to keep track of how many Panels' data has
                # been added to the arrays
                panel_num = 0

                # Unravel the Wing's Panel matrix and iterate through it
                panels = np.ravel(wing.panels)
                for panel in panels:
                    # Stack this Panel's vertices and faces
                    panel_vertices_to_add = np.vstack(
                        (
                            panel.Flpp_G_Cg,
                            panel.Frpp_G_Cg,
                            panel.Brpp_G_Cg,
                            panel.Blpp_G_Cg,
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

                    # Stack this Panel's vertices and faces with the array of all
                    # vertices and faces
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

    @property
    def num_panels(self) -> int:
        """Sets a property for the total number of Panels across all Wings.

        :return: The total number of Panels.
        """
        return sum(wing.num_panels for wing in self.wings)

    def validate_first_airplane_constraints(self) -> None:
        """Validates that the first Airplane in a simulation has Cg_E_CgP1 set to zeros.

        This method should be called by SteadyProblem or UnsteadyProblem.

        :return: None
        """
        if not np.allclose(self.Cg_E_CgP1, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The first Airplane in a simulation must have Cg_E_CgP1 set to"
                "(0.0, 0.0, 0.0) by definition."
            )

    # TEST: Add unit tests for this method.
    def compute_T_pas_G_Cg_to_GP1_CgP1(self, first_airplane: "Airplane") -> np.ndarray:
        """Computes the passive transformation matrix from this Airplane's geometry
        axes, relative to this Airplane's CG to the first Airplane's geometry axes,
        relative to the first Airplane's CG.

        Computes the transformation chain: G_Cg > B_Cg > E_CgP1 > BP1_CgP1 > GP1_CgP1.
        This transformation matrix is used to position Airplanes relative to one
        another, in problems with more than one Airplane. If this Airplane is the first
        Airplane (where Cg_E_CgP1 = [0, 0, 0]), it returns an identity transformation.

        :param first_airplane: The first Airplane in a problem.
        :return: A (4,4) ndarray of floats representing the passive transformation
            matrix from this Airplane's geometry axes, relative to its CG to the first
            Airplane's geometry axes, relative to its CG.
        """
        # Step 1: G_Cg > B_Cg (geometry axes to body axes: 180-degree rotation about
        # y-axis)
        # This transforms from geometry axes (aft/right/up) to body axes (
        # forward/right/down)
        T_pas_G_Cg_to_B_Cg = _transformations.generate_rot_T(
            angles=np.array([0.0, 180.0, 0.0], dtype=float),
            passive=True,
            intrinsic=True,
            order="xyz",
        )

        # Step 2: B_Cg > E_Cg (body axes to Earth axes using this Airplane's
        # angles_E_to_B_izyx)
        # We need to invert the E_Cg > B_Cg transformation to get B_Cg > E_Cg
        T_pas_E_Cg_to_B_Cg = _transformations.generate_rot_T(
            angles=self.angles_E_to_B_izyx,
            passive=True,
            intrinsic=True,
            order="zyx",
        )
        T_pas_B_Cg_to_E_Cg = _transformations.invert_T_pas(T_pas_E_Cg_to_B_Cg)

        # Step 3: E_Cg > E_CgP1
        # This offsets by Cg_E_CgP1 (the position of this Airplane's CG relative to the
        # first Airplane's CG)
        T_pas_E_Cg_to_E_CgP1 = _transformations.generate_trans_T(
            translations=self.Cg_E_CgP1, passive=True
        )

        # Step 4: E_CgP1 > BP1_CgP1 (Earth axes to the first Airplane's body axes)
        T_pas_E_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
            angles=first_airplane.angles_E_to_B_izyx,
            passive=True,
            intrinsic=True,
            order="zyx",
        )

        # Step 5: BP1_CgP1 > GP1_CgP1 (the first Airplane's body axes to geometry
        # axes: 180° rotation about the y-axis)
        T_pas_BP1_CgP1_to_GP1_CgP1 = _transformations.generate_rot_T(
            angles=np.array([0.0, 180.0, 0.0], dtype=float),
            passive=True,
            intrinsic=True,
            order="xyz",
        )

        # Compose the full passive transformation matrix
        T_pas_G_Cg_to_GP1_CgP1 = _transformations.compose_T_pas(
            T_pas_G_Cg_to_B_Cg,
            T_pas_B_Cg_to_E_Cg,
            T_pas_E_Cg_to_E_CgP1,
            T_pas_E_CgP1_to_BP1_CgP1,
            T_pas_BP1_CgP1_to_GP1_CgP1,
        )

        return T_pas_G_Cg_to_GP1_CgP1

    # TEST: Add more thorough unit tests for this method.
    @staticmethod
    def process_wing_symmetry(wing: wing_mod.Wing) -> list[wing_mod.Wing]:
        """Processes a Wing to determine what type of symmetry it has. If necessary, it
        then modifies the Wing. If type 5 symmetry is detected, it also creates a second
        reflected Wing. Finally, it returns a list of Wings.

        :param wing: The Wing to process for symmetry analysis and potential
            modification.
        :return: The list of processed Wings. For types 1-4 symmetry it contains only
            the one modified Wing, but for type 5 symmetry it contains the modified Wing
            followed by the new reflected Wing. Before returning them, it also calls
            each Wing's generate_mesh method, preparing them for use simulation.
        """
        # Determine if the symmetry plane is coincident with the wing axes' xz-plane.
        # If symmetryNormal_G or symmetryPoint_G_Cg is None, then there is no
        # symmetry and the symmetry plane doesn't exist. Otherwise, the symmetry
        # plane is coincident to the wing axes' xz-plane if Ler_Gs_Cgs lies on the
        # symmetry plane, and if symmetryNormal_G is parallel with WnY_G. We don't
        # need to check types, values, or normalize because this is done in Wing's
        # init method.
        coincident_symmetry_plane = True
        if wing.symmetryPoint_G_Cg is None or wing.symmetryNormal_G is None:
            coincident_symmetry_plane = False
        else:
            # If the symmetry plane exists, we first need to check if its normal
            # vector is parallel with the wing axes' y-axis vector.

            # Actively transform geometry axes' second basis vector (in geometry
            # axes) to this Wing's axes' second basis vector (in geometry axes). We
            # can skip the translation step (step 2) as we are only transforming a
            # direction vector, not a position vector.
            GY_G = np.array([0.0, 1.0, 0.0], dtype=float)
            GsY_G = _transformations.apply_T_to_vectors(
                _transformations.generate_reflect_T(
                    plane_point_A_a=wing.symmetryPoint_G_Cg,
                    plane_normal_A=wing.symmetryNormal_G,
                    passive=False,
                ),
                GY_G,
                has_point=False,
            )
            WnY_G = _transformations.apply_T_to_vectors(
                _transformations.generate_rot_T(
                    wing.angles_Gs_to_Wn_ixyz,
                    passive=False,
                    intrinsic=True,
                    order="xyz",
                ),
                GsY_G,
                has_point=False,
            )

            # If symmetryNormal_G is parallel with WnY_G, their cross product will be
            # the zero vector.
            is_parallel = np.allclose(
                np.cross(wing.symmetryNormal_G, WnY_G),
                np.array([0.0, 0.0, 0.0], dtype=float),
            )

            if not is_parallel:
                coincident_symmetry_plane = False
            else:
                # If the symmetry plane's normal vector and the wing axes y-axis
                # vector are parallel, then the last check for a coincident symmetry
                # plane is to check if the Ler is on the symmetry plane.

                # To do this, we first find the symmetry plane's normal vector (in
                # geometry axes after accounting for symmetry) and the symmetry
                # plane's point (in geometry axes after accounting for symmetry,
                # relative to the CG after accounting for symmetry). As the symmetry
                # plane is defined using these quantities, they don't change after
                # reflection.
                symmetryPoint_Gs_Cgs = wing.symmetryPoint_G_Cg
                symmetryNormal_Gs_Cgs = wing.symmetryNormal_G

                # The leading edge root point is on the symmetry plane if the
                # distance between it and the symmetry plane is zero.
                Ler_on_plane = np.allclose(
                    np.dot(
                        symmetryNormal_Gs_Cgs, (wing.Ler_Gs_Cgs - symmetryPoint_Gs_Cgs)
                    ),
                    0.0,
                )

                if not Ler_on_plane:
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

                reflected_airfoil = airfoil_mod.Airfoil(
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
                    wing_cross_section_mod.WingCrossSection(
                        airfoil=reflected_airfoil,
                        num_spanwise_panels=wing_cross_section.num_spanwise_panels,
                        chord=wing_cross_section.chord,
                        Lp_Wcsp_Lpp=np.copy(wing_cross_section.Lp_Wcsp_Lpp),
                        angles_Wcsp_to_Wcs_ixyz=np.copy(
                            wing_cross_section.angles_Wcsp_to_Wcs_ixyz
                        ),
                        control_surface_symmetry_type=None,
                        control_surface_hinge_point=wing_cross_section.control_surface_hinge_point,
                        control_surface_deflection=reflected_control_surface_deflection,
                        spanwise_spacing=wing_cross_section.spanwise_spacing,
                    )
                )

            reflected_wing = wing_mod.Wing(
                wing_cross_sections=reflected_wing_cross_sections,
                name=f"Reflected {wing.name}",
                Ler_Gs_Cgs=np.copy(wing.Ler_Gs_Cgs),
                angles_Gs_to_Wn_ixyz=np.copy(wing.angles_Gs_to_Wn_ixyz),
                symmetric=False,
                mirror_only=True,
                symmetryNormal_G=np.copy(wing.symmetryNormal_G),
                symmetryPoint_G_Cg=np.copy(wing.symmetryPoint_G_Cg),
                num_chordwise_panels=wing.num_chordwise_panels,
                chordwise_spacing=wing.chordwise_spacing,
            )

            wing.symmetric = False
            wing.mirror_only = False
            wing.symmetryNormal_G = None
            wing.symmetryPoint_G_Cg = None

            # Reset control_surface_symmetry_type to None for Type 1 symmetry.
            for wing_cross_section in wing.wing_cross_sections:
                wing_cross_section.control_surface_symmetry_type = None

            wing.generate_mesh(symmetry_type=1)
            reflected_wing.generate_mesh(symmetry_type=3)
            return [wing, reflected_wing]
