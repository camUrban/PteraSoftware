"""Contains the Wing class.

**Contains the following classes:**

Wing: A class used to contain wings of an Airplane.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pyvista as pv

from .. import _panel, _parameter_validation, _transformations
from . import _meshing
from . import wing_cross_section as wing_cross_section_mod


class Wing:
    """A class used to contain the wings of an Airplane.

    **Contains the following methods:**

    generate_mesh: Generates this Wing's mesh, which finishes the process of preparing
    the Wing to be used in a simulation. It is called by the Wing's parent Airplane,
    after it's determined its symmetry type.

    get_plottable_data: Returns plottable data for this Wing's Airfoils' outlines and
    mean camber lines.

    T_pas_G_Cg_to_Wn_Ler: The passive transformation matrix which maps in homogeneous
    coordinates from geometry axes relative to the CG to wing axes relative to the
    leading edge root point. It is None if the Wing's symmetry type hasn't been defined
    yet.

    T_pas_Wn_Ler_to_G_Cg: The passive transformation matrix which maps in homogeneous
    coordinates from wing axes relative to the leading edge root point to geometry axes
    relative to the CG point. It is None if the Wing's symmetry type hasn't been defined
    yet.

    WnX_G: The wing axes' first basis vector (in geometry axes).

    WnY_G: The wing axes' second basis vector (in geometry axes).

    WnZ_G: The wing axes' third basis vector (in geometry axes).

    children_T_pas_Wn_Ler_to_Wcs_Lp: A list of passive transformation matrices which map
    in homogeneous coordinates from wing axes, relative to the leading edge root point,
    to each of this Wing's WingCrossSection's axes, relative to their respective leading
    points.

    children_T_pas_Wcs_Lp_to_Wn_Ler: A list of passive transformation matrices which map
    in homogeneous coordinates from each of this Wing's WingCrossSection's axes,
    relative to their respective leading points, to wing axes, relative to the leading
    edge root point.

    children_T_pas_G_Cg_to_Wcs_Lp: A list of passive transformation matrices which map
    in homogeneous coordinates from geometry axes, relative to the CG, to each of this
    Wing's WingCrossSection's axes, relative to their respective leading points.

    children_T_pas_Wcs_Lp_to_G_Cg: A list of passive transformation matrices which map
    in homogeneous coordinates from each of this Wing's WingCrossSection's axes,
    relative to their respective leading points, to geometry axes, relative to the CG.

    projected_area: The area of the Wing projected onto the plane defined by the wing
    axes' xy plane.

    wetted_area: The Wing's wetted area.

    average_panel_aspect_ratio: The average aspect ratio of the Wing's Panels.

    span: The Wing's span.

    standard_mean_chord: The Wing's standard mean chord.

    mean_aerodynamic_chord: The Wing's mean aerodynamic chord.

    **Notes:**

    Every Wing has its own axis system, known as wing axes. The user sets the
    relationship between these axes and geometry axes with the Ler_Gs_Cgs and
    angles_Gs_to_Wn_ixyz parameters. However, the steps for transforming a vector from
    geometry axes to wing axes, and the interpretation of the wing axes orientation and
    position relative to an Airplane's geometry axes, also depend on the parameters
    symmetric, mirror_only, symmetryNormal_G, and symmetryPoint_G_Cg. In all cases, the
    order of transformations from geometry axes to wing axes is reflection (if
    applicable), translation, and then rotation.

    There are five symmetry types. Type 1: symmetric=False, mirror_only=False, and the
    symmetry plane must be undefined. Type 2: symmetric=False, mirror_only=True, and the
    symmetry plane is coincident with the wing axes' xz plane. Type 3: symmetric=False,
    mirror_only=True, and the symmetry plane is not coincident with the wing axes' xz
    plane. Type 4: symmetric=True, mirror_only=False, and the symmetry plane is
    coincident with the wing axes' xz plane. Type 5: symmetric=True, mirror_only=False,
    and the symmetry plane is not coincident with the wing axes' xz plane.

    **Citation:**

    Adapted from: geometry.Wing in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/24/2020
    """

    def __init__(
        self,
        wing_cross_sections: list[wing_cross_section_mod.WingCrossSection],
        name: str = "Untitled Wing",
        Ler_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        angles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        symmetric: bool | np.bool_ = False,
        mirror_only: bool | np.bool_ = False,
        symmetryNormal_G: None | np.ndarray | Sequence[float | int] = None,
        symmetryPoint_G_Cg: None | np.ndarray | Sequence[float | int] = None,
        num_chordwise_panels: int = 8,
        chordwise_spacing: str = "cosine",
    ) -> None:
        """The initialization method.

        :param wing_cross_sections: A list of WingCrossSections representing the wing's
            cross sections in order from root to tip. It must contain at least two
            WingCrossSections.
        :param name: A sensible name for the Wing. The default is "Untitled Wing".
        :param Ler_Gs_Cgs: An array-like object of 3 numbers (int or float) representing
            the position of the origin of this Wing's axes (in geometry axes after
            accounting for symmetry, relative to the CG after accounting for symmetry).
            Can be a tuple, list, or ndarray. Values are converted to floats internally.
            The units are meters. The default is (0.0, 0.0, 0.0).
        :param angles_Gs_to_Wn_ixyz: An array-like object of 3 numbers (int or float)
            representing the angle vector that defines the orientation of this Wing's
            axes relative to the geometry axes (after accounting for symmetry). Can be a
            tuple, list, or ndarray. Values are converted to floats internally. All
            angles must be in the range [-90, 90] degrees. Rotations are intrinsic, and
            proceed in the xy'z" order. The units are degrees. The default is (0.0, 0.0,
            0.0).
        :param symmetric: Set this to True if the Wing's geometry should be mirrored
            across the symmetry plane while retaining the non mirrored side. If
            mirror_only is True, symmetric must be False. If symmetric is True, then
            neither symmetryNormal_G nor symmetryPoint_G_Cg can be None. If the symmetry
            plane is coincident with this Wing's axes' xz plane, the mirrored and non
            mirrored geometry will be meshed as a single wing. If not, this Wing's
            Airplane will automatically create another Wing with the mirrored geometry,
            modify both Wings' parameters, and add the reflected Wing to its list of
            Wings immediately following this one. For more details on that process, and
            how this parameter interacts with symmetryNormal_G, symmetryPoint_G_Cg, and
            mirror_only, see the class docstring. Can be a bool or a numpy bool and will
            be converted internally to a bool. The default is False.
        :param mirror_only: Set this to True if the Wing's geometry should be reflected
            about the symmetry plane without retaining the non reflected geometry. If
            symmetric is True, mirror_only must be False. If mirror_only is True, then
            neither symmetryNormal_G nor symmetryPoint_G_Cg can be None. For more
            details on how this parameter interacts with symmetryNormal_G,
            symmetryPoint_G_Cg, and symmetric, see the class docstring. Can be a bool or
            a numpy bool and will be converted internally to a bool. The default is
            False.
        :param symmetryNormal_G: None, or an array-like of 3 numbers (int or float)
            representing the unit normal vector (in geometry axes) that, together with
            symmetryPoint_G_Cg, defines the plane used for symmetry or mirroring. Can be
            None, or a tuple, list, or ndarray. If not None, values are converted to
            floats and normalized internally. Note that reversing the normal direction
            (using the antiparallel vector) defines the same plane and produces the same
            result. This value must be None if both symmetric and mirror_only are False,
            and cannot be None if either are True. For more details on how this
            parameter interacts with symmetryPoint_G_Cg, symmetric, and mirror_only, see
            the class docstring. The default is None.
        :param symmetryPoint_G_Cg: None or an array-like object of 3 numbers (int or
            float) representing a point (in geometry axes, relative to the CG) that,
            along with symmetryNormal_G, defines the location of the plane about which
            symmetry or mirroring is applied. Can be None, or a list, tuple, or ndarray.
            If not None, values are converted to floats internally. This value must be
            None if both symmetric and mirror_only are False, and cannot be None if
            either are True. For more details on how this parameter interacts with
            symmetryNormal_G, symmetric, and mirror_only, see the class docstring. The
            units are meters. The default is None.
        :param num_chordwise_panels: The number of chordwise panels to be used on this
            Wing, which must be set to a positive integer. The default is 8.
        :param chordwise_spacing: The type of spacing between the Wing's chordwise
            panels. Can be "cosine" or "uniform". Using cosine spacing is highly
            recommended for steady simulations and uniform spacing is highly recommended
            for unsteady simulations. The default is "cosine".
        :return: None
        """
        # Validate wing_cross_sections.
        wing_cross_sections = _parameter_validation.non_empty_list_return_list(
            wing_cross_sections, "wing_cross_sections"
        )
        num_wing_cross_sections = len(wing_cross_sections)
        if num_wing_cross_sections < 2:
            raise ValueError("wing_cross_sections must contain at least two elements.")
        for wing_cross_section_id, wing_cross_section in enumerate(wing_cross_sections):
            if not isinstance(
                wing_cross_section, wing_cross_section_mod.WingCrossSection
            ):
                raise TypeError(
                    "Every element in wing_cross_sections must be a WingCrossSection."
                )
            if wing_cross_section_id == 0:
                # Validate root WingCrossSection constraints.
                wing_cross_section.validate_root_constraints()
            elif wing_cross_section_id == num_wing_cross_sections - 1:
                # Validate tip WingCrossSection constraints.
                wing_cross_section.validate_tip_constraints()
            else:
                wing_cross_section.validate_mid_constraints()
            # Set the validated flag for this WingCrossSection.
            wing_cross_section.validated = True
        self.wing_cross_sections: list[wing_cross_section_mod.WingCrossSection] = (
            wing_cross_sections
        )

        # Validate name and Ler_Gs_Cgs.
        self.name = _parameter_validation.str_return_str(name, "name")
        self.Ler_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            Ler_Gs_Cgs, "Ler_Gs_Cgs"
        )

        # Validate angles_Gs_to_Wn_ixyz.
        angles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                angles_Gs_to_Wn_ixyz, "angles_Gs_to_Wn_ixyz"
            )
        )
        if not np.all((-90.0 <= angles_Gs_to_Wn_ixyz) & (angles_Gs_to_Wn_ixyz <= 90.0)):
            raise ValueError(
                "All elements of angles_Gs_to_Wn_ixyz must lie in the range [-90, "
                "90] degrees."
            )
        self.angles_Gs_to_Wn_ixyz = angles_Gs_to_Wn_ixyz

        # Validate symmetric and mirror_only.
        symmetric = _parameter_validation.boolLike_return_bool(symmetric, "symmetric")
        mirror_only = _parameter_validation.boolLike_return_bool(
            mirror_only, "mirror_only"
        )
        if symmetric and mirror_only:
            raise ValueError("symmetric and mirror_only cannot both be True.")
        self.symmetric = symmetric
        self.mirror_only = mirror_only

        # Validate symmetryNormal_G and symmetryPoint_G_Cg.
        if self.symmetric or self.mirror_only:
            if symmetryNormal_G is None:
                raise ValueError(
                    "symmetryNormal_G cannot be None when symmetric or mirror_only is "
                    "True."
                )
            symmetryNormal_G = (
                _parameter_validation.threeD_number_vectorLike_return_float_unit_vector(
                    symmetryNormal_G, "symmetryNormal_G"
                )
            )
            if symmetryPoint_G_Cg is None:
                raise ValueError(
                    "symmetryPoint_G_Cg cannot be None when symmetric or mirror_only "
                    "is True."
                )
            symmetryPoint_G_Cg = (
                _parameter_validation.threeD_number_vectorLike_return_float(
                    symmetryPoint_G_Cg, "symmetryPoint_G_Cg"
                )
            )
        else:
            if symmetryNormal_G is not None:
                raise ValueError(
                    "symmetryNormal_G must be None when both symmetric and "
                    "mirror_only are False."
                )
            if symmetryPoint_G_Cg is not None:
                raise ValueError(
                    "symmetryPoint_G_Cg must be None when both symmetric and "
                    "mirror_only are False."
                )
        self.symmetryNormal_G = symmetryNormal_G
        self.symmetryPoint_G_Cg = symmetryPoint_G_Cg

        # Validate num_chordwise_panels and chordwise_spacing.
        self.num_chordwise_panels = _parameter_validation.int_in_range_return_int(
            num_chordwise_panels,
            "num_chordwise_panels",
            min_val=1,
            min_inclusive=True,
        )
        if chordwise_spacing not in ["cosine", "uniform"]:
            raise ValueError('chordwise_spacing must be "cosine" or "uniform".')
        self.chordwise_spacing = chordwise_spacing

        # These attributes will be initialized or populated once this Wing's parent
        # Airplane calls generate_mesh.
        self.symmetry_type: int | None = None
        self.num_spanwise_panels: int | None = None
        self.num_panels: int | None = None
        self.panels: np.ndarray | None = None
        self.gridWrvp_GP1_CgP1: np.ndarray | None = None
        self.wake_ring_vortices: np.ndarray | None = None

    def generate_mesh(self, symmetry_type: int) -> None:
        """Generates this Wing's mesh, which finishes the process of preparing the Wing
        to be used in a simulation. It is called by the Wing's parent Airplane, after
        it's determined its symmetry type.

        :param symmetry_type: The symmetry type of this Wing as an int from 1-4. See the
            class docstring for details on how to interpret the symmetry types.
        :return: None
        """
        # Validate and apply symmetry_type. 5 isn't a valid symmetry type, because
        # the parent Airplane should have modified a Wing that initially had type 5
        # symmetry to have type 1 symmetry, and then made a new reflected Wing with
        # type 3 symmetry.
        self.symmetry_type = _parameter_validation.int_in_range_return_int(
            symmetry_type,
            "symmetry_type",
            min_val=1,
            min_inclusive=True,
            max_val=4,
            max_inclusive=True,
        )

        # Set this Wing's children WingCrossSections' symmetry type parameters.
        for wing_cross_section in self.wing_cross_sections:
            wing_cross_section.symmetry_type = self.symmetry_type

        # Find the number of spanwise panels on the wing by adding each cross
        # section's number of spanwise panels. Exclude the last cross section's
        # number of spanwise panels as this is irrelevant. If the wing has type 4
        # symmetry multiply the summation by two.
        self.num_spanwise_panels = 0
        for wing_cross_section in self.wing_cross_sections[:-1]:
            assert wing_cross_section.num_spanwise_panels is not None
            self.num_spanwise_panels += wing_cross_section.num_spanwise_panels
        if self.symmetry_type == 4:
            self.num_spanwise_panels *= 2

        # Calculate the number of panels on this wing.
        self.num_panels = self.num_spanwise_panels * self.num_chordwise_panels

        # Initialize empty arrays to hold this wing's wake RingVortices and its wake
        # RingVortex points.
        self.wake_ring_vortices = np.zeros((0, self.num_spanwise_panels), dtype=object)
        self.gridWrvp_GP1_CgP1 = np.empty(
            (0, self.num_spanwise_panels + 1, 3), dtype=float
        )

        # Generate the wing's mesh, which populates the Panels attribute.
        _meshing.mesh_wing(self)

    # TEST: Consider adding unit tests for this method.
    def get_plottable_data(
        self, show: bool | np.bool_ = False
    ) -> list[list[np.ndarray]] | None:
        """Returns plottable data for this Wing's Airfoils' outlines and mean camber
        lines.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: None if the Wing's symmetry type hasn't been set yet, or if show is
            True. Otherwise, returns a list of two lists, each containing one ndarray
            for every one of this Wing's Airfoils. These ndarrays represent points on
            each Airfoil's outline and mean camber lines, respectively. The points are
            in wing axes, relative to the leading edge root point. The units are in
            meters.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        # If this Wing hasn't had its symmetry type set, return None.
        if self.symmetry_type is None:
            return None

        airfoilOutlines_Wn_ler = []
        airfoilMcls_Wn_ler = []
        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections
        ):
            plottable_data = wing_cross_section.get_plottable_data(show=False)
            assert plottable_data is not None

            [airfoilOutline_Wcs_lp, airfoilMcl_Wcs_lp] = plottable_data

            T_pas_Wcs_Lp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[
                wing_cross_section_id
            ]

            airfoilOutline_Wn_ler = _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_Wn_Ler, airfoilOutline_Wcs_lp, has_point=True
            )
            airfoilMcl_Wn_ler = _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_Wn_Ler, airfoilMcl_Wcs_lp, has_point=True
            )

            airfoilOutlines_Wn_ler.append(airfoilOutline_Wn_ler)
            airfoilMcls_Wn_ler.append(airfoilMcl_Wn_ler)

        if not show:
            return [airfoilOutlines_Wn_ler, airfoilMcls_Wn_ler]

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

        plotter.add_actor(AxesGCg)  # type: ignore[arg-type]

        _T_pas_G_Cg_to_Wn_Ler = self.T_pas_G_Cg_to_Wn_Ler
        assert _T_pas_G_Cg_to_Wn_Ler is not None

        AxesWLerWcs1Lp1_G_Cg = pv.AxesAssembly(
            x_label="WX@Ler/Wcs1XLp1",
            y_label="WY@Ler/Wcs1YLp1",
            z_label="WZ@Ler/Wcs1ZLp1",
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
            user_matrix=np.linalg.inv(_T_pas_G_Cg_to_Wn_Ler),
            # user_matrix=wingAxes_T_act,
            name="W/Wcs1",
            shaft_type="cylinder",
            shaft_radius=0.025,
            shaft_length=(0.8, 0.8, 0.8),
            tip_type="cone",
            tip_radius=0.1,
            tip_length=(0.2, 0.2, 0.2),
            symmetric_bounds=False,
        )

        plotter.add_actor(AxesWLerWcs1Lp1_G_Cg)  # type: ignore[arg-type]

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections
        ):
            airfoilOutline_Wn_ler = airfoilOutlines_Wn_ler[wing_cross_section_id]
            airfoilMcl_Wn_ler = airfoilMcls_Wn_ler[wing_cross_section_id]

            airfoilOutline_G_Cg = _transformations.apply_T_to_vectors(
                _T_pas_Wn_Ler_to_G_Cg, airfoilOutline_Wn_ler, has_point=True
            )
            airfoilMcl_G_Cg = _transformations.apply_T_to_vectors(
                _T_pas_Wn_Ler_to_G_Cg, airfoilMcl_Wn_ler, has_point=True
            )

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
                    x_label=f"Wcs{wing_cross_section_num}X@Lp{wing_cross_section_num}",
                    y_label=f"Wcs{wing_cross_section_num}Y@Lp{wing_cross_section_num}",
                    z_label=f"Wcs{wing_cross_section_num}Z@Lp{wing_cross_section_num}",
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
                        self.children_T_pas_G_Cg_to_Wcs_Lp[wing_cross_section_id]
                    ),
                    name=f"Wcs{wing_cross_section_id}",
                    shaft_type="cylinder",
                    shaft_radius=0.025,
                    shaft_length=(0.8, 0.8, 0.8),
                    tip_type="cone",
                    tip_radius=0.1,
                    tip_length=(0.2, 0.2, 0.2),
                    symmetric_bounds=False,
                )

                plotter.add_actor(AxesWcsLp_G_Cg)  # type: ignore[arg-type]

        if self.panels is not None:
            # Initialize empty arrays to hold the Panels' vertices and faces.
            panel_vertices = np.empty((0, 3), dtype=float)
            panel_faces = np.empty(0, dtype=int)

            # Initialize a variable to keep track of how many Panels' data has been
            # added to the arrays.
            panel_num = 0

            # Unravel the Wing's Panel matrix and iterate through it.
            panels = np.ravel(self.panels)
            for panel in panels:
                # Stack this Panel's vertices and faces.
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
                # vertices and faces.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # Update the number of previous Panels.
                panel_num += 1

                # Convert the Panel vertices and faces to PolyData.
                panel_surfaces = pv.PolyData(panel_vertices, panel_faces)

                # Add the Panels to the plotter.
                plotter.add_mesh(
                    panel_surfaces,
                    show_edges=True,
                    color="chartreuse",
                    smooth_shading=False,
                )

        plotter.enable_parallel_projection()  # type: ignore[call-arg]

        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )

        return None

    @property
    def T_pas_G_Cg_to_Wn_Ler(self) -> None | np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        geometry axes relative to the CG to wing axes relative to the leading edge root
        point. Is None if the Wing's symmetry type hasn't been defined yet.

        :return: A (4,4) ndarray of floats representing the transformation matrix or
            None if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Step 1: Create T_reflect_pas_G_Cg_to_Gs_Cgs, which maps from which maps in
        # homogeneous coordinates from geometry axes relative to the CG to reflected
        # geometry axes (after accounting for symmetry) relative to the CG (after
        # accounting for symmetry). This is the reflection step. Only apply reflection
        # for mirror-only Wings (types 2 and 3), not for symmetric Wings (type 4).
        if self.symmetry_type in (2, 3):
            assert self.symmetryPoint_G_Cg is not None
            assert self.symmetryNormal_G is not None
            T_reflect_pas_G_Cg_to_Gs_Cgs = _transformations.generate_reflect_T(
                plane_point_A_a=self.symmetryPoint_G_Cg,
                plane_normal_A=self.symmetryNormal_G,
                passive=True,
            )
        else:
            T_reflect_pas_G_Cg_to_Gs_Cgs = np.eye(4, dtype=float)

        # Step 2: Create T_trans_pas_Gs_Cgs_to_Gs_Ler, which maps in homogeneous
        # coordinates from geometry axes (after accounting for symmetry) relative to
        # the CG (after accounting for symmetry) to geometry axes (after accounting
        # for symmetry) relative to the leading edge root point. This is the
        # translation step.
        T_trans_pas_Gs_Cgs_to_Gs_Ler = _transformations.generate_trans_T(
            self.Ler_Gs_Cgs, passive=True
        )

        # Step 3: Create T_rot_pas_Gs_to_Wn, which maps in homogeneous coordinates
        # from geometry axes (after accounting for symmetry) to wing axes. This is
        # the rotation step.
        T_rot_pas_Gs_to_Wn = _transformations.generate_rot_T(
            self.angles_Gs_to_Wn_ixyz, passive=True, intrinsic=True, order="xyz"
        )

        return _transformations.compose_T_pas(
            T_reflect_pas_G_Cg_to_Gs_Cgs,
            T_trans_pas_Gs_Cgs_to_Gs_Ler,
            T_rot_pas_Gs_to_Wn,
        )

    @property
    def T_pas_Wn_Ler_to_G_Cg(self) -> None | np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        wing axes relative to the leading edge root point to geometry axes relative to
        the CG. Is None if the Wing's symmetry type hasn't been defined yet.

        :return: A (4,4) ndarray of floats representing the transformation matrix or
            None if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        _T_pas_G_Cg_to_Wn_Ler = self.T_pas_G_Cg_to_Wn_Ler
        assert _T_pas_G_Cg_to_Wn_Ler is not None

        return _transformations.invert_T_pas(_T_pas_G_Cg_to_Wn_Ler)

    @property
    def WnX_G(self) -> None | np.ndarray:
        """The wing axes' first basis vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the wing axes' first basis vector
            (in geometry axes) or None if the Wing's symmetry type hasn't been defined
            yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        WnX_Wn = np.array([1.0, 0.0, 0.0])

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        return _transformations.apply_T_to_vectors(
            _T_pas_Wn_Ler_to_G_Cg, WnX_Wn, has_point=False
        )

    @property
    def WnY_G(self) -> None | np.ndarray:
        """The wing axes' second basis vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the wing axes' second basis
            vector (in geometry axes) or None if the Wing's symmetry type hasn't been
            defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        WnY_Wn = np.array([0.0, 1.0, 0.0])

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        return _transformations.apply_T_to_vectors(
            _T_pas_Wn_Ler_to_G_Cg, WnY_Wn, has_point=False
        )

    @property
    def WnZ_G(self) -> None | np.ndarray:
        """The wing axes' third basis vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the wing axes' third basis vector
            (in geometry axes) or None if the Wing's symmetry type hasn't been defined
            yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        WnZ_Wn = np.array([0.0, 0.0, 1.0])

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        return _transformations.apply_T_to_vectors(
            _T_pas_Wn_Ler_to_G_Cg, WnZ_Wn, has_point=False
        )

    @property
    def children_T_pas_Wn_Ler_to_Wcs_Lp(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from wing axes, relative to the leading edge root point, to each of
        this Wing's WingCrossSection's axes, relative to their respective leading
        points.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """

        return [
            _transformations.compose_T_pas(
                *(
                    _assert_T_not_none(wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp)
                    for wing_cross_section in self.wing_cross_sections[: i + 1]
                )
            )
            for i in range(len(self.wing_cross_sections))
        ]

    @property
    def children_T_pas_Wcs_Lp_to_Wn_Ler(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from each of this Wing's WingCrossSection's axes, relative to their
        respective leading points, to wing axes, relative to the leading edge root
        point.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """

        return [
            _transformations.invert_T_pas(self.children_T_pas_Wn_Ler_to_Wcs_Lp[i])
            for i in range(len(self.wing_cross_sections))
        ]

    @property
    def children_T_pas_G_Cg_to_Wcs_Lp(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from geometry axes, relative to the CG, to each of this Wing's
        WingCrossSection's axes, relative to their respective leading points.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """

        _T_pas_G_Cg_to_Wn_Ler = self.T_pas_G_Cg_to_Wn_Ler
        assert _T_pas_G_Cg_to_Wn_Ler is not None

        return [
            _transformations.compose_T_pas(
                _T_pas_G_Cg_to_Wn_Ler, self.children_T_pas_Wn_Ler_to_Wcs_Lp[i]
            )
            for i in range(len(self.wing_cross_sections))
        ]

    @property
    def children_T_pas_Wcs_Lp_to_G_Cg(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from each of this Wing's WingCrossSection's axes, relative to their
        respective leading points, to geometry axes, relative to the CG.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """

        return [
            _transformations.invert_T_pas(self.children_T_pas_G_Cg_to_Wcs_Lp[i])
            for i in range(len(self.wing_cross_sections))
        ]

    @property
    def projected_area(self) -> None | float:
        """The area of the Wing projected onto the plane defined by the wing axes' xy
        plane.

        **Notes:**

        If the Wing is symmetric and continuous, the area of the mirrored half is
        included.

        :return: The projected area of the Wing. It has units of square meters. If the
            Wing hasn't been meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self.panels is None:
            return None

        projected_area = 0.0

        # Iterate through the chordwise and spanwise indices of the Panels and add
        # their area to the total projected area.
        assert self.num_spanwise_panels is not None
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                this_panel: _panel.Panel = self.panels[
                    chordwise_location, spanwise_location
                ]

                thisWnZ_G = self.WnZ_G
                assert thisWnZ_G is not None

                projected_area += this_panel.calculate_projected_area(thisWnZ_G)

        return projected_area

    @property
    def wetted_area(self) -> None | float:
        """The Wing's wetted area.

        **Notes:**

        If the Wing is symmetric and continuous, the area of the mirrored half is
        included.

        :return: The wetted area of the Wing. It has units of square meters. If the Wing
            hasn't been meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self.panels is None:
            return None

        wetted_area = 0.0

        # Iterate through the chordwise and spanwise indices of the panels and add
        # their area to the total wetted area.
        assert self.num_spanwise_panels is not None
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                this_panel: _panel.Panel = self.panels[
                    chordwise_location, spanwise_location
                ]
                wetted_area += this_panel.area

        return wetted_area

    # TEST: Consider adding unit tests for this method.
    @property
    def average_panel_aspect_ratio(self) -> None | float:
        """The average aspect ratio of the Wing's Panels.

        :return: The average aspect ratio of the Wing's Panels. If the Wing hasn't been
            meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self.panels is None:
            return None

        aspect_ratio_sum = 0.0

        # Iterate through the chordwise and spanwise indices of the Panels and sum
        # all the Panels' aspect ratios.
        assert self.num_spanwise_panels is not None
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                this_panel: _panel.Panel = self.panels[
                    chordwise_location, spanwise_location
                ]
                aspect_ratio_sum += this_panel.aspect_ratio

        assert self.num_panels is not None
        average_aspect_ratio = aspect_ratio_sum / self.num_panels

        return average_aspect_ratio

    @property
    def span(self) -> None | float:
        """The Wing's span.

        **Notes:**

        The span is derived by first finding the vector connecting the leading edges of
        the root and tip WingCrossSections. Then, this vector is projected onto the wing
        axes' second basis vector. The span is defined as the magnitude of this
        projection.

        If the Wing is symmetric and continuous, this method includes the span of the
        mirrored half.

        :return: The Wing's span. It has units of meters. None is returned if the Wing's
            symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        tipLp_Wcsp_Lpp = self.wing_cross_sections[-1].Lp_Wcsp_Lpp

        tip_T_pas_Wcsp_Lpp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[-2]

        tipLp_Wn_Ler = _transformations.apply_T_to_vectors(
            tip_T_pas_Wcsp_Lpp_to_Wn_Ler, tipLp_Wcsp_Lpp, has_point=True
        )

        # Project the tip position onto the wing axes' y direction (spanwise direction)
        projectedTipLp_Wn_Ler = np.dot(
            tipLp_Wn_Ler, np.array([0.0, 1.0, 0.0])
        ) * np.array([0.0, 1.0, 0.0])

        span = float(np.linalg.norm(projectedTipLp_Wn_Ler))

        # If the wing is symmetric and continuous, multiply the span by two.
        if self.symmetry_type == 4:
            span *= 2

        return span

    @property
    def standard_mean_chord(self) -> None | float:
        """The Wing's standard mean chord.

        **Notes:**

        The standard mean chord is defined as the projected area divided by the span.
        See their respective methods for the definitions of span and projected area.

        :return: The standard mean chord of the Wing. It has units of meters. None is
            returned if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        _projected_area = self.projected_area
        assert _projected_area is not None

        _span = self.span
        assert _span is not None

        return _projected_area / _span

    @property
    def mean_aerodynamic_chord(self) -> None | float:
        """The Wing's mean aerodynamic chord.

        :return: The mean aerodynamic chord of the Wing. It has units of meters. None is
            returned if the Wing's symmetry type hasn't been defined yet.
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
        integral = 0.0

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
            Lp_Wcs_Lp = np.array([0.0, 0.0, 0.0])

            T_pas_Wcs_Lp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[
                wing_cross_section_id
            ]

            Lp_Wn_Ler = _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_Wn_Ler, Lp_Wcs_Lp, has_point=True
            )

            # Calculate next WingCrossSection's position in wing axes
            nextLp_nextWcs_nextLp = np.array([0.0, 0.0, 0.0])

            T_pas_nextWcs_nextLp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[
                wing_cross_section_id + 1
            ]

            nextLp_Wn_Ler = _transformations.apply_T_to_vectors(
                T_pas_nextWcs_nextLp_to_Wn_Ler, nextLp_nextWcs_nextLp, has_point=True
            )

            # Find the section vector and project it onto spanwise direction (wing axes y direction)
            nextLp_Wn_Lp = nextLp_Wn_Ler - Lp_Wn_Ler

            nextLpProj_Wn_Lp = np.dot(
                nextLp_Wn_Lp, np.array([0.0, 1.0, 0.0])
            ) * np.array([0.0, 1.0, 0.0])

            section_span = float(np.linalg.norm(nextLpProj_Wn_Lp))

            # Each Wing section is, by definition, trapezoidal (at least when
            # projected on to the wing axes' xy plane). For a trapezoid,
            # the integral from the cited equation can be shown to evaluate to the
            # following.
            integral += (
                section_span * (chord**2 + chord * next_chord + next_chord**2) / 3
            )

        _projected_area = self.projected_area
        assert _projected_area is not None

        # Multiply the integral's value by the coefficients from the cited equation.
        # Double if the wing is symmetric and continuous.
        if self.symmetry_type == 4:
            return 2 * integral / _projected_area
        return integral / _projected_area


def _assert_T_not_none(T: np.ndarray | None) -> np.ndarray:
    """Assert that a transformation matrix is not None and return it.

    :param T: None, or a (4,4) ndarray of floats representing the transformation matrix.
    :return: A (4,4) ndarray of floats representing the transformation matrix
        (guaranteed not to be None).
    """
    assert T is not None
    return T
