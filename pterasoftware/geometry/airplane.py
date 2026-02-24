"""Contains the Airplane class.

**Contains the following classes:**

Airplane: A class used to contain airplanes.

**Contains the following functions:**

None
"""

from __future__ import annotations

import copy
import time
from collections.abc import Sequence
from typing import Any, cast

import numpy as np
import pyvista as pv
import webp

from .. import _parameter_validation, _transformations
from . import wing as wing_mod
from . import wing_cross_section as wing_cross_section_mod


class Airplane:
    """A class used to contain airplanes.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this Airplane, preserving mesh geometry but
    resetting solver state.

    deep_copy_with_Cg_GP1_CgP1: Creates a deep copy of this Airplane with a different
    Cg_GP1_CgP1 position.

    num_panels: The total number of Panels across all Wings.

    T_pas_G_Cg_to_GP1_CgP1: The passive transformation matrix from this Airplane's
    geometry axes, relative to this Airplane's CG to the first Airplane's geometry axes,
    relative to the first Airplane's CG.

    draw: Draws the 3D geometry of this Airplane.

    get_plottable_data: Returns plottable data for this Airplane's Airfoils' outlines
    and mean camber lines.

    validate_first_airplane_constraints: Validates that the first Airplane in a
    simulation has Cg_GP1_CgP1 set to zeros.

    process_wing_symmetry: Processes a Wing to determine what type of symmetry it has.
    If necessary, it then modifies the Wing. If type 5 symmetry is detected, it also
    creates a second reflected Wing. Finally, it returns a list of Wings.

    **Notes:**

    The Airplane class is responsible for: (1) Defining the local body axes and geometry
    axes, (2) managing Wings and their coordinate transformations, (3) processing
    symmetric Wings and converting them to separate wings when the symmetry plane is not
    coincident with the Wing's axes xz plane (type 5 symmetry), and (4) providing
    reference dimensions for aerodynamic calculations.

    Every Airplane has a body axis system, where +x points forward along fuselage, +y
    points to the right (starboard direction), and +z points downward (completing a
    right-handed system).

    Every Airplane also has a geometry axis system, where +x points aft along fuselage,
    +y points to the right (starboard direction), and +z points upward (completing a
    right-handed system).

    Immutable attributes (wings, name, Cg_GP1_CgP1, weight, s_ref, c_ref, and b_ref) are
    set during initialization and cannot be modified afterward. The numpy array
    Cg_GP1_CgP1 is made read only to prevent in place mutation. The wings attribute is
    stored as a tuple to prevent external mutation.

    Derived properties (num_panels and T_pas_G_Cg_to_GP1_CgP1) are lazily evaluated and
    cached since they depend only on immutable attributes.

    The forces_W, forceCoefficients_W, moments_W_CgP1, and momentCoefficients_W_CgP1
    attributes remain mutable as they are set by the solver during simulation.

    **Citation:**

    Adapted from: geometry.Airplane in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/23/2020
    """

    def __init__(
        self,
        wings: list[wing_mod.Wing],
        name: str = "Untitled Airplane",
        Cg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        weight: float | int = 0.0,
        s_ref: float | int | None = None,
        c_ref: float | int | None = None,
        b_ref: float | int | None = None,
    ) -> None:
        """The initialization method.

        :param wings: A list of the airplane's wings defined as Wings. It must contain
            at least one Wing. Wings with symmetric=True and non coincident symmetry
            planes will be automatically processed into separate Wings during
            initialization (type 5 symmetry).
        :param name: A sensible name for your airplane. The default is "Untitled
            Airplane".
        :param Cg_GP1_CgP1: An array-like object of 3 numbers representing the position
            of this Airplane's CG (in the first Airplane's geometry axes, relative to
            the first Airplane's CG). Can be a list, tuple, or ndarray. Values are
            converted to floats internally. For the first Airplane in a simulation, this
            must be equivalent to (0.0, 0.0, 0.0) by definition. The units are in
            meters. The default is (0.0, 0.0, 0.0).
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
        # Initialize the immutable attributes. Set those that are numpy arrays to be
        # read only. Store wings as a tuple to prevent external mutation.
        wings = _parameter_validation.non_empty_list_return_list(wings, "wings")
        processed_wings: list[wing_mod.Wing] = []
        for wing in wings:
            if not isinstance(wing, wing_mod.Wing):
                raise TypeError("Every element in wings must be a Wing")
            processed_wings.extend(self.process_wing_symmetry(wing))
        self._wings = tuple(processed_wings)

        self._name = _parameter_validation.str_return_str(name, "name")

        self._Cg_GP1_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            Cg_GP1_CgP1, "Cg_GP1_CgP1"
        )
        self._Cg_GP1_CgP1.flags.writeable = False

        self._weight = _parameter_validation.number_in_range_return_float(
            weight,
            "weight",
            min_val=0.0,
            min_inclusive=True,
        )

        # If any of the passed reference dimensions are None, set them to first Wing's
        # corresponding reference. Otherwise, set them to the passed dimension after
        # checking that it is valid.
        if s_ref is None:
            first_wing_projected_area = self._wings[0].projected_area
            if first_wing_projected_area is None:
                raise ValueError(
                    "s_ref was not provided and the first Wing's projected_area is "
                    "None. Either provide an explicit s_ref or ensure the first Wing "
                    "is meshed."
                )
            self._s_ref = first_wing_projected_area
        else:
            self._s_ref = _parameter_validation.number_in_range_return_float(
                s_ref, "s_ref", min_val=0.0, min_inclusive=False
            )
        if c_ref is None:
            first_wing_mean_aerodynamic_chord = self._wings[0].mean_aerodynamic_chord
            if first_wing_mean_aerodynamic_chord is None:
                raise ValueError(
                    "c_ref was not provided and the first Wing's "
                    "mean_aerodynamic_chord is None. Either provide an explicit c_ref "
                    "or ensure the first Wing is meshed."
                )
            self._c_ref = first_wing_mean_aerodynamic_chord
        else:
            self._c_ref = _parameter_validation.number_in_range_return_float(
                c_ref, "c_ref", min_val=0.0, min_inclusive=False
            )
        if b_ref is None:
            first_wing_span = self._wings[0].span
            if first_wing_span is None:
                raise ValueError(
                    "b_ref was not provided and the first Wing's span is None. Either "
                    "provide an explicit b_ref or ensure the first Wing is meshed."
                )
            self._b_ref = first_wing_span
        else:
            self._b_ref = _parameter_validation.number_in_range_return_float(
                b_ref, "b_ref", min_val=0.0, min_inclusive=False
            )

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._num_panels: int | None = None
        self._T_pas_G_Cg_to_GP1_CgP1: np.ndarray | None = None

        # Initialize mutable attributes to hold the forces, moments, force
        # coefficients, and moment coefficients this Airplane experiences.
        self.forces_W: np.ndarray | None = None
        self.forceCoefficients_W: np.ndarray | None = None
        self.moments_W_CgP1: np.ndarray | None = None
        self.momentCoefficients_W_CgP1: np.ndarray | None = None

    # --- Deep copy methods ---
    def __deepcopy__(self, memo: dict) -> Airplane:
        """Creates a deep copy of this Airplane, preserving mesh geometry but resetting
        solver state.

        The copy preserves: (1) Wings tuple (each Wing is deep copied, preserving mesh
        and Panels) (2) Airplane parameters (name, Cg_GP1_CgP1, weight, reference
        dimensions), and (3) cached derived properties (num_panels,
        T_pas_G_Cg_to_GP1_CgP1).

        The copy resets to None: (1) loads and load coefficients.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new Airplane with preserved mesh geometry and reset solver state.
        """
        # Create a new Airplane instance without calling __init__ to avoid redundant
        # validation and Wing symmetry processing.
        new_airplane = object.__new__(Airplane)

        # Store this Airplane in memo to handle potential circular references.
        memo[id(self)] = new_airplane

        # Deep copy the Wings into a new tuple.
        new_airplane._wings = tuple(copy.deepcopy(wing, memo) for wing in self._wings)

        # Copy immutable attributes. For those that are numpy arrays, make the copies
        # read only.
        new_airplane._name = self._name
        new_airplane._Cg_GP1_CgP1 = self._Cg_GP1_CgP1.copy()
        new_airplane._Cg_GP1_CgP1.flags.writeable = False
        new_airplane._weight = self._weight
        new_airplane._s_ref = self._s_ref
        new_airplane._c_ref = self._c_ref
        new_airplane._b_ref = self._b_ref

        # Copy cached derived properties. For those that are numpy arrays, make the
        # copies read only.
        new_airplane._num_panels = self._num_panels
        if self._T_pas_G_Cg_to_GP1_CgP1 is not None:
            new_airplane._T_pas_G_Cg_to_GP1_CgP1 = self._T_pas_G_Cg_to_GP1_CgP1.copy()
            new_airplane._T_pas_G_Cg_to_GP1_CgP1.flags.writeable = False
        else:
            new_airplane._T_pas_G_Cg_to_GP1_CgP1 = None

        # Reset loads and load coefficients to None (solver will compute these).
        new_airplane.forces_W = None
        new_airplane.forceCoefficients_W = None
        new_airplane.moments_W_CgP1 = None
        new_airplane.momentCoefficients_W_CgP1 = None

        return new_airplane

    def deep_copy_with_Cg_GP1_CgP1(
        self, new_Cg_GP1_CgP1: np.ndarray | Sequence[float | int]
    ) -> Airplane:
        """Creates a deep copy of this Airplane with a different Cg_GP1_CgP1 position.

        This method is used by AirplaneMovement to create Airplanes at different time
        steps that share the same geometry but have different positions in the
        formation. It maintains immutability by returning a new Airplane rather than
        modifying the existing one.

        Only Cg_GP1_CgP1 and its derived cache (_T_pas_G_Cg_to_GP1_CgP1) need to differ
        from a standard deep copy because (1) Wing geometry (Ler_Gs_Cgs, panels, etc.)
        is defined relative to this Airplane's own CG, not the formation position, so it
        remains valid, (2) Panel local coordinates (_G_Cg) are independent of formation
        position (global coordinates (_GP1_CgP1) are reset to None by Panel's
        __deepcopy__ and will be recomputed by the Problem using the new transformation
        matrix), and (3) all other child objects (WingCrossSections, Airfoils, vortices)
        have no dependency on Cg_GP1_CgP1.

        :param new_Cg_GP1_CgP1: An array-like object of 3 numbers representing the
            position of the new Airplane's CG (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). Can be a list, tuple, or ndarray.
            Values are converted to floats internally. The units are in meters.
        :return: A new Airplane with the specified position and deep copied geometry.
        """
        # Validate the new position.
        validated_Cg_GP1_CgP1 = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                new_Cg_GP1_CgP1, "new_Cg_GP1_CgP1"
            )
        )
        validated_Cg_GP1_CgP1.flags.writeable = False

        # Create a new Airplane instance without calling __init__ to avoid redundant
        # validation and Wing symmetry processing.
        new_airplane = object.__new__(Airplane)

        # Deep copy the Wings into a new tuple.
        memo: dict = {id(self): new_airplane}
        new_airplane._wings = tuple(copy.deepcopy(wing, memo) for wing in self._wings)

        # Copy immutable attributes, using the new Cg_GP1_CgP1.
        new_airplane._name = self._name
        new_airplane._Cg_GP1_CgP1 = validated_Cg_GP1_CgP1
        new_airplane._weight = self._weight
        new_airplane._s_ref = self._s_ref
        new_airplane._c_ref = self._c_ref
        new_airplane._b_ref = self._b_ref

        # Copy _num_panels cache (depends only on Wings, not position).
        # Reset _T_pas_G_Cg_to_GP1_CgP1 to None (depends on Cg_GP1_CgP1).
        new_airplane._num_panels = self._num_panels
        new_airplane._T_pas_G_Cg_to_GP1_CgP1 = None

        # Reset loads and load coefficients to None (solver will compute these).
        new_airplane.forces_W = None
        new_airplane.forceCoefficients_W = None
        new_airplane.moments_W_CgP1 = None
        new_airplane.momentCoefficients_W_CgP1 = None

        return new_airplane

    # --- Immutable: read only properties ---
    @property
    def wings(self) -> tuple[wing_mod.Wing, ...]:
        return self._wings

    @property
    def name(self) -> str:
        return self._name

    @property
    def Cg_GP1_CgP1(self) -> np.ndarray:
        return self._Cg_GP1_CgP1

    @property
    def weight(self) -> float:
        return self._weight

    @property
    def s_ref(self) -> float:
        return self._s_ref

    @property
    def c_ref(self) -> float:
        return self._c_ref

    @property
    def b_ref(self) -> float:
        return self._b_ref

    # --- Immutable derived: manual lazy caching ---
    @property
    def num_panels(self) -> int:
        """The total number of Panels across all Wings.

        :return: The total number of Panels.
        """
        if self._num_panels is None:
            self._num_panels = sum(
                wing.num_panels if wing.num_panels is not None else 0
                for wing in self._wings
            )
        return self._num_panels

    @property
    def T_pas_G_Cg_to_GP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix from this Airplane's geometry axes,
        relative to this Airplane's CG to the first Airplane's geometry axes, relative
        to the first Airplane's CG.

        Computes the transformation chain: G_Cg > GP1_CgP1. This transformation matrix
        is used to position Airplanes relative to one another, in problems with more
        than one Airplane. If this Airplane is the first Airplane (where Cg_GP1_CgP1 =
        [0, 0, 0]), it returns an identity transformation.

        :return: A (4,4) ndarray of floats representing the passive transformation
            matrix from this Airplane's geometry axes, relative to its CG to the first
            Airplane's geometry axes, relative to its CG.
        """
        if self._T_pas_G_Cg_to_GP1_CgP1 is None:
            # generate_trans_T with passive=True expects the `translations` parameter
            # to be the position of the target reference point (CgP1) relative to the
            # source reference point (Cg). Using the notation from
            # AXES_POINTS_AND_FRAMES.md: translations=CgP1_G_Cg. However, we have
            # Cg_GP1_CgP1 (position of Cg, in GP1 axes, relative to CgP1). Since
            # geometry axes G and GP1 are parallel (pure translation, no rotation):
            # CgP1_G_Cg=-Cg_GP1_CgP1.
            self._T_pas_G_Cg_to_GP1_CgP1 = _transformations.generate_trans_T(
                translations=-self._Cg_GP1_CgP1, passive=True
            )
            self._T_pas_G_Cg_to_GP1_CgP1.flags.writeable = False
        return self._T_pas_G_Cg_to_GP1_CgP1

    # --- Other methods ---
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

        # Define visualization constants.
        panel_color = "chartreuse"
        plotter_background_color = "black"
        window_size = [1024, 768]
        quality = 75

        # Initialize the plotter and set it to use parallel projection.
        plotter = pv.Plotter(window_size=window_size, lighting=None)
        plotter.enable_parallel_projection()  # type: ignore[call-arg]

        # Initialize empty arrays to hold the Panels' vertices and faces.
        panel_vertices = np.empty((0, 3), dtype=float)
        panel_faces = np.empty(0, dtype=int)

        # Initialize a variable to keep track of how many Panels' data has been added
        # to the ndarrays.
        panel_num = 0

        # Iterate through this Airplane's Wings.
        for wing in self._wings:
            # Unravel the Wing's Panel matrix and iterate through it.
            if wing.panels is None:
                continue
            panels = np.ravel(wing.panels)
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

        # Add the Panels to the plotter
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
            smooth_shading=False,
        )

        # Set the plotter's background color
        plotter.set_background(color=plotter_background_color)  # type: ignore[call-arg]

        if not testing:
            # Show the plotter so the user can adjust the camera position and window
            plotter.show(
                title=f"Airplane: {self._name}",
                cpos=(-1, -1, 1),
                full_screen=False,
                auto_close=False,
            )
        else:
            # Show the plotter for 1 second, then proceed automatically (for testing)
            plotter.show(
                title=f"Airplane: {self._name}",
                cpos=(-1, -1, 1),
                full_screen=False,
                interactive=False,
                auto_close=False,
            )
            time.sleep(1)

        # If the user wants to save the image, take a screenshot and save as WebP
        if save:
            screenshot = plotter.screenshot(
                filename=None,
                transparent_background=True,
                return_img=True,
            )
            image = webp.Image.fromarray(
                cast(np.ndarray[Any, Any], screenshot),
            )
            webp.save_image(
                img=image,
                file_path=f"{self._name}_geometry.webp",
                lossless=False,
                quality=quality,
            )

        # Close all the plotters
        pv.close_all()

    def get_plottable_data(
        self, show: bool | np.bool_ = False
    ) -> list[list[list[np.ndarray]]] | None:
        """Returns plottable data for this Airplane's Airfoils' outlines and mean camber
        lines.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: If show is True, returns None. If show is False, returns a list of sub
            lists (one sub list for each of this Airplane's Wings). Each sub list
            contains sub sub lists (one for each of this Wing's WingCrossSections). Each
            sub sub list contains two ndarrays. The first ndarray contains points on
            that WingCrossSection's Airfoil's outline and the second contains points on
            its mean camber line. The points are in geometry axes, relative to the CG.
            The units are in meters.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        airfoilOutlines_G_Cg = []
        airfoilMcls_G_Cg = []
        for wing_id, wing in enumerate(self._wings):
            plottable_data = wing.get_plottable_data(show=False)

            assert plottable_data is not None
            [airfoilOutlines_Wn_ler, airfoilMcls_Wn_ler] = plottable_data

            these_airfoilOutlines_G_Cg = []
            these_airfoilMcls_G_Cg = []
            for airfoil_id in range(len(airfoilOutlines_Wn_ler)):
                airfoilOutline_Wn_ler = airfoilOutlines_Wn_ler[airfoil_id]
                airfoilMcl_Wn_ler = airfoilMcls_Wn_ler[airfoil_id]

                assert wing.T_pas_Wn_Ler_to_G_Cg is not None
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

        plotter.add_actor(AxesGCg)  # type: ignore[arg-type]

        for wing_id, wing in enumerate(self._wings):
            wing_num = wing_id + 1

            assert wing.T_pas_G_Cg_to_Wn_Ler is not None
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

            plotter.add_actor(AxesWLerWcs1Lp1_G_Cg)  # type: ignore[arg-type]

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
                        x_label=f"Wcs{wing_cross_section_num}Wn{wing_num}X@Lp"
                        f"{wing_cross_section_num}Wn{wing_num}",
                        y_label=f"Wcs{wing_cross_section_num}Wn{wing_num}Y@Lp"
                        f"{wing_cross_section_num}Wn{wing_num}",
                        z_label=f"Wcs{wing_cross_section_num}Wn{wing_num}Z@Lp"
                        f"{wing_cross_section_num}Wn{wing_num}",
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

                    plotter.add_actor(AxesWcsLp_G_Cg)  # type: ignore[arg-type]

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
                    # vertices and faces.
                    panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                    panel_faces = np.hstack((panel_faces, panel_face_to_add))

                    # Update the number of previous Panels.
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

        plotter.enable_parallel_projection()  # type: ignore[call-arg]

        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )

        return None

    def validate_first_airplane_constraints(self) -> None:
        """Validates that the first Airplane in a simulation has Cg_GP1_CgP1 set to
        zeros.

        This method should be called by SteadyProblem or UnsteadyProblem.

        :return: None
        """
        if not np.allclose(self._Cg_GP1_CgP1, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The first Airplane in a simulation must have Cg_GP1_CgP1 set to ("
                "0.0, 0.0, 0.0) by definition."
            )

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
        # Determine if the symmetry plane is coincident with the wing axes' xz plane.
        # If symmetryNormal_G or symmetryPoint_G_Cg is None, then there is no
        # symmetry and the symmetry plane doesn't exist. Otherwise, the symmetry
        # plane is coincident to the wing axes' xz plane if Ler_Gs_Cgs lies on the
        # symmetry plane, and if symmetryNormal_G is parallel with WnY_G. We don't
        # need to check types, values, or normalize because this is done in Wing's
        # init method.
        coincident_symmetry_plane = True
        if wing.symmetryPoint_G_Cg is None or wing.symmetryNormal_G is None:
            coincident_symmetry_plane = False
        else:
            # If the symmetry plane exists, we first need to check if its normal
            # vector is parallel with the wing axes' y axis vector.

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
                # If the symmetry plane's normal vector and the wing axes y axis
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

                reflected_airfoil = copy.deepcopy(airfoil)

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

            assert wing.symmetryNormal_G is not None
            assert wing.symmetryPoint_G_Cg is not None
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
