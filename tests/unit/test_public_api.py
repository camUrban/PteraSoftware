"""This module contains tests that pin Ptera Software's public API.

The public API is the names in the package's __all__, the public members of each public
class, and the parameter names and kinds of each public callable. Each is pinned below
as a deliberate second record, so a change to the public API fails these tests until the
matching pin is updated. Deprecated members and parameters are not part of the public
API, so they are left out of the pins, as found by the same source scan that keeps them
out of the API reference.
"""

import ast
import builtins
import inspect
import sys
import textwrap
import types
import unittest
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pterasoftware as ps
from pterasoftware import _deprecation

DEPRECATED_MEMBERS, DEPRECATED_PARAMETERS = _deprecation.find_deprecations(
    Path(ps.__file__).parent
)

# This is the pinned copy of the package's __all__, sorted.
PUBLIC_NAMES: tuple[str, ...] = (
    "AeroelasticAirplaneMovement",
    "AeroelasticMovement",
    "AeroelasticUnsteadyProblem",
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver",
    "AeroelasticWingCrossSectionMovement",
    "AeroelasticWingMovement",
    "Airfoil",
    "Airplane",
    "AirplaneMovement",
    "FreeFlightMovement",
    "FreeFlightOperatingPointMovement",
    "FreeFlightUnsteadyProblem",
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver",
    "Movement",
    "OperatingPoint",
    "OperatingPointMovement",
    "SteadyHorseshoeVortexLatticeMethodSolver",
    "SteadyProblem",
    "SteadyRingVortexLatticeMethodSolver",
    "UnsteadyProblem",
    "UnsteadyRingVortexLatticeMethodSolver",
    "Wing",
    "WingCrossSection",
    "WingCrossSectionMovement",
    "WingMovement",
    "analyze_steady_convergence",
    "analyze_steady_trim",
    "analyze_unsteady_convergence",
    "analyze_unsteady_trim",
    "animate",
    "draw",
    "load",
    "log_results",
    "plot_results_versus_time",
    "save",
    "set_up_logging",
)

# These are the pinned public members of each public class, sorted. They include the
# members inherited from internal parents.
PUBLIC_MEMBERS: dict[str, tuple[str, ...]] = {
    "AeroelasticAirplaneMovement": (
        "all_periods",
        "ampCg_GP1_CgP1",
        "base_airplane",
        "max_period",
        "periodCg_GP1_CgP1",
        "phaseCg_GP1_CgP1",
        "spacingCg_GP1_CgP1",
        "wing_movements",
    ),
    "AeroelasticMovement": (
        "airplane_movements",
        "delta_time",
        "lcm_period",
        "max_period",
        "max_wake_rows",
        "min_period",
        "num_steps",
        "operating_point_movement",
        "operating_points",
        "static",
    ),
    "AeroelasticUnsteadyProblem": (
        "damping_constant_rad",
        "delta_time",
        "diagram",
        "finalCrosswindForceCoefficients_W",
        "finalCrosswindForces_W",
        "finalForceCoefficients_G",
        "finalForceCoefficients_W",
        "finalForces_G",
        "finalForces_W",
        "finalInducedDragCoefficients_W",
        "finalInducedDrags_W",
        "finalLiftCoefficients_W",
        "finalLifts_W",
        "finalMeanCrosswindForceCoefficients_W",
        "finalMeanCrosswindForces_W",
        "finalMeanForceCoefficients_G",
        "finalMeanForceCoefficients_W",
        "finalMeanForces_G",
        "finalMeanForces_W",
        "finalMeanInducedDragCoefficients_W",
        "finalMeanInducedDrags_W",
        "finalMeanLiftCoefficients_W",
        "finalMeanLifts_W",
        "finalMeanMomentCoefficients_G_Cg",
        "finalMeanMomentCoefficients_W_Cg",
        "finalMeanMomentCoefficients_W_CgP1",
        "finalMeanMoments_G_Cg",
        "finalMeanMoments_W_Cg",
        "finalMeanMoments_W_CgP1",
        "finalMeanPitchingMomentCoefficients_W_Cg",
        "finalMeanPitchingMoments_W_Cg",
        "finalMeanRollingMomentCoefficients_W_Cg",
        "finalMeanRollingMoments_W_Cg",
        "finalMeanYawingMomentCoefficients_W_Cg",
        "finalMeanYawingMoments_W_Cg",
        "finalMomentCoefficients_G_Cg",
        "finalMomentCoefficients_W_Cg",
        "finalMomentCoefficients_W_CgP1",
        "finalMoments_G_Cg",
        "finalMoments_W_Cg",
        "finalMoments_W_CgP1",
        "finalPitchingMomentCoefficients_W_Cg",
        "finalPitchingMoments_W_Cg",
        "finalRmsForceCoefficients_G",
        "finalRmsForceCoefficients_W",
        "finalRmsForces_G",
        "finalRmsForces_W",
        "finalRmsMomentCoefficients_G_Cg",
        "finalRmsMomentCoefficients_W_Cg",
        "finalRmsMomentCoefficients_W_CgP1",
        "finalRmsMoments_G_Cg",
        "finalRmsMoments_W_Cg",
        "finalRmsMoments_W_CgP1",
        "finalRollingMomentCoefficients_W_Cg",
        "finalRollingMoments_W_Cg",
        "finalYawingMomentCoefficients_W_Cg",
        "finalYawingMoments_W_Cg",
        "first_averaging_step",
        "first_results_step",
        "listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz",
        "max_wake_rows",
        "movement",
        "num_steps",
        "only_final_results",
        "spring_constant_rad",
        "steady_problems",
        "step_discards",
        "wing_density",
    ),
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver": (
        "current_airplanes",
        "current_operating_point",
        "delta_time",
        "diagram",
        "first_results_step",
        "force_method",
        "gridStreamlinePoints_GP1_CgP1",
        "listStackBlwrvp_GP1_CgP1",
        "listStackBrwrvp_GP1_CgP1",
        "listStackFlwrvp_GP1_CgP1",
        "listStackFrwrvp_GP1_CgP1",
        "list_num_wake_vortices",
        "moments_GP1_Slep",
        "num_airplanes",
        "num_panels",
        "num_steps",
        "panel_areas",
        "panel_is_leading_edge",
        "panel_is_left_edge",
        "panel_is_right_edge",
        "panel_is_trailing_edge",
        "panels",
        "prescribed_wake",
        "ran",
        "run",
        "stackBbrv_GP1",
        "stackBlbrvp_GP1_CgP1",
        "stackBrbrvp_GP1_CgP1",
        "stackCblvpb_GP1_CgP1",
        "stackCblvpf_GP1_CgP1",
        "stackCblvpl_GP1_CgP1",
        "stackCblvpr_GP1_CgP1",
        "stackCpp_GP1_CgP1",
        "stackFbrv_GP1",
        "stackFlbrvp_GP1_CgP1",
        "stackFrbrvp_GP1_CgP1",
        "stackLbrv_GP1",
        "stackRbrv_GP1",
        "stackSeedPoints_GP1_CgP1",
        "stackUnitNormals_GP1",
        "steady_problems",
        "unsteady_problem",
    ),
    "AeroelasticWingCrossSectionMovement": (
        "all_periods",
        "ampAngles_Wcsp_to_Wcs_ixyz",
        "ampLp_Wcsp_Lpp",
        "base_wing_cross_section",
        "max_period",
        "periodAngles_Wcsp_to_Wcs_ixyz",
        "periodLp_Wcsp_Lpp",
        "phaseAngles_Wcsp_to_Wcs_ixyz",
        "phaseLp_Wcsp_Lpp",
        "spacingAngles_Wcsp_to_Wcs_ixyz",
        "spacingLp_Wcsp_Lpp",
    ),
    "AeroelasticWingMovement": (
        "all_periods",
        "ampAngles_Gs_to_Wn_ixyz",
        "ampLer_Gs_Cgs",
        "base_wing",
        "max_period",
        "periodAngles_Gs_to_Wn_ixyz",
        "periodLer_Gs_Cgs",
        "phaseAngles_Gs_to_Wn_ixyz",
        "phaseLer_Gs_Cgs",
        "rotationPointOffset_Gs_Ler",
        "spacingAnglesSecondDerivative_Gs_to_Wn_ixyz",
        "spacingAngles_Gs_to_Wn_ixyz",
        "spacingLer_Gs_Cgs",
        "wing_cross_section_movements",
    ),
    "Airfoil": (
        "diagram",
        "mcl_A_Lp",
        "n_points_per_side",
        "name",
        "outline_A_Lp",
        "resample",
    ),
    "Airplane": (
        "Cg_GP1_CgP1",
        "T_pas_G_Cg_to_GP1_CgP1",
        "b_ref",
        "c_ref",
        "crosswindForceCoefficient_W",
        "crosswindForce_W",
        "deep_copy_with_Cg_GP1_CgP1",
        "diagram",
        "forceCoefficients_G",
        "forceCoefficients_W",
        "forces_G",
        "forces_W",
        "inducedDragCoefficient_W",
        "inducedDrag_W",
        "liftCoefficient_W",
        "lift_W",
        "momentCoefficients_G_Cg",
        "momentCoefficients_W_Cg",
        "momentCoefficients_W_CgP1",
        "moments_G_Cg",
        "moments_W_Cg",
        "moments_W_CgP1",
        "name",
        "num_panels",
        "pitchingMomentCoefficient_W_Cg",
        "pitchingMoment_W_Cg",
        "process_wing_symmetry",
        "rollingMomentCoefficient_W_Cg",
        "rollingMoment_W_Cg",
        "s_ref",
        "validate_first_airplane_constraints",
        "weight",
        "wings",
        "yawingMomentCoefficient_W_Cg",
        "yawingMoment_W_Cg",
    ),
    "AirplaneMovement": (
        "all_periods",
        "ampCg_GP1_CgP1",
        "base_airplane",
        "max_period",
        "periodCg_GP1_CgP1",
        "phaseCg_GP1_CgP1",
        "spacingCg_GP1_CgP1",
        "wing_movements",
    ),
    "FreeFlightMovement": (
        "airplane_movements",
        "airplanes",
        "delta_time",
        "free_num_steps",
        "lcm_period",
        "max_period",
        "max_wake_rows",
        "min_period",
        "num_steps",
        "operating_point_movement",
        "prescribed_num_steps",
        "static",
    ),
    "FreeFlightOperatingPointMovement": (
        "ampVCg__E",
        "base_operating_point",
        "max_period",
        "operating_points",
        "periodVCg__E",
        "phaseVCg__E",
        "spacingVCg__E",
    ),
    "FreeFlightUnsteadyProblem": (
        "I_BP1_CgP1",
        "delta_time",
        "diagram",
        "external_loads_fn",
        "finalCrosswindForceCoefficients_W",
        "finalCrosswindForces_W",
        "finalForceCoefficients_G",
        "finalForceCoefficients_W",
        "finalForces_G",
        "finalForces_W",
        "finalInducedDragCoefficients_W",
        "finalInducedDrags_W",
        "finalLiftCoefficients_W",
        "finalLifts_W",
        "finalMeanCrosswindForceCoefficients_W",
        "finalMeanCrosswindForces_W",
        "finalMeanForceCoefficients_G",
        "finalMeanForceCoefficients_W",
        "finalMeanForces_G",
        "finalMeanForces_W",
        "finalMeanInducedDragCoefficients_W",
        "finalMeanInducedDrags_W",
        "finalMeanLiftCoefficients_W",
        "finalMeanLifts_W",
        "finalMeanMomentCoefficients_G_Cg",
        "finalMeanMomentCoefficients_W_Cg",
        "finalMeanMomentCoefficients_W_CgP1",
        "finalMeanMoments_G_Cg",
        "finalMeanMoments_W_Cg",
        "finalMeanMoments_W_CgP1",
        "finalMeanPitchingMomentCoefficients_W_Cg",
        "finalMeanPitchingMoments_W_Cg",
        "finalMeanRollingMomentCoefficients_W_Cg",
        "finalMeanRollingMoments_W_Cg",
        "finalMeanYawingMomentCoefficients_W_Cg",
        "finalMeanYawingMoments_W_Cg",
        "finalMomentCoefficients_G_Cg",
        "finalMomentCoefficients_W_Cg",
        "finalMomentCoefficients_W_CgP1",
        "finalMoments_G_Cg",
        "finalMoments_W_Cg",
        "finalMoments_W_CgP1",
        "finalPitchingMomentCoefficients_W_Cg",
        "finalPitchingMoments_W_Cg",
        "finalRmsForceCoefficients_G",
        "finalRmsForceCoefficients_W",
        "finalRmsForces_G",
        "finalRmsForces_W",
        "finalRmsMomentCoefficients_G_Cg",
        "finalRmsMomentCoefficients_W_Cg",
        "finalRmsMomentCoefficients_W_CgP1",
        "finalRmsMoments_G_Cg",
        "finalRmsMoments_W_Cg",
        "finalRmsMoments_W_CgP1",
        "finalRollingMomentCoefficients_W_Cg",
        "finalRollingMoments_W_Cg",
        "finalYawingMomentCoefficients_W_Cg",
        "finalYawingMoments_W_Cg",
        "first_averaging_step",
        "first_results_step",
        "k_max",
        "mass",
        "max_wake_rows",
        "movement",
        "num_steps",
        "only_final_results",
        "steady_problems",
    ),
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver": (
        "current_airplanes",
        "current_operating_point",
        "delta_time",
        "diagram",
        "first_results_step",
        "force_method",
        "gridStreamlinePoints_GP1_CgP1",
        "listStackBlwrvp_GP1_CgP1",
        "listStackBrwrvp_GP1_CgP1",
        "listStackFlwrvp_GP1_CgP1",
        "listStackFrwrvp_GP1_CgP1",
        "list_num_wake_vortices",
        "num_airplanes",
        "num_panels",
        "num_steps",
        "panel_areas",
        "panel_is_leading_edge",
        "panel_is_left_edge",
        "panel_is_right_edge",
        "panel_is_trailing_edge",
        "panels",
        "prescribed_wake",
        "ran",
        "run",
        "stackBbrv_GP1",
        "stackBlbrvp_GP1_CgP1",
        "stackBrbrvp_GP1_CgP1",
        "stackCblvpb_GP1_CgP1",
        "stackCblvpf_GP1_CgP1",
        "stackCblvpl_GP1_CgP1",
        "stackCblvpr_GP1_CgP1",
        "stackCpp_GP1_CgP1",
        "stackFbrv_GP1",
        "stackFlbrvp_GP1_CgP1",
        "stackFrbrvp_GP1_CgP1",
        "stackLbrv_GP1",
        "stackRbrv_GP1",
        "stackSeedPoints_GP1_CgP1",
        "stackUnitNormals_GP1",
        "steady_problems",
        "unsteady_problem",
    ),
    "Movement": (
        "airplane_movements",
        "airplanes",
        "delta_time",
        "lcm_period",
        "max_period",
        "max_wake_chords",
        "max_wake_cycles",
        "max_wake_rows",
        "min_period",
        "num_chords",
        "num_cycles",
        "num_steps",
        "operating_point_movement",
        "operating_points",
        "static",
    ),
    "OperatingPoint": (
        "CgP1_E_Eo",
        "T_pas_BP1_CgP1_to_E_CgP1",
        "T_pas_BP1_CgP1_to_GP1_CgP1",
        "T_pas_BP1_CgP1_to_W_CgP1",
        "T_pas_E_CgP1_to_BP1_CgP1",
        "T_pas_E_CgP1_to_GP1_CgP1",
        "T_pas_E_CgP1_to_W_CgP1",
        "T_pas_GP1_CgP1_to_BP1_CgP1",
        "T_pas_GP1_CgP1_to_E_CgP1",
        "T_pas_GP1_CgP1_to_W_CgP1",
        "T_pas_W_CgP1_to_BP1_CgP1",
        "T_pas_W_CgP1_to_E_CgP1",
        "T_pas_W_CgP1_to_GP1_CgP1",
        "alpha",
        "angles_E_to_BP1_izyx",
        "beta",
        "externalFX_W",
        "g_E",
        "nu",
        "omegas_BP1__E",
        "qInf__E",
        "rho",
        "surfaceNormal_E",
        "surfaceNormal_GP1",
        "surfacePoint_E_Eo",
        "surfacePoint_GP1_CgP1",
        "surfaceReflect_T_act_GP1_CgP1",
        "vCg__E",
        "vInfHat_GP1__E",
        "vInf_GP1__E",
    ),
    "OperatingPointMovement": (
        "ampVCg__E",
        "base_operating_point",
        "max_period",
        "periodVCg__E",
        "phaseVCg__E",
        "spacingVCg__E",
    ),
    "SteadyHorseshoeVortexLatticeMethodSolver": (
        "airplanes",
        "diagram",
        "gridStreamlinePoints_GP1_CgP1",
        "num_airplanes",
        "num_panels",
        "operating_point",
        "panels",
        "ran",
        "reynolds_numbers",
        "run",
        "stackFreestreamWingInfluences__E",
        "stackSeedPoints_GP1_CgP1",
        "stackUnitNormals_GP1",
        "vInf_GP1__E",
    ),
    "SteadyProblem": (
        "airplanes",
        "diagram",
        "operating_point",
        "reynolds_numbers",
    ),
    "SteadyRingVortexLatticeMethodSolver": (
        "airplanes",
        "diagram",
        "gridStreamlinePoints_GP1_CgP1",
        "num_airplanes",
        "num_panels",
        "num_trailing_edge_panels",
        "operating_point",
        "panel_areas",
        "panel_is_leading_edge",
        "panel_is_left_edge",
        "panel_is_right_edge",
        "panel_is_trailing_edge",
        "panels",
        "ran",
        "reynolds_numbers",
        "run",
        "stackBbrv_GP1",
        "stackBlbrvp_GP1_CgP1",
        "stackBrbrvp_GP1_CgP1",
        "stackCblvpb_GP1_CgP1",
        "stackCblvpf_GP1_CgP1",
        "stackCblvpl_GP1_CgP1",
        "stackCblvpr_GP1_CgP1",
        "stackCpp_GP1_CgP1",
        "stackFbrv_GP1",
        "stackFlbrvp_GP1_CgP1",
        "stackFrbrvp_GP1_CgP1",
        "stackFreestreamWingInfluences__E",
        "stackLbrv_GP1",
        "stackRbrv_GP1",
        "stackSeedPoints_GP1_CgP1",
        "stackUnitNormals_GP1",
        "vInf_GP1__E",
    ),
    "UnsteadyProblem": (
        "delta_time",
        "diagram",
        "finalCrosswindForceCoefficients_W",
        "finalCrosswindForces_W",
        "finalForceCoefficients_G",
        "finalForceCoefficients_W",
        "finalForces_G",
        "finalForces_W",
        "finalInducedDragCoefficients_W",
        "finalInducedDrags_W",
        "finalLiftCoefficients_W",
        "finalLifts_W",
        "finalMeanCrosswindForceCoefficients_W",
        "finalMeanCrosswindForces_W",
        "finalMeanForceCoefficients_G",
        "finalMeanForceCoefficients_W",
        "finalMeanForces_G",
        "finalMeanForces_W",
        "finalMeanInducedDragCoefficients_W",
        "finalMeanInducedDrags_W",
        "finalMeanLiftCoefficients_W",
        "finalMeanLifts_W",
        "finalMeanMomentCoefficients_G_Cg",
        "finalMeanMomentCoefficients_W_Cg",
        "finalMeanMomentCoefficients_W_CgP1",
        "finalMeanMoments_G_Cg",
        "finalMeanMoments_W_Cg",
        "finalMeanMoments_W_CgP1",
        "finalMeanPitchingMomentCoefficients_W_Cg",
        "finalMeanPitchingMoments_W_Cg",
        "finalMeanRollingMomentCoefficients_W_Cg",
        "finalMeanRollingMoments_W_Cg",
        "finalMeanYawingMomentCoefficients_W_Cg",
        "finalMeanYawingMoments_W_Cg",
        "finalMomentCoefficients_G_Cg",
        "finalMomentCoefficients_W_Cg",
        "finalMomentCoefficients_W_CgP1",
        "finalMoments_G_Cg",
        "finalMoments_W_Cg",
        "finalMoments_W_CgP1",
        "finalPitchingMomentCoefficients_W_Cg",
        "finalPitchingMoments_W_Cg",
        "finalRmsForceCoefficients_G",
        "finalRmsForceCoefficients_W",
        "finalRmsForces_G",
        "finalRmsForces_W",
        "finalRmsMomentCoefficients_G_Cg",
        "finalRmsMomentCoefficients_W_Cg",
        "finalRmsMomentCoefficients_W_CgP1",
        "finalRmsMoments_G_Cg",
        "finalRmsMoments_W_Cg",
        "finalRmsMoments_W_CgP1",
        "finalRollingMomentCoefficients_W_Cg",
        "finalRollingMoments_W_Cg",
        "finalYawingMomentCoefficients_W_Cg",
        "finalYawingMoments_W_Cg",
        "first_averaging_step",
        "first_results_step",
        "max_wake_rows",
        "movement",
        "num_steps",
        "only_final_results",
        "steady_problems",
    ),
    "UnsteadyRingVortexLatticeMethodSolver": (
        "current_airplanes",
        "current_operating_point",
        "delta_time",
        "diagram",
        "first_results_step",
        "force_method",
        "gridStreamlinePoints_GP1_CgP1",
        "listStackBlwrvp_GP1_CgP1",
        "listStackBrwrvp_GP1_CgP1",
        "listStackFlwrvp_GP1_CgP1",
        "listStackFrwrvp_GP1_CgP1",
        "list_num_wake_vortices",
        "num_airplanes",
        "num_panels",
        "num_steps",
        "panel_areas",
        "panel_is_leading_edge",
        "panel_is_left_edge",
        "panel_is_right_edge",
        "panel_is_trailing_edge",
        "panels",
        "prescribed_wake",
        "ran",
        "run",
        "stackBbrv_GP1",
        "stackBlbrvp_GP1_CgP1",
        "stackBrbrvp_GP1_CgP1",
        "stackCblvpb_GP1_CgP1",
        "stackCblvpf_GP1_CgP1",
        "stackCblvpl_GP1_CgP1",
        "stackCblvpr_GP1_CgP1",
        "stackCpp_GP1_CgP1",
        "stackFbrv_GP1",
        "stackFlbrvp_GP1_CgP1",
        "stackFrbrvp_GP1_CgP1",
        "stackLbrv_GP1",
        "stackRbrv_GP1",
        "stackSeedPoints_GP1_CgP1",
        "stackUnitNormals_GP1",
        "steady_problems",
        "unsteady_problem",
    ),
    "Wing": (
        "Ler_Gs_Cgs",
        "T_pas_G_Cg_to_Wn_Ler",
        "T_pas_Wn_Ler_to_G_Cg",
        "WnX_G",
        "WnY_G",
        "WnZ_G",
        "angles_Gs_to_Wn_ixyz",
        "average_panel_aspect_ratio",
        "children_T_pas_G_Cg_to_Wcs_Lp",
        "children_T_pas_Wcs_Lp_to_G_Cg",
        "children_T_pas_Wcs_Lp_to_Wn_Ler",
        "children_T_pas_Wn_Ler_to_Wcs_Lp",
        "chordwise_spacing",
        "diagram",
        "from_edge_points",
        "generate_mesh",
        "gridWrvp_GP1_CgP1",
        "leadingEdgePoints_Wn_Ler",
        "mean_aerodynamic_chord",
        "mirror_only",
        "name",
        "num_chordwise_panels",
        "num_panels",
        "num_spanwise_panels",
        "panels",
        "projected_area",
        "span",
        "spanwise_mesh",
        "standard_mean_chord",
        "symmetric",
        "symmetryNormal_G",
        "symmetryPoint_G_Cg",
        "symmetry_type",
        "tip_trim_fraction",
        "trailingEdgePoints_Wn_Ler",
        "wing_cross_sections",
    ),
    "WingCrossSection": (
        "Lp_Wcsp_Lpp",
        "T_pas_Wcs_Lp_to_Wcsp_Lpp",
        "T_pas_Wcsp_Lpp_to_Wcs_Lp",
        "airfoil",
        "angles_Wcsp_to_Wcs_ixyz",
        "chord",
        "control_surface_deflection",
        "control_surface_hinge_point",
        "control_surface_symmetry_type",
        "diagram",
        "num_spanwise_panels",
        "spanwise_spacing",
        "symmetry_type",
    ),
    "WingCrossSectionMovement": (
        "all_periods",
        "ampAngles_Wcsp_to_Wcs_ixyz",
        "ampLp_Wcsp_Lpp",
        "base_wing_cross_section",
        "max_period",
        "periodAngles_Wcsp_to_Wcs_ixyz",
        "periodLp_Wcsp_Lpp",
        "phaseAngles_Wcsp_to_Wcs_ixyz",
        "phaseLp_Wcsp_Lpp",
        "spacingAngles_Wcsp_to_Wcs_ixyz",
        "spacingLp_Wcsp_Lpp",
    ),
    "WingMovement": (
        "all_periods",
        "ampAngles_Gs_to_Wn_ixyz",
        "ampLer_Gs_Cgs",
        "base_wing",
        "max_period",
        "periodAngles_Gs_to_Wn_ixyz",
        "periodLer_Gs_Cgs",
        "phaseAngles_Gs_to_Wn_ixyz",
        "phaseLer_Gs_Cgs",
        "rotationPointOffset_Gs_Ler",
        "spacingAngles_Gs_to_Wn_ixyz",
        "spacingLer_Gs_Cgs",
        "wing_cross_section_movements",
    ),
}

# These are the pinned parameters of each public callable, keyed by the public
# function's name, or by the public class's name followed by its method's name. They
# include the methods inherited from internal parents. Each parameter is listed by name
# in order, with a "/" after the positional-only parameters and a "*" before the
# keyword-only ones, as in a Python signature. A method's self and a class method's cls
# are left out.
PUBLIC_PARAMETERS: dict[str, tuple[str, ...]] = {
    "AeroelasticAirplaneMovement.__init__": (
        "*",
        "base_airplane",
        "wing_movements",
        "ampCg_GP1_CgP1",
        "periodCg_GP1_CgP1",
        "spacingCg_GP1_CgP1",
        "phaseCg_GP1_CgP1",
    ),
    "AeroelasticMovement.__init__": (
        "*",
        "airplane_movements",
        "operating_point_movement",
        "delta_time",
        "num_steps",
        "max_wake_rows",
    ),
    "AeroelasticUnsteadyProblem.__init__": (
        "*",
        "movement",
        "wing_density",
        "spring_constant_rad",
        "damping_constant_rad",
        "step_discards",
    ),
    "AeroelasticUnsteadyProblem.diagram": (
        "*",
        "step",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver.__init__": (
        "*",
        "aeroelastic_unsteady_problem",
    ),
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver.diagram": (
        "*",
        "step",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "simplify_vortices",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver.run": (
        "*",
        "prescribed_wake",
        "calculate_streamlines",
        "show_progress",
        "force_method",
    ),
    "AeroelasticWingCrossSectionMovement.__init__": (
        "*",
        "base_wing_cross_section",
        "ampLp_Wcsp_Lpp",
        "periodLp_Wcsp_Lpp",
        "spacingLp_Wcsp_Lpp",
        "phaseLp_Wcsp_Lpp",
        "ampAngles_Wcsp_to_Wcs_ixyz",
        "periodAngles_Wcsp_to_Wcs_ixyz",
        "spacingAngles_Wcsp_to_Wcs_ixyz",
        "phaseAngles_Wcsp_to_Wcs_ixyz",
    ),
    "AeroelasticWingMovement.__init__": (
        "*",
        "base_wing",
        "wing_cross_section_movements",
        "ampLer_Gs_Cgs",
        "periodLer_Gs_Cgs",
        "spacingLer_Gs_Cgs",
        "phaseLer_Gs_Cgs",
        "ampAngles_Gs_to_Wn_ixyz",
        "periodAngles_Gs_to_Wn_ixyz",
        "spacingAngles_Gs_to_Wn_ixyz",
        "phaseAngles_Gs_to_Wn_ixyz",
        "rotationPointOffset_Gs_Ler",
        "spacingAnglesSecondDerivative_Gs_to_Wn_ixyz",
    ),
    "Airfoil.__init__": ("*", "name", "outline_A_Lp", "resample", "n_points_per_side"),
    "Airfoil.diagram": ("*", "show_mcl", "math_labels", "save", "path", "quality"),
    "Airplane.__init__": (
        "*",
        "wings",
        "name",
        "Cg_GP1_CgP1",
        "weight",
        "s_ref",
        "c_ref",
        "b_ref",
    ),
    "Airplane.deep_copy_with_Cg_GP1_CgP1": ("new_Cg_GP1_CgP1",),
    "Airplane.diagram": (
        "*",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "Airplane.process_wing_symmetry": ("wing",),
    "Airplane.validate_first_airplane_constraints": (),
    "AirplaneMovement.__init__": (
        "*",
        "base_airplane",
        "wing_movements",
        "ampCg_GP1_CgP1",
        "periodCg_GP1_CgP1",
        "spacingCg_GP1_CgP1",
        "phaseCg_GP1_CgP1",
    ),
    "FreeFlightMovement.__init__": (
        "*",
        "airplane_movements",
        "operating_point_movement",
        "delta_time",
        "prescribed_num_steps",
        "free_num_steps",
        "max_wake_rows",
    ),
    "FreeFlightOperatingPointMovement.__init__": ("*", "base_operating_point"),
    "FreeFlightUnsteadyProblem.__init__": (
        "*",
        "movement",
        "mass",
        "I_BP1_CgP1",
        "external_loads_fn",
        "k_max",
        "integrator",
        "extra_xml",
        "mujoco_assets",
    ),
    "FreeFlightUnsteadyProblem.diagram": (
        "*",
        "step",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver.__init__": (
        "*",
        "free_flight_unsteady_problem",
    ),
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver.diagram": (
        "*",
        "step",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "simplify_vortices",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver.run": (
        "*",
        "prescribed_wake",
        "calculate_streamlines",
        "show_progress",
        "force_method",
    ),
    "Movement.__init__": (
        "*",
        "airplane_movements",
        "operating_point_movement",
        "delta_time",
        "num_cycles",
        "num_chords",
        "num_steps",
        "max_wake_rows",
        "max_wake_chords",
        "max_wake_cycles",
    ),
    "OperatingPoint.__init__": (
        "*",
        "rho",
        "vCg__E",
        "alpha",
        "beta",
        "angles_E_to_BP1_izyx",
        "CgP1_E_Eo",
        "surfaceNormal_E",
        "surfacePoint_E_Eo",
        "externalFX_W",
        "nu",
        "g_E",
        "omegas_BP1__E",
    ),
    "OperatingPointMovement.__init__": (
        "*",
        "base_operating_point",
        "ampVCg__E",
        "periodVCg__E",
        "spacingVCg__E",
        "phaseVCg__E",
    ),
    "SteadyHorseshoeVortexLatticeMethodSolver.__init__": ("*", "steady_problem"),
    "SteadyHorseshoeVortexLatticeMethodSolver.diagram": (
        "*",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "simplify_vortices",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "SteadyHorseshoeVortexLatticeMethodSolver.run": ("*", "calculate_streamlines"),
    "SteadyProblem.__init__": ("*", "airplanes", "operating_point"),
    "SteadyProblem.diagram": (
        "*",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "SteadyRingVortexLatticeMethodSolver.__init__": ("*", "steady_problem"),
    "SteadyRingVortexLatticeMethodSolver.diagram": (
        "*",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "simplify_vortices",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "SteadyRingVortexLatticeMethodSolver.run": ("*", "calculate_streamlines"),
    "UnsteadyProblem.__init__": ("*", "movement", "only_final_results"),
    "UnsteadyProblem.diagram": (
        "*",
        "step",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "UnsteadyRingVortexLatticeMethodSolver.__init__": ("*", "unsteady_problem"),
    "UnsteadyRingVortexLatticeMethodSolver.diagram": (
        "*",
        "step",
        "show_airplane_axes_and_points",
        "show_wing_axes_and_points",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "simplify_vortices",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "UnsteadyRingVortexLatticeMethodSolver.run": (
        "*",
        "prescribed_wake",
        "calculate_streamlines",
        "show_progress",
        "force_method",
    ),
    "Wing.__init__": (
        "*",
        "wing_cross_sections",
        "name",
        "Ler_Gs_Cgs",
        "angles_Gs_to_Wn_ixyz",
        "symmetric",
        "mirror_only",
        "symmetryNormal_G",
        "symmetryPoint_G_Cg",
        "explode_into_strips",
        "num_chordwise_panels",
        "chordwise_spacing",
    ),
    "Wing.diagram": (
        "*",
        "show_wing_cross_section_axes_and_points",
        "show_airfoil_axes_and_points",
        "show_airfoils",
        "show_mcls",
        "show_collocation_points",
        "label_collocation_points",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "Wing.from_edge_points": (
        "*",
        "leadingEdgePoints_Wn_Ler",
        "trailingEdgePoints_Wn_Ler",
        "num_wing_cross_sections",
        "airfoil",
        "name",
        "Ler_Gs_Cgs",
        "angles_Gs_to_Wn_ixyz",
        "symmetric",
        "mirror_only",
        "symmetryNormal_G",
        "symmetryPoint_G_Cg",
        "num_chordwise_panels",
        "chordwise_spacing",
        "tip_trim_fraction",
    ),
    "Wing.generate_mesh": ("symmetry_type",),
    "WingCrossSection.__init__": (
        "*",
        "airfoil",
        "num_spanwise_panels",
        "chord",
        "Lp_Wcsp_Lpp",
        "angles_Wcsp_to_Wcs_ixyz",
        "control_surface_symmetry_type",
        "control_surface_hinge_point",
        "control_surface_deflection",
        "spanwise_spacing",
    ),
    "WingCrossSection.diagram": (
        "*",
        "show_airfoil_axes_and_points",
        "show_airfoil",
        "show_mcl",
        "math_labels",
        "save",
        "path",
        "quality",
    ),
    "WingCrossSectionMovement.__init__": (
        "*",
        "base_wing_cross_section",
        "ampLp_Wcsp_Lpp",
        "periodLp_Wcsp_Lpp",
        "spacingLp_Wcsp_Lpp",
        "phaseLp_Wcsp_Lpp",
        "ampAngles_Wcsp_to_Wcs_ixyz",
        "periodAngles_Wcsp_to_Wcs_ixyz",
        "spacingAngles_Wcsp_to_Wcs_ixyz",
        "phaseAngles_Wcsp_to_Wcs_ixyz",
    ),
    "WingMovement.__init__": (
        "*",
        "base_wing",
        "wing_cross_section_movements",
        "ampLer_Gs_Cgs",
        "periodLer_Gs_Cgs",
        "spacingLer_Gs_Cgs",
        "phaseLer_Gs_Cgs",
        "ampAngles_Gs_to_Wn_ixyz",
        "periodAngles_Gs_to_Wn_ixyz",
        "spacingAngles_Gs_to_Wn_ixyz",
        "phaseAngles_Gs_to_Wn_ixyz",
        "rotationPointOffset_Gs_Ler",
    ),
    "analyze_steady_convergence": (
        "ref_problem",
        "/",
        "*",
        "solver_type",
        "panel_aspect_ratio_bounds",
        "num_chordwise_panels_bounds",
        "rtol",
        "atol",
        "coefficient_mask",
        "resolve_converged_solver",
        "cache_path",
    ),
    "analyze_steady_trim": (
        "problem",
        "/",
        "*",
        "solver_type",
        "boundsVCg__E",
        "alpha_bounds",
        "beta_bounds",
        "boundsExternalFX_W",
        "objective_cut_off",
        "num_calls",
    ),
    "analyze_unsteady_convergence": (
        "ref_problem",
        "/",
        "*",
        "prescribed_wake",
        "free_wake",
        "num_cycles_bounds",
        "num_chords_bounds",
        "panel_aspect_ratio_bounds",
        "num_chordwise_panels_bounds",
        "rtol",
        "atol",
        "coefficient_mask",
        "show_solver_progress",
        "resolve_converged_solver",
        "cache_path",
        "force_method",
    ),
    "analyze_unsteady_trim": (
        "problem",
        "/",
        "*",
        "boundsVCg__E",
        "alpha_bounds",
        "beta_bounds",
        "boundsExternalFX_W",
        "objective_cut_off",
        "num_calls",
        "show_solver_progress",
        "force_method",
    ),
    "animate": (
        "unsteady_solver",
        "/",
        "*",
        "scalar_type",
        "show_wake_vortices",
        "show_mujoco_geometry",
        "window_size",
        "save",
        "path",
        "quality",
        "speed",
        "testing",
    ),
    "draw": (
        "solver",
        "/",
        "*",
        "scalar_type",
        "show_streamlines",
        "show_wake_vortices",
        "show_mujoco_geometry",
        "window_size",
        "save",
        "path",
        "quality",
        "testing",
    ),
    "load": ("path", "/", "*", "max_size", "callables"),
    "log_results": ("solver", "/"),
    "plot_results_versus_time": (
        "unsteady_solver",
        "/",
        "*",
        "show",
        "figure_size_in",
        "save",
        "save_csv",
        "directory",
        "prefix",
        "resolution_dpi",
        "file_format",
        "show_titles",
        "font_size",
        "text_color",
        "line_width",
    ),
    "save": ("path", "obj", "/"),
    "set_up_logging": ("*", "level", "handler", "format_string"),
}

# These are the public surfaces allowed to name an internal type. A callable's or
# property's surface is keyed by its label, as in PUBLIC_PARAMETERS, and the parameter
# or "return" whose annotation names the type. A slot attribute's surface is keyed by
# its class's name and the attribute's name. The base unsteady solver's constructor
# accepts the shared parent problem type so the derived solvers can pass their own
# problem types through it, until the core solver split removes the need.
INTERNAL_TYPE_ALLOWLIST = {
    ("UnsteadyRingVortexLatticeMethodSolver.__init__", "unsteady_problem"),
}


def get_public_classes() -> dict[str, type]:
    """Returns the public classes, keyed by their public names.

    :return: A dict mapping each public class's name to the class.
    """
    return {
        name: getattr(ps, name)
        for name in ps.__all__
        if inspect.isclass(getattr(ps, name))
    }


def get_public_members(public_class: type) -> dict[str, tuple[type, Any]]:
    """Returns the public members of a public class.

    A member is public when its name is clean, a class in the package defines it, and it
    is not deprecated. Members inherited from the standard library are left out. A
    member defined by several classes in the MRO is read from the first one, as
    attribute lookup does.

    :param public_class: The public class.
    :return: A dict mapping each public member's name to the class that defines it and
        the member's raw value in that class's namespace.
    """
    members: dict[str, tuple[type, Any]] = {}
    for klass in public_class.__mro__:
        if not klass.__module__.startswith("pterasoftware."):
            continue
        for name, value in vars(klass).items():
            if not name.startswith("_") and name not in members:
                members[name] = (klass, value)
    return {
        name: (klass, value)
        for name, (klass, value) in members.items()
        if f"{klass.__module__}.{klass.__qualname__}.{name}" not in DEPRECATED_MEMBERS
    }


def get_public_callables() -> dict[str, tuple[Callable[..., Any], str, bool]]:
    """Returns the public callables.

    These are the public functions, each public class's constructor, and each public
    class's public methods, class methods, and static methods.

    :return: A dict mapping each public callable's label to the callable as users reach
        it, its fully qualified name as the deprecation scan writes it, and whether
        users reach it with a leading self parameter still in its signature.
    """
    callables: dict[str, tuple[Callable[..., Any], str, bool]] = {}
    for name in ps.__all__:
        obj = getattr(ps, name)
        if not inspect.isclass(obj):
            callables[name] = (obj, f"{obj.__module__}.{obj.__qualname__}", False)
    for class_name, public_class in get_public_classes().items():
        init_class = next(
            klass for klass in public_class.__mro__ if "__init__" in vars(klass)
        )
        callables[f"{class_name}.__init__"] = (
            vars(init_class)["__init__"],
            f"{init_class.__module__}.{init_class.__qualname__}",
            True,
        )
        for name, (klass, value) in get_public_members(public_class).items():
            if not (
                isinstance(value, (classmethod, staticmethod))
                or inspect.isfunction(value)
            ):
                continue
            callables[f"{class_name}.{name}"] = (
                getattr(public_class, name),
                f"{klass.__module__}.{klass.__qualname__}.{name}",
                inspect.isfunction(value),
            )
    return callables


def describe_parameters(
    callable_: Callable[..., Any], has_self: bool, deprecated: list[str]
) -> tuple[str, ...]:
    """Returns a callable's parameters in the form that PUBLIC_PARAMETERS pins.

    :param callable_: The callable, as users reach it.
    :param has_self: Whether the callable's signature starts with a self parameter,
        which is left out.
    :param deprecated: The names of the callable's deprecated parameters, which are left
        out.
    :return: The parameter names in order, with a "/" after the positional-only
        parameters and a "*" before the keyword-only ones.
    """
    parameters = list(inspect.signature(callable_).parameters.values())
    if has_self:
        parameters = parameters[1:]
    described: list[str] = []
    previous_kind: Any = None
    for parameter in parameters:
        if parameter.name in deprecated:
            continue
        kind = parameter.kind
        if (
            previous_kind is inspect.Parameter.POSITIONAL_ONLY
            and kind is not inspect.Parameter.POSITIONAL_ONLY
        ):
            described.append("/")
        if kind is inspect.Parameter.KEYWORD_ONLY and previous_kind not in (
            inspect.Parameter.KEYWORD_ONLY,
            inspect.Parameter.VAR_POSITIONAL,
        ):
            described.append("*")
        if kind is inspect.Parameter.VAR_POSITIONAL:
            described.append(f"*{parameter.name}")
        elif kind is inspect.Parameter.VAR_KEYWORD:
            described.append(f"**{parameter.name}")
        else:
            described.append(parameter.name)
        previous_kind = kind
    if previous_kind is inspect.Parameter.POSITIONAL_ONLY:
        described.append("/")
    return tuple(described)


def find_internal_types(annotation: str, namespace: dict[str, Any]) -> list[str]:
    """Returns each dotted name in an annotation that names an internal class or does
    not resolve.

    Each dotted name is resolved against the namespace of the module where the
    annotation is written, and then against the builtins. A class is internal when the
    package defines it and does not list it in __all__. A name that does not resolve at
    runtime, such as one imported only for type checking, is returned too, since its
    type cannot be checked.

    :param annotation: The annotation's source text.
    :param namespace: The namespace of the module where the annotation is written.
    :return: The offending dotted names, each followed by " (unresolved)" when it does
        not resolve.
    """
    chains: list[str] = []

    def collect(node: ast.AST) -> None:
        """Collects the maximal dotted names in an annotation's expression, skipping
        constants such as Literal values.

        :param node: The node to collect from.
        :return: None
        """
        parts: list[str] = []
        current = node
        while isinstance(current, ast.Attribute):
            parts.append(current.attr)
            current = current.value
        if isinstance(current, ast.Name):
            chains.append(".".join([current.id, *reversed(parts)]))
            return
        for child in ast.iter_child_nodes(node):
            collect(child)

    collect(ast.parse(annotation, mode="eval"))

    offending: list[str] = []
    for chain in chains:
        head, *rest = chain.split(".")
        if head in namespace:
            value = namespace[head]
        elif hasattr(builtins, head):
            value = getattr(builtins, head)
        else:
            offending.append(f"{chain} (unresolved)")
            continue
        try:
            for part in rest:
                value = getattr(value, part)
        except AttributeError:
            offending.append(f"{chain} (unresolved)")
            continue
        if not (
            inspect.isclass(value) and value.__module__.startswith("pterasoftware.")
        ):
            continue
        if value.__name__ not in ps.__all__ or getattr(ps, value.__name__) is not value:
            offending.append(chain)
    return offending


def get_attribute_annotations(klass: type) -> dict[str, list[str]]:
    """Returns the annotations that a class's source gives its instance attributes.

    These are the annotated assignments to self attributes anywhere in the class's body,
    such as an annotated assignment in its constructor.

    :param klass: The class.
    :return: A dict mapping each annotated attribute's name to the source text of each
        of its annotations.
    """
    annotations: dict[str, list[str]] = {}
    for node in ast.walk(ast.parse(textwrap.dedent(inspect.getsource(klass)))):
        if (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Attribute)
            and isinstance(node.target.value, ast.Name)
            and node.target.value.id == "self"
        ):
            annotations.setdefault(node.target.attr, []).append(
                ast.unparse(node.annotation)
            )
    return annotations


def find_internal_type_uses() -> dict[tuple[str, str], list[str]]:
    """Returns every public surface that names an internal type, ignoring the allowlist.

    The surfaces are the parameter and return annotations of every public callable, the
    getter and setter annotations of every public property, and the annotations of every
    public slot attribute, read from the self assignments in the source of the classes
    that define them. Deprecated parameters are skipped.

    :return: A dict mapping each offending surface's label and its parameter, attribute,
        or "return" to the offending dotted names.
    """
    uses: dict[tuple[str, str], list[str]] = {}

    def check_annotations(
        this_label: str, function: Callable[..., Any], skipped: list[str]
    ) -> None:
        """Records the internal types named by a function's annotations.

        :param this_label: The label of the surface the function belongs to.
        :param function: The function, which is unwrapped to read its own annotations.
        :param skipped: The names of the parameters to skip.
        :return: None
        """
        function = inspect.unwrap(function)
        namespace = vars(sys.modules[function.__module__])
        for this_name, this_annotation in inspect.get_annotations(function).items():
            if this_name in skipped:
                continue
            if not isinstance(this_annotation, str):
                uses[(this_label, this_name)] = [f"{this_annotation!r} (not postponed)"]
                continue
            this_offending = find_internal_types(this_annotation, namespace)
            if this_offending:
                uses[(this_label, this_name)] = this_offending

    for label, (callable_, qualified_name, _) in get_public_callables().items():
        check_annotations(
            label,
            getattr(callable_, "__func__", callable_),
            DEPRECATED_PARAMETERS.get(qualified_name, []),
        )

    attribute_annotations: dict[type, dict[str, list[str]]] = {}
    for class_name, public_class in get_public_classes().items():
        for name, (_, value) in get_public_members(public_class).items():
            if isinstance(value, property):
                for accessor in (value.fget, value.fset):
                    if accessor is not None:
                        check_annotations(f"{class_name}.{name}", accessor, [])
            elif isinstance(value, types.MemberDescriptorType):
                offending = []
                for klass in public_class.__mro__:
                    if not klass.__module__.startswith("pterasoftware."):
                        continue
                    if klass not in attribute_annotations:
                        attribute_annotations[klass] = get_attribute_annotations(klass)
                    for annotation in attribute_annotations[klass].get(name, []):
                        offending.extend(
                            find_internal_types(
                                annotation, vars(sys.modules[klass.__module__])
                            )
                        )
                if offending:
                    uses[(class_name, name)] = offending
    return uses


class TestPublicNames(unittest.TestCase):
    """Tests that pin the package's public names."""

    def test_all_matches_the_pinned_names(self) -> None:
        """Test that the package's __all__ lists exactly the pinned public names.

        :return: None
        """
        self.assertEqual(
            tuple(sorted(ps.__all__)),
            PUBLIC_NAMES,
            msg=(
                "The package's __all__ no longer matches PUBLIC_NAMES. If the public API "
                "change is intended, update the pins in this module."
            ),
        )


class TestPublicMembers(unittest.TestCase):
    """Tests that pin each public class's public members."""

    def test_every_public_class_has_a_pinned_member_list(self) -> None:
        """Test that the pinned member lists cover exactly the public classes.

        :return: None
        """
        self.assertEqual(sorted(get_public_classes()), sorted(PUBLIC_MEMBERS))

    def test_every_public_class_has_exactly_its_pinned_members(self) -> None:
        """Test that each public class's public members, including those inherited from
        internal parents, are exactly its pinned members.

        :return: None
        """
        for name, public_class in get_public_classes().items():
            with self.subTest(name=name):
                self.assertEqual(
                    tuple(sorted(get_public_members(public_class))),
                    PUBLIC_MEMBERS.get(name),
                    msg=(
                        f"The public members of {name} no longer match "
                        "PUBLIC_MEMBERS. If the public API change is intended, update "
                        "the pins in this module."
                    ),
                )


class TestPublicParameters(unittest.TestCase):
    """Tests that pin each public callable's parameter names and kinds."""

    def test_every_public_callable_has_a_pinned_parameter_list(self) -> None:
        """Test that the pinned parameter lists cover exactly the public callables.

        :return: None
        """
        self.assertEqual(sorted(get_public_callables()), sorted(PUBLIC_PARAMETERS))

    def test_every_public_callable_has_exactly_its_pinned_parameters(self) -> None:
        """Test that each public callable's parameter names and kinds, as
        inspect.signature reports them, are exactly its pinned parameters.

        :return: None
        """
        for label, (
            callable_,
            qualified_name,
            has_self,
        ) in get_public_callables().items():
            with self.subTest(label=label):
                self.assertEqual(
                    describe_parameters(
                        callable_,
                        has_self,
                        DEPRECATED_PARAMETERS.get(qualified_name, []),
                    ),
                    PUBLIC_PARAMETERS.get(label),
                    msg=(
                        f"The parameters of {label} no longer match "
                        "PUBLIC_PARAMETERS. If the public API change is intended, "
                        "update the pins in this module."
                    ),
                )


class TestNoInternalTypesInPublicSignatures(unittest.TestCase):
    """Tests that no public surface names an internal type, other than the allowlisted
    ones."""

    internal_type_uses: dict[tuple[str, str], list[str]]

    @classmethod
    def setUpClass(cls) -> None:
        """Find every public surface that names an internal type."""
        cls.internal_type_uses = find_internal_type_uses()

    def test_no_public_surface_names_an_internal_type(self) -> None:
        """Test that every public surface naming an internal type is allowlisted.

        :return: None
        """
        self.assertEqual(
            {
                surface: names
                for surface, names in self.internal_type_uses.items()
                if surface not in INTERNAL_TYPE_ALLOWLIST
            },
            {},
            msg=(
                "These public surfaces name internal types. Annotate them with public "
                "types instead."
            ),
        )

    def test_every_allowlist_entry_still_names_an_internal_type(self) -> None:
        """Test that the allowlist holds no stale entries.

        :return: None
        """
        self.assertEqual(
            INTERNAL_TYPE_ALLOWLIST - set(self.internal_type_uses),
            set(),
            msg=(
                "These allowlisted surfaces no longer name internal types, so remove "
                "them from INTERNAL_TYPE_ALLOWLIST."
            ),
        )
