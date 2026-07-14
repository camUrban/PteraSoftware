"""Contains functions for calculating induced velocities."""

from __future__ import annotations

import math
import warnings

import numba
import numpy as np
from numba import njit, prange

from . import _logging

_logger = _logging.get_logger("_aerodynamics_functions")

# Squire's parameter relates to the size of the vortex cores and the rate at which they
# grow. The value of this parameter is slightly controversial. It dramatically affects
# the stability of the result. I'm using this value, as cited for use in flapping-wing
# vehicles in "Role of Filament Strain in the Free-Vortex Modeling of Rotor Wakes"
# (Ananthan and Leishman, 2004). It is unitless.
_squire = 1.0e-4

# Lamb's constant relates to the size of the vortex cores and the rate at which they
# grow. The value of this parameter is well agreed upon, and published in "Extended
# Unsteady Vortex-Lattice Method for Insect Flapping Wings" (Nguyen et al., 2016). It is
# unitless.
_lamb = 1.25643

# The local machine error is used to detect degenerate (zero length) line vortices.
_eps = np.finfo(float).eps

# The relative tolerance for scale invariant singularity checks. It provides two orders
# of magnitude of safety margin before catastrophic cancellation in 1.0 - cos(theta)
# begins at theta ~ 2.1e-8.
_tol = 1.0e-10

# Pre compute 4 * pi and 4.0 * _lamb as they used repeatedly.
_four_pi = 4.0 * math.pi
_four_lamb = 4.0 * _lamb

# The smallest number of evaluations worth handing one thread of a parallel kernel
# launch, where one evaluation is computing one line vortex's induced velocity at one
# point. Sized so a thread's share of the work takes several times longer than the cost
# of rallying the thread pool: waking its threads, splitting the work, and waiting for
# the slowest member to finish. A launch too small to clear one grain would spend more
# time rallying than computing, so it runs on a single thread. The grain is expressed as
# work per thread rather than as a core count, so it transfers across machines of
# different sizes without per-machine calibration.
_GRAIN = 10_000


def _ceiling() -> int:
    """Returns the largest share of the Numba thread pool that a kernel launch may use.

    The share is three quarters of the pool's width, rounded down, never below 1. A team
    that claims every logical processor stalls on its closing barrier whenever anything
    else on the machine wants CPU time, because the scheduler delays one member and the
    barrier waits for that laggard. Three quarters leaves enough slack to absorb
    realistic background bursts, and it costs nothing where it binds, having measured as
    the fastest width in every non-idle condition tested. The ceiling applies uniformly,
    including to thread masks set at or above it. Masks below the ceiling are honored
    exactly.

    The pool's width is read here rather than at import, because a module-level read
    would do arithmetic on a Numba attribute while the module is being imported, which
    fails wherever Numba is not fully importable (Sphinx's autodoc mocks it, for one).
    Reading it is cheap and, unlike numba.get_num_threads, does not launch the thread
    pool.

    :return: A positive int representing the most threads a single kernel launch may
        use.
    """
    # Numba populates its config module's attributes dynamically, so mypy sees this one
    # as untyped. Convert it explicitly rather than returning an untyped value.
    pool_width = int(numba.config.NUMBA_NUM_THREADS)
    return max((3 * pool_width) // 4, 1)


def _threads_for_launch(num_points: int, num_vortices: int) -> int:
    """Returns the work-proportional thread count for a parallel kernel launch.

    :param num_points: A non negative int representing the number of points at which the
        launch computes induced velocities.
    :param num_vortices: A non negative int representing the number of line vortices
        whose induced velocities the launch computes at each point.
    :return: A positive int representing the number of threads the launch's work
        justifies: one thread per whole grain of evaluations, floored at 1 and capped at
        the launch's point count.
    """
    evaluations = num_points * num_vortices

    # The kernels parallelize over points alone, summing each point's vortices in a
    # serial inner loop, so a launch can never use more threads than it has points.
    # Handing it more would leave the surplus threads with no iteration to run while
    # they still pay the parallel region's entry and its closing barrier, which is the
    # very overhead this dispatch exists to avoid. The cap binds when a launch carries
    # more than one grain of vortices and fewer points than the ceiling, as an unsteady
    # run's wake influences do on a small mesh once the wake grows, and as streamline
    # launches do generally.
    return min(max(evaluations // _GRAIN, 1), max(num_points, 1))


def report_thread_settings() -> None:
    """Reports the thread dispatch settings that a solver run's parallel kernel launches
    will operate under.

    Logged at the debug level: the active Numba threading layer, the Numba thread pool's
    width, the current thread mask, the kernel thread ceiling, and the upper bound those
    put on a launch's thread count. A launch whose work falls below the grain uses fewer
    threads than that bound.

    Raises a UserWarning, rather than logging one, if Numba is running on its workqueue
    threading layer. That layer slows every simulation, users land on it without being
    told, and the fix is a change to the user's environment rather than to anything
    Ptera Software can do, which is the case the warnings module exists for. The
    warnings module also suppresses the repeats itself, which matters because a trim or
    convergence analysis runs the solver many times over and Numba fixes the threading
    layer for the life of the process.

    :return: None
    """
    # Read the thread mask before the threading layer. Reading the mask launches the
    # thread pool if it isn't already running, and launching the pool initializes the
    # threading layer, so the layer is only queryable afterwards. The records below are
    # emitted in the opposite order, which the reads do not constrain.
    external_cap = numba.get_num_threads()
    threading_layer = numba.threading_layer()
    ceiling = _ceiling()
    kernel_cap = min(external_cap, ceiling)

    if threading_layer == "workqueue":
        warnings.warn(
            "Numba is running on its workqueue threading layer, which it falls "
            "back to when neither TBB nor OpenMP can be loaded. That layer charges "
            "a large fixed cost to every parallel kernel launch, no matter how "
            "small the launch or how few threads it uses, which slows every "
            "simulation. Install tbb, or a system OpenMP runtime, to let Numba use "
            "a faster layer. Numba also falls back to workqueue when an installed "
            "tbb cannot be loaded, so having tbb installed does not rule this out."
        )

    _logger.debug(_logging.indent() + f"Numba threading layer: {threading_layer}")
    _logger.debug(
        _logging.indent()
        + f"Numba thread pool width: {numba.config.NUMBA_NUM_THREADS}, thread mask: "
        + f"{external_cap}, ceiling: {ceiling} (three quarters of the pool)"
    )
    _logger.debug(
        _logging.indent()
        + f"Kernel launches will use at most {kernel_cap} threads, and fewer when a "
        + "launch's work is below the grain"
    )


def collapsed_velocities_from_ring_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrrvp_GP1_CgP1: np.ndarray,
    stackFrrvp_GP1_CgP1: np.ndarray,
    stackFlrvp_GP1_CgP1: np.ndarray,
    stackBlrvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of ring vortices and
    finds the cumulative induced velocity at every point.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations. Each kernel
    launch runs with a work-proportional thread count, which never exceeds the launch's
    point count, the current Numba thread mask, or the module's kernel thread ceiling.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M ring vortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' front left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackBlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' back left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M ring
        vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        ring vortices. Based on results from Ramasamy and Leishman (2007), a reasonable
        value that works across scales is 3.0% the chord length of each line vortex's
        parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across all four legs.
    :param ages: For bound ring vortices, this must be None. For ring vortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M ring vortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrrvp_GP1_CgP1,
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
        stackBrrvp_GP1_CgP1,
    ]

    stackVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], 3), dtype=float)

    # Read the current thread mask so the dispatch honors any cap the user has set, and
    # restore it once the kernel launches finish.
    external_cap = numba.get_num_threads()
    numba.set_num_threads(
        min(
            _threads_for_launch(stackP_GP1_CgP1.shape[0], strengths.shape[0]),
            external_cap,
            _ceiling(),
        )
    )

    # Get the velocity induced by each leg of each ring vortex (in the first Airplane's
    # geometry axes, observed from the Earth frame).
    try:
        for i in range(4):
            stackVInd_GP1__E += _collapsed_velocities_from_line_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
                stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
                strengths=strengths,
                r_c0s=r_c0s,
                singularity_counts=singularity_counts,
                ages=ages,
                nu=nu,
            )
    finally:
        numba.set_num_threads(external_cap)
    return stackVInd_GP1__E


def collapsed_velocities_from_ring_vortices_chordwise_segments(
    stackP_GP1_CgP1: np.ndarray,
    stackBrrvp_GP1_CgP1: np.ndarray,
    stackFrrvp_GP1_CgP1: np.ndarray,
    stackFlrvp_GP1_CgP1: np.ndarray,
    stackBlrvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of ring vortices and
    finds the cumulative induced velocity at every point due to the ring vortices' left
    and right line vortex legs.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations. Each kernel
    launch runs with a work-proportional thread count, which never exceeds the launch's
    point count, the current Numba thread mask, or the module's kernel thread ceiling.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M ring vortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' front left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackBlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' back left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M ring
        vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        ring vortices. Based on results from Ramasamy and Leishman (2007), a reasonable
        value that works across scales is 3.0% the chord length of each line vortex's
        parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across both legs.
    :param ages: For bound ring vortices, this must be None. For ring vortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M ring vortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame) due to the M ring vortices' left and right line vortex legs. The units
        are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
    ]

    stackVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], 3), dtype=float)

    # Read the current thread mask so the dispatch honors any cap the user has set, and
    # restore it once the kernel launches finish.
    external_cap = numba.get_num_threads()
    numba.set_num_threads(
        min(
            _threads_for_launch(stackP_GP1_CgP1.shape[0], strengths.shape[0]),
            external_cap,
            _ceiling(),
        )
    )

    # Get the velocity induced by the left and right legs of each ring vortex (in the
    # first Airplane's geometry axes, observed from the Earth frame).
    try:
        for i in range(2):
            stackVInd_GP1__E += _collapsed_velocities_from_line_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
                stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
                strengths=strengths,
                r_c0s=r_c0s,
                singularity_counts=singularity_counts,
                ages=ages,
                nu=nu,
            )
    finally:
        numba.set_num_threads(external_cap)
    return stackVInd_GP1__E


def expanded_velocities_from_ring_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrrvp_GP1_CgP1: np.ndarray,
    stackFrrvp_GP1_CgP1: np.ndarray,
    stackFlrvp_GP1_CgP1: np.ndarray,
    stackBlrvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of ring vortices and
    finds the induced velocity at every point due to each ring vortex.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations. Each kernel
    launch runs with a work-proportional thread count, which never exceeds the launch's
    point count, the current Numba thread mask, or the module's kernel thread ceiling.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M ring vortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' front left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackBlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M ring vortices' back left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M ring
        vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        ring vortices. Based on results from Ramasamy and Leishman (2007), a reasonable
        value that works across scales is 3.0% the chord length of each line vortex's
        parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across all four legs.
    :param ages: For bound ring vortices, this must be None. For ring vortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M ring vortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M ring vortices. The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrrvp_GP1_CgP1,
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
        stackBrrvp_GP1_CgP1,
    ]

    gridVInd_GP1__E = np.zeros(
        (stackP_GP1_CgP1.shape[0], strengths.shape[0], 3), dtype=float
    )

    # Read the current thread mask so the dispatch honors any cap the user has set, and
    # restore it once the kernel launches finish.
    external_cap = numba.get_num_threads()
    numba.set_num_threads(
        min(
            _threads_for_launch(stackP_GP1_CgP1.shape[0], strengths.shape[0]),
            external_cap,
            _ceiling(),
        )
    )

    # Get the velocity induced by each leg of each ring vortex (in the first Airplane's
    # geometry axes, observed from the Earth frame).
    try:
        for i in range(4):
            gridVInd_GP1__E += _expanded_velocities_from_line_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
                stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
                strengths=strengths,
                r_c0s=r_c0s,
                singularity_counts=singularity_counts,
                ages=ages,
                nu=nu,
            )
    finally:
        numba.set_num_threads(external_cap)
    return gridVInd_GP1__E


# TODO: Remove the ability to specify horseshoe vortex ages, since they are never used
#  in unsteady simulations.
def collapsed_velocities_from_horseshoe_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrhvp_GP1_CgP1: np.ndarray,
    stackFrhvp_GP1_CgP1: np.ndarray,
    stackFlhvp_GP1_CgP1: np.ndarray,
    stackBlhvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of horseshoe vortices
    and finds the cumulative induced velocity at every point.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations. Each kernel
    launch runs with a work-proportional thread count, which never exceeds the launch's
    point count, the current Numba thread mask, or the module's kernel thread ceiling.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M horseshoe vortices' back right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M horseshoe vortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M horseshoe vortices' front left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackBlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M horseshoe vortices' back left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        horseshoe vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        horseshoe vortices. Based on results from Ramasamy and Leishman (2007), a
        reasonable value that works across scales is 3.0% the chord length of each line
        vortex's parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across all three legs.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrhvp_GP1_CgP1,
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
        stackBlhvp_GP1_CgP1,
    ]

    stackVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], 3), dtype=float)

    # Read the current thread mask so the dispatch honors any cap the user has set, and
    # restore it once the kernel launches finish.
    external_cap = numba.get_num_threads()
    numba.set_num_threads(
        min(
            _threads_for_launch(stackP_GP1_CgP1.shape[0], strengths.shape[0]),
            external_cap,
            _ceiling(),
        )
    )

    # Get the velocity induced by each leg of each horseshoe vortex (in the first
    # Airplane's geometry axes, observed from the Earth frame).
    try:
        for i in range(3):
            stackVInd_GP1__E += _collapsed_velocities_from_line_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
                stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
                strengths=strengths,
                r_c0s=r_c0s,
                singularity_counts=singularity_counts,
                ages=None,
                nu=nu,
            )
    finally:
        numba.set_num_threads(external_cap)
    return stackVInd_GP1__E


# TODO: Remove the ability to specify horseshoe vortex ages, since they are never used
#  in unsteady simulations.
def expanded_velocities_from_horseshoe_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrhvp_GP1_CgP1: np.ndarray,
    stackFrhvp_GP1_CgP1: np.ndarray,
    stackFlhvp_GP1_CgP1: np.ndarray,
    stackBlhvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of horseshoe vortices
    and finds the induced velocity at every point due to each horseshoe vortex.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations. Each kernel
    launch runs with a work-proportional thread count, which never exceeds the launch's
    point count, the current Numba thread mask, or the module's kernel thread ceiling.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M horseshoe vortices' back right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M horseshoe vortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M horseshoe vortices' front left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackBlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M horseshoe vortices' back left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of M horseshoe
        vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        horseshoe vortices. Based on results from Ramasamy and Leishman (2007), a
        reasonable value that works across scales is 3.0% the chord length of each line
        vortex's parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across all three legs.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M horseshoe vortices. The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrhvp_GP1_CgP1,
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
        stackBlhvp_GP1_CgP1,
    ]

    gridVInd_GP1__E = np.zeros(
        (stackP_GP1_CgP1.shape[0], strengths.shape[0], 3), dtype=float
    )

    # Read the current thread mask so the dispatch honors any cap the user has set, and
    # restore it once the kernel launches finish.
    external_cap = numba.get_num_threads()
    numba.set_num_threads(
        min(
            _threads_for_launch(stackP_GP1_CgP1.shape[0], strengths.shape[0]),
            external_cap,
            _ceiling(),
        )
    )

    # Get the velocity induced by each leg of each horseshoe vortex (in the first
    # Airplane's geometry axes, observed from the Earth frame).
    try:
        for i in range(3):
            gridVInd_GP1__E += _expanded_velocities_from_line_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
                stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
                strengths=strengths,
                r_c0s=r_c0s,
                singularity_counts=singularity_counts,
                ages=None,
                nu=nu,
            )
    finally:
        numba.set_num_threads(external_cap)
    return gridVInd_GP1__E


@njit(cache=True, fastmath=False, parallel=True)
def _collapsed_velocities_from_line_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackSlvp_GP1_CgP1: np.ndarray,
    stackElvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of line vortices and
    finds the cumulative induced velocity at every point.

    This function uses a modified version of the Biot-Savart law to create a smooth
    induced velocity decay based on a line vortex's core radius. The core radius grows
    from an initial value based on the line vortex's age.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    **Citation:**

    Core radius equation adapted from Eq. 3 of: "A Reynolds Number-Based Blade Tip
    Vortex Model"

    Authors: Manikandan Ramasamy and J. Gordon Leishman

    Biot-Savart equation adapted from: "Extended Unsteady Vortex-Lattice Method for
    Insect Flapping Wings"

    Authors: Anh Tuan Nguyen, Joong-Kwan Kim, Jong-Seob Han, and Jae-Hung Han

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackSlvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M line vortices' starting vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackElvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M line vortices' ending vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M line
        vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        line vortices. Based on results from Ramasamy and Leishman (2007), a reasonable
        value that works across scales is 3.0% the chord length of each line vortex's
        parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across calls.
    :param ages: For bound line vortices, this must be None. For line vortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M line vortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    num_vortices = stackSlvp_GP1_CgP1.shape[0]
    num_points = stackP_GP1_CgP1.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in the
    # first Airplane's geometry axes, observed from the Earth frame).
    stackVInd_GP1__E = np.zeros((num_points, 3))

    # If the user didn't specify any ages, set the age of each line vortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    # Pre compute per line vortex quantities in a serial pass. This hoists invariant
    # work out of the parallel point loop and correctly counts degenerate filaments,
    # which are a per line vortex property, not per point.
    vortex_valid = np.empty(num_vortices, dtype=np.bool_)
    vortex_c1 = np.empty(num_vortices)
    vortex_c2 = np.empty(num_vortices)
    vortex_r0_times_tol = np.empty(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
        Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]

        # The r0_GP1 vector goes from the line vortex's start point to its end point (in
        # the first Airplane's geometry axes).
        r0X_GP1 = Elvp_GP1_CgP1[0] - Slvp_GP1_CgP1[0]
        r0Y_GP1 = Elvp_GP1_CgP1[1] - Slvp_GP1_CgP1[1]
        r0Z_GP1 = Elvp_GP1_CgP1[2] - Slvp_GP1_CgP1[2]

        # Find r0_GP1's length.
        r0 = math.sqrt(r0X_GP1**2.0 + r0Y_GP1**2.0 + r0Z_GP1**2.0)

        # Skip degenerate filaments where the start and end points coincide.
        if r0 < _eps:
            singularity_counts[0] += 1
            vortex_valid[vortex_id] = False
            continue

        vortex_valid[vortex_id] = True

        strength = strengths[vortex_id]
        age = ages[vortex_id]
        r_c0 = r_c0s[vortex_id]

        # Pre compute r0 * _tol outside the parallel loop.
        vortex_r0_times_tol[vortex_id] = r0 * _tol

        # Calculate the radius of the line vortex's core squared. The initial core
        # radius ensures nonzero regularization even for bound vortices with zero age.
        r_c_sq = r_c0**2.0 + _four_lamb * (nu + _squire * abs(strength)) * age

        vortex_c1[vortex_id] = strength / _four_pi
        vortex_c2[vortex_id] = r0**2.0 * r_c_sq

    # Use per point singularity counts to avoid write races in the parallel loop. Index
    # mapping: [1] vertex start proximity, [2] vertex end proximity, [3] collinearity.
    # Index [0] (degenerate filament) is handled in the serial pre pass above.
    point_singularity_counts = np.zeros((num_points, 4), dtype=np.int64)

    # Numba does not annotate prange.__new__, so mypy cannot tell that calling prange
    # returns an iterable. Once this is fixed upstream, warn_unused_ignores will flag
    # the ignore below as unused.
    for point_id in prange(num_points):  # type: ignore[attr-defined]
        P_GP1_CgP1 = stackP_GP1_CgP1[point_id]

        for vortex_id in range(num_vortices):
            if not vortex_valid[vortex_id]:
                continue

            Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
            Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]

            # The r1_GP1 vector goes from P_GP1_CgP1 to the line vortex's start point
            # (in the first Airplane's geometry axes).
            r1X_GP1 = Slvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r1Y_GP1 = Slvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r1Z_GP1 = Slvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r2_GP1 vector goes from P_GP1_CgP1 to the line vortex's end point (in
            # the first Airplane's geometry axes).
            r2X_GP1 = Elvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r2Y_GP1 = Elvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r2Z_GP1 = Elvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r3_GP1 vector is the cross product of r1_GP1 and r2_GP1 (in the first
            # Airplane's geometry axes).
            r3X_GP1 = r1Y_GP1 * r2Z_GP1 - r1Z_GP1 * r2Y_GP1
            r3Y_GP1 = r1Z_GP1 * r2X_GP1 - r1X_GP1 * r2Z_GP1
            r3Z_GP1 = r1X_GP1 * r2Y_GP1 - r1Y_GP1 * r2X_GP1

            # Find the lengths of r1_GP1 and r2_GP1.
            r1 = math.sqrt(r1X_GP1**2.0 + r1Y_GP1**2.0 + r1Z_GP1**2.0)
            r2 = math.sqrt(r2X_GP1**2.0 + r2Y_GP1**2.0 + r2Z_GP1**2.0)

            # Check for singularities using scale invariant criteria. The vertex
            # proximity checks (r1/r0 and r2/r0 but refactored below to use
            # multiplication instead of slower division) guard 1/r singularities.
            r0_times_tol = vortex_r0_times_tol[vortex_id]
            if r1 < r0_times_tol:
                point_singularity_counts[point_id, 1] += 1
                continue
            if r2 < r0_times_tol:
                point_singularity_counts[point_id, 2] += 1
                continue

            # Cache squared length of r3_GP1 as it is used in the c_4 calculation.
            r3_sq = r3X_GP1**2.0 + r3Y_GP1**2.0 + r3Z_GP1**2.0

            # Find the length of r3_GP1.
            r3 = math.sqrt(r3_sq)

            # Cache r1 * r2 as it is used for the collinearity check and twice in the
            # c_4 calculation.
            r1_times_r2 = r1 * r2

            c_3 = r1X_GP1 * r2X_GP1 + r1Y_GP1 * r2Y_GP1 + r1Z_GP1 * r2Z_GP1

            # The collinearity check (r3/(r1*r2) = |sin(theta)| but with the same
            # multiplication instead of division refactor) guards catastrophic
            # cancellation in 1-cos(theta).
            if r3 < (_tol * r1_times_r2):
                # Collinearity can indicate one of two things. If the point is collinear
                # and between the filament's vertices, it is a true singularity (the
                # Biot-Savart equation diverges), so we exclude the contribution as it
                # is the influence of the filament on itself. If the point is collinear
                # and off to one side of the filament, it isn't a true singularity, as
                # the Biot-Savart equation (if calculated with infinite precision)
                # correctly returns zero induced velocity. However, we still run into
                # the catastrophic cancellation issue, so we again manually return zero
                # induced velocity contribution. These two situations are distinguished
                # by the sign of the c_3 (the dot product of r1 and r2).
                if c_3 < 0.0:
                    point_singularity_counts[point_id, 3] += 1
                continue

            c_4 = (
                vortex_c1[vortex_id]
                * (r1 + r2)
                * (r1_times_r2 - c_3)
                / (r1_times_r2 * (r3_sq + vortex_c2[vortex_id]))
            )
            stackVInd_GP1__E[point_id, 0] += c_4 * r3X_GP1
            stackVInd_GP1__E[point_id, 1] += c_4 * r3Y_GP1
            stackVInd_GP1__E[point_id, 2] += c_4 * r3Z_GP1

    # Aggregate per point singularity counts into the output array.
    for k in range(1, 4):
        for p in range(num_points):
            singularity_counts[k] += point_singularity_counts[p, k]

    return stackVInd_GP1__E


@njit(cache=True, fastmath=False, parallel=True)
def _expanded_velocities_from_line_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackSlvp_GP1_CgP1: np.ndarray,
    stackElvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of line vortices and
    finds the induced velocity at every point due to each line vortex.

    This function uses a modified version of the Biot-Savart law to create a smooth
    induced velocity decay based on a line vortex's core radius. The core radius grows
    from an initial value based on the line vortex's age.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    **Citation:**

    Core radius equation adapted from Eq. 3 of: "A Reynolds Number-Based Blade Tip
    Vortex Model"

    Authors: Manikandan Ramasamy and J. Gordon Leishman

    Biot-Savart equation adapted from: "Extended Unsteady Vortex-Lattice Method for
    Insect Flapping Wings"

    Authors: Anh Tuan Nguyen, Joong-Kwan Kim, Jong-Seob Han, and Jae-Hung Han

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackSlvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of M
        line vortices' starting vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackElvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M line vortices' ending vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M line
        vortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        line vortices. Based on results from Ramasamy and Leishman (2007), a reasonable
        value that works across scales is 3.0% the chord length of each line vortex's
        parent Wing. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across calls.
    :param ages: For bound line vortices, this must be None. For line vortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M line vortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M line vortices. The units are in meters per second.
    """
    num_vortices = stackSlvp_GP1_CgP1.shape[0]
    num_points = stackP_GP1_CgP1.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in the
    # first Airplane's geometry axes, observed from the Earth frame).
    gridVInd_GP1__E = np.zeros((num_points, num_vortices, 3))

    # If the user didn't specify any ages, set the age of each line vortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    # Pre compute per line vortex quantities in a serial pass. This hoists invariant
    # work out of the parallel point loop and correctly counts degenerate filaments,
    # which are a per line vortex property, not per point.
    vortex_valid = np.empty(num_vortices, dtype=np.bool_)
    vortex_c1 = np.empty(num_vortices)
    vortex_c2 = np.empty(num_vortices)
    vortex_r0_times_tol = np.empty(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
        Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]

        # The r0_GP1 vector goes from the line vortex's start point to its end point (in
        # the first Airplane's geometry axes).
        r0X_GP1 = Elvp_GP1_CgP1[0] - Slvp_GP1_CgP1[0]
        r0Y_GP1 = Elvp_GP1_CgP1[1] - Slvp_GP1_CgP1[1]
        r0Z_GP1 = Elvp_GP1_CgP1[2] - Slvp_GP1_CgP1[2]

        # Find r0_GP1's length.
        r0 = math.sqrt(r0X_GP1**2.0 + r0Y_GP1**2.0 + r0Z_GP1**2.0)

        # Skip degenerate filaments where the start and end points coincide.
        if r0 < _eps:
            singularity_counts[0] += 1
            vortex_valid[vortex_id] = False
            continue

        vortex_valid[vortex_id] = True

        strength = strengths[vortex_id]
        age = ages[vortex_id]
        r_c0 = r_c0s[vortex_id]

        # Pre compute r0 * _tol outside the parallel loop.
        vortex_r0_times_tol[vortex_id] = r0 * _tol

        # Calculate the radius of the line vortex's core squared. The initial core
        # radius ensures nonzero regularization even for bound vortices with zero age.
        r_c_sq = r_c0**2.0 + _four_lamb * (nu + _squire * abs(strength)) * age

        vortex_c1[vortex_id] = strength / _four_pi
        vortex_c2[vortex_id] = r0**2.0 * r_c_sq

    # Use per point singularity counts to avoid write races in the parallel loop. Index
    # mapping: [1] vertex start proximity, [2] vertex end proximity, [3] collinearity.
    # Index [0] (degenerate filament) is handled in the serial pre pass above.
    point_singularity_counts = np.zeros((num_points, 4), dtype=np.int64)

    # Numba does not annotate prange.__new__, so mypy cannot tell that calling prange
    # returns an iterable. Once this is fixed upstream, warn_unused_ignores will flag
    # the ignore below as unused.
    for point_id in prange(num_points):  # type: ignore[attr-defined]
        P_GP1_CgP1 = stackP_GP1_CgP1[point_id]

        for vortex_id in range(num_vortices):
            if not vortex_valid[vortex_id]:
                continue

            Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
            Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]

            # The r1_GP1 vector goes from P_GP1_CgP1 to the line vortex's start point
            # (in the first Airplane's geometry axes).
            r1X_GP1 = Slvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r1Y_GP1 = Slvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r1Z_GP1 = Slvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r2_GP1 vector goes from P_GP1_CgP1 to the line vortex's end point (in
            # the first Airplane's geometry axes).
            r2X_GP1 = Elvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r2Y_GP1 = Elvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r2Z_GP1 = Elvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r3_GP1 vector is the cross product of r1_GP1 and r2_GP1 (in the first
            # Airplane's geometry axes).
            r3X_GP1 = r1Y_GP1 * r2Z_GP1 - r1Z_GP1 * r2Y_GP1
            r3Y_GP1 = r1Z_GP1 * r2X_GP1 - r1X_GP1 * r2Z_GP1
            r3Z_GP1 = r1X_GP1 * r2Y_GP1 - r1Y_GP1 * r2X_GP1

            # Find the lengths of r1_GP1 and r2_GP1.
            r1 = math.sqrt(r1X_GP1**2.0 + r1Y_GP1**2.0 + r1Z_GP1**2.0)
            r2 = math.sqrt(r2X_GP1**2.0 + r2Y_GP1**2.0 + r2Z_GP1**2.0)

            # Check for singularities using scale invariant criteria. The vertex
            # proximity checks (r1/r0 and r2/r0 but refactored below to use
            # multiplication instead of slower division) guard 1/r singularities.
            r0_times_tol = vortex_r0_times_tol[vortex_id]
            if r1 < r0_times_tol:
                point_singularity_counts[point_id, 1] += 1
                continue
            if r2 < r0_times_tol:
                point_singularity_counts[point_id, 2] += 1
                continue

            # Cache squared length of r3_GP1 as it is used in the c_4 calculation.
            r3_sq = r3X_GP1**2.0 + r3Y_GP1**2.0 + r3Z_GP1**2.0

            # Find the length of r3_GP1.
            r3 = math.sqrt(r3_sq)

            # Cache r1 * r2 as it is used for the collinearity check and twice in the
            # c_4 calculation.
            r1_times_r2 = r1 * r2

            c_3 = r1X_GP1 * r2X_GP1 + r1Y_GP1 * r2Y_GP1 + r1Z_GP1 * r2Z_GP1

            # The collinearity check (r3/(r1*r2) = |sin(theta)| but with the same
            # multiplication instead of division refactor) guards catastrophic
            # cancellation in 1-cos(theta).
            if r3 < (_tol * r1_times_r2):
                # Collinearity can indicate one of two things. If the point is collinear
                # and between the filament's vertices, it is a true singularity (the
                # Biot-Savart equation diverges), so we exclude the contribution as it
                # is the influence of the filament on itself. If the point is collinear
                # and off to one side of the filament, it isn't a true singularity, as
                # the Biot-Savart equation (if calculated with infinite precision)
                # correctly returns zero induced velocity. However, we still run into
                # the catastrophic cancellation issue, so we again manually return zero
                # induced velocity contribution. These two situations are distinguished
                # by the sign of the c_3 (the dot product of r1 and r2).
                if c_3 < 0.0:
                    point_singularity_counts[point_id, 3] += 1
                continue

            c_4 = (
                vortex_c1[vortex_id]
                * (r1 + r2)
                * (r1_times_r2 - c_3)
                / (r1_times_r2 * (r3_sq + vortex_c2[vortex_id]))
            )
            gridVInd_GP1__E[point_id, vortex_id, 0] = c_4 * r3X_GP1
            gridVInd_GP1__E[point_id, vortex_id, 1] = c_4 * r3Y_GP1
            gridVInd_GP1__E[point_id, vortex_id, 2] = c_4 * r3Z_GP1

    # Aggregate per point singularity counts into the output array.
    for k in range(1, 4):
        for p in range(num_points):
            singularity_counts[k] += point_singularity_counts[p, k]

    return gridVInd_GP1__E
