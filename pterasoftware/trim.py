"""This module contains functions to analyze the trim conditions of steady and unsteady solvers.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    analyze_steady_trim: This function attempts to calculate a trim condition of a steady solver by varying the
    angles of attack and sideslip until the steady pitching and yawing moments are zero. If a trim condition can be
    found, it returns the angles of attack and sideslip. Otherwise, it returns NaN angles and logs the failure.

    analyze_unsteady_trim: This function attempts to calculate a cycle-averaged trim condition of an unsteady solver
    by varying the angles of attack and sideslip until the unsteady, cycle-averaged pitching and yawing moments are
    zero. If a trim condition can be found, it returns the angles of attack and sideslip. Otherwise, it returns NaN
    angles and logs the failure.
"""


# ToDo: Document this function.
def analyze_steady_trim():
    """This function attempts to calculate a trim condition of a steady solver by varying the angles of attack and
    sideslip until the steady pitching and yawing moments are zero. If a trim condition can be found, it returns the
    angles of attack and sideslip. Otherwise, it returns NaN angles and logs the failure.

    :return:
    """
    pass


# ToDo: Document this function.
def analyze_unsteady_trim():
    """This function attempts to calculate a cycle-averaged trim condition of an unsteady solver by varying the
    angles of attack and sideslip until the unsteady, cycle-averaged pitching and yawing moments are zero. If a trim
    condition can be found, it returns the angles of attack and sideslip. Otherwise, it returns NaN angles and logs
    the failure.

    :return:
    """
    pass
