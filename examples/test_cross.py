import numpy as np
import timeit
from numba import njit

num_elem = 1000

a = np.random.random(size=num_elem * 3).reshape(num_elem, 3)
b = np.random.random(size=num_elem * 3).reshape(num_elem, 3)


@njit
def nb_explicit_cross():
    e = np.zeros_like(a)
    e[:, 0] = a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1]
    e[:, 1] = a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2]
    e[:, 2] = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    return e


@njit
def nb_np_cross():
    return np.cross(a, b)


def explicit_cross():
    e = np.zeros_like(a)
    e[:, 0] = a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1]
    e[:, 1] = a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2]
    e[:, 2] = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    return e


def np_cross():
    return np.cross(a, b)


nb_explicit_cross()
nb_np_cross()
true_1 = np.allclose(nb_explicit_cross(), nb_np_cross())
true_2 = np.allclose(explicit_cross(), np_cross())
true_3 = np.allclose(nb_explicit_cross(), explicit_cross())

num_exec = 100000
num_rep = 3

print("Numba:")
print("\tExplicit:")
print(
    "\t\t"
    + str(
        round(min(timeit.repeat(nb_explicit_cross, number=num_exec, repeat=num_rep)), 3)
    )
    + " s"
)
print("\tNumpy:")
print(
    "\t\t"
    + str(round(min(timeit.repeat(nb_np_cross, number=num_exec, repeat=num_rep)), 3))
    + " s"
)
print("Python:")
print("\tExplicit:")
print(
    "\t\t"
    + str(round(min(timeit.repeat(explicit_cross, number=num_exec, repeat=num_rep)), 3))
    + " s"
)
print("\tNumpy:")
print(
    "\t\t"
    + str(round(min(timeit.repeat(np_cross, number=num_exec, repeat=num_rep)), 3))
    + " s"
)
print("Check:")
print("\t\t" + str(all([true_1, true_2, true_3])))
