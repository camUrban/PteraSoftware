import numpy as np
import timeit
from numba import njit

num_elem = 1000

a = np.random.random(size=num_elem * 3).reshape(num_elem, 3)


@njit
def nb_explicit_norm():
    return np.sqrt((a[:, 0]) ** 2 + (a[:, 1]) ** 2 + (a[:, 2]) ** 2)


def explicit_norm():
    return np.sqrt((a[:, 0]) ** 2 + (a[:, 1]) ** 2 + (a[:, 2]) ** 2)


def np_norm():
    return np.linalg.norm(a, axis=-1)


nb_explicit_norm()
true_1 = np.allclose(nb_explicit_norm(), explicit_norm())
true_2 = np.allclose(explicit_norm(), np_norm())

num_exec = 100000
num_rep = 3

print("Numba:")
print("\tExplicit:")
print(
    "\t\t"
    + str(
        round(min(timeit.repeat(nb_explicit_norm, number=num_exec, repeat=num_rep)), 3)
    )
    + " s"
)
print("Python:")
print("\tExplicit:")
print(
    "\t\t"
    + str(round(min(timeit.repeat(explicit_norm, number=num_exec, repeat=num_rep)), 3))
    + " s"
)
print("\tNumpy:")
print(
    "\t\t"
    + str(round(min(timeit.repeat(np_norm, number=num_exec, repeat=num_rep)), 3))
    + " s"
)
print("Check:")
print("\t\t" + str(all([true_1, true_2])))
