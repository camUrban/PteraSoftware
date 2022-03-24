import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opto


def obj(arguments):
    x, y = arguments
    return np.sin(np.sqrt(x**2 + y**2))


x_init = 5
y_init = 5
z_init = obj((x_init, y_init))
initial_guess = (x_init, y_init)
result = opto.minimize(obj, initial_guess, tol=0.01)

x_opto, y_opto = result.x
z_opto = result.fun

x_vec = np.linspace(-5, 5, 50)
y_vec = np.linspace(-5, 5, 50)

x_arr, y_arr = np.meshgrid(x_vec, y_vec)
z_arr = obj((x_arr, y_arr))

fig = plt.figure()
ax = plt.axes(projection="3d")

ax.plot_surface(x_arr, y_arr, z_arr, cmap="plasma", alpha=0.5, edgecolor="none")
plt.xlabel("x")
plt.ylabel("y")

ax.scatter(x_init, y_init, z_init, s=50, c="red", label="Guess")
ax.scatter(x_opto, y_opto, z_opto, s=50, c="blue", label="Solution")
plt.legend()
plt.show()

print(result)
