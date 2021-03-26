import numpy as np


# Cartesian frame unit vectors
x_hat = np.array([1, 0, 0])
y_hat = np.array([0, 1, 0])
z_hat = np.array([0, 0, 1])

# Spaceflight constants
mu_earth = 3.98600441e14  # m^3/s^2

# Programming constants
eps64 = np.finfo(np.float64).eps
