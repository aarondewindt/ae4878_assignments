from math import acos, atan2, cos, sin, sqrt, pi, tan, atan
from typing import Sequence, Optional

import numpy as np
from numpy.linalg import norm
from scipy.optimize import root_scalar

from constants import z_hat, mu_earth


def cartesian_to_kepler(r: Sequence[float], v: Sequence[float], mu: float=mu_earth):
    """
    Convert Cartesian components to Kepler elements.

    Source:
        Wertz, James Richard. "Mission geometry; orbit and constellation
        design and management." Space Technology Library (2001).

    :param r: Position vector.
    :param v: Velocity vector.
    :param mu: Optional, Standard gravitational parameter, by defaults it's Earth's (3.98600441e14).
    :return: Tuple with the semi-major axis (a), eccentricity (e), inclination (i),
             right ascension of the ascending node (Ω), argument of periapsis (ω), true anomaly (θ).
             eccentric anomaly (E) and the mean anomaly (M).
    """
    # Make sure the arguments are numpy arrays.
    r = np.asarray(r)
    v = np.asarray(v)

    r_norm = norm(r)

    # Orbit angular momentum
    h = np.cross(r, v)

    # Ascending node
    n = np.cross(z_hat, h / norm(h))
    n_hat = n / norm(n)

    # Eccentricity vector
    e_vec = np.cross(v, h) / mu - r / r_norm

    # Semi major axis
    a = 1 / (2 / r_norm - norm(v)**2/mu)

    # Eccentricity
    e = norm(e_vec)

    # Inclination
    i = acos(h[2] / norm(h))

    # Right Ascension of the Ascending Node
    raan = atan2(n[1], n[0])
    raan = (raan + 2 * np.pi) % (2 * np.pi)  # Force range between 0 and 2pi.

    # Argument of perigee
    omega = acos(e_vec @ n_hat / norm(e_vec))

    # True anomaly
    sign = 1 if np.cross(e_vec, r) @ h > 0 else -1
    theta = sign * acos(r @ e_vec / (e * r_norm))
    theta = (theta + 2 * np.pi) % (2 * np.pi)  # Force range between 0 and 2pi.

    # Eccentric anomaly
    e_ano = atan2(sqrt(1 - e**2) * sin(theta), e + cos(theta))
    e_ano = (e_ano + 2 * np.pi) % (2 * np.pi)  # Force range between 0 and 2pi.

    # Mean anomaly
    m = e_ano - e * sin(e_ano)
    m = (m + 2 * np.pi) % (2 * np.pi)  # Force range between 0 and 2pi.

    return a, e, i, raan, omega, theta, e_ano, m


def kepler_to_cartesian(a: float, e: float, i: float, raan: float, omega: float,
                        theta: Optional[float]=None, e_ano: Optional[float]=None, m: Optional[float]=None, 
                        mu: float=mu_earth):
    """
    Convert to Kepler elements to Cartesian components.

    Source:
        Wertz, James Richard. "Mission geometry; orbit and constellation
        design and management." Space Technology Library (2001).

    :param a: Semi-major axis.
    :param e: Eccentricity.
    :param i: Inclination.
    :param raan: Right ascension of the ascending node (Omega).
    :param omega: Argument of periapsis.
    :param theta: True anomaly.
    :param e_ano: Eccentric anomaly.
    :param m: Mean anomaly.
    :param mu: Optional, Standard gravitational parameter, by defaults it's Earth's (3.98600441e14).
    :return: Tuple with the cartesian position (r) and velocity (v).
    """

    # We need the true anomaly for the calculations, check if it was given.
    if theta is None:
        # If not we can calculate it from the eccentric anomaly, check if it was given.
        if e_ano is None:
            # If not we can calculate the eccentric anomaly from the mean anomaly.
            # Check if it was given.
            if m is None:
                # We need at least one of the three anomalies. Raise error.
                raise ValueError('At least one of "True anomaly", "Eccentric anomaly" or "Mean anomaly" '
                                 'is required')
            else:
                e_ano = eccentric_anomaly_from_mean_anomaly(e, m)

        # We should have a value for the eccentric anomaly at this point.
        theta = true_anomaly_from_eccentric_anomaly(e, e_ano)

    # Semiparameter
    p = a * (1 - e**2)

    # Position in the perifocal coordinate system (pf)
    r_pf = np.zeros(3)
    r_pf[0] = p * cos(theta) / (1 + e * cos(theta))
    r_pf[1] = p * sin(theta) / (1 + e * cos(theta))

    # Velocity in the perifocal coordinate system
    v_pf = np.zeros(3)
    v_pf[0] = -sqrt(mu / p) * sin(theta)
    v_pf[1] = sqrt(mu / p) * (e + cos(theta))

    # Tranformation matrix from the perifocal to inertial coordinate system.
    cr = cos(raan)
    co = cos(omega)
    ci = cos(i)
    sr = sin(raan)
    so = sin(omega)
    si = sin(i)

    c_pf = np.array([[cr*co - sr*so*ci, -cr*so - sr*co*ci, sr*si],
                     [sr*co + cr*so*ci, -sr*so + cr*co*ci, -cr*si],
                     [so*si, co*si, ci]])

    # Position and velocity in the inertial frame
    r = c_pf @ r_pf
    v = c_pf @ v_pf

    return r, v


def eccentric_anomaly_from_mean_anomaly(e: float, m: float):
    """
    Calculates the eccentric anomaly from the mean anomaly using a
    numerical rootfinding method.
    
    :param e: Eccentricity.
    :param m: Mean anomaly.
    :return: Mean anomaly
    """
    # Equation to solve. This is the equation to calculate the mean
    # anomaly from the eccentric anomaly minus the mean anomaly.
    # Thus the solusion will be at the root.
    def f(e_ano):
        return e_ano - e * sin(e_ano) - m
    
    # Solve it numerically, we're letting scipy choose a method for us here.
    solution = root_scalar(f, x0=m, bracket=[0, 2*pi])
    return solution.root


def true_anomaly_from_eccentric_anomaly(e: float, e_ano):
    """
    Calculates the true anomaly from the eccentric anomaly
    
    :param e: Eccentricity.
    :param m: Eccentric anomaly.
    :return: True anomaly
    """
    return 2 * atan(sqrt((1+e)/(1-e)) * tan(e_ano/2))
    