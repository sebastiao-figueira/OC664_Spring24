# Entrained sand load
import numpy as np


# D_star (Soulsby 1997)
def dimensionless_grainsize(d: float = 0.2, g: float = 9.81, rho_s: float = 2650.0,
                            rho: float = 1000.0, nu: float = 0.00001):
    """Calculates dimensionless grain size, D*.
    Parameters
    ----------
    d : float
        Grain diameter (default is 0.2mm)
    g : float
        Acceleration due to gravity (default is 9.81)
    rho_s : float
        Density of sediment (default is 2650.0)
    rho : float
        Density of water (default is 1000.0)
    nu : float
        Kinematic viscosity of water (default is 0.00001)

    Returns
    -------
    d_star : float
        Dimensionless grain size
    """
    d = d * 0.001
    d_star = d * (g * ((rho_s / rho) - 1) / nu ** 2) ** (1 / 3)
    return d_star


# Critical Shields (Soulsby 1997)
def critical_shields(d_star):
    """Calculates critical Shields number.
    Parameters
    ----------
    d_star : float
        Dimensionless grain size

    Returns
    -------
    theta_cr : float
        Critical Shields number
    """

    theta_cr = 0.30 / (1 + 1.2 * d_star) + 0.055 * (1 - np.exp(-0.020) * d_star)
    return theta_cr


# Calculate sand load for half-cycles: omega_c or omega_t
# If shields_c and shields_t are actually arrays, will need to loop through each value,
# but this is the correct logic
def sandload_crest(theta_c, theta_cr, m: float = 11.0, n: float = 1.2):
    """Calculates sand load for crest half-cycle.
    Parameters
    ----------
    theta_c
        Shields Number for crest half cycle
    theta_cr : float
        Critical Shields number
    m : float
        Proportionality constant (default is 11.0 from van der A et al. 2013)
    n : float
        Power constant (default is 1.2 from van der A et al. 2013)
    Returns
    -------
    theta_cr : float
        Critical Shields number
    """

    omega_c = np.empty(len(theta_c))
    for i in range(len(theta_c)):
        if theta_c[i] > theta_cr:
            omega_c_single = m * (theta_c[i] - theta_cr) ** n
        else:
            omega_c_single = 0

        omega_c[i] = omega_c_single

    return omega_c


def sandload_trough(theta_t, theta_cr, m: float = 11.0, n: float = 1.2):
    """Calculates sand load for crest half-cycle.
    Parameters
    ----------
    theta_t
        Shields Number for trough half cycle
    theta_cr : float
        Critical Shields number
    m : float
        Proportionality constant (default is 11.0 from van der A et al. 2013)
    n : float
        Power constant (default is 1.2 from van der A et al. 2013)
    Returns
    -------
    theta_cr : float
        Critical Shields number
    """
    omega_t = np.empty(len(theta_t))
    for i in range(len(theta_t)):
        if theta_t[i] > theta_cr:
            omega_t_single = m * (theta_t[i] - theta_cr) ** n
        else:
            omega_t_single = 0
        omega_t[i] = omega_t_single

    return omega_t

# then omega_c and omega_t go into Phase lag formulas

# TEST
dstar = dimensionless_grainsize()
theta_critical = critical_shields(dstar)

input_thetac = np.random.rand(100, 1)
input_thetat = np.random.rand(100, 1)

sandload_c = sandload_crest(input_thetac, theta_critical)
sandload_t = sandload_trough(input_thetat, theta_critical)
