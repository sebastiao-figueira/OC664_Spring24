# -*- coding: utf-8 -*-
"""
Created on Mon May 13 21:22:03 2024

@author: sebas
"""

import numpy as np

# omega = []
# wave_period = []
# shields = []


def VDA_13_sediment_transport(omega, wave_period, shields, rho=1000, rho_s=2650, d_50=0.2, g=9.81):
    """
    Calculate the net transport rate induced by non-breaking waves and currents
    based on the method set forth by van der A (2013).

    Args:
        omega (array): Sand load entrainment and transport rates during different wave phases.
        wave_period (array): Durations of wave phases.
        shields (array): Magnitudes of Shields parameters and bed shear stresses.
        rho (float): Water density (default: 1000 kg/m^3).
        rho_s (float): Sediment density (default: 2650 kg/m^3).
        d_50 (float): Median sediment grain size (default: 0.2 mm).
        g (float): Acceleration due to gravity (default: 9.81 m/s^2).

    Returns:
        float: Net sediment transport rate.
    """
    
    # Check input array sizes
    if len(omega) != 4 or len(wave_period) != 5 or len(shields) != 4:
        raise ValueError("Input arrays must have sizes [4], [5], and [4] respectively.")

    # Convert d_50 to meters
    d_50 = d_50 / 1000
    
    # Calculate s parameter
    s = (rho_s - rho) / rho 
    
    # Extract values from input arrays
    omega_cc, omega_ct, omega_tt, omega_tc = omega
    T, T_c, T_cu, T_t, T_tu = wave_period
    shields_c, shields_t, shields_cx, shields_tx = shields
    
    # Calculate transport rates
    q_c = np.sqrt(shields_c) * T_c * (omega_cc + T_c / (2 * T_cu) * omega_tc) * shields_cx / shields_c
    q_t = np.sqrt(shields_t) * T_t * (omega_tt + T_t / (2 * T_tu) * omega_ct) * shields_tx / shields_t
    
    # Calculate net sediment transport rate
    q_s = (q_c + q_t) * np.sqrt((s - 1) * g * d_50 ** 3) / T
    
    # Sum sediment transport over selected time period
    q_sum = sum(q_s)
    
    return q_s, q_sum
