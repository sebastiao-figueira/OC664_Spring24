# Entrained sand load
import numpy as np

# Constants needed:
rho_s = 2650
rho = 1000
nu = 0.00001
d = d50
g = 9.81
m = 11.0  # proportionality constant, van der A 2013
n = 1.2  # power of the excess Shields parameter, van der A 2013

# D_star (Soulsby 1997)
D_star = d * (g * ((rho_s / rho) - 1) / n**2)**(1/3)
# Critical Shields (Soulsby 1997)
critical_shields = 0.30/(1 + 1.2 * D_star) + 0.055 * (1 - np.exp(-0.020)*D_star)

# Variables from other formulas:
# (magnitude of) Shields parameter at crest, empty array for now
shields_c = np.empty(5, dtype=float)
# (magnitude of) Shields parameter at trough, empty array for now
shields_t = np.empty(5, dtype=float)

# Calculate sand load for half-cycles: omega_c or omega_t
# If shields_c and shields_t are actually arrays, will need to loop through each value,
# but this is the correct logic
if shields_c > critical_shields:
    omega_c = m * (shields_c - critical_shields)**n
else:
    omega_c = 0
if shields_t > critical_shields:
    omega_t = m * (shields_t - critical_shields)**n
else:
    omega_c = 0

# then omega_c and omega_t go into Phase lag formulas
