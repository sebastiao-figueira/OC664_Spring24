# Entrained sand load
import numpy as np

# Constants needed:
critical_shields = 0.05
m = 11.0  # proportionality constant, van der A 2013
n = 1.2  # power of the excess Shields parameter, van der A 2013

# Variables from other formulas:
# (magnitude of) Shields parameter at crest, empty array for now
shields_c = np.empty(5, dtype=float)
# (magnitude of) Shields parameter at trough, empty array for now
shields_t = np.empty(5, dtype=float)

# Calculate sand load for half-cycles: omega_c or omega_t
# If shields_c and shields_t are actually arrays, will need to loop through each value,
# but this is the correct logic
if shields_c > critical_shields:
    omega_c = m * (shields_c - critical_shields) ** np
else:
    omega_c = 0
if shields_t > critical_shields:
    omega_t = m * (shields_t - critical_shields) ** np
else:
    omega_c = 0

# then omega_c and omega_t go into Phase lag formulas
