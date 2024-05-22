
# NOTE: translated from matlab

# Define necessary constants and inputs
# ripple model equations
# required defined factors: 
    # ripple length (lambda)
    # ripple height (eta)
    # time-averaged absolute shields stress (named shields_aa here)
        # defined in category C1
    # d50, d90: grain size (grain size, static values assumed)
    # mobility number psi 
# INPUTS: 
    # CONSTANTS:
        # d50
        # d90
        # s
        # g
    # Defined elsewhere
        # ahat (category A)
        # shields_aa (category C1, can rename to fit)
    # Define here (for outputs)
        # eta (ripple height)
        # lambda (ripple length) 
        # psimax (maximum mobility number for ripple dimension eqns)
 # OUTPUTS
    # mu (fine factor adjustment for sheet flow)
    # ksdelta (current-related bed roughness)
    # ksw (wave-related bed roughness)

import numpy as np 

def ripples(d50, d90, s, g, ahat, shields_aa):
    # Maximum mobility number
    psimax = (max(u_hat))**2 / ((s - 1) * g * d50)

    # Lambda multipliers m
    if d50 <= 0.22:
        mlambda = 0.73
    elif d50 <= 0.3:
        mlambda = 0.73 + ((0.27 * (d50 - 0.22)) / (0.3 - 0.22))
    else:
        mlambda = 1

    # Eta multipliers
    if d50 <= 0.22:
        meta = 0.55
    elif d50 <= 0.3:
        meta = 0.55 + ((0.45 * (d50 - 0.22)) / (0.3 - 0.22))
    else:
        meta = 1

    # Multipliers nlambda and neta for smooth transition
    if psimax <= 190:
        n = 1
    elif psimax <= 240:
        n = 0.5 * (1 + np.cos(np.pi * ((psimax - 190) / (240 - 190))))
    else:
        n = 0

    # Ripple height (eta) equation
    eta = ahat * meta * n * (0.275 - ((0.022 * psimax)**0.42))

    # Ripple length (lambda) equation
    lambda_ = ahat * mlambda * n * (1.97 - ((0.44 * psimax)**0.21))

    # Factor mu for fine sand adjustment
    if d50 <= 0.15:
        mu = 6
    elif d50 <= 0.2:
        mu = 6 - ((5 * (d50 - 0.15)) / (0.2 - 0.15))
    else:
        mu = 1

    # Current-related bed roughness
    rough = d50 * (mu + 6 * (shields_aa - 1))
    ksdelta = max((3 * d90), rough) + ((0.4 * eta**2) / lambda_)

    # Wave-related bed roughness
    ksw = max(d50, rough) + ((0.4 * eta**2) / lambda_)

    return (mu, ksdelta, ksw)
