################################
# This is the combined functions for wave friction factor and ripple model
# It was necessary to combine the two because an iterative calculation takes place
# This code reflects the work of both Emily An and Colin Arnowil. Their individual contributions are more fully commented and documented in separate files.
################################

def combined_wavefric_ripples(u, T_cu, T_c, T_tu, T_t, ahat, u_hat, u_deltamag, delta, d50, d90, s, g):
    
    import numpy as np
    
    # Emily section, with the iterative part moved down below
    
    # Maximum mobility number
    psimax = (max(u))**2 / ((s - 1) * g * d50)

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

    
    # Iterative part
    epsilon = 0.000000001 # tolerance to determine convergence
    new_calc = 1 # arbitrary initial guess
    last_calc = 0 # arbitrary arbitrary starting point for comparison
    ksdelta = 0 # arbitrary starting point for comparison

    while abs(new_calc - last_calc) >= epsilon: # testing convergence
        last_calc = new_calc
        ksdelta = ksdelta + 0.000001 # change to make iteration go
        fdelta = 2 * ((0.4 / ( np.log( 30 * delta / ksdelta )))**2) # from Maggie section; current friction factor for boundary layer
        shields_aa = ((1/2)*fdelta*(u_deltamag**2))/((s-1)*g*d50) + ((1/4)*f_w*(u_hat**2))/((s-1)*g*d50) # time averaged absolute shields stress
        rough = d50 * (mu + 6 * (shields_aa - 1)) # current-related bed roughness
        ksdelta = max((3 * d90), rough) + ((0.4 * eta**2) / lambda_) # also current-related bed roughness
        new_calc = ksdelta
    
    
    # Wave-related bed roughness
    ksw = max(d50, rough) + ((0.4 * eta**2) / lambda_)
    
    
    
    ### Colin section 
    
    c1 = 2.6 # noted in VDA13, empirically derived
    
    # Swart Equation, which goes into the time-averaged absolute Shields stress below
    if ahat/ksw > 1.587:
        f_w = 0.00251*np.exp(5.21*((ahat/ksw)**(-0.19))) # equation A.4; the Swart Equation. In the absence of acceleration skewness f_wc/f_wt below reduce to this
    else:
        f_w = 0.3    
    
    # Wave friction factor, separated out into crest and trough half cycles to account for acceleration skewness
    if ahat/ksw > 1.587:
        f_wc = 0.00251*np.exp(5.21*((((2*T_cu/T_c)**c1)*ahat)/ksw)**(-0.19)) # equation 21 for crest
        f_wt = 0.00251*np.exp(5.21*((((2*T_tu/T_t)**c1)*ahat)/ksw)**(-0.19)) # equation 21 for trough
    else:
       f_wc = 0.3
       f_wt = 0.3
     
    return shields_aa, f_wc, f_wt, ksdelta, ksw, mu