'''
Phase lag effect equations for the Sediment Transport group term project

By: Luis D. Perez-Squeo
Updated: May 13, 2024
'''

import numpy as np

def D_star(s, d50):
    
    '''
    Inputs:
    s           -> specific gravity of the sediment (considering quartz) 
    d50         -> median grain diameter [m]
    g           -> gravitational acceleration [m/s^2]   
    nu          -> kinematic viscosity of water [m^2/s]

    Output:
    D_star      -> non-dimensional grain size 
    '''
    
    g      = 9.81
    nu     = 2e-6
    D_star = ((g*(s-1)/(nu**2))**(1/3)) * d50
    
    return D_star

def w_s(D_star, d50):
    
    '''
    Inputs:
    D_star      -> non-dimensional grain size 
    d50         -> median grain diameter [m]  
    
    Output:
    w_s         -> settling velocity [m/s]
    '''
    
    nu = 2e-6
    
    if D_star**3 <= 16.187:
        w_s = nu*(D_star**3)/(18*d50)
    elif 16.187 < D_star**3 <= 16187: 
        w_s = (10*nu/d50)*(np.sqrt(1 + 0.01*(D_star**3))-1)
    elif D_star**3 > 16187:
        w_s = 1.1*nu*(D_star**1.5)/d50
    
    return w_s

def P_c(eta, u_c, c_w, T_c, T_cu, delta_sc, w_s, alpha=8.2, xi=1.7):
    
    '''
    Inputs:
    eta         -> ripple height [m]    
    u_c         -> peak crest orbital velocity [m/s] 
    c_w         -> wave speed [m/s]
    T_c         -> period of the crest half cycle [s]
    T_cu        -> period of accelerating flow within the crest half cycle [s]
    delta_sc    -> sheet flow layer thickness for the crest half cycle [m]
    w_s         -> particle settling velocity [m/s]
    alpha       -> calibration coefficient (alpha=8.2 in VDA13)
    xi          -> calibration factor (xi=1.7 in VDA13)    
    
    Output:
    P_c         -> phase lag parameter for the crest half cycle
    '''
    
    if eta > 0:    P_c = alpha * ((1-xi*u_c)/c_w) * (eta/(2*(T_c-T_cu)*w_s))
    elif eta == 0: P_c = alpha * ((1-xi*u_c)/c_w) * (delta_sc/(2*(T_c-T_cu)*w_s))
    
    return P_c

def P_t(eta, u_t, c_w, T_t, T_tu, delta_st, w_s, alpha=8.2, xi=1.7):
    
    '''
    Inputs:
    eta         -> ripple height [m]       
    u_t         -> peak trough orbital velocity [m/s]
    c_w         -> wave speed [m/s]
    T_t         -> period of the trough half cycle [s]
    T_tu        -> period of accelerating flow within the trough half cycle [s]
    delta_st    -> sheet flow layer thickness for the trough half cycle [m]
    w_s         -> particle settling velocity [m/s]
    alpha       -> calibration coefficient (alpha=8.2 in VDA13)
    xi          -> calibration factor (xi=1.7 in VDA13)    
    
    Output:
    P_t         -> phase lag parameter for the trough half cycle
    '''
    
    if eta > 0:    P_t = alpha * ((1-xi*u_t)/c_w) * (eta/(2*(T_t-T_tu)*w_s))
    elif eta == 0: P_t = alpha * ((1-xi*u_t)/c_w) * (delta_st/(2*(T_t-T_tu)*w_s))
    
    return P_t

def omega_cc(P_c, omega_c):
    
    '''
    Inputs:
    P_c         -> phase lag parameter for the crest half cycle    
    omega_c     -> sand load entrained during the wave crest period   
    
    Output:
    omega_cc    -> sand load entrained during the wave crest period and transported during the crest period
    '''
    
    if P_c <= 1:  omega_cc = omega_c
    elif P_c > 1: omega_cc = omega_c/P_c
    
    return omega_cc

def omega_ct(P_c, omega_c):
    
    '''
    Inputs:
    P_c         -> phase lag parameter for the crest half cycle    
    omega_c     -> sand load entrained during the wave crest period       
    
    Output:
    omega_ct    -> sand load entrained during the wave crest period and transported during the trough period
    '''
    
    if P_c <= 1:  omega_ct = 0
    elif P_c > 1: omega_ct = (1-1/P_c)*omega_c
    
    return omega_ct

def omega_tt(P_t, omega_t):
    
    '''
    Inputs:
    P_t         -> phase lag parameter for the trough half cycle    
    omega_t     -> sand load entrained during the wave trough period       
    
    Output:
    omega_tt    -> sand load entrained during the wave trough period and transported during the trough period
    '''
    
    if P_t <= 1:  omega_tt = omega_t
    elif P_t > 1: omega_tt = omega_t/P_t
    
    return omega_tt

def omega_tc(P_t, omega_t):
    
    '''
    Inputs:
    P_t         -> phase lag parameter for the trough half cycle    
    omega_t     -> sand load entrained during the wave trough period    
    
    Output:
    omega_tc    -> sand load entrained during the wave trough period and transported during the crest period
    '''
    
    if P_t <= 1:  omega_tc = 0
    elif P_t > 1: omega_tc = (1-1/P_t)*omega_t
    
    return omega_tc