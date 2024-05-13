'''
Phase lag effect variables for the Sediment Transport group term project

By: Luis D. Perez-Squeo
Updated: May 13, 2024
'''

def D_star():
    
    '''
    Inputs:
    eta         -> ripple height    
    u_c         -> peak crest orbital velocity    
    c_w         -> wave speed
    T_c         -> period of the crest (positive) half cycle
    T_cu        -> period of accelerating flow within the crest half cycle
    delta_sc    -> sheet flow layer thickness for the crest half cycle
    w_s         -> particle settling velocity
    alpha       -> calibration coefficient (alpha=8.2 in VDA13)
    xi          -> calibration factor (xi=1.7 in VDA13)    
    
    Output:
    P_c         -> phase lag parameter for the crest half cycle
    '''
    
    Dstar    = ((g*(s-1)/(nu**2))**(1/3)) * d50             # Non-dimensional grain size 
    
    return w_s

def w_s():
    
    '''
    Inputs:
    eta         -> ripple height    
    u_c         -> peak crest orbital velocity    
    c_w         -> wave speed
    T_c         -> period of the crest (positive) half cycle
    T_cu        -> period of accelerating flow within the crest half cycle
    delta_sc    -> sheet flow layer thickness for the crest half cycle
    w_s         -> particle settling velocity
    alpha       -> calibration coefficient (alpha=8.2 in VDA13)
    xi          -> calibration factor (xi=1.7 in VDA13)    
    
    Output:
    P_c         -> phase lag parameter for the crest half cycle
    '''
    
    # Settling velocity [m/s] 
    if Dstar**3 <= 16.187:
        ws = nu*(Dstar**3)/(18*d50)
    elif 16.187 < Dstar**3 <= 16187: 
        ws = (10*nu/d50)*(np.sqrt(1 + 0.01*(Dstar**3))-1)
    elif Dstar**3 > 16187:
        ws = 1.1*nu*(Dstar**1.5)/d50
    
    return w_s

def P_c(eta, u_c, c_w, T_c, T_cu, delta_sc, w_s, alpha=8.2, xi=1.7):
    
    '''
    Inputs:
    eta         -> ripple height    
    u_c         -> peak crest orbital velocity    
    c_w         -> wave speed
    T_c         -> period of the crest (positive) half cycle
    T_cu        -> period of accelerating flow within the crest half cycle
    delta_sc    -> sheet flow layer thickness for the crest half cycle
    w_s         -> particle settling velocity
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
    eta         -> ripple height         
    u_t         -> peak trough orbital velocity    
    c_w         -> wave speed
    T_t         -> period of the trough (negative) half cycle
    T_tu        -> period of accelerating flow within the trough half cycle
    delta_st    -> sheet flow layer thickness for the trough half cycle
    w_s         -> particle settling velocity
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