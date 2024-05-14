import numpy as np

# Mobile-bed effects in oscillatory flow sheet flow
# Dohmen-Janssen 2001
# Code by Carson Williams - Git: CarsonW503

def shields(fw, ui, s, g, d50):
    
    '''
    function for maximum shields parameter in oscillatory 
    flow - dimensionless
    
    inputs:
    fw: wave-current friction factor --> dimensionless
    ui: crest trough velocity amplitude velocity
    s: relative density
    g: gravity acceleration
    d50: median grain size [m]
    
    output:
    shields parameter
    '''
    
    shield = ((1/2)*fw*(ui**2)) / ((s-1)*g*d50)
    return np.array(shield)

def sfl_thickness(shield, d50):
    
    '''
    function for Sheet Flow Layer Thickness
    *dependant on grain size*
    
    inputs:
    shield: shields parameter --> dimensionless
    d50: median grain size [m] (NEEDS TO BE IN METERS)
    
    output:
    Sheet flow layer thickness [m]
    '''
    
    d50 = d50*1000             #conversion from m to mm for calculation
    
    sfl = []
    if d50 <= 0.15:
        sfl.append(25*shield)
    elif 0.15 <= d50 <= 0.20:
        sfl.append(25 - (12*(d50-0.15))/(0.20-0.15))
    elif d50 >= 0.20:
        sfl.append(13*shield)
    else:
        sfl.append(np.nan)
        
    return np.array(sfl).flatten()*(d50/1000)    #conversion back to m