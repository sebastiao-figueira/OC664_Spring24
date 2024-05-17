import numpy as np

# Mobile-bed effects in oscillatory flow sheet flow
# Dohmen-Janssen 2001
# Code by Carson Williams - Git: CarsonW503

def sfl_thickness(shield, d50):
    
    '''
    function for Sheet Flow Layer Thickness
    *dependant on grain size*
    
    inputs:
    shield: shields parameter --> dimensionless (FROM MAGGIE EQ.17)
    d50: median grain size [m] (NEEDS TO BE IN METERS)
    
    output:
    Sheet flow layer thickness [m] as an array
    '''
    
    d50 = d50*1000      #conversion from m to mm for calculation
    
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