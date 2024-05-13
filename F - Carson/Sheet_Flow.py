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
    d50: median grain size [mm]
    
    output:
    shields parameter
    '''
    
    shield = ((1/2)*fw*(ui**2)) / ((s-1)*g*d50)
    return shield

def sfl_thickness(shield, d50):
    
    '''
    function for Sheet Flow Layer Thickness
    *dependant on grain size*
    
    inputs:
    shield: shields parameter --> dimensionless
    d50: median grain size [mm]
    
    output:
    Sheet flow layer thickness
    '''
    
    sfl = []
    if d50 <= 0.00015:
        sfl.append(25*shield)
    elif 0.00015 <= d50 <= 0.00020:
        sfl.append(25 - ())