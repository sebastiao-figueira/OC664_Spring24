# -*- coding: utf-8 -*-
"""
Created on Fri May 10 16:46:17 2024

@author: libbym
"""

import numpy as np

#function to calculate wave-current friction factor and related outputs
#eqn numbers are in reference to VDA13

#INPUTS(collaborators' names in brackets)
    #udelta_mag=magnitude of current velocity udelta (time-mean of velocity time series) [JACOB]
    #uhat=representative orbital velocity amplitude [JACOB] (eqn 8)
    #ksdelta=current-related roughness [EMILY] (eqn A.1)
    #fw_c=wave friction factor at crest [COLIN] (eqn 21)
    #fw_t=wave friction factor at trough [COLIN] (eqn 21)
    #rho=density of seawater in kg/m^3 [SEBAS]
    #cw=wave speed where cw=L/T, L=wavelength, T=wave period, from LWT [LUIS]
    #fw_Delt=full-cycle wave friction factor [COLIN] (eqn A.4)
    #uc_r=representative half-cycle wave-current velocity at crest [JACOB] (eqn 12)
    #ut_r=representative half-cycle wave-current velocity at trough [JACOB] (eqn 13)
    #s=specific gravity [SEBAS] (section 2 of VDA13)
    #g=constant of gravitational acceleration [SEBAS]
    #d50=median grain size of sediment [SEBAS]

def currentfric(u_deltamag,uhat,ksdelta,fw_c,fw_t,rho,fw_Delt,cw,uc_r,ut_r,s,g,d50):
    #assumed wave boundary layer thickness (see section 6 VDA13)
    delta=0.2 #[m]
    #compare current vs wave orbital velocity(eqn 19)
    alpha=np.divide(u_deltamag,(u_deltamag+uhat)) #[-]
    #current-related friction factor (eqn 20)
    fdelta=2*(np.divide(0.4,(np.log(np.divide(30*delta,ksdelta)))))**2 #[-]
    #wave-current friction factor at crest (fwdelt_c) and trough (fwdelt_t) (eqn 18)
    fwdelt_c=np.multiply(alpha,fdelta)+np.multiply((1-alpha),fw_c) #[-]
    fwdelt_t=np.multiply(alpha,fdelta)+np.multiply((1-alpha),fw_t) #[-]
    #wave Reynolds stress TwRe (eqn 22)
    alpha_w=0.424 #dimensionless scale factor
    TwRe=np.multiply(np.divide((rho*fw_Delt),2*cw),alpha_w*uhat**3)
    #Shields magnitude at crest (theta_cmag) and trough (theta_tmag) (eqn 17)
    theta_cmag=np.divide(np.multiply(0.5*fwdelt_c,(np.abs(uc_r))**2),((s-1)*g*d50))
    theta_tmag=np.divide(np.multiply(0.5*fwdelt_t,(np.abs(ut_r))**2),((s-1)*g*d50))
    #Bed shear stress (Shields vector) (eqn 15)
    theta_cx=np.multiply(theta_cmag,np.divide(uc_r,np.abs(uc_r)))+np.divide(TwRe,((s-1)*g*d50))
    theta_tx=np.multiply(theta_cmag,np.divide(uc_r,np.abs(uc_r)))+np.divide(TwRe,((s-1)*g*d50))
    
    return fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx
    
    