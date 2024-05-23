##########################
# Wave friction factor code for OC 664 - Sediment Transport group project
# Author: Colin Arnowil
# Source material: "Practical sand transport formula for non-breaking waves and currents". Van der A (2013); "VDA13"
# Equation/page #'s are in reference to this article
##########################

# My individual contribution alculates the wave related friction factor and the time-averaged absolute Shields stress. By necessity in the complete product it was combined with code from Emily An and Maggie Libby, for ripples and current related friction factor respectively. This readme reflects only my separate portion:

# necessary inputs:

# T_cu, T_c, T_tu, T_t -- varying crest vs trough half cycle periods due to acceleration skewness [Jacob]
# a_hat [Jacob - Equation 9]
# u_hat -- a representative orbital velocity [Jacob - Equation 8]
# u_deltamag -- magnitude of current velocity udelta (time-mean of velocity time series) [Jacob]
# k_sw [Emily - Equation A.5]
# f_delta -- current-related friction factor [Maggie - Equation 20]
# d50 -- a median grain size for the sediment [Sebas]
# s -- specific gravity of sediment [Sebas]
# g -- gravity [let's just guess this one]

# Equation 21 contains an exponent, c1, which should be c1 = 2.6
# "Optimisation of c1 against the measurements of bed shear stress by Van der A et al. (2011) for a range of acceleration-skewed oscillatory flows resulted in c1 = 2.6." VDA13 p. 29

Rough code blocks:

c1 = 2.6 # noted above
    
# Equation A.4: Swart Equation, which goes into the time-averaged absolute Shields stress below. In the absence of acceleration skewness f_wc/f_wt below reduces to this.
if a_hat/k_sw > 1.587:
    f_w = 0.00251*np.exp(5.21*((a_hat/k_sw)**(-0.19)))
else:
    f_w = 0.3
    

# Time-Averaged absolute Shields Stress
theta_mag_timeavg = ((1/2)*f_delta*(u_deltamag**2))/((s-1)*g*d50) + ((1/4)*f_w*(u_hat**2))/((s-1)*g*d50) # equation A.3
###### Note: the time-averaged absolute shields stress depends on f_delta, which depends on k_sdelta, which depends on the time-averaged absolute shields stress
# the above equation must necessarily be iterated with those others in order to arrive at a solution to a given desired accuracy
# this iteration (i.e. the integration of these equations) will be performed later in the full model
    
# Wave friction factor, separated out into crest and trough half cycles to account for acceleration skewness
if a_hat/k_sw > 1.587:
    f_wc = 0.00251*np.exp(5.21*((((2*T_cu/T_c)**c1)*a_hat)/k_sw)**(-0.19)) # equation 21 for crest
    f_wt = 0.00251*np.exp(5.21*((((2*T_tu/T_t)**c1)*a_hat)/k_sw)**(-0.19)) # equation 21 for trough
else:
   f_wc = 0.3
   f_wt = 0.3
     
return theta_mag_timeavg, f_wc, f_wt
