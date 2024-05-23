The file current_friction_factor_v2.py contains a script to calculate wave-current friction factor and related outputs
Eqn numbers are in reference to VDA13

#INPUTS(collaborators' names in brackets)
    #udelta=current velocity (time-mean of velocity time series) [JACOB]
    #uhat=representative orbital velocity amplitude [JACOB] (eqn 8)
    #ksdelta=current-related roughness [EMILY] (eqn A.1)
    #fw_c=wave friction factor at crest [COLIN] (eqn 21)
    #fw_t=wave friction factor at trough [COLIN] (eqn 21)
    #rho=density of seawater in kg/m^3 [SEBAS]
    #cw=wave speed where cw=L/T, L=wavelength, T=wave period, from LWT [LUIS]
    #fw_swart=Swart's friction factor from appendix A.4 [COLIN] (eqn A.4)
    #uc_r=representative half-cycle wave-current velocity at crest [JACOB] (eqn 12)
    #ut_r=representative half-cycle wave-current velocity at trough [JACOB] (eqn 13)
    #s=specific gravity [SEBAS] (section 2 of VDA13)
    #g=constant of gravitational acceleration [SEBAS]
    #d50=median grain size of sediment [SEBAS]
	
#OUTPUTS
	#fwdelt_c=wave friction factor at crest (eqn 18)
	#fwdelt_t=wave friction factor at trough (eqn 18)
	#TwRe=wave Reynolds stress (eqn 22)
	#theta_cmag=magnitude of the Shields parameter at crest (eqn 17)
	#theta_tmag=magnitude of the Shields parameter at trough (eqn 17)
	#theta_cx=x component of the Shields parameter at crest (eqn 15)
	#theta_tx=x component of the Shields parameter at trough (eqn 15)
	

	
