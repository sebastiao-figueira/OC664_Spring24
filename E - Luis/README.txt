The function 'phaseLagSingle' calculates phase lag related parameters for a single wave. The function 'phaseLag' uses 'phaseLagSingle' to calculate phase lag related parameters for more than one wave. 

        Inputs:
        s           -> specific gravity of the sediment (considering quartz)
        d50         -> median grain diameter [m]
        eta         -> ripple height [m]
        u_hat_c     -> peak crest orbital velocity [m/s]
        u_hat_t     -> peak trough orbital velocity [m/s]
        c_w         -> wave speed [m/s]
        T_c         -> period of the crest half cycle [s]
        T_cu        -> period of accelerating flow within the crest half cycle [s]
        delta_sc    -> sheet flow layer thickness for the crest half cycle [m]
        T_t         -> period of the trough half cycle [s]
        T_tu        -> period of accelerating flow within the trough half cycle [s]
        delta_st    -> sheet flow layer thickness for the trough half cycle [m]
        omega_c     -> sand load entrained during the wave crest period
        omega_t     -> sand load entrained during the wave trough period
        alpha       -> calibration coefficient (alpha=8.2 in VDA13)
        xi          -> calibration factor (xi=1.7 in VDA13)
        g           -> gravitational acceleration [m/s^2] (default 9.81)
        nu          -> kinematic viscosity of water [m^2/s] (default 2e-6)
        
        Outputs:
        omega_cc    -> sand load entrained during the wave crest period and transported during the crest period
        omega_ct    -> sand load entrained during the wave crest period and transported during the trough period
        omega_tt    -> sand load entrained during the wave trough period and transported during the trough period
	omega_tc    -> sand load entrained during the wave trough period and transported during the crest period	

	
