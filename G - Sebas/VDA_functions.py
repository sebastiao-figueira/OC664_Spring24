import scipy
import numpy as np


class modelfunctions:

    # %% Colin and Emily's function
    @staticmethod
    def combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat, u_delta, delta, d50, d90, s, g):
        # Maximum mobility number
        psimax = (1.27 * u_hat) ** 2 / ((s - 1) * g * d50)

        # diamteter in mm for mlambda and meta
        d50_mm = 1000 * d50
        
        # Equation B.4 in VDA13
        if d50_mm <= 0.22:
            mlambda = 0.73
        elif 0.22 < d50_mm < 0.3:
            mlambda = 0.73 + (0.27 * (d50_mm - 0.22)) / (0.3 - 0.22)
        else:
            mlambda = 1

        # Equation B.3 in VDA13
        if d50_mm <= 0.22:
            meta = 0.55
        elif 0.22 < d50_mm < 0.3:
            meta = 0.55 + (0.45 * (d50_mm - 0.22)) / (0.3 - 0.22)
        else:
            meta = 1

        # Equation A.2 from VDA13
        if d50_mm <= 0.15:
            mu = 6
        elif 0.15 < d50_mm < 0.2:
            mu = 6 - (5 * (d50_mm - 0.15)) / (0.2 - 0.15)
        else:
            mu = 1
        
        # Equation B.5 in VDA13
        n = np.ones(np.shape(psimax)) # preallocating empty array
            
        for i in range(0, len(psimax)):
            if psimax[i] <= 190:
                n[i] = 1
            elif 190 < psimax[i] < 240:
                n[i] = 0.5 * (1 + np.cos(np.pi * (psimax[i] - 190) / 50))
            else:
                n[i] = 0
    
        # Equation B.1 in VDA13
        eta = np.ones(np.shape(psimax)) # preallocating empty array
            
        for i in range(0, len(psimax)):
            eta[i] = meta * n[i] * a_hat[i] * (0.275 - 0.022 * psimax[i]**0.42)
            # eta[i] = 0
    
        # Equation B.2 in VDA13
        lambda_ = np.ones(np.shape(psimax)) # preallocating empty array
    
        for i in range(0, len(psimax)): # Ripple length (lambda) equation
            lambda_[i] = mlambda * a_hat[i] * n[i] * (1.97 - 0.44 * psimax[i]**0.21)
            # lambda_[i] = 0
            
        c1 = 2.6 # noted in VDA13, empirically derived        
    
        ksdelta = np.ones(np.shape(psimax))*0.5 # preallocating empty arrays; used 0.5 so as not to interfere with dummy variables below
        fdelta = np.ones(np.shape(psimax))*0.5
        shields_aa = np.ones(np.shape(psimax))*0.5
        rough = np.ones(np.shape(psimax))*0.5
        ksw = np.ones(np.shape(psimax))*0.5
        f_w = np.ones(np.shape(psimax))*0.5
    
        # Iterative part
        epsilon = 0.000000001 # tolerance to determine convergence
        
    
        for i in range(0, len(psimax)):
            new_calc = 1 # arbitrary initial guess
            last_calc = 0 # arbitrary arbitrary starting point for comparison

            while abs(new_calc - last_calc) >= epsilon: # testing convergence
                last_calc = new_calc
                
                ksdelta[i] = ksdelta[i] + 0.000001 # change to make iteration go
                # Equation 20 from VDA13
                fdelta[i] = 2 * (0.4 / np.log(30 * delta / ksdelta[i]))**2 # from Maggie section; current friction factor for boundary layer
                # Equation A.3 from VDA13
                shields_aa[i] = 0.5 * np.multiply(fdelta[i], u_delta[i]**2) / ((s - 1) * g * d50) + 0.25 * np.multiply(f_w[i], u_hat[i]**2) / ((s - 1) * g * d50) # time averaged absolute shields stress
                # Part of equation A.1 and A.2 from VDA13
                rough[i] = d50 * (mu + 6 * (shields_aa[i] - 1)) # current-related bed roughness
                # Equation A.5 from VDA13
                if lambda_[i] == 0:
                    ksw[i] = max(d50, rough[i])
        
                else:
                    ksw[i] = max(d50, rough[i]) + (0.4 * eta[i]**2) / lambda_[i] # Wave-related bed roughness
                    
                # Equation A.4 from VDA13
                if np.divide(a_hat[i],ksw[i]) > 1.587: # Swart Equation, which goes into the time-averaged absolute Shields stress
                    f_w[i] = 0.00251*np.exp(5.21*((np.divide(a_hat[i],ksw[i]))**(-0.19))) # equation A.4; the Swart Equation. In the absence of acceleration skewness f_wc/f_wt below reduce to this
                    # f_w[i] = 0.007
                else:
                    f_w[i] = 0.3
                    # f_w[i] = 0.007
                # Equation A.1 from VDA13
                if lambda_[i] == 0:
                    ksdelta[i] = max((3 * d90), rough[i])
                else:
                    ksdelta[i] = max((3 * d90), rough[i]) + (0.4 * eta[i]**2) / lambda_[i] # also current-related bed roughness
                new_calc = ksdelta[i]
            
    
        f_wc = np.ones(np.shape(psimax)) # preallocating empty array
        f_wt = np.ones(np.shape(psimax))
    
            # Wave friction factor, separated out into crest and trough half cycles to account for acceleration skewness
            
        # Equation 21 from VDA13
        for i in range(0, len(psimax)):
            if a_hat[i]/ksw[i] > 1.587:
                f_wc[i] = 0.00251*np.exp(5.21 * (a_hat[i]*((2*T_cu[i] / T_c[i])**c1)/ksw[i])**(-0.19))
                f_wt[i] = 0.00251*np.exp(5.21 * (a_hat[i]*((2*T_tu[i] / T_t[i])**c1)/ksw[i])**(-0.19))
            else:
                f_wc[i] = 0.3
                f_wt[i] = 0.3

        return f_wc, f_wt, fdelta, f_w, eta, shields_aa, ksdelta, ksw, lambda_

    # %% Maggie's function
    @staticmethod
    def currentfric(u_delta, uhat, fw_c, fw_t, rho, fw_swart, cw, uc_r, ut_r, s, g, d50, fdelta):
        #assumed wave boundary layer thickness (see section 6 VDA13)
        
        #compare current vs wave orbital velocity(eqn 19)
        alpha = np.abs(u_delta)/(np.abs(u_delta) + uhat)  #[-]
        
        #wave-current friction factor at crest (fwdelt_c) and trough (fwdelt_t) (eqn 18)
        fwdelt_c = alpha * fdelta + (1 - alpha) * fw_c  #[-]
        fwdelt_t = alpha * fdelta + (1 - alpha) * fw_t  #[-]
        
        #wave Reynolds stress TwRe (eqn 22)
        alpha_w = 0.424  #dimensionless scale factor
        fw_Delt = alpha * fdelta + (1 - alpha) * fw_swart  #[-] full-cycle wave-current friction factor
        TwRe    = fw_Delt * alpha_w * ((uhat)**3)/(2*cw)
        
        #Shields magnitude at crest (theta_cmag) and trough (theta_tmag) (eqn 17)
        theta_cmag = 0.5 * fwdelt_c * ((np.abs(uc_r)) ** 2) / ((s - 1) * g * d50)
        theta_tmag = 0.5 * fwdelt_t * ((np.abs(ut_r)) ** 2) / ((s - 1) * g * d50)
        
        #Bed shear stress (Shields vector) at crest (theta_cx) and trough (theta_tx) (eqn 15)
        theta_cx = theta_cmag * uc_r / np.abs(uc_r) + TwRe / ((s - 1) * g * d50)
        theta_tx = theta_tmag * ut_r / np.abs(ut_r) + TwRe / ((s - 1) * g * d50)

        return fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx, fw_Delt

    # %% Carly's functions
    # D_star (Soulsby 1997)
    @staticmethod
    def dimensionless_grainsize(d: float = 0.0002, g: float = 9.81, rho_s: float = 2650.0,
                                rho: float = 1025.0, nu: float = 2e-6):
        """Calculates dimensionless grain size, D*.
        Parameters
        ----------
        d : float
            Grain diameter (default is 0.2mm, 0.0002m)
        g : float
            Acceleration due to gravity (default is 9.81)
        rho_s : float
            Density of sediment (default is 2650.0)
        rho : float
            Density of water (default is 1025.0)
        nu : float
            Kinematic viscosity of water (default is 2e-6)

        Returns
        -------
        d_star : float
            Dimensionless grain size
        """

        d_star = d * (g * ((rho_s / rho) - 1) / nu ** 2) ** (1 / 3)
        return d_star

    # Critical Shields (Soulsby 1997)
    @staticmethod
    def critical_shields(d_star):
        """Calculates critical Shields number.
        Parameters
        ----------
        d_star : float
            Dimensionless grain size

        Returns
        -------
        theta_cr : float
            Critical Shields number
        """

        theta_cr = 0.30 / (1 + 1.2 * d_star) + 0.055 * (1 - np.exp(-0.020 * d_star))
        return theta_cr

    # Calculate sand load for half-cycles: omega_c or omega_t
    # If shields_c and shields_t are actually arrays, will need to loop through each value,
    # but this is the correct logic
    @staticmethod
    def sandload(theta_c, theta_t, theta_cr, m: float = 11.0, n: float = 1.2):
        """Calculates sand load for crest half-cycle.
        Parameters
        ----------
        theta_c
            Shields Number for crest half cycle
        theta_t
            Shields Number for trough half cycle
        theta_cr : float
            Critical Shields number
        m : float
            Proportionality constant (default is 11.0 from van der A et al. 2013)
        n : float
            Power constant (default is 1.2 from van der A et al. 2013)
        Returns
        -------
        theta_c : np.array
            Array of sand load under crest half cycle
        theta_t : np.array
            Array of sand load under trough half cycle
        """

        omega_c = []
        for val in theta_c:
            if val > theta_cr:
                omega_c_single = m * (val - theta_cr) ** n
            else:
                omega_c_single = 0
            omega_c.append(omega_c_single)
        omega_c = np.array(omega_c)

        omega_t = []
        for val in theta_t:
            if val > theta_cr:
                omega_t_single = m * (val - theta_cr) ** n
            else:
                omega_t_single = 0
            omega_t.append(omega_t_single)
        omega_t = np.array(omega_t)

        return omega_c, omega_t

    # %% Carson's function
    @staticmethod
    def shields_hat(fwdelt_c, fwdelt_t, u_hat_c, u_hat_t, rho_s=2065, rho=1025, d50=0.0002, g=9.81):
        """Calculates Shields parameter based on crest/trough velocity amplitude.
        Parameters
        ----------
        fwdelt_c
            Wave-current friction factor - crest
        fwdelt_t
            Wave-current friction factor - trough
        u_hat_c
            Crest velocity amplitude
        u_hat_t
            Trough velocity amplitude
        rho_s
            Sediment density (default 2650)
        rho
            Water density (default 1025)
        g
            Acceleration due to gravity (default 9.81)
        d50
            Median grain size (default 0.0002 m)
        Returns
        -------
        theta_hat_c : np.array
            Shields parameter - crest
        theta_hat_t : np.array
            Shields parameter - trough
        """
        
        # calculate s parameter from VDA13
        s = rho_s/rho
        
        # Equation C.2 from VDA13
        theta_hat_c = (0.5 * fwdelt_c * u_hat_c ** 2) / ((s - 1) * g * d50)
        theta_hat_t = (0.5 * fwdelt_t * u_hat_t ** 2) / ((s - 1) * g * d50)
        return theta_hat_c, theta_hat_t

    @staticmethod
    def sfl_thickness(shield, d50):
        '''
        function for Sheet Flow Layer Thickness
        *dependant on grain size*

        inputs:
        shield: shields parameter --> dimensionless (FROM MAGGIE EQ.17)
        d50: median grain size [m] (NEEDS TO BE IN METERS)

        output:
        Sheet flow layer thickness [m] as an array (FOR LUIS EQ.27 and EQ.28)
        '''

        d50 = d50 * 1000  # conversion from m to mm for calculation

        # Equation C.1 from VDA13
        sfl = []
        if d50 <= 0.15:
            sfl.append(25 * shield)
        elif 0.15 < d50 < 0.20:
            sfl.append(25 - ((12 * (d50 - 0.15)) / (0.20 - 0.15)))
        else:
            sfl.append(13 * shield)

        return np.array(sfl).flatten() * (d50 / 1000)  # conversion back to m

    # %% Luis' function
    @staticmethod
    def phaseLagSingle(rho, rho_s, d50, eta, u_hat_c, u_hat_t, c_w, T_c, T_cu, delta_sc, T_t, T_tu, delta_st, omega_c,
                       omega_t, alpha=8.2, xi=1.7, g=9.81, nu=2e-6):
        '''
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
        '''

        # Calculate specific gravity
        s = rho_s / rho

        # Non-dimensional grain size
        D_star = ((g * (s - 1) / (nu ** 2)) ** (1 / 3)) * d50

        # Settling velocity [m/s] (following van Rijn (1984))
        if D_star ** 3 <= 16.187:
            w_s = nu * (D_star ** 3) / (18 * d50)
        elif 16.187 < D_star ** 3 <= 16187:
            w_s = (10 * nu / d50) * (np.sqrt(1 + 0.01 * (D_star ** 3)) - 1)
        elif D_star ** 3 > 16187:
            w_s = 1.1 * nu * (D_star ** 1.5) / d50

        # Phase lag parameters for the crest half cycle (P_c) and trough half cycle (P_t)
        if eta > 0:  # Ripple regime
            P_c = alpha * (1 - (xi * u_hat_c / c_w)) * (eta / (2 * (T_c - T_cu) * w_s))
            P_t = alpha * (1 + (xi * u_hat_t / c_w)) * (eta / (2 * (T_t - T_tu) * w_s))

        elif eta == 0:  # Sheet flow regime
            P_c = alpha * (1 - (xi * u_hat_c / c_w)) * (delta_sc / (2 * (T_c - T_cu) * w_s))
            P_t = alpha * (1 + (xi * u_hat_t / c_w)) * (delta_st / (2 * (T_t - T_tu) * w_s))
        
        # Sand load entrainment during wave cycles
        if P_c <= 1:
            omega_cc = omega_c
            omega_ct = 0
        elif P_c > 1:
            omega_cc = omega_c / P_c
            omega_ct = (1 - 1 / P_c) * omega_c

        if P_t <= 1:
            omega_tt = omega_t
            omega_tc = 0
        elif P_t > 1:
            omega_tt = omega_t / P_t
            omega_tc = (1 - 1 / P_t) * omega_t
            
        print(P_c)
        print(P_t)

        return omega_cc, omega_ct, omega_tt, omega_tc

    @staticmethod
    def phaseLag(rho, rho_s, d50, eta, u_hat_c, u_hat_t, c_w, T_c, T_cu, delta_sc, T_t, T_tu, delta_st, omega_c,
                 omega_t, alpha=8.2, xi=1.7, g=9.81, nu=2e-6):
        '''
        Same as phaseLag above but works with vectors!
        '''

        omega_cc = []
        omega_ct = []
        omega_tt = []
        omega_tc = []

        for i in range(len(eta)):
            omega_cc_temp, omega_ct_temp, omega_tt_temp, omega_tc_temp = modelfunctions.phaseLagSingle(rho, rho_s, d50,
                                                                                                       eta[i],
                                                                                                       u_hat_c[i],
                                                                                                       u_hat_t[i],
                                                                                                       c_w[i], T_c[i],
                                                                                                       T_cu[i],
                                                                                                       delta_sc[i],
                                                                                                       T_t[i], T_tu[i],
                                                                                                       delta_st[i],
                                                                                                       omega_c[i],
                                                                                                       omega_t[i],
                                                                                                       alpha=8.2,
                                                                                                       xi=1.7)

            omega_cc.append(omega_cc_temp)
            omega_ct.append(omega_ct_temp)
            omega_tt.append(omega_tt_temp)
            omega_tc.append(omega_tc_temp)

        omega_cc = np.array(omega_cc)
        omega_ct = np.array(omega_ct)
        omega_tt = np.array(omega_tt)
        omega_tc = np.array(omega_tc)

        return omega_cc, omega_ct, omega_tt, omega_tc

    @staticmethod
    def sediment_transport(omega, wave_period, shields, rho=1000, rho_s=2650, d_50=0.0002, g=9.81):
        """
        Calculate the net transport rate induced by non-breaking waves and currents
        based on the method set forth by van der A (2013).
    
        Args:
            omega (array): Sand load entrainment and transport rates during different wave phases.
            wave_period (array): Durations of wave phases.
            shields (array): Magnitudes of Shields parameters and bed shear stresses.
            rho (float): Water density (default: 1000 kg/m^3).
            rho_s (float): Sediment density (default: 2650 kg/m^3).
            d_50 (float): Median sediment grain size (default: 0.0002 m).
            g (float): Acceleration due to gravity (default: 9.81 m/s^2).
    
        Returns:
            float: Net sediment transport rate.
        """

        # Check input array sizes
        if len(omega) != 4 or len(wave_period) != 5 or len(shields) != 4:
            raise ValueError("Input arrays must have sizes [4], [5], and [4] respectively.")

        # Calculate s parameter
        s = rho_s/rho

        # Extract values from input arrays
        omega_cc, omega_ct, omega_tt, omega_tc = omega
        T, T_c, T_cu, T_t, T_tu = wave_period
        shields_c, shields_t, shields_cx, shields_tx = shields

        # Calculate transport rates
        q_c = np.sqrt(shields_c) * T_c * (omega_cc + T_c / (2 * T_cu) * omega_tc) * shields_cx / shields_c
        q_t = np.sqrt(shields_t) * T_t * (omega_tt + T_t / (2 * T_tu) * omega_ct) * shields_tx / shields_t

        # Calculate net sediment transport rate
        q_s = (q_c + q_t) * np.sqrt((s - 1) * g * d_50 ** 3) / T
        
        # Represenative statistic to report (could be np.mean as well)
        Q_sum = np.median(q_s)
        Q_sum = Q_sum * 86400 # [m^2/day]

        return q_s, Q_sum
