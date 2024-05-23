import scipy
import numpy as np


class modelfunctions:

    # %% Colin and Emily's function
    @staticmethod
    def combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat, u_delta, delta, d50, d90, s, g):
        # Maximum mobility number
        psimax = (1.27 * u_hat) ** 2 / ((s - 1) * g * d50)

        # Lambda multipliers m
        mlambda = 0.73  # given our d50 this value is correct

        # Eta multipliers
        meta = 0.55  # given our d50 this value is correct

        # Factor mu for fine sand adjustment
        mu = 1  # given our d50 this value is correct

        n = np.ones(np.shape(psimax))  # preallocating empty array

        for i in range(0, len(psimax)):
            if psimax[i] <= 190:
                n[i] = 1
            elif psimax[i] <= 240:
                n[i] = 0.5 * (1 + np.cos(np.pi * (psimax[i] - 190) / 50))
            elif psimax[i] > 240:
                n[i] = 0

        eta = np.ones(np.shape(psimax))  # preallocating empty array

        for i in range(0, len(psimax)):
            eta[i] = meta * n[i] * a_hat[i] * (0.275 - (0.022 * (psimax[i] ** 0.42)))

        lambda_ = np.ones(np.shape(psimax))  # preallocating empty array

        for i in range(0, len(psimax)):  # Ripple length (lambda) equation
            lambda_[i] = mlambda * a_hat[i] * n[i] * (1.97 - 0.44 * (psimax[i] ** 0.21))

        c1 = 2.6  # noted in VDA13, empirically derived

        ksdelta = np.ones(np.shape(
            psimax)) * 0.5  # preallocating empty arrays; used 0.5 so as not to interfere with dummy variables below
        fdelta = np.ones(np.shape(psimax)) * 0.5
        shields_aa = np.ones(np.shape(psimax)) * 0.5
        rough = np.ones(np.shape(psimax)) * 0.5
        ksw = np.ones(np.shape(psimax)) * 0.5
        f_w = np.ones(np.shape(psimax)) * 0.5

        # Iterative part
        epsilon = 0.000000001  # tolerance to determine convergence

        for i in range(0, len(psimax)):
            new_calc = 1  # arbitrary initial guess
            last_calc = 0  # arbitrary arbitrary starting point for comparison
            while abs(new_calc - last_calc) >= epsilon:  # testing convergence
                last_calc = new_calc
                ksdelta[i] = ksdelta[i] + 0.000001  # change to make iteration go
                fdelta[i] = 2 * ((0.4 / (np.log(np.divide((30 * delta), ksdelta[
                    i])))) ** 2)  # from Maggie section; current friction factor for boundary layer
                shields_aa[i] = ((1 / 2) * np.multiply(fdelta[i], (u_delta[i] ** 2))) / ((s - 1) * g * d50) + (
                        (1 / 4) * np.multiply(f_w[i], (u_hat[i] ** 2)) / (
                        (s - 1) * g * d50))  # time averaged absolute shields stress
                rough[i] = d50 * (mu + 6 * (shields_aa[i] - 1))  # current-related bed roughness
                if lambda_[i] == 0:
                    ksw[i] = max(d50, rough[i])
                else:
                    ksw[i] = max(d50, rough[i]) + ((0.4 * eta[i] ** 2) / lambda_[i])  # Wave-related bed roughness
                if np.divide(a_hat[i],
                             ksw[
                                 i]) > 1.587:  # Swart Equation, which goes into the time-averaged absolute Shields stress
                    f_w[i] = 0.00251 * np.exp(5.21 * ((np.divide(a_hat[i], ksw[i])) ** (
                        -0.19)))  # equation A.4; the Swart Equation. In the absence of acceleration skewness f_wc/f_wt below reduce to this
                else:
                    f_w[i] = 0.3
                if lambda_[i] == 0:
                    ksdelta[i] = max((3 * d90), rough[i])
                else:
                    ksdelta[i] = max((3 * d90), rough[i]) + np.divide((0.4 * eta[i] * 2),
                                                                      lambda_[i])  # also current-related bed roughness
                new_calc = ksdelta[i]

        f_wc = np.ones(np.shape(psimax))  # preallocating empty array
        f_wt = np.ones(np.shape(psimax))

        # Wave friction factor, separated out into crest and trough half cycles to account for acceleration skewness

        for i in range(0, len(psimax)):
            if a_hat[i] / ksw[i] > 1.587:
                f_wc[i] = 0.00251 * np.exp(5.21 * ((((2 * T_cu[i] / T_c[i]) ** c1) * a_hat[i]) / ksw[i]) ** (-0.19))
                f_wt[i] = 0.00251 * np.exp(5.21 * ((((2 * T_tu[i] / T_t[i]) ** c1) * a_hat[i]) / ksw[i]) ** (-0.19))
            else:
                f_wc[i] = 0.3
                f_wt[i] = 0.3

        return shields_aa, f_wc, f_wt, ksdelta, ksw, fdelta, f_w, eta

    # %% Maggie's function
    @staticmethod
    def currentfric(u_delta, uhat, fw_c, fw_t, rho, fw_swart, cw, uc_r, ut_r, s, g, d50, fdelta):
        #assumed wave boundary layer thickness (see section 6 VDA13)
        #compare current vs wave orbital velocity(eqn 19)
        alpha = np.divide(np.abs(u_delta), (np.abs(u_delta) + uhat))  #[-]
        #wave-current friction factor at crest (fwdelt_c) and trough (fwdelt_t) (eqn 18)
        fwdelt_c = np.multiply(alpha, fdelta) + np.multiply((1 - alpha), fw_c)  #[-]
        fwdelt_t = np.multiply(alpha, fdelta) + np.multiply((1 - alpha), fw_t)  #[-]
        #wave Reynolds stress TwRe (eqn 22)
        alpha_w = 0.424  #dimensionless scale factor
        fw_Delt = np.multiply(alpha, fdelta) + np.multiply((1 - alpha),
                                                           fw_swart)  #[-] full-cycle wave-current friction factor
        TwRe = np.multiply(np.divide((rho * fw_Delt), (2 * cw)), (alpha_w * uhat ** 3))
        #Shields magnitude at crest (theta_cmag) and trough (theta_tmag) (eqn 17)
        theta_cmag = np.divide(np.multiply((0.5 * fwdelt_c), ((np.abs(uc_r)) ** 2)), ((s - 1) * g * d50))
        theta_tmag = np.divide(np.multiply((0.5 * fwdelt_t), ((np.abs(ut_r)) ** 2)), ((s - 1) * g * d50))
        #Bed shear stress (Shields vector) at crest (theta_cx) and trough (theta_tx) (eqn 15)
        theta_cx = np.multiply(theta_cmag, np.divide(uc_r, np.abs(uc_r))) + np.divide(TwRe, ((s - 1) * g * d50))
        theta_tx = np.multiply(theta_tmag, np.divide(ut_r, np.abs(ut_r))) + np.divide(TwRe, ((s - 1) * g * d50))

        return fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx

    # %% Carly's functions
    # D_star (Soulsby 1997)
    @staticmethod
    def dimensionless_grainsize(d: float = 0.0002, g: float = 9.81, rho_s: float = 2650.0,
                                rho: float = 1000.0, nu: float = 0.00001):
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
            Density of water (default is 1000.0)
        nu : float
            Kinematic viscosity of water (default is 0.00001)

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

        theta_cr = 0.30 / (1 + 1.2 * d_star) + 0.055 * (1 - np.exp(-0.020) * d_star)
        return theta_cr

    # Calculate sand load for half-cycles: omega_c or omega_t
    # If shields_c and shields_t are actually arrays, will need to loop through each value,
    # but this is the correct logic
    @staticmethod
    def sandload_crest(theta_c, theta_cr, m: float = 11.0, n: float = 1.2):
        """Calculates sand load for crest half-cycle.
        Parameters
        ----------
        theta_c
            Shields Number for crest half cycle
        theta_cr : float
            Critical Shields number
        m : float
            Proportionality constant (default is 11.0 from van der A et al. 2013)
        n : float
            Power constant (default is 1.2 from van der A et al. 2013)
        Returns
        -------
        theta_cr : float
            Critical Shields number
        """

        omega_c = []
        for val in theta_c:
            if val > theta_cr:
                omega_c_single = m * (val - theta_cr) ** n
            else:
                omega_c_single = 0
            omega_c.append(omega_c_single)
        omega_c = np.array(omega_c)
        return omega_c

    @staticmethod
    def sandload_trough(theta_t, theta_cr, m: float = 11.0, n: float = 1.2):
        """Calculates sand load for crest half-cycle.
        Parameters
        ----------
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
        theta_cr : float
            Critical Shields number
        """
        omega_t = []
        for val in theta_t:
            if val > theta_cr:
                omega_t_single = m * (val - theta_cr) ** n
            else:
                omega_t_single = 0
            omega_t.append(omega_t_single)
        omega_t = np.array(omega_t)
        return omega_t

    # %% Carson's function
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

        sfl = []
        if d50 <= 0.15:
            sfl.append(25 * shield)
        elif 0.15 < d50 < 0.20:
            sfl.append(25 - (12 * (d50 - 0.15)) / (0.20 - 0.15))
        elif d50 >= 0.20:
            sfl.append(13 * shield)
        else:
            sfl.append(np.nan)

        return np.array(sfl).flatten() * (d50 / 1000)  # conversion back to m
    
    # %% Luis' function
    @staticmethod
    def phaseLagSingle(s, d50, eta, u_hat_c, u_hat_t, c_w, T_c, T_cu, delta_sc, T_t, T_tu, delta_st, omega_c, omega_t, alpha=8.2, xi=1.7, g=9.81, nu=2e-6):
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

        # Non-dimensional grain size
        D_star = ((g * (s - 1) / (nu**2))**(1 / 3)) * d50

        # Settling velocity [m/s] (following van Rijn (1984))
        if D_star**3 <= 16.187:
            w_s = nu * (D_star**3) / (18 * d50)
        elif 16.187 < D_star**3 <= 16187:
            w_s = (10 * nu / d50) * (np.sqrt(1 + 0.01 * (D_star**3)) - 1)
        elif D_star**3 > 16187:
            w_s = 1.1 * nu * (D_star**1.5) / d50

        # Phase lag parameters for the crest half cycle (P_c) and trough half cycle (P_t)
        if eta > 0:    # Ripple regime
            P_c = alpha * ((1 - xi * u_hat_c) / c_w) * (eta / (2 * (T_c - T_cu) * w_s))
            P_t = alpha * ((1 + xi * u_hat_t) / c_w) * (eta / (2 * (T_t - T_tu) * w_s))
            
        elif eta == 0: # Sheet flow regime
            P_c = alpha * ((1 - xi * u_hat_c) / c_w) * (delta_sc / (2 * (T_c - T_cu) * w_s))
            P_t = alpha * ((1 + xi * u_hat_t) / c_w) * (delta_st / (2 * (T_t - T_tu) * w_s))

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

        return omega_cc, omega_ct, omega_tt, omega_tc
    
    @staticmethod
    def phaseLag(s, d50, eta, u_hat_c, u_hat_t, c_w, T_c, T_cu, delta_sc, T_t, T_tu, delta_st, omega_c, omega_t, alpha=8.2, xi=1.7, g=9.81, nu=2e-6):
        '''
        Same as phaseLag above but works with vectors!
        '''
    
        omega_cc = []
        omega_ct = []
        omega_tt = []
        omega_tc = []
        
        for i in range(len(eta)):
            omega_cc_temp, omega_ct_temp, omega_tt_temp, omega_tc_temp = modelfunctions.phaseLagSingle(s, d50, eta[i], u_hat_c[i], u_hat_t[i], c_w[i], T_c[i], T_cu[i], delta_sc[i], T_t[i], T_tu[i], delta_st[i], omega_c[i], omega_t[i], alpha=8.2, xi=1.7)
        
            omega_cc.append(omega_cc_temp)
            omega_ct.append(omega_ct_temp)
            omega_tt.append(omega_tt_temp)
            omega_tc.append(omega_tc_temp)
        
        omega_cc = np.array(omega_cc)
        omega_ct = np.array(omega_ct)
        omega_tt = np.array(omega_tt)
        omega_tc = np.array(omega_tc)
        
        return omega_cc, omega_ct, omega_tt, omega_tc