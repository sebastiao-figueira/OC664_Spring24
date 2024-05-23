'''
Testing Phase Lag

By: Luis D. Perez-Squeo
Updated: May 22, 2024
'''
#%% PREPARING PYTHON ENVIRONMENT

# Importing Python packages
import scipy
import numpy as np
import VDA_functions as vda

#%% Important paths
data_path  = r'D:\Github\OC664_Spring24\A - Jacob'

#%% Loading data
mat = scipy.io.loadmat(data_path + '\A_intrawave_velocity_timeseries_output_file.mat')

T_c        = mat['T_c'][0]
T_cu       = mat['T_cu'][0]
T_t        = mat['T_t'][0]
T_tu       = mat['T_tu'][0]
u_hat_c    = mat['u_hat_c'][0]
u_hat_t    = mat['u_hat_t'][0]
c_w        = mat['c_w'][0]

a_hat      = mat['a_hat'][0]
T          = mat['T'][0]
u_crx      = mat['u_crx'][0]
u_delta    = mat['u_delta'][0]
u_hat      = mat['u_hat'][0]
u_tilda_cr = mat['u_tilda_cr'][0]
u_tilda_tr = mat['u_tilda_tr'][0]
u_trx      = mat['u_trx'][0]
u_w        = mat['u_w'][0]
u_x        = mat['u_x'][0]

#%% Constants 

g          = 9.81            # Gravitational acceleration [m/s^2]
nu         = 2e-6            # Kinematic viscosity of water [m^2/s]
d50        = 0.2/1000        # Median grain diameter [m]
d90        = 0.00025         # [mm] value taken from: "Uncertainty in Nearshore Sand Bar Migration" which had d90 of 0.00024. We rounded up slightly given that our d50 is also larger than that paper's d50.
delta      = 0.2             # [m]
rho        = 1025            # Seawater density [kg/m^3]
rho_s      = 2650            # Sediment density [kg/m^3]
s          = rho_s / rho     # Specific gravity of the sediment (considering quartz)

#%%



# Colin/Emily output:
shields_aa, f_wc, f_wt, ksdelta, ksw, fdelta, f_w, eta = vda.modelfunctions.combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat,
                                                                              u_delta, delta, d50, d90, s, g)
# Maggie output:
fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx = vda.modelfunctions.currentfric(u_delta, u_hat, f_wc, f_wt, rho, f_w,
                                                                                   c_w, u_crx, u_trx, s, g, d50, fdelta)

# Carly output:
dstar      = vda.modelfunctions.dimensionless_grainsize(d50)
shields_cr = vda.modelfunctions.critical_shields(dstar)
omega_c    = vda.modelfunctions.sandload_crest(theta_cx, shields_cr)
omega_t    = vda.modelfunctions.sandload_trough(theta_tx, shields_cr)

# Carson output:
delta_sc = vda.modelfunctions.sfl_thickness(theta_cmag, d50)
delta_st = vda.modelfunctions.sfl_thickness(theta_tmag, d50)

# Luis output:
def phaseLag(s, d50, eta, u_hat_c, u_hat_t, c_w, T_c, T_cu, delta_sc, T_t, T_tu, delta_st, omega_c, omega_t, alpha=8.2, xi=1.7):
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


#%% Outputs

# Tell Carly to change these outputs as numpy arrays
omega_c = np.array(omega_c)
omega_t = np.array(omega_t)
delta_sc = delta_sc[0]
delta_st = delta_st[0]

omega_cc = []
omega_ct = []
omega_tt = []
omega_tc = []

for i in range(len(eta)):
    omega_cc_temp, omega_ct_temp, omega_tt_temp, omega_tc_temp = phaseLag(s, d50, eta[i], u_hat_c[i], u_hat_t[i], c_w[i], T_c[i], T_cu[i], delta_sc, T_t[i], T_tu[i], delta_st, omega_c[i], omega_t[i], alpha=8.2, xi=1.7)

    omega_cc.append(omega_cc_temp)
    omega_ct.append(omega_ct_temp)
    omega_tt.append(omega_tt_temp)
    omega_tc.append(omega_tc_temp)

omega_cc = np.array(omega_cc)
omega_ct = np.array(omega_ct)
omega_tt = np.array(omega_tt)
omega_tc = np.array(omega_tc)