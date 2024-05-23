import numpy as np
import scipy
import VDA_functions as vda


# TEST MODEL HERE
# %% general inputs

delta = 0.2  # m
d50 = 0.0002  # mm; this was given in our midterm so we've followed the same assumptions
d90 = 0.00025  # mm; value taken from: "Uncertainty in Nearshore Sand Bar Migration" which had d90 of 0.00024. We rounded up slightly given that our d50 is also larger than that paper's d50.
g = 9.81
# check this
rho = 1025
rho_s = 2650
# check this
s = rho_s / rho

# Load Jacob's output
mat = scipy.io.loadmat('A_intrawave_velocity_timeseries_output_file.mat')

a_hat = mat['a_hat'][0]
c_w = mat['c_w'][0]
T = mat['T'][0]
T_c = mat['T_c'][0]
T_cu = mat['T_cu'][0]
T_t = mat['T_t'][0]
T_tu = mat['T_tu'][0]
u_crx = mat['u_crx'][0]
u_delta = mat['u_delta'][0]
u_hat = mat['u_hat'][0]
u_hat_c = mat['u_hat_c'][0]
u_hat_t = mat['u_hat_t'][0]
u_tilda_cr = mat['u_tilda_cr'][0]
u_tilda_tr = mat['u_tilda_tr'][0]
u_trx = mat['u_trx'][0]
u_w = mat['u_w'][0]
u_x = mat['u_x'][0]

# Colin/Emily output:
shields_aa, f_wc, f_wt, ksdelta, ksw, fdelta, f_w, eta = vda.modelfunctions.combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat,
                                                                              u_delta, delta, d50, d90, s, g)
# Maggie output:
fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx = vda.modelfunctions.currentfric(u_delta, u_hat, f_wc, f_wt, rho, f_w,
                                                                                   c_w, u_crx, u_trx, s, g, d50, fdelta)

# Carly output:
dstar = vda.modelfunctions.dimensionless_grainsize(d50)
shields_cr = vda.modelfunctions.critical_shields(dstar)
omega_c = vda.modelfunctions.sandload_crest(theta_cx, shields_cr)
omega_t = vda.modelfunctions.sandload_trough(theta_tx, shields_cr)

# Carson output:
sheetflow_thickness_c = vda.modelfunctions.sfl_thickness(theta_cmag, d50)
sheetflow_thickness_t = vda.modelfunctions.sfl_thickness(theta_tmag, d50)

# Luis output:
omega_cc, omega_ct, omega_tt, omega_tc = vda.modelfunctions.phaseLag(s, d50, eta, u_hat_c, u_hat_t, c_w, T_c, T_cu, sheetflow_thickness_c, T_t, T_tu, sheetflow_thickness_t, omega_c, omega_t, alpha=8.2, xi=1.7, g=9.81, nu=2e-6)