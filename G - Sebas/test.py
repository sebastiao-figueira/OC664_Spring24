import numpy as np
import scipy.io
import VDA_functions as vda
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# TEST MODEL HERE
# %% general inputs

delta = 0.2  # m
d50   = 0.0002  # m; this was given in our midterm so we've followed the same assumptions
d90   = 0.00025  # m; value taken from: "Uncertainty in Nearshore Sand Bar Migration" which had d90 of 0.00024. We rounded up slightly given that our d50 is also larger than that paper's d50.
g     = 9.81
rho   = 1025
rho_s = 2650
s     = rho_s/rho
# s = (rho_s - rho) / rho

# Load Jacob's output
mat = scipy.io.loadmat('A_intrawave_velocity_timeseries_output_file.mat')

a_hat      = mat['a_hat'][0]
c_w        = mat['c_w'][0]
T          = mat['T'][0]
T_c        = mat['T_c'][0]
T_cu       = mat['T_cu'][0]
T_t        = mat['T_t'][0]
T_tu       = mat['T_tu'][0]
u_crx      = mat['u_crx'][0]
u_delta    = mat['u_delta'][0]
u_hat      = mat['u_hat'][0]
u_hat_c    = mat['u_hat_c'][0]
u_hat_t    = mat['u_hat_t'][0]
u_tilda_cr = mat['u_tilda_cr'][0]
u_tilda_tr = mat['u_tilda_tr'][0]
u_trx      = mat['u_trx'][0]
u_w        = mat['u_w'][0]
u_x        = mat['u_x'][0]

# Colin/Emily output:
f_wc, f_wt, fdelta, f_w, eta, shields_aa, ksdelta, ksw, lambda_ = vda.modelfunctions.combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat,
                                                                            u_delta, delta, d50, d90, s, g)

# Maggie output:
fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx, fw_Delt = vda.modelfunctions.currentfric(u_delta, u_hat,
                                                                                                      f_wc, f_wt, rho,
                                                                                                      f_w,
                                                                                                      c_w, u_crx, u_trx,
                                                                                                      s, g, d50, fdelta)

# Carly output:
dstar            = vda.modelfunctions.dimensionless_grainsize()
shields_cr       = vda.modelfunctions.critical_shields(dstar)
omega_c, omega_t = vda.modelfunctions.sandload(theta_cmag, theta_tmag, shields_cr)

# Carson output:
shields_hat_c, shields_hat_t = vda.modelfunctions.shields_hat(fwdelt_c, fwdelt_t, u_hat_c, u_hat_t)
sheetflow_thickness_c        = vda.modelfunctions.sfl_thickness(shields_hat_c, d50)
sheetflow_thickness_t        = vda.modelfunctions.sfl_thickness(shields_hat_t, d50)

# Luis output:
omega_cc, omega_ct, omega_tt, omega_tc = vda.modelfunctions.phaseLag(rho, rho_s, d50, eta, u_hat_c, u_hat_t, c_w, T_c,
                                                                     T_cu, sheetflow_thickness_c, T_t, T_tu,
                                                                     sheetflow_thickness_t, omega_c, omega_t, alpha=8.2,
                                                                     xi=1.7, g=9.81, nu=2e-6)

# Sebas output
omega       = [omega_cc, omega_ct, omega_tt, omega_tc]
wave_period = [T, T_c, T_cu, T_t, T_tu]
shields     = [theta_cmag, theta_tmag, theta_cx, theta_tx]

q_s, Q_sum = vda.modelfunctions.sediment_transport(omega, wave_period, shields, rho, rho_s, d50, g)

print("for calculated theta avg q_s = ", np.mean(q_s))

# Define the bins with a width of 0.001
bins = np.arange(0, 0.1 + 0.001, 0.001)

# Create the histogram
plt.hist(q_s, bins=bins, color='crimson')

# Set the x and y axis labels
plt.xlabel('Net Transport Rate (m$^2$/s)')
plt.ylabel('Number of Waves')

# Set the x-axis limit
plt.xlim(0, 0.1)

# Display the plot
plt.show()

# %% general inputs
# Edited theta_cx and theta_tx
theta_cx = np.ones(np.shape(T)) 
theta_tx = np.ones(np.shape(T))

# Carly output:
dstar            = vda.modelfunctions.dimensionless_grainsize()
shields_cr       = vda.modelfunctions.critical_shields(dstar)
omega_c, omega_t = vda.modelfunctions.sandload(theta_cmag, theta_tmag, shields_cr)

# Carson output:
shields_hat_c, shields_hat_t = vda.modelfunctions.shields_hat(fwdelt_c, fwdelt_t, u_hat_c, u_hat_t)
sheetflow_thickness_c        = vda.modelfunctions.sfl_thickness(shields_hat_c, d50)
sheetflow_thickness_t        = vda.modelfunctions.sfl_thickness(shields_hat_t, d50)

# Luis output:
omega_cc, omega_ct, omega_tt, omega_tc = vda.modelfunctions.phaseLag(rho, rho_s, d50, eta, u_hat_c, u_hat_t, c_w, T_c,
                                                                     T_cu, sheetflow_thickness_c, T_t, T_tu,
                                                                     sheetflow_thickness_t, omega_c, omega_t, alpha=8.2,
                                                                     xi=1.7, g=9.81, nu=2e-6)

# Sebas output
omega       = [omega_cc, omega_ct, omega_tt, omega_tc]
wave_period = [T, T_c, T_cu, T_t, T_tu]
shields     = [theta_cmag, theta_tmag, theta_cx, theta_tx]

q_s_edited, Q_sum_edited = vda.modelfunctions.sediment_transport(omega, wave_period, shields, rho, rho_s, d50, g)

print("for constant theta avg q_s = ", np.mean(q_s_edited))

# Define the bins with a width of 0.000001
bins = np.arange(0, 0.00015 + 0.000001, 0.000001)

# Create the histogram with forest green color
plt.hist(q_s_edited, bins=bins, color='forestgreen')

# Set the x and y axis labels with superscript units
plt.xlabel('Net Transport Rate (m$^2$/s)')
plt.ylabel('Number of Waves')

# Set the x-axis limit
plt.xlim(0, 0.00015)

# Use scientific notation for the x-axis
plt.gca().xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# Display the plot
plt.show()