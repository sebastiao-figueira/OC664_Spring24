#%% PREPARING PYTHON ENVIRONMENT

# Importing Python packages

import numpy as np
import scipy.io
import VDA_functions as vda
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Font settings for figures
font_label = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 20,
        }

font = {'family' : 'serif',
'weight' : 'bold',
'size' : 16}
plt.rc('font', **font)

SMALL_SIZE = 18
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('axes', titlesize=MEDIUM_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)  
plt.rc('xtick', labelsize=SMALL_SIZE)   
plt.rc('ytick', labelsize=SMALL_SIZE)    
plt.rc('legend', fontsize=SMALL_SIZE)   
plt.rc('figure', titlesize=BIGGER_SIZE)

#%% GENERAL INPUTS

delta = 0.2     # m
d50   = 0.0002  # m; this was given in our midterm so we've followed the same assumptions
d90   = 0.00025 # m; value taken from: "Uncertainty in Nearshore Sand Bar Migration" which had d90 of 0.00024. We rounded up slightly given that our d50 is also larger than that paper's d50.
g     = 9.81
rho   = 1025
rho_s = 2650
s     = rho_s/rho

#%% LOADING DATA

# A - intra wave data extraction
# Load Jacob's output from SandyDuck '97 data
mat = scipy.io.loadmat('A_intrawave_velocity_timeseries_output_file.mat')

# Representative Wave
# a_hat      = mat['a_hat'][0][4:5]
# c_w        = mat['c_w'][0][4:5]
# T          = mat['T'][0][4:5]
# T_c        = mat['T_c'][0][4:5]
# T_cu       = mat['T_cu'][0][4:5]
# T_t        = mat['T_t'][0][4:5]
# T_tu       = mat['T_tu'][0][4:5]
# u_crx      = mat['u_crx'][0][4:5]
# u_delta    = mat['u_delta'][0][4:5]
# u_hat      = mat['u_hat'][0][4:5]
# u_hat_c    = mat['u_hat_c'][0][4:5]
# u_hat_t    = mat['u_hat_t'][0][4:5]
# u_tilda_cr = mat['u_tilda_cr'][0][4:5]
# u_tilda_tr = mat['u_tilda_tr'][0][4:5]
# u_trx      = mat['u_trx'][0][4:5]
# u_w        = mat['u_w'][0][4:5]
# u_x        = mat['u_x'][0][4:5]

# Wave-by-wave
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

#%% APPLYING VDA13 EQS

# C1 and D - Wave friction factor and ripple model
# Colin/Emily output:
f_wc, f_wt, fdelta, f_w, eta, shields_aa, ksdelta, ksw, lambda_ = vda.modelfunctions.combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat,
                                                                            u_delta, delta, d50, d90, s, g)

# C2 - Current friction factor
# Maggie output:
fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx, fw_Delt = vda.modelfunctions.currentfric(u_delta, u_hat,
                                                                                                      f_wc, f_wt, rho,
                                                                                                      f_w,
                                                                                                      c_w, u_crx, u_trx,
                                                                                                      s, g, d50, fdelta)
# B - Entrained Sand Load
# Carly output:
dstar            = vda.modelfunctions.dimensionless_grainsize()
shields_cr       = vda.modelfunctions.critical_shields(dstar)
omega_c, omega_t = vda.modelfunctions.sandload(theta_cmag, theta_tmag, shields_cr)

# F - Sheet flow thickness
# Carson output:
shields_hat_c, shields_hat_t = vda.modelfunctions.shields_hat(fwdelt_c, fwdelt_t, u_hat_c, u_hat_t)
sheetflow_thickness_c        = vda.modelfunctions.sfl_thickness(shields_hat_c, d50)
sheetflow_thickness_t        = vda.modelfunctions.sfl_thickness(shields_hat_t, d50)

# E - Phase Lag
# Luis output:
omega_cc, omega_ct, omega_tt, omega_tc, P_c, P_t = vda.modelfunctions.phaseLag(rho, rho_s, d50, eta, u_hat_c, u_hat_t, c_w, T_c,
                                                                     T_cu, sheetflow_thickness_c, T_t, T_tu,
                                                                     sheetflow_thickness_t, omega_c, omega_t, alpha=8.2, xi=1.7, g=9.81, nu=2e-6)

# G - Main Function
# Sebas output
omega       = [omega_cc, omega_ct, omega_tt, omega_tc]
wave_period = [T, T_c, T_cu, T_t, T_tu]
shields     = [theta_cmag, theta_tmag, theta_cx, theta_tx]

q_s, Q_sum = vda.modelfunctions.sediment_transport(omega, wave_period, shields, rho, rho_s, d50, g)

# q_s   [m^2/s]
# Q_sum [m^2/day]

#%% q_s histograms

# Set up the figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 10))

# First histogram
step1 = 1e-4
bins1 = np.arange(0, max(q_s) + step1, step1)
ax1.hist(q_s, bins=bins1, density=False, color='r', edgecolor='black', zorder=2)
ax1.set_xlabel('Net Transport Rate (m$^2$/s)', fontdict=font_label)
ax1.set_ylabel('Number of Waves', fontdict=font_label)
ax1.grid()
ax1.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax1.set_xlim(0, 5e-3)

# Second histogram
step2 = 1e-5
bins2 = np.arange(0, max(q_s) + step2, step2)
ax2.hist(q_s, bins=bins2, density=False, color='b', edgecolor='black', zorder=2)
ax2.set_xlabel('Net Transport Rate (m$^2$/s)', fontdict=font_label)
ax2.set_ylabel('Number of Waves', fontdict=font_label)
ax2.grid()
ax2.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax2.set_xlim(0, 1e-3)
ax2.set_ylim(0, 249)

# Save the combined figure
fig_name = r'qs_Histogram_Combined.png'
plt.savefig(fig_name, dpi=600, bbox_inches='tight')
plt.show()