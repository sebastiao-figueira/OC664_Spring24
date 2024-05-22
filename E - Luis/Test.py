'''
Testing Phase Lag

By: Luis D. Perez-Squeo
Updated: May 22, 2024
'''
#%% PREPARING PYTHON ENVIRONMENT

# Importing Python packages
import scipy

#%% Important paths
data_path  = r'D:\Github\OC664_Spring24\A - Jacob'

#%% Loading data
mat = scipy.io.loadmat(data_path + '\A_intrawave_velocity_timeseries_output_file.mat')

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