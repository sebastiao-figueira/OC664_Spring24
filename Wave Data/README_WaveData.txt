Read-Me for "A_intrawave_velocity_timeseries_matlab_v3.m"

This matlab code takes data collected from October 22, 1997 from Duck, NC and extracts various intra-wave data.  Some of the intra wave data is applied to various equations from VanderA2013.  



Inputs:

dataset_19971022_0700EST.mat - this file contains almost 3 hours worth of data for 8 different cross shore locations. It includes:

u 	= cross-shore velocity (m/s, 2Hz) 
p 	= water surface elevation (m, 2Hz) relative to zp
x 	= cross-shore position of sensor
y 	= alongshore position of sensor
zu 	= vertical position for u-sensor
zp 	= vertical position for p-sensor
zbed 	= vertical position of seabed (from sonar altimeter)

The code applies only data from the first gauge, which is column 1 of the data file, dataset_19971022_0700EST.mat.



Outputs:

T 	= intra-wave period in s, 1x1299
T_c	= intra-wave duration of the crest (positive) half cycle in s, 1x1299 matrix
T_t	= intra-wave duration of the trough (negative) half cycle in s, 1x1299 matrix
T_cu	= intra-wave duration of accelerating flow within the crest half cycle in s, 1x1299 matrix
T_tu	= intra-wave duration of accelerating flow within the trough half cycle in s, 1x1299 matrix

	  Note: for the above 5 output parameters, refer to Fig. 1 in VanderA2013.

u_w	= time varying freestream orbital velocities in m/s, 1x1299 cell 
	  (each cell inside has varying number of rows x 1 column)
u_delta = input for eqns (4), (A.3), and (19), steady state current velocity in m/s, 1x1299 matrix
u_x	= eqn (4) output, time varying orbital velocity in x-direction in m/s, 1x1299 cell 
	  (each cell inside has varying number of rows x 1 column)
u_hat	= eqn (8) output, representative orbital velocity amplitude in m/s, 1x1299
u_hat_c	= eqn (8a) output, representative orbital velocity amplutude (with T_c as an input) in m/s, 1x1299
u_hat_t	= eqn (8b) output, representative orbital velocity amplutude (with T_t as an input) in m/s, 1x1299
a_hat	= eqn (9) output, representative orbital excursion amplitude in m, 1x1299
u_tilda_cr = eqn (10) output, representative half-cycle orbital velocity for the wave crest in m/s, 1x1299
u_tilda_tr = eqn (11) output, representative half-cycle orbital velocity for the wave trough in m/s, 1x1299 
u_crx	= eqn (12) output, representative combined wave current velocity vector for crest halfcycle in x-direction in m/s, 1x1299
u_trx	= eqn (13) output, representative combined wave current velocity vector for trough halfcycle in x-direction in m/s, 1x1299
c_w	= celerity for each intra-wave in m/s, 1x1299
















