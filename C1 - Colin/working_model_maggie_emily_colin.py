import scipy
import numpy as np

# %% general inputs

delta = 0.2 # m
d50 = 0.0002 # mm; this was given in our midterm so we've followed the same assumptions
d90 = 0.00025 # mm; value taken from: "Uncertainty in Nearshore Sand Bar Migration" which had d90 of 0.00024. We rounded up slightly given that our d50 is also larger than that paper's d50.
s = 2650/1000
g = 9.81
rho = 1025

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

# %% Colin and Emily's function

def combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat, u_delta, delta, d50, d90, s, g):

        # Maximum mobility number
    psimax = (1.27 * u_hat)**2 / ((s - 1) * g * d50)

        # Lambda multipliers m
    mlambda = 0.73 # given our d50 this value is correct

        # Eta multipliers
    meta = 0.55 # given our d50 this value is correct

        # Factor mu for fine sand adjustment
    mu = 1 # given our d50 this value is correct

    n = np.ones(np.shape(psimax)) # preallocating empty array
            
    for i in range(0, len(psimax)):
        if psimax[i] <= 190:
            n[i] = 1
        elif psimax[i] <= 240:
            n[i] = 0.5 * (1 + np.cos(np.pi * (psimax[i] - 190) / 50))
        elif psimax[i] > 240:
            n[i] = 0

    eta = np.ones(np.shape(psimax)) # preallocating empty array
        
    for i in range(0, len(psimax)):
        eta[i] = meta * n[i] * a_hat[i] * (0.275 - (0.022 * (psimax[i]**0.42)))

    lambda_ = np.ones(np.shape(psimax)) # preallocating empty array

    for i in range(0, len(psimax)): # Ripple length (lambda) equation
        lambda_[i] = mlambda * a_hat[i] * n[i] * (1.97 - 0.44 * (psimax[i]**0.21))
        
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
            fdelta[i] = 2 * ((0.4 / ( np.log(np.divide((30 * delta), ksdelta[i]))))**2) # from Maggie section; current friction factor for boundary layer
            shields_aa[i] = ((1/2)*np.multiply(fdelta[i],(u_delta[i]**2)))/((s-1)*g*d50) + ((1/4)*np.multiply(f_w[i],(u_hat[i]**2))/((s-1)*g*d50)) # time averaged absolute shields stress
            rough[i] = d50 * (mu + 6 * (shields_aa[i] - 1)) # current-related bed roughness
            if lambda_[i] == 0:
                ksw[i] = max(d50, rough[i])
            else:
                ksw[i] = max(d50, rough[i]) + ((0.4 * eta[i]**2) / lambda_[i]) # Wave-related bed roughness
            if np.divide(a_hat[i],ksw[i]) > 1.587: # Swart Equation, which goes into the time-averaged absolute Shields stress
                f_w[i] = 0.00251*np.exp(5.21*((np.divide(a_hat[i],ksw[i]))**(-0.19))) # equation A.4; the Swart Equation. In the absence of acceleration skewness f_wc/f_wt below reduce to this
            else:
                f_w[i] = 0.3
            if lambda_[i] == 0:
                ksdelta[i] = max((3 * d90), rough[i])
            else:
                ksdelta[i] = max((3 * d90), rough[i]) + np.divide((0.4 * eta[i]*2), lambda_[i]) # also current-related bed roughness
            new_calc = ksdelta[i]
        

    f_wc = np.ones(np.shape(psimax)) # preallocating empty array
    f_wt = np.ones(np.shape(psimax))

        # Wave friction factor, separated out into crest and trough half cycles to account for acceleration skewness
        
    for i in range(0, len(psimax)):
        if a_hat[i]/ksw[i] > 1.587:
            f_wc[i] = 0.00251*np.exp(5.21 * ((((2*T_cu[i] / T_c[i])**c1)*a_hat[i])/ksw[i])**(-0.19))
            f_wt[i] = 0.00251*np.exp(5.21 * ((((2*T_tu[i] / T_t[i])**c1)*a_hat[i])/ksw[i])**(-0.19))
        else:
            f_wc[i] = 0.3
            f_wt[i] = 0.3
            
            
    return shields_aa, f_wc, f_wt, ksdelta, ksw, fdelta, f_w

shields_aa, f_wc, f_wt, ksdelta, ksw, fdelta, f_w = combined_wavefric_ripples(T_cu, T_c, T_tu, T_t, a_hat, u_hat, u_delta, delta, d50, d90, s, g)

# %% Maggie's function

def currentfric(u_delta,uhat,fw_c,fw_t,rho,fw_swart,cw,uc_r,ut_r,s,g,d50,fdelta):
    #assumed wave boundary layer thickness (see section 6 VDA13)
    #compare current vs wave orbital velocity(eqn 19)
    alpha=np.divide(np.abs(u_delta),(np.abs(u_delta)+uhat)) #[-]
    #wave-current friction factor at crest (fwdelt_c) and trough (fwdelt_t) (eqn 18)
    fwdelt_c=np.multiply(alpha,fdelta)+np.multiply((1-alpha),fw_c) #[-]
    fwdelt_t=np.multiply(alpha,fdelta)+np.multiply((1-alpha),fw_t) #[-]
    #wave Reynolds stress TwRe (eqn 22)
    alpha_w=0.424 #dimensionless scale factor
    fw_Delt=np.multiply(alpha,fdelta)+np.multiply((1-alpha),fw_swart) #[-] full-cycle wave-current friction factor
    TwRe=np.multiply(np.divide((rho*fw_Delt),(2*cw)),(alpha_w*uhat**3))
    #Shields magnitude at crest (theta_cmag) and trough (theta_tmag) (eqn 17)
    theta_cmag=np.divide(np.multiply((0.5*fwdelt_c),((np.abs(uc_r))**2)),((s-1)*g*d50))
    theta_tmag=np.divide(np.multiply((0.5*fwdelt_t),((np.abs(ut_r))**2)),((s-1)*g*d50))
    #Bed shear stress (Shields vector) at crest (theta_cx) and trough (theta_tx) (eqn 15)
    theta_cx=np.multiply(theta_cmag,np.divide(uc_r,np.abs(uc_r)))+np.divide(TwRe,((s-1)*g*d50))
    theta_tx=np.multiply(theta_tmag,np.divide(ut_r,np.abs(ut_r)))+np.divide(TwRe,((s-1)*g*d50))
    
    return fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx

fwdelt_c, fwdelt_t, TwRe, theta_cmag, theta_tmag, theta_cx, theta_tx = currentfric(u_delta,u_hat,f_wc,f_wt,rho,f_w,c_w,u_crx,u_trx,s,g,d50,fdelta)
# %%
