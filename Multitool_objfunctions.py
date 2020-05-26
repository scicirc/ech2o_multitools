#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''

Python scripts for ECH2O

-------
Routine: Subroutines objective functions
-------
Author: S. Kuppel
Created on 05/2020
-------------------------------------------------
'''

from __future__ import division
from scipy.special import gamma
import numpy as np
import pandas as pd
import math
# import sys
import spotpy_forked.spotpy as spotpy
import sys
# ==============================================================================
# -- Multi-objective function(s)


def MultiObj(obs, sim, Data, Opti, w=True):

    like = 0
    w_tot = 0
    Ltot = []
    
    for i in range(Data.nobs):

        oname = Data.names[i]

        # Have obervation and simulations matching the same time period
        # obs: already pre-processed
        tobs = pd.to_datetime(obs[oname]['Date'].values)
        o = np.asanyarray(obs[oname]['value'].values)
        # sim: trim sim to obs timespan
        # + only keep dates with obs (even if nan)
        # print(sim)
        # print(sim.shape)
        if Data.nobs == 1:
            s = np.asanyarray([sim[j] for j in range(Data.lsimEff) if
                              Data.simt[j] in tobs])
        else:
            s = np.asanyarray([sim[i][j] for j in range(Data.lsimEff) if
                               Data.simt[j] in tobs])

        # Remove nan
        tmp = s*o
        s = np.asanyarray([s[k] for k in range(len(tmp)) if not
                           np.isnan(tmp[k])])
        o = np.asanyarray([o[j] for j in range(len(tmp)) if not
                           np.isnan(tmp[j])])

        # For GWD, we center using mean
        if oname.split('_')[0] in ['GWD', 'WTD', 'GWL', 'WTL']:
            # print(oname, 'remove mean')
            # s -= np.mean(s)
            # o -= np.mean(o)
            # Use RMSE for water table
            # L = -1*RMSE(o, s)
            L = - spotpy.objectivefunctions.rmse(o, s) / np.mean(o)
        else:
            # MAE for streamflow
            # L = -1*MAE(o, s)
            L = - spotpy.objectivefunctions.mae(o, s) / np.mean(o)

        Ltot += [L]
        # Now use your favorite metrics for each obs type
        # for now, Schoups & Vrugt
        # L = SchoupsVrugt_GL(o, s, oname, max_B=0.5)

        # Weight by the number of obs ?
        if w is True:
            like += len(s) * L
            w_tot += len(s)
        else:
            like += L
        # print(oname, mean_o, n, d, 1-n/d)
    if w is True:
        like /= w_tot
    print(np.round(Ltot, 2))
    # print('OJ: ', like)
    return like
# ==============================================================================
# -- Likelihood: Mean absolute error


def MAE(obs, sim):

    # Normalize
    sd = np.std(obs)
    obs /= sd
    sim /= sd
    # MAE
    mae = np.mean(np.abs(obs-sim))

    return mae
# ==============================================================================
# -- Likelihood: RMSE


def RMSE(obs, sim):

    # Normalize
    sd = np.std(obs)
    obs /= sd
    sim /= sd
    # Rmse
    rmse = np.sqrt(np.mean((obs-sim)**2))

    return rmse

# ==============================================================================
# -- Likelihood: Nash-Sutcliffe


def NashSutcliffe(obs, sim):

    # And finally...
    mean_o = np.mean(obs)  # , axis=axis)
    n = ((obs - sim) ** 2).sum()  # axis=axis)
    d = ((obs - mean_o) ** 2).sum()  # axis=axis)

    return 1 - n/d
# ==============================================================================
# -- Generalized likelihood from Schoups & Vrught, 2010


def SchoupsVrugt_GL(obs, sim, oname, max_B):

    nobs = len(obs)  # Equal to length of sim (see preprocessing)

    # Normalize the data
    sd = np.std(obs)
    obs /= sd
    sim /= sd

    base_xi = 1
    base_beta = 0.5
    base_sig0 = 0.5
    base_sig1 = 0
    base_phi = 0.5*max_B
    sd_mult = 0.1

    # Algorithm based on Johnson
    for i in range(3000):

        # Skewness parameter xi
        xi = np.random.normal(base_xi, sd_mult*(10 - 0.1))
        # Reflexive boundaries
        if xi > 10:
            xi = 10 + (10-xi)
        if xi < 0.1:
            xi = 0.1 + abs(xi-0.1)

        # Kurtosis parameter beta
        beta = np.random.normal(base_beta, sd_mult*(1 - 0))
        if beta > 1:
            beta = 1 + (1-beta)
        if beta < 0:
            beta = abs(beta)

        # Residual error: fixed part
        sig0 = np.random.normal(base_sig0, sd_mult*(1 - 0))
        if sig0 > 1:
            sig0 = 1 + (1-sig0)
        if sig0 < 0:
            sig0 = 0 + abs(sig0)
        # Residual error: heteroscedastic part
        if oname.split('_')[0] in ['GWD', 'WTD', 'GWL', 'WTL']:
            # Assumed homoscedastic for water table depth
            sig1 = 0
        else:
            sig1 = np.random.normal(base_sig1, sd_mult*(2 - 0))
            if sig1 > 2:
                sig1 = 2 + (2-sig1)
            if sig1 < 0:
                sig1 = abs(sig1)

        # Autocorrelation coefficient
        phi = np.random.normal(base_phi, sd_mult*(max_B - 0))
        if phi > max_B:
            phi = max_B + (max_B - phi)
        if phi < 0:
            phi = abs(phi)

        # Remove autocorrelation from residuals
        Res = sim - obs
        Res_Norm = np.zeros(len(Res))
        Res_Norm_Hom = np.zeros(len(Res))

        for c in range(len(Res)):
            Res_Norm[c] = Res[c] - phi*Res[c-1]

        # Remove heteroscedasticity
        Res_Norm_Hom[1] = 0
        sigma_t = np.zeros((len(Res_Norm)))
        for c in range(len(Res_Norm)):
            sigma_t[c] = sig0 + sig1*abs(sim[c])
            Res_Norm_Hom[c] = Res_Norm[c]/sigma_t[c]

        # Compute Likelihood
        M1 = gamma(1 + beta)/(math.sqrt(gamma(3*(1+beta)/2)) *
                              math.sqrt(gamma((1+beta)/2)))
        M2 = 1

        mu_xi = M1*(xi - 1/xi)
        sig_xi = math.sqrt((M2 - M1**2)*(xi**2 + xi**-2) +
                           2*(M1**2) - M2)

        kurto1 = gamma(3*(1 + beta)/2)
        kurto2 = gamma((1 + beta)/2)
        wb = math.sqrt(kurto1)/((1+beta)*(kurto2**(3/2)))
        cb = (kurto1/kurto2)**(1/(1 + beta))

        alpha_wt = (mu_xi + sig_xi*Res_Norm_Hom) * \
            (xi**np.sign(mu_xi + sig_xi*Res_Norm_Hom))

        ln_sigma_t = sigma_t
        for c in range(len(sigma_t)):
            ln_sigma_t[c] = math.log(sigma_t[c])

        abs_alpha_wt = alpha_wt
        for c in range(len(alpha_wt)):
            abs_alpha_wt[c] = abs(alpha_wt[c])

        new_L = nobs*math.log((2*sig_xi*wb)/(xi + xi**-1)) - \
            sum(ln_sigma_t) - cb*sum(abs_alpha_wt**(2/(1 + beta)))
        if np.isnan(new_L):
            new_L = -100000

        # Initialize
        if i == 0:
            max_L = new_L

        # Update if it yields a higher log-likelihood
        if new_L > max_L:
            max_L = new_L
            base_xi = xi
            base_beta = beta
            base_sig0 = sig0
            base_sig1 = sig1
            base_phi = phi

    # print(oname, 'xi:', np.round(xi, 3), 'beta:', np.round(beta, 3),
    #       'phi:', np.round(phi, 3),
    #       'sig0', np.round(sig0, 3), 'sig1', np.round(sig1, 3))
    return max_L
    # Results[0] = max_L
    # Results[1] = base_xi
    # Results[2] = base_beta
    # Results[3] = base_sig0
    # Results[4] = base_sig1
    # Results[5] = phi

    # return Results
