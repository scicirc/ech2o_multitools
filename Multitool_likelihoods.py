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
import math

# ==============================================================================
# -- global likelihood from Schoups & Vrught, 2010


def Multi_SchoupsVrugtGL(obs, sim, Data, Opti):

    Results = np.zeros(6)

    # remove flagged observed data (-999)
    rawcount = 0
    cleancount = 0
    SimObs_clean = []
    for i in SimObs:
        if SimObs[rawcount,1] > -999 and SimObs[rawcount,0] > -1000:
            SimObs_clean.append(SimObs[rawcount])
            cleancount = cleancount + 1
        rawcount = rawcount + 1

    SimObs_clean = np.asarray(SimObs_clean)  
    
    #normalize the data
    sd = np.std(SimObs_clean[:,1]);
    SimObs_clean[:,0] = SimObs_clean[:,0]/sd
    SimObs_clean[:,1] = SimObs_clean[:,1]/sd
    
    current_epsilon = 1
    current_beta = 0.5
    current_sigma_0 = 0.5
    current_sigma_1 = 0
    current_phi = 0.5*max_B
    sd_mult = 0.1
  
    for i in range(evals):
        epsilon = np.random.normal(current_epsilon, sd_mult*(10 - 0.5))
        beta = np.random.normal(current_beta, sd_mult*(1 - 0))
        sigma_0 = np.random.normal(current_sigma_0, sd_mult*(1 - 0))
        sigma_1 = np.random.normal(current_sigma_1, sd_mult*(2 - 0))
        phi = np.random.normal(current_phi, sd_mult*(max_B - 0))
      
        if epsilon > 10:
            epsilon = -1*epsilon + 20
        if epsilon < 0.5:
            epsilon = -1*epsilon + 1
        if beta > 1:
            beta = -1*beta + 2
        if beta < 0:
            beta = abs(beta)
        if sigma_0 > 1:
            sigma_0 = -1*sigma_0 + 2
        if sigma_0 < 0:
            sigma_0 = -1*sigma_0 + 0
        if sigma_1 > 2:
            sigma_1 = -1*sigma_1 + 4
        if sigma_1 < 0:
            sigma_1 = -1*sigma_1 + 0
        if phi > max_B:
            phi = -1*phi + 2*max_B
        if phi < 0:
            phi = -1*phi + 0
        
        #Remove autocorrelation from residuals
        Res = SimObs_clean[:,0] - SimObs_clean[:,1]
        Res_Norm = np.zeros(len(Res))
        Res_Norm_Hom = np.zeros(len(Res))
        
        for c in range(len(Res)):
            Res_Norm[c] = Res[c] - phi*Res[c-1]

        #Remove heteroscedasticity
        Res_Norm_Hom[1] = 0
        sigma_t = np.zeros((len(Res_Norm)))
        for c in range(len(Res_Norm)):
            sigma_t[c] = sigma_0 + sigma_1*abs(SimObs_clean[c,0])
            Res_Norm_Hom[c] = Res_Norm[c]/sigma_t[c]
        
        #Compute Likelihood
        n = len(SimObs_clean)
        M1 = gamma(1 + beta)/(math.sqrt(gamma(3*(1+beta)/2))*math.sqrt(gamma((1 + beta)/2)))
        M2 = 1

        mew_e = M1*(epsilon - epsilon**-1)
        sigma_e = math.sqrt((M2 - M1**2)*(epsilon**2 + epsilon**-2) + 2*(M1**2) - M2)

        A1_mcmc = gamma(3*(1 + beta)/2)
        A2_mcmc = gamma((1 + beta)/2)
        wb = math.sqrt(A1_mcmc)/((1+beta)*(A2_mcmc**(3/2)))
        cb = (A1_mcmc/A2_mcmc)**(1/(1 + beta))

        alpha_wt = (mew_e + sigma_e*Res_Norm_Hom)*(epsilon**np.sign(mew_e + sigma_e*Res_Norm_Hom))
        
        ln_sigma_t = sigma_t
        for c in range(len(sigma_t)):
            ln_sigma_t[c] = math.log(sigma_t[c])

        abs_alpha_wt = alpha_wt
        for c in range(len(alpha_wt)):
            abs_alpha_wt[c] = abs(alpha_wt[c])
        
        new_l = n*math.log((2*sigma_e*wb)/(epsilon + epsilon**-1)) - sum(ln_sigma_t) - cb*sum(abs_alpha_wt**(2/(1 + beta)))
       
        if math.isnan(new_l):
            new_l = -100000
      
        if i == 0:
            max_l = new_l
         
        if new_l > max_l:
            max_l = new_l
            current_epsilon = epsilon
            current_beta = beta
            current_sigma_0 = sigma_0
            current_sigma_1 = sigma_1
            current_phi = phi

        Results[0] = max_l
        Results[1] = current_epsilon
        Results[2] = current_beta
        Results[3] = current_sigma_0
        Results[4] = current_sigma_1
        Results[5] = phi
   
    return(Results)
