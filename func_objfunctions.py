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

import numpy as np
import os
import pandas as pd
# import math
import sys
import spotpy_forked.spotpy as spotpy

# ==============================================================================
# -- Multi-objective function(s)


def MultiObj(obs, sim, Data, Opti, w=False):

    like = 0
    Ltot = []

    # Sanity check: did the simulation run?
    # (necessary to avoid error in the subsequent postprocessing)
    if sim is None:
        if Data.nobs > 1:
            like = [np.nan for i in range(Data.nobs+1)]
            # Data.nobs+1 because of total objfunc + indiv ones
        else:
            like = np.nan

    else:    
        for i in range(Data.nobs):

            oname = Data.names[i]

            # Have obervation and simulations matching the same time period
            # obs: already pre-processed
            tobs = pd.to_datetime(obs[oname]['Date'].values)
            o = np.asanyarray(obs[oname]['value'].values)
            
            # First step for sim: trim sim to obs timespan
            # + only keep dates with obs (even if nan)
            # print(sim)
            # print(sim.shape)
            if Data.nobs == 1:
                s = np.asanyarray([sim[j] for j in range(Data.lsimEff) if
                                   Data.simt[j] in tobs])
            else:
                s = np.asanyarray([sim[i][j] for j in range(Data.lsimEff) if
                                   Data.simt[j] in tobs])

            # Second step (both o and s): remove nan due to gaps in obs
            tmp = s*o
            s = np.asanyarray([s[k] for k in range(len(tmp)) if not
                               np.isnan(tmp[k])])
            o = np.asanyarray([o[j] for j in range(len(tmp)) if not
                               np.isnan(tmp[j])])
        
            # Another sanity check: if there any data/sim left after nan screening?
            if s.__len__() == 0 or o.__len__() == 0:
                L = np.nan

            else:
                # Now use your favorite likelihood estimator for each obs type

                # Specific treatment for different obs types?
                if oname.split('_')[0] in ['GWD', 'WTD', 'GWL', 'WTL']:
                    # print(oname, 'remove mean')
                    # s -= np.mean(s)
                    # o -= np.mean(o)
                    # Log Gaussian likelihood without error (temporary)
                    # L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s) * \
                    #     1 / (o.__len__() * np.mean(o))
                    # Log Gaussian with error set using std
                    # LogLikehood as used in Vrugt et al. (2016)
                    # Using standard deviation as indicator of error
                    o_err = np.repeat(np.std(o)*0.25, o.__len__())
                    L = spotpy.likelihoods.logLikelihood(o, s, measerror=o_err)
                else:
                    # Log Gaussian likelihood without error (temporary)
                    # L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s) * \
                    #    1 / (o.__len__() * np.mean(o))
                    # Log Gaussian with error set using std
                    # LogLikehood as used in Vrugt et al. (2016)
                    # Using standard deviation as indicator of error...
                    o_err = np.std(o) + 0.1*o
                    L = spotpy.likelihoods.logLikelihood(o, s, measerror=o_err)
                
                # Normalize by data length
                L /= o.__len__()

            Ltot += [L]  # list of all likelihoods
            like += L  # "main" likelihood (used by algorithm)

        # Several datasets: list multi-objective function and individual ones
        # Only the first (multi-obj) one will be used by the algorithm
        if Data.nobs > 1:
            like = [like] + Ltot

    # print(np.round(like, 2))

    return like
