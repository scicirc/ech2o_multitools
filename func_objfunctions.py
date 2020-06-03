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
import pandas as pd
# import math
import sys
import spotpy_forked.spotpy as spotpy

# ==============================================================================
# -- Multi-objective function(s)


def MultiObj(obs, sim, Data, Opti, w=False):

    like = 0
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

        # Now use your favorite likelihood estimator for each obs type

        # For GWD, we center using mean
        if oname.split('_')[0] in ['GWD', 'WTD', 'GWL', 'WTL']:
            # print(oname, 'remove mean')
            # s -= np.mean(s)
            # o -= np.mean(o)
            # Use RMSE for water table
            # L = - spotpy.objectivefunctions.rmse(o, s) / np.mean(o)
            # LogLikehood as used in Vrugt et al. (2016)
            # o_err = np.repeat(0.3, len(o)).tolist()
            # L = spotpy.likelihoods.logLikelihood(o, s,
            #                                     measerror=o_err)
            # Log Gaussian likelihood without error (temporary)
            L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s)

        else:
            # MAE for streamflow
            # L = - spotpy.objectivefunctions.mae(o, s) / np.mean(o)
            # LogL gaussian likehood as used in Vrugt et al. (2016)
            # L = spotpy.likelihoods.logLikelihood(o, s)
            # Log Gaussian likelihood without error (temporary)
            L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s)

        Ltot += [L]  # list of all likelihoods
        like += L  # "final" likelihood

    # Several datasets: list multi-objective function and individual ones
    # Only the first (multi-obj) one will be used by DREAM
    if Data.nobs > 1:
        like = [like] + Ltot

    # print(np.round(like, 2))

    return like
