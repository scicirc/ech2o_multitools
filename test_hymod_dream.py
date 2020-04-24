#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************

import spotpy
import cmf
from datetime import timedelta, datetime  # Standard in Python to work with data that have date and time stamp
import Load_Data as loader # Used to import meterological data and the evaluation data from csv-files.
import numpy as np

# This chapter shows you, how to link the external hydrological model HYMOD
# with SPOTPY (works only on Windows systems).

# We use the hydrological model HYMOD as an example, to calibrate it with the
# Differential Evolution Adaptive Metropolis (DREAM) algorithm. For detailed
# information about the underlying theorie, have a look at the Vrugt (2016).
# The SPOTPY package comes with an example which is desgined to help you to
# set up your own research project.

# First, we need to setup the model within a spot setup class.
# The model needs some meteorological input data and five parameters to estimate discharge:

## Connect HYMOD with SPOTPY
# Here we use to _init_ function, to initialize the parameter for our model.


class spot_setup(object):
    def __init__(self):

        self.params = [spotpy.parameter.Uniform('x1',low=1.0 , high=500,  optguess = 412.33),
                       spotpy.parameter.Uniform('x2',low=0.1 , high=2.0,  optguess = 0.1725),
                       spotpy.parameter.Uniform('x3',low=0.1 , high=0.99, optguess = 0.8127),
                       spotpy.parameter.Uniform('x4',low=0.0 , high=0.10, optguess=0.0404),
                       spotpy.parameter.Uniform('x5',low=0.1 , high=0.99, optguess=0.5592)
                       ]

        self.curdir = os.getcwd()
        self.owd = os.path.realpath(__file__)+os.sep+'..'
        self.evals = list(np.genfromtxt(self.owd+os.sep+'hymod'+os.sep+'bound.txt',skip_header=65)[:,3])[:730]
        self.Factor = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24)
        print(len(self.evals))

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    # We use the simulation function to write one random parameter set into a
    # parameter file, like it is needed for the HYMOD model, start the model
    # and read the model discharge output data:
    def simulation(self, x):
        os.chdir(self.owd+os.sep+'hymod')
        if sys.version_info.major == 2:
            params = file('Param.in', 'w')
        elif sys.version_info.major == 3:
            params = open('Param.in', 'w')
        else:
            raise Exception("Your python is too old for this example")
        for i in range(len(x)):
            if i == len(x):
                params.write(str(round(x[i], 5)))
            else:
                params.write(str(round(x[i], 5))+' ')
        params.close()
        os.system('HYMODsilent.exe')

        # try:
        if sys.version_info.major == 2:
            SimRR = file('Q.out', 'r')
        elif sys.version_info.major == 3:
            SimRR = open('Q.out', 'r')
        else:
            raise Exception("Your python is too old for this example")
        simulations = []
        for i in range(64):
            SimRR.readline()
        for i in range(730):
            val = SimRR.readline()
            simulations.append(float(val)*self.Factor)
        # except:#Assign bad values - model might have crashed
        #    SimRR = 795 * [np.nan]
        os.chdir(self.curdir)

        return simulations

    # And in a last step, we compare the observed and the simulated data.
    # Here we choose one of the implemented Likelihood functions in the SPOTPY
    # package. Please mind that the selection of the Likelihood highly
    # influences the results gained with this algorithm:
    def evaluation(self):
        return self.evals

    def objectivefunction(self, simulation, evaluation, params=None):
        like = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(evaluation,
                                                                 simulation)
        return like

# -------------------------------------------------------------------------
# Now we can initialize the Hymod example:
spot_setup = spot_setup()

# Create the Dream sampler of spotpy,
# alt_objfun is set to None to force SPOTPY to jump into the def
# objectivefunction in the spot_setup class (default is
# spotpy.objectivefunctions.log_p).
# Results are saved in a DREAM_hymod.csv file:
sampler = spotpy.algorithms.dream(spot_setup, dbname='DREAM_hymod',
                                  dbformat='csv', alt_objfun=None)

# Select number of maximum repetitions, the number of chains used by DREAM
# (default = 5) and set the Gelman-Rubin convergence limit (default 1.2).
# We further allow 100 runs after convergence is achieved:
nChains = 4
convergence_limit = 1.2
runs_after_convergence = 100

# We start the sampler and collect the gained r_hat convergence
# values after the sampling:
r_hat = sampler.sample(rep, nChains=nChains,
                       convergence_limit=convergence_limit,
                       runs_after_convergence=runs_after_convergence)

