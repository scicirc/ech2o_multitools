#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Python routines for ECH2O
#
# -------
# Routine: Subroutines for variable sampling
# -------
# Author: S. Kuppel
# Created on 01/2020
# -------------------------------------------------

import time
import os
import glob
import sys
import copy
from pyDOE import *
import random

import spotpy
import Multitool_outputs
import Multitool_params


class spot_setup(object):

    def __init__(self, Config):

        # Parameters
        self.params = []
        for i in range(Config.npar):
            self.params += [spotpy.parameter.Uniform(Config.names[i],
                                                     low=Config.min[i],
                                                     high=Opti.max[i],
                                                     optguess=Opti.guess[i])]

        # Evaluation data
        self.evals = Opti.obs

        self.curdir = blabla
        self.owd = blabla

    # Retrieve parameters
    def parameters(self):
        return spotpy.parameter.generate(self.params)

    # We use the simulation function to write one random parameter set into a
    # parameter file, like it is needed for the HYMOD model, start the model
    # and read the model discharge output data:
    def simulation(self, Opti, Config):
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

        # Create the inputs for ECH2O
        Multitool_params.create_inputs(Opti, Paras, Site, Config, it)

        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('--> running ECH2O')
        start = time.time()
        os.system(Config.cmde_ech2o+' > ech2o.log')
        print('    run time:', time.time() - start,
              'seconds (limit at '+Config.tlimit+')')

        # Store outputs
        simulations = Multitool_outputs.manage_outputs()
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
