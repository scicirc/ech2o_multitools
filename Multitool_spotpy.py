'''
Python routines for ECH2O

Routine: Subroutines for variable sampling
-------
:author: Sylvain Kuppel

Adapted from spot_setup_hymod_unix.py and tutorial_dream_hymod.py
(available at github.com/thouska/spotpy/tree/master/spotpy/examples),
created by Tobias Houska
'''

import time
import os
# import sys
import glob
import spotpy
import numpy as np

import Multitool_outputs as outputs
import Multitool_params as params
import Multitool_likelihoods as likelihoods

# # ==============================================================================


# def spotpy_dream(Config, Opti):

#     # Initialize
#     spot_setup = spot_setup(Config, Opti, _used_algorithm='dream')

#     sampler = spotpy.algorithms.dream(spot_setup, dbname='DREAM_ech2o',
#                                       dbformat='csv')
#     r_hat = sampler.sample(Opti.rep, nChains=Opti.nChains,
#                            convergence_limit=Opti.conv_lim,
#                            runs_after_convergence=Opti.run_after_conv)
# ==============================================================================


class spot_setup(object):

    def __init__(self, Config, Opti, parallel='seq',
                 _used_algorithm='default'):

        # Parameters
        sdmult = 0.1

        self.params = []
        for i in range(Opti.nvar):
            minp = Opti.min[i]
            maxp = Opti.max[i]
            print(Opti.names[i], minp, maxp)
            # Logarithmic sampling, will be de-logged during Ech2O's
            # inputs writing
            if Opti.log[i] == 1:
                self.params += [spotpy.parameter.Uniform(Opti.names[i],
                                                         low=np.log10(minp),
                                                         high=np.log10(maxp))]
                # tmp = spotpy.parameter.Normal(name=Opti.names[i],
                #                               low=np.log10(minp),
                #                               high=np.log10(maxp),
                #                               stddev=sdmult*np.log10(maxp/minp))
            else:
                self.params += \
                    [spotpy.parameter.Uniform(name=Opti.names[i],
                                              low=minp, high=maxp)]
                # [spotpy.parameter.Normal(name=Opti.names[i],
                #                          low=minp, high=maxp,
                #                          stddev=sdmult*(maxp-minp))]

        print(self.params)

        # Evaluation data
        self.evals = Opti.obs

        self.curdir = Config.PATH_MAIN
        self.owd = Config.PATH_EXEC

    # Retrieve parameters
    def parameters(self):
        return spotpy.parameter.generate(self.params)

    # We use the simulation function to write one random parameter set into a
    # parameter file, like it is needed for the HYMOD model, start the model
    # and read the model discharge output data:
    def simulation(self, Opti, Config, Paras, Data, Site):

        # Create the inputs for ECH2O
        params.sim_inputs(Opti, Paras, Site, Config, it)

        # To store simulations
        simulations = np.full((Data.nobs, Config.trimL), np.nan)

        # Create or clean exec directory
        if len(glob.glob(Config.PATH_EXEC)) == 0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('--> running ECH2O')

        try:
            start = time.time()
            os.system(Config.cmde_ech2o+' > '+Config.PATH_EXEC+'/ech2o.log')
            print('    run time:', time.time() - start, 'seconds (limit at ',
                  Config.tlimit, ')')

            # Store outputs: for now restricted to time series
            for i in range(Data.nobs):
                oname = Data.names[i]
                simulations[i, :] = outputs.readsim(Config, Data, oname)

        except('Model has failed'):
            print('prout')
            # Report param config that failed
            # f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
            # if len(glob.glob(f_failpar)) == 0:
            #     with open(f_failpar, 'w') as f_in:
            #         f_in.write('Sample,'+','.join(Opti.names)+'\n')
            #         f_in.write(str(it+1)+','+','.join([str(x) for x in
            #                                            Opti.x])+'\n')
            # else:
            #     with open(f_failpar, 'a') as f_in:
            #         f_in.write(str(it+1)+','+','.join([str(x) for x in
            #                                            Opti.x])+'\n')
            # # Save the failed screen outputs
            # os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
            #           Config.PATH_OUT + '/ech2o_'+it+'.log')

        os.chdir(self.curdir)

        return simulations

    # And in a last step, we compare the observed and the simulated data.
    # Please bear in mind that the selection of the Likelihood highly
    # influences the results gained with this algorithm:
    def evaluation(self):
        return self.evals

    def objectivefunction(self, simulation, evaluation, params=None):
        # Use multi-objective GL from Shoups and Vrught, as in
        # Knighton et al., 2017, 2020
        like = likelihoods.ShoupsVrught_GL(evaluation, simulation)
        return like
