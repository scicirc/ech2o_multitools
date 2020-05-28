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
import copy
import sys
# import glob
import spotpy_forked.spotpy as spotpy
import numpy as np

import func_outputs as outputs
import func_params as params
# import func_likelihoods as likelihoods
import func_objfunctions as objfunc

from distutils.dir_util import copy_tree, remove_tree, mkpath
from distutils.file_util import copy_file

# ==============================================================================


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

    def __init__(self, Config, Opti, Paras, Data, Site, parallel='seq',
                 _used_algorithm='default', dbname=None):

        # Pass on classes
        self.config = copy.copy(Config)
        self.opti = Opti
        self.par = Paras
        self.site = Site
        self.data = Data

        # Initial parameters sampling
        # In case of logarithmic sampling, min/max/guess/sd are
        # are already "converted", and are de-logged during Ech2O's
        # inputs writing (see func_params.sim_inputs)
        self.params = []
        if Opti.initSample == 'normal':
            for i in range(Opti.nvar):
                if np.isnan(Opti.std[i]*Opti.guess[i]):
                    sys.exit('Error: no std and/or guess to sample with',
                             Opti.names[i])
                else:
                    self.params += \
                        [spotpy.parameter.Normal(Opti.names[i],
                                                 stddev=Opti.std[i],
                                                 mean=Opti.guess[i])]
                # print(Opti.names[i], Opti.std[i], Opti.guess[i],
                # self.params[i])
        else:
            for i in range(Opti.nvar):
                self.params += \
                    [spotpy.parameter.Uniform(name=Opti.names[i],
                                              low=Opti.min[i],
                                              high=Opti.max[i])]
                # print(Opti.names[i], Opti.min[i], Opti.max[i],
                # self.params[i])
        # print(self.params)

        # Evaluation data
        self.evals = Opti.obs

        self.curdir = Config.PATH_OUT
        self.owd = Config.PATH_EXEC
        self.parallel = parallel

        # Initialize custom database output.
        # More than one simulation type: save all likelihood
        # (DREAM only use the first one, i.e here the multi-objective)
        # chain, and use explicit simulations name in header
        if dbname is not None:
            self.db_headers = ['Likelihood']
            if self.data.nobs > 1:
                self.db_headers.extend(['like.' + a for a in self.data.names])
            self.db_headers.extend(['par.{0}'.format(x) for x in
                                    self.opti.names])
            self.db_headers.append('chain')
            if self.data.nobs > 1:
                for i in range(self.data.nobs):
                    for j in range(self.data.lsimEff):
                        self.db_headers.extend(['sim.'+self.data.names[i] +
                                                '.'+str(j+1)])
            else:
                for i in range(self.data.lsimEff):
                    self.db_headers.extend(['sim' + '.' + str(i)])
            self.database = open(dbname, 'w')
            self.database.write(",".join(self.db_headers) + "\n")

    # Retrieve parameters
    def parameters(self):
        return spotpy.parameter.generate(self.params)

    # We use the simulation function to write one random parameter set into a
    # parameter file, like it is needed for the HYMOD model, start the model
    # and read the model discharge output data:
    def simulation(self, x):

        os.chdir(self.curdir)

        if self.parallel == 'seq':
            self.PATH_SPA = copy.copy(self.config.PATH_SPA)
            self.PATH_EXEC = copy.copy(self.config.PATH_EXEC)
            self.cfg_ech2o = 'config.ini'
            mkpath(self.PATH_EXEC)
        elif self.parallel == 'mpi':
            # Running n parallel on a unix system.
            # Check the ID of the current computer core
            call = str(int(os.environ['OMPI_COMM_WORLD_RANK'])+2)
        elif self.parallel == 'mpc':
            # Running n parallel on a single (Windows) computer.
            # ID of the current computer core
            call = str(os.getpid())
        else:
            raise 'No call variable was assigned'

        # For parallel computing, a few more preparations
        if self.parallel in ['mpi', 'mpc']:
            # Generate a new input folder with all underlying files
            self.PATH_SPA = self.config.PATH_SPA+call
            copy_tree(self.config.PATH_SPA, self.PATH_SPA)
            # An execution directory next to the "template" one
            self.PATH_EXEC = self.config.PATH_EXEC+call
            mkpath(self.PATH_EXEC)
            # Copy ech2o config file define input maps directory
            self.cfg_ech2o = self.PATH_EXEC+'/config.ini'
            copy_file(self.curdir+'/config.ini', self.cfg_ech2o)
            # Update the "new" config file with output and map locations
            with open(self.cfg_ech2o, 'a') as fw:
                fw.write('Maps_Folder = '+self.PATH_SPA+'\n')
                fw.write('Output_Folder = '+self.PATH_EXEC+'\n')

        # To store simulations (will remain np.nan if simulation fails)
        simulations = np.full((self.data.nobs, self.data.lsimEff),
                              np.nan).tolist()

        try:
            # Create the inputs for ECH2O's run
            params.sim_inputs(self.opti, self.par, self.site,
                              self.PATH_SPA,
                              0, mode='spotpy', paramcur=self.parameters)

            # Run ECH2O
            print('|| running ECH2O...', end='\r')
            start = time.time()
            os.system(self.config.cmde_ech2o + ' ' + self.cfg_ech2o +
                      ' > '+self.PATH_EXEC + '/ech2o.log')
            if self.parallel in ['mpi', 'mpc']:
                print('|| EcH2O run done for core ID#'+call+', using',
                      str(self.config.ncpu), 'cpu(s). Run time:',
                      np.round(time.time()-start, 3), 'seconds')
            else:
                print('|| EcH2O run done using',
                      str(self.config.ncpu), 'cpu(s). Run time:',
                      np.round(time.time()-start, 3), 'seconds')
            # (limit at',self.config.tlimit, ')')

            os.chdir(self.PATH_EXEC)
            # Store outputs: for now restricted to time series
            if self.data.nobs == 1:
                simulations = outputs.read_sim(self.config, self.data,
                                               self.data.names[0]).tolist()
                # print(type(simulations))
            else:
                for i in range(self.data.nobs):
                    oname = self.data.names[i]
                    simulations[i][:] = outputs.read_sim(self.config,
                                                         self.data,
                                                         oname).tolist()

        except():  # 'Model has failed'):
            print('Something went wrong, this run is useless')
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

        # Clean up
        if self.parallel in ['mpc', 'mpi']:
            remove_tree(self.PATH_EXEC)

        return simulations

    # And in a last step, we compare the observed and the simulated data.
    # Please bear in mind that the selection of the Likelihood highly
    # influences the results gained with this algorithm:
    def evaluation(self):
        return self.evals

    def objectivefunction(self, simulation, evaluation, params=None):

        # Use multi-objective function
        like = objfunc.MultiObj(evaluation, simulation, self.data, self.opti)
        return like

    # Custom txt-formatted database (see __init__)
    def save(self, objectivefunctions, parameter, simulations,
             chain, *args, **kwargs):

        if self.data.nobs > 1:
            like_str = ",".join((str(l) for l in objectivefunctions))
            sim_str = []
            for i in range(self.data.nobs):
                sim_str += simulations[i]
                # print([','.join([str(s) for s in simulations[i]])])
            sim_str = ','.join([str(s) for s in sim_str])
        else:
            like_str = objectivefunctions
            sim_str = ",".join([str(s) for s in simulations])

        param_str = ",".join((str(p) for p in parameter))

        line = ",".join([like_str, param_str, str(int(chain)), sim_str]) + '\n'
        self.database.write(line)
