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
import multiprocessing

from distutils.dir_util import copy_tree, remove_tree, mkpath
from distutils.file_util import copy_file
from contextlib import ExitStack

# ==============================================================================


class spot_setup(object):

    def __init__(self, Config, Opti, Paras, Obs, Site, parallel='seq',
                 _used_algorithm='default', dbname=None):

        # Pass on classes
        self.config = copy.copy(Config)
        self.opti = Opti
        self.par = Paras
        self.site = Site
        self.obs = Obs

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

        # Evaluation data
        self.evals = copy.copy(Opti.obs)

        self.parallel = copy.copy(parallel)

        # Initialize custom database output.
        # More than one simulation type: save all likelihood
        # (DREAM only use the first one, i.e here the multi-objective)
        # chain, and use explicit simulations name in header
        if dbname is not None:
            ## -- Likelihood, parameters and time series
            # If multi-objecive calibration, we'll write one file per
            # simulation type containing "individual" likelihood and
            # times series, plus a "general" file with parameters
            # and likelihoods
            if self.obs.nobs > 1:
                self.dbheaders = {}
                self.filenames = [dbname+'_LikeParGeneral.txt'] + \
                    [dbname+'.'+a+'.txt' for a in self.obs.names]
                self.outnames = ['general'] + self.obs.names
                # General file
                self.dbheaders['general'] = ['Iteration']
                self.dbheaders['general'].append('Likelihood')
                self.dbheaders['general'].extend(['like.' + a
                                                  for a in self.obs.names])
                self.dbheaders['general'].extend(['par.{0}'.format(x) for x in
                                                  self.opti.names])
                self.dbheaders['general'].append('chain')
                # Simulation-specific files
                for oname in self.obs.names:
                    self.dbheaders[oname] = ['Iteration']
                    self.dbheaders[oname].append('like.'+oname)
                    self.dbheaders[oname].extend(['par.{0}'.format(x) for x in
                                                 self.opti.names])
                    for j in range(self.obs.saveL):
                        self.dbheaders[oname].extend(['sim.t'+str(j+1)])
                    self.dbheaders[oname].append('chain')
                # Write headers
                with ExitStack() as stack:
                    # Open all files
                    self.database = [stack.enter_context(open(fname, 'w')) for
                                     fname in self.filenames]
                    # All opened files will automatically be closed at the end
                    # ofthe with statement, even if attempts to open files
                    # later in the list raise an exception
                    for i in range(self.obs.nobs+1):
                        self.database[i].write(",".
                                               join(self.dbheaders[
                                                   self.outnames[i]])+"\n")
            else:
                # Only one calibration datasets: some tweaks from the
                # standard .csv style database
                self.filename = dbname+'.txt'
                self.dbheaders = ['Iteration']
                self.dbheaders.append('Likelihood')
                self.dbheaders.extend(['par.{0}'.format(x) for x in
                                       self.opti.names])
                for i in range(self.obs.saveL):
                    self.dbheaders.extend(['sim' + '.t' + str(i)])
                self.dbheaders.append('chain')
                # Write
                self.database = open(self.filename, 'w')
                self.database.write(",".join(self.dbheaders) + "\n")

            ## -- Diagnostics from the calibration algorithm (for now
            # mostly suited for DREAM)
            if self.opti.SPOTalgo == 'DREAM':
                self.f_diagnostics = dbname + '_Diagnostics.txt'
                headers = ['Iteration','BestLike']
                headers.extend(['PctAccept.ch'+str(c) for c in 
                                range(1,self.opti.nChains+1)])
                headers.extend(['ConvRates.{0}'.format(x) for x in self.opti.names])
                with open(self.f_diagnostics, 'w') as f:
                    f.write(','.join(headers)+'\n')
            else:
                print('Diagnostics for algorithm progression: file output',
                      'under construction...')
            

    # Retrieve parameters
    def parameters(self):
        return spotpy.parameter.generate(self.params)

    # We use the simulation function to write one random parameter set into a
    # parameter file, like it is needed for the HYMOD model, start the model
    # and read the model discharge output data:
    def simulation(self, x):

        os.chdir(self.config.PATH_OUT)

        if self.parallel == 'seq':
            # Sequential
            self.PATH_SPA = copy.copy(self.config.PATH_SPA)
            self.PATH_EXEC = copy.copy(self.config.PATH_EXEC)
            self.cfg_ech2o = 'config.ini'
            mkpath(self.PATH_EXEC)
        elif self.parallel == 'mpi':
            # Running n parallel on a unix system.
            # Check the ID of the current mpi task using mpi4py
            # call = str(int(self.rank))
            if 'OMPI_COMM_WORLD_RANK' in os.environ.keys():
                call = str(int(os.environ['OMPI_COMM_WORLD_RANK']))
            elif 'PMI_RANK' in os.environ.keys():
                call = str(int(os.environ['PMI_RANK']))
            else:
                sys.exit('The ID of this task could not be found...')
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
            os.makedirs(self.PATH_EXEC, exist_ok=True)
            # REmove ech2o log file if it exists
            if os.path.isfile(self.PATH_EXEC+'/ech2o.log'):
                os.remove(self.PATH_EXEC+'/ech2o.log')
            # Copy ech2o config file define input maps directory
            self.cfg_ech2o = self.PATH_EXEC+'/config.ini'
            copy_file(self.config.FILE_CFGdest, self.cfg_ech2o)
            # Update the "new" config file with output and map locations
            with open(self.cfg_ech2o, 'a') as fw:
                fw.write('Maps_Folder = '+self.PATH_SPA+'\n')
                fw.write('Output_Folder = '+self.PATH_EXEC+'\n')

        # Create the inputs for ECH2O's run
        # Here, x is the current set of sampled parameter values.
        # It is ordered as in Opti.names etc., so it can be used as is            
        # print('rank', call, ': generating maps & veg params...')
        params.sim_inputs(self.config, self.opti, self.par, self.site,
                          self.PATH_SPA, 0, mode='spotpy', paramcur=x)

        # nan values so that simulations can be returned even if the run fails
        if self.obs.nobs > 1:
            simulations = np.full((self.obs.nobs, self.obs.saveL),
                                  np.nan).tolist()
        else:
            simulation = [np.nan] * self.obs.saveL

        # Run EcH2O
        print('|| running ECH2O...', end='\r')
        start = time.time()

        try:
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

            # Store outputs: for now restricted to time series
            os.chdir(self.PATH_EXEC)
        
            if self.obs.nobs == 1:
                simulations = outputs.read_sim(self.config, self.obs,
                                               self.obs.names[0]).tolist()
                if(len(simulations) < self.obs.saveL):
                    print('sim length:', len(simulations),
                          ', expected:', self.obs.saveL)
                    simulations = [np.nan] * self.obs.saveL
            else:
                # To store simulations
                simulations = np.full((self.obs.nobs, self.obs.saveL),
                                      np.nan).tolist()
                for i in range(self.obs.nobs):
                    oname = self.obs.names[i]
                    simulations[i][:] = outputs.read_sim(self.config,
                                                         self.obs,
                                                         oname).tolist()
                    if(len(simulations[i]) < self.obs.saveL):
                        print('sim length:', len(simulations[i]),
                              ', expected:', self.obs.saveL)
                        simulations[i][:] = [np.nan] * self.obs.saveL
               
        except:
            #'Model has failed'
            print('Something went wrong, this run is useless')        
            # simulations = [np.nan] * self.obs.saveL
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

        # sys.exit()
        os.chdir(self.config.PATH_OUT)

        # Clean up
        # if self.parallel in ['mpc', 'mpi']:
        #    remove_tree(self.PATH_EXEC)

        return simulations

    # And in a last step, we compare the observed and the simulated data.
    # Please bear in mind that the selection of the Likelihood highly
    # influences the results gained with this algorithm:
    def evaluation(self):
        return self.evals

    def objectivefunction(self, simulation, evaluation, params=None):

        # Use multi-objective function
        like = objfunc.MultiObj(evaluation, simulation, self.obs, self.opti)
        # Check if there was any issue with the simulation outputs
        if np.isnan(like).any():
            self.runFail = True
        else:
            self.runFail = False

        return like

    # Custom txt-formatted database (see __init__)
    def save(self, objfuncs, parameter, simulations,
             chains=1, rep=1, *args, **kwargs):

        param_str = ",".join([str(p) for p in parameter])

        if self.obs.nobs > 1:
            with ExitStack() as stack:
                # Open all files for append-writing
                self.database = [stack.enter_context(open(fname, 'a')) for
                                 fname in self.filenames]
                # "general" file
                like_str = ','.join([str(l) for l in objfuncs])
                self.database[0].write(",".join([str(rep+1), like_str, param_str,
                                                 str(int(chains+1))])+'\n')
                # Sim-specific files
                for i in range(self.obs.nobs):
                    sim_str = ','.join([str(s) for s in simulations[i]])
                    line = ','.join([str(rep+1), str(objfuncs[i+1]), param_str, 
                                     sim_str, str(int(chains+1))]) + '\n'
                    self.database[i+1].write(line)
        else:
            # One calibration datasets: one file
            sim_str = ",".join([str(s) for s in simulations])
            self.database.write(",".join([str(rep+1), objfuncs, param_str,
                                          str(int(chains+1)), sim_str]) + '\n')

    # Custom txt-formatted database (see __init__)
    def save_diagnostics(self, it, objfunc_ref, accepted, conv_rates):

        ## -- Diagnostics from the calibration algorithm (for now
        # mostly suited for DREAM)
        if self.opti.SPOTalgo == 'DREAM':
            line = ','.join([str(it), str(np.round(objfunc_ref,2)),
                             ','.join([str(s) for s in np.round(accepted,2)]),
                             ','.join([str(s) for s in np.round(conv_rates,4)])])
            # print(line)
            print('')
            with open(self.f_diagnostics, 'a') as f:
                f.write(line+ '\n')
     
