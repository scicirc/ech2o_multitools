#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
Multi-purpose routines ECH2O-iso:
sensitivity analysis, calibration, and ensemble runs in general

-------
Routine: Main program
-------
Author: S. Kuppel
Created on 10/2016
-------------------------------------------------
'''

import os, sys, glob
from optparse import OptionParser
# from itertools import chain

# --- Subroutines
import func_init as init
import func_runs as runs
import func_spotsetup as spot_setup
import spotpy_forked.spotpy as spotpy
# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# =============================================================================
# --- Definition ---

'''
What do you want to do ? ==> Specified by the "mode" in Config (def file)
 'calib_MCsampling': generate brute force Monte Carlo ensemble
                   of parameters sets for calibration
 'calib_MCruns': runs the model for all MC-sampled parameters
 'calib_DREAM': calibration using the Differential Evolution Adaptative
                Metropolis (DREAM) algorithm, using the Spotpy package
 'forward_runs': runs the model for an ensemble of runs,
                 (usually the best configurations from the calibration.
                 Allows to look at observations not used in calibration)
 'sensi_morris': performs a Morris sensitivity analysis
'''

# Other options: see subroutine in func_init.py
# Configuration file
parser.add_option("--file", dest="file", metavar="FILE",
                  help="Name of the file that defines the calibration " +
                  "to perform")
# Output directory
parser.add_option("--outdir", dest="outdir", metavar="outdir",
                  help="Output directory")

# MPI parallel computing activated
parser.add_option("--mpi", dest="mpi", metavar="mpi",
                  help="MPI parallel option (0 or 1)")
# Number of CPUs used per EcH2O run (useful in options in mpirun mode)
parser.add_option("--ncpu", dest="ncpu", metavar="ncpu",
                  help="Number of CPUs for each EcH2O run")
# Use scratch ? (saves EcH2O tmp outputs on scratch, then saves on users after
# each run)
parser.add_option("--scratch", dest="scratch", metavar="scratch",
                  help="Uses of localscratch (1: fastest storage) " +
                  "or shared (2: fast)")

# == Options for specific routines modes ==

# -- If mode == 'calib_MCsampling'
# Sampling method
parser.add_option("--sampling", dest="sampling", metavar="sampling",
                  help="Sampling method: LHS (uniform), LHS_r (min corr), " +
                  "LHS_m (max distance)")

# -- If mode != 'calib_MCsampling' (i.e., EcH2O is actually run)
# Name of the ECH2O config tracking file (if mode!='calib_MCsampling', and
# tracking activated in simulation)
parser.add_option("--cfgTrck", dest="cfgTrck", metavar="cfgTrck",
                  help="Name of the ECH2O configTrck file")
# Name of the reference ECH2O config file
parser.add_option("--cfg", dest="cfg", metavar="cfg",
                  help="Name of the ECH2O config file")
# Time limit for the runs
parser.add_option("--tlimit", dest="tlimit", metavar="tlimit",
                  help="Time limit of one ECH2O run (in seconds)")

# -- If mode == 'calib_MCruns'
# Spinup ? (if post optim, mode = 2 ) ?
parser.add_option("--spinup", dest="spinup", metavar="spinup",
                  help="Spinup switch (0) or length (days)")
# Restart ? if calibration stopped for some reasons
# --> 0 if not (by default), or 1: the code will take the antelast that worked
parser.add_option("--restart", dest="restart", metavar="restart",
                  help="Restart (1) or not (0)")

# -- If mode == 'forward_runs'
# Input dir (if post optim, mode = 2) ?
parser.add_option("--inEns", dest="inEns", metavar="inEns",
                  help="Ensenble input prefix")
# Size of the ensemble (if post optim, mode = 2 ) ?
parser.add_option("--nEns", dest="nEns", metavar="nEns",
                  help="Size of the ensemble")
# Average maps ? (saves a lot of disk space)
parser.add_option("--MapAv", dest="MapAv", metavar="MapAv",
                  help="0 or 1, to temporally average spatial outputs)")
parser.add_option("--MapAvT", dest="MapAvT", metavar="MapAvT",
                  help="Averaging period (week, month, season)")

parser.add_option("--OMP_it", dest="OMP_it", metavar="OMP_it",
                  help="Ensemble iteration number")

# If mode == 'sensi_morris'
# Only generating trajectories ?
parser.add_option("--MSinit", dest="MSinit", metavar="MSinit",
                  help="Switch for Morris sensitivity: " +
                  "0=trajectories generation, 1=runs")
parser.add_option("--MSspace", dest="MSspace", metavar="MSspace",
                  help="Walk in the parameter space for Morris: " +
                  "'trajectory' or 'radial'")

# Read the options
(options, args) = parser.parse_args()

# For MPI
# rank = -1
# if options.mpi == '1':
#     options.mpi = int(options.mpi)
#     from mpi4py import MPI
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     size = comm.Get_size()
#     print('rank/size', rank,'/', size, 'before initialiation routines')

# =============================================================================
# == Initialization

# Configuration --------------------------------
(Config, Opti, Data, Paras, Site) = init.config(options)

# Parameters: from definition to all values ------
init.parameters(Config, Opti, Paras, Site, options)

# Whenever there are EcH2O runs: a few others initializations
if Config.runECH2O == 1:
    # Runs' properties etc.
    init.runs(Config, Opti, Data, Paras, Site, options)
    # Initialize observation names
    # (and read datasets if in DREAM calibration mode)
    init.observations(Config, Opti, Data)

# Files and verbose: only do it once
if options.mpi != 1 or \
       (options.mpi == 1):  # and rank == 0):

    # print(options.mpi+1)
    # print('rank', str(rank), 'in initialization routines')

    # ==== Introductory verbose
    # if Opti.parallel == False or \
    #   (Opti.DREAMpar != 'mpi' or Config.MPIrank == 0):
    print('')
    print('******************************************************************')
    print('The user provided definition file is: '+options.file)
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    if Config.mode.split('_')[0] == 'calib':
        print('CALIBRATION with EcH2O: \n')
    elif Config.mode == 'forward_runs':
        print('ENSEMBLE RUNS with EcH2O')
    elif Config.mode == 'sensi_morris':
        print('MORRIS SENSITIVITY with EcH2O: ')
        print('- construction of the trajectories')
        print('- forward runs')
        print('- storage of outputs and info for posterior analysis : ' +
              'elementary effects, etc.')
        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
        print('')
    print('')
    print('Original/template maps & parameters:\n', Config.PATH_SPA_REF)
    print('Climate data:\n', Config.PATH_CLIM)
    print('Final outputs:\n', Config.PATH_OUT)
    if Opti.parallel is False:
        print('Maps & parameters:', Config.PATH_SPA)
        print('Raw output from simulations:', Config.PATH_EXEC)
    print('The user provided CFG file is: '+Config.cfg_ech2o)
    if Site.isTrck == 1:
        print('The user provided CFGtrck file is: ' + Config.cfgTrck_ech2o)
    print('-----------------------------------------')
    # -- How many variables ?
    print('Total number of parameters :', Opti.nvar)
    print('')

    # Copy, edit directories / files
    # In MPI parallel mode, care as to be taken to write only once
    init.files(Config, Opti, Paras, Site)

#     if options.mpi == 1:  # i.e., rank == 0:
#         for i in range(1,size):
#             print('Rank 0: rank',i, ', wait for me!')
#             req = comm.isend('waited for you, chillax...', dest=i, tag=1)
#             req.wait()

# elif options.mpi == 1:  # i.e., rank > 0
#     # sys.path.insert(0, cwd_tmp)
#     req = comm.irecv(source=0, tag=1)
#     data = req.wait()
#     print('rank', str(rank), data)
    

    # -- Import classes and setup from the def file
#     Config = None #__import__(file_py).Config
#     Opti = None #__import__(file_py).Opti
#     Data = None #__import__(file_py).Data
#     Paras = None # __import__(file_py).Paras
#     Site = None # __import__(file_py).Site

# #print(Config.__dict__.keys())
# #if options.mpi == 1:
#     # In MPI mode, broadcast the class to all processes
#     # Config = comm.bcast(Config, root=0)
#     Opti = comm.bcast(Opti, root=0)
#     Data = comm.bcast(Data, root=0)
#     Paras = comm.bcast(Paras, root=0)
#     Site = comm.bcast(Site, root=0)

# print(Config.__dict__.keys())

# sys.exit()

# === Runs ========================================
# -------------------
if Config.mode == 'calib_MCruns':
    # Calibration loop
    runs.calibMC_runs(Config, Opti, Data, Paras, Site)

elif Config.mode == 'calib_DREAM':

    # Initialize
    spot_setup = \
        spot_setup.spot_setup(Config, Opti, Paras,
                              Data, Site,
                              parallel=Opti.DREAMpar,
                              _used_algorithm='dream',
                              # file extension added automatically
                              dbname=Config.PATH_OUT+'/DREAMech2o')
    sampler = \
        spotpy.algorithms.dream(spot_setup, parallel=Opti.DREAMpar,
                                dbformat='custom')
    # dbname=Config.PATH_OUT+'/DREAM_ech2o',
    # dbformat='csv')

    # Run DREAM
    r_hat = sampler.sample(Opti.rep, nChains=Opti.nChains,
                           convergence_limit=Opti.convergence_limit,
                           runs_after_convergence=Opti.runs_after_conv)

    # Close the created txt file (if nobs=1, otherwise it's automatic)
    if Data.nobs == 1:
        spot_setup.database.close()

elif Config.mode == 'forward_runs':
    # Ensemble "forward" runs
    runs.forward_runs(Config, Opti, Data, Paras, Site, options)

elif Config.mode == 'sensi_morris' and Config.MSinit == 0:
    # Simulations when varying the parameters, Morris's one-at-a-time
    runs.morris_runs(Config, Opti, Data, Paras, Site)


# END MAIN
