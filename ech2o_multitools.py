#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
Multi-purpose routines ECH2O-iso:
sensitivity analysis, calibration, and ensemble runs in general

-------
Main program
-------
Author: S. Kuppel
Created on 10/2016
-------------------------------------------------
'''

# import os
import sys
# import glob
from optparse import OptionParser
# from itertools import chain

# --- Subroutines
import functions as func
import spotsetup
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
 'calib_SPOTPY': calibration using one of the algorithms available in SPOTPY
                 and adapted here (Opti.SPOTalgo in def file).
                 For now, it's only the Differential Evolution Adaptative
                 Metropolis (DREAM) algorithm
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
# Output directory
parser.add_option("--task", dest="task", metavar="task",
                  help="Identifier of subtasks in jobs array for MonteCarlo simulations")

# MPI parallel computing activated
parser.add_option("--mpi", dest="mpi", metavar="mpi",
                  help="MPI parallel option (0 or 1)")
# Number of CPUs used per EcH2O run (useful in options in mpirun mode)
# overrides the Config.ncpu that may be in definition file
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

# Read the options
(options, args) = parser.parse_args()

# =============================================================================
# == Initialization

# Configuration --------------------------------
(Config, Opti, Obs, Paras, Site) = func.config_init(options)

# Whenever there are EcH2O runs: a few others initializations
if Config.runECH2O == 1:
    # Runs' properties etc.
    func.runs_init(Config, Opti, Obs, Paras, Site, options)
    # Initialize observation names
    # (and read datasets if in DREAM calibration mode)
    func.obs_init(Config, Opti, Obs)

# Parameters: from definition to all values ------
func.param_init(Config, Opti, Paras, Site, options)

# Files and verbose: only do it once
# print(options.mpi+1)
# print('rank', str(rank), 'in initialization routines')
func.files_init(Config, Opti, Paras, Site)

# === Runs ========================================
# -------------------
if Config.mode == 'calib_MCruns':
    # Calibration loop
    func.calibMC_runs(Config, Opti, Obs, Paras, Site)

elif Config.mode == 'calib_SPOTPY':

    if Opti.SPOTalgo in ['DREAM']:
        # Initialize
        spot_setup = spotsetup.spot_setup(Config, Opti, Paras,
                                          Obs, Site,
                                          parallel=Opti.SPOTpar,
                                          _used_algorithm=Opti.SPOTalgo.lower(),
                                          # file extension added automatically
                                          dbname=Opti.dbname)
        if Opti.SPOTalgo == 'DREAM':
            sampler = spotpy.algorithms.dream(spot_setup,
                                              parallel=Opti.SPOTpar,
                                              dbname=Opti.dbname2,
                                              dbformat=Opti.dbformat)
    else:
        sys.exit('Error: other SPOTPY algorithms are to be integrated in' +
                 'ech2o_multitools...')

    # Run DREAM
    r_hat = sampler.sample(Opti.rep, nChains=Opti.nChains,
                           convergence_limit=Opti.convergence_limit,
                           runs_after_convergence=Opti.runs_after_conv)

    # Close the created txt file (if nobs=1, otherwise it's automatic)
    if Obs.nobs == 1:
        spot_setup.database.close()

elif Config.mode == 'forward_runs':
    # Ensemble "forward" runs
    func.forward_runs(Config, Opti, Obs, Paras, Site, options)

elif Config.mode == 'sensi_morris':
    # Simulations when varying the parameters, Morris's one-at-a-time
    func.morris_runs(Config, Opti, Obs, Paras, Site)


# END MAIN
