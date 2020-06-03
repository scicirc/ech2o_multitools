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

import os
# import sys
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

# =============================================================================

# Current working directory (temporary, just to import the def file)
cwd_tmp = os.getcwd()+'/'

# =============================================================================
# == Initialization

# Configuration --------------------------------
(Config, Opti, Data, Paras, Site) = init.config(cwd_tmp, options)

# Parameters: from definition to all values ------
init.parameters(Config, Opti, Paras, Site, options)

# Whenever there are EcH2O runs: a few others initializations
if Config.runECH2O == 1:
    # Runs' properties etc.
    init.runs(Config, Opti, Data, Paras, Site, options)
    # Initialize observation names
    # (and read datasets if in DREAM calibration mode)
    init.observations(Config, Opti, Data)

# -- How many variables ?
print('Total number of parameters :', Opti.nvar)

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
