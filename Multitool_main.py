#!/usr/bin/python3
# -*- coding: utf-8 -*-
# *************************************************
#
# Multi-purpose routines ECH2O-iso:
# sensitivity analysis, calibration, and ensemble runs in general
#
# -------
# Routine: Main program
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import os
from optparse import OptionParser
# from itertools import chain

# --- Subroutines
import Multitool_init as init
import Multitool_runs as runs

# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# =============================================================================
# --- Definition ---

# What do you want to do ?
# 'calib_sampling': generate ensemble of parameters sets for calibration
# 'calib_runs': runs the model for all sampled parameters
# 'forward_runs': runs the model for an ensemble of runs, usually the best
#     configurations from the calibration. Allows to look at observations not
#     used in calibration, e.g. maps
# 'sensi_morris': performs a Morris sensitivity analysis
parser.add_option("--mode", dest="mode", metavar="mode",
                  help="Switch ('calib_sampling','calib_runs'," +
                  "'forwards_runs','sensi_morris')")

# Other options: see subroutine in Multitool_init.py
# Configuration file
parser.add_option("--file", dest="file", metavar="FILE",
                  help="Name of the file that defines the calibration " +
                  "to perform")
# Output directory
parser.add_option("--outdir", dest="outdir", metavar="outdir",
                  help="Output directory")
# Number of CPUs used
parser.add_option("--ncpu", dest="ncpu", metavar="ncpu",
                  help="Number of CPUs used")
# Use scratch ? (saves EcH2O tmp outputs on scratch, then saves on users after
# each run)
parser.add_option("--scratch", dest="scratch", metavar="scratch",
                  help="Uses of localscratch (1: fastest storage) " +
                  "or shared (2: fast)")
# Resolution of EcH2O
parser.add_option("--Resol", dest="Resol", metavar="Resol",
                  help="resolution (in m) for simulation, by default 30m)")
# Report BasinSummary.txt ? (by default = 0)
parser.add_option("--BSum", dest="BSum", metavar="BSum",
                  help="report BasinSummary.txt files for each simulation")

# == Options for specific routines modes ==

# -- If mode == 'calib_sampling'
# Sampling method
parser.add_option("--sampling", dest="sampling", metavar="sampling",
                  help="Sampling method: LHS (uniform), LHS_r (min corr), " +
                  "LHS_m (max distance)")

# -- If mode != 'calib_sampling' (i.e., EcH2O is actually run)
# Name of the ECH2O config tracking file (if mode!='calib_sampling', and
# tracking activated in simulation)
parser.add_option("--cfgTrck", dest="cfgTrck", metavar="cfgTrck",
                  help="Name of the ECH2O configTrck file")
# Name of the ECH2O executable
parser.add_option("--exe", dest="exe", metavar="exe",
                  help="Name of the ECH2O exec file")
# Name of the reference ECH2O config file
parser.add_option("--cfg", dest="cfg", metavar="cfg",
                  help="Name of the ECH2O config file")
# Flux tracking activated ?
parser.add_option("--isTrck", dest="isTrck", metavar="isTrck",
                  help="Switch of water tracking")
# Time limit for the runs
parser.add_option("--tlimit", dest="tlimit", metavar="tlimit",
                  help="Time limit of one ECH2O run (in seconds)")
# Trim outputs ? If so, the beginning and/or ending of outputs should not be
# saved (e.g. transient state)
parser.add_option("--trimB", dest="trimB", metavar="trimB",
                  help="Drop the beginning of outputs: 0 if not, " +
                  "length otherwise")
parser.add_option("--trimBmap", dest="trimBmap", metavar="trimBmap",
                  help="Drop the beginning of map outputs: 0 if not, " +
                  "length otherwise")
parser.add_option("--trimL", dest="trimL", metavar="trimL",
                  help="Length of the trim (if trimB>0): full trim by " +
                  "default, integer otherwise (has to be larger than " +
                  "total length - trimB)")

# -- If mode == 'calib_runs'
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


class Config:
    # set working directory
    PATH_MAIN = os.getcwd()+'/'
    cfgdir = PATH_MAIN+'Input_Configs/'


# =============================================================================
# == Initialization

# Configuration --------------------------------
(Opti, Data, Paras, Site) = init.config(Config, options)

# Parameters: from definition to all values ------
init.parameters(Config, Opti, Paras, Site, options)

# Whenever there are EcH2O runs: a few others initializations
if Config.runECH2O == 1:
    init.runs(Config, Opti, Data, Paras, Site)

# -- How many variables ?
print('Total number of variables :', Opti.nvar)

# === Runs ========================================
# -------------------
if Config.mode == 'calib_runs':
    # Calibration loop
    runs.calib_runs(Config, Opti, Data, Paras, Site)

elif Config.mode == 'forward_runs':
    # Calibration loop
    runs.forward_runs(Config, Opti, Data, Paras, Site, options)

elif Config.mode == 'sensi_morris' and Config.MSinit == 0:
    # Simulations when varying the parameters, Morris's one-at-a-time
    runs.morris_runs(Config, Opti, Data, Paras, Site)


# END MAIN
