# **************************************************************************
# CONFIGURATION FILE
# **************************************************************************

import os
# import numpy as np
from datetime import datetime
# , timedelta

# ==========================================================================
#  General configuration
#
# --------------------------------------------------------------------------


class Config:

    # -- Name of the root working directory
    # if no name is provided, the default path will be where
    # the script is launched
    PATH_MAIN = os.getcwd()+'/'
    # Script mode
    mode = 'calib_SPOTPY'
    # EcH2O exectuable name
    exe = 'ech2o_iso'
    # path to input/output directories
    PATH_CLIM = PATH_MAIN+'/Input_Climate'  # Clim inputs
    PATH_CFG = PATH_MAIN+'/Input_Configs'  # Ech2O's configs template dir
    PATH_SPA_REF = PATH_MAIN+'/Input_Maps_90m'  # Inputs maps / veg param
    PATH_OBS = PATH_MAIN+'/Calibration_Datasets'  # Observation for calib
    # Number of CPUs used (=1 if not specified)
    ncpu = 2 # can also be given as option in job 
    # (e.g., using $SLURM_CPUS_PER_TASK on clusters using slurms)
    # EcH2O's Config files (not used for now, still specified in job file)
    # cfg = 'config_90m_nospin_Ts.ini'  # Main one
    # cfgTrck = 'configTrck_Ts.ini'     # Tracking one
    # outdir = 'Calib_90m_nospin'   # Output directory

# ==========================================================================
#  Optimization/calibration characteristics
#
# --------------------------------------------------------------------------


class Opti:

    SPOTalgo = 'DREAM'
    # Parameters for DREAM calibration
    rep = 2000  # Maximum number of repetitions
    nChains = 6  # Number of chains
    convergence_limit = 1.2  # Gelman-Rubin convegence limit criterion
    runs_after_conv = 50  # To construct posterior distrib

    SPOTpar = 'mpi'
    # Sampling sequence: 'seq' (default)-> normal iterations
    # 'mpc': multiprocessing on a single windows pc
    # 'mpi': parallel computing for unix (mac/linux/HPC)

    # Initial parameter sampling
    # if uniform, min and max have to be given
    # if normal, guess is needed, and stddev (taken as 0.1*(max-min) otherwise)
    initSample = 'uniform'

# ==========================================================================
# Site conceptualization
# ------------------------------------------------------------------------


class Site:

    # Resolution, used to locate input directories (='50m' if not specified)
    Resol = '90m'
    # Enable water tracking (isotopes, ages, ...)
    isTrck = 0
    # -- Soil
    soils = ['VrtS', 'VrtD', 'FrlD', 'FrlS']
    ns = len(soils)
    sfiles = ['unit.soil_' + s + '.map' for s in soils]
    # sfiles = ['unit.map']
    # -- Vegetation
    vegs = ['Tree', 'Grass']
    nv = len(vegs)
    vfile = 'SpeciesParams.tab'
    # Take into account bare rock? (limiting soil evap)
    simRock = 0

# ==============================================================================
# Define the measurements on which assimilation is performed
# ------------------------------------------------------------------------


class Data:

    # -- Simulations outputs
    # Starting date
    lsim = 1274  # length of the simulation bit to look at (days)

    # length of spinup (to be discarded WITHIN lsim, days)
    lspin = 909  # ONLY for calibration
    # for forward_runs, use the Config.trimB, trimL keywords
    # (TODO: merge those)

    # Starting date after spinup
    simbeg = datetime(2006, 1, 1) # going until 2011-12-31

    # -- Observations used
    # Number of obsrvations points
    nts = 10
    # In case the tsmask.map is not ordered as Ech2O does:
    # ascending order left->right, then top->bottom
    sim_order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    obs = {}
    # Hydrometric data
    obs['Streamflow'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 5,
                         'sim_conv': 1,
                         'obs_file': 'Discharge_grouped_daily.csv',
p                         #'fit_beg': datetime(2006, 1, 1),
                         #'fit_end': datetime(2006, 12, 31),
                         'obs_col': 2, 'obs_conv': 1,
                         'type': 'Ts'}
    # Groundwater level
    obs['GWD_P1'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 5,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_grouped_daily.csv',
                     'fit_beg': 32, 'obs_col': 2,
                     'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P3'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 1,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 3,
                     'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P5'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 10,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 4,
                     'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P6'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 7,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 5,
                     'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P7'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 6,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 6,
                     'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P9'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 4,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 7,
                     'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P10'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 4,
                      'sim_conv': 1,
                      'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 8,
                      'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P12'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 9,
                      'sim_conv': 1,
                      'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 9,
                      'obs_conv': -1, 'type': 'Ts'}
    obs['GWD_P13'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 8,
                      'sim_conv': 1,
                      'obs_file': 'Piezo_grouped_daily.csv', 'obs_col': 10,
                      'obs_conv': -1, 'type': 'Ts'}

    nobs = len(obs)

# ==============================================================================
# Parameters that are optimized
# ------------------------------------------------------------------------------


class Paras:
    # -- Parameters to optimize: dimensions
    ref = {}

    # ref['Manning'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'chanmanningn'}
    # - soil-dependent
    ref['DepthL2'] = {'soil': 0, 'veg': 0, 'log': 0,
                      'min': 5, 'max': 30, 'file': 'soildepth.L2'}
    ref['Porosity0'] = {'soil': 0, 'veg': 0, 'log': 0,
                        'min': 0.25, 'max': 0.4, 'file': 'poros'}
    # ref['kPorosity'] = {'soil': 1, 'veg': 0, 'log': 0, 'file': 'kporos'}
    ref['PorosL2'] = {'soil': 0, 'veg': 0, 'log': 0,
                      'min': 0.1, 'max': 0.25, 'file': 'poros.L2'}
    ref['PorosL3'] = {'soil': 0, 'veg': 0, 'log': 0,
                      'min': 0.01, 'max': 0.05, 'file': 'poros.L3'}

    ref['Khoriz0'] = {'soil': 0, 'veg': 0, 'log': 1,
                      'min': 1e-5, 'max': 1e-3, 'file': 'Khsat'}
    ref['Anisotropy'] = {'soil': 0, 'veg': 0, 'log': 1,
                         'min': 0.001, 'max': 0.1, 'file': 'KvKh'}
    # ref['kKhoriz'] = {'soil': 0, 'veg': 0, 'log': 0,
    #                   'min': 5, 'max': 30, 'file': 'kKhsat'}
    # Vegetation
    # ref['p_Tree'] = {'soil': 0, 'veg': 0, 'log': 0,
    #                  'min': 5, 'max': 30, 'file': 'p_0'}
    # ref['p_Grass'] = {'soil': 0, 'veg': 0, 'log': 0,
    #                   'min': 5, 'max': 30, 'file': 'p_1'}
    ref['Kroot'] = {'soil': 0, 'veg': 1, 'log': 1,
                    'min': [0.01, 0.5], 'max': [0.5, 10]}
    # water use
    ref['Gs_max'] = {'soil': 0, 'veg': 1, 'log': 1,
                     'min': [0.001, 0.001], 'max': [0.1, 0.1]}
