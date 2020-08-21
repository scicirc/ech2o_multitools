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
    ncpu = 8
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
    SPOTdb = 'DREAM'  # Output mode for SPOTPY, custom offers more
    # flexibility (separate files in multi-objective, add SPOTPY diagnostics...)
    # standard options: 'csv', 'hdf5', 'ram', 'sql', 'noData'
    SPOTpar = 'seq'
    # Sampling sequence: 'seq' (default)-> normal iterations
    # 'mpc': multiprocessing on a single windows pc
    # 'mpi': parallel computing for unix (mac/linux/HPC)

    # Parameters for DREAM calibration
    rep = 5000  # Maximum number of repetitions
    nChains = 20  # Number of chains
    convergence_limit = 1.2  # Gelman-Rubin convegence limit criterion
    runs_after_conv = 100  # To construct posterior distrib

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
# Define the observations charactreristics, including sim outputs and 
# the measurements on which assimilation is performed
# ------------------------------------------------------------------------


class Obs:

    # -- Simulations outputs

    # Starting date
    lsim = 1274  # length of the whole simulation (days)
    lspin = 909  # length of spinup (to be discarded WITHIN lsim, days)
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
                         'obs_file': 'Discharge_outlet.csv',
                         #'fit_beg': datetime(2006, 1, 1),
                         #'fit_end': datetime(2006, 12, 31),
                         'obs_col': 2, 'obs_conv': 1,
                         'type': 'Ts'}
    # Groundwater depth
    obs['GWD_P1'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 5,
                     'sim_conv': 1,
                     'obs_file': 'Piezo_well1.csv',
                     'obs_col': 2,
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
    ref['DepthL2'] = {'min': 10, 'max': 30, 'guess': 20, 
                      'file': 'soildepth.L2'}
    ref['Porosity0'] = {'min': 0.2, 'max': 0.4, 'guess': 0.3,
                        'file': 'poros'}
    # ref['kPorosity'] = {'soil': 1, 'veg': 0, 'log': 0, 'file': 'kporos'}
    ref['PorosL2'] = {'min': 0.1, 'max': 0.25, 'guess': 0.2,
                      'file': 'poros.L2'}
    ref['PorosL3'] = {'min': 0.005, 'max': 0.05, 'guess': 0.01,
                      'file': 'poros.L3'}
    ref['Khoriz0'] = {'min': 1e-6, 'max': 1e-3, 'log': 0, 'guess':1e-4,
                      'file': 'Khsat'}
    ref['Anisotropy'] = {'min': 0.001, 'max': 1, 'guess':0.1,
                         'file': 'KvKh'}
    ref['BClambda'] = {'min': 2, 'max': 10, 'guess': 8.8,
                       'file': 'BClambda'}
    ref['PsiAE'] = {'min': 0.05, 'max': 0.8, 'guess': 0.23,
                    'file': 'psi_ae'}
    # ref['kKhoriz'] = {'soil': 0, 'veg': 0, 'log': 0,
    #                   'min': 5, 'max': 30, 'file': 'kKhsat'}
    # Vegetation
    # ref['p_Tree'] = {'soil': 0, 'veg': 0, 'log': 0,
    #                  'min': 5, 'max': 30, 'file': 'p_0'}
    # ref['p_Grass'] = {'soil': 0, 'veg': 0, 'log': 0,
    #                   'min': 5, 'max': 30, 'file': 'p_1'}
    ref['Kroot'] = {'veg': [1, 0], 'log': 0, # Only Trees
                    'min': 0.05, 'max': 0.5, 'guess': 0.1}
    # water use
    ref['Gs_max'] = {'veg': [1, 0], 'log': 0, # Only Trees
                     'min': 0.005, 'max': 0.08, 'guess': 0.01}
