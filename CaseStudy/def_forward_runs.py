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
    # Script mode (see ech2o_multitools.py for details)
    mode = 'forwards_runs'
    # EcH2O exectuable name
    exe = 'ech2o_iso'
    # Number of CPUs used (=1 if not specified)
    ncpu = 8
    # EcH2O's Config files (specified in job file)
    # cfg = 'config_90m_nospin_Ts.ini'  # Main one
    # cfgTrck = 'configTrck_Ts.ini'     # Tracking one
    # outdir = 'Calib_90m_nospin'   # Output directory

    # Trim outputs ? If so, the beginning and/or ending of outputs
    # should not be saved (e.g. transient state)
    trimB = 1  # Initial cutoff for time series (default=1)
    trimBmap = 1  # Initial cutoff for maps (default=1)
    # trimL =  # Length of saved outputs (default: full length from trimB)

    repBS = 1  # Report BasinSummary_runX.txt? (and tracking summaries)
    # (default =0)

# ==========================================================================
#  Optimization/calibration characteristics
#
# --------------------------------------------------------------------------


class Opti:

    # Parameters for DREAM calibration
    rep = 2000  # Maximum number of repetitions
    nChains = 6  # Number of chains
    convergence_limit = 1.2  # Gelman-Rubin convegence limit criterion
    runs_after_conv = 50  # To construct posterior distrib

    DREAMpar = 'seq'
    # Sampling sequence: 'seq' (default)-> normal iterations
    # 'mpc': multiprocessing on a single windows pc
    # 'mpi': parallel computing for unix (mac/linux/HPC)

# ==========================================================================
# Site conceptualization
# ------------------------------------------------------------------------


class Site:

    # Resolution, used to locate input directories (='50m' if not specified)
    Resol = '90m'
    # Enable water tracking (isotopes, ages, ...)
    isTrck = 1
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
    lsim = 10586  # length of the whole simulation (spinup included, days)

    # length of spinup (to be discarded WITHIN lsim, days)
    lspin = 909  # ONLY for calibration
    # for forward_runs, use the Config.trimB, trimL keywords
    # (TODO: merge those)

    # Starting date after spinup
    simbeg = datetime(2004, 1, 1)

    # -- Observations used
    # Number of obsrvations points
    nts = 10
    # In case the tsmask.map is not ordered as Ech2O does:
    # ascending order left->right, then top->bottom
    sim_order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    obs = {}

    # Hydrometric data
    obs['Streamflow'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 5,
                         'conv': 1, 'type': 'Ts'}
    # Groundwater level
    obs['GWD_P1'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':5,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P2'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':2,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P3'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':1,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P5'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':10,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P6'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':7,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P7'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':6,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P8'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':4,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P9'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':4,
                     'conv':1, 'type': 'Ts'}
    obs['GWD_P10'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':4,
                      'conv':1, 'type': 'Ts'}
    obs['GWD_P12'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':9,
                      'conv':1, 'type': 'Ts'}
    obs['GWD_P13'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts':8,
                      'conv':1, 'type': 'Ts'}
    
    # == Catchment-scale values: From BasinSummary.txt
    # & BasinAgeSummary.txt ------------
    
    # Water budget
    obs['P_tot'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 1,
                    'conv': 1, 'type': 'Total'}
    obs['S_tot.SWE'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 2,
                        'conv': 1, 'type': 'Total'}
    obs['S_tot.Veg'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 3,
                        'conv': 1, 'type': 'Total'}
    obs['S_tot.Srf'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 4,
                        'conv': 1, 'type': 'Total'}
    obs['S_tot.Soil'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 5,
                         'conv': 1, 'type': 'Total'}
    obs['S_tot.L1'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 6,
                       'conv': 1, 'type': 'Total'}
    obs['S_tot.L2'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 7,
                       'conv': 1, 'type': 'Total'}
    obs['S_tot.RZ'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 9,
                       'conv': 1, 'type': 'Total'}
    obs['S_tot.GW'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 10,
                       'conv': 1, 'type': 'Total'}
    obs['Q_tot.E'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 11,
                      'conv': 1, 'type': 'Total'}
    obs['Q_tot.Es'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 12,
                       'conv': 1, 'type': 'Total'}
    obs['Q_tot.Ei'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 13,
                       'conv': 1, 'type': 'Total'}
    obs['Q_tot.Et'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 14,
                       'conv': 1, 'type': 'Total'}
    obs['Q_tot.Srf'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 16,
                        'conv': 1, 'type': 'Total'}
    obs['Q_tot.GW'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 17,
                       'conv': 1, 'type': 'Total'}
    obs['Q_tot.srfQ'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 18,
                         'conv': 1, 'type': 'Total'}
    obs['Q_tot.GWQ'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 19,
                        'conv': 1, 'type': 'Total'}
    obs['Q_tot.Rch'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 20,
                        'conv': 1, 'type': 'Total'}
    # Lumped ages
    obs['Age_tot.S'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 1,
                        'conv': 1, 'type': 'Total'}
    obs['Age_tot.Soil'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 5,
                           'conv': 1, 'type': 'Total'}
    obs['Age_tot.L1'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 6,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.L2'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 7,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.RZ'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 9,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.GW'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 10,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.E'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 11,
                        'conv': 1, 'type': 'Total'}
    obs['Age_tot.Es'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 12,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.Ei'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 13,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.Et'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 14,
                         'conv': 1, 'type': 'Total'}
    obs['Age_tot.Q'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 16,
                        'conv': 1, 'type': 'Total'}
    obs['Age_tot.Qgw'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 17,
                          'conv': 1, 'type': 'Total'}
    obs['Age_tot.Out'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 18,
                          'conv': 1, 'type': 'Total'}
    obs['Age_tot.srfQ'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 19,
                           'conv': 1, 'type': 'Total'}
    obs['Age_tot.GWQ'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 20,
                          'conv': 1, 'type': 'Total'}
    obs['Age_tot.Rch'] = {'sim_file': 'BasinAgeSummary.txt', 'sim_pts': 21,
                          'conv': 1, 'type': 'Total'}

    nobs = len(obs)

# ==============================================================================
# Parameters that are optimized
# ------------------------------------------------------------------------------


class Paras:
    # -- Parameters to optimize: dimensions
    ref = {}

    ref['Manning'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'chanmanningn'}
    # - soil-dependent
    ref['DepthL2'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'soildepth.L2'}
    ref['Porosity0'] = {'soil': 1, 'veg': 0, 'log': 0, 'file': 'poros'}
    # ref['kPorosity'] = {'soil': 1, 'veg': 0, 'log': 0, 'file': 'kporos'}
    ref['PorosL2'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'poros.L2'}
    ref['PorosL3'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'poros.L3'}

    ref['Khoriz0'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'Khsat'}
    ref['Anisotropy'] = {'soil': 1, 'veg': 0, 'log': 0, 'file': 'KvKh'}
    ref['kKhoriz'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'kKhsat'}
    # Vegetation
    ref['p_Tree'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'p_0'}
    ref['p_Grass'] = {'soil': 0, 'veg': 0, 'log': 0, 'file': 'p_1'}
    ref['Kroot'] = {'soil': 0, 'veg': 1, 'log': 1}
    # water use
    ref['Gs_max'] = {'soil': 0, 'veg': 1, 'log': 1}
