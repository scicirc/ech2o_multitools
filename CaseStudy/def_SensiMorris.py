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
    mode = 'sensi_morris'
    # EcH2O exectuable name
    exe = 'ech2o-iso'
    # path to input/output directories
    PATH_CLIM = PATH_MAIN+'/Input_Climate'  # Clim inputs
    PATH_CFG = PATH_MAIN+'/Input_Configs'  # Ech2O's configs template dir
    PATH_SPA_REF = PATH_MAIN+'/Input_Maps_P5'  # Inputs maps / veg param
    PATH_PAR = PATH_MAIN+'/Sensitivity_ParameterSets'  # Inputs maps / veg param
    # Number of CPUs used (=1 if not specified)
    # ncpu = 6 # given as option in job using $SLURM_CPUS_PER_TASK
    # EcH2O's Config files
    # cfg = 'config_90m_nospin_Ts.ini'  # Main one
    # cfgTrck = 'configTrck_Ts.ini'     # Tracking one
    # outdir = 'Calib_90m_nospin'   # Output directory
    PATH_SCRATCH = '/workdir/kups'

# ==========================================================================
#  Optimization/calibration characteristics
#
# --------------------------------------------------------------------------


class Opti:

    # Type of metrics used for elementary effects
    # (normalised) rmse, bias, or differenc in std
    EEs = ['rmse','bias','dstd']
    # Walk in the parameter space for Morris
    MSspace = 'radial'#trajectory' 
    # Number of trajectories
    nr = 20

# ==========================================================================
# Site conceptualization
# ------------------------------------------------------------------------


class Site:

    # Enable water tracking (isotopes, ages, ...)
    isTrck = 1
    # -- Soil
    soils = []#'black', 'red']
    ns = len(soils)
    sfiles = ['unit.soil_' + s + '.map' for s in soils]
    # sfiles = ['unit.map']
    # -- Vegetation
    vegs = ['ATT'] #, 'ATTmid','ATTtop','Shorea','Swamp']
    #vegs = ['Tree', 'Grass']
    nv = len(vegs)
    vfile = 'SpeciesParams.tab'
    # Take into account bare rock? (limiting soil evap)
    simRock = 0
    # Initial soil moisture
    # absolute value (m3 m-3)
    # SWC1init = 0.2
    # SWC2init = 0.15 
    # SWC3init = 0.1 
    # proportion of porosity
    frac_SWC1 = 0.4 
    frac_SWC2 = 0.4 
    frac_SWC3 = 0.7   

# ==============================================================================
# Define the observations charactreristics, including sim outputs
# and the measurements on which assimilation is performed
# ------------------------------------------------------------------------


class Obs:

    # -- Simulations outputs

    # Starting date
    lsim = 7483
    # length of the whole simulation (days, here 2003-7-7 : 2023-12-31)

    # Sim outputs can be limited trimmed to only keep the desired
    # time span (e.g. shaving off spinup / transient periods before
    # for calibration, ensemble runs, etc.)
    saveB = 1  # Timestep to store times series (see class Obs, default = 1)
    saveBmap = 1  # Timestep to stroage maps (default=1)
    # saveL =  # Length of saved outputs (default: full length from saveB on)

    # Starting date of time series after trimming (trimB)
    simbeg = datetime(2003, 7, 7) # going until 2023-12-24

    # -- Observations used
    # Number of obsrvations points
    nts = 1
    # In case the tsmask.map is not ordered as Ech2O does:
    # ascending order left->right, then top->bottom
    sim_order = [1]

    # -- Observations used
    obs = {}

    # Groundwater level
    obs['GWD'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 1,
                  'sim_conv': 1, 'type': 'Ts'}
    obs['GWD.L3'] = {'sim_file': 'WaterTableDepth_L3.tab', 'sim_pts': 1,
                     'sim_conv': 1, 'type': 'Ts'}
    obs['Evap'] = {'sim_file': 'Evap.tab', 'sim_pts': 1,
                  'sim_conv': 1, 'type': 'Ts'}
    obs['EvapS'] = {'sim_file': 'EvapS.tab', 'sim_pts': 1,
                  'sim_conv': 1, 'type': 'Ts'}
    obs['EvapT'] = {'sim_file': 'EvapT.tab', 'sim_pts': 1,
                  'sim_conv': 1, 'type': 'Ts'}

    # Chloride concentration in stream and groundwater
    obs['cCl.WT'] = {'sim_file': 'cCl_watertable.tab', 'sim_pts': 1,
                     'sim_conv': 1, 'type': 'Ts'}
    obs['cCl.GW'] = {'sim_file': 'cCl_groundwater.tab', 'sim_pts': 1,
                     'sim_conv': 1, 'type': 'Ts'}
    obs['cCl.L3'] = {'sim_file': 'cCl_soilL3.tab', 'sim_pts': 1,
                     'sim_conv': 1, 'type': 'Ts'}



    nobs = len(obs)

# ==============================================================================
# Parameters that are optimized
# ------------------------------------------------------------------------------


class Paras:

    # -- Parameters to optimize: dimensions
    ref = {}

    # - soil-dependent
    # ref['Depth']     = {'min':40, 'max': 50, 'file':'soildepth'}
    # ref['DepthL1']   = {'soil':1, 'min': [2, 1, 1, 0.2], 'max': [6, 3, 3, 1], 'file': 'soildepth.L1'}
    # ref['DepthL2']   = {'min': 10, 'max': 35,  'file': 'soildepth.L2'}

    ref['PorosL1'] = {'min': 0.3, 'max': 0.5, 'file': 'poros'}
    # ref['PorosL2'] = {'min': 0.05, 'max': 0.5, 'file': 'poros.L2'}
    # ref['PorosL3'] = {'min': 0.005, 'max': 0.2, 'file': 'poros.L3'}

    ref['KhorizL1']  = {'min': 1e-7, 'max': 1e-2, 'log': 1, 'file': 'Khsat'}
    # ref['KhorizL2']  = {'min': 1e-9, 'max': 1e-5, 'log': 1, 'file': 'Khsat.L2'}
    # ref['KhorizL3']  = {'min': 1e-9, 'max': 1e-5, 'log': 1, 'file': 'Khsat.L3'}
    # #ref['Leakance']  = {'min': 1e-5, 'max': 1, 'log': 1, 'file':'leakance'}
    # ref['AnisoL1'] = {'min': 0.01, 'max': 100, 'log': 1, 'file': 'KvKh'}
    # ref['AnisoL2'] = {'min': 0.01, 'max': 100, 'log': 1, 'file': 'KvKh.L2'}
    # ref['AnisoL3'] = {'min': 0.01, 'max': 100, 'log': 1, 'file': 'KvKh.L3'}

    # ref['BClambdaL1'] = {'min': 2., 'max': 12.,  'file': 'BClambda'}
    # ref['BClambdaL2'] = {'min': 2., 'max': 12., 'file': 'BClambda.L2'}
    # ref['BClambdaL3'] = {'min': 2., 'max': 12., 'file': 'BClambda.L3'}
    # ref['PsiAEL1']    = {'min': 0.02, 'max': 0.8, 'file': 'psi_ae'}
    # ref['PsiAEL2']    = {'min': 0.02, 'max': 0.8, 'file': 'psi_ae.L2'}
    # ref['PsiAEL3']    = {'min': 0.02, 'max': 0.8, 'file': 'psi_ae.L3'}
    # ref['SMresL1']    = {'min': 0.005, 'max': 0.15, 'file':'theta_r'}
    # ref['SMresL2']    = {'min': 0.005, 'max': 0.1, 'file':'theta_r.L2'}
    # ref['SMresL3']    = {'min': 0.005, 'max': 0.05, 'file':'theta_r.L3'}

    # #ref['Manning']   = {'min': 5, 'max': 25, 'file': 'chanmanningn'}
    # ref['AlbedoSoil'] = {'min': 0.1, 'max': 0.3, 'file': 'albedo.soil'}
    
    # Vegetation (mostly ATTlow/mid/top, also Shorea params)
    #ref['NPP/GPP'] = {'veg':1, 'min': 0.4, 'max': 0.6}
    ref['Gs_max'] = {'veg': 1, 'min': 0.005, 'max': 0.03}
    #ref['CnpQEff']   = {'veg':1, 'min': 1e-6,'max': 5e-6}
    
    # ref['Topt']      = {'veg': 1, 'min': 15., 'max': 30.}
    # ref['Tmax']      = {'veg': 1, 'min': 31., 'max': 45.}
    # ref['Tmin']      = {'veg': 1, 'min': 0.,'max': 9.}
    # # ref['AllocLeaf_a']      = {'veg': 1, 'min': 0.,'max': 9.}
    # # ref['AllocLeaf_b']      = {'veg': 1, 'min': 0.,'max': 9.}
    # # ref['AllocStem_a']      = {'veg': 1, 'min': 0.,'max': 9.}
    # # ref['AllocStem_b']      = {'veg': 1, 'min': 0.,'max': 9.}
    # ref['Gs_light'] = {'veg': 1, 'min': 0.,'max': 300.}
    # ref['Gs_vpd'] = {'veg': 1, 'log': 1, 'min': 1e-6,'max': 3e-3}
    # ref['LWP_low'] = {'veg': 1, 'min': 100.,'max': 600.}
    # ref['LWP_high'] = {'veg': 1, 'min': 3.5, 'max': 80}
    # #ref['WiltPnt'] = {'veg':1, 'min': 120., 'max': 250.}
    # ref['CWS_max'] = {'veg': 1, 'min': 1e-5, 'max': 5e-4}

    # ref['KBeers'] = {'veg': 1, 'min': 0.3,'max': 0.7}
    # #ref['WUE_cnp'] = {'veg':1, 'min': 900.,'max': 1500.}
    # ref['Albedo'] = {'veg': 1, 'min': 0.05, 'max': 0.3}
    # # ref['Kroot'] = {'veg': 1, 'log': 1, 'min': 0.05, 'max': 0.3}     
