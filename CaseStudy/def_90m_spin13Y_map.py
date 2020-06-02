# **************************************************************************
# CONFIGURATION FILE
# **************************************************************************

# import numpy as np
from datetime import datetime
# , timedelta

# ==========================================================================
#  Optimization characteristics
#
# --------------------------------------------------------------------------


class Opti:

    main_path = '/users/s01ks6/Data/'
    # -- Name of the output directory for the current assimilations
    #    if no name is provided, the default path will be that of this def file
    PATH_EXEC = 'setup_all'

    # Number of iterations
    nit = 1000

    # Take into account bare rock? (limiting soil evap)
    simRock = 0

# ==========================================================================
# Site conceptualization
# ------------------------------------------------------------------------


class Site:
    # -- Soil
    soils = ['VrtS', 'VrtD', 'FrlD', 'FrlS']
    ns = len(soils)
    sfiles = ['unit.soil_' + s + '.map' for s in soils]
    # sfiles = ['unit.map']
    # -- Vegetation
    vegs = ['Tree', 'Grass']
    nv = len(vegs)
    vfile = 'SpeciesParams.tab'

# ==============================================================================
# Define the measurements on which assimilation is performed as well as
# the error structures on the observations
# ------------------------------------------------------------------------


class Data:

    # obsdir = '/users/s01ks6/Data/BB/'

    # -- Observations used
    obs = {}

    # == Maps for spatialized values -----------------------------------------
    # Water budget
    obs['S_map.Srf'] = {'sim_file': 'Ponding_', 'conv': 1, 'type': 'mapTs'}
    obs['SWC_map.L1'] = {'sim_file': 'SWC1_', 'conv': 1, 'type': 'mapTs'}
    obs['SWC_map.L2'] = {'sim_file': 'SWC2_', 'conv': 1, 'type': 'mapTs'}
    obs['SWC_map.L3'] = {'sim_file': 'SWC3_', 'conv': 1, 'type': 'mapTs'}
    obs['SWC_map.Soil'] = {'sim_file': 'SWCav', 'conv': 1, 'type': 'mapTs'}
    obs['GWD_map'] = {'sim_file': 'WTD_', 'conv': 1, 'type': 'mapTs'}
    # Upward fluxes
    obs['Q_map.Et'] = {'sim_file': 'EvapT', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.EtL1'] = {'sim_file': 'EvTL1', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.EtL2'] = {'sim_file': 'EvTL2', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.EtL3'] = {'sim_file': 'EvTL3', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.Es'] = {'sim_file': 'EvapS', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.L1toSrf'] = {'sim_file': 'RSrf', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.L2toL1'] = {'sim_file': 'RL1', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.L3toL2'] = {'sim_file': 'RL2', 'conv': 1, 'type': 'mapTs'}
    # Downward fluxes
    obs['Q_map.SrftoL1'] = {'sim_file': 'Inf', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.L1toL2'] = {'sim_file': 'PrcL2', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.L2toL3'] = {'sim_file': 'PrcL3', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.Rech'] = {'sim_file': 'Rchg', 'conv': 1, 'type': 'mapTs'}
    # Lateral out
    obs['Q_map.ChnOut'] = {'sim_file': 'LChno', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.SrfOut'] = {'sim_file': 'LSrfo', 'conv': 1, 'type': 'mapTs'}
    obs['Q_map.GWOut'] = {'sim_file': 'LGWo', 'conv': 1, 'type': 'mapTs'}

    # Water age
    # obs['Age_map.L1'] = {'sim_file': 'AgesL1', 'conv': 1, 'type': 'mapTs'}
    # obs['Age_map.L2'] = {'sim_file': 'AgesL2', 'conv': 1, 'type': 'mapTs'}
    # obs['Age_map.L3'] = {'sim_file': 'AgesL3', 'conv': 1, 'type': 'mapTs'}
    obs['Age_map.Soil'] = {'sim_file': 'AgesAv', 'conv': 1/365.25,
                           'type': 'mapTs'}
    obs['Age_map.Et'] = {'sim_file': 'AgeeT', 'conv': 1/365.25,
                         'type': 'mapTs'}
    obs['Age_map.Et0'] = {'sim_file': 'AgEt0_', 'conv': 1/365.25,
                          'type': 'mapTs'}
    obs['Age_map.Et1'] = {'sim_file': 'AgEt1_', 'conv': 1/365.25,
                          'type': 'mapTs'}
    obs['Age_map.Es'] = {'sim_file': 'AgeeS', 'conv': 1/365.25,
                         'type': 'mapTs'}
    obs['Age_map.GW'] = {'sim_file': 'AgesL3', 'conv': 1/365.25,
                         'type': 'mapTs'}

    nobs = len(obs)

    # -- Simulations outputs

    # Starting date
    simbeg = datetime(2004,1,1)
    lsim = 10586
    # simt = [simbeg+timedelta(days=x) for x in range(lsim)]

    # Number of obsrvations points
    nts = 8
    # In case the tsmask.map is not ordered like Ech2O does:
    # ascending order left->right, then top->bottom
    sim_order = [1, 2, 3, 4, 5, 6, 7, 8]

# ==============================================================================
# Parameters that are optimized
# ------------------------------------------------------------------------------


class Paras:
    # -- Parameters to optimize: dimensions
    ref = {}

    ref['Manning']   = {'soil':0, 'veg':0, 'log':0, 'file':'chanmanningn'}
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
