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
    exe = 'ech2o-iso_tmp'
    # path to input/output directories
    PATH_CLIM = PATH_MAIN+'/Input_Climate'  # Clim inputs
    PATH_CFG = PATH_MAIN+'/Input_Configs'  # Ech2O's configs template dir
    PATH_SPA_REF = PATH_MAIN+'/Input_Maps_300m'  # Inputs maps / veg param
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

    # Walk in the parameter space for Morris
    MSspace = 'radial'#trajectory' 
    # Number of trajectories
    nr = 10
    
# ==========================================================================
# Site conceptualization
# ------------------------------------------------------------------------


class Site:

    # Enable water tracking (isotopes, ages, ...)
    isTrck = 1
    # -- Soil
    maps1 = ['loam', 'sandyloam']
    nmap1 = len(maps1)
    map1files = ['soil_' + s + '.map' for s in maps1]
    # -- Geol
    maps2 = ['mancoshale', 'igneous','sandstone',
             'maroon', 'deposits']             
    nmap2 = len(maps2)
    map2files = ['geol_' + g + '.map' for g in maps2]
    # -- Land cover classes
    maps3 = ['grassripa', 'forest','barren']
    nmap3 = len(maps3)
    map3files = ['lc_' + l + '.map' for l in maps3]
    # -- Vegetation
    vegs = ['Grass','Decidous','Needleleaf','Riparian']
    #vegs = ['Tree', 'Grass']
    nv = len(vegs)
    vfile = 'SpeciesParams_LAI_Eqw.tab'
    # Take into account bare rock? (limiting soil evap)
    simRock = 0
    # Initial soil moisture
    # absolute value (m3 m-3)
    #SWC1init = 0.4
    #SWC2init = 0.2 
    #SWC3init = 0.1 
    # proportion of porosity
    frac_SWC1 = 0.9 
    frac_SWC2 = 0.9 
    frac_SWC3 = 0.9

    

# ==============================================================================
# Define the observations charactreristics, including sim outputs
# and the measurements on which assimilation is performed
# ------------------------------------------------------------------------


class Obs:

    # -- Simulations outputs

    # Starting date
    lsim = 2921#4382  # length of the whole simulation (spinup included, days)

    # Sim outputs can be limited trimmed to only keep the desired
    # time span (e.g. shaving off spinup / transient periods before
    # for calibration, ensemble runs, etc.)
    saveB = 731#2192  # Timestep to store times series (see class Obs, default = 1)
    saveBmap = 1  # Timestep to store maps (default=1)
    # saveL =  # Length of saved outputs (default: full length from saveB on)
    # Starting date of time series after trimming (trimB)
    simbeg = datetime(2014, 10, 1) # going until 2020-09-29

    # -- Observations used
    # Number of obsrvations points
    nts = 29
    # In case the tsmask.map is not ordered as Ech2O does:
    # ascending order left->right, then top->bottom
    sim_order = [i for i in range(1,nts+1)]

    # -- Observations used
    obs = {}

    # Snow water equivalent
    obs['SWE_Schofield'] = {'sim_file': 'SWE.tab', 'sim_pts': 1,
                            'sim_conv': 1, 'type': 'Ts'}
    # CAREFUL !! the map name should include spinup time steps (+730 for 2Y)
    obs['SWEmap_20160404'] = {'sim_file': 'SWE00001.282', 'type': 'mapStep'}
    obs['SWEmap_20180330'] = {'sim_file': 'SWE00002.007', 'type': 'mapStep'}
    obs['SWEmap_20180524'] = {'sim_file': 'SWE00002.062', 'type': 'mapStep'}
    obs['SWEmap_20190407'] = {'sim_file': 'SWE00002.380', 'type': 'mapStep'}
    # Basin-scale, using BasinSummary
    obs['fSnow'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 29,
                    'type': 'Total'}
    # + convert from m3 to mm using area in m2 (~86km2, see DEM.map)
    # (cumulative->.d-1 done in the script)
    obs['ET_total'] = {'sim_file': 'BasinSummary.txt', 'sim_pts': 14,
                       'sim_conv': 1000./86019299, 'type': 'Total'}
    # Stream discharge
    obs['Q_ERaboveQuigley'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 2,
                               'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Quigley'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 3,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Rustlers'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 4,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Bradley'] = {'sim_file': 'Streamflow.tab', 'sim_pts':5,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Rock'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 6,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Gothic'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 7,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Marmot'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 9,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_Copper'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 13,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['Q_ERbelowCopper'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 17,
                              'sim_conv': 1, 'type': 'Ts'}
    obs['Q_PumpHouse'] = {'sim_file': 'Streamflow.tab', 'sim_pts': 27,
                       'sim_conv': 1, 'type': 'Ts'}
    # Groundwater level
    obs['VWC.L1_SMN1-PHS1'] = {'sim_file': 'SoilMoistureL1.tab', 'sim_pts': 29,
                               'sim_conv': 1, 'type': 'Ts'}
    obs['VWC.L1_SMN3-4-5-PHS2'] = {'sim_file': 'SoilMoistureL1.tab', 'sim_pts': 27,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['VWC.L1_PHS4'] = {'sim_file': 'SoilMoistureL1.tab', 'sim_pts': 28,
                              'sim_conv': 1, 'type': 'Ts'}
    obs['VWC.L2_SMN1-PHS1'] = {'sim_file': 'SoilMoistureL2.tab', 'sim_pts': 29,
                          'sim_conv': 1, 'type': 'Ts'}
    obs['VWC.L2_SMN3-4-5-PHS2'] = {'sim_file': 'SoilMoistureL2.tab', 'sim_pts': 27,
                              'sim_conv': 1, 'type': 'Ts'}
    obs['VWC.L2_PHS4'] = {'sim_file': 'SoilMoistureL2.tab', 'sim_pts': 28,
                          'sim_conv': 1, 'type': 'Ts'}
    # Groundwater level
    obs['GWD_CPA-DOE-DOW-PLM6'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 27,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['GWD_PLM1'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 29,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['GWD_UPE-UPW'] = {'sim_file': 'WaterTableDepth.tab', 'sim_pts': 24,
                          'sim_conv': 1, 'type': 'Ts'}

    # d18O (stream)
    obs['d18O.Q_ERaboveQuigley'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 2,
                               'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_Quigley'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 3,
                       'sim_conv': 1, 'type': 'Ts'}
    # obs['d18O.Q_Rustlers'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 4,
    #                    'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_Bradley'] = {'sim_file': 'd18O_channel.tab', 'sim_pts':5,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_Rock'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 6,
                       'sim_conv': 1, 'type': 'Ts'}
    # obs['d18O.Q_Gothic'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 7,
    #                    'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_Marmot'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 9,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_Copper'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 13,
                       'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_ERbelowCopper'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 17,
                              'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Q_PumpHouse'] = {'sim_file': 'd18O_channel.tab', 'sim_pts': 27,
                       'sim_conv': 1, 'type': 'Ts'}
    # d18O (soil)
    obs['d18O.L1_Snodgrass1'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 19,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Snodgrass2'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 20,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Snodgrass3'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 21,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Snodgrass4'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 22,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Snodgrass5'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 23,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Snodgrass6'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 25,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Copper1'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 14,
                             'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L1_Copper2'] = {'sim_file': 'd18O_soilL1.tab', 'sim_pts': 16,
                             'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.L2_Copper2'] = {'sim_file': 'd18O_soilL2.tab', 'sim_pts': 16,
                             'sim_conv': 1, 'type': 'Ts'}
    # d18O (groundwater)
    obs['d18O.GW_GLS'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 8,
                    'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_Billys'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 11,
                             'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_RMBL'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 12,
                           'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_GUM'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 15,
                           'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_Inouye'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 16,
                             'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_Shumway-Tuttle'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 18,
                                     'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_PLM6-8'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 27,
                             'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.GW_PLM1-7'] = {'sim_file': 'd18O_groundwater.tab', 'sim_pts': 29,
                             'sim_conv': 1, 'type': 'Ts'}
    # d18O (xylem)
    obs['d18O.Deci_Snodgrass1'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 19,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Snodgrass1'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 19,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Snodgrass2'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 20,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Snodgrass2'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 20,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Snodgrass3'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 21,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Snodgrass3'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 21,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Snodgrass4'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 22,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Snodgrass4'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 22,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Snodgrass5'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 23,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Snodgrass5'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 23,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Snodgrass6'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 25,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Snodgrass6'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 25,
                                   'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Copper1'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 14,
                               'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Copper1'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 14,
                               'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Deci_Copper2'] = {'sim_file': 'd18OevapT_1.tab', 'sim_pts': 16,
                               'sim_conv': 1, 'type': 'Ts'}
    obs['d18O.Ever_Copper2'] = {'sim_file': 'd18OevapT_2.tab', 'sim_pts': 16,
                               'sim_conv': 1, 'type': 'Ts'}

    # # d2H (stream)
    # obs['d2H.Q_ERaboveQuigley'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 2,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_Quigley'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 3,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['d2H.Q_Rustlers'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 4,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_Bradley'] = {'sim_file': 'd2H_channel.tab', 'sim_pts':5,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_Rock'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 6,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['d2H.Q_Gothic'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 7,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_Marmot'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 9,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_Copper'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 13,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_ERbelowCopper'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 17,
    #                           'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Q_PumpHouse'] = {'sim_file': 'd2H_channel.tab', 'sim_pts': 27,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # d2H (soil)
    # obs['d2H.L1_Snodgrass1'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 19,
    #                              'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Snodgrass2'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 20,
    #                              'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Snodgrass3'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 21,
    #                              'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Snodgrass4'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 22,
    #                              'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Snodgrass5'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 23,
    #                              'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Snodgrass6'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 25,
    #                              'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Copper1'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 14,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L1_Copper2'] = {'sim_file': 'd2H_soilL1.tab', 'sim_pts': 16,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.L2_Copper2'] = {'sim_file': 'd2H_soilL2.tab', 'sim_pts': 16,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # # d2H (groundwater)
    # obs['d2H.GW_GLS'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 8,
    #                 'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_Billys'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 11,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_RMBL'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 12,
    #                        'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_GUM'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 15,
    #                        'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_Inouye'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 16,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_Shumway-Tuttle'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 18,
    #                                  'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_PLM6-8'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 27,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.GW_PLM1-7'] = {'sim_file': 'd2H_groundwater.tab', 'sim_pts': 29,
    #                          'sim_conv': 1, 'type': 'Ts'}
    # # d2H (xylem)
    # obs['d2H.Deci_Snodgrass1'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 19,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Snodgrass1'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 19,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Snodgrass2'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 20,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Snodgrass2'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 20,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Snodgrass3'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 21,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Snodgrass3'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 21,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Snodgrass4'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 22,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Snodgrass4'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 22,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Snodgrass5'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 23,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Snodgrass5'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 23,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Snodgrass6'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 25,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Snodgrass6'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 25,
    #                                'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Copper1'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 14,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Copper1'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 14,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Deci_Copper2'] = {'sim_file': 'd2HevapT_1.tab', 'sim_pts': 16,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['d2H.Ever_Copper2'] = {'sim_file': 'd2HevapT_2.tab', 'sim_pts': 16,
    #                            'sim_conv': 1, 'type': 'Ts'}
    
    # Energy balance
    obs['ET_PHFluxTower'] = {'sim_file': 'Evap.tab', 'sim_pts': 27,
                             'sim_conv': 8.64e7, 'type': 'Ts'}
    obs['NetRad_PHFluxTower'] = {'sim_file': 'NetRad_tot.tab', 'sim_pts': 27,
                                 'sim_conv': 1, 'type': 'Ts'}
    obs['LatHeat_PHFluxTower'] = {'sim_file': 'Evap.tab', 'sim_pts': 27,
                                  'sim_conv': -1, 'type': 'Ts'}
    obs['SensHeat_PHFluxTower'] = {'sim_file': 'Evap.tab', 'sim_pts': 27,
                                  'sim_conv': -1, 'type': 'Ts'}
    obs['Tsoil_SMN1-PHS1'] = {'sim_file': 'SoilTemp.tab', 'sim_pts': 29,
                              'sim_conv': 1, 'type': 'Ts'}
    obs['Tsoil_SMN3-4-5-PHS2'] = {'sim_file': 'SoilTemp.tab', 'sim_pts': 27,
                                  'sim_conv': 1, 'type': 'Ts'}
    obs['Tsoil_PHS4'] = {'sim_file': 'SoilTemp.tab', 'sim_pts': 28,
                         'sim_conv': 1, 'type': 'Ts'}
    obs['NetRad_PHS1'] = {'sim_file': 'NetRad_tot.tab', 'sim_pts': 29,
                          'sim_conv': 1, 'type': 'Ts'}
    obs['NetRad_PHS2'] = {'sim_file': 'NetRad_tot.tab', 'sim_pts': 27,
                          'sim_conv': 1, 'type': 'Ts'}
    obs['NetRad_PHS4'] = {'sim_file': 'NetRad_tot.tab', 'sim_pts': 28,
                          'sim_conv': 1, 'type': 'Ts'}
    obs['GrndHeat_PHS1'] = {'sim_file': 'GrndHeat.tab', 'sim_pts': 29,
                            'sim_conv': -1, 'type': 'Ts'}
    obs['GrndHeat_PHS2'] = {'sim_file': 'GrndHeat.tab', 'sim_pts': 27,
                            'sim_conv': -1, 'type': 'Ts'}
    obs['GrndHeat_PHS4'] = {'sim_file': 'GrndHeat.tab', 'sim_pts': 28,
                            'sim_conv': -1, 'type': 'Ts'}

    # Age (stream)
    # obs['Age.Q_ERaboveQuigley'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 2,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_Quigley'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 3,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['Age.Q_Rustlers'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 4,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_Bradley'] = {'sim_file': 'Age_channel.tab', 'sim_pts':5,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_Rock'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 6,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['Age.Q_Gothic'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 7,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_Marmot'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 9,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_Copper'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 13,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_ERbelowCopper'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 17,
    #                           'sim_conv': 1, 'type': 'Ts'}
    # obs['Age.Q_PumpHouse'] = {'sim_file': 'Age_channel.tab', 'sim_pts': 27,
    #                    'sim_conv': 1, 'type': 'Ts'}

    # # fMelt (stream)
    # obs['fMelt.Q_ERaboveQuigley'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 2,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_Quigley'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 3,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['fMelt.Q_Rustlers'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 4,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_Bradley'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts':5,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_Rock'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 6,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['fMelt.Q_Gothic'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 7,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_Marmot'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 9,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_Copper'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 13,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_ERbelowCopper'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 17,
    #                           'sim_conv': 1, 'type': 'Ts'}
    # obs['fMelt.Q_PumpHouse'] = {'sim_file': 'fMelt_channel.tab', 'sim_pts': 27,
    #                    'sim_conv': 1, 'type': 'Ts'}

    # # fGWlat (stream)
    # obs['fGWlat.Q_ERaboveQuigley'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 2,
    #                            'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_Quigley'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 3,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['fGWlat.Q_Rustlers'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 4,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_Bradley'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts':5,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_Rock'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 6,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # # obs['fGWlat.Q_Gothic'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 7,
    # #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_Marmot'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 9,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_Copper'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 13,
    #                    'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_ERbelowCopper'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 17,
    #                           'sim_conv': 1, 'type': 'Ts'}
    # obs['fGWlat.Q_PumpHouse'] = {'sim_file': 'fGWlat_channel.tab', 'sim_pts': 27,
    #                    'sim_conv': 1, 'type': 'Ts'}
    
    nobs = len(obs)

# ==============================================================================
# Parameters that are optimized
# ------------------------------------------------------------------------------


class Paras:

    # -- Parameters to optimize: dimensions
    ref = {}

    # - soil-dependent
    # ref['Depth']     = {'min':40, 'max': 50, 'file':'soildepth'}
    ref['DepthL1']   = {'map1': 0, 'file': 'soildepth.L1',
                        'min': 0.05, 'max': 0.3}
    ref['Depth']   = {'map2':0, 'file': 'soildepth',
                      'min': 30., 'max': 150.}

    ref['PorosL1-L2'] = {'map1': 1,  'file': 'poros',
                      'min': [0.3, 0.3], 'max': [0.6, 0.6]}
    # ref['PorosL2'] = {'map1': 1, 'file': 'poros.L2',
    #                   'min': [0.05, 0.05], 'max': [0.5, 0.5]}
    ref['PorosL3'] = {'map2': 1,'file': 'poros.L3',
                      'min': [0.01, 0.01, 0.01, 0.01, 0.01],
                      'max': [0.3, 0.3, 0.3, 0.3, 0.3]}
                      
    ref['KhorizL1-L2']  = {'map1': 1, 'file': 'Khsat',
                           'min': [1e-7, 1e-7], 'max': [0.01, 0.01], 'log': 1}
    # ref['KhorizL2']  = {'map1': 1, 'file': 'Khsat.L2',
    #                     'min': [1e-7, 1e-7], 'max': [0.01, 0.01], 'log': 1}
    ref['KhorizL3']  = {'map2': 1, 'file': 'Khsat.L3', 'log': 1,
                        'min': [1e-9, 1e-9, 1e-9, 1e-9, 1e-9],
                        'max':[1e-4, 1e-4, 1e-4, 1e-4, 1e-4]}
    # ref['Leakance']  = {'map2': 1, 'file': 'leakance', 'log':1,
    #                     'min': [1e-12, 1e-12, 1e-12, 1e-12, 1e-12],
    #                     'max': [1e-3, 1e-3, 1e-3, 1e-3, 1e-3]}
    
    ref['AnisoL1-L2'] = {'map1':1 , 'file': 'KvKh', 'log': 1,
                         'min': [0.01, 0.01], 'max': [100, 100]}
    # ref['AnisotropyL2'] = {'map1':1 , 'file': 'KvKh',
    #                        'min': [0.1, 0.1], 'max': [10, 10]}
    # ref['AnisotropyL3'] = {'map2':1 , 'file': 'KvKh.L3',
    #                        'min': [0.01, 0.01, 0.01, 0.01, 0.01],
    #                        'max': [100, 100, 100, 100, 100]}

    ref['BClambdaL1-L2'] = {'map1': 1, 'file': 'BClambda',
                            'min': [2., 2.], 'max': [12., 12.]}
    # ref['BClambdaL2'] = {'map1': 1, 'file': 'BClambda.L2',
    #                      'min': [2., 2.], 'max': [12., 12.]}
    ref['BClambdaL3'] = {'map2': 1, 'file': 'BClambda.L3',
                         'min': [2., 2., 2., 2., 2.], 'max': [12.,12.,12.,12.,12.]}

    ref['PsiAEL1-L2']    = {'map1': 1, 'file': 'psi_ae',
                         'min': [0.1,0.1], 'max': [0.8, 0.8]}
    # ref['PsiAEL2']    = {'map1': 1, 'file': 'psi_ae.L2',
    #                      'min': [0.1,0.1], 'max': [0.8, 0.8]}
    ref['PsiAEL3']    = {'map2': 1, 'file': 'psi_ae.L3',
                         'min': [0.1,0.1,0.1,0.1,0.1], 'max': [0.8,0.8,0.8,0.8,0.8]}
    
    ref['SMresL1-L2']    = {'map1': 1, 'file': 'theta_r',
                            'min':[0.03,0.03], 'max': [0.15,0.15]}
    # ref['SMresL2']    = {'map1': 1, 'file': 'theta_r.L2',
    #                      'min':[0.03,0.03], 'max': [0.15,0.15]}
    ref['SMresL3']    = {'map2': 1, 'file':'theta_r.L3',
                         'min': [0.005,0.005,0.005,0.005,0.005],
                         'max': [0.1,0.1,0.1,0.1,0.1]}

    ref['GWseepL1-L2'] = {'file':'chanparam', 'log':1,
                       'min': 1e-5, 'max': 10}
    ref['GWseepL3'] = {'file':'chanparam.L3', 'log':1,
                       'min': 1e-5, 'max': 10}
    
    ref['Manning']   = {'file': 'chanmanningn', 'min': 5, 'max': 25}
    ref['AlbedoSoil'] = {'map1': 0, 'file': 'albedo.soil',
                         'min': 0.1, 'max': 0.3}
    ref['SnowmeltCoeff'] = {'map3': 1, 'file': 'snowmeltCoeff', 'log': 1,
                            'min': [1e-9,1e-9,1e-9], 'max': [1e-6,1e-6,1e-6]}
    ref['SnowRainTemp'] = {'map3':1, 'file': 'snowrainTemp',
                           'min': [-1., -1., -1.], 'max': [3., 3., 3.]}
    
    # -- Vegetation parameters
    # Reminder: species0 = grass, species1= deciduous (aspen)
    # species2 = evergreen, species3 = ~deciduous (riparian shrub)
    # ref['NPP/GPP'] = {'veg':1, 'min':[0.35, 0.35, 0.35, 0.35],
    #                   'max':[0.65, 0.65, 0.65, 0.65]}
    ref['Gs_max'] = {'veg': 1, 'min': [0.001, 0.001, 0.001, 0.001],
                     'max': [0.035, 0.035, 0.035, 0.035]}
    # WUE fixed to 1000 gC.m-1 -> quantum efficiency sensitivity
    # test canopy efficiency (quantum*water use), with unit
    # ref['CnpQEff']   = {'veg':1, 'min': [1.5e-6, 1.5e-6, 1.5e-6, 1.5e-6],
    #                     'max': [1.5e-5, 1.5e-5, 1.5e-5, 1.5e-5]}
    ref['Topt']   = {'veg': 1, 'min': [15., 15., 15., 15.],
                     'max': [30., 30., 30., 30.]}
    ref['Tmax']   = {'veg': 1, 'min': [31., 31. , 31., 31.],
                     'max': [45., 45., 45., 45.]}
    ref['Tmin']   = {'veg': 1, 'min': [0., 0., 0., 0.],
                     'max': [10., 10., 10., 10.]}
    # AllocLeaf_a: coeff for NPP to foliage (evergreen), not used (grass)
    # minimum root alloc factor (veg=2)
    # ref['AllocLeaf_a'] = {'veg': [0, 1, 2, 3],
    #                       'min': [0.4, 1.5, 0.4],'max': [0.9, 2.5, 0.9]}
    # AllocLeaf_b: coeff for NPP to foliage (evergreen), not used (grass)
    # minimum stem alloc factor (veg=2)
    # ref['AllocLeaf_b']      = {'veg': [0, 1, 2 ,3],
    #                            'min': [0.01, 0.004, 0.01],'max': [0.15, 0.03, 0.15]}
    # AllocStem_a: coeff for NPP to stem (evergreen), not used (grass)
    # water/temperature modulator of allocation (veg=2)
    # ref['AllocStem_a']      = {'veg': [0, 1, 2, 3],
    #                            'min': [0.5, 2., 0.5],'max': [0.9, 3.5, 0.9]}
    # AllocStem_b: coeff for NPP to foliage (evergreen), not used (grass & deciduous)
    # ref['AllocStem_b']      = {'veg': [0, 0, 1, 0], 'min': 0.03,'max': 0.2}
    
    ref['Gs_light'] = {'veg': 1, 'min': [0., 0., 0., 0.],
                       'max': [300., 300., 300., 300.]}
    ref['Gs_vpd'] = {'veg': 1, 'min': [1e-5, 1e-5, 1e-5, 1e-5],
                     'max': [2e-3, 2e-3, 2e-3, 2e-3], 'log': 1}
    ref['LWP_low'] = {'veg': 1, 'min': [200.,200., 200., 200.],
                      'max': [600., 600., 600., 600.]}
    ref['LWP_high'] = {'veg':1, 'min': [50., 50., 50., 50.],
                       'max': [150., 150., 150., 150.]}
    # Wilting point taken from LWP_low
    # ref['WiltPnt'] = {'veg':1, 'min': 120., 'max': 250.} 
    # ref['SLA'] = {'veg': 1, 'min': [0.015, 0.015, 0.002, 0.003],
    #               'max': [0.07, 0.07, 0.02, 0.03]}
    # ref['TurnovL'] = {'veg': [1, 2, 0, 3], 'min': [1.5e-8, 1.5e-8, 1.5e-8],
    #                   'max': [6.5e-8, 6.5e-8, 6.5e-8]}
    # ref['TurnovR'] = {'veg': [1, 2, 2, 2],
    #                   'min': [4e-9, 4e-9], 'max': [1e-8, 1e-8]}
    ref['CWS_max'] = {'veg': 1, 'min': [1e-5, 1e-5, 1e-5, 1e-5],
                      'max': [4e-4, 4e-4, 4e-4, 4e-4]}
    ref['Albedo'] = {'veg': 1, 'min': [0.1, 0.1, 0.1, 0.1],
                     'max': [0.3, 0.3, 0.3, 0.3]}
    ref['KBeers'] = {'veg': 1, 'min': [0.3, 0.3, 0.3, 0.3],
                     'max': [0.8, 0.8, 0.8, 0.8]}
    # ref['WUE_cnp'] = {'veg':1, 'min': 900.,'max': 1500.}
    # ref['Kroot'] = {'veg': 1, 'min': [0.3, 0.2, 0.05, 0.05, 0.3],
    #                 'max': [1, 0.5, 0.3, 0.5, 10.], 'log': 1} 



    
    
