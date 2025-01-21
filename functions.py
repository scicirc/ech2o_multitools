#!/usr/bin/python3
# -*- coding: utf-8 -*-
# *************************************************
#
# Functions for ech2o_multitools.py
# -------
# Author: S. Kuppel
# Created on 04/2020
# -------------------------------------------------

import pcraster as pcr
import numpy as np
import pandas as pd
import scipy.io as spio

import os
import sys
import glob
import copy
import time
import csv
import pyDOE

from datetime import timedelta
from datetime import datetime
from distutils.dir_util import mkpath
from distutils.file_util import copy_file

import GOFs
import spotpy_forked.spotpy as spotpy

# ==================================================================================
# Initialization functions


def config_init(options):

    # Current working directory (temporary, just to import the def file)
    cwd_tmp = os.getcwd()+'/'

    # General definition file for the script
    if options.file is not None:
        [file_py, ext] = os.path.splitext(options.file)
        if len(glob.glob(options.file)) == 0:
            sys.exit('# STOP. The file that define the script configuration' +
                     'does not exist : \n   ...'+file_py+'.py')
    else:
        sys.exit('Please specify a file defining the cscript configuration')

    # -- Import classes and setup from the def file
    sys.path.insert(0, cwd_tmp)
    Config = __import__(file_py).Config
    Opti = __import__(file_py).Opti
    Obs = __import__(file_py).Obs
    Paras = __import__(file_py).Paras
    Site = __import__(file_py).Site

    # Store name of def file in Config
    Config.file = copy.copy(options.file)
    frep = os.path.dirname(Config.file)
    if frep == '':
        Config.file = os.path.join(cwd_tmp, Config.file)

    # -- Main mode
    if Config.mode not in ['calib_MCsampling', 'calib_MCruns',
                           'calib_SPOTPY', 'forward_runs',
                           'sensi_morris']:
        sys.exit("Please choose a valid script mode!")

    # Working directory specified in the def file
    if not hasattr(Config, 'PATH_MAIN'):
        Config.PATH_MAIN = cwd_tmp

    # -- MS init: many things don't happen if it is the case
    Config.SA_init = 0
    if Config.mode == 'sensi_morris':
        if options.SA_init is None:
            sys.exit("Please state if you're initializating the MS sampling")
        else:
            Config.SA_init = int(options.SA_init)
            
        if not hasattr(Opti, 'MSspace'):
            Opti.MSspace = 'trajectory'
        elif Opti.MSspace not in ['trajectory', 'radial']:
            sys.exit('Wrong specification of morris walking mode in the' +
                     "parameter space ('trajectory' or 'radial')")
        # Number of trajectories -> even number determined form ncpu
        if not hasattr(Opti, 'nr'):
            Opti.nr = 10

    # -- Output directory
    # Base name
    if options.outdir is None:
        if Config.mode != 'forward_runs':
            Config.outdir = os.path.splitext(options.file)[0]
        else:
            tmp = options.cfg.split('_')[1].split('.')[0].split('-')
            if(len(tmp) == 2):
                [pref1, pref2] = tmp
                Config.outdir = 'Res.ens'+options.nEns+'_'+pref2 + \
                                '.'+pref1+'.'+options.inEns
            if(len(tmp) == 1):
                Config.outdir = 'Res.ens'+options.nEns+'.'+tmp[0]+'.' +\
                                options.inEns
            if(len(tmp) > 2):
                sys.exit('Error: incorrect config file format.')
    else:
        Config.outdir = copy.copy(options.outdir)

    if Config.mode == 'calib_MCruns' or (Config.mode=='sensi_morris' and Config.SA_init==0):
        if options.task is not None:
            tmp = options.task.split('.')
            print(tmp)
            if len(tmp) == 2:
                Config.indir = tmp[0]
                Config.tasknum = '%03i' % int(tmp[1])
                Config.tasknum2 = tmp[1]
            else:
                tmp = options.task.split('_')
                if len(tmp)==2:
                    Config.indir = tmp[0]
                    Config.tasknum = '%03i' % int(tmp[1])
                    Config.tasknum2 = tmp[1]
                else:
                    sys.exit('Please specify an task ID like XX.tasknumber or +',
                             'XX_tasknumber as a format')
        else:
            sys.exit('Please specify an subtask like XX.number or +',
                     'XX_tasknumber as a format')


    # Absolute location
    if Config.mode == 'forward_runs' and options.OMP_it is not None:
        Config.OMP_it = int(options.OMP_it)
        Config.PATH_OUTmain = \
            os.path.abspath(os.path.join(Config.PATH_MAIN, Config.outdir))
        if len(glob.glob(Config.PATH_OUTmain)) == 0:
            mkpath(Config.PATH_OUTmain)
        Config.PATH_OUT = \
            os.path.abspath(os.path.join(Config.PATH_OUTmain,
                                         'EnsembleRun_' + str(Config.OMP_it)))
    elif Config.mode == 'calib_MCruns':
        Config.PATH_OUTmain = \
            os.path.abspath(os.path.join(Config.PATH_MAIN, Config.outdir))
        if len(glob.glob(Config.PATH_OUTmain)) == 0:
            mkpath(Config.PATH_OUTmain)
        Config.PATH_OUT = \
            os.path.abspath(os.path.join(Config.PATH_OUTmain,
                                         'task' + Config.tasknum))
    elif Config.mode == 'sensi_morris' and Config.SA_init == 0:
        Config.PATH_OUTmain = \
            os.path.abspath(os.path.join(Config.PATH_MAIN, Config.outdir))
        if len(glob.glob(Config.PATH_OUTmain)) == 0:
            mkpath(Config.PATH_OUTmain)
        Config.PATH_OUT = \
            os.path.abspath(os.path.join(Config.PATH_OUTmain,
                                         'task' + Config.tasknum))
    else:
        Config.PATH_OUT = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                       Config.outdir))

    # -- Restart: do not start from first iteration
    if options.restart is None:
        Config.restart = 0
    else:
        Config.restart = int(options.restart)
        if Config.restart > 1:
            sys.exit('Wrong value for restart')

    # -- Run ECH2O?
    Config.runECH2O = 1
    if Config.mode == 'calib_MCsampling' or \
       (Config.mode == 'sensi_morris' and Config.SA_init == 1):
        Config.runECH2O = 0
    # Number of CPUs used in multi-threading: first check if it's given
    # in the options (may be the case for Slurm-type MPI mode), otherwise
    # in Config, default = 1
    if options.ncpu is not None:
        Config.ncpu = copy.copy(options.ncpu)
    elif not hasattr(Config, 'ncpu'):
        Config.ncpu = 1

    # Determine if a parallel computing mode is activated
    Opti.parallel = False
    if not hasattr(Opti, 'SPOTpar'):
        # Even outside SPOTpy mode
        Opti.SPOTpar = 'seq'  # By default, sequential runs
    elif Opti.SPOTpar in ['mpc', 'mpi']:
        Opti.parallel = True
        # Determine with process we're on, to avoid multiple print/file
        # edits later on
        if Opti.SPOTpar == 'mpi':
            if 'OMPI_COMM_WORLD_RANK' in os.environ.keys():
                Opti.rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
            elif 'PMI_RANK' in os.environ.keys():
                Opti.rank = int(os.environ['PMI_RANK'])

    # Database output for SPOTPY
    if Config.mode == 'calib_SPOTPY':
        if hasattr(Opti, 'SPOTdb'):
            # - Custom case, format txt with separate files per obs, one
            # general file and a algo diagnostic file
            # Used by spot_setup
            Opti.dbname = Config.PATH_OUT + '/' + Opti.SPOTalgo + 'ech2o'
            Opti.dbformat = 'custom'
            # Used by sampler.sample
            Opti.dbname2 = None
        else:
            # - Default case, format csv and db call in sampler.sample
            Opti.dbname = None  # Used by spot_setup
            Opti.dbformat = 'csv'
            # Used by sampler.sample
            Opti.dbname2 = Config.PATH_OUT + '/' + Opti.SPOTalgo + 'ech2o'

    # print(options.outdir)
    # -- Calibration: all parameter path (and datasets, if needed)
    if Config.mode in ['calib_MCsampling', 'calib_MCruns']:
        if not hasattr(Config, 'PATH_PAR'):
            print('Warning: path to parameter samples not specified, ' +
                  'set to default')
            Config.PATH_PAR = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                           'Calibration_Samples'))

        if Config.mode == 'calib_MCsampling':
            Config.FILE_PAR = Config.PATH_PAR+'/'+Config.outdir + '_parameters.'

        if Config.mode == 'calib_MCruns':
            Config.FILE_PAR = Config.PATH_PAR+'/'+Config.indir + '_parameters.' + \
                              Config.tasknum2+'.txt'

        # -- Creation of output directory
        if len(glob.glob(Config.PATH_PAR)) == 0:
            mkpath(Config.PATH_PAR)
        # -- Some verbose
        # print('')
        # print("Parameter samples file:\n", Config.FILE_PAR)
        # print('')

    if Config.mode.split('_')[0] == 'calib':
        if not hasattr(Config, 'PATH_OBS'):
            print('Warning: path to calibration datasets not specified, ' +
                  'set to default')
            # Observations (for now only needed in DREAM mode)
            Config.PATH_OBS =\
                os.path.abspath(os.path.join(Config.PATH_MAIN,
                                             'Calibration_Datasets'))

            # -- Output of elementary effects
            # if(Config.SA_init == 0):
            # Config.PATH_EE = os.path.abspath(os.path.join(Config.PATH_MAIN,
            #                                               'ElementaryEffects'))
        # print("Elementary effects directory:          ", Config.PATH_EE)
        # Config.FILE_EE = Config.PATH_OUT+'/'+options.outdir.split('.')[0]

    # print('')

    # -- Morris sensitivity: parameter sets to use
    if Config.mode == 'sensi_morris':
        if not hasattr(Config, 'PATH_PAR'):
            print('Warning: path to parameter samples / trajectories not specified, ' +
                  'set to default')
            Config.PATH_PAR = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                           'Sensitivity_ParameterSets'))

        if Config.SA_init == 1:
            if Opti.MSspace == 'trajectory':
                Config.FILE_PAR = Config.PATH_PAR+'/'+Config.outdir + '_parameters.traj'
            if Opti.MSspace == 'radial':
                Config.FILE_PAR = Config.PATH_PAR+'/'+Config.outdir + '_parameters.rad'
        else :
            if Opti.MSspace == 'trajectory':
                Config.FILE_PAR = Config.PATH_PAR+'/'+Config.indir + '_parameters.traj' + \
                                  Config.tasknum2+'_values.txt'
            if Opti.MSspace == 'radial':
                Config.FILE_PAR = Config.PATH_PAR+'/'+Config.indir + '_parameters.rad' + \
                                  Config.tasknum2+'_values.txt'

        # -- Creation of output directory
        if len(glob.glob(Config.PATH_PAR)) == 0:
            mkpath(Config.PATH_PAR)

    # ---------------------------------------------------------------------------
    # Return the classes read in the definition file
    return (Config, Opti, Obs, Paras, Site)


def param_init(Config, Opti, Paras, Site, options):

    Paras.names = sorted(Paras.ref.keys())
    Paras.n = len(Paras.names)
    # Read dictionary to get all param setup
    Opti.min = []
    Opti.max = []
    Opti.guess = []
    Opti.std = []
    Opti.log = []
    Opti.names = []
    Opti.ind = []
    Paras.comp = {}
    Paras.ind = {}
    ipar = 0
    ipar2 = 0
    Paras.Veg = 0 # report if at least one param is veg-dependent
    Paras.spa1 = 0 # report if at least one param is soil-unit-dependent
    Paras.spa2 = 0 # report if at least one param is geol-unit-dependent
    Paras.spa3 = 0 # report if at least one param is geol-unit-dependent

    # Initial Sampling type (for DREAM)
    if not hasattr(Opti, 'initSample'):
        Opti.initSample = 'uniform'

    # -- Construct parameter-realted vectors used in calibration,
    # -- sensitivity or ensemble runs
    for par in Paras.names:

        # print(par)

        # Default
        if 'map1' not in Paras.ref[par].keys():
            Paras.ref[par]['map1'] = 0
        if 'map2' not in Paras.ref[par].keys():
            Paras.ref[par]['map2'] = 0
        if 'map3' not in Paras.ref[par].keys():
            Paras.ref[par]['map3'] = 0
        if 'veg' not in Paras.ref[par].keys():
            Paras.ref[par]['veg'] = 0

        # Number of components for this parameters : depends on map1, map2 or map3 &
        # veg dependence (veg can be not all species present), default 1
        # Opti.names: unique name using par + component (none, soil/geol/veg type)
        # Opti.comp: which of the component are calibrated? (none, all, some)
        if Paras.ref[par]['map1'] == 0 and Paras.ref[par]['map2'] == 0 and \
           Paras.ref[par]['map3'] == 0 and Paras.ref[par]['veg'] == 0:

            # No mapped nor vegetation dependance
            nr = 1
            Opti.names = Opti.names + [par]
            Paras.comp[par] = 0
            Paras.ref[par]['map'] = 0

        elif Paras.ref[par]['map1'] != 0:

            Paras.ref[par]['map'] = 1
            Paras.spa1 = 1

            # Case where mapped parameters are distributed
            # over all the map1 units listed in the def file
            if Paras.ref[par]['map1'] ==1:
                nr = Site.nmap1
                Opti.names = Opti.names + [par + '_' + m for m in Site.maps1]
                Paras.comp[par] = range(Site.nmap1)

            # Case where mapped parameters either:
            # - don't have component calibrated in some map1 units
            # [..,0,...] for that component and/or
            # - some component are shared across various 1 or more map1 units
            # e.g. [0,1,2,2,0] mean that out of 5 possible map1 components,
            # 1st and 5th are not calibrated, 2nd is one component,
            # 3st and 4st is one 2nd component shared across 3st and 4st map1 units
            elif type(Paras.ref[par]['map1']) == list:
                # The length of the list must match that of map1 units
                if len(Paras.ref[par]['map2']) != Site.nmap1:
                    sys.exit('Invalid spatalized (#1) dependence for parameter '+par)
                # number of components
                nr = sum(i>0 for i in np.unique(Paras.ref[par]['map1']))
                # name
                for i in range(1,nr+1):
                    # param_map1Unit if one component
                    if sum(j==i for j in Paras.ref[par]['map1'])==1:
                        Opti.names = Opti.names + [par + '_' + Site.maps1[j] for
                                                   j in range(Site.nmap1) if
                                                   Paras.ref[par]['map1'][j] == i]
                    # param_map1Unit1+map2Unit2+... if more than one
                    if sum(j==i for j in Paras.ref[par]['map1'])>1:
                        Opti.names = Opti.names + [par + '_' +
                                                   '+'.join(Site.maps1[j] for
                                                            j in range(Site.nmap1) if
                                                            Paras.ref[par]['map1'][j] == i)]
                # component, here, rows indices
                Paras.comp[par] = [i for i in range(Site.nmap1) if
                                   Paras.ref[par]['map1'][i] > 0                                ]
            else:
                sys.exit('Invalid spatialized (#1) dependence for parameter '+par)


        elif Paras.ref[par]['map2'] != 0:

            Paras.ref[par]['map'] = 1
            Paras.spa2 = 1

            # Case where mapped parameters are distributed
            # over all the map2 units listed in the def file
            if Paras.ref[par]['map2'] ==1:
                nr = Site.nmap2
                Opti.names = Opti.names + [par + '_' + m for m in Site.maps2]
                Paras.comp[par] = range(Site.nmap2)

            # Case where mapped parameters either:
            # - don't have component calibrated in some map2 units
            # [..,0,...] for that component and/or
            # - some component are shared across various 1 or more map2 units
            # e.g. [0,1,2,2,0] mean that out of 5 possible map2 components,
            # 1st and 5th are not calibrated, 2nd is one component,
            # 3st and 4st is one 2nd component shared across 3st and 4st map2 units
            elif type(Paras.ref[par]['map2']) == list:
                # The length of the list must match that of map2 units
                if len(Paras.ref[par]['map2']) != Site.nmap2:
                    sys.exit('Invalid spatalized (#2) dependence for parameter '+par)
                # number of components
                nr = sum(i>0 for i in np.unique(Paras.ref[par]['map2']))
                # name
                for i in range(1,nr+1):
                    # param_map2Unit if one component
                    if sum(j==i for j in Paras.ref[par]['map2'])==1:
                        Opti.names = Opti.names + [par + '_' + Site.maps2[j] for
                                                   j in range(Site.nmap2) if
                                                   Paras.ref[par]['map2'][j] == i]
                    # param_map2Unit1+map2Unit2+... if more than one
                    if sum(j==i for j in Paras.ref[par]['map2'])>1:
                        Opti.names = Opti.names + [par + '_' +
                                                   '+'.join(Site.maps2[j] for
                                                            j in range(Site.nmap2) if
                                                            Paras.ref[par]['map2'][j] == i)]
                # component, here, rows indices
                Paras.comp[par] = [i for i in range(Site.nmap2) if
                                   Paras.ref[par]['map2'][i] > 0                                ]
            else:
                sys.exit('Invalid spatialized (#2) dependence for parameter '+par)

        elif Paras.ref[par]['map3'] != 0:

            Paras.ref[par]['map'] = 1
            Paras.spa3 = 1

            # Case where mapped parameters are distributed
            # over all the map3 units listed in the def file
            if Paras.ref[par]['map3'] ==1:
                nr = Site.nmap3
                Opti.names = Opti.names + [par + '_' + m for m in Site.maps3]
                Paras.comp[par] = range(Site.nmap3)

            # Case where mapped parameters either:
            # - don't have component calibrated in some map3 units
            # [..,0,...] for that component and/or
            # - some component are shared across various 1 or more map3 units
            # e.g. [0,1,2,2,0] mean that out of 5 possible map3 components,
            # 1st and 5th are not calibrated, 2nd is one component,
            # 3st and 4st is one 2nd component shared across 3st and 4st map3 units
            elif type(Paras.ref[par]['map3']) == list:
                # The length of the list must match that of map3 units
                if len(Paras.ref[par]['map3']) != Site.nmap3:
                    sys.exit('Invalid spatalized (#3) dependence for parameter '+par)
                # number of components
                nr = sum(i>0 for i in np.unique(Paras.ref[par]['map3']))
                # name
                for i in range(1,nr+1):
                    # param_map3Unit if one component
                    if sum(j==i for j in Paras.ref[par]['map3'])==1:
                        Opti.names = Opti.names + [par + '_' + Site.maps3[j] for
                                                   j in range(Site.nmap3) if
                                                   Paras.ref[par]['map3'][j] == i]
                    # param_map3Unit1+map2Unit2+... if more than one
                    if sum(j==i for j in Paras.ref[par]['map3'])>1:
                        Opti.names = Opti.names + [par + '_' +
                                                   '+'.join(Site.maps3[j] for
                                                            j in range(Site.nmap3) if
                                                            Paras.ref[par]['map3'][j] == i)]
                # component, here, rows indices
                Paras.comp[par] = [i for i in range(Site.nmap3) if
                                   Paras.ref[par]['map3'][i] > 0                                ]
            else:
                sys.exit('Invalid spatialized (#3) dependence for parameter '+par)


        elif Paras.ref[par]['veg'] != 0:

            Paras.ref[par]['map'] = 0
            Paras.Veg = 1

            # Case where there will be as many vegetation components
            # for this parameter as there are vegetation types
            if Paras.ref[par]['veg'] == 1:
                nr = Site.nv
                Opti.names = Opti.names + [par + '_' + s for s in Site.vegs]
                Paras.comp[par] = range(Site.nv)

            # Case where some vegetation type don't have the parameter calibrated,
            # ([..,1 or more,...]>0) while others do [..,0,...])
            # all components > 0 with same integer share the associated parameter
            # e.g. [0,1,2,2,0] mean that out of 5 veg types, 1st and 5th are not
            # calibrated, 2nd has one component, and 3st and 4st share common parameter
            elif type(Paras.ref[par]['veg']) == list:
                # The length of the list must match that of vegetation types
                if len(Paras.ref[par]['veg']) != Site.nv:
                    sys.exit('Invalid veg dependence for parameter '+par)
                # number of components
                nr = sum(i>0 for i in np.unique(Paras.ref[par]['veg']))
                # name
                for i in range(1,nr+1):
                    # param_VegName if one component
                    if sum(j==i for j in Paras.ref[par]['veg'])==1:
                        Opti.names = Opti.names + [par + '_' + Site.vegs[j] for
                                                   j in range(Site.nv) if
                                                   Paras.ref[par]['veg'][j] == i]
                    # param_VegName1+VegName2+... if more than one
                    if sum(j==i for j in Paras.ref[par]['veg'])>1:
                        Opti.names = Opti.names + [par + '_' +
                                                   '+'.join(Site.vegs[j] for
                                                            j in range(Site.nv) if
                                                            Paras.ref[par]['veg'][j] == i)]
                # component, here, rows indices in reference SpeciesParams.tab param file
                Paras.comp[par] = [i for i in range(Site.nv) if
                                   Paras.ref[par]['veg'][i] > 0                                ]
            else:
                sys.exit('Invalid veg dependence for parameter '+par)
        else:
            sys.exit('Invalid spatial/veg dependence for parameter '+par)

        # Link between Paras and Opti indices:
        # 1. Which paras entry is covered in each Opti.* position?
        Opti.ind += list(np.repeat(ipar, nr))
        # Vice versa: which Opti.* indice(s) correspond to a given par?
        if nr>1 or type(Paras.comp[par]) == list :

            # If there are as many component as there are different soil/geol/species calibrated
            if len(Paras.comp[par]) == nr:
                Paras.ind[par] = list(np.arange(ipar2, ipar2+nr, 1))

            # Else if there are one component covering several veg types
            if len(Paras.comp[par]) > nr:

                if type(Paras.ref[par]['map1']) == list:
                    # soil parameter case
                    tmp = [i for i in Paras.ref[par]['map1'] if i > 0]
                elif type(Paras.ref[par]['map2']) == list:
                    # geol parameter case
                    tmp = [i for i in Paras.ref[par]['map2'] if i > 0]
                elif type(Paras.ref[par]['map3']) == list:
                    # geol parameter case
                    tmp = [i for i in Paras.ref[par]['map3'] if i > 0]
                elif type(Paras.ref[par]['veg']) == list:
                    # vegetation parameter case
                    tmp = [i for i in Paras.ref[par]['veg'] if i > 0]
                else:
                    sys.exit('Error when mapping parameter '+par+' to optimization vector')

                Paras.ind[par] = []
                for i in range(nr):
                    Paras.ind[par] += [i+ipar2 for j in tmp if j==np.unique(tmp)[i]]

        else:
            Paras.ind[par] = [ipar2]

        # Build vectors used in the optimisation
        if Config.mode not in ['calib_MCruns', 'forward_runs'] or Config.SA_init == 1 :

            # Log-sampling ?
            if 'log' not in Paras.ref[par].keys():
                Opti.log += list(np.repeat(0, nr))
            else:
                Opti.log += list(np.repeat(Paras.ref[par]['log'], nr))

            # For sampling (uniform -and normal if no sd provided-)
            if 'min' in Paras.ref[par].keys():
                if type(Paras.ref[par]['min']) in [float, int]:
                    Opti.min += [Paras.ref[par]['min']]
                    Opti.max += [Paras.ref[par]['max']]
                elif len(Paras.ref[par]['min']) == nr:
                    Opti.min += Paras.ref[par]['min']
                    Opti.max += Paras.ref[par]['max']
                else:
                    sys.exit('Wrong "min" or "max" format for '+par+
                             ',maybe check soil/geol/veg units in def file?')
            else:
                Opti.min += list(np.repeat(np.nan, nr))
                Opti.max += list(np.repeat(np.nan, nr))
            # In log sampling case with spotpy, log-transform boundaries too
            if Config.mode == 'calib_SPOTPY' and Opti.log[ipar2] == 1:
                Opti.min[ipar2:ipar2+nr+1] = \
                    np.log10(Opti.min[ipar2:ipar2+nr+1]).tolist()
                Opti.max[ipar2:ipar2+nr+1] = \
                    np.log10(Opti.max[ipar2:ipar2+nr+1]).tolist()

            # For normal sampling (or uniform if guess provided)
            if 'guess' in Paras.ref[par].keys():
                if type(Paras.ref[par]['guess']) in [float, int]:
                    Opti.guess += [Paras.ref[par]['guess']]
                elif len(Paras.ref[par]['guess']) == nr:
                    Opti.guess += Paras.ref[par]['guess']
                else:
                    sys.exit('Wrong "guess" format for', par)
            else:
                Opti.guess += list(np.repeat(np.nan, nr))

            # For normal sampling: 0.1 of min/max range if not specified
            if 'std' in Paras.ref[par].keys():
                if type(Paras.ref[par]['std']) in [float, int]:
                    Opti.std += [Paras.ref[par]['std']]
                elif len(Paras.ref[par]['std']) == nr:
                    Opti.std += Paras.ref[par]['std']
                else:
                    sys.exit('Wrong "std" format for', par)
            else:
                # Parameters
                sdmult = 0.1
                # (Log sampling case taken into account above)
                Opti.std += [sdmult*(Opti.max[ipar2+j]-Opti.min[ipar2+j])
                             for j in range(nr)]

        # Increment indices
        ipar += 1
        ipar2 += nr

    # Total number of variables
    Opti.nvar = len(Opti.names)

    # -- Calibration sampling: generate calibration parameters samples
    if Config.mode == 'calib_MCsampling':

        if options.sampling is None:
            Config.sampling = 'uniform'
        else:
            Config.sampling = options.sampling

        # Total number of samples
        Opti.nsamptot = Opti.nit*int(Config.ncpu)
        print('')
        print(str(Opti.nsamptot)+' sets of '+str(Opti.nvar) +
              ' parameters will be generated...')

        # Generate enough parameters sets and save them
        MC_sample(Opti, Config)
        # Done
        print('Parameters sample generation done.')

    # -- MC calibration or ensemble runs: retrieve previously stored
    # -- samples or ensemble
    if Config.mode in ['calib_MCruns', 'forward_runs']:
        # print('Get parameters samples for this job...')
        param_get(Opti, Config, options)

    # -- Sensitivity analysis: generate morris trajectories
    # (SA_init = 1) or read them from files (SA_init = 0)
    if Config.mode == 'sensi_morris':

        if Config.SA_init == 1 :
            
            trajs(Config, Opti)
            print('Parameters trajectory generation done.')

        else:

            # Normalized step: plus-minus 0.5 of the normalized range of
            # each parameter
            Opti.stepN = np.zeros((Opti.nvar), np.float64) + 0.5   
            
            # Total number of runs
            Opti.nruns = Opti.nvar+1
         
            # Read parameter trajectories / sets from files
            print('Get parameters trajectory for this job...')
            param_get(Opti, Config, options)

            
def runs_init(Config, Opti, Obs, Paras, Site, options):

    # -- Directories for runs
    # Forcings and reference maps
    if not hasattr(Config, 'PATH_SPA_REF'):
        print('Warning: path to EcH2O input maps not specified, ' +
              'set to default')
        Config.PATH_SPA_REF = os.path.join(Config.PATH_MAIN, 'Input_Maps')

    if not hasattr(Config, 'PATH_CLIM'):
        print('Warning: path to EcH2O climate input not specified, ' +
              'set to default')
        Config.PATH_CLIM = os.path.join(Config.PATH_MAIN, 'Input_Climate')

    # Creation of inputs directory (PATH_SPA will use parameter sampling)
    # (in case of parallel runs as in DREAM+mpc/mpi, it will serve as
    # -yet another- template)
    # print(os.environ.keys())
    # if Opti.parallel == False:
    Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,
                                                   'Spatial'))
    # else:
    #     # In parallel computing case, one PATH_SPA per task
    #     if Opti.SPOTpar  == 'mpi':
    #         # Running n parallel on a unix system (open MPI type).
    #         # Check the ID of the current mpi task
    #         if 'OMPI_COMM_WORLD_RANK' in os.environ.keys():
    #             call = str(int(os.environ['OMPI_COMM_WORLD_RANK'])+1)
    #         elif 'PMI_RANK' in os.environ.keys():
    #             call = str(int(os.environ['PMI_RANK'])+1)
    #         else:
    #             sys.exit('The ID of this task could not be found...')
    #     if Opti.SPOTpar == 'mpc':
    #         # Running n parallel on a single (Windows) computer.
    #         # ID of the current computer core
    #         call = str(os.getpid())
    #     Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,
    #                                                    'Spatial_'+call))

    # # Copy of reference input data
    # # print(Config.PATH_SPA_REF)
    # # print(Config.PATH_SPA)
    # mkpath(Config.PATH_SPA)
    # # copy_tree(Config.PATH_SPA_REF, Config.PATH_SPA)
    # try:
    #     os.system('cp -f '+Config.PATH_SPA_REF+'/*.map '+
    #               Config.PATH_SPA)
    # except(FileExistsError):
    #     print('')
    # try:
    #     os.system('cp -f '+Config.PATH_SPA_REF+'/SpeciesParams.tab '+
    #               Config.PATH_SPA)
    # except(FileExistsError):
    #     print('')
    # copy_file(Config.PATH_SPA_REF+'/SpeciesParams.tab', Config.PATH_SPA,
    #         update=1)

    # -- About running EcH2O
    # Executable
    if hasattr(Config, 'exe'):
        if len(glob.glob(Config.exe)) == 0:
            sys.exit('The user provided EXEC file was not found: ' +
                     Config.exe)
            print('The user provided EXEC file is: '+Config.exe)
    elif Config.runECH2O == 1:
        sys.exit('The EcH2O EXEC file needs to be specified')

    # Remove if existing (serves as marker for parallel process, see init.files)
    if not (Opti.SPOTpar == 'mpi' and Opti.rank !=0):
        try:
            os.remove(os.path.join(Config.PATH_OUT, Config.exe))
        except:
            print("Could not delete", Config.exe, "(doesn't exist)")
    # try:
    #     os.system('rm -f '+os.path.join(Config.PATH_OUT, Config.exe))
    # except(FileNotFound):
    #     print('No symlink found, no need to remove')
    # # Symbolic link to executable
    # try:
    #     os.symlink(os.path.join(Config.PATH_MAIN, Config.exe),
    #                os.path.join(Config.PATH_OUT, Config.exe))
    # except(FileExistsError):
    #     print('No symlink created, the exec file is alredy there')

    # Time wall for ECH2O execution
    if options.tlimit is None:
        Config.tlimit = '250'
    else:
        Config.tlimit = options.tlimit
    # Time limit
    Config.tcmd = 'ulimit -t ' + \
        str(int(Config.tlimit)*int(Config.ncpu))+' ;'
    # Execution command with Open MP use
    Config.cmde_ech2o = ' '.join([Config.tcmd, 'OMP_NUM_THREADS=' +
                                  str(Config.ncpu), './' + Config.exe])
    # Scratch disk?
    Config.scratch = 0
    if options.scratch is not None:
        if int(options.scratch) == 1:
            Config.scratch = 1
    # Execution directory (in case of parallel runs
    # with DREAM, the acutal directory will be next to this one)
    if Config.scratch == 1:
        print('-----------------------------------------')
        print("Scratch storage activated! Hopefully that will speed " +
              "things up...")
        if hasattr(Config, 'PATH_SCRATCH'):
            if Config.mode == 'calib_MCruns':
                #Config.PATH_EXEC = '/scratch/sylvain.kuppel/MCruns.'+Config.outdir+ \
                Config.PATH_EXEC = Config.PATH_SCRATCH+'/MCruns.'+Config.outdir+ \
                                   '_tmp'+Config.tasknum
            elif Config.mode == 'sensi_morris':
                #Config.PATH_EXEC = '/scratch/sylvain.kuppel/MCruns.'+Config.outdir+ \
                Config.PATH_EXEC = Config.PATH_SCRATCH+'/SAruns.'+Config.outdir+ \
                                   '_tmp'+Config.tasknum
            else:
                #Config.PATH_EXEC = '/scratch/sylvain.kuppel/'+Config.outdir
                Config.PATH_EXEC = Config.PATH_SCRATCH+'/'+Config.outdir
        else:
            sys.exit('Error: no scratch/workdir storage path specified.')
        # if Config.scratch == 2:
        #     Config.PATH_EXEC = '/nobackup/users/s08sk8/'+options.outdir
    else:
        if Config.mode in ['calib_MCruns','sensi_morris']:
            Config.PATH_EXEC = copy.copy(Config.PATH_OUT)
        else:
            Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT, 'tmp'))

    # # In parellel computing case, one execution folder per task
    # if Opti.parallel == True:
    #     # In parallel computing case, one PATH_SPA per task
    #     if Opti.SPOTpar  == 'mpi':
    #         # Running n parallel on a unix system (open MPI type).
    #         # Check the ID of the current mpi task
    #         if 'OMPI_COMM_WORLD_RANK' in os.environ.keys():
    #             call = str(int(os.environ['OMPI_COMM_WORLD_RANK'])+1)
    #         elif 'PMI_RANK' in os.environ.keys():
    #             call = str(int(os.environ['PMI_RANK'])+1)
    #         else:
    #             sys.exit('The ID of this task could not be found...')
    #     elif Opti.SPOTpar == 'mpc':
    #         # Running n parallel on a single (Windows) computer.
    #         # ID of the current computer core
    #         call = str(os.getpid())

    #     Config.PATH_EXEC += '_' + call
    #     print('path exec 1:', Config.PATH_EXEC)

    # -- Ech2O run config file
    # Path for EcH2O config file
    if not hasattr(Config, 'PATH_CFG'):
        print('Warning: path to EcH2O configs not specified, set to default')
        Config.PATH_CFG = os.path.join(Config.PATH_MAIN, 'Input_Configs')

    if options.cfg is not None:
        Config.cfg_ech2o = options.cfg+'.ini'
        # Full path
        Config.FILE_CFG = os.path.join(Config.PATH_CFG, Config.cfg_ech2o)
        if len(glob.glob(Config.FILE_CFG)) == 0:
            sys.exit('The user provided CFG file was not found: ' +
                     Config.FILE_CFG)
        # Where the config file will be (first) copied
        Config.cfg2_ech2o = options.cfg+'.ini'
        Config.FILE_CFGdest = os.path.join(Config.PATH_OUT, Config.cfg2_ech2o)

    else:
        sys.exit('Error: the script need a template ech2o config file!')

    # (if needed) get the parallel job number, based on the output dir name
    # Config.tasknum = options.outdir.split('.')[::-1][0]

    # -- Tracking age and/or tracers?
    if not hasattr(Site, 'isTrck'):
        Site.isTrck = 0
    if Site.isTrck == 1:
        if options.cfgTrck is not None:
            Config.cfgTrck_ech2o = options.cfgTrck+'.ini'
        else:
            Config.cfgTrck_ech2o = options.cfg.split('_')[0] + \
                'Trck_'+options.cfg.split('_')[1]+'.ini'
        # Full path
        Config.FILE_CFGtrck = os.path.join(Config.PATH_CFG,
                                           Config.cfgTrck_ech2o)
        if len(glob.glob(Config.FILE_CFG)) == 0:
            sys.exit('The user provided CFGtrck file was not found: ' +
                     Config.FILE_CFGtrck)

    # -- Forward runs: parameter sets to use
    if Config.mode == 'forward_runs':
        if options.inEns is not None:
            if not hasattr(Config, 'PATH_PAR'):
                print('Warning: path to ensemble parameters not ' +
                      'specified, set to default')
                Config.PATH_PAR = os.path.join(Config.PATH_MAIN,
                                               'Input_Params')
            Config.FILE_PAR = Config.PATH_PAR + '/' + options.inEns+'.txt'
            # Config.FILE_PAR = Config.PATH_MAIN+'Input_Params/'+\
            # options.inEns+'.'+
            # options.nEns+'bestParams.txt'
            if len(glob.glob(Config.FILE_PAR)) == 0:
                sys.exit('The param file (ensemble set) was not found: ' +
                         Config.FILE_PAR)
        else:
            sys.exit('The param file (ensemble set) needs to be ' +
                     'specified (--inEns)')
        # Number of runs is the size of the parameters file, unless specified
        nv = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                           unpack=True)[1::].shape[0]
        if options.nEns is not None:
            Config.nEns = min(int(options.nEns), nv)
        else:
            Config.nEns = nv

        # -- Some verbose
        # print('')
        # print("Parameter samples file:\n", Config.FILE_PAR)
        # print('')
            
    # Maps averaging
    if options.MapAv is not None:
        Config.MapAv = int(options.MapAv)
        if Config.MapAv == 1:
            if options.MapAvT in ['week', 'month', 'season']:
                Config.MapAvT = options.MapAvT
            else:
                sys.exit('Wrong maps averaging option!')
    else:
        Config.MapAv = 0

# ==================================================================================
# -- Read measurements for calibration


def obs_init(Config, Opti, Obs):

    # Set the observations types collected from runs (sim outputs)
    # (and compared to measurements if there is calibration)
    Obs.names = sorted(Obs.obs.keys())

    # Just a flag to make sure when the first map-time files is actually written
    Obs.firstMapTs = {}
    for oname in Obs.names:
        Obs.firstMapTs[oname] = 1

    # Time resolution in simulations, in seconds
    if not hasattr(Obs, 'simres'):
        Obs.simres = 86400  # Daily, by default

    # --- Reporting stuff
    # Trim: only saves the time steps within the trim
    if not hasattr(Obs, 'saveB'):
        Obs.saveB = 1  # Initial time step to report time series
    if not hasattr(Obs, 'saveBmap'):
        Obs.saveBmap = 1  # Initial time step to report maps
    # Length of saved outputs
    if not hasattr(Obs, 'saveL'):
        Obs.saveL = Obs.lsim - Obs.saveB + 1
    else:
        if Obs.saveL > Obs.lsim - Obs.saveB + 1:
            sys.exit('Error: the specified time series slicing start+length ' +
                     'goes beyond simulation time!')
    if not hasattr(Obs, 'saveLmap'):
        Obs.saveLmap = Obs.lsim - Obs.saveBmap + 1
    else:
        if Obs.saveLmap > Obs.lsim - Obs.saveBmap + 1:
            sys.exit('Error: the specified map slicing start+length ' +
                     'goes beyond simulation time!')
    # Report ech2o.log files
    if not hasattr(Config, 'replog'):
        Config.replog = 0
    # Report BasinSummary.txt
    if not hasattr(Obs, 'repBS'):
        Obs.repBS = 0
    # Report BasinSummary.txt for lines at saveB and last line
    if not hasattr(Opti, 'repBS_interval'):
        Opti.repBS_interval = 0

    # Date vector of the simulations timesteps
    Obs.simt = [Obs.simbeg + timedelta(seconds=Obs.simres*x) for x in range(Obs.saveL)]
    Obs.simtmap = [Obs.simbeg + timedelta(seconds=Obs.simres*x) for x in range(Obs.saveLmap)]

    # Default values obs & sim converter = 1
    for oname in Obs.names:
        if not 'obs_conv' in Obs.obs[oname]:
            Obs.obs[oname]['obs_conv'] = 1
        if not 'sim_conv' in Obs.obs[oname]:
            Obs.obs[oname]['sim_conv'] = 1
            
    if Config.mode in ['calib_MCruns','calib_SPOTPY']:
        # print('Reading measured datasets for calibration...')

        Opti.obs = {}  # np.full((Obs.nobs, Obs.saveL), np.nan)
        Opti.obs2 = {}  # in case there's a second model-data fit window given
        Opti.calib2 = {}  # in case there's a second model-data fit window given
        Opti.iscalib2 = False
        Obs.nobs2 = 0

        # Use derivate-based GOFs ?
        Opti.gof_d1 = False
        if Config.mode == 'calib_MCruns':
            for gof in Opti.GOFs:
                if 'd1' in gof:
                    Opti.gof_d1 = True

        for oname in Obs.names:

            # print(oname)
            
            # -- Get the obs
            f_obs = Config.PATH_OBS + '/' + Obs.obs[oname]['obs_file']

            # Time series (point or domain-scale) ---------------------------
            if Obs.obs[oname]['type'] == 'Ts' or Obs.obs[oname]['type'] == 'Total':

                # Read the file, keeping only the data and relevant TS columns
                tmp = pd.read_csv(f_obs, sep=';').iloc[
                    :, [0, Obs.obs[oname]['obs_col']-1]]
                tmp.columns = ['Date', 'value']
                # Convert date column to datetime
                tmp['Date'] = pd.to_datetime(tmp['Date'])
                # Convert date column to datetime
                tmp['value'] = tmp['value'] * Obs.obs[oname]['obs_conv']

                # -- Calibration period:
                # Check if specified, otherwise use the whole simulation
                # in any case, remove the spinup (it will be removed from
                # simulation outputs in post-processing)
                # (no need of Obs.saveB because simt starts at saveB)
                if 'fit_beg' not in Obs.obs[oname].keys() or \
                   type(Obs.obs[oname]['fit_beg']) is not datetime.date:
                    fitbeg = Obs.simt[0]
                else:
                    fitbeg = max(Obs.obs[oname]['fit_beg'], Obs.simt[0])
                if 'fit_end' not in Obs.obs[oname].keys() or \
                   type(Obs.obs[oname]['fit_end']) is not datetime.date:
                    fitend = Obs.simt[Obs.saveL-1]
                else:
                    fitend = min(Obs.obs[oname]['fit_end'],
                                 Obs.simt[Obs.saveL-1])
                # Crop obs between desired time frame
                tmp2 = tmp.loc[(tmp.Date >= fitbeg) & (tmp.Date <= fitend)]
                Opti.obs[oname] = tmp2.dropna(how='any')

                # -- Second calibration period?
                if 'fit_beg2' in Obs.obs[oname].keys() and \
                   'fit_end2' in Obs.obs[oname].keys() :
                    #print(type(Obs.obs[oname]['fit_beg2']))
                    #print(type(Obs.obs[oname]['fit_end2']))
                    if type(Obs.obs[oname]['fit_beg2']) is datetime and \
                       type(Obs.obs[oname]['fit_end2']) is datetime:
                        fitbeg2 = max(Obs.obs[oname]['fit_beg2'], Obs.simt[0])
                        fitend2 = min(Obs.obs[oname]['fit_end2'], Obs.simt[Obs.saveL-1])
                        # Crop obs between desired time frame
                        tmp2 = tmp.loc[(tmp.Date >= fitbeg2) & (tmp.Date <= fitend2)]
                        Opti.obs2[oname] = tmp2.dropna(how='any')
                        Opti.calib2[oname] = True
                        Opti.iscalib2 = True
                        Obs.nobs2 += 1
                    else:
                        Opti.calib2[oname] = False
                else:
                    Opti.calib2[oname] = False

            # Snapshot maps (reading one time step) --------------------------
            elif Obs.obs[oname]['type'] == 'mapStep':

                # Read map, convert to array using a missing value (-9999)
                try:
                    tmp = pcr.pcr2numpy(pcr.readmap(f_obs),
                                        np.nan)*Obs.obs[oname]['sim_conv']
                    Opti.obs[oname] = tmp[~np.isnan(tmp)]
                    #print(tmp)
                except RuntimeError:
                    print('Warning: RuntimeError - could not read ',f_obs)

            # # -- Get derivative if needed for GOFs
            # if Opti.gof_d1 is True:
            #     # between available data points
            #     # with can be longer than simulation resolution
            #     tmp3 = Opti.obs[oname].diff()[1::] # interval & delta value
            #     # tmp3.rename(columns={'Date':'Dt', 'value':'DY'}, inplace=True)
            #     tmp3.columns = ['Dt', 'DY']
            #     # starting date for Dt
            #     tmp3['tbeg'] = Opti.obs[oname]['Date'].loc[tmp3.index] - tmp3['Dt']
            #     # actual rate in (ObsUnit) per second
            #     tmp3['value'] = tmp3['DY'] / tmp3['Dt'].dt.total_seconds()
            #     Opti.obs[oname+'_d1'] = copy.copy(tmp3[['tbeg','value','Dt']])
            #
            #     # Same for second calibration period, if needed
            #     if Opti.calib2[oname] is True:
            #         tmp3 = Opti.obs2[oname].diff()[1::] # interval & delta value
            #         tmp3.columns = ['Dt', 'DY']
            #         # starting date for Dt
            #         tmp3['tbeg'] = Opti.obs2[oname]['Date'].loc[tmp3.index] - tmp3['Dt']
            #         # actual rate in (ObsUnit) per second
            #         tmp3['value'] = tmp3['DY'] / tmp3['Dt'].dt.total_seconds()
            #         Opti.obs2[oname+'_d1'] = copy.copy(tmp3[['tbeg','value','Dt']])



    if Config.mode == 'calib_MCruns':
        # -- Group the output files in one across simulations,
        #    separating by observations points and veg type where it applies
        GOFref = ['NSE','KGE','KGE2012','RMSE','MAE', # classic
                'KGEc','KGE2012c','RMSEc','MAEc', # mean-centered
                'RMSE_d1','MAE_d1', # using derivatives
                  'gauL','gauLerr','logL','corr','rstd','rmu','bias','dstd','sprmn']
        # Check that use-defined GOFs are in the reference list
        tmp = []
        for gof in Opti.GOFs:
            if gof in GOFref:
                tmp += [gof]
        if len(tmp) == 0:
            sys.exit('ERROR: None of that user-defined GOFs are in the reference list!')
        else:
            Opti.GOFs = tmp
            Opti.nGOF = len(Opti.GOFs)
            Opti.GOFfiles = {}
            if Opti.iscalib2 is True:
                Opti.GOFfiles2 = {}

            for gof in Opti.GOFs:
                # Historic time series file names
                Opti.GOFfiles[gof] = Config.PATH_OUTmain+'/'+gof+'.task'+Config.tasknum+'.tab'
                # Header of files
                if Config.restart == 0:
                    #print('first period:',fitbeg,'to',fitend)
                    with open(Opti.GOFfiles[gof], 'w') as f_out:
                        f_out.write('Sample,'+','.join(Obs.names)+'\n')
                # Second calibration period
                if Opti.iscalib2 is True:
                    # Historic time series file names
                    Opti.GOFfiles2[gof] = Config.PATH_OUTmain+'/'+gof+'_p2.task'+ \
                                          Config.tasknum+'.tab'
                    # Header of files
                    if Config.restart == 0:
                        #print('second period:',fitbeg2,'to',fitend2)
                        with open(Opti.GOFfiles2[gof], 'w') as f_out:
                            f_out.write('Sample,'+
                                        ','.join([Obs.names[i]
                                                  for i in range(Obs.nobs)
                                                  if Obs.obs[Obs.names[i]]['type'] != 'mapStep'])+'\n')


        # -- Optional: group the BasinSummary files in one across simulations
        if Opti.repBS_interval == 1:
            # Interval start (saveB)
            Opti.BSfile_tb = Config.PATH_OUTmain+'/BasinBudget_tb.task'+Config.tasknum+'.tab'
            # Interval start (saveB+lsim)
            Opti.BSfile_te = Config.PATH_OUTmain+'/BasinBudget_te.task'+Config.tasknum+'.tab'


# ==========================================================================


def files_init(Config, Opti, Paras, Site):

    # Verbose, copying and editing directories and files.
    # Only do it once, in mpi mode the Opti.rank variable variable is used
    if not (Opti.SPOTpar == 'mpi' and Opti.rank !=0):

        # ==== Introductory verbose
        # if Opti.parallel == False or \
        # (Opti.DREAMpar != 'mpi' or Config.MPIrank == 0):
        print('')
        print('**************************************************************')
        print('The user provided definition file is:\n'+Config.file)
        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
        if Config.mode.split('_')[0] == 'calib':
            print('CALIBRATION with EcH2O: \n')
        elif Config.mode == 'forward_runs':
            print('ENSEMBLE RUNS with EcH2O')
        elif Config.mode == 'sensi_morris':
            print('MORRIS SENSITIVITY with EcH2O: ')
            if(Config.SA_init == 1):
                print('- construction and storage of the parameter trajectories')
            else:
                print('- forward runs')
                print('- storage of outputs and info for posterior analysis : ' +
                      'elementary effects, etc.')
            print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
            print('')
        print('')
        if Config.runECH2O == 1:
            print('Original/template maps & parameters:\n', Config.PATH_SPA_REF)
            print('Climate data:\n', Config.PATH_CLIM)
            if Config.mode.split('_')[0] == 'calib':
                print('Calibration datasets:\n', Config.PATH_OBS)
            if Config.mode == 'calib_MCruns':
                print('Calibration sample file: \n'+Config.FILE_PAR)
            if Opti.parallel is False:
                print('Maps & parameters:\n', Config.PATH_SPA)
                print('Raw output from simulations:\n', Config.PATH_EXEC)
            print('The user provided CFG file is: \n'+Config.FILE_CFG)
            if Config.mode == 'calib_MCruns':
                print('-> for this subtask moved/edited to: \n'+Config.FILE_CFGdest)
            if Site.isTrck == 1:
                print('The user provided CFGtrck file is:\n' + Config.FILE_CFGtrck)
            if Config.mode == 'calib_MCruns' or \
               (Config.mode == 'sensi_morris' and Config.SA_init == 0) or \
               (Config.mode == 'forward_runs' and hasattr(Config, 'OMP_it')):
                print('Final outputs:\n', Config.PATH_OUTmain)
            else:
                print('Final outputs:\n', Config.PATH_OUT)
        elif Config.mode == 'calib_MCsampling' or Config.SA_init == 1:
            print('Calibration samples in :\n', Config.FILE_PAR+'* files')

        print('-----------------------------------------')
        # -- How many variables ?
        print('Total number of parameters used :', Opti.nvar)
        print('')

        if Config.runECH2O == 1:
            # Main output directory: create if needed
            mkpath(Config.PATH_OUT)
            # copy definition file there
            copy_file(Config.file, Config.PATH_OUT)
        # If ensemble runs, copy parameter file
        if(Config.mode == 'forward_runs'):
            copy_file(Config.FILE_PAR, Config.PATH_OUT)
        # EcH2O config file: copy from template and edit
        if Config.runECH2O == 1:
            copy_file(Config.FILE_CFG, Config.FILE_CFGdest)
                
        if Config.runECH2O == 1:
            with open(Config.FILE_CFGdest, 'a') as fw:
                fw.write('\n\n\n#Simulation-specific folder section\n#\n\n')
                fw.write('Clim_Maps_Folder = '+Config.PATH_CLIM+'\n')
                # Further edit regarding tracking
                if Site.isTrck == 1:
                    fw.write('Tracking = 1\n')
                    fw.write('TrackingConfig = configTrck.ini\n')
                else:
                    fw.write('Tracking = 0\n')
                # Input maps directory defined here (unless you need parallel runs
                # with DREAM, in which case the directory will change on the fly,
                # see spot_setup.simulation)
                if Opti.parallel is False:
                    fw.write('Maps_Folder = '+Config.PATH_SPA+'\n')
                    fw.write('Output_Folder = '+Config.PATH_EXEC+'\n')

                # Vegetation parameter files to be used
                fw.write('Species_Parameters = '+Site.vfile+'\n')

            # If tracking, copy template configTrck file for EcH2O
            if Site.isTrck == 1:
                copy_file(Config.FILE_CFGtrck,
                          os.path.join(Config.PATH_OUT, 'configTrck.ini'))
                # Where the configTrck file will be (first) copied
                Config.FILE_CFGdest_trck = os.path.join(Config.PATH_OUT,
                                                        'configTrck.ini')


            # Copy of reference input parameters
            # remove_tree(Config.PATH_SPA)
            mkpath(Config.PATH_SPA)
            os.system('cp -f '+Config.PATH_SPA_REF+'/*.map ' + Config.PATH_SPA)
            # Remove the default map files of calibrated param in the inputs dir.
            # --> helps checking early on if there is an improper map update
            # for pname in Paras.names:
            #     if Paras.ref[pname]['veg'] == 0:
            #     os.remove(Config.PATH_SPA+'/' + Paras.ref[pname]['file']+'.map')
            # os.system('cp -p '+Config.PATH_SPA_REF+'/SpeciesParams*.tab '+
            #           Config.PATH_SPA)
            # Copy reference tables (if Species_State_Variable_Input_Method = tables in ech2o cobfig file)
            os.system('cp -f '+Config.PATH_SPA_REF+'/Spec*.tab ' + Config.PATH_SPA)
            # Copy species parameters reference file (may be redundant if Site.vfile starts with Spec...)
            #copy_file(Config.PATH_SPA_REF+'/SpeciesParams.tab', Config.PATH_SPA)
            copy_file(Config.PATH_SPA_REF+'/'+Site.vfile, Config.PATH_SPA)

            # Keep it as last action in this if !
            # EcH2O executable file: clean up / update old symlink
            print('copy executable from'+os.path.join(Config.PATH_MAIN, Config.exe),'to', Config.PATH_OUT)
            copy_file(os.path.join(Config.PATH_MAIN, Config.exe),
                      Config.PATH_OUT)

    else:
        # In MPI mode, make other processes wait until EcH2O has been copied
        # which means the rank=0 process finished all the preps
        while len(glob.glob(os.path.join(Config.PATH_OUT, Config.exe))) == 0:
            time.sleep(1)

    # === Preparing inputs maps/files for site geometry etc.
    # === (for all processes)

    if Config.runECH2O == 1:
        # -- Soils / geologies / units maps
        # Initialization by cloning base map
        Config.cloneMap = pcr.boolean(pcr.readmap(Config.PATH_SPA+'/base.map'))
        pcr.setclone(Config.PATH_SPA+'/base.map')
        Site.bmaps = {}
        if(Paras.spa1 == 1):
            for im in range(Site.nmap1):
                Site.bmaps[Site.maps1[im]] = pcr.readmap(Config.PATH_SPA+'/' +
                                                         Site.map1files[im])
        if(Paras.spa2 == 1):
            for im in range(Site.nmap2):
                Site.bmaps[Site.maps2[im]] = pcr.readmap(Config.PATH_SPA+'/' +
                                                         Site.map2files[im])
        if(Paras.spa3 == 1):
            for im in range(Site.nmap3):
                Site.bmaps[Site.maps3[im]] = pcr.readmap(Config.PATH_SPA+'/' +
                                                        Site.map3files[im])

        Site.bmaps['unit'] = pcr.readmap(Config.PATH_SPA+'/unit.map')
        # Stream network
        Site.bmaps['chanmask'] = pcr.readmap(Config.PATH_SPA+'/chanmask.map')
        # Site.bmaps['chanmask_NaN'] = pcr.readmap(Config.PATH_SPA +
        #                                          '/chanmask_NaN.map')
        # Bare rock patches
        if hasattr(Site, 'simRock'):
            if Site.simRock == 1:
                # Site.bmaps['nolowK'] = readmap(Config.PATH_SPA+'/unit.nolowK.map')
                Site.bmaps['rock'] = pcr.readmap(Config.PATH_SPA+'/unit.rock.map')

        # -- For initial soil moisture maps generation at each run of EcH2O,
        # get the value of keyword in config file (porosity profile type and
        # name of maps)
        if not hasattr(Site, 'frac_SWC1'):
            Site.frac_SWC1 = 0.9
        if not hasattr(Site, 'frac_SWC2'):
            Site.frac_SWC2 = 0.9
        if not hasattr(Site, 'frac_SWC3'):
            Site.frac_SWC3 = 0.9
        # print('initial conditions 1')
        # print(Config.FILE_CFGdest)
        with open(Config.FILE_CFGdest, 'r') as f:
            datafile = f.readlines()

        # Search and destroy, er, get the values
        for line in datafile:
            # Check which porosity profile mode is on in config file
            if 'Porosity_profile =' in line:
                Site.poros_mode = int(line.split('=')[1].strip())
            # Files names needed in any case
            if 'Porosity =' in line:
                Site.f_poros = line.split('=')[1].strip()
            if 'Soil_moisture_1' in line:
                Site.f_initSWC1 = line.split('=')[1].strip()
            if 'Soil_moisture_2' in line:
                Site.f_initSWC2 = line.split('=')[1].strip()
            if 'Soil_moisture_3' in line:
                Site.f_initSWC3 = line.split('=')[1].strip()
        if not any(s in Site.__dict__.keys() for s in
                   ['poros_mode', 'f_poros',
                    'f_initSWC1', 'f_initSWC2', 'f_initSWC3']):
            sys.exit('Error: file names for poros and init SWC not found')

        if Site.poros_mode == 1:
            # Exponential: profile coeff and depths file names are needed
            for line in datafile:
                if 'Porosity_Profile_coeff =' in line:
                    Site.f_kporos = line.split('=')[1].strip()
                if 'Depth_soil_layer_1 =' in line:
                    Site.f_dL1 = line.split('=')[1].strip()
                if 'Depth_soil_layer_2 =' in line:
                    Site.f_dL2 = line.split('=')[1].strip()
                if 'Soil_depth =' in line:
                    Site.f_dTot = line.split('=')[1].strip()
            if not any(s in Site.__dict__.keys() for s in
                       ['f_kporos', 'f_dL1', 'f_dL2', 'f_dTot']):
                sys.exit('Error: file names for poros profile and depths not found')

        elif Site.poros_mode == 2:
            # Porosity map given for each layer
            for line in datafile:
                if 'Porosity_Layer2 =' in line:
                    Site.f_porosL2 = line.split('=')[1].strip()
                if 'Porosity_Layer3 =' in line:
                    Site.f_porosL3 = line.split('=')[1].strip()
            if not any(s in Site.__dict__.keys() for s in
                       ['f_porosL2', 'f_porosL3']):
                sys.exit('Error: files names for L2 & L3 poros not found')

        # -- Vegetation inputs file: reference dictionary
        Opti.vref = {}
        # Read template file
        with open(Config.PATH_SPA_REF+'/' + Site.vfile, 'r') as csvfile:
            paramread = list(csv.reader(csvfile, delimiter='\t'))
        # "Head": number of species and of params
        Opti.vref['header'] = paramread[0][0:len(paramread[0])]
        # Check that the number of species per params in def files matches
        if Site.nv != len(paramread)-3:
            sys.exit('ERROR: the number of species in def files ('+ str(Site.nv)+
                     ') differs from that in template params file ('+
                     str(len(paramread)-3)+ ')')
        # All parameters values (keep strings!)
        for iv in range(Site.nv):
            Opti.vref[iv] = paramread[iv+1][0:len(paramread[iv+1])]
        # "Footers" : name of head1 and of parameters
        Opti.vref['footer'] = paramread[Site.nv+1][0:len(paramread[Site.nv+1])]
        Opti.vref['name'] = paramread[Site.nv+2][0:len(paramread[Site.nv+2])]

        # -- For initial coyping (or not) initial tracking maps,
        # get the value of keyword in configTrck file
        # if Site.isTrck == 1:
        #     with open(Config.FILE_CFGdest_trck, 'r') as f:
        #         datafile = f.readlines()
        #     # Search and destroy, er, get the trck values
        #     for line in datafile:
        #         # 2H tracking ?
        #         if 'sw_2H =' in line:
        #             Site.sw_2H = int(line.split('=')[1].strip())
        #         # 2H tracking ?
        #         if 'sw_18O =' in line:
        #             Site.sw_18O = int(line.split('=')[1].strip())
        #         # 2H tracking ?
        #         if 'sw_Cl =' in line:
        #             Site.sw_Cl = int(line.split('=')[1].strip())
        #         # 2H tracking ?
        #         if 'sw_Age =' in line:
        #             Site.sw_Age = int(line.split('=')[1].strip())

        # if Config.isTrck == 1:
        #     os.system('cp '+Config.PATH_SPA+'/d2H_snowpack.map '+Config.PATH_SPA+
        # '/d2H.snowpack.map')
        #     os.system('cp '+Config.PATH_SPA+'/d2H_surface.map '+Config.PATH_SPA+
        # '/d2H.surface.map')
        #     os.system('cp '+Config.PATH_SPA+'/d2H_soil1.map '+Config.PATH_SPA+
        # '/d2H.L1.map')
        #     os.system('cp '+Config.PATH_SPA+'/d2H_soil2.map '+Config.PATH_SPA+
        # '/d2H.L2.map')
        #     os.system('cp '+Config.PATH_SPA+'/d2H_soil3.map '+Config.PATH_SPA+
        # '/d2H.L3.map')
        #     os.system('cp '+Config.PATH_SPA+'/d2H_groundwater.map '+
        # Config.PATH_SPA+'/d2H.GW.map')

        #     os.system('cp '+Config.PATH_SPA+'/d18O_snowpack.map '+
        # Config.PATH_SPA+'/d18O.snowpack.map')
        #     os.system('cp '+Config.PATH_SPA+'/d18O_surface.map '+
        # Config.PATH_SPA+'/d18O.surface.map')
        #     os.system('cp '+Config.PATH_SPA+'/d18O_soil1.map '+Config.PATH_SPA+
        # '/d18O.L1.map')
        #     os.system('cp '+Config.PATH_SPA+'/d18O_soil2.map '+Config.PATH_SPA+
        # '/d18O.L2.map')
        #     os.system('cp '+Config.PATH_SPA+'/d18O_soil3.map '+Config.PATH_SPA+
        # '/d18O.L3.map')
        #     os.system('cp '+Config.PATH_SPA+'/d18O_groundwater.map '+
        # Config.PATH_SPA+'/d18O.GW.map')

        #     os.system('cp '+Config.PATH_SPA+'/Age_snowpack.map '+
        # Config.PATH_SPA+'/Age.snowpack.map')
        #     os.system('cp '+Config.PATH_SPA+'/Age_surface.map '+
        # Config.PATH_SPA+'/Age.surface.map')
        #     os.system('cp '+Config.PATH_SPA+'/Age_soil1.map '+Config.PATH_SPA+
        # '/Age.L1.map')
        #     os.system('cp '+Config.PATH_SPA+'/Age_soil2.map '+Config.PATH_SPA+
        # '/Age.L2.map')
        #     os.system('cp '+Config.PATH_SPA+'/Age_soil3.map '+Config.PATH_SPA+
        # '/Age.L3.map')
        #     os.system('cp '+Config.PATH_SPA+'/Age_groundwater.map '+
        # Config.PATH_SPA+'/Age.GW.map')

# ============================================================================
# Routine: Runs functions
# -------------------------------------------------


def calibMC_runs(Config, Opti, Obs, Paras, Site):
    # Calibration runs
    # -------------------

    print('Number of iterations      :', Opti.nit)
    if Config.restart == 1:
        restart(Config, Opti, Obs)
        #print('...but directly restarting from iter. ', Config.itres)
    print('')
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print(' Entering the optimisation loop...')
    print('')

    # Output initialization
    Config.initpar = 0
    Config.initobs = 0
    Opti.begfail = 0

    # Initial clean up, if not restarting
    f_failpar = Config.PATH_OUTmain+'/Parameters_fail.task'+Config.tasknum+'.txt'
    if len(glob.glob(f_failpar)) != 0 and Config.restart == 0:
        os.system('rm -f '+f_failpar)


    if Config.restart == 1:
        it0 = Config.itres-1
    else:
        it0 = 0

    for it in range(it0, Opti.nit):

        Opti.itout = '%04i' % int(it+1)
        print('Iteration ', Opti.itout, ' of ', Opti.nit)

        # Create (if needed) the run outputs directory
        mkpath(Config.PATH_EXEC)
        # Create the inputs for ECH2O
        sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...', end='\r')
        start = time.time()
        os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time()-start,
              'seconds (limit at '+Config.tlimit+')')

        if runOK(Obs, Opti, Config, mode='verbose') == 0:
            # If the run fails, let's give it one more chance!
            time.sleep(1)
            # os.system('rm -f '+Config.PATH_EXEC+'/*')
            print('|| running ECH2O (2nd try)...', end='\r')
            start = time.time()
            os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                      Config.PATH_EXEC+'/ech2o.log')
            print('run time:', time.time()-start,
                  'seconds (limit at '+Config.tlimit+')')
            os.chdir(Config.PATH_EXEC)
            # Still not running properly? Report and move on
            if runOK(Obs, Opti, Config, mode='verbose') == 0:
                f_failpar = Config.PATH_OUTmain+'/Parameters_fail.task'+Config.tasknum+'.txt'
                if len(glob.glob(f_failpar)) == 0:
                    with open(f_failpar, 'w') as f_in:
                        f_in.write('Sample,'+','.join(Opti.names)+'\n')
                with open(f_failpar, 'a') as f_in:
                    f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')
                # If it is the very first iteration, record it for later use
                if it == 0:
                    Opti.begfail = 1

                os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
                          Config.PATH_OUT + '/ech2o.run'+Opti.itout+'.log')



        # Store goodness of fit across outputs (or NA/blank if the run crashed)
        os.chdir(Config.PATH_EXEC)
        # Write goodness of fit (even if it's nan)
        store_GOF(Obs, Opti, Config, Site, it)
        # Write parameters values for this sequence
        # params.store(Opti, Config, it)

        # Store BasinSummary.txt at intervals if asked (and if the simulation worked)
        if Opti.repBS_interval == 1:
            store_BS_interval(Obs, Opti, Config, it)

        # os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
        #           Config.PATH_EXEC + '/ech2o.run'+Opti.itout+'.log')

        os.chdir(Config.PATH_OUT)
        # sys.exit()
        # Clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*.tab')
        os.system('rm -f '+Config.PATH_EXEC+'/Basin*.txt')

    # Final cleanup
    os.system('rm -f '+Config.PATH_EXEC+'/*')
    os.system('rm -fr '+Config.PATH_EXEC)


def forward_runs(Config, Opti, Obs, Paras, Site, options):

    # Non-calibration runs
    # -------------------

    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print('')

    if options.OMP_it is None:
        # Mininum between params file size and prescribed runs
        nruns = Config.nEns
        # print(nruns)
    else:
        nruns = 1

    for it in range(nruns):

        print('Iteration '+str(it+1)+' of '+str(nruns))

        # Create / clean up the run outputs directory
        if len(glob.glob(Config.PATH_EXEC)) == 0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        time.sleep(1)

        # Create the inputs for ECH2O
        sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...', end='\r')
        start = time.time()
        os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time()-start,
              'seconds (limit at', Config.tlimit, ')')
        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        # Group outputs
        store_sim(Obs, Opti, Config, Site, it)

        # if runOK(Obs, Opti, Config, mode='silent') == 0:
          # print('Warning: something came up during the run;' +
            #      ' some outputs might be truncated or missing')
            # Group outputs
            # outputs.store_sim(Obs, Opti, Config, Site, it)

        # Store log, if requested
        if Config.replog == 1:
            os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
                      Config.PATH_OUT + '/ech2o_run'+str(it+1)+'.log')

        # Intermediate clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*.tab')
        os.system('rm -f '+Config.PATH_EXEC+'/Basin*.txt')

    # Clean up
    os.system('rm -f '+Config.PATH_EXEC+'/*')
    os.system('rm -fr '+Config.PATH_EXEC)


def morris_runs(Config, Opti, Obs, Paras, Site):

    # Simulations when varying the parameters, Morris's one-at-a-time
    # ---------------------------------------------------------------
    print('Total Number of iterations:', Opti.nruns)
    print('')
    print('-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print('')

    Config.initobs = 0
    Opti.begfail = 0

    # irun_tot = 0

    # Initial clean up
    f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
    if len(glob.glob(f_failpar)) != 0:
        os.system('rm -f '+f_failpar)

    print('======================================')
    print('## Runs along trajectory #', str(Config.tasknum2))
    print('--------------------------------------')

    # There are npara+1 runs for each trajectory
    for irun in range(Opti.nvar+1):

        # runnb = '%02i' % int(irun+1)
        print('Run '+str(irun+1)+' out of '+str('%02i' % int(Opti.nvar+1)))

        # Create / clean up the run outputs directory
        if len(glob.glob(Config.PATH_EXEC)) == 0:
            os.system('mkdir '+Config.PATH_EXEC)
        #else:
        #    os.system('rm -f '+Config.PATH_EXEC+'/*')

        # print
        # print('|- Creating parameter maps / table for this run...'

        # Create the inputs for ECH2O
        sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=irun)

        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...')
        start = time.time()
        os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time() - start,
              'seconds (limit at ', Config.tlimit, ')')

        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)

        # Store_simulations for elementary effects calculations at the end
        store_sim(Obs, Opti, Config, Site, irun) #_tot)

        # if runOK(Obs, Opti, Config) == 1:
        # Group outputs
        # outputs.store_sim(Obs, Opti, Config, Site, irun)
        # os.system('rm -f *.tab')

        if runOK(Obs, Opti, Config, mode='silent') == 0:
            # else:  # Not running properly? Report

            f_failpar = Config.PATH_OUTmain+'/Parameters_fail.task'+Config.tasknum+'.txt'
            if len(glob.glob(f_failpar)) == 0:
                with open(f_failpar, 'w') as f_in:
                    f_in.write('Trajectory/RadPoint,Sample,'+','.join(Opti.names)+'\n')
            with open(f_failpar, 'a') as f_in:
                f_in.write(str(irun+1)+','+','.join([str(x) for x in Opti.x])+'\n')
            # If it is the very first iteration, record it for later storage
            if irun == 0:
                Opti.begfail = 1

            os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT +
                      '/ech2o.run'+str(irun+1)+'.log')


        # Intermediate clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*.tab')
        os.system('rm -f '+Config.PATH_EXEC+'/Basin*.txt')

        # irun_tot += 1

    # print(Obs.obs)
    # Only for debugging ------------------------------------------------------
    # for oname in Obs.names:
    #     if Obs.obs[oname]['type'] == 'Ts' or Obs.obs[oname]['type'] == 'Total':
    #         Obs.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'
    # ----------------------------------------------------------------------------
    
    # Calculate and output the elementary effects
    ee(Config, Obs, Opti)

    # Clean up after trajectory/radial point
    os.system('rm -f '+Config.PATH_EXEC+'/*')
    os.system('rm -fr '+Config.PATH_EXEC)
                    

def runOK(Obs, Opti, Config, mode='silent'):
    # -- Check if ECH2O ran properly

    isOK = 1

    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0 or \
       os.stat(Config.PATH_EXEC+'/BasinSummary.txt').st_size == 0:
        if(mode == 'verbose'):
            print("Something went wrong, BasinSummary.txt is missing/empty...")
        isOK = 0

    else:
        # for oname in Obs.names:
        #     # print oname
        #     if Obs.obs[oname]['type'] != 'mapStep':
        #         if Obs.obs[oname]['type'] != 'map':
        #             f_test = Config.PATH_EXEC+'/'+Obs.obs[oname]['sim_file']
        #         else:
        #             f_test = Config.PATH_EXEC+'/' + \
        #                 Obs.obs[oname]['sim_file'] + '.map'

        #         if len(glob.glob(f_test)) == 0:
        #             if(mode == 'verbose'):
        #                 print("Something went wrong, no output for "+oname+" !!")
        #                 print('(i.e., '+f_test+' is missing...)')
        #             isOK = 0
        #             # print 'Outputs are there...'

        # 2. Check it ran until the end
        f_test = Config.PATH_EXEC+'/BasinSummary.txt'
        try:
            tmp = np.genfromtxt(f_test, skip_header=1, unpack=True)[0]
        except ValueError:
            isOK = 0
        else:
            # print str(len(tmp2))
            # print str(Obs.lsim)
            if type(tmp) == float or type(tmp) == np.float64 or \
               type(tmp) == np.float32:
                isOK = 0
                if mode == 'verbose':
                    print("Something went wrong, output of length 1 !")
            elif len(tmp) != Obs.lsim:
                # & len(tmp2)==Obs.lsim:
                # print 'It ran until the end !'
                isOK = 0
                if mode == 'verbose':
                    print("Something went wrong, output does not match the " +
                          "supposed sim length!")
                    print('Output: '+str(len(tmp))+' , supposed to be: ' +
                          str(Obs.lsim))
    return isOK


def restart(Config, Opti, Obs):
    # Restart: trim outputs files to match the specified restart iteration

    # -- Get the last iteration that worked
    # use the GOF files

    # File names for grouped simulations
    it = []
    for gof in Opti.GOFs:
        # Get last referenced index (i.e. run number) that worked (without NaN)
        tmp = pd.read_csv(Opti.GOFfiles[gof]).set_index('Sample')
        tmp2 = tmp.dropna()
        if len(tmp2)==0:
            it += [0]
        else:
            it += [tmp2.index[-1]+1]
    # in case this index is different, take minimum
    Config.itres = min(it)

    if Config.itres > 0 :

        # There's something to keep for previous jobs
        print('...but directly restarting from iter. ', Config.itres)

        # Rewrite the GOF files: to evenize between variable last it or
        # remove the final NaN lines (but not those in between "good" runs)
        for gof in Opti.GOFs:
            #print(gof)
            # Get last referenced index (i.e. run number) that worked (without NaN)
            tmp = pd.read_csv(Opti.GOFfiles[gof]).set_index('Sample').loc[1:Config.itres-1]
            tmp.to_csv(Opti.GOFfiles[gof], na_rep="nan")
            if Opti.iscalib2 is True:
                tmp2 = pd.read_csv(Opti.GOFfiles2[gof]).set_index('Sample').loc[1:Config.itres-1]
                tmp2.to_csv(Opti.GOFfiles2[gof], na_rep="nan")

        # Rewrite the BSfile (if needed): to evenize between variable last it or
        # remove the final NaN lines (but not those in between "good" runs)
        if Opti.repBS_interval == 1:
            Config.restart2 = 2
            # Get last referenced index (i.e. run number) that worked (without NaN)
            tmp = pd.read_csv(Opti.BSfile_tb).set_index('Sample').loc[1:Config.itres-1]
            tmp.to_csv(Opti.BSfile_tb, na_rep="nan")
            tmp = pd.read_csv(Opti.BSfile_te).set_index('Sample').loc[1:Config.itres-1]
            tmp.to_csv(Opti.BSfile_te, na_rep="nan")



    else:
        # Nothing worth saving from previous job(s): reinitialize files
        for gof in Opti.GOFs:
            with open(Opti.GOFfiles[gof], 'w') as f_out:
                f_out.write('Sample,'+','.join(Obs.names)+'\n')
            if Opti.iscalib2 is True:
                with open(Opti.GOFfiles2[gof], 'w') as f_out:
                    f_out.write('Sample,'+','.join(Obs.names)+'\n')


        # For RepBS files, it's done elsewhere
        if Opti.repBS_interval == 1:
            Config.restart2 = 1

    # # -- Some cleaning for parameters
    # Config.f_par = Config.PATH_OUT+'/Parameters.txt'
    # # Read grouped simulations
    # tmp = np.genfromtxt(Config.f_par, delimiter=',', skip_header=1,
    #                     max_rows=mxRow)[::, 1::]
    # # Take out the iteration from/beyond the restarting point
    # # #tmp = tmp[0:Config.itres,1::]

    # # Write
    # with open(Config.f_par, 'w') as f_out:
    #     f_out.write('Iteration,'+','.join(Opti.names)+'\n')
    # with open(Config.f_par, 'a') as f_out:
    #     for i in range(mxRow):
    #         f_out.write(str(idx[i])+','+','.join([str(x) for x in
    #                                               list(tmp[i])]) + '\n')

# ==================================================================================
# Parameters-realted functions


def param_get(Opti, Config, options):

    if Config.mode == 'calib_MCruns':

        # Read parameters sample for this instance of array task
        # Open one file for all samples
        Opti.xpar = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                                  unpack=True)[1::]

        tmp = list(pd.read_csv(Config.FILE_PAR, header=None).loc[:, 0])
        if tmp != Opti.names:
            print("ERROR: the full list of parameters are not the same"+
                " in the MCrun definition file and in the input parameter file")
            print('== Dictionary from the MCrun definition file:')
            print(Opti.names)
            print('== Dictionary from stored parameter sets (from MCsampling):')
            print(tmp)
            sys.exit("Aborting calibration")

        # if(Opti.xpar.shape[1] != len(Opti.names)):
        #     sys.exit("The definition file and input parameter file ain't " +
        #              "matching!")
        # print(Opti.xpar.shape)
        Opti.nit = Opti.xpar.shape[0]

    # -- Forward ensemble runs: read directly the params from "best params"
    #    file
    elif Config.mode == 'forward_runs':

        # Sanity check
        tmp = list(pd.read_csv(Config.FILE_PAR, header=None).loc[:, 0])
        if tmp != Opti.names:
            print('== Dictionary from the definition file:')
            print(Opti.names)
            print('== Dictionary from file-stored parameter sets:')
            print(tmp)
            sys.exit("The definition file and input parameter file ain't " +
                     "matching!")

        if(options.OMP_it is None):
            Opti.xpar = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                                      unpack=True)[1::]
        else:
            Opti.xpar = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                                      unpack=True)[1::][None,
                                                        Config.OMP_it-1, :]

    # -- Sensitivity analysis: read directly the params from "best params"
    #    SA_init = 0 from where param_get is called
    elif Config.mode == 'sensi_morris':

        # Sanity check
        tmp = list(pd.read_csv(Config.FILE_PAR, header=None).loc[:, 0])
        if tmp != Opti.names:
            print('== Dictionary from the definition file:')
            print(Opti.names)
            print('== Dictionary from file-stored parameter sets:')
            print(tmp)
            sys.exit("The definition file and input parameter file ain't " +
                     "matching!")

        Opti.xpar = np.genfromtxt(Config.FILE_PAR, delimiter=',',unpack=True)[1::]


# ----------------------------------------------------------------------------
# -- Write parameters values file


def param_store(Opti, Config, it):

    # Open one file for all samples
    if Config.initpar == 0:
        if Config.mode == 'calib_MCruns':
            Config.f_par = Config.PATH_OUTmain+'/Parameters.task' + \
                Config.tasknum+'.txt'
        else:
            Config.f_par = Config.PATH_OUT+'/Parameters.txt'

        if Config.restart == 0:
            with open(Config.f_par, 'w') as f_in:
                f_in.write('Sample,'+','.join(Opti.names)+'\n')
        Config.initpar = 1

    with open(Config.f_par, 'a') as f_in:
        f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')

# ----------------------------------------------------------------------------
# -- Creating/updating inputs for ECH2O


def sim_inputs(Config, Opti, Paras, Site, path_spa, it=0, mode='no_spotpy',
               paramcur=None):

    # Small switch not to write the vegetation params file if
    # there is no need
    sw_veg = 0

    # -- Get the parameter values
    if mode == 'spotpy':
        # In spotpy mode there is no Opti.xpar array with all samples,
        # values are read directly from the spot_setup class
        # In addition, parameter with log10 variation should be "de-log10ned"
        # here!
        Opti.x = []
        for i in range(Opti.nvar):
            if Opti.log[i] == 1:
                # Opti.x += [10**(paramcur()[i][0])]
                Opti.x += [10**paramcur[i]]
            else:
                # Opti.x += [paramcur()[i][0]]
                Opti.x += [paramcur[i]]

    else:
        # Otherwise, just get the samples of the current iteration
        if type(Opti.xpar[it]) is not np.float64:
            Opti.x = np.array(Opti.xpar[it], dtype=np.float64)
        else:
            Opti.x = [Opti.xpar[it]]

    # print(Opti.x)
    # sys.exit()

    for pname in Paras.names:

        # print(pname)
        # print(Paras.ref[pname])

        # - Mapped parameters
        if Paras.ref[pname]['map'] == 1:

            # Map1 unit dependence
            if Paras.ref[pname]['map1'] != 0:
                # print 'Map1 dependent !!'

                # Full map1 dependence
                if Paras.ref[pname]['map1'] == 1:
                    # Start from 0 map
                    outmap = Site.bmaps['unit']*0
                    # Read each soil map unit and apply param value
                    for im in range(Site.nmap1):
                        outmap += Site.bmaps[Site.maps1[im]]*Opti.x[Paras.ind[pname][im]]

                # Heterogeneous map1 unit dependence
                elif type(Paras.ref[pname]['map1']) == list:
                    # print 'Soil dependent !!'

                    # Import reference map (default values)
                    outmap = pcr.readmap(Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
                    # Change the value based on param name correspondance
                    im2 = 0
                    # Only pick the calibrated map/soil units
                    for im in Paras.comp[pname]:
                        outmap = pcr.ifthenelse(Site.bmaps['unit'] == Site.bmaps[Site.maps1[im]],
                                                Opti.x[Paras.ind[pname][im2]], outmap)
                        #outmap += Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im2]]
                        im2 += 1

            # Map2 unit dependence
            elif Paras.ref[pname]['map2'] != 0:
                # print 'Map2 dependent !!'

                # Full map2 dependence
                if Paras.ref[pname]['map2'] == 1:
                    # Start from 0 map
                    outmap = Site.bmaps['unit']*0
                    # Read each soil map unit and apply param value
                    for im in range(Site.nmap2):
                        # print(pname)
                        # print(im)
                        # print(Site.maps2[im])
                        # print(Opti.x[Paras.ind[pname][im]])
                        outmap += Site.bmaps[Site.maps2[im]]*Opti.x[Paras.ind[pname][im]]

                # Heterogeneous map2 unit dependence
                elif type(Paras.ref[pname]['map2']) == list:
                    # print 'Soil dependent !!'

                    # Import reference map (default values)
                    outmap = pcr.readmap(Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
                    # Change the value based on param name correspondance
                    im2 = 0
                    # Only pick the calibrated map/soil units
                    for im in Paras.comp[pname]:
                        outmap = pcr.ifthenelse(Site.bmaps['unit'] == Site.bmaps[Site.maps2[im]],
                                                Opti.x[Paras.ind[pname][im2]], outmap)
                        #outmap += Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im2]]
                        im2 += 1

            # Map3 unit dependence
            elif Paras.ref[pname]['map3'] != 0:
                # print 'Map3 dependent !!'
                
                # Full map3 dependence
                if Paras.ref[pname]['map3'] == 1:
                    # Start from 0 map
                    outmap = Site.bmaps['unit']*0
                    # Read each map3 map unit and apply param value
                    for im in range(Site.nmap3):
                        outmap += Site.bmaps[Site.maps3[im]]*Opti.x[Paras.ind[pname][im]]

                # Heterogeneous map3 unit dependence
                elif type(Paras.ref[pname]['map3']) == list:
                    # print 'Soil dependent !!'

                    # Import reference map (default values)
                    outmap = pcr.readmap(Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
                    # Change the value based on param name correspondance
                    im2 = 0
                    # Only pick the calibrated map/soil units
                    for im in Paras.comp[pname]:
                        outmap = pcr.ifthenelse(Site.bmaps['unit'] == Site.bmaps[Site.maps3[im]],
                                                Opti.x[Paras.ind[pname][im2]], outmap)
                        #outmap += Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im2]]
                        im2 += 1

            # if Site.simRock == 1:
            #     # Taking into account rock/scree: micro-topsoil,
            #     low poros and fixed anisotropy
            #     if pname == 'HLayer1':
            #         outmap = \
            #             outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
            #             Site.bmaps['rock']*0.001
            #     # if pname=='Porosity':
            #     #    outmap = \
            #         # outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
            #         # Site.bmaps['rock']*0.25
            #     if pname == 'Khoriz':
            #         outmap = \
            #             outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
            #             Site.bmaps['rock']*0.000001
            #     if pname == 'Anisotropy':
            #         outmap = \
            #             outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
            #             Site.bmaps['rock']*0.1

            # No spatial
            else:
                # print 'Not map dependent !!'
                # Channel stuff (redundant?)
                # if Paras.ref[pname]['file'] in \
                #    ['chanwidth', 'chanmanningn', 'chanparam']:
                #     outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname][0]]
                # else:
                    # print(type(Opti.x[Paras.ind[pname][0]]))
                outmap = Site.bmaps['unit']*Opti.x[Paras.ind[pname][0]]

            # Create map
            pcr.report(outmap,
                       path_spa+'/'+Paras.ref[pname]['file']+'.map')

            # print('rank', it, ': map of',pname,
            #       'in', path_spa+'/'+Paras.ref[pname]['file']+'.map')

        # - Vegetation parameters
        elif Paras.ref[pname]['veg'] != 0:
            sw_veg = 1
            # Load reference dictionary with lines etc.
            vegnew = copy.copy(Opti.vref)
            # Change the value based on param name correspondance
            # print Opti.vref
            iv2 = 0
            # Only pick the calibrated veg species
            for iv in Paras.comp[pname]:
                vegnew[iv][vegnew['name'].index(pname)] = \
                        str(Opti.x[Paras.ind[pname][iv2]])
                iv2 += 1

        else:
            sys.exit('Error: invalid soil/geol/veg flags to update parameter '+pname)

    # - Write the vegetation parameterization
    if sw_veg == 1:
        # Equalize leaf turnover and additional turnover due to water and/or
        # temperature stress
        # for iv in range(Site.nv):
        #    vegnew[iv][vegnew['name'].index('TurnovL_MWS')] = \
        # copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
        #    vegnew[iv][vegnew['name'].index('TurnovL_MTS')] = \
        # copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
        # Write the vegetation params file (if there's at least one veg param)
        vegfile = open(path_spa+'/'+Site.vfile, 'w')
        vegfile.write('\t'.join(Opti.vref['header'])+'\n')
        for iv in range(Site.nv):
            vegfile.write('\t'.join(vegnew[iv])+'\n')
        vegfile.write('\t'.join(Opti.vref['footer'])+'\n')
        vegfile.write('\t'.join(Opti.vref['name'])+'\n')
        vegfile.close()
    # sys.exit()

    # ------------------------------------------------------------------------------
    # Finalizing the preps....

    # Initial soil water content: needs porosity profile and depth to have
    # consistent initial soil moisture for each layer

    # In any case, get base porosity
    poros = pcr.readmap(path_spa+'/' + Site.f_poros)
    # use it to get number of non-nan pixels
    tmp = pcr.pcr2numpy(poros, np.nan)
    Site.npix = len(list(tmp[~np.isnan(tmp)]))

    # Depending on porosity profile mode, different ways to each layer's porosity
    if Site.poros_mode == 0:
        # Constant profile
        porosL1 = poros
        porosL2 = poros
        porosL3 = poros
    elif Site.poros_mode == 1:
        # Exponential: profile coeff and depths are needed
        kporos = pcr.readmap(path_spa+'/' + Site.f_kporos)
        dL1 = pcr.readmap(path_spa+'/' + Site.f_dL1)
        dL2 = pcr.readmap(path_spa+'/' + Site.f_dL2)
        dTot = pcr.readmap(path_spa+'/' + Site.f_dTot)
        # Layer-integrated values from profile
        porosL1 = kporos*poros*(1-pcr.exp(-dL1/kporos))/dL1
        porosL2 = kporos*poros*(pcr.exp(-dL1/kporos) -
                                pcr.exp(-(dL1+dL2)/kporos))/dL2
        porosL3 = kporos*poros*(pcr.exp(-(dL1+dL2)/kporos) -
                                pcr.exp(-dTot/kporos))/(dTot-dL1-dL2)
    elif Site.poros_mode == 2:
        # Porosity map given for each layer
        porosL1 = poros
        porosL2 = pcr.readmap(path_spa+'/' + Site.f_porosL2)
        porosL3 = pcr.readmap(path_spa+'/' + Site.f_porosL3)

    # Finally, prescribe initial soil moisture depending on user's choice
    if hasattr(Site, 'SWC1init'):
        pcr.report(Site.SWC1init*Site.bmaps['unit'], path_spa+'/'+Site.f_initSWC1)
    else:
        pcr.report(porosL1*Site.frac_SWC1, path_spa+'/'+Site.f_initSWC1)

    if hasattr(Site, 'SWC2init'):
        pcr.report(Site.SWC2init*Site.bmaps['unit'], path_spa+'/'+Site.f_initSWC2)
    else:
        pcr.report(porosL2*Site.frac_SWC2, path_spa+'/'+Site.f_initSWC1)
    if hasattr(Site, 'SWC3init'):
        pcr.report(Site.SWC3init*Site.bmaps['unit'], path_spa+'/'+Site.f_initSWC3)
    else:
        pcr.report(porosL3*Site.frac_SWC3, path_spa + '/' + Site.f_initSWC3)


# ==================================================================================
# Functions for outputs management


def read_sim(Config, Obs, oname, it=0):
    # -- Read a given simulation output (time series only)

    #print(oname)

    # Time series ---------------------------------------------------
    if Obs.obs[oname]['type'] == 'Ts' or Obs.obs[oname]['type'] == 'Total':
        
        # Point-scale
        if Obs.obs[oname]['type'] == 'Ts':
            # print(Obs.obs[oname]['sim_file'])
            # HEader in EcH2O files
            hskip = Obs.nts+3
            # column, factoring in the "custom point sorting" of EcH2O, using
            # the sim_order vector
            idx = np.argsort(np.array(Obs.sim_order))[Obs.obs[oname]['sim_pts']-1]
            # Read
            sim = pd.read_table(Obs.obs[oname]['sim_file'],#error_bad_lines=False,
                            skiprows=hskip, header=None).set_index(0)
            # Check if it ran properly (sometimes some time steps are skipped in *tab...)
            # If steps are missing but the series is long enough, let's replace with nan
            if len(sim) < Obs.lsim:  # and len(sim) > 365:
                idx2 = np.arange(1,Obs.lsim+1)
                sim = sim.reindex(idx2)
                print("Warning: some steps were missing in", Obs.obs[oname]['sim_file'],
                      '(',','.join([str(x) for x in list(pd.isnull(sim[1]).nonzero()[0]+1)]))
                copy_file(Obs.obs[oname]['sim_file'],
                          Config.PATH_OUT+'/'+Obs.obs[oname]['sim_file']+
                          '.run'+str(it+1)+'.txt')


        # Domain-scale (in Basin*Summary.txt) -----------------
        if Obs.obs[oname]['type'] == 'Total':
            # column
            idx = Obs.obs[oname]['sim_pts']-1
            # Read
            sim = pd.read_table(Obs.obs[oname]['sim_file'])  #, error_bad_lines=False)
            # Basin*Summary.txt files don't have an index
            sim = sim.set_axis([str(i)
                                for i in np.arange(1, sim.shape[0]+1)])

        # Get observation column
        sim = sim.iloc[:,idx] * Obs.obs[oname]['sim_conv']

        # If flux in BasinSummary.txt, convert from cumulative to "instantaneous"
        if Obs.obs[oname]['type'] == 'Total' and Obs.obs[oname]['sim_file'] == 'BasinSummary.txt' and \
           Obs.obs[oname]['sim_pts'] in [1, 2, 5, 7, 14, 15, 16, 17, 18, 19, 20,
                                                 21, 22, 23, 24, 25, 26, 27]:
            sim1 = sim.diff() # differences
            sim1.iloc[0] = sim.iloc[0] # restore first value
            sim = copy.copy(sim1)
            
        # print(len(sim),'(2)')
        # Trim (spinup, transient state, etc.)
        sim = sim[Obs.saveB-1:Obs.saveB+Obs.saveL-1]
        #sim = sim.loc[Obs.saveB:Obs.saveB+Obs.saveL-1]
        # print(len(sim),'(3', Obs.saveB, Obs.saveB+Obs.saveL-1,')')
        
        if len(sim) != Obs.saveL:
            print("Warning: problem with "+oname+" trim: read length is " +
                  str(len(sim)) + ' instead of '+str(Obs.saveL))
            sim = [np.nan] * Obs.saveL

        # if Obs.saveB > 1 or Obs.saveL < Obs.lsim:
        #     sim = sim[Obs.saveB-1:Obs.saveB-1+Obs.saveL]
        #     # print(len(sim))
        # #print(sim)

        
    # Time-varying maps (reading one time step) --------------------------
    elif Obs.obs[oname]['type'] == 'mapInit' or Obs.obs[oname]['type'] == 'mapStep' :
        
        # Read map, convert to vector and remove nan pixels
        try:
            sim = pcr.pcr2numpy(pcr.readmap(Obs.obs[oname]['sim_file']),
                                np.nan)*Obs.obs[oname]['sim_conv']
            sim = sim[~np.isnan(sim)]
            #print(sim)
        except RuntimeError:
            print('Warning: RuntimeError - could not read ',
                    Obs.obs[oname]['sim_file'])


    return list(sim)


def store_sim(Obs, Opti, Config, Site, it):
    # -- Store in files for later use

    mode = 'verbose'

    # -- For time series, group the output files in one across simulations,
    #    separating by observations points and veg type where it applies
    for oname in Obs.names:
        if Obs.obs[oname]['type'] != 'mapTs' and \
           (it == 0 or (Opti.begfail == 1 and Config.mode != 'sensi_morris')):
            
            # Historic time series file names
            if Config.mode == 'sensi_morris':
                Obs.obs[oname]['sim_hist'] = Config.PATH_OUTmain+'/'+oname+'_traj'+str(Config.tasknum2)+'.tab'
            else:
                Obs.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'

            # Header of files
            if Obs.obs[oname]['type'] == 'Ts' or Obs.obs[oname]['type'] == 'Total' :
                with open(Obs.obs[oname]['sim_hist'], 'w') as f_out:
                    f_out.write('Sample,'+','.join([str(i+1) for i in
                                                    range(Obs.saveL)])+'\n')
            elif Config.mode == 'sensi_morris':
                # For initial or snapshot maps : only for morris (to get elementary effects)
                with open(Obs.obs[oname]['sim_hist'], 'w') as f_out:
                    f_out.write('Sample,'+','.join([str(i+1) for i in
                                                    range(Site.npix)])+'\n')
                
                
    # Reinit begfail (otherwise will never write all!)
    Opti.begfail = 0

    # Save current run outputs (and delete the file to relieve Maxwell...)
    for oname in Obs.names:

        #print(oname)
        firstMapTs = 1

        if(oname != Obs.names[0]):
            mode = 'silent'

        # Time series ---------------------------------------------------------
        # Integrated variables (in Basin*Summary.txt)

        if Obs.obs[oname]['type'] == 'Total':

            if runOK(Obs, Opti, Config, mode) == 1:

                idx = Obs.obs[oname]['sim_pts']-1
                # Read
                sim = pd.read_table(Obs.obs[oname]['sim_file'],
                                    error_bad_lines=False)
                # Basin*Summary.txt files don't have an index
                sim = sim.set_axis([str(i)
                                    for i in np.arange(1, sim.shape[0]+1)])
                # Get observation column
                sim = sim.iloc[:, idx] * Obs.obs[oname]['sim_conv']

                # If flux, convert from cumulative to "instantaneous"
                if Obs.obs[oname]['sim_file'] == 'BasinSummary.txt' and \
                   Obs.obs[oname]['sim_pts'] in [1, 2, 5, 7, 14, 15, 16, 17, 18, 19, 20,
                                                 21, 22, 23, 24, 25, 26, 27]:
                    sim1 = sim.diff() # differences
                    sim1.iloc[0] = sim.iloc[0] # restore first value
                    sim = copy.copy(sim1) 
                
                # Trim (spinup, transient state, etc.)
                # Index starts at 0
                if Obs.obs[oname]['sim_file'] in ['BasinSummary.txt',
                                                  'BasinVegSummary.txt',
                                                  'BasinAgeSummary.txt',
                                                  'BasincClSummary.txt',
                                                  'Basind2HSummary.txt',
                                                  'Basind18OSummary.txt',
                                                  'BasinfMeltSummary.txt',
                                                  'BasinfGWlat1Summary.txt',
                                                  'BasinfGWlat2Summary.txt',
                                                  'BasinfGWlat3Summary.txt',
                                                  'BasinfGWlatSummary.txt',
                                                  'BasinfGWseep1Summary.txt',
                                                  'BasinfGWseep2Summary.txt',
                                                  'BasinfGWseep3Summary.txt',
                                                  'BasinfGWseepSummary.txt']:
                    sim = sim[Obs.saveB-1:Obs.saveB+Obs.saveL-1]
                else:
                    #sim = sim.loc[Obs.saveB-1:Obs.saveB+Obs.saveL-2]
                    sys.exit('Error: a "Total" obs type should be associated with Basin*Summary file')

                if len(sim) != Obs.saveL:
                    sys.exit("Warning: problem with " + oname +
                             " trim: we've got" +
                             str(len(sim)) + ' instead of '+str(Obs.saveL))

                with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                    f_out.write(str(it+1)+','+','.join([str(j) for j in
                                                        list(sim)])+'\n')
            else:
                # If run failed, write nan line
                with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                    f_out.write(str(it+1)+','+
                                ','.join(list(np.repeat('nan', Obs.saveL)))+
                                '\n')

        # Time series ---------------------------------------------------------
        # Pixel-scale variables (in *tab)

        if Obs.obs[oname]['type'] == 'Ts':

            if runOK(Obs, Opti, Config, mode) == 1:

                hskip = Obs.nts+3
                idx = np.argsort(np.array(Obs.sim_order))[Obs.obs[oname]
                                                          ['sim_pts']-1]

                # Read
                sim = pd.read_table(Obs.obs[oname]['sim_file'],#error_bad_lines=False,
                                    skiprows=hskip, header=None).set_index(0)


                # Check if it ran properly (sometimes some time steps are skipped in *tab...)
                # If steps are missing but the series is long enough, let's replace with nan
                if len(sim) < Obs.lsim:  # and len(sim) > 365:
                    idx2 = np.arange(1,Obs.lsim+1)
                    sim = sim.reindex(idx2)
                    print("Warning: some steps were missing in", Obs.obs[oname]['sim_file'],
                          '(',','.join([str(x) for x in list(pd.isnull(sim[1]).nonzero()[0]+1)]))
                    copy_file(Obs.obs[oname]['sim_file'],
                              Config.PATH_OUT+'/'+Obs.obs[oname]['sim_file']+
                              '.run'+str(it+1)+'.txt')

                # Get observation column
                sim = sim.iloc[:,idx] * Obs.obs[oname]['sim_conv']
                # Trim (spinup, transient state, etc.)
                # Index starts at 1
                sim = sim.loc[Obs.saveB:Obs.saveB+Obs.saveL-1]


                if len(sim) != Obs.saveL:
                    print("Warning: problem with "+oname+" trim: we've got" +
                          str(len(sim)) + ' instead of '+str(Obs.saveL))

                with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                    f_out.write(str(it+1)+','+
                                ','.join([str(j) for j in list(sim)])+'\n')
            else:
                # If run failed, write nan line
                with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                    f_out.write(str(it+1)+','+
                                ','.join(list(np.repeat('nan',Obs.saveL)))+'\n')

        # Fixed-value (initial-value) maps ------------------------------------
        if Obs.obs[oname]['type'] == 'mapInit' or Obs.obs[oname]['type'] == 'mapStep' :

            fail = 0
            
            # Output map name
            if Obs.obs[oname]['type'] == 'mapInit':
                f_m = Config.PATH_EXEC+'/'+Obs.obs[oname]['sim_file']+'.map'
            else:
                f_m = Config.PATH_EXEC+'/'+Obs.obs[oname]['sim_file']

            # Read map
            if len(glob.glob(f_m)) != 0 and runOK(Obs, Opti, Config, mode) == 1 :
                # Now that we have what we need, read the PCraster map...
                var_val = pcr.pcr2numpy(pcr.readmap(f_m), np.nan)*Obs.obs[oname]['sim_conv']
            else :
                print("Warning: the map " + f_m +
                      " seems to be missing from the EcH2O outputs...")
                var_val = np.nan
                fail = 1

            # In forward mode, store netCDF
            if Config.mode == 'forward_runs':
                
                ncFile = Config.PATH_OUT+'/'+oname+'_all.nc'                
                # Missing vaue for PCraster to numpy conversion
                MV = -9999.
                
                # If first run, create file
                if(it == 0):
                    # print('Creating '+ncFile+'...')
                    ncFile = Config.PATH_OUT+'/'+oname+'_all.nc'
                    rootgrp = spio.netcdf_file(ncFile, 'w')
                    rootgrp.createDimension('time', 0)
                    var_y = pcr.pcr2numpy(pcr.ycoordinate(
                        Config.cloneMap), MV)[:, 1]
                    var_x = pcr.pcr2numpy(pcr.xcoordinate(
                        Config.cloneMap), MV)[1, :]
                    rootgrp.createDimension('latitude', len(var_y))
                    rootgrp.createDimension('longitude', len(var_x))
                    rootgrp.createDimension('ensemble', Config.nEns)
                    lat = rootgrp.createVariable('latitude', 'f4', ('latitude',))
                    lat.standard_name = 'Latitude'
                    lat.long_name = 'Latitude cell centres'
                    lon = rootgrp.createVariable('longitude', 'f4', ('longitude',))
                    lon.standard_name = 'Longitude'
                    lon.long_name = 'Longitude cell centres'
                    ens = rootgrp.createVariable('ensemble', 'i', ('ensemble',))
                    ens.standard_name = 'Ensemble'
                    ens.long_name = 'Ensembles of runs'
                    # -assign lat and lon to variables
                    lat[:] = var_y
                    lon[:] = var_x
                    ens[:] = np.arange(Config.nEns)+1
                    # -set netCDF attribute
                    rootgrp.title = 'Maps of '+oname
                    rootgrp.institution = 'Gaia'
                    rootgrp.author = 'P. Camenzind'
                    rootgrp.history = 'Created on %s' % (datetime.now())
                    varStructure = ('latitude', 'longitude', 'ensemble')
                    ncVariable = rootgrp.createVariable(oname, 'f4', varStructure)
                    ncVariable.standard_name = oname
                    # -write to file
                    rootgrp.sync()
                    rootgrp.close()
                                                
                # print('Appending to '+ncFile+'...')
                # Write the actual values for this run
                rootgrp = spio.netcdf_file(ncFile, 'a')
                # - write data
                ncVariable = rootgrp.variables[oname]
                ncVariable[:, :, it] = var_val
                # -update file and close
                rootgrp.sync()
                rootgrp.close()

                
            else: # If Morris mode, store in text file just for elementary effects

                if fail == 0:
                    var_val = var_val[~np.isnan(var_val)]
                    if len(list(var_val)) == Site.npix:
                        with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                            f_out.write(str(it+1)+','+
                                        ','.join([str(j) for j in list(var_val)])+'\n')
                else:
                    # If run failed, write nan line
                    with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                        f_out.write(str(it+1)+','+
                                    ','.join(list(np.repeat('nan', Site.npix)))+'\n')


        # Time-varying maps ---------------------------------------------------
        if Obs.obs[oname]['type'] == 'mapTs':

            # print(oname)

            # Missing value for PCraster to numpy conversion
            MV = -9999.
            lensuf = 8 - len(Obs.obs[oname]['sim_file'])

            MapNames = []
            itOK = []
            itNotOK = []

            for it2 in range(1, Obs.lsim+1):

                # Only save files beyond the spinup/transient period (if any)
                if it2 >= Obs.saveBmap and \
                   it2 < Obs.saveBmap+Obs.saveLmap:

                    # The basic format of maps outputs is XXXXXXXX.xxx
                    # where xxx is the iteration number below 1000.
                    # XXXXXXX (8 characters) concatenate the "base" map name
                    # and the number above 1000

                    suf = ''.join(list(np.repeat('0', lensuf)))+'.' +\
                        format(it2, '03')

                    if(it2 >= 1000):
                        n = np.max([lensuf-1, 0])
                        suf = ''.join(list(np.repeat('0', n))) +\
                            str(it2//1000)+'.'+format(it2 % 1000, '03')

                    if(it2 >= 10000):
                        n = np.max([lensuf-2, 0])
                        suf = ''.join(list(np.repeat('0', n))) +\
                            str(it2//1000)+'.'+format(it2 % 1000, '03')

                    if(it2 >= 100000):
                        n = np.max([lensuf-3, 0])
                        suf = ''.join(list(np.repeat('0', n))) +\
                            str(it2//1000)+'.'+format(it2 % 1000, '03')

                    if(it2 >= 100000):
                        n = np.max([lensuf-4, 0])
                        suf = ''.join(list(np.repeat('0', n))) +\
                            str(it2//1000)+'.'+format(it2 % 1000, '03')

                    # Store names and it2 index
                    # If the length of iteration number conflicts with the
                    # space taken by the "base" name of the map, then the
                    # latter is "eaten" in EcH2O outputting
                    if len(Obs.obs[oname]['sim_file'])+len(suf) > 12:
                        excess = len(Obs.obs[oname]['sim_file'])+len(suf) - 12
                        if len(suf) > 12:
                            print("Warning: the time step is too high to " +
                                  "sort the maps!")
                        f_m = Config.PATH_EXEC+'/' +\
                            Obs.obs[oname]['sim_file'][:-excess]+suf
                    else:
                        f_m = Config.PATH_EXEC+'/' +\
                            Obs.obs[oname]['sim_file'] + suf

                    if len(glob.glob(f_m)) == 0:
                        itNotOK += [it2]
                    else:
                        MapNames += [f_m]
                        itOK += [it2]
                        # print(f_m)

            # Time values for netCDF output
            print(len(Obs.simt))
            var_t = np.array([(Obs.simtmap[x-Obs.saveBmap] -
                               datetime(1901, 1, 1, 0, 0)).days for x in itOK])

            if(len(itOK) == 0):
                print("Nothing to be read from the "+oname+" maps...")
                print("it could be due to an incorrect template map name, or" +
                      "a saveBmap value beyond the last simulation timestep?")
                continue

            if(len(itNotOK) > 0 and firstMapTs == 1):
                firstMapTs = 0
                if(len(itOK) > 0):
                    print("Warning: some of the demanded "+oname +
                          " maps are missing") #, before t=", itOK[0])
                    print("(that's normal if maps output does not start " +
                          "at t=saveBmap or \nif map interval>86400 "+
                          "in EcH2O config file)")
                else:
                    print("Warning: all of the demanded "+oname +
                          " maps are missing!")

            # Second now that we have what we need...
            for it2 in range(len(itOK)):
                # Read map at first time step of interest, convert to array
                # using a missing value,
                # and add an extra 3rd dimension (empty) for later appending
                if(it2 == 0):
                    try:
                        var_val = pcr.pcr2numpy(
                            pcr.readmap(MapNames[it2]),
                            MV)[None, ...]*Obs.obs[oname]['sim_conv']
                    except RuntimeError:
                        print('Warning: RuntimeError - could not read',MapNames[it2])
                # Read subsequent map, same procedure and then append
                else:
                    try:
                        var_val = np.append(var_val,
                                            pcr.pcr2numpy(pcr.readmap(MapNames[it2]),
                                                          MV)[None, ...], axis=0)
                    except RuntimeError:
                        print('Warning: RuntimeError - could not read',MapNames[it2])
                # Clean up
                os.system('rm -f '+MapNames[it2])

            # Write output NCDF file
            ncFile = Config.PATH_OUT+'/'+oname+'_all.nc'
            # print(ncFile)
            # -open nc dataset
            # If first run (that works!), create file
            if(Obs.firstMapTs[oname] == 1):
                ncFile = Config.PATH_OUT+'/'+oname+'_all.nc'
                rootgrp = spio.netcdf_file(ncFile, 'w')
                rootgrp.createDimension('time', 0)
                var_y = pcr.pcr2numpy(
                    pcr.ycoordinate(Config.cloneMap), MV)[:, 1]
                var_x = pcr.pcr2numpy(
                    pcr.xcoordinate(Config.cloneMap), MV)[1, :]
                rootgrp.createDimension('latitude', len(var_y))
                rootgrp.createDimension('longitude', len(var_x))
                rootgrp.createDimension('ensemble', Config.nEns)
                date_time = rootgrp.createVariable('time', 'f8', ('time',))
                date_time.standard_name = 'time'
                date_time.long_name = 'Days since 1901-01-01 00:00:00.0'
                date_time.units = 'Days since 1901-01-01 00:00:00.0'
                date_time.calendar = 'gregorian'
                lat = rootgrp.createVariable('latitude', 'f4', ('latitude',))
                lat.standard_name = 'Latitude'
                lat.long_name = 'Latitude cell centres'
                lon = rootgrp.createVariable('longitude', 'f4', ('longitude',))
                lon.standard_name = 'Longitude'
                lon.long_name = 'Longitude cell centres'
                ens = rootgrp.createVariable('ensemble', 'i', ('ensemble',))
                ens.standard_name = 'Ensemble'
                ens.long_name = 'Ensembles of runs'
                # -assign lat, lon and t to variables
                lat[:] = var_y
                lon[:] = var_x
                date_time[:] = var_t

                ens[:] = np.arange(Config.nEns)+1

                # -set netCDF attribute
                rootgrp.title = 'Maps of '+oname
                rootgrp.institution = 'Gaia'
                rootgrp.author = 'P. Camenzind'
                rootgrp.history = 'Created on %s' % (datetime.now())
                varStructure = ('time', 'latitude', 'longitude', 'ensemble')
                ncVariable = rootgrp.createVariable(oname, 'f4', varStructure)
                ncVariable.standard_name = oname
                # -write to file
                rootgrp.sync()
                rootgrp.close()

                Obs.firstMapTs[oname] = 0

            # Write the actual values for this run
            rootgrp = spio.netcdf_file(ncFile, 'a')
            # - write data
            ncVariable = rootgrp.variables[oname]
            ncVariable[:, :, :, it] = var_val
            # -update file and close
            rootgrp.sync()
            rootgrp.close()

    # print
    # -- Report the full BasinSummary.txt files?
    if Obs.repBS == 1 and Config.mode != 'calib_MCruns':
        os.system('cp '+Config.PATH_EXEC+'/BasinSummary.txt ' +
                  Config.PATH_OUT+'/BasinSummary_run'+str(it+1)+'.txt')

        # Vegetation summary
        if len(glob.glob(Config.PATH_EXEC+'/BasinVegSummary.txt')) != 0:
            os.system('cp '+Config.PATH_EXEC + '/BasinVegSummary.txt ' +
                      Config.PATH_OUT+'/BasinVegSummary_run' +
                      str(it+1)+'.txt')

        if Site.isTrck == 1:
            # deuterium summary
            if len(glob.glob(Config.PATH_EXEC+'/Basind2HSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/Basind2HSummary.txt ' +
                          Config.PATH_OUT+'/Basind2HSummary_run' +
                          str(it+1)+'.txt')
            # oxygen 18 summary
            if len(glob.glob(Config.PATH_EXEC+'/Basind18OSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/Basind18OSummary.txt ' +
                          Config.PATH_OUT+'/Basind18OSummary_run' +
                          str(it+1)+'.txt')
            # chloride summary
            if len(glob.glob(Config.PATH_EXEC + '/BasincClSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasincClSummary.txt ' +
                          Config.PATH_OUT+'/BasincClSummary_run' +
                          str(it+1)+'.txt')
            # age summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinAgeSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinAgeSummary.txt ' +
                          Config.PATH_OUT+'/BasinAgeSummary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfMeltSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfMeltSummary.txt ' +
                          Config.PATH_OUT+'/BasinfMeltSummary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWlat1Summary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWlat1Summary.txt ' +
                          Config.PATH_OUT+'/BasinfGWlat1Summary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWlat2Summary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWlat2Summary.txt ' +
                          Config.PATH_OUT+'/BasinfGWlat2Summary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWlat3Summary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWlat3Summary.txt ' +
                          Config.PATH_OUT+'/BasinfGWlat3Summary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWlatSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWlatSummary.txt ' +
                          Config.PATH_OUT+'/BasinfGWlatSummary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWseep1Summary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWseep1Summary.txt ' +
                          Config.PATH_OUT+'/BasinfGWseep1Summary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWseep2Summary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWseep2Summary.txt ' +
                          Config.PATH_OUT+'/BasinfGWseep2Summary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWseep3Summary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWseep3Summary.txt ' +
                          Config.PATH_OUT+'/BasinfGWseep3Summary_run' +
                          str(it+1)+'.txt')
            # fMelt summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinfGWseepSummary.txt')) != 0:
                os.system('cp '+Config.PATH_EXEC + '/BasinfGWseepSummary.txt ' +
                          Config.PATH_OUT+'/BasinfGWseepSummary_run' +
                          str(it+1)+'.txt')

def store_BS_interval(Obs, Opti, Config, it):
    # -- Store in files for later use
    # Integrated variables (in BasinSummary.txt)

    # Read
    sim = pd.read_table(Config.PATH_EXEC+'/BasinSummary.txt',
                        error_bad_lines=False).loc[:, 'Precip':'MBErr']
    #print(sim)
    # Basin*Summary.txt files don't have an index
    sim = sim.set_axis([str(i) for i in np.arange(1, sim.shape[0]+1)])

    # Header
    if (it == 0 or Opti.begfail == 1) and (Config.restart == 0 or
                                           (Config.restart==1 and Config.restart2==1)):
        with open(Opti.BSfile_tb, 'w') as f_out:
            f_out.write('Sample,'+','.join([j for j in sim.columns])+'\n')
        with open(Opti.BSfile_te, 'w') as f_out:
            f_out.write('Sample,'+','.join([j for j in sim.columns])+'\n')


    # Write values at given intervals
    if runOK(Obs, Opti, Config, 'silent') == 1:

        # Get observation lines, at beginning and end
        sim_tb = sim.iloc[Obs.saveB-1, :]
        sim_te = sim.iloc[Obs.lsim-1, :]
        #print(sim_tb)
        #print(sim_te)

        # Simulation, or nan
        with open(Opti.BSfile_tb, 'a') as f_out:
            f_out.write(str(it+1)+','+
                        ','.join([str(j) for j in list(sim_tb)])+'\n')
        with open(Opti.BSfile_te, 'a') as f_out:
            f_out.write(str(it+1)+','+
                        ','.join([str(j) for j in list(sim_te)])+'\n')
    else:
        # If run failed, write nan line
        with open(Opti.BSfile_tb, 'a') as f_out:
            f_out.write(str(it+1)+','+
                        ','.join(list(np.repeat('nan', 22)))+'\n')
        with open(Opti.BSfile_te, 'a') as f_out:
            f_out.write(str(it+1)+','+
                        ','.join(list(np.repeat('nan', 22)))+'\n')
# -----------------------------------------------------------------------


def store_GOF(Obs, Opti, Config, Site, it):
    # -- Store goodness-of-fit using using several metrics:
    #    NSE, KGE, RMSE, MAE...

    # Did it run OK?
    if runOK(Obs, Opti, Config, mode='verbose') == 1:
        # # Read outputs
        # if Obs.nobs == 1:
        #     simulations = read_sim(Config, Obs, Obs.names[0], it)
        #     # if(len(simulations) < Obs.saveL):
        #     #     print('sim length:', len(simulations),
        #     # ', expected:', Obs.saveL)
        #     #     simulations = [np.nan] * Obs.saveL
        # else:
        #     simulations = np.full((Obs.nobs, Obs.saveL), np.nan).tolist()
        #     for i in range(Obs.nobs):
        #         oname = Obs.names[i]
        #         simulations[i][:] = read_sim(Config, Obs, oname, it)
        #         # if(len(simulations[i]) < Obs.saveL):
        #         #     print('sim length:', len(simulations[i]),
        #         # ', expected:', Obs.saveL)
        #         #     simulations[i][:] = [np.nan] * Obs.saveL

        # Opti.tmp = copy.copy(simulations)

        for i in range(Obs.nobs):

            oname = Obs.names[i]
            # Read outputs
            simulation = read_sim(Config, Obs, oname, it)

            ncalib = 1
            
            if Obs.obs[oname]['type'] != 'mapStep':
                # Check if there is one or two calibration periods (except for snapshot maps)
                if Opti.calib2[oname] is True:
                    ncalib = 2

            for ic in range(ncalib):

                # Have obervation and simulations matching the same time period
                # obs: already pre-processed
                if ic == 0:
                    if Obs.obs[oname]['type'] != 'mapStep':
                        tobs = pd.to_datetime(Opti.obs[oname]['Date'].values)
                        # o = np.asanyarray(Opti.obs[oname]['value'].values)
                        o = Opti.obs[oname].reset_index(drop=True)
                        # if Opti.gof_d1 is True: # Derivative
                        #     o_d1 = np.asanyarray(Opti.obs[oname+'_d1']['value'].values)
                    else:
                        o = pd.Series(Opti.obs[oname]).reset_index(drop=True)
                if ic == 1: # ic always =0 for type=mapStep
                    tobs = pd.to_datetime(Opti.obs2[oname]['Date'].values)
                    # o = np.asanyarray(Opti.obs2[oname]['value'].values)
                    o = Opti.obs2[oname].reset_index(drop=True)
                    # if Opti.gof_d1 is True: # Derivative
                    #     o_d1 = np.asanyarray(Opti.obs2[oname+'_d1']['value'].values)

                # First step for sim: trim sim to obs timespan
                # + only keep dates with obs (even if nan)
                # print(sim)
                # print(sim.shape)
                if Obs.obs[oname]['type'] != 'mapStep': # Time series
                    s = pd.Series([simulation[j] for j in range(Obs.saveL)
                                   if Obs.simt[j] in tobs]).reset_index(drop=True)
                else: # Map (already with right number of non-NA pixels
                    s = pd.Series(list(simulation)).reset_index(drop=True)

                if Opti.gof_d1 is True: # Derivative
                    dt = tobs.to_series().reset_index(drop=True).diff().dt.total_seconds()[1::]
                    o_d1 = o['value'].diff()[1::] / dt
                    s_d1 = s.diff()[1::] / dt
                    
                # print(s)
                # print(o)
                # print(s_d1)
                # #sys.exit()

                # Second step (both o and s): remove nan due to gaps in obs
                # (or missing steps in sims...)
                if Obs.obs[oname]['type'] != 'mapStep':
                    tmp = s.values*o['value'].values
                    s = np.asanyarray([s.values[k] for k in range(len(tmp)) if not
                                       np.isnan(tmp[k])])
                    o = np.asanyarray([o['value'].values[j] for j in range(len(tmp)) if not
                                       np.isnan(tmp[j])])
                else:
                    tmp = s.values*o.values
                    s = np.asanyarray([s.values[k] for k in range(len(tmp)) if not
                                       np.isnan(tmp[k])])
                    o = np.asanyarray([o.values[j] for j in range(len(tmp)) if not
                                       np.isnan(tmp[j])])

                # print(oname)
                # print('obs')
                # print(o)
                # print('sim')
                # print(s)
                
                # Prepare lists of GOFs
                if i == 0 and ic == 0:
                    gofs = {}
                    for gof in Opti.GOFs:
                        gofs[gof] = []
                if i == 0 and ic == 1:
                    gofs2 = {}
                    for gof in Opti.GOFs:
                        gofs2[gof] = []

                # Another sanity check: any data/sim left after nan screening?
                if s.__len__() == 0 or o.__len__() == 0:
                    print('Warning: nothing to compare to after trimming for '+oname+'!')
                    # Add nan
                    for gof in Opti.GOFs:
                        if ic == 0:
                            gofs[gof] += [np.nan]
                        if ic == 1:
                            gofs2[gof] += [np.nan]

                else:
                    # Now use your favorite likelihood estimator for each obs type
                    for gof in Opti.GOFs:
                        if gof == 'corr':  # pearson correlation coefficient
                            tmp = GOFs.corr(s, o)
                        elif gof == 'sprmn': # separman's rank correlation
                            tmp = np.array(pd.DataFrame(data={'x':s, 'y':o}).corr(method='spearman'))[0,1]
                        elif gof == 'rstd':  # ratio of standard deviations
                            tmp = GOFs.rstd(s, o)
                        elif gof == 'dstd':  # difference of standard deviations
                            tmp = GOFs.dstd(s, o)
                        elif gof == 'rmu':  # ratio of mean
                            tmp = GOFs.rmu(s, o)
                        elif gof == 'bias':  # bias
                            tmp = GOFs.bias(s, o)
                        elif gof == 'NSE':  # NSE
                            tmp = GOFs.nash_sutcliffe(s, o)
                        elif gof == 'KGE':  # KGE 2009
                            tmp = GOFs.kling_gupta(s, o)
                        elif gof == 'KGE2012':  # KGE 2012
                            tmp = GOFs.kling_gupta(s, o, method='2012')
                        elif gof == 'RMSE':  # RMSE
                            tmp = GOFs.rmse(s, o)
                        elif gof == 'MAE':  # MAE
                            tmp = GOFs.meanabs(s, o)
                        elif gof == 'KGEc':  # mean-centered KGE
                            tmp = GOFs.kling_gupta(s-np.mean(s), o-np.mean(o))
                        elif gof == 'KGE2012c':  # mean-centered KGE2012
                            tmp = GOFs.kling_gupta(s-np.mean(s),o-np.mean(o), method='2012')
                        elif gof == 'RMSEc':  # mean-centered RMSE
                            tmp = GOFs.rmse(s-np.mean(s), o-np.mean(o))
                        elif gof == 'MAEc':  # mean-centered MAE
                            tmp = GOFs.meanabs(s-np.mean(s), o-np.mean(o))
                        elif gof == 'RMSE_d1':  # RMSE
                            tmp = GOFs.rmse(s_d1, o_d1)
                        elif gof == 'MAE_d1':  # RMSE
                            tmp = GOFs.meanabs(s_d1, o_d1)
                        elif gof == 'gauL':  # Gaussian likelihood, measurement error out
                            tmp = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s)
                        elif gof == 'logL':  # Log likelihood
                            tmp = spotpy.likelihoods.logLikelihood(o, s)
                        elif gof == 'gauLerr':  # Gaussian likelihood
                            tmp = spotpy.likelihoods.gaussianLikelihoodHomoHeteroDataError(o, s)
                        else:
                            sys.exit('Error, this GOF is not treated (yet)')
                        # Add to the gof, depending on calibration period
                        if ic == 0:
                            gofs[gof] += [tmp]
                        if ic == 1:
                            gofs2[gof] += [tmp]

        # Store goodnesses of fit, one files per GOF
        for gof in Opti.GOFs:
            # GOF
            with open(Opti.GOFfiles[gof], 'a') as f_out:
                f_out.write(str(it+1)+',' +
                            ','.join([str(x) for x in gofs[gof]])+'\n')
            # Where applicable, GOF second period
            if Opti.iscalib2 is True:
                with open(Opti.GOFfiles2[gof], 'a') as f_out:
                    f_out.write(str(it+1)+',' +
                                ','.join([str(x) for x in gofs2[gof]])+'\n')

    else:
        # Write NaN
        for gof in Opti.GOFs:
            with open(Opti.GOFfiles[gof], 'a') as f_out:
                f_out.write(str(it+1)+',' +
                            ','.join(list(np.repeat('nan', Obs.nobs)))+'\n')
            # Where applicable, nan is GOF second period
            if Opti.iscalib2 is True:
                with open(Opti.GOFfiles2[gof], 'a') as f_out:
                    f_out.write(str(it+1)+',' +
                                ','.join(list(np.repeat('nan', Obs.nobs2)))+'\n')


# ==================================================================================
# Functions for variable sampling


def MC_sample(Opti, Config):
    # -- Generate random set of parameters values of write it

    # Opti.xtot = np.arange(1,Opti.nsamptot+1)

    # Latin Hypercube sampling
    if Config.sampling in ['LHS', 'LHS_m', 'LHS_r', 'LHS_c']:

        print('...using a latin hypercube sampling...')

        # 'Normal' LHS
        if Config.sampling == 'LHS':
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(pyDOE.lhs(Opti.nvar, samples=Opti.nsamptot))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling == 'LHS_m':
            print('...with maximin criterion -- it will take longer and a ' +
                  'lot of memory')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(pyDOE.lhs(Opti.nvar, samples=Opti.nsamptot,
                                         criterion='m'))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling == 'LHS_r':
            print('...with correlation criterion -- it will ' +
                  'take longer and a lot of memory')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(pyDOE.lhs(Opti.nvar, samples=Opti.nsamptot,
                                         criterion='corr'))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling == 'LHS_c':
            print('...with center criterion')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(pyDOE.lhs(Opti.nvar, samples=Opti.nsamptot,
                                         criterion='center'))

        print('...LHS matrix generated...')

        # Second, scale with the actual range
        for i in range(Opti.nvar):
            # Log transform where needed
            if Opti.log[i] == 1:
                tmp = 10**(mat[i]*np.log10(Opti.max[i]/Opti.min[i]) +
                           np.log10(Opti.min[i]))
            else:
                tmp = mat[i]*(Opti.max[i]-Opti.min[i]) + Opti.min[i]

            if i == 0:
                Opti.xtot = copy.copy(tmp)
            else:
                Opti.xtot = np.vstack((Opti.xtot, tmp))

    # -- Uniform distribution
    elif Config.sampling == 'uniform':

        print('...using uniform distributions...')

        for i in range(Opti.nvar):
            # Log transform where needed
            if Opti.log[i] == 1:
                if i == 0:
                    Opti.xtot = 10**(np.random.uniform(np.log10(Opti.min[i]),
                                                       np.log10(Opti.max[i]),
                                                       Opti.nsamptot))
                else:
                    Opti.xtot = \
                        np.vstack((Opti.xtot,
                                   10**(np.random.uniform(
                                       np.log10(Opti.min[i]),
                                       np.log10(Opti.max[i]),
                                       Opti.nsamptot))))
            else:
                if i == 0:
                    Opti.xtot = np.random.uniform(Opti.min[i], Opti.max[i],
                                                  Opti.nsamptot)
                else:
                    Opti.xtot = np.vstack((Opti.xtot,
                                           np.random.uniform(Opti.min[i],
                                                             Opti.max[i],
                                                             Opti.nsamptot)))
        # print i, Opti.names[i], Opti.x[i]

    else:
        sys.exit('No proper sampling method selected!')

    print('Parameters sampling done!')

    # -- Reshape for output
    # Opti.xtot = np.transpose(Opti.xtot)#.shape((Opti.nvar,Opti)
    # print(Opti.xtot.shape)

    print
    print('Writing in '+Config.ncpu+' files...('+str(Opti.nit)+' sets each)')

    # -- Write one file per parallel job
    for i in range(int(Config.ncpu)):
        f_out = Config.FILE_PAR+str(i+1)+'.txt'
        k = i*Opti.nit
        # print str(k)+','+str(k+Opti.nit)
        with open(f_out, 'w') as fw:
            for j in range(Opti.nvar):
                # print i
                # print j
                # print k
                tmp = [a for a in Opti.xtot[j][k:k+Opti.nit]]
                # print len(tmp)
                fw.write(Opti.names[j]+',' +
                         ','.join([str(a) for a in tmp])+'\n')

    # Write on file giving parameters range, log...(for later plots)
    f_out = Config.FILE_PAR+'char.txt'
    with open(f_out, 'w') as fw:
        # fw.write('Sample,'+','.join(Opti.names)+'\n')
        fw.write('Names,Min,Max,Log\n')
        for i in range(Opti.nvar):
            fw.write(','.join([Opti.names[i],str(Opti.min[i]),str(Opti.max[i]),str(Opti.log[i])])+'\n')

    print('')

    # np.savetxt(f_out,np.transpose(Opti.xtot))#,header=

# ==================================================================================
# Functions for Morris sensitivity analysis


def trajs(Config, Opti):
    # -- Creation of the Morris trajectories
    # Follows methodology and recommandations of
    # Sohier, Farges, and Piet-Lahanier (2014), Improvements of the
    # representativity of the Morris
    # method for air-launch-to-orbit separation, Proc. 19th congress of IFAC.

    # Number of levels -> equal number of trajectories/radial points
    # This derives from a recommendation of Sohier et al. (IFAC, 2014)
    Opti.nlev = copy.copy(Opti.nr)

    # vals = {}

    # Step: plus-minus 0.5 of the range of each parameter
    Opti.step = np.zeros((Opti.nvar), np.float32)

    # Construct B* matrix, for each repetition
    Opti.Bnorm = np.zeros((Opti.nvar, Opti.nvar+1, Opti.nr),
                          np.float32)  # the normalized
    Opti.Bstar = np.zeros((Opti.nvar, Opti.nvar+1, Opti.nr),
                          np.float32)  # the final used in runs

    # Starting point: latin hypercube sampling, maximizing 'distance'
    # between point, and centered in the nvar intervals
    Opti.Bnorm[:, 0, :] = np.transpose(pyDOE.lhs(Opti.nvar,
                                                 samples=Opti.nr,
                                                 criterion='cm'))
    # print(Opti.nvar, Opti.nvar+1, Opti.nr)
    # print(Opti.Bnorm[:, 0, :])
    # print('-----------------------------------------')

    # Construct samples
    for ir in range(Opti.nr):
        for iv in range(Opti.nvar):

            # Mode 1 : trajectories
            # Mode 2 : radial points
            # In both cases the other of one-at-a-time change is fixed:
            # 1st param changes,
            # then 2nd param, etc.
            # the randomness is assured by the initial LHS + the fixed step of
            # +-0.5

            if(Opti.MSspace == 'trajectory'):
                # copy previous location
                Opti.Bnorm[:, iv+1, ir] = copy.copy(Opti.Bnorm[:, iv, ir])
            elif(Opti.MSspace == 'radial'):
                # alway start from initial
                Opti.Bnorm[:, iv+1, ir] = copy.copy(Opti.Bnorm[:, 0, ir])
            else:
                sys.exit('Wrong option for the MS parameter space definition!')

            # Successive changes with + (resp. -) 0.5, depending on if they
            # are above (resp. below) the mid-interval
            if Opti.Bnorm[iv, iv+1, ir] < 0.5:
                Opti.Bnorm[iv, iv+1, ir] += 0.5
            elif Opti.Bnorm[iv, iv+1, ir] >= 0.5:
                Opti.Bnorm[iv, iv+1, ir] -= 0.5

            # Check for error
            if Opti.Bnorm[iv, iv+1, ir] > 1 or \
               Opti.Bnorm[iv, iv+1, ir] < 0:
                print('Error in the incrementation of the parameter',
                      Opti.names[iv])
                print(1/(2*Opti.nlev), Opti.Bnorm[iv, iv, ir],
                      Opti.Bnorm[iv, iv+1, ir], 1-1/(2*Opti.nlev))
                sys.exit()

    # Construct the actual Bstar, with non-normalized values
    for iv in range(Opti.nvar):
        if Opti.log[iv] == 1:
            Opti.Bstar[iv, :, :] = 10**(Opti.Bnorm[iv, :, :]*np.log10(
                Opti.max[iv]/Opti.min[iv])+np.log10(Opti.min[iv]))
            Opti.step[iv] = \
                0.5 * (np.log10(Opti.max[iv])-np.log10(Opti.min[iv]))
        else:
            Opti.Bstar[iv, :, :] = \
                Opti.Bnorm[iv, :, :]*(Opti.max[iv]-Opti.min[iv]) + Opti.min[iv]
            Opti.step[iv] = 0.5 * (Opti.max[iv]*Opti.min[iv])

    # # Reconstruct step (+- 0.5)
    # if(Opti.MSspace == 'trajectory'):
    #     Opti.dx = np.diff(Opti.Bstar, axis=1)
    #     #Opti.dx = np.diff(Opti.Bnorm, axis=1)
    # elif(Opti.MSspace == 'radial'):
    #     Opti.dx = Opti.Bstar[:, 1::, :] - Opti.Bstar[: ,0 , :][:, None, :]
    #     #Opti.dx = Opti.Bnorm[:, 1::, :] - Opti.Bnorm[:, 0, :][:, None, :]            
    # Opti.dx[Opti.dx > 0] = 0.5
    # Opti.dx[Opti.dx < 0] = -0.5

    # if(np.ptp(Opti.dx) != 1.0 or np.min(Opti.dx) != -0.5 or
    #    np.max(Opti.dx) != 0.5):
    #     print(np.ptp(Opti.dx), np.min(Opti.dx), np.max(Opti.dx))
    #     sys.exit('Error: Bnorm has a problem...')
        
    # Outputs of the parameter trajectory / radial points
    #print(Opti.Bstar.shape)
    print('Writing in '+str(Opti.nr)+' trajectory files...('+str(Opti.nvar+1)+' sets each)')
    
    #Bstart2 = pd.DataFrame(Opti.Bstar)
    #print(Opti.Bstar2)
    for ir in range(Opti.nr):
        trajnb = str(ir+1)  # '%02i' % int(ir+1)
        # print trajnb

        ## Write normalized values [0, 1] of sampling interval
        with open(Config.FILE_PAR+trajnb+'_norm.txt', 'w') as fw:
            for j in range(Opti.nvar):
                tmp = [a for a in Opti.Bnorm[j, :, ir]]
                #print(tmp)
                fw.write(Opti.names[j]+','+','.join([str(a) for a in tmp])+'\n')


        ## Write actual values within the sampling interval
        with open(Config.FILE_PAR+trajnb+'_values.txt', 'w') as fw:
            for j in range(Opti.nvar):
                tmp = [a for a in Opti.Bstar[j, :, ir]]
                #print(tmp)
                fw.write(Opti.names[j]+','+','.join([str(a) for a in tmp])+'\n')

    # Write on file giving parameters range, log...(for later plots)
    f_out = Config.FILE_PAR+'_attributes.txt'
    with open(f_out, 'w') as fw:
        fw.write('Names,Min,Max,StepN,Log\n')
        for i in range(Opti.nvar):
            fw.write(','.join([Opti.names[i],str(Opti.min[i]),str(Opti.max[i]),
                               str(Opti.step[i]),str(Opti.log[i])])+'\n')

    print('')

# ----------------------------------------------------------------------------
# -- Calculation of elementary effects for Morris Sensitivity Analysis
# Applied to one trajectory / radial points with same origin
# Since time series are the outputs, metrics for d(outputs) are bias and RMSE

# !! Pandas unavailable so using numpy (unhandy) numpy for I/O+data-manip
# EDIT : pandas available, has to be reverted to it
# (in theory just uncomment previous code, but testing would be necessary)

# NEW: it's used with summary time series/vectorized maps which are specific for
# one trajcory / set of radial points (and not over all runs *nr), hence itraj
# does not matter

def ee(Config, Obs, Opti):

    firstObs = 0
    numObs = 0
    outObs = []#['Parameter']

    for oname in Obs.names:

        # Time series or snapshot maps
        # print(oname)
        
        # Read file        
        f_in = Obs.obs[oname]['sim_hist']
        df_sim = pd.read_csv(f_in).set_index('Sample')#.iloc[itraj*(Opti.nvar+1)::, ]

        #print(df_sim)
        
        # Diff between sims
        if Opti.MSspace == 'trajectory':
            df_diff = df_sim.diff().iloc[1::, ]
        if Opti.MSspace == 'radial':
            df_diff = df_sim.iloc[1::]-df_sim.iloc[0]

        # print(df_diff.shape)
        # Normalization factor: 80%-spread of "reference" simulations
        qt90 = np.quantile(df_sim.iloc[0], 0.9)
        qt10 = np.quantile(df_sim.iloc[0], 0.1)
        ynorm = qt90-qt10
        if ynorm == 0:
            if qt90 != 0:
                ynorm = qt90
            else:
                ynorm = 1
            
        bias = df_diff.mean(axis=1) / ynorm
        # print(bias.shape)
        # print(Opti.stepN.shape)
        # Get RMSE, normalized by 80% spread of the (origin) time series
        # qt90 = df_sim.set_index('Sample').iloc[0].quantile(0.9)
        # qt10 = df_sim.set_index('Sample').iloc[0].quantile(0.1)
        # qt90 = df_sim.set_index('Sample').quantile(0.9, axis=1)
        # qt10 = df_sim.set_index('Sample').quantile(0.1, axis=1)
        RMSE = np.sqrt((df_diff**2).mean(axis=1)) / ynorm
        # Get corresponding (normalized) elementary effect
        bias_ee = bias / Opti.stepN
        RMSE_ee = RMSE / Opti.stepN

        # Associate the corresponding parameter being tested
        # (the order if the basic parameter order)
        bias_ee.index = Opti.names
        RMSE_ee.index = Opti.names
        
        # sim = np.genfromtxt(f_in, delimiter=',', skip_header=1,
        #                     unpack=True)[1:Config.trimL+1, :]
        
        # # Take into account accumulated fluxes

        # # Diff between sims
        # if(Opti.MSspace == 'trajectory'):
        #     simd = np.diff(sim)
        # elif(Opti.MSspace == 'radial'):
        #     simd = sim[:, 1::] - sim[:, 0][..., None]
        
        # # Elementary effect (keep direction of param change)
        # bias_ee = np.zeros((Opti.nvar), np.float32)*np.nan
        # RMSE_ee = np.zeros((Opti.nvar), np.float32)*np.nan
        # for i in range(Opti.nvar):
        #     bias_ee[i] = np.mean(simd[:, i]) / Opti.dx[i, i]
        #     RMSE_ee[i] = np.sqrt(np.mean(simd[:, i]**2)) / Opti.dx[i, i]
        # # bias_ee = bias / Opti.BnormstepN
        # # RMSE_ee = RMSE / Opti.stepN
        
        # Build the overall data frame
        if(firstObs == 0):
            # bias_ee_tot = bias_ee[..., None]  # Creates a (..,1) dimension
            # RMSE_ee_tot = RMSE_ee[..., None]  # Creates a (..,1) dimension
            bias_ee_tot = pd.DataFrame(bias_ee) #.assign(oname=bias_ee)
            RMSE_ee_tot = pd.DataFrame(RMSE_ee) #.assign(oname=RMSE_ee)
        else:
            # bias_ee_tot = np.append(bias_ee_tot, bias_ee[..., None], 1)
            # RMSE_ee_tot = np.append(RMSE_ee_tot, RMSE_ee[..., None], 1)
            bias_ee_tot[str(numObs)] = bias_ee
            RMSE_ee_tot[str(numObs)] = RMSE_ee

        # Update
        firstObs = 1
        # Increment number of obs actually evaluated
        numObs += 1
        # Append name of obs actually evaluated
        outObs = outObs + [oname]

    # Write outputs -----------------------------------------------------------
    bias_ee_tot.columns = outObs
    RMSE_ee_tot.columns = outObs

    if(Opti.MSspace == 'trajectory'):
        bias_ee_tot.to_csv(Config.PATH_OUTmain+'/EE.Traj'+Config.tasknum2+'.bias.txt')
        RMSE_ee_tot.to_csv(Config.PATH_OUTmain+'/EE.Traj'+Config.tasknum2+'.RMSE.txt')
        # with open(Config.FILE_EE+'.EE.Traj'+str(itraj+1) +
        #           '.bias.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(bias_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

        # with open(Config.FILE_EE+'.EE.Traj'+str(itraj+1) +
        #           '.RMSE.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(RMSE_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

    if(Opti.MSspace == 'radial'):
        bias_ee_tot.to_csv(Config.PATH_OUTmain+'/EE.RadP'+Config.tasknum2+'.bias.txt')
        RMSE_ee_tot.to_csv(Config.PATH_OUTmain+'/EE.RadP'+Config.tasknum2+'.RMSE.txt')
        # with open(Config.FILE_EE+'.EE.RadP'+str(itraj+1) +
        #           '.bias.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(bias_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

        # with open(Config.FILE_EE+'.EE.RadP'+str(itraj+1) +
        #           '.RMSE.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(RMSE_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

# ==================================================================================

# Function(s) for multi-objective calibration with spotpy


def MultiObj(obs, sim, Obs, Opti, w=False):

    like = 0
    Ltot = []

    # Sanity check: did the simulation run?
    # (necessary to avoid error in the subsequent postprocessing)
    if sim is None:
        if Obs.nobs > 1:
            like = [np.nan for i in range(Obs.nobs+1)]
            # Obs.nobs+1 because of total objfunc + indiv ones
        else:
            like = np.nan

    else:
        for i in range(Obs.nobs):

            oname = Obs.names[i]

            # Have obervation and simulations matching the same time period
            # obs: already pre-processed
            tobs = pd.to_datetime(obs[oname]['Date'].values)
            o = np.asanyarray(obs[oname]['value'].values)

            # First step for sim: trim sim to obs timespan
            # + only keep dates with obs (even if nan)
            # print(sim)
            # print(sim.shape)
            if Obs.nobs == 1:
                s = np.asanyarray([sim[j] for j in range(Obs.saveL) if
                                   Obs.simt[j] in tobs])
            else:
                s = np.asanyarray([sim[i][j] for j in range(Obs.saveL) if
                                   Obs.simt[j] in tobs])

            # Second step (both o and s): remove nan due to gaps in obs
            tmp = s*o
            s = np.asanyarray([s[k] for k in range(len(tmp)) if not
                               np.isnan(tmp[k])])
            o = np.asanyarray([o[j] for j in range(len(tmp)) if not
                               np.isnan(tmp[j])])

            # Another sanity check: any data/sim left after nan screening?
            if s.__len__() == 0 or o.__len__() == 0:
                L = np.nan

            else:
                # Now use your favorite likelihood estimator for each obs type

                if obs[oname]['gof'] == 'corr':  # pearson correlation coefficient
                    L = GOFs.corr(s, o)
                elif obs[oname]['gof'] == 'sprmn': # spearman's rank correlation
                    L = np.array(pd.DataFrame(data={'x':s, 'y':o}).corr(method='spearman'))[0,1]
                elif obs[oname]['gof'] == 'rstd':  # ratio of standard deviations
                    L = GOFs.rstd(s, o)
                elif obs[oname]['gof'] == 'rmu':  # ratio of mean
                    L = GOFs.rmu(s, o)
                elif obs[oname]['gof'] == 'bias':  # bias
                    L = GOFs.bias(s, o)
                elif obs[oname]['gof'] == 'NSE':  # NSE
                    L = GOFs.nash_sutcliffe(s, o)
                elif obs[oname]['gof'] == 'KGE':  # KGE 2009
                    L = GOFs.kling_gupta(s, o)
                elif obs[oname]['gof'] == 'KGE2012':  # KGE 2012
                    L = GOFs.kling_gupta(s, o, method='2012')
                elif obs[oname]['gof'] == 'RMSE':  # RMSE
                    L = GOFs.rmse(s, o)
                elif obs[oname]['gof'] == 'MAE':  # MAE
                    L = GOFs.meanabs(s, o)
                elif obs[oname]['gof'] == 'KGEc':  # mean-centered KGE
                    L = GOFs.kling_gupta(s-np.mean(s), o-np.mean(o))
                elif obs[oname]['gof'] == 'KGE2012c':  # mean-centered KGE2012
                    L = GOFs.kling_gupta(s-np.mean(s),o-np.mean(o), method='2012')
                elif obs[oname]['gof'] == 'RMSEc':  # mean-centered RMSE
                    L = GOFs.rmse(s-np.mean(s), o-np.mean(o))
                elif obs[oname]['gof'] == 'MAEc':  # mean-centered MAE
                    L = GOFs.meanabs(s-np.mean(s), o-np.mean(o))
                elif obs[oname]['gof'] == 'gauL':  # Gaussian likelihood, measurement error out
                    L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s)
                elif obs[oname]['gof'] == 'logL':  # Log likelihood
                    L = spotpy.likelihoods.logLikelihood(o, s)
                elif obs[oname]['gof'] == 'gauLerr':  # Gaussian likelihood
                    L = spotpy.likelihoods.gaussianLikelihoodHomoHeteroDataError(o, s)
                else:
                    sys.exit('ERROR: unvalid gof for ', oname)

                # # Specific treatment for different obs types?
                # if oname.split('_')[0] in ['GWD', 'WTD', 'GWL', 'WTL']:
                #     # print(oname, 'remove mean')
                #     # s -= np.mean(s)
                #     # o -= np.mean(o)
                #     # Kling-gupta
                #     # L = spotpy.objectivefunctions.kge(o, s)
                #     # Log Gaussian likelihood without error
                #     L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s)
                #     #   *  1 / np.mean(o)
                #     # Log Gaussian with error set using std
                #     # LogLikehood as used in Vrugt et al. (2016)
                #     # Using standard deviation as indicator of error
                #     # o_err = 0.5*np.std(o) + 0.1*o
                #     # res = o - s
                #     # cor = np.corrcoef(res[:-1],res[1:])[1,0]
                #     # o_err = np.repeat(np.std(res)*np.sqrt((1+cor)/(1-cor)),
                #     #                   res.__len__())
                #     # L = spotpy.likelihoods.logLikelihood(o, s, measerror=o_err)
                # else:
                #     # Mean absolute error
                #     L = spotpy.objectivefunctions.mae(o, s) / mean(o)
                #     # Log Gaussian likelihood without error
                #     # L = spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(o, s) * \
                #     #     1 / np.mean(o)
                #     # Log Gaussian with error set using std
                #     # LogLikehood as used in Vrugt et al. (2016)
                #     # Using standard deviation as indicator of error...
                #     # o_err = np.std(o) + 0.25*o
                #     # res = o - s
                #     # cor = np.corrcoef(res[:-1],res[1:])[1,0]
                #     # o_err = np.repeat(np.std(res)*np.sqrt((1+cor)/(1-cor)),
                #     #                   res.__len__())
                #     # L = spotpy.likelihoods.logLikelihood(o, s, measerror=o_err)
                # # Normalize by data length
                # # L /= o.__len__()

            Ltot += [L]  # list of all likelihoods
            like += L  # "main" likelihood (used by algorithm)

        # Several datasets: list multi-objective function and individual ones
        # Only the first (multi-obj) one will be used by the algorithm
        if Obs.nobs > 1:
            like = [like] + Ltot

    # print(np.round(like, 2))

    return like
