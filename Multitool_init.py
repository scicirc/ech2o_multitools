#!/usr/bin/python3
# -*- coding: utf-8 -*-
# *************************************************
#
# Multi-purpose routines ECH2O-iso:
# sensitivity analysis, calibration, and ensemble runs in general
#
# -------
# Routine: Initialization routines
# -------
# Author: S. Kuppel
# Created on 04/2020
# -------------------------------------------------

import pcraster as pcr
import numpy as np
import os
import sys
import glob
import copy
from datetime import timedelta

import Multitool_sampling as sampling
import Multitool_params as params
import Multitool_morris as morris
import pandas as pd
import csv
# ----------------------------------------------------------------------------


def config(Config, options):

    # -- General config file (where parameters info is taken from)
    if options.file is not None:
        [file_py, ext] = os.path.splitext(options.file)
        if len(glob.glob(options.file)) == 0:
            sys.exit('# STOP. The file that define the assimilation ' +
                     'characteristics does not exist : \n   ...'+file_py+'.py')
        file = options.file
        frep = os.path.dirname(file)
        if frep == '':
            file = os.path.join(Config.PATH_MAIN, file)
        print('The user provided definition file is: '+options.file)

    # -- Number of CPUs used in parallel
    if options.ncpu is None:
        options.ncpu = 1
    Config.ncpu = copy.copy(options.ncpu)

    # -- Main mode
    if options.mode is None:
        sys.exit("Please choose which a script mode!")
    else:
        Config.mode = copy.copy(options.mode)

    # -- Restart: do not start from first iteration
    if options.restart is None:
        Config.restart = 0
    else:
        Config.restart = int(options.restart)
        if Config.restart > 1:
            sys.exit('Wrong value for restart')

    # -- Output path
    if options.outdir is None:
        if Config.mode != 'forward_runs':
            options.outdir = os.path.splitext(options.file)[0]
        else:
            tmp = options.cfg.split('_')[1].split('.')[0].split('-')
            if(len(tmp) == 2):
                [pref1, pref2] = tmp
                options.outdir = 'Res.ens'+options.nEns+'_'+pref2 + \
                    '.'+pref1+'.'+options.inEns
            if(len(tmp) == 1):
                options.outdir = 'Res.ens'+options.nEns+'.'+tmp[0]+'.' +\
                    options.inEns
            if(len(tmp) > 2):
                sys.exit('Error: incorrect config file format.')

    # -- scratch?
    if options.scratch is not None:
        if int(options.scratch) == 1:
            Config.scratch = 1
        else:
            Config.scratch = 0
    else:
        Config.scratch = 0

    # -- Time wall for ECH2O execution
    if options.tlimit is None:
        Config.tlimit = '100000'
    else:
        Config.tlimit = options.tlimit
    Config.tcmd = 'ulimit -t '+str(int(Config.tlimit)*int(Config.ncpu))+' ;'

    # -- Maps averaging
    if options.MapAv is not None:
        Config.MapAv = int(options.MapAv)
        if Config.MapAv == 1:
            if options.MapAvT in ['week', 'month', 'season']:
                Config.MapAvT = options.MapAvT
            else:
                sys.exit('Wrong maps averaging option!')
    else:
        Config.MapAv = 0

    # -- Resolution (default 50m, used with reference Spatial folder)
    if options.Resol is not None:
        Config.Resol = copy.copy(options.Resol)
    else:
        Config.Resol = '50m'

    # -- Report BasinSummary.txt
    if options.BSum is not None:
        Config.repBS = int(options.BSum)
    else:
        Config.repBS = 0

    # -- MS init: many things don't happen if it is the case
    if Config.mode == 'sensi_morris':
        if options.MSinit is None:
            Config.MSinit = 1
            sys.exit("Please state if you're initializating the MS sampling")
        else:
            Config.MSinit = int(options.MSinit)

        if options.MSspace is None:
            Config.MSspace = 'trajectory'
        else:
            Config.MSspace = copy.copy(options.MSspace)

    # -- Run ECH2O?
    Config.runECH2O = 1
    if Config.mode == 'calib_sampling' or (Config.mode == 'sensi_morris' and
                                           Config.MSinit == 1):
        Config.runECH2O = 0

    # print(options.outdir)

    print('')
    print('******************************************************************')
    if Config.mode.split('_')[0] == 'calib':
        print('CALIBRATION with EcH2O: ')
    elif Config.mode == 'forward_runs':
        print('ENSEMBLE RUNS with EcH2O')
    elif Config.mode == 'sensi_morris':
        print('MORRIS SENSITIVITY with EcH2O: ')
        print('- construction of the trajectories')
        print('- forward runs')
        print('- storage of outputs and info for posterior analysis : ' +
              'elementary effects, etc.')
        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
        print('')

    # -- Calibration: all parameter path
    if Config.mode.split('_')[0] == 'calib':
        Config.PATH_PAR = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                       'Parameters_samples'))
        Config.FILE_PAR = Config.PATH_PAR+'/'+options.outdir.split('.')[0] + \
            '_parameters.'
        # -- Creation of output directory
        if len(glob.glob(Config.PATH_PAR)) == 0:
            os.system('mkdir ' + Config.PATH_PAR)
        # -- Some verbose
        print('')
        print("Parameter samples' directory:          ", Config.PATH_PAR)
        print('')

    # -- Sensitivity: all parameter path
    if Config.mode == 'sensi_morris':

        print('')
        if(Config.MSspace == 'trajectory'):
            Config.PATH_TRAJ = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                            'Trajectories'))
            print("Trajectories directory:          ", Config.PATH_TRAJ)
        if(Config.MSspace == 'radial'):
            Config.PATH_TRAJ = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                            'RadialPoints'))
            print("Radial points directory:         ", Config.PATH_TRAJ)

        Config.FILE_TRAJ = Config.PATH_TRAJ+'/'+options.outdir.split('.')[0]
        # -- Creation of output directory
        if len(glob.glob(Config.PATH_TRAJ)) == 0:
            os.system('mkdir ' + Config.PATH_TRAJ)

        # -- Output of elementary effects
        if(Config.MSinit == 0):
            Config.PATH_EE = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                          'ElementaryEffects'))
            print("Elementary effects directory:          ", Config.PATH_EE)
            Config.FILE_EE = Config.PATH_EE+'/'+options.outdir.split('.')[0]

    print('')

    # -- For runs
    if Config.runECH2O == 1:
        # -- Execution command
        if options.exe is not None:
            if len(glob.glob(options.exe)) == 0:
                sys.exit('The user provided EXEC file was not found: ' +
                         options.exe)
            exe_ech2o = options.exe
            print('The user provided EXEC file is: '+options.exe)

        # Ech2O run config file
        if options.cfg is not None:
            cfg_ech2o = options.cfg+'.ini'
            if Config.mode == 1 and \
               len(glob.glob(Config.cfgdir+cfg_ech2o)) == 0:
                sys.exit('The user provided CFG file was not found: ' +
                         cfg_ech2o)
            if Config.mode == 2 and \
               len(glob.glob(Config.cfgdir+cfg_ech2o)) == 0:
                sys.exit('The user provided CFG file was not found: ' +
                         cfg_ech2o)
            print('The user provided CFG file is: '+cfg_ech2o)

        # Time limit and OMP use
        Config.cmde_ech2o = ' '.join([Config.tcmd, 'OMP_NUM_THREADS=' +
                                      str(Config.ncpu), './' + exe_ech2o,
                                      cfg_ech2o])

        # (if needed) get the parallel job number, based on the output dir name
        Config.numsim = options.outdir.split('.')[::-1][0]

        # Tracking age and/or tracers?
        if options.isTrck is not None:
            Config.isTrck = int(options.isTrck)
        else:
            Config.isTrck = 0

        if Config.isTrck == 1:
            if options.cfgTrck is not None:
                cfgTrck_ech2o = options.cfgTrck+'.ini'
            else:
                cfgTrck_ech2o = options.cfg.split('_')[0] + \
                    'Trck_'+options.cfg.split('_')[1]+'.ini'

            if len(glob.glob(Config.cfgdir+cfgTrck_ech2o)) == 0:
                sys.exit('The user provided CFGtrck file was not found: ' +
                         cfgTrck_ech2o)
            else:
                print('The user provided CFGtrck file is: '+cfgTrck_ech2o)

        # --- Defining the various PATHs
        if options.OMP_it is not None:
            Config.OMP_it = int(options.OMP_it)
            Config.PATH_OUTmain = \
                os.path.abspath(os.path.join(Config.PATH_MAIN,
                                             options.outdir))
            if len(glob.glob(Config.PATH_OUTmain)) == 0:
                os.system('mkdir ' + Config.PATH_OUTmain)
            Config.PATH_OUT = \
                os.path.abspath(os.path.join(Config.PATH_OUTmain,
                                             'EnsembleRun_' +
                                             str(Config.OMP_it)))
        else:
            Config.PATH_OUT = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                           options.outdir))

        Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,
                                                       'Spatial'))

        # -- Define execution directoy, depends on scratch options
        if Config.scratch > 0:
            if Config.scratch == 1:
                Config.PATH_EXEC = '/scratch/users/s08sk8/'+options.outdir
            if Config.scratch == 2:
                Config.PATH_EXEC = '/nobackup/users/s08sk8/'+options.outdir
        else:
            Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,
                                                            'tmp'))
        # -- Forward runs: parameter sets to use
        if Config.mode == 'forward_runs':

            if options.inEns is not None:
                Config.nEns = int(options.nEns)
                Config.FILE_PAR = Config.PATH_MAIN+'Input_Params/' +\
                    options.inEns+'.txt'
                # Config.FILE_PAR = Config.PATH_MAIN+'Input_Params/'+\
                # options.inEns+'.'+
                # options.nEns+'bestParams.txt'
                if len(glob.glob(Config.FILE_PAR)) == 0:
                    sys.exit('The param file (ensemble set) was not found: ' +
                             Config.FILE_PAR)

                print('')
                print('The ensemble param file is : '+Config.FILE_PAR)
            else:
                sys.exit('The size of ensemble simulations needs to be ' +
                         'specified (--nEns)')

    # -- Runs directories
    if Config.runECH2O == 1:
        # Forcings and reference maps
        Config.PATH_SPA_REF = \
            os.path.abspath(os.path.join(Config.PATH_MAIN,
                                         'Input_Maps_' + Config.Resol))
        Config.PATH_CLIM_REF = \
            os.path.abspath(os.path.join(Config.PATH_MAIN, 'Input_Climate'))

        # Output directory
        if len(glob.glob(Config.PATH_OUT)) == 0:
            os.system('mkdir ' + Config.PATH_OUT)

        # Creation of inputs directory (PATH_SPA will use parameter sampling)
        if len(glob.glob(Config.PATH_SPA)) == 0:
            os.system('mkdir ' + Config.PATH_SPA)
        # Copy of reference data
        os.system('cp -p '+Config.PATH_SPA_REF+'/*.map '+Config.PATH_SPA)
        os.system('cp -p '+Config.PATH_SPA_REF+'/SpeciesParams*.tab ' +
                  Config.PATH_SPA)

        # Execution directory
        if len(glob.glob(Config.PATH_EXEC)) == 0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        # -- Some verbose
        print('')
        print('Original/template maps & parameters:' + Config.PATH_SPA_REF)
        print('Original climate data:' + Config.PATH_CLIM_REF)
        if Config.scratch > 0:
            print('-----------------------------------------')
            print("Scratch storage activated! Hopefully that will speed " +
                  "things up...")
            print('Temporary maps & parameters:' + Config.PATH_SPA)
            # print 'Temporary climate data:', Config.PATH_CLIMS
            print('Temporary outputs:' + Config.PATH_EXEC)
            print('-----------------------------------------')
        else:
            print('Maps & parameters:' + Config.PATH_SPA)
        # print 'Climatic forcing: ', Config.PATH_CLIM
        print('Final outputs:          ' + Config.PATH_OUT)

        print('')

        # Symbolic link to executable
        if len(glob.glob(os.path.join(Config.PATH_OUT, exe_ech2o))) != 0:
            os.system('rm '+os.path.join(Config.PATH_OUT, exe_ech2o))
        os.symlink(os.path.join(Config.PATH_MAIN, exe_ech2o),
                   os.path.join(Config.PATH_OUT, exe_ech2o))

        # Copy config files (on /users even if there's scratch, for debugging)
        os.system('cp '+os.path.join(Config.PATH_MAIN, Config.cfgdir,
                                     cfg_ech2o) +
                  ' ' + os.path.join(Config.PATH_OUT, cfg_ech2o))
        if Config.isTrck == 1:
            os.system('cp -p '+os.path.join(Config.PATH_MAIN, Config.cfgdir,
                                            cfgTrck_ech2o) + ' ' +
                      os.path.join(Config.PATH_OUT, cfgTrck_ech2o))

        # -- Modify template config file with config-specific info
        with open(os.path.join(Config.PATH_OUT, cfg_ech2o), 'a') as fw:
            fw.write('\n\n\n#Simulation-specific folder section\n#\n\n')
            fw.write('Maps_Folder = '+Config.PATH_SPA+'\n')
            fw.write('Clim_Maps_Folder = '+Config.PATH_CLIM_REF+'\n')
            fw.write('Output_Folder = '+Config.PATH_EXEC+'\n')
            fw.write('ClimateZones = ClimZones_'+Config.Resol+'.map\n')
            fw.write('Isohyet_map = isohyet_'+Config.Resol+'.map\n')
            if Config.isTrck == 1:
                fw.write('Tracking = 1\n')
                fw.write('TrackingConfig = '+os.path.join(Config.PATH_OUT,
                                                          cfgTrck_ech2o)+'\n')
            else:
                fw.write('Tracking = 0')
            fw.write('\n\n')

        # -- Copy def file in the output
        os.system('cp '+file+' '+Config.PATH_OUT)

    # -- Import classes and setup from the def file
    sys.path.insert(0, Config.PATH_MAIN)
    Opti = __import__(file_py).Opti
    Data = __import__(file_py).Data
    Paras = __import__(file_py).Paras
    Site = __import__(file_py).Site

    # -- Trim: only saves the time steps within the trim
    if options.trimB is None:
        Config.trimB = 1
    else:
        Config.trimB = int(options.trimB)

    if options.trimBmap is None:
        Config.trimBmap = 1
    else:
        Config.trimBmap = int(options.trimBmap)

    if options.trimL is None:
        Config.trimL = Data.lsim - Config.trimB + 1
    else:
        if int(options.trimL) <= Data.lsim - Config.trimB+1:
            Config.trimL = int(options.trimL)
        else:
            sys.exit('Error: the specified output slicing start+length ' +
                     'goes beyond simulation time!')

    # Bare rock simulation
    if Opti.simRock is None:
        Opti.simRock = 0

    # -- Used for simulations
    Config.treal = [Data.simbeg+timedelta(days=x) for x in range(Config.trimL)]
    # sys.exit()

    # ---------------------------------------------------------------------------
    # Return the classes read in the definition file
    return (Opti, Data, Paras, Site)
# ==================================================================================


def parameters(Config, Opti, Paras, Site, options):

    Paras.names = sorted(Paras.ref.keys())
    Paras.n = len(Paras.names)
    # Read dictionary to get all param setup
    Opti.min = []
    Opti.max = []
    Opti.log = []
    Opti.names = []
    Opti.ind = []
    Opti.comp = []
    Paras.ind = {}
    ipar = 0
    ipar2 = 0
    Paras.isveg = 0

    for par in Paras.names:

        # Dimensions (soil or veg or 1)
        nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg'] *\
            (Site.nv-1)+1

        # Build vectors used in the optimisation
        if Config.mode != 'forward_runs':
            if type(Paras.ref[par]['min']) == float or \
               type(Paras.ref[par]['min']) == int:
                Opti.min = Opti.min + [Paras.ref[par]['min']]
                Opti.max = Opti.max + [Paras.ref[par]['max']]
            else:
                Opti.min = Opti.min + Paras.ref[par]['min']
                Opti.max = Opti.max + Paras.ref[par]['max']

            Opti.log = Opti.log + list(np.repeat(Paras.ref[par]['log'], nr))

        # Link betwen params and all variables
        Opti.ind = Opti.ind+list(np.repeat(ipar, nr))
        ipar += 1
        if nr > 1:
            Paras.ind[par] = list(np.arange(ipar2, ipar2+nr, 1))
        if nr == 1:
            Paras.ind[par] = ipar2
        ipar2 += nr
        # For outputs
        if Paras.ref[par]['soil'] == 1:
            Opti.names = Opti.names + [par + '_' + s for s in Site.soils]
            Opti.comp = Opti.comp + [i for i in range(Site.ns)]
        elif Paras.ref[par]['veg'] == 1:
            Opti.names = Opti.names + [par + '_' + s for s in Site.vegs]
            Opti.comp = Opti.comp + [i for i in range(Site.nv)]
            Paras.isveg += 1
        else:
            Opti.names = Opti.names + [par]
            Opti.comp = Opti.comp + [0]

    # Total number of variables
    Opti.nvar = len(Opti.names)

    # -- Calibration sampling: generate calibration parameters samples
    if Config.mode == 'calib_sampling':

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
        sampling.MC_sample(Opti, Config)
        # Done
        print('Parameters sample generation done.')

    # -- Calibration runs: retrieve previously generated samples
    if Config.mode == 'calib_runs':

        print('Get parameters samples for this job...')
        # Read parameters sample for this parallel run
        params.get(Opti, Config)

    # -- Forward ensemble runs: read directly the params from "best params"
    #    file
    if Config.mode == 'forward_runs':

        print('Get parameter set for these ensemble runs...')
        # tmp = list(np.genfromtxt(Config.FILE_PAR,delimiter=',',dtype= '|S30',
        # unpack=True)[0])
        tmp = list(pd.read_csv(Config.FILE_PAR, header=None).loc[:, 0])
        # print(Opti.names)
        # print(tmp)

    if tmp != Opti.names:
        sys.exit("The definition file and input parameter file ain't " +
                 "matching!")

    # pnames = ','.join([str(tmp[id.x]) for idx in range(Opti.nvar)])
    # print(tmp#pnames
    # print(Opti.names
    if(options.OMP_it is None):
        Opti.xpar = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                                  unpack=True)[1::]
    else:
        Opti.xpar = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                                  unpack=True)[1::][None, Config.OMP_it-1, :]

    # -- Sensitivity analysis: generate morris trajectories
    if Config.mode == 'sensi_morris':

        # Normalized step: plus-minus 0.5 of the normalized range of
        # each parameter
        Opti.stepN = np.zeros((Opti.nvar), np.float64) + 0.5

        if Config.MSinit == 1:
            morris.trajs(Config, Opti)
            print('Parameters trajectory generation done.')

        else:
            # Get the trajectory
            f_in = Config.PATH_TRAJ+'/'+options.outdir.split('.')[0] + \
                '.Bstar_traj' + Config.numsim+'.txt'
            # print(f_in
            Opti.xpar = np.genfromtxt(f_in, delimiter=',', skip_header=1)
            # print(Opti.xpar.shape
            # Reconstruct step (+- 0.5)
            if(Config.MSspace == 'trajectory'):
                Opti.dx = np.diff(Opti.xpar, axis=0)
            elif(Config.MSspace == 'radial'):
                Opti.dx = Opti.xpar[1::, :] - Opti.xpar[0, :]
            Opti.dx[Opti.dx != 0] = Opti.dx[Opti.dx != 0] / \
                np.abs(Opti.dx[Opti.dx != 0]) * 0.5
            if(np.ptp(Opti.dx) != 1.0 or np.min(Opti.dx) != -0.5 or
               np.max(Opti.dx) != 0.5):
                sys.exit('Error: The fetched Bnorm has a problem...')

        # Total number of runs
        Opti.nruns = (Opti.nvar+1) * Opti.nr

    # Calibration or SA runs: Check that the same parameters are used
    if Config.mode == 'calib_runs' or \
       (Config.mode == 'sensi_morris' and Config.MSinit == 0):
        if(len(Opti.xpar[0]) != len(Opti.names)):
            sys.exit("The definition file and input parameter file ain't " +
                     "matching!")
# ==================================================================================


def runs(Config, Opti, Data, Paras, Site):

    # Get the observations used for optimisation
    Data.names = sorted(Data.obs.keys())

    # Remove the default map files of calibrated param in the inputs directory
    # --> helps checking early on if there is an improper map update
    for pname in Paras.names:
        if Paras.ref[pname]['veg'] == 0:
            os.system('rm -f '+Config.PATH_SPA+'/' +
                      Paras.ref[pname]['file']+'.map')

    # Soils / units maps
    Config.cloneMap = pcr.boolean(pcr.readmap(Config.PATH_SPA+'/base.map'))
    pcr.setclone(Config.PATH_SPA+'/base.map')

    Site.bmaps = {}
    for im in range(Site.ns):
        Site.bmaps[Site.soils[im]] = pcr.readmap(Config.PATH_SPA+'/' +
                                                 Site.sfiles[im])
    Site.bmaps['unit'] = pcr.readmap(Config.PATH_SPA+'/unit.map')
    if Opti.simRock == 1:
        # Site.bmaps['nolowK'] = readmap(Config.PATH_SPA+'/unit.nolowK.map')
        Site.bmaps['rock'] = pcr.readmap(Config.PATH_SPA+'/unit.rock.map')

    Site.bmaps['chanmask'] = pcr.readmap(Config.PATH_SPA+'/chanmask.map')
    Site.bmaps['chanmask_NaN'] = pcr.readmap(Config.PATH_SPA +
                                             '/chanmask_NaN.map')

    # Reference dictionary for vegetation inputs file
    Opti.vref = {}
    with open(Config.PATH_SPA_REF+'/' + Site.vfile, 'r') as csvfile:
        paramread = list(csv.reader(csvfile, delimiter='\t'))
    exit
    # "Head": number of species and of params
    Opti.vref['header'] = paramread[0][0:len(paramread[0])]
    # All parameters values (keep strings!)
    for iv in range(Site.nv):
        Opti.vref[iv] = paramread[iv+1][0:len(paramread[iv+1])]
    # "Footers" : name of head1 and of parameters
    Opti.vref['footer'] = paramread[Site.nv+1][0:len(paramread[Site.nv+1])]
    Opti.vref['name'] = paramread[Site.nv+2][0:len(paramread[Site.nv+2])]

    # if Config.isTrck == 1:
    #     os.system('cp '+Config.PATH_SPA+'/dD_snowpack.map '+Config.PATH_SPA+
    # '/dD.snowpack.map')
    #     os.system('cp '+Config.PATH_SPA+'/dD_surface.map '+Config.PATH_SPA+
    # '/dD.surface.map')
    #     os.system('cp '+Config.PATH_SPA+'/dD_soil1.map '+Config.PATH_SPA+
    # '/dD.L1.map')
    #     os.system('cp '+Config.PATH_SPA+'/dD_soil2.map '+Config.PATH_SPA+
    # '/dD.L2.map')
    #     os.system('cp '+Config.PATH_SPA+'/dD_soil3.map '+Config.PATH_SPA+
    # '/dD.L3.map')
    #     os.system('cp '+Config.PATH_SPA+'/dD_groundwater.map '+
    # Config.PATH_SPA+'/dD.GW.map')

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
