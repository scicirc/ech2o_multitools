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
from datetime import datetime


import func_sampling as sampling
import func_params as params
import func_morris as morris
import pandas as pd
import csv

# import multiprocessing as mp


from distutils.dir_util import mkpath, copy_tree, remove_tree
from distutils.file_util import copy_file

# ----------------------------------------------------------------------------


def config(options):


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
    Data = __import__(file_py).Data
    Paras = __import__(file_py).Paras
    Site = __import__(file_py).Site

    # Store name of def file in Config
    Config.file = copy.copy(options.file)
    frep = os.path.dirname(Config.file)
    if frep == '':
        Config.file = os.path.join(cwd_tmp, Config.file)

    # -- Main mode
    if Config.mode not in ['calib_MCsampling', 'calib_MCruns',
                           'calib_DREAM', 'forward_runs',
                           'sensi_morris']:
        sys.exit("Please choose a valid script mode!")

    # Working directory specified in the def file
    if Config.PATH_MAIN is None:
        Config.PATH_MAIN = cwd_tmp

    # -- Output directory
    # Base name
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
    # Absolute location
    if Config.mode == 'forward_runs' and options.OMP_it is not None:
        Config.OMP_it = int(options.OMP_it)
        Config.PATH_OUTmain = \
            os.path.abspath(os.path.join(Config.PATH_MAIN,
                                         options.outdir))
        if len(glob.glob(Config.PATH_OUTmain)) == 0:
            mkpath(Config.PATH_OUTmain)
        Config.PATH_OUT = \
            os.path.abspath(os.path.join(Config.PATH_OUTmain,
                                         'EnsembleRun_' +
                                         str(Config.OMP_it)))
    else:
        Config.PATH_OUT = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                       options.outdir))

    # -- Restart: do not start from first iteration
    if options.restart is None:
        Config.restart = 0
    else:
        Config.restart = int(options.restart)
        if Config.restart > 1:
            sys.exit('Wrong value for restart')

    # -- Resolution (default 50m, used with reference Spatial folder)
    if Site.Resol is None:
        Site.Resol = '50m'

    # -- MS init: many things don't happen if it is the case
    if Config.mode == 'sensi_morris':
        # if options.MSinit is None:
        #     Config.MSinit = 1
        #     sys.exit("Please state if you're initializating the MS sampling")
        # else:
        #     Config.MSinit = int(options.MSinit)

        if options.MSspace is None:
            Config.MSspace = 'trajectory'
        else:
            Config.MSspace = copy.copy(options.MSspace)

    # -- Run ECH2O?
    Config.runECH2O = 1
    if Config.mode == 'calib_MCsampling':  # or
        # (Config.mode == 'sensi_morris' and Config.MSinit == 1):
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
    if Config.mode == 'calib_DREAM':
        if Opti.DREAMpar is None:  # By default, sequential runs
            Opti.DREAMpar = 'seq'
        elif Opti.DREAMpar in ['mpc', 'mpi']:
            Opti.parallel = True
            #     from mpi4py import MPI
            #     Config.MPIcomm = MPI.COMM_WORLD
            #     Config.MPIsize = comm.Get_size()
            #     Config.MPIrank = comm.Get_rank()
            # Determine with process we're on, to avoid multiple prints
            
            # If in parallel mode, use as many CPUs as possible unless
            # the number of chains is larger
            # NOPE; not fit for cluster where cpu_count() doesn't work?
            # just choose cpus wisely in def + job files
            # Config.ncpu = max(1, mp.cpu_count() // Opti.nChains)

    # 
    if (Opti.DREAMpar == 'mpi' and options.mpi != 1) or \
       (Opti.DREAMpar != 'mpi' and options.mpi == 1):
        sys.exit('Inconsitent MPI flags between Opti and options')

    # print(options.outdir)
    # -- Calibration: all parameter path (and datasets, if needed)
    if Config.mode in ['calib_MCsampling', 'calib_MCruns']:
        Config.PATH_PAR = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                       'Calibration_Samples'))
        Config.FILE_PAR = Config.PATH_PAR+'/'+options.outdir.split('.')[0] + \
            '_parameters.'
        # -- Creation of output directory
        if len(glob.glob(Config.PATH_PAR)) == 0:
            mkdir(Config.PATH_PAR)
        # -- Some verbose
        print('')
        print("Parameter samples' directory:          ", Config.PATH_PAR)
        print('')

    if Config.mode.split('_')[0] == 'calib':
        # Observations (for now only needed in DREAM mode)
        Config.PATH_OBS = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                       'Calibration_Datasets'))

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
            mkdir(Config.PATH_TRAJ)

        # -- Output of elementary effects
        # if(Config.MSinit == 0):
        Config.PATH_EE = os.path.abspath(os.path.join(Config.PATH_MAIN,
                                                      'ElementaryEffects'))
        print("Elementary effects directory:          ", Config.PATH_EE)
        Config.FILE_EE = Config.PATH_EE+'/'+options.outdir.split('.')[0]

    # print('')

    # ---------------------------------------------------------------------------
    # Return the classes read in the definition file
    return (Config, Opti, Data, Paras, Site)
# ==================================================================================


def parameters(Config, Opti, Paras, Site, options):

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
    Opti.comp = []
    Paras.ind = {}
    ipar = 0
    ipar2 = 0
    Paras.isveg = 0

    # Initial Sampling type (for DREAM)
    if not hasattr(Opti, 'initSample'):
        Opti.initSample = 'uniform'

    for par in Paras.names:

        # Dimensions (soil or veg or 1)
        nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg'] *\
            (Site.nv-1)+1

        # Build vectors used in the optimisation
        if Config.mode != 'forward_runs':

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
                    sys.exit('Wrong "min" or "max" format for', par)
            else:
                Opti.min += list(np.repeat(np.nan, nr))
                Opti.max += list(np.repeat(np.nan, nr))
            # In log sampling case with spotpy, log-transform boundaries too
            if Config.mode == 'calib_DREAM' and Opti.log[ipar2] == 1:
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
                for j in range(nr):
                    # (Log sampling case taken into account above)
                    Opti.std += [sdmult*(Opti.max[ipar2+j]-Opti.min[ipar2+j])]

        # Link betwen params and all variables
        Opti.ind += list(np.repeat(ipar, nr))
        ipar += 1
        if nr > 1:
            Paras.ind[par] = list(np.arange(ipar2, ipar2+nr, 1))
        if nr == 1:
            Paras.ind[par] = ipar2
        ipar2 += nr
        # For outputs
        if 'soil' not in Paras.ref[par].keys():
            Paras.ref[par]['soil'] = 0
        if 'veg' not in Paras.ref[par].keys():
            Paras.ref[par]['veg'] = 0

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
        sampling.MC_sample(Opti, Config)
        # Done
        print('Parameters sample generation done.')

    # -- Calibration runs: retrieve previously generated samples
    if Config.mode == 'calib_MCruns':
        print('Get parameters samples for this job...')
        # Read parameters sample for this instance of array task
        params.get(Opti, Config)

    # -- Forward ensemble runs: read directly the params from "best params"
    #    file
    if Config.mode == 'forward_runs':

        print('Get parameter set for these ensemble runs...')
        # tmp = list(np.genfromtxt(Config.FILE_PAR,delimiter=',',dtype= '|S30',
        # unpack=True)[0])
        tmp = list(pd.read_csv(Config.FILE_PAR, header=None).loc[:, 0])

        if tmp != Opti.names:
            print(Opti.names)
            print(tmp)
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
                                      unpack=True)[1::][None,
                                                        Config.OMP_it-1, :]

    # -- Sensitivity analysis: generate morris trajectories
    if Config.mode == 'sensi_morris':

        # Normalized step: plus-minus 0.5 of the normalized range of
        # each parameter
        Opti.stepN = np.zeros((Opti.nvar), np.float64) + 0.5

        # if Config.MSinit == 1:
        morris.trajs(Config, Opti)
        print('Parameters trajectory generation done.')

        # else:
        #     # Get the trajectory
        #     f_in = Config.PATH_TRAJ+'/'+options.outdir.split('.')[0] + \
        #         '.Bstar_traj' + Config.numsim+'.txt'
        #     # print(f_in
        #     Opti.xpar = np.genfromtxt(f_in, delimiter=',', skip_header=1)
        #     # print(Opti.xpar.shape
        #     # Reconstruct step (+- 0.5)
        #     if(Config.MSspace == 'trajectory'):
        #         Opti.dx = np.diff(Opti.xpar, axis=0)
        #     elif(Config.MSspace == 'radial'):
        #         Opti.dx = Opti.xpar[1::, :] - Opti.xpar[0, :]
        #     Opti.dx[Opti.dx != 0] = Opti.dx[Opti.dx != 0] / \
        #         np.abs(Opti.dx[Opti.dx != 0]) * 0.5
        #     if(np.ptp(Opti.dx) != 1.0 or np.min(Opti.dx) != -0.5 or
        #        np.max(Opti.dx) != 0.5):
        #         sys.exit('Error: The fetched Bnorm has a problem...')
        Opti.xpar = Opti.Bstar
        # Reconstruct step (+- 0.5)
        if(Config.MSspace == 'trajectory'):
            Opti.dx = np.diff(Opti.xpar, axis=0)
        elif(Config.MSspace == 'radial'):
            Opti.dx = Opti.xpar[1::, :, :] - Opti.xpar[0, :, :]
        Opti.dx[Opti.dx != 0] = Opti.dx[Opti.dx != 0] / \
            np.abs(Opti.dx[Opti.dx != 0]) * 0.5
        if(np.ptp(Opti.dx) != 1.0 or np.min(Opti.dx) != -0.5 or
           np.max(Opti.dx) != 0.5):
            sys.exit('Error: Bnorm has a problem...')

        # Total number of runs
        Opti.nruns = (Opti.nvar+1) * Opti.nr

    # Calibration or SA runs: Check that the same parameters are used
    if Config.mode == 'calib_MCruns' or Config.mode == 'sensi_morris':
        # (Config.mode == 'sensi_morris' and Config.MSinit == 0):
        if(len(Opti.xpar[0]) != len(Opti.names)):
            sys.exit("The definition file and input parameter file ain't " +
                     "matching!")
# ==================================================================================


def runs(Config, Opti, Data, Paras, Site, options):

    # -- Runs directories
    # Forcings and reference maps
    Config.PATH_SPA_REF = \
        os.path.abspath(os.path.join(Config.PATH_MAIN,
                                     'Input_Maps_' + Site.Resol))
    Config.PATH_CLIM = \
        os.path.abspath(os.path.join(Config.PATH_MAIN, 'Input_Climate'))

    # Creation of inputs directory (PATH_SPA will use parameter sampling)
    # (in case of parallel runs as in DREAM+mpc/mpi, it will serve as
    # -yet another- template)
    # print(os.environ.keys())
    # if Opti.parallel == False:
    Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,
                                                   'Spatial'))
    # else:
    #     # In parallel computing case, one PATH_SPA per task
    #     if Opti.DREAMpar  == 'mpi':
    #         # Running n parallel on a unix system (open MPI type).
    #         # Check the ID of the current mpi task
    #         if 'OMPI_COMM_WORLD_RANK' in os.environ.keys():
    #             call = str(int(os.environ['OMPI_COMM_WORLD_RANK'])+1)
    #         elif 'PMI_RANK' in os.environ.keys():
    #             call = str(int(os.environ['PMI_RANK'])+1)
    #         else:
    #             sys.exit('The ID of this task could not be found...')
    #     if Opti.DREAMpar == 'mpc':
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
    if Config.exe is not None:
        if len(glob.glob(Config.exe)) == 0:
            sys.exit('The user provided EXEC file was not found: ' +
                     Config.exe)
            print('The user provided EXEC file is: '+Config.exe)
    
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
        Config.tlimit = '2000'
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
    if Config.scratch > 0:
        if Config.scratch == 1:
            Config.PATH_EXEC = '/scratch/sylvain.kuppel/'+options.outdir
        # if Config.scratch == 2:
        #     Config.PATH_EXEC = '/nobackup/users/s08sk8/'+options.outdir
        print('-----------------------------------------')
        print("Scratch storage activated! Hopefully that will speed " +
              "things up...")
    else:
        Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT, 'tmp'))

    # # In parellel computing case, one execution folder per task
    # if Opti.parallel == True:
    #     # In parallel computing case, one PATH_SPA per task
    #     if Opti.DREAMpar  == 'mpi':
    #         # Running n parallel on a unix system (open MPI type).
    #         # Check the ID of the current mpi task
    #         if 'OMPI_COMM_WORLD_RANK' in os.environ.keys():
    #             call = str(int(os.environ['OMPI_COMM_WORLD_RANK'])+1)
    #         elif 'PMI_RANK' in os.environ.keys():
    #             call = str(int(os.environ['PMI_RANK'])+1)
    #         else:
    #             sys.exit('The ID of this task could not be found...')
    #     elif Opti.DREAMpar == 'mpc':
    #         # Running n parallel on a single (Windows) computer.
    #         # ID of the current computer core
    #         call = str(os.getpid())

    #     Config.PATH_EXEC += '_' + call
    #     print('path exec 1:', Config.PATH_EXEC)

    # -- Ech2O run config file
    # Path for EcH2O config file
    Config.cfgdir = Config.PATH_MAIN+'Input_Configs/'
    if options.cfg is not None:
        Config.cfg_ech2o = options.cfg+'.ini'
        if len(glob.glob(Config.cfgdir+Config.cfg_ech2o)) == 0:
            sys.exit('The user provided CFG file was not found: ' +
                     Config.cfg_ech2o)
    else:
        sys.exit('Error: the script need a template ech2o config file!')

    # (if needed) get the parallel job number, based on the output dir name
    Config.numsim = options.outdir.split('.')[::-1][0]

    # -- Tracking age and/or tracers?
    if Site.isTrck is None:
        Site.isTrck = 0
    if Site.isTrck == 1:
        if options.cfgTrck is not None:
            cfgTrck_ech2o = options.cfgTrck+'.ini'
        else:
            cfgTrck_ech2o = options.cfg.split('_')[0] + \
                'Trck_'+options.cfg.split('_')[1]+'.ini'

        if len(glob.glob(Config.cfgdir+cfgTrck_ech2o)) == 0:
            sys.exit('The user provided CFGtrck file was not found: ' +
                     cfgTrck_ech2o)

    # -- Forward runs: parameter sets to use
    if Config.mode == 'forward_runs':
        if options.inEns is not None:
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
            sys.exit('The param file (ensemble set) needs to be ' +
                     'specified (--inEns)')
        # Number of runs is the size of the parameters file, unless specified
        nv = np.genfromtxt(Config.FILE_PAR, delimiter=',',
                           unpack=True)[1::].shape[0]
        if options.nEns is not None:
            Config.nEns = min(int(options.nEns), nv)
        else:
            Config.nEns = nv

    # --- Reporting stuff
    # Trim: only saves the time steps within the trim
    if not hasattr(Config, 'trimB'):
        Config.trimB = 1

    # Initial cutoff for map reporting
    if not hasattr(Config, 'trimBmap'):
        Config.trimBmap = 1

    # Length of saved outputs
    if not hasattr(Config, 'trimL'):
        Config.trimL = Data.lsim - Config.trimB + 1
    else:
        if Config.trimL > Data.lsim - Config.trimB+1:
            sys.exit('Error: the specified output slicing start+length ' +
                     'goes beyond simulation time!')

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

    # Report BasinSummary.txt
    if hasattr(Config, 'repBS'):
        Config.repBS = 0

# ==================================================================================
# -- Read measurements for calibration


def observations(Config, Opti, Data):

    # Set the observations types collected from runs (sim outputs)
    # (and compared to measurements if there is calibration)
    Data.names = sorted(Data.obs.keys())

    # Date of the simulations (used for calibration periods, mostly)
    Data.lsimEff = Data.lsim - Data.lspin
    Data.simt = [Data.simbeg + timedelta(days=x)
                 for x in range(Data.lsimEff)]
    # Same thing in froward mode (needs to be merged...)
    Config.treal = [Data.simbeg+timedelta(days=x) for x in range(Config.trimL)]

    if Config.mode == 'calib_DREAM':

        # print('Reading measured datasets for calibration...')

        Opti.obs = {}  # np.full((Data.nobs, Config.trimL), np.nan)

        for oname in Data.names:

            # print(oname)
            # -- Get the obs
            f_obs = Config.PATH_OBS + '/' + Data.obs[oname]['obs_file']

            # Read the file, keeping only the data and relevant TS columns
            tmp = pd.read_csv(f_obs, sep=';').iloc[
                :, [0, Data.obs[oname]['obs_col']-1]]
            tmp.columns.values[:] = ['Date', 'value']
            # Convert date column to datetime
            tmp['Date'] = pd.to_datetime(tmp['Date'], format='%Y-%m-%d')
            # Convert date column to datetime
            tmp['value'] = tmp['value'] * Data.obs[oname]['obs_conv']

            # Calibration period:
            # Check if specified, otherwise use the whole simulation
            # in any case, remove the spinup (it will be removed from
            # simulation outputs in post-processing)
            if 'fit_beg' not in Data.obs[oname].keys() or \
               type(Data.obs[oname]['fit_beg']) is not datetime.date:
                fitbeg = Data.simt[0]
            else:
                fitbeg = max(Data.obs[oname]['fit_beg'], Data.simt[0])
            if 'fit_end' not in Data.obs[oname].keys() or \
               type(Data.obs[oname]['fit_end']) is not datetime.date:
                fitend = Data.simt[Data.lsimEff-1]
            else:
                fitend = min(Data.obs[oname]['fit_end'],
                             Data.simt[Data.lsimEff - 1])

            # Crop obs between desired time frame
            tmp = tmp.loc[(tmp.Date >= fitbeg) & (tmp.Date <= fitend)]

            Opti.obs[oname] = tmp.dropna(how='all')
# ==========================================================================


def files(Config, Opti, Paras, Site):

    # Copying and editing directories and files.

    # Main output directory: create if needed
    mkpath(Config.PATH_OUT)
    # copy definition file there
    copy_file(Config.file, Config.PATH_OUT)

    # Copy of reference input parameters
    # remove_tree(Config.PATH_SPA)
    mkpath(Config.PATH_SPA)
    os.system('cp -f '+Config.PATH_SPA_REF+'/*.map '+ Config.PATH_SPA)
    copy_file(Config.PATH_SPA_REF+'/SpeciesParams.tab', Config.PATH_SPA)
    
    # EcH2O executable file: clean up / update old symlink
    copy_file(os.path.join(Config.PATH_MAIN, Config.exe),
              Config.PATH_OUT)

    # EcH2O config file: copy from template and edit 
    copy_file(os.path.join(Config.PATH_MAIN, Config.cfgdir, Config.cfg_ech2o),
              os.path.join(Config.PATH_OUT, 'config.ini'))
    with open(os.path.join(Config.PATH_OUT, 'config.ini'), 'a') as fw:
        fw.write('\n\n\n#Simulation-specific folder section\n#\n\n')
        fw.write('Clim_Maps_Folder = '+Config.PATH_CLIM+'\n')
        fw.write('ClimateZones = ClimZones_'+Site.Resol+'.map\n')
        fw.write('Isohyet_map = isohyet_'+Site.Resol+'.map\n')
        # Input maps directory defined here (unless you need parallel runs
        # with DREAM, in which case the directory will change on the fly,
        # see spot_setup.simulation)
        if Opti.parallel is False:
            fw.write('Maps_Folder = '+Config.PATH_SPA+'\n')
            fw.write('Output_Folder = '+Config.PATH_EXEC+'\n')
        # Further edit regarding tracking
        if Site.isTrck == 1 :
            fw.write('Tracking = 1\n')
            fw.write('TrackingConfig = configTrck.ini\n')
        else:
            fw.write('Tracking = 0\n')

    # If tracking, copy template configTrck file for EcH2O
    if Site.isTrck == 1:
        copy_file(os.path.join(Config.PATH_MAIN, Config.cfgdir,cfgTrck_ech2o),
                  os.path.join(Config.PATH_OUT, 'configTrck.ini'))


    # -- Preparing inputs maps/files for site geometry etc.

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
    # Stream network
    Site.bmaps['chanmask'] = pcr.readmap(Config.PATH_SPA+'/chanmask.map')
    Site.bmaps['chanmask_NaN'] = pcr.readmap(Config.PATH_SPA +
                                             '/chanmask_NaN.map')
    # Bare rock simulation
    if Site.simRock is not None and Site.simRock == 1:
        # Site.bmaps['nolowK'] = readmap(Config.PATH_SPA+'/unit.nolowK.map')
        Site.bmaps['rock'] = pcr.readmap(Config.PATH_SPA+'/unit.rock.map')

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
