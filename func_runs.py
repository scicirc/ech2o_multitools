#!/usr/bin/python3
# -*- coding: utf-8 -*-
# *************************************************
#
# Multi-purpose routines ECH2O-iso:
# sensitivity analysis, calibration, and ensemble runs in general
#
# -------
# Routine: Runs routines
# -------
# Author: S. Kuppel
# Created on 04/2020
# -------------------------------------------------

import os
# import sys
import time
import glob
import numpy as np
import func_params as params
import func_outputs as outputs
import func_morris as morris
# ----------------------------------------------------------------------------


def calibMC_runs(Config, Opti, Obs, Paras, Site):

    # Calibration runs
    # -------------------

    print('Number of iterations      :', Opti.nit)
    if Config.restart == 1:
        restart(Config, Opti, Obs)
        print('...but directly restarting from iter. ', Config.itres)
    print('')
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print(' Entering the optimisation loop...')
    print('')

    # Output initialization
    Config.initpar = 0
    Config.initobs = 0
    Opti.begfail = 0

    if Config.restart == 1:
        it0 = Config.itres-1
    else:
        it0 = 0

    for it in range(it0, Opti.nit):

        Opti.itout = '%05i' % int(it+1)
        print('Iteration ', Opti.itout, ' of ', Opti.nit)

        # Create / clean up the run outputs directory
        if len(glob.glob(Config.PATH_EXEC)) == 0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        # Create the inputs for ECH2O
        params.sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...', end='\r')
        start = time.time()
        os.system(Config.cmde_ech2o+' config.ini > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time()-start,
              'seconds (limit at '+Config.tlimit+')')

        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        if runOK(Obs, Opti, Config) == 0:
            # If the run fails, let's give it one more chance!
            os.chdir(Config.PATH_OUT)
            os.system('rm -f '+Config.PATH_EXEC+'/*')
            print('|| running ECH2O...', end='\r')
            start = time.time()
            os.system(Config.cmde_ech2o+' config.ini > ' +
                      Config.PATH_EXEC+'/ech2o.log')
            print('run time:', time.time()-start,
                  'seconds (limit at '+Config.tlimit+')')
            os.chdir(Config.PATH_EXEC)
            # Still not running properly? Report and move on
            if runOK(Obs, Opti, Config) == 0:
                f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
                if len(glob.glob(f_failpar)) == 0:
                    with open(f_failpar, 'w') as f_in:
                        f_in.write('Sample,'+','.join(Opti.names)+'\n')
                        f_in.write(str(it+1)+','+','.join([str(x) for x in
                                                           Opti.x])+'\n')
                else:
                    with open(f_failpar, 'a') as f_in:
                        f_in.write(str(it+1)+','+','.join([str(x) for x in
                                                           Opti.x])+'\n')
                os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
                          Config.PATH_OUT + '/ech2o_'+Opti.itout+'.log')
                # If it is the very first iteration, record it for
                # main_outputs later
                if it == 0:
                    Opti.begfail = 1

        # If it worked...
        if runOK(Obs, Opti, Config) == 1:
            # Write parameters values for this sequence
            params.store(Opti, Config, it)
            # Group sampling outputs
            outputs.store_sim(Obs, Opti, Config, it)

        os.chdir(Config.PATH_OUT)
        # sys.exit()
        # Clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*')
# ============================================================================


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

        # Create the inputs for ECH2O
        params.sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...', end='\r')
        start = time.time()
        os.system(Config.cmde_ech2o+' config.ini > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time()-start,
              'seconds (limit at', Config.tlimit, ')')
        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        if runOK(Obs, Opti, Config) == 1:

            # Group outputs
            outputs.store_sim(Obs, Opti, Config, it)

            # Clean up
            os.system('rm -f '+Config.PATH_EXEC+'/*')
# ===========================================================================


def morris_runs(Config, Opti, Obs, Paras, Site):

    # Simulations when varying the parameters, Morris's one-at-a-time
    # ---------------------------------------------------------------
    if Config.mode == 'sensi_morris' and Config.MSinit == 0:

        print('Total Number of iterations : '+Opti.nruns)
    print('')
    print('-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print('')

    Config.initobs = 0
    Opti.begfail = 0

    print('')
    print('======================================')
    print('## Runs along trajectory #'+Config.numsim)
    print('--------------------------------------')

    # There are npara+1 runs for each trajectory
    for irun in range(Opti.nvar+1):

        # runnb = '%02i' % int(irun+1)
        print('Run '+str(irun+1)+' out of '+str('%02i' % int(Opti.nvar+1)))

        # Create / clean up the run outputs directory
        if len(glob.glob(Config.PATH_EXEC)) == 0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        # print
        # print('|- Creating parameter maps / table for this run...'

        # Create the inputs for ECH2O
        params.sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=irun)

        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...')
        start = time.time()
        os.system(Config.cmde_ech2o+' config.ini > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time() - start,
              'seconds (limit at ', Config.tlimit, ')')

        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        if runOK(Obs, Opti, Config) == 1:
            # Group outputs
            outputs.store_sim(Obs, Opti, Config, irun)
                # os.system('rm -f *.tab')
        # Not running properly? Report
        else:
            f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
            if len(glob.glob(f_failpar)) == 0:
                with open(f_failpar, 'w') as f_in:
                    f_in.write('Sample,'+','.join(Opti.names)+'\n')
                    f_in.write(str(irun+1)+','+','.join([str(x) for x in
                                                         Opti.x])+'\n')
            else:
                with open(f_failpar, 'a') as f_in:
                    f_in.write(str(irun+1)+','+','.join([str(x) for x in
                                                         Opti.x])+'\n')
            # If it is the very first iteration, record it for later storage
            if irun == 0:
                Opti.begfail = 1

            os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT +
                      '/ech2o_'+str(irun+1)+'.log')

        # Clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*')

    # print(Obs.obs)
    # Only for debugging ------------------------------------------------------
    for oname in Obs.names:
        if Obs.obs[oname]['type'] != 'map':
            Obs.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'
    # ----------------------------------------------------------------------------

    # Calculate and output the elementary effects
    morris.ee(Config, Obs, Opti)
# ===========================================================================


def runOK(Obs, Opti, Config):
    # -- Check if ECH2O ran properly

    isOK = 1

    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0 or \
       os.stat(Config.PATH_EXEC+'/BasinSummary.txt').st_size == 0:
        print("Something went wrong, BasinSummary.txt is missing/empty...")
        isOK = 0

    else:
        for oname in Obs.names:
            # print oname
            if Obs.obs[oname]['type'] != 'mapTs':
                if Obs.obs[oname]['type'] != 'map':
                    f_test = Config.PATH_EXEC+'/'+Obs.obs[oname]['sim_file']
                else:
                    f_test = Config.PATH_EXEC+'/' + \
                        Obs.obs[oname]['sim_file'] + '.map'

                if len(glob.glob(f_test)) == 0:
                    print("Something went wrong, no output for "+oname+" !!")
                    print('(i.e., '+f_test+' is missing...)')
                    isOK = 0
                    # print 'Outputs are there...'

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
                print("Something went wrong, output of length 1 !")
            elif len(tmp) != Obs.lsim:
                # & len(tmp2)==Obs.lsim:
                # print 'It ran until the end !'
                isOK = 0
                print("Something went wrong, output does not match the " +
                      "supposed sim length!")
                print('Output: '+str(len(tmp))+' , supposed to be: ' +
                      str(Obs.lsim))

    return isOK
# ==================================================================================


def restart(Config, Opti, Obs):
    # Restart: trim outputs files to match the specified restart iteration

    # -- Get the last iteration that worked

    # File names for grouped simulations
    for oname in Obs.names:
        Obs.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'

    # Read one and get the last iteration written (take one before,
    # just to be sure writing was not ongoing there when it stopped)
    idx = np.genfromtxt(Obs.obs[Obs.names[0]]['sim_hist'],
                        delimiter=',', skip_header=1, unpack=True, usecols=(0))

    Config.itres = int(idx[::-1][0])
    # Removed the "-1" since want to re-run whichever was the last run being
    # written / completed. In case some iterations failed the max number of
    # rows to read is smaller than Config.itres-1!!

    mxRow = len(idx)-1
    # Still keep -1 here since want to discard the last run, since this is
    # being re-run
    # print Config.itres-1
    # print mxRow

    # Some cleaning of the outputs
    for oname in Obs.names:
        # -- Read grouped simulations
        tmp = np.genfromtxt(Obs.obs[oname]['sim_hist'], delimiter=',',
                            skip_header=1, max_rows=mxRow)[::, 1::]
        # print tmp.shape
        # -- Take out the iterations from/beyond restarting point
        # tmp = tmp[::,1::]
        # print Config.itres
        # print Config.trim
        # tmp = tmp[0:Config.itres-1,Config.trim+1::] ## TEMPORARY!!
        # print tmp.shape
        # -- Overwrite
        # Header
        with open(Obs.obs[oname]['sim_hist'], 'w') as f_out:
            f_out.write('Sample,'+','.join([str(i+1) for i in
                                            range(len(tmp[0]))])+'\n')
        # Iterations before starting point
        with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
            for i in range(mxRow):
                f_out.write(str(idx[i])+','+','.join([str(j) for j in
                                                      list(tmp[i])]) + '\n')

    # -- Some cleaning for parameters
    Config.f_par = Config.PATH_OUT+'/Parameters.txt'
    # Read grouped simulations
    tmp = np.genfromtxt(Config.f_par, delimiter=',', skip_header=1,
                        max_rows=mxRow)[::, 1::]
    # Take out the iteration from/beyond the restarting point
    # #tmp = tmp[0:Config.itres,1::]

    # Write
    with open(Config.f_par, 'w') as f_out:
        f_out.write('Iteration,'+','.join(Opti.names)+'\n')
    with open(Config.f_par, 'a') as f_out:
        for i in range(mxRow):
            f_out.write(str(idx[i])+','+','.join([str(x) for x in
                                                  list(tmp[i])]) + '\n')
