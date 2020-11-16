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

import pandas as pd

from distutils.dir_util import mkpath
# ----------------------------------------------------------------------------


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
        params.sim_inputs(Config, Opti, Paras, Site, Config.PATH_SPA, it=it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print('|| running ECH2O...', end='\r')
        start = time.time()
        os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time()-start,
              'seconds (limit at '+Config.tlimit+')')

        # Store goodness of fit across outputs (or NA/blank if the run crashed)
        os.chdir(Config.PATH_EXEC)
        # Write goodness of fit (even if it's nan)
        outputs.store_GOF(Obs, Opti, Config, Site, it)
        # Write parameters values for this sequence
        #params.store(Opti, Config, it)

        if runOK(Obs, Opti, Config, mode='silent') == 0:
            # # If the run fails, let's give it one more chance!
            # os.chdir(Config.PATH_OUT)
            # os.system('rm -f '+Config.PATH_EXEC+'/*')
            # print('|| running ECH2O...', end='\r')
            # start = time.time()
            # os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
            #           Config.PATH_EXEC+'/ech2o.log')
            # print('run time:', time.time()-start,
            #       'seconds (limit at '+Config.tlimit+')')
            # os.chdir(Config.PATH_EXEC)
            # # Still not running properly? Report and move on
            # if runOK(Obs, Opti, Config, mode='verbose') == 0:
            f_failpar = Config.PATH_OUTmain+'/Parameters_fail.task'+Config.outnum+'.txt'
            if len(glob.glob(f_failpar)) == 0:
                with open(f_failpar, 'w') as f_in:
                    f_in.write('Sample,'+','.join(Opti.names)+'\n')
            with open(f_failpar, 'a') as f_in:
                f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')
            # If it is the very first iteration, record it for later use
            if it == 0:
                Opti.begfail = 1

            os.system('cp -p '+Config.PATH_EXEC+'/ech2o.log ' +
                      Config.PATH_OUT + '/ech2o.run'+Opti.itout+'.log')

        # # If it worked...
        # if runOK(Obs, Opti, Config, mode='silent') == 1:
        #     # Group sampling outputs
        #     outputs.store_sim(Obs, Opti, Config, Site, it)
        
        os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
                  Config.PATH_EXEC + '/ech2o.run'+Opti.itout+'.log')

        os.chdir(Config.PATH_OUT)
        # sys.exit()
        # Clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*.tab')
        os.system('rm -f '+Config.PATH_EXEC+'/Basin*.txt')
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
        os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                  Config.PATH_EXEC+'/ech2o.log')
        print('run time:', time.time()-start,
              'seconds (limit at', Config.tlimit, ')')
        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        # Group outputs
        outputs.store_sim(Obs, Opti, Config, Site, it)

        # if runOK(Obs, Opti, Config, mode='silent') == 0:
          # print('Warning: something came up during the run;' +
            #      ' some outputs might be truncated or missing')
            # Group outputs
            # outputs.store_sim(Obs, Opti, Config, Site, it)

        # Store log (in case)
        os.system('mv '+Config.PATH_EXEC+'/ech2o.log ' +
                  Config.PATH_OUT + '/ech2o_run'+str(it+1)+'.log')
        # Clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*')
# ===========================================================================


def morris_runs(Config, Opti, Obs, Paras, Site):

    # Simulations when varying the parameters, Morris's one-at-a-time
    # ---------------------------------------------------------------
    print('Total Number of iterations:', Opti.nruns)
    print('')
    print('-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    print('')

    Config.initobs = 0
    Opti.begfail = 0

    irun_tot = 0

    for itraj in range(Opti.nr):

        print('======================================')
        print('## Runs along trajectory #', itraj+1)
        print('--------------------------------------')

        # Array of parameters for this trajectory
        Opti.xpar = np.transpose(Opti.Bstar[:, :, itraj]) 

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
            os.system(Config.cmde_ech2o+' '+Config.cfg2_ech2o+' > ' +
                      Config.PATH_EXEC+'/ech2o.log')
            print('run time:', time.time() - start,
                  'seconds (limit at ', Config.tlimit, ')')

            # Check if it ran properly
            os.chdir(Config.PATH_EXEC)
            outputs.store_sim(Obs, Opti, Config, Site, irun_tot)

            # if runOK(Obs, Opti, Config) == 1:
                # Group outputs
                # outputs.store_sim(Obs, Opti, Config, Site, irun)
                # os.system('rm -f *.tab')

            if runOK(Obs, Opti, Config, mode='silent') == 0:
                # else:  # Not running properly? Report
                f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
                if len(glob.glob(f_failpar)) == 0:
                    with open(f_failpar, 'w') as f_in:
                        f_in.write('Trajectory/RadPoint,Sample,'+','.join(Opti.names)+'\n')
                with open(f_failpar, 'a') as f_in:
                    f_in.write(str(itraj+1)+','+str(irun+1)+','+','.join([str(x) for x in
                                                                          Opti.x])+'\n')
                # If it is the very first iteration, record it for later storage
                if irun == 0:
                    Opti.begfail = 1

                os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT +
                          '/ech2o_traj'+str(itraj+1)+'_run'+str(irun+1)+'.log')

            # Clean up
            os.system('rm -f '+Config.PATH_EXEC+'/*')

            irun_tot += 1

        # print(Obs.obs)
        # Only for debugging ------------------------------------------------------
        for oname in Obs.names:
            if Obs.obs[oname]['type'] != 'map':
                Obs.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'
        # ----------------------------------------------------------------------------

        # Calculate and output the elementary effects
        morris.ee(Config, Obs, Opti, itraj)

# ===========================================================================


def runOK(Obs, Opti, Config, mode='silent'):
    # -- Check if ECH2O ran properly

    isOK = 1

    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0 or \
       os.stat(Config.PATH_EXEC+'/BasinSummary.txt').st_size == 0:
        #if(mode == 'verbose'):
        print("Something went wrong, BasinSummary.txt is missing/empty...")
        isOK = 0

    else:
        # for oname in Obs.names:
        #     # print oname
        #     if Obs.obs[oname]['type'] != 'mapTs':
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
    # use the GOF files

    # File names for grouped simulations
    it = []
    for gof in Opti.GOFs:
        # Get last referenced index (i.e. run number) that worked (without NaN)
        tmp = pd.read_csv(Opti.GOFfiles[gof]).set_index('Sample')
        tmp2 = tmp.dropna()
        it += [tmp2.index[-1]+1]
    # in case this index is different, take minimum
    Config.itres = min(it)
    print('...but directly restarting from iter. ', Config.itres)

    # Rewrite the GOF files: to evenize between variable last it or 
    # remove the final NaN lines (but not those in between "good" runs)
    for gof in Opti.GOFs:
        print(gof)
        # Get last referenced index (i.e. run number) that worked (without NaN)
        tmp = pd.read_csv(Opti.GOFfiles[gof]).set_index('Sample').loc[1:Config.itres+1]
        tmp.to_csv(Opti.GOFfiles[gof], na_rep="nan")
            
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
