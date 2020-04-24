#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Python scripts for ECH2O
#
# -------
# Routine: Subroutines, miscellaneous
# -------
# Contributors: S. Kuppel, A.J. Neill
# Created on 10/2016
# -------------------------------------------------

import time, os, glob, sys, copy
import numpy as np

# ----------------------------------------------------------------------------
# -- Check if ECH2O ran properly

def runOK(Data, Opti, Config):

    isOK = 1

    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0 or os.stat(Config.PATH_EXEC+'/BasinSummary.txt').st_size == 0 :
        print("Something went wrong, BasinSummary.txt is missing/empty...")
        isOK = 0

    else:
        for oname in Data.names:
            # print oname
            if Data.obs[oname]['type']!='mapTs':
                if Data.obs[oname]['type']!='map':
                    f_test = Config.PATH_EXEC+'/'+Data.obs[oname]['sim_file']
                else:
                    f_test = Config.PATH_EXEC+'/'+Data.obs[oname]['sim_file']+'.map'
                
                if len(glob.glob(f_test)) == 0:
                    print("Something went wrong, no output for "+oname+" !!")
                    print('(i.e., '+f_test+' is missing...)')
                    isOK = 0    
                    # print 'Outputs are there...'
    
        # 2. Check it ran until the end
        f_test = Config.PATH_EXEC+'/BasinSummary.txt'
        try:
            tmp = np.genfromtxt(f_test,skip_header=1,unpack=True)[0]
        except ValueError:
            isOK = 0
        else:
            # print str(len(tmp2))
            # print str(Data.lsim)
            if type(tmp)==float or type(tmp)==np.float64 or type(tmp)==np.float32:
                isOK = 0
                print("Something went wrong, output of length 1 !")
            elif len(tmp)!=Data.lsim: # & len(tmp2)==Data.lsim:
                # print 'It ran until the end !'
                isOK = 0
                print("Something went wrong, output does not match the supposed sim length!")
                print('Output: '+str(len(tmp))+' , supposed to be: '+str(Data.lsim))

    return isOK
        
# ----------------------------------------------------------------------------
# -- Restart: trim outputs files to match the specified restart iteration 

def restart(Config, Opti, Data):

    # -- Get the last iteration that worked

    # File names for grouped simulations
    for oname in Data.names:
        Data.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'

    # Read one and get the last iteration written
    # (take one before, just to be sure writing was not ongoing there when it stopped)
    idx = np.genfromtxt(Data.obs[Data.names[0]]['sim_hist'],delimiter=',',skip_header=1,
                        unpack=True,usecols=(0))
    
    Config.itres = int(idx[::-1][0])      # Removed the "-1" since want to re-run whichever was the last run being written / completed
    # In case some iterations failed the max number of rows to read is smaller than Config.itres-1!!
    
    mxRow = len(idx)-1      # Still keep -1 here since want to discard the last run, since this is being re-run
    #print Config.itres-1
    #print mxRow

    # Some cleaning of the outputs
    for oname in Data.names:
        # -- Read grouped simulations
        tmp = np.genfromtxt(Data.obs[oname]['sim_hist'],delimiter=',',skip_header=1,
                            max_rows=mxRow)[::,1::]
        #print tmp.shape
        # -- Take out the iterations from/beyond restarting point
        #tmp = tmp[::,1::]
        #print Config.itres
        #print Config.trim
        #tmp = tmp[0:Config.itres-1,Config.trim+1::] ## TEMPORARY!!
        #print tmp.shape
        # -- Overwrite
        # Header
        with open(Data.obs[oname]['sim_hist'],'w') as f_out:
            f_out.write('Sample,'+','.join([str(i+1) for i in range(len(tmp[0]))])+'\n')
        # Iterations before starting point
        with open(Data.obs[oname]['sim_hist'],'a') as f_out:
            for i in range(mxRow):
                f_out.write(str(idx[i])+','+','.join([str(j) for j in list(tmp[i])])+'\n')

    # -- Some cleaning for parameters
    Config.f_par = Config.PATH_OUT+'/Parameters.txt'
    # Read grouped simulations
    tmp = np.genfromtxt(Config.f_par,delimiter=',',skip_header=1,max_rows=mxRow)[::,1::]
    # Take out the iteration from/beyond the restarting point
    ###tmp = tmp[0:Config.itres,1::]
    
    # Write
    with open(Config.f_par,'w') as f_out:
        f_out.write('Iteration,'+','.join(Opti.names)+'\n')
    with open(Config.f_par,'a') as f_out:
        for i in range(mxRow):
            f_out.write(str(idx[i])+','+','.join([str(x) for x in list(tmp[i])])+'\n')

