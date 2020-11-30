#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Python scripts for ECH2O
#
# -------
# Routine: Subroutines for outputs management
# -------
# Contributors: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import scipy.io as spio
import pcraster as pcr
import os
import glob
import sys
import numpy as np
from datetime import datetime
import pandas as pd

import func_runs as runs
import func_GOFs as GOFs

from distutils.dir_util import mkpath
from distutils.file_util import copy_file

# ----------------------------------------------------------------------------
# -- Read a given simulation output (time series only)


def read_sim(Config, Obs, oname, it=0):

    #print(oname)
    # Point-scale time series -------------------------------------
    if Obs.obs[oname]['type'] == 'Ts':
        #print(Obs.obs[oname]['sim_file'])
        # HEader in EcH2O files
        hskip = Obs.nts+3
        idx = np.argsort(np.array(Obs.sim_order))[Obs.obs[oname]['sim_pts']-1]
        # Read
        sim = pd.read_table(Obs.obs[oname]['sim_file'],error_bad_lines=False,
                            skiprows=hskip, header=None).set_index(0)
        # Check if it ran properly (sometimes some time steps are skipped in *tab...)
        # If steps are missing but the series is long enough, let's replace with nan
        if len(sim) < Obs.lsim and len(sim) > 365:
            idx2 = np.arange(1,Obs.lsim+1)
            sim = sim.reindex(idx2)
            print("Warning: some steps were missing in", Obs.obs[oname]['sim_file'],
                  '(',','.join([str(x) for x in list(pd.isnull(sim[1]).nonzero()[0]+1)]))
            copy_file(Obs.obs[oname]['sim_file'],
                      Config.PATH_OUT+'/'+Obs.obs[oname]['sim_file']+ 
                      '.run'+str(it+1)+'.txt')

    # Integrated variables (in BasinSummary.txt) -----------------
    elif Obs.obs[oname]['type'] == 'Total':
        idx = Obs.obs[oname]['sim_pts']-1
        # Read
        sim = pd.read_table(Obs.obs[oname]['sim_file'], error_bad_lines=False)
        sim = sim.set_axis([str(i) for i in np.arange(1,sim.shape[0]+1)])


    # Get observation column
    sim = sim.iloc[:,idx] * Obs.obs[oname]['sim_conv']
    #print(len(sim),'(2)')
    # Trim (spinup, transient state, etc.)
    sim = sim.loc[Obs.saveB:Obs.saveB+Obs.saveL-1]
    #print(len(sim),'(3', Obs.saveB, Obs.saveB+Obs.saveL-1,')')

    if len(sim) != Obs.saveL:
        print("Warning: problem with "+oname+" trim: read length is " + 
              str(len(sim)) + ' instead of '+str(Obs.saveL))
        sim = [np.nan] * Obs.saveL

    # if Obs.saveB > 1 or Obs.saveL < Obs.lsim:
    #     sim = sim[Obs.saveB-1:Obs.saveB-1+Obs.saveL]
    #     # print(len(sim))
    # #print(sim)
    
    return list(sim)
# ----------------------------------------------------------------------------
# -- Store in files for later use


def store_sim(Obs, Opti, Config, Site, it):

    mode = 'verbose'
    
    # -- Group the output files in one across simulations,
    #    separating by observations points and veg type where it applies
    for oname in Obs.names:
        if Obs.obs[oname]['type'] != 'map' and \
           Obs.obs[oname]['type'] != 'mapTs' and \
                                     (it == 0 or \
                                      (Opti.begfail == 1 and Config.mode!='sensi_morris')):
            # Historic time series file names
            Obs.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'
            # Header of files
            with open(Obs.obs[oname]['sim_hist'], 'w') as f_out:
                f_out.write('Sample,'+','.join([str(i+1) for i in
                                                range(Obs.saveL)])+'\n')

    # Reinit begfail (otherwise will never write all!)
    Opti.begfail = 0

    # Save current run outputs (and delete the file to relieve Maxwell...)
    for oname in Obs.names:

        #print(oname)
        firstMapTs = 1

        if(oname != Obs.names[0]):
            mode = 'silent'

        # Time series ---------------------------------------------------------
        # Integrated variables (in BasinSummary.txt)

        if Obs.obs[oname]['type'] == 'Total':

            if runs.runOK(Obs, Opti, Config, mode) == 1:

                idx = Obs.obs[oname]['sim_pts']-1
                # Read
                sim = pd.read_table(Obs.obs[oname]['sim_file'], error_bad_lines=False)
                # Basin*Summary.txt files don't have an index
                sim = sim.set_axis([str(i) for i in np.arange(1,sim.shape[0]+1)])
                # Get observation column
                sim = sim.iloc[:, idx] * Obs.obs[oname]['sim_conv']
                # Trim (spinup, transient state, etc.)
                # Index starts at 0
                if Obs.obs[oname]['sim_file'] == 'BasinSummary.txt':
                    sim = sim[Obs.saveB-1:Obs.saveB+Obs.saveL-1]
                else:
                    sim = sim.loc[Obs.saveB-1:Obs.saveB+Obs.saveL-2]

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

            if runs.runOK(Obs, Opti, Config, mode) == 1:
                hskip = Obs.nts+3
                idx = np.argsort(np.array(Obs.sim_order))[Obs.obs[oname]
                                                          ['sim_pts']-1]

                # Read
                sim = pd.read_table(Obs.obs[oname]['sim_file'],error_bad_lines=False,
                                    skiprows=hskip, header=None).set_index(0)
                # Check if it ran properly (sometimes some time steps are skipped in *tab...)
                # If steps are missing but the series is long enough, let's replace with nan
                if len(sim) < Obs.lsim and len(sim) > 365:
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
                    f_out.write(str(it+1)+','+','.join([str(j) for j in
                                                        list(sim)])+'\n')
            else:
                # If run failed, write nan line
                with open(Obs.obs[oname]['sim_hist'], 'a') as f_out:
                    f_out.write(str(it+1)+','+','.join(list(np.repeat('nan',
                                                                      Obs.saveL)))+'\n')

        # Fixed-value (initial-value) maps ------------------------------------
        if Obs.obs[oname]['type'] == 'map':

            # Missing vaue for PCraster to numpy conversion
            MV = -9999.

            f_m = Config.PATH_EXEC+'/'+Obs.obs[oname]['sim_file']+'.map'
            if(len(glob.glob(f_m)) == 0):
                print("Warning: the map " + f_m +
                      " seems to be missing from the EcH2O outputs...")
                continue

            # Now that we have what we need, read the PCraster map...
            var_val = \
                pcr.pcr2numpy(pcr.readmap(f_m), MV)*Obs.obs[oname]['sim_conv']

            # Write output NCDF file
            ncFile = Config.PATH_OUT+'/'+oname+'_all.nc'
            # -open nc dataset

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

        # Time-varying maps ---------------------------------------------------
        if Obs.obs[oname]['type'] == 'mapTs':

            # print(oname)

            # Missing vaue for PCraster to numpy conversion
            MV = -9999.
            lensuf = 8 - len(Obs.obs[oname]['sim_file'])

            MapNames = []
            itOK = []
            itNotOK = []

            for it2 in range(1, Obs.lsim+1):

                # Only save files beyond the spinup/transient period (if any)
                if it2 >= Obs.saveBmap and \
                   it2 < Obs.saveBmap+Obs.saveL:

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
            var_t = np.array([(Config.treal[x-Obs.saveBmap] -
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
                          " maps are missing, before t=", itOK[0])
                    print("(that's normal if maps output does not start " +
                          "at t=saveBmap in EcH2O config file)")
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
        os.system('mv '+Config.PATH_EXEC+'/BasinSummary.txt ' +
                  Config.PATH_OUT+'/BasinSummary_run'+str(it+1)+'.txt')

        if Site.isTrck == 1:
            # deuterium summary
            if len(glob.glob(Config.PATH_EXEC+'/Basind2HSummary.txt')) != 0:
                os.system('mv '+Config.PATH_EXEC + '/Basind2HSummary.txt ' +
                          Config.PATH_OUT+'/Basind2HSummary_run' +
                          str(it+1)+'.txt')
            # oxygen 18 summary
            if len(glob.glob(Config.PATH_EXEC+'/Basind18OSummary.txt')) != 0:
                os.system('mv '+Config.PATH_EXEC + '/Basind18OSummary.txt ' +
                          Config.PATH_OUT+'/Basind18OSummary_run' +
                          str(it+1)+'.txt')
            # age summary
            if len(glob.glob(Config.PATH_EXEC + '/BasinAgeSummary.txt')) != 0:
                os.system('mv '+Config.PATH_EXEC + '/BasinAgeSummary.txt ' +
                          Config.PATH_OUT+'/BasinAgeSummary_run' +
                          str(it+1)+'.txt')

# ==============================================================================
# -- Store goodness-of-fit using using several metrics:
#    NSE, KGE, RMSE, MAE


def store_GOF(Obs, Opti, Config, Site, it):
        
    # Did it run OK?
    if runs.runOK(Obs, Opti, Config, mode='verbose') == 1:
        # Read outputs
        if Obs.nobs == 1:
            simulations = read_sim(Config, Obs, Obs.names[0], it)
            # if(len(simulations) < Obs.saveL):
            #     print('sim length:', len(simulations), ', expected:', Obs.saveL)
            #     simulations = [np.nan] * Obs.saveL
        else:
            simulations = np.full((Obs.nobs, Obs.saveL), np.nan).tolist()
            for i in range(Obs.nobs):
                oname = Obs.names[i]
                simulations[i][:] = read_sim(Config, Obs, oname, it)
                # if(len(simulations[i]) < Obs.saveL):
                #     print('sim length:', len(simulations[i]), ', expected:', Obs.saveL)
                #     simulations[i][:] = [np.nan] * Obs.saveL

        for i in range(Obs.nobs):

            oname = Obs.names[i]

            # Have obervation and simulations matching the same time period
            # obs: already pre-processed
            tobs = pd.to_datetime(Opti.obs[oname]['Date'].values)
            o = np.asanyarray(Opti.obs[oname]['value'].values)
            
            # First step for sim: trim sim to obs timespan
            # + only keep dates with obs (even if nan)
            # print(sim)
            # print(sim.shape)
            if Obs.nobs == 1:
                s = np.asanyarray([simulations[j] for j in range(Obs.saveL) if
                                   Obs.simt[j] in tobs])
            else:
                s = np.asanyarray([simulations[i][j] for j in range(Obs.saveL) if
                                   Obs.simt[j] in tobs])
            # Second step (both o and s): remove nan due to gaps in obs
            # (or missing steps in sims...)
            tmp = s*o
            s = np.asanyarray([s[k] for k in range(len(tmp)) if not
                               np.isnan(tmp[k])])
            o = np.asanyarray([o[j] for j in range(len(tmp)) if not
                               np.isnan(tmp[j])])


            # Prepare lists of GOFs
            if i==0:
                gofs = {}
                for gof in Opti.GOFs:
                    gofs[gof] = []
        
            # Another sanity check: if there any data/sim left after nan screening?
            if s.__len__() == 0 or o.__len__() == 0:
                print('Warning: nothing to compare to after date trimming!')
                # Add nan
                for gof in Opti.GOFs:
                    gofs[gof] += [np.nan]
            else:
                # Now use your favorite likelihood estimator for each obs type
                for gof in Opti.GOFs:
                    if gof == 'NSE':
                        gofs[gof] += [GOFs.nash_sutcliffe(s,o)] # NSE
                    if gof == 'KGE':
                        gofs[gof] += [GOFs.kling_gupta(s,o)] # KGE 2009
                    if gof == 'KGE2012':
                        gofs[gof] += [GOFs.kling_gupta(s,o,method='2012')] # KGE 2012
                    if gof == 'RMSE':
                        gofs[gof] += [GOFs.rmse(s,o)] # RMSE
                    if gof == 'MAE':
                        gofs[gof] += [GOFs.meanabs(s,o)] # MAE
                    if gof == 'RMSEc':
                        gofs[gof] += [GOFs.rmse(s-np.mean(s),o-np.mean(o))] # RMSE
                    if gof == 'MAEc':
                        gofs[gof] += [GOFs.meanabs(s-np.mean(s),o-np.mean(o))] # MAE

        # Store goodnesses of fit, one files per GOF
        for gof in Opti.GOFs:
            # GOF
            with open(Opti.GOFfiles[gof], 'a') as f_out:
                f_out.write(str(it+1)+','+','.join([str(x) for x in gofs[gof]])+'\n')

    else:
        # Write NaN
        for gof in Opti.GOFs:
            with open(Opti.GOFfiles[gof], 'a') as f_out:
                f_out.write(str(it+1)+','+','.join(list(np.repeat('nan',Obs.nobs)))+'\n')
    
