#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Python scripts for ECH2O
#
# -------
# Routine: Subroutines for outputs management
# -------
# Contributors: S. Kuppel, A.J. Neill
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

# ----------------------------------------------------------------------------
# -- Read a given simulation output (time series only)


def read_sim(Config, Data, oname):

    # Point-scale time series -------------------------------------
    if Data.obs[oname]['type'] == 'Ts':
        # HEader in EcH2O files
        hskip = Data.nts+3
        idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]
                                                   ['sim_pts']-1]+1
        # Read
        # tmp = np.genfromtxt(Data.obs[oname]['sim_file'],
        #                     delimiter='\t', skip_header=hskip,
        #                     unpack=True)[idx] * \
        #     Data.obs[oname]['sim_conv']
        tmp = pd.read_table(Data.obs[oname]['sim_file'],
                            skiprows=hskip, header=None).iloc[:, idx] * \
            Data.obs[oname]['sim_conv']
        # Shave off the transient part
        # if Config.trimB > 1:
        #    sim = tmp[Config.trimB-1:Config.trimB-1+Config.trimL]
        if Data.lspin > 1:
            sim = tmp[Data.lspin-1:Data.lsim-1]
        # print(sim)

    # Integrated variables (in BasinSummary.txt) -----------------
    if Data.obs[oname]['type'] == 'Total':
        idx = Data.obs[oname]['sim_pts']-1
        tmp = np.genfromtxt(Data.obs[oname]['sim_file'],
                            delimiter='\t', skip_header=1,
                            unpack=True)[idx] * \
            Data.obs[oname]['sim_conv']
        # Shave off the transient part
        # if Config.trimB > 1:
        #     sim = tmp[Config.trimB-1:Config.trimB-1+Config.trimL]
        if Data.lspin > 1:
            sim = tmp[Data.lspin-1:Data.lsim-1]

    return sim
# ----------------------------------------------------------------------------
# -- Store in files for later use


def store_sim(Data, Opti, Config, it):

    # -- Report the full BasinSummary.txt files?
    # if Config.repBS == 1:
    #    os.system('mv '+Config.PATH_EXEC+'/BasinSummary.txt '+
    # Config.PATH_OUT+'/BasinSummary_run'+str(it+1)+'.txt')

    # -- Group the output files in one across simulations,
    #    separating by observations points and veg type where it applies
    for oname in Data.names:
        if Data.obs[oname]['type'] != 'map' and \
           Data.obs[oname]['type'] != 'mapTs' and (it == 0 or
                                                   Opti.begfail == 1):
            # Historic time series file names
            Data.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'
            # Header of files
            with open(Data.obs[oname]['sim_hist'], 'w') as f_out:
                f_out.write('Sample,'+','.join([str(i+1) for i in
                                                range(Config.trimL)])+'\n')

    # Reinit begfail (otherwise will never write all!)
    Opti.begfail = 0

    # Save current run outputs (and delete the file to relieve Maxwell...)
    for oname in Data.names:
        # print(oname)

        # Integrated variables (in BasinSummary.txt) --------------------------
        if Data.obs[oname]['type'] == 'Total':
            idx = Data.obs[oname]['sim_pts']-1

            tmp = np.genfromtxt(Data.obs[oname]['sim_file'], delimiter='\t',
                                skip_header=1,
                                unpack=True)[idx]*Data.obs[oname]['sim_conv']

            # Shave off the transient part (if any)
            if Config.trimB > 1:
                tmp = tmp[Config.trimB-1:Config.trimB-1+Config.trimL]
                if len(tmp) != Config.trimL:
                    sys.exit("ERROR -> Problem with output trim: we've got" +
                             str(len(tmp)) + ' instead of '+str(Config.trimL))

            with open(Data.obs[oname]['sim_hist'], 'a') as f_out:
                f_out.write(str(it+1)+','+','.join([str(j) for j in
                                                    list(tmp)])+'\n')

        # Time series ---------------------------------------------------------
        if Data.obs[oname]['type'] == 'Ts':
            hskip = Data.nts+3
            idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]
                                                       ['sim_pts']-1]+1

            tmp = np.genfromtxt(
                Data.obs[oname]['sim_file'], delimiter='\t',
                skip_header=hskip, unpack=True)[idx]*Data.obs[oname]['sim_conv']

            # Shave off the transient part (if any)
            if Config.trimB > 1:
                tmp = tmp[Config.trimB-1:Config.trimB-1+Config.trimL]

            with open(Data.obs[oname]['sim_hist'], 'a') as f_out:
                f_out.write(str(it+1)+','+','.join([str(j) for j in
                                                    list(tmp)])+'\n')

        # Fixed-value (initial-value) maps ------------------------------------
        if Data.obs[oname]['type'] == 'map':

            # Missing vaue for PCraster to numpy conversion
            MV = -9999.

            f_m = Config.PATH_EXEC+'/'+Data.obs[oname]['sim_file']+'.map'
            if(len(glob.glob(f_m)) == 0):
                print("Warning: the map " + f_m +
                      " seems to be missing from the EcH2O outputs...")
                continue

            # Now that we have what we need, read the PCraster map...
            var_val = \
                pcr.pcr2numpy(pcr.readmap(f_m), MV)*Data.obs[oname]['sim_conv']

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
        if Data.obs[oname]['type'] == 'mapTs':

            # print(oname)

            # Missing vaue for PCraster to numpy conversion
            MV = -9999.
            lensuf = 8 - len(Data.obs[oname]['sim_file'])

            MapNames = []
            itOK = []
            itNotOK = []

            for it2 in range(1, Data.lsim+1):

                # Only save files beyond the spinup/transient period (if any)
                if it2 >= Config.trimBmap and \
                   it2 < Config.trimBmap+Config.trimL:

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
                    if len(Data.obs[oname]['sim_file'])+len(suf) > 12:
                        excess = len(Data.obs[oname]['sim_file'])+len(suf) - 12
                        if len(suf) > 12:
                            print("Warning: the time step is too high to " +
                                  "sort the maps!")
                        f_m = Config.PATH_EXEC+'/' +\
                            Data.obs[oname]['sim_file'][:-excess]+suf
                    else:
                        f_m = Config.PATH_EXEC+'/' +\
                            Data.obs[oname]['sim_file'] + suf

                    if len(glob.glob(f_m)) == 0:
                        itNotOK += [it2]
                    else:
                        MapNames += [f_m]
                        itOK += [it2]
                        # print(f_m)

            # Time values for netCDF output
            var_t = np.array([(Config.treal[x-Config.trimBmap] -
                               datetime(1901, 1, 1, 0, 0)).days for x in itOK])

            if(len(itOK) == 0):
                print("Nothing to be read from the "+oname+" maps...")
                print("it could be due to an incorrect template map name, or" +
                      "a trimBmap value beyond the last simulation timestep?")
                continue

            if(len(itNotOK) > 0):
                if(len(itOK) > 0):
                    print("Warning: some of the demanded "+oname +
                          " maps are missing, before t=", itOK[0])
                    print("(that's normal if maps output does not start " +
                          "at t=trimBmap in EcH2O config file)")
                else:
                    print("Warning: all of the demanded "+oname +
                          " maps are missing!")

            # Second now that we have what we need...
            for it2 in range(len(itOK)):
                # Read map at first time step of interest, convert to array
                # using a missing value,
                # and add an extra 3rd dimension (empty) for later appending
                if(it2 == 0):
                    var_val = pcr.pcr2numpy(
                        pcr.readmap(MapNames[it2]),
                        MV)[None, ...]*Data.obs[oname]['sim_conv']
                # Read subsequent map, same procedure and then append
                else:
                    var_val = \
                        np.append(var_val,
                                  pcr.pcr2numpy(pcr.readmap(MapNames[it2]),
                                                MV)[None, ...], axis=0)
                # Clean up
                os.system('rm -f '+MapNames[it2])

            # Write output NCDF file
            ncFile = Config.PATH_OUT+'/'+oname+'_all.nc'
            # print(ncFile)
            # -open nc dataset
            # If first run, create file
            if(it == 0):
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

            # Write the actual values for this run
            rootgrp = spio.netcdf_file(ncFile, 'a')
            # - write data
            ncVariable = rootgrp.variables[oname]
            ncVariable[:, :, :, it] = var_val
            # -update file and close
            rootgrp.sync()
            rootgrp.close()

            # print
