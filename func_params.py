#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Python scripts for ECH2O
#
# -------
# Routine: Subroutines for parameters
# -------
# Contributors: S. Kuppel, AJ Neill
# Created on 10/2016
# -------------------------------------------------

import pcraster as pcr
import copy
import numpy as np
import re

# ----------------------------------------------------------------------------
# -- Write parameters values file


def get(Opti, Config):

    # Open one file for all samples
    f_in = Config.FILE_PAR+Config.numsim+'.txt'
    print(f_in)
    Opti.xpar = np.genfromtxt(f_in, delimiter=',', unpack=True)[1::]
    print(Opti.xpar.shape)

# ----------------------------------------------------------------------------
# -- Write parameters values file


def store(Opti, Config, it):

    # Open one file for all samples
    if Config.initpar == 0:
        Config.f_par = Config.PATH_OUT+'/Parameters.txt'
        if Config.restart == 0:
            with open(Config.f_par, 'w') as f_in:
                f_in.write('Iteration,'+','.join(Opti.names)+'\n')
        Config.initpar = 1

    with open(Config.f_par, 'a') as f_in:
        f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')

# ----------------------------------------------------------------------------
# -- Creating/updating inputs for ECH2O


def sim_inputs(Opti, Paras, Site, path_spa, it=0, mode='no_spotpy',
               paramcur=None):

    # Small switch not to read the vegetation params file every time
    # readveg = 1

    # -- Get the parameter values
    if mode == 'spotpy':
        # In spotpy mode there is no Opti.xpar array with all samples,
        # values are read directly from the spot_setup class
        # In addition, parameter with log10 variation should be "de-log10ned"
        # here!
        # print('Parameter value:')

        Opti.x = []
        for i in range(Opti.nvar):
            if Opti.log[i] == 1:
                # Opti.x += [10**(paramcur()[i][0])]
                Opti.x += [10**paramcur[i]]
            else:
                # Opti.x += [paramcur()[i][0]]
                Opti.x += [paramcur[i]]

            # print('|',Opti.names[i],':',Opti.x[i],end='\r')
    else:
        # Therwise, just get the samples of the current iteration
        Opti.x = Opti.xpar[it]
    

    for pname in Paras.names:

        # print pname

        # - Mapped parameters
        if Paras.ref[pname]['veg'] == 0:

            # Soil unit dependence
            if Paras.ref[pname]['soil'] == 1:
                # print 'Soil dependent !!'
                outmap = Site.bmaps['unit']*0
                # Read each soil map unit and apply param value
                for im in range(Site.ns):
                    outmap += \
                        Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im]]

                if Site.simRock == 1:
                    # Taking into account rock/scree: micro-topsoil, low poros
                    # and fixed anisotropy
                    if pname == 'HLayer1':
                        outmap = \
                            outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
                            Site.bmaps['rock']*0.001
                    # if pname=='Porosity':
                    #    outmap = \
                        # outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
                        # Site.bmaps['rock']*0.25
                    if pname == 'Khoriz':
                        outmap = \
                            outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
                            Site.bmaps['rock']*0.000001
                    if pname == 'Anisotropy':
                        outmap = \
                            outmap*(Site.bmaps['unit']-Site.bmaps['rock']) +\
                            Site.bmaps['rock']*0.1

                pcr.report(outmap,
                           path_spa+'/'+Paras.ref[pname]['file']+'.map')

            # No spatial/veg dependence, but channel stuff
            else:
                # print 'Not dependent !!'
                if Paras.ref[pname]['file'] in ['chanwidth', 'chanmanningn', 'chanparam']:
                    outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname]]
                else:
                    outmap = Site.bmaps['unit']*Opti.x[Paras.ind[pname]]

                pcr.report(outmap,
                           path_spa+'/'+Paras.ref[pname]['file']+'.map')

            print('rank',it, ': map of',pname,
                  'in',path_spa+'/'+Paras.ref[pname]['file']+'.map')
        # - Vegetation parameters
        else:
            # Change the value based on param name correspondance
            vegnew = copy.copy(Opti.vref)
            # print Opti.vref
            for iv in range(Site.nv):
                vegnew[iv][vegnew['name'].index(pname)] = \
                    str(Opti.x[Paras.ind[pname][iv]])

    # ------------------------------------------------------------------------------
    # Finalizing the preps....
    # Check for initial condition/other parameter dependence
    # Back up relevant maps

    # Initial soil water content: needs porosity profile and depth to have
    # consistent initial soil moisture for each layer
    
    # Check which porosity profile mode is on in config file
    pattern = re.compile('Porosity_profile')
    with open (Config.FILE_CFGdest, 'rt') as myfile:    
        for line in myfile:
            if pattern.search(line) != None:      # If a match is found 
                poros_mode = int(line.strip('=')[1].strip())
                break
    print('Porosity profile:',poros_mode)
    # Get base porosity
    pattern = re.compile('Top-of-profile_Porosity')
    with open (Config.FILE_CFGdest, 'rt') as myfile:    
        for line in myfile:
            if pattern.search(line) != None:      # If a match is found 
                poros_file = line.strip('=')[1].strip()
                break
    print('Porosity (base) map file:',poros_file)
    poros = pcr.readmap(path_spa+'/' + poros_file)
    
    # Depending on porosity profile mode, different ways to each layer's porosity
    if poros_mode == 0:
        # Constant profile
        porosL1 = poros
        porosL2 = poros
        porosL3 = poros

    elif poros_mode == 1:
        # Exponential: profile coeff and depths are needed
        pattern1 = re.compile('Porosity_Profile_coeff')
        pattern2 = re.compile('Depth_soil_layer_1')
        pattern3 = re.compile('Depth_soil_layer_2')
        pattern4 = re.compile('Soil_depth')
        count = 0 
        with open (Config.FILE_CFGdest, 'rt') as myfile:    
            for line in myfile:
                if pattern1.search(line) != None:      # If a match is found 
                    kporos_file = line.strip('=')[1].strip()
                    count += 1
                if pattern2.search(line) != None:      # If a match is found 
                    dL1_file = line.strip('=')[1].strip()
                    count += 1
                if pattern3.search(line) != None:      # If a match is found 
                    dL2_file = line.strip('=')[1].strip()
                    count += 1
                if pattern4.search(line) != None:      # If a match is found 
                    dsoil_file = line.strip('=')[1].strip()
                    count += 1
                if count == 4:
                    break
        print('Porosity coeff map file:', kporos_file)
        # Read maps
        kporos = pcr.readmap(path_spa+'/' + kporos_file)
        dL1 = pcr.readmap(path_spa+'/' + dL1_file)
        dL2 = pcr.readmap(path_spa+'/' + dL2_file)
        dsoil = pcr.readmap(path_spa+'/' + dsoil_file)
        # Layer-integrated values from profile
        porosL1 = kporos*poros*(1-pcr.exp(-dL1/kporos))/dL1
        porosL2 = kporos*poros*(pcr.exp(-dL1/kporos) -
                                pcr.exp(-(dL1+dL2)/kporos))/dL2
        porosL3 = kporos*poros*(pcr.exp(-(dL1+dL2)/kporos) -
                                pcr.exp(-dsoil/kporos))/(dsoil-dL1-dL2)
  
    elif poros_mode == 2:
        # Porosity map given for each layer
        # L2
        pattern1 = re.compile('Porosity_Layer2')
        pattern2 = re.compile('Porosity_Layer3')
        count = 0
        with open (Config.FILE_CFGdest, 'rt') as myfile:    
            for line in myfile:
                if pattern1.search(line) != None:      # If a match is found 
                    porosL2_file = line.strip('=')[1].strip()
                    count += 1
                if pattern2.search(line) != None:      # If a match is found 
                    porosL3_file = line.strip('=')[1].strip()
                    count += 1
                if count == 2:
                    break
        # Assign
        porosL1 = poros
        porosL2 = pcr.readmap(path_spa+'/' + porosL2_file)
        porosL3 = pcr.readmap(path_spa+'/' + porosL3_file)
    
    # -- Use a fraction of these porosities as initial soil moisture
    pattern1 = re.compile('Soil_moisture_1')
    pattern2 = re.compile('Soil_moisture_2')
    pattern3 = re.compile('Soil_moisture_3')
    count = 0
    with open (Config.FILE_CFGdest, 'rt') as myfile:    
        for line in myfile:
            if pattern1.search(line) != None:      # If a match is found 
                initSWC1_file = line.strip('=')[1].strip()
                print('key SWC1 found !')
                count += 1
            if pattern2.search(line) != None:      # If a match is found 
                initSWC2_file = line.strip('=')[1].strip()
                print('key SWC2 found !')
                count += 1
            if pattern3.search(line) != None:      # If a match is found 
                initSWC3_file = line.strip('=')[1].strip()
                print('key SWC3 found !')
                count += 1
            if count == 3:
                break
    # Write initial moisture maps
    pcr.report(porosL1*0.9, path_spa + '/' + initSWC1_file)
    pcr.report(porosL2*0.9, path_spa + '/' + initSWC2_file)
    pcr.report(porosL3*0.9, path_spa + '/' + initSWC3_file)

    # - Updating the vegetation parameterization
    if Paras.isveg > 0:
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

    sys.exit()
