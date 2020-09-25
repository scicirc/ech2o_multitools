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
import pandas as pd
import sys

# ----------------------------------------------------------------------------
# -- Write parameters values file


def get(Opti, Config, options):

    if Config.mode == 'calib_MCruns':

        # Read parameters sample for this instance of array task
        # Open one file for all samples
        Opti.xpar = np.genfromtxt(Config.FILE_PAR+Config.numsim+'.txt',
                                  delimiter=',', unpack=True)[1::]
        if(len(Opti.xpar[0]) != len(Opti.names)):
            sys.exit("The definition file and input parameter file ain't " +
                     "matching!")
        # print(Opti.xpar.shape)

    # -- Forward ensemble runs: read directly the params from "best params"
    #    file
    elif Config.mode == 'forward_runs':

        # Sanity check
        tmp = list(pd.read_csv(Config.FILE_PAR, header=None).loc[:, 0])
        if tmp != Opti.names:
            print(Opti.names)
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

        # - Mapped parameters
        if Paras.ref[pname]['map'] == 1:

            # Soil unit dependence
            if Paras.ref[pname]['soil'] == 1:
                # print 'Soil dependent !!'

                outmap = Site.bmaps['unit']*0

                # Read each soil map unit and apply param value
                for im in range(Site.ns):
                    outmap += \
                        Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im]]

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
                # print 'Not dependent !!'
                # Channel stuff (redundant?)
                if Paras.ref[pname]['file'] in \
                   ['chanwidth', 'chanmanningn', 'chanparam']:
                    outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname][0]]
                else:
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
            sys.exit('Error: invalid soil/veg flags to update parameter '+pname)

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
    # Finally, use a fraction of these porosities as initial soil moisture
    pcr.report(porosL1*0.9, path_spa + '/' + Site.f_initSWC1)
    pcr.report(porosL2*0.9, path_spa + '/' + Site.f_initSWC2)
    pcr.report(porosL3*0.9, path_spa + '/' + Site.f_initSWC3)
