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
# -- Updating inputs for ECH2O


def sim_inputs(Opti, Paras, Site, Config, it):

    # Small switch not to read the vegetation params file every time
    # readveg = 1

    # -- Get the parameter sample of the current iteration
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

                if Opti.simRock == 1:
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
                           Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

            # No spatial/veg dependence, but channel stuff
            else:
                # print 'Not dependent !!'
                if Paras.ref[pname]['file'] in ['chanwidth', 'chanmanningn']:
                    outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname]]
                elif Paras.ref[pname]['file'] == 'chanparam':
                    outmap = \
                        Site.bmaps['chanmask_NaN']*Opti.x[Paras.ind[pname]]
                else:
                    outmap = Site.bmaps['unit']*Opti.x[Paras.ind[pname]]

                pcr.report(outmap,
                           Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

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
    if 'Porosity0' in Paras.names:
        poros = pcr.readmap(Config.PATH_SPA+'/' +
                            Paras.ref['Porosity0']['file']+'.map')
    elif 'Porosity' in Paras.names:
        poros = pcr.readmap(Config.PATH_SPA+'/'+Paras.ref['Porosity']['file'] +
                            '.map')
    else:
        poros = pcr.readmap(Config.PATH_SPA+'/poros0.map')

    if 'kPorosity' in Paras.names:

        kporos = pcr.readmap(Config.PATH_SPA+'/' +
                             Paras.ref['kPorosity']['file']+'.map')

        if 'HLayer1' in Paras.names:
            dL1 = pcr.readmap(Config.PATH_SPA+'/' +
                              Paras.ref['HLayer1']['file']+'.map')
        else:
            dL1 = pcr.readmap(Config.PATH_SPA+'/soildepth.L1.map')
        if 'HLayer2' in Paras.names:
            dL2 = pcr.readmap(Config.PATH_SPA+'/' +
                              Paras.ref['HLayer2']['file']+'.map')
        else:
            dL2 = pcr.readmap(Config.PATH_SPA+'/soildepth.L2.map')
        if 'Depth' in Paras.names:
            dTot = pcr.readmap(Config.PATH_SPA+'/' +
                               Paras.ref['Depth']['file']+'.map')
        else:
            dTot = pcr.readmap(Config.PATH_SPA+'/soildepth.map')

        # Layer-integrated values from profile
        porosL1 = kporos*poros*(1-pcr.exp(-dL1/kporos))/dL1
        porosL2 = kporos*poros*(pcr.exp(-dL1/kporos) -
                                pcr.exp(-(dL1+dL2)/kporos))/dL2
        porosL3 = kporos*poros*(pcr.exp(-(dL1+dL2)/kporos) -
                                pcr.exp(-dTot/kporos))/(dTot-dL1-dL2)

    else:

        porosL1 = poros
        porosL2 = poros
        porosL3 = poros

    # -- Use a fraction of these porosities as initial soil moisture
    pcr.report(porosL1*0.4, Config.PATH_SPA+'/Init_SWC.L1.map')
    pcr.report(porosL2*0.4, Config.PATH_SPA+'/Init_SWC.L2.map')
    pcr.report(porosL3*0.4, Config.PATH_SPA+'/Init_SWC.L3.map')

    # - Finalizing soil parameterization
    # Check that initial soil moisture is not smaller residual soil
    # tbd... for now just pick the porosity and thetar reange wisely enough

    # - Finalizing the vegetation parameterization
    if Paras.isveg > 0:
        # Equalize leaf turnover and additional turnover due to water and/or
        # temperature stress
        # for iv in range(Site.nv):
        #    vegnew[iv][vegnew['name'].index('TurnovL_MWS')] = \
        # copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
        #    vegnew[iv][vegnew['name'].index('TurnovL_MTS')] = \
        # copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
        # Write the vegetation params file (if there's at least one veg param)
        vegfile = open(Config.PATH_SPA+'/'+Site.vfile, 'w')
        vegfile.write('\t'.join(Opti.vref['header'])+'\n')
        for iv in range(Site.nv):
            vegfile.write('\t'.join(vegnew[iv])+'\n')
        vegfile.write('\t'.join(Opti.vref['footer'])+'\n')
        vegfile.write('\t'.join(Opti.vref['name'])+'\n')
        vegfile.close()
