#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Subroutines for spinup (!depreciated now!)
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import time, os, glob, sys, copy

# ----------------------------------------------------------------------------
# -- Spinup + clean up/edit after

def spinup(Config):

    # Run spinup (1Y-run)
    os.system(Config.spin_ech2o+' > ech2o_spin.log')

    # Check if it ran properly
    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0:
        sys.exit("Something went wrong in the spinup, BasinSummary.txt is missing...")
    # 2. Check it ran until the end
    # get the last day number
    if Config.leap == 1:
        lspin = 366*Config.spinup
    else:
        lspin = 365*Config.spinup
    #format for ECH2O map outputs
    if lspin < 1000:
        espin = '0.'+format(lspin,'03')
    if lspin>=1000 :
        espin = '1.'+format(lspin-1000,'03')
    if lspin>=2000 :
        espin = '2.'+format(lspin-2000,'03')
    if lspin>=3000 :
        espin = '3.'+format(lspin-3000,'03')
    if lspin>=4000 :
        espin = '4.'+format(lspin-4000,'03')
    if lspin>=5000 :
        espin = '4.'+format(lspin-5000,'03')
    if lspin>=6000 :
        espin = '4.'+format(lspin-6000,'03')
    if lspin>=7000 :
        espin = '4.'+format(lspin-7000,'03')

    if len(glob.glob(Config.PATH_EXEC+'/Q000000'+espin)) == 0:
        sys.exit("Something went wrong in the spinup, the last-step-map Q000000"+espin+" is missing...")

    # Copy last-step maps to Spatial directory
    os.system('cp Q000000'+espin+' '+Config.PATH_SPA+'/Q.map')      # Initial streamflow
    os.system('cp SWE0000'+espin+' '+Config.PATH_SPA+'/SWE.map')    # Initial snow water equiv
    os.system('cp SWC1_00'+espin+' '+Config.PATH_SPA+'/SWC.L1.map') # Initial soil moisture L1
    os.system('cp SWC2_00'+espin+' '+Config.PATH_SPA+'/SWC.L2.map') # Initial soil moisture L2
    os.system('cp SWC3_00'+espin+' '+Config.PATH_SPA+'/SWC.L3.map') # Initial Soil moisture L3
    os.system('cp Ts00000'+espin+' '+Config.PATH_SPA+'/Ts.map')     # Initial soil temperature

    if Config.isTrck == 1:
        # -- For initial isotopes, use measurement-derived knowledge as initial value
        # Snowpack, approximate generic value, may not have huge influence
        os.system('cp '+Config.PATH_SPA+'/dD_snowpack.map '+Config.PATH_SPA+'/dD.snowpack.map') 
        os.system('cp '+Config.PATH_SPA+'/d18O_snowpack.map '+Config.PATH_SPA+'/d18O.snowpack.map') 
        # Stream : extrapolated outlet value
        os.system('cp '+Config.PATH_SPA+'/dD.stream.20130221.map '+Config.PATH_SPA+'/dD.surface.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.stream.20130221.map '+Config.PATH_SPA+'/d18O.surface.map')
        # Soil : Josie's measurements (peat/gley/podzol)
        # Podzol extrapolated to ranker
        # L3 is taken as L2 (measurements depth)
        os.system('cp '+Config.PATH_SPA+'/dD.s10cm.20130221.map '+Config.PATH_SPA+'/dD.L1.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.s10cm.20130221.map '+Config.PATH_SPA+'/d18O.L1.map')
        os.system('cp '+Config.PATH_SPA+'/dD.s20cm.20130221.map '+Config.PATH_SPA+'/dD.L2.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.s20cm.20130221.map '+Config.PATH_SPA+'/d18O.L2.map')
        os.system('cp '+Config.PATH_SPA+'/dD.s40cm.20130221.map '+Config.PATH_SPA+'/dD.L3.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.s40cm.20130221.map '+Config.PATH_SPA+'/d18O.L3.map')
        # Groundwater: extrapolation of Bernhard's deep wells (DW1-4)
        # (although it's not the same year, as it is very stable)
        os.system('cp '+Config.PATH_SPA+'/dD.GW.DW201603.map '+Config.PATH_SPA+'/dD.GW.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.GW.DW201603.map '+Config.PATH_SPA+'/d18O.GW.map')

        # -- For initial water age, use spinup
        os.system('cp Agesnw0'+espin+' '+Config.PATH_SPA+'/Age.snowpack.map')
        os.system('cp Agesrf0'+espin+' '+Config.PATH_SPA+'/Age.surface.map')
        os.system('cp AgesL10'+espin+' '+Config.PATH_SPA+'/Age.L1.map')
        os.system('cp AgesL20'+espin+' '+Config.PATH_SPA+'/Age.L2.map')
        os.system('cp AgesL30'+espin+' '+Config.PATH_SPA+'/Age.L3.map')
        os.system('cp Agegw00'+espin+' '+Config.PATH_SPA+'/Age.GW.map') 

    # Clean up
    os.system('rm -f *.txt *.tab *'+espin)

