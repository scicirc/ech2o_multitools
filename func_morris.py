#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Subroutines for parameters & obs manip & outputs
# -------
# Contributors: S. Kuppel, AJ Neill
# Created on 10/2016
# -------------------------------------------------

import os
import glob
import sys
import copy
import csv
import numpy as np
import pyDOE

# --------------------------------------------------------------------------------------------
# -- Creation of the Morris trajectories
# Follows methodology and recommandations of
# Sohier, Farges, and Piet-Lahanier (2014), Improvements of the
# representativity of the Morris
# method for air-launch-to-orbit separation, Proc. 19th congress of IFAC.


def trajs(Config, Opti):

    # Number of levels -> equal number of trajectories/radial points
    # This derives from a recommendation of Sohier et al. (IFAC, 2014)
    Opti.nlev = copy.copy(Opti.nr)

    # vals = {}

    # Step: plus-minus 0.5 of the range of each parameter
    Opti.step = np.zeros((Opti.nvar), np.float32)

    # Construct B* matrix, for each repetition
    Opti.Bnorm = np.zeros((Opti.nvar, Opti.nvar+1, Opti.nr),
                          np.float32)  # the normalized
    Opti.Bstar = np.zeros((Opti.nvar, Opti.nvar+1, Opti.nr),
                          np.float32)  # the final used in runs

    # Starting point: latin hypercube sampling, maximizing 'distance'
    # between point, and centered in the nvar intervals
    Opti.Bnorm[:, 0, :] = np.transpose(pyDOE.lhs(Opti.nvar,
                                                 samples=Opti.nr,
                                                 criterion='cm'))

    # Construct samples
    for ir in range(Opti.nr):
        for iv in range(Opti.nvar):

            # Mode 1 : trajectories
            # Mode 2 : radial points
            # In both cases the other of one-at-a-time change is fixed:
            # 1st param changes,
            # then 2nd param, etc.
            # the randomness is assured by the initial LHS + the fixed step of
            # +-0.5

            if(Opti.MSspace == 'trajectory'):
                # copy previous location
                Opti.Bnorm[:, iv+1, ir] = copy.copy(Opti.Bnorm[:, iv, ir])
            elif(Opti.MSspace == 'radial'):
                # alway start from initial
                Opti.Bnorm[:, iv+1, ir] = copy.copy(Opti.Bnorm[:, 0, ir])
            else:
                sys.exit('Wrong option for the MS parameter space definition!')

            # Successive changes with + (resp. -) 0.5, depending on if they
            # are above (resp. below) the mid-interval
            if Opti.Bnorm[iv, iv+1, ir] < 0.5:
                Opti.Bnorm[iv, iv+1, ir] += 0.5
            elif Opti.Bnorm[iv, iv+1, ir] >= 0.5:
                Opti.Bnorm[iv, iv+1, ir] -= 0.5

            # Check for error
            if Opti.Bnorm[iv, iv+1, ir] > 1 or \
               Opti.Bnorm[iv, iv+1, ir] <=0 :
                print('Error in the incrementation of the parameter',
                      Opti.names[iv])
                print(1/(2*Opti.nlev),Opti.Bnorm[iv, iv, ir], 
                      Opti.Bnorm[iv, iv+1, ir], 1-1/(2*Opti.nlev))
                sys.exit()

    # Construct the actual Bstar, with non-normalized values
    for iv in range(Opti.nvar):
        if Opti.log[iv] == 1:
            Opti.Bstar[iv, :, :] = 10**(Opti.Bnorm[iv, :, :]*np.log10(
                Opti.max[iv]/Opti.min[iv])+np.log10(Opti.min[iv]))
            Opti.step[iv] = \
                0.5 * (np.log10(Opti.max[iv])-np.log10(Opti.min[iv]))
        else:
            Opti.Bstar[iv, :, :] = \
                Opti.Bnorm[iv, :, :]*(Opti.max[iv]-Opti.min[iv]) + Opti.min[iv]
            Opti.step[iv] = 0.5 * (Opti.max[iv]*Opti.min[iv])

    # Check if outputs directory exists
    #if len(glob.glob(Config.PATH_TRAJ)) == 0:
    #    os.system('mkdir ' + Config.PATH_TRAJ)

    # Write Bstar for each trajectory
    # for ir in range(Opti.nr):
    #     trajnb = str(ir+1)  # '%02i' % int(ir+1)
    #     # print trajnb
    #     with open(Config.FILE_TRAJ+'.Bstar_traj'+trajnb+'.txt', 'wb') as fw:
    #         csv_writer = csv.writer(fw)
    #         csv_writer.writerow(Opti.names)
    #         for irun in range(Opti.nvar+1):
    #             csv_writer.writerow(Opti.Bstar[:, irun, ir])
    #     exit

# ----------------------------------------------------------------------------
# -- Calculation of elementary effects for Morris Sensitivity Analysis
# Applied to one trajectory / radial points with same origin
# Since time series are the outputs, metrics for d(outputs) are bias and RMSE

# !! Pandas unavailable so using numpy (unhandy) numpy for I/O+data-manip
# EDIT : pandas available, has to be reverted to it
# (in theory just uncomment previous code, but testing would be necessary)


def ee(Config, Obs, Opti, itraj):

    firstObs = 0
    numObs = 0
    outObs = []

    trajnb = str(itraj+1)

    for oname in Obs.names:

        # Only look into time series
        if Obs.obs[oname]['type'] != 'map' and \
           Obs.obs[oname]['type'] != 'mapTs':

            # Read file
            f_in = Obs.obs[oname]['sim_hist']
            df_sim = pd.read_csv(f_in,skiprows=1)

            # Diff between two sims
            df_diff = df_sim.set_index('Sample').diff().loc[2:Opti.nvar+1,]

            # Get bias
            bias = df_diff.mean(axis=1)
            # Get RMSE
            RMSE = np.sqrt((df.diff**2).mean(axis=1))
            # Get corresponding elementary effect
            bias_ee = bias / Opti.stepN
            RMSE_ee = RMSE / Opti.stepN
            # Write the files
            bias_ee.index = Opti.names
            RMSE_ee.index = Opti.names

            # sim = np.genfromtxt(f_in, delimiter=',', skip_header=1,
            #                     unpack=True)[1:Config.trimL+1, :]

            # # Take into account accumulated fluxes
            # if Obs.obs[oname]['type'] == 'Total' and \
            #    Obs.obs[oname]['sim_pts'] in [1, 11, 12, 13, 14, 15, 16,
            #                                   17, 18, 19, 20]:
            #     sim = np.diff(sim, axis=0)

            # # Diff between sims
            # if(Opti.MSspace == 'trajectory'):
            #     simd = np.diff(sim)
            # elif(Opti.MSspace == 'radial'):
            #     simd = sim[:, 1::] - sim[:, 0][..., None]

            # # Elementary effect (keep direction of param change)
            # bias_ee = np.zeros((Opti.nvar), np.float32)*np.nan
            # RMSE_ee = np.zeros((Opti.nvar), np.float32)*np.nan
            # for i in range(Opti.nvar):
            #     bias_ee[i] = np.mean(simd[:, i]) / Opti.dx[i, i]
            #     RMSE_ee[i] = np.sqrt(np.mean(simd[:, i]**2)) / Opti.dx[i, i]
            # # bias_ee = bias / Opti.BnormstepN
            # # RMSE_ee = RMSE / Opti.stepN

            # Add the overall data frame
            if(firstObs == 0):
                # bias_ee_tot = bias_ee[..., None]  # Creates a (..,1) dimension
                # RMSE_ee_tot = RMSE_ee[..., None]  # Creates a (..,1) dimension
                bias_ee_tot = pd.DataFrame(bias_ee).assign(oname=bias_ee)
                RMSE_ee_tot = pd.DataFrame(RMSE_ee).assign(oname=RMSE_ee)
            else:
                # bias_ee_tot = np.append(bias_ee_tot, bias_ee[..., None], 1)
                # RMSE_ee_tot = np.append(RMSE_ee_tot, RMSE_ee[..., None], 1)
                bias_ee_tot = bias_ee_tot.assign(oname=bias_ee)
                RMSE_ee_tot = RMSE_ee_tot.assign(oname=RMSE_ee)

            # Update
            firstObs = 1
            # Increment number of obs actually evaluated
            numObs += 1
            # Append name of obs actually evaluated
            outObs = outObs + [oname]

    # Ugly: drop the first column (named 0) that had to be left
    # (duplicate of 2nd column)
    bias_ee_tot.drop(0,axis=1,inplace=True)
    RMSE_ee_tot.drop(0,axis=1,inplace=True)
    # Write outputs -----------------------------------------------------------

    # Check if directory exists
    #if len(glob.glob(Config.PATH_EE)) == 0:
    #    os.system('mkdir ' + Config.PATH_EE)

    if(Opti.MSspace == 'trajectory'):
        bias_ee_tot.to_csv(Config.PATH_OUT+'/EE.Traj'+trajnb+'.bias.txt')
        RMSE_ee_tot.to_csv(Config.PATH_OUT+'/EE.Traj'+trajnb+'.RMSE.txt')
        # with open(Config.FILE_EE+'.EE.Traj'+str(itraj+1) + 
        #           '.bias.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(bias_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

        # with open(Config.FILE_EE+'.EE.Traj'+str(itraj+1) +
        #           '.RMSE.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(RMSE_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

    if(Opti.MSspace == 'radial'):
        bias_ee_tot.to_csv(Config.PATH_OUT+'/EE.RadP'+trajnb+'.bias.txt')
        RMSE_ee_tot.to_csv(Config.PATH_OUT+'/EE.RadP'+trajnb+'.RMSE.txt')
        # with open(Config.FILE_EE+'.EE.RadP'+str(itraj+1) +
        #           '.bias.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(bias_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')

        # with open(Config.FILE_EE+'.EE.RadP'+str(itraj+1) +
        #           '.RMSE.txt', 'w') as f_out:
        #     f_out.write('Parameter'+','+','.join([outObs[j] for j in
        #                                           range(numObs)])+'\n')
        #     for i in range(Opti.nvar):
        #         f_out.write(Opti.names[i]+',' +
        #                     ','.join([str(RMSE_ee_tot[i, j]) for j in
        #                               range(numObs)])+'\n')
