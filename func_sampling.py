#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Python routines for ECH2O
#
# -------
# Routine: Subroutines for variable sampling
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import sys
import copy
import numpy as np
import pyDOE
# import random

# ----------------------------------------------------------------------------
# -- Generate random set of parameters values of write it


def MC_sample(Opti, Config):

    # Opti.xtot = np.arange(1,Opti.nsamptot+1)

    # Latin Hypercube sampling
    if Config.sampling in ['LHS', 'LHS_m', 'LHS_r']:

        print('...using a latin hypercube sampling...')

        # 'Normal' LHS
        if Config.sampling == 'LHS':
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(pyDOE.lhs(Opti.nvar, samples=Opti.nsamptot))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling == 'LHS_m':
            print('...with maximin criterion -- it will take longer and a ' +
                  'lot of memory')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(pyDOE.lhs(Opti.nvar, samples=Opti.nsamptot,
                                         criterion='m'))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling == 'LHS_r':
            print('...with correlation criterion -- it will take longer and a lot of memory')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(lhs(Opti.nvar,samples=Opti.nsamptot,criterion='corr'))


        print('...LHS matrix generated...')
        
        # Second, scale with the actual range
        for i in range(Opti.nvar):
            # Log transform where needed
            if Opti.log[i]==1:                
                tmp = 10**(mat[i]*np.log10(Opti.max[i]/Opti.min[i])+np.log10(Opti.min[i]))
            else:
                tmp = mat[i]*(Opti.max[i]-Opti.min[i]) + Opti.min[i]
                
            
            if i==0:
                Opti.xtot = copy.copy(tmp)
            else:
                Opti.xtot = np.vstack((Opti.xtot,tmp))

    #-- Uniform distribution
    elif Config.sampling=='uniform':

        print('...using uniform distributions...')

        for i in range(Opti.nvar):
            # Log transform where needed
            if Opti.log[i]==1:
                if i==0:
                    Opti.xtot = 10**(np.random.uniform(np.log10(Opti.min[i]),np.log10(Opti.max[i]),Opti.nsamptot))
                else:
                    Opti.xtot = np.vstack((Opti.xtot,10**(np.random.uniform(np.log10(Opti.min[i]),
                                                                            np.log10(Opti.max[i]),
                                                                            Opti.nsamptot))))
            else:
                if i==0:
                    Opti.xtot = np.random.uniform(Opti.min[i],Opti.max[i],Opti.nsamptot)
                else:
                    Opti.xtot = np.vstack((Opti.xtot,np.random.uniform(Opti.min[i],Opti.max[i],Opti.nsamptot)))
        #print i, Opti.names[i], Opti.x[i]

    else:
        sys.exit('No proper sampling method selected!')


    print('Parameters sampling done!')

    # -- Reshape for output
    #Opti.xtot = np.transpose(Opti.xtot)#.shape((Opti.nvar,Opti)
    print(Opti.xtot.shape)

    print
    print('Writing in '+Config.ncpu+' files...('+str(Opti.nit)+' sets each)')

    # -- Write one file per parallel job
    for i in range(int(Config.ncpu)):
        f_out = Config.FILE_PAR+str(i+1)+'.txt'
        k = i*Opti.nit
        #print str(k)+','+str(k+Opti.nit)
        with open(f_out,'w') as fw:
            for j in range(Opti.nvar):
                #print i
                #print j
                #print k
                tmp=[a for a in Opti.xtot[j][k:k+Opti.nit]]
                #print len(tmp)
                fw.write(Opti.names[j]+','+','.join([str(a) for a in tmp])+'\n')

    # Write on file giving parameters range, log...(for later plots)
    f_out = Config.FILE_PAR+'char.txt'
    with open(f_out,'w') as fw:
        # fw.write('Sample,'+','.join(Opti.names)+'\n')
        fw.write('Names,Min,Max,Log\n')
        for i in range(Opti.nvar):
            fw.write(','.join([Opti.names[i],str(Opti.min[i]),str(Opti.max[i]),str(Opti.log[i])])+'\n')

    print('')
    
    #np.savetxt(f_out,np.transpose(Opti.xtot))#,header= 

