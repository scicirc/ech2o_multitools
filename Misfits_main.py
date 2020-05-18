#    !/usr/bin/env python3
# -*- coding: utf-8 -*-
# *************************************************
#
# Model-data misfit after Monte-Carlo sampling +
# assemble parameters samples
#
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import numpy as np
import scipy as sp
import pandas as pd
from scipy.interpolate import interp1d, UnivariateSpline, splrep, splint
from datetime import datetime, timedelta
import os, time, sys, glob, copy
import json
from optparse import OptionParser

# --- Subroutine(s)
import Misfits_funcs as metrics

parser = OptionParser()

# Switch : fetching obs metrics or corresponding parameters / best runs ?
# 0 = metrics + par, 1 = best runs
parser.add_option("--switch",dest="switch",metavar="SWITCH")
parser.add_option("--ext",dest="ext",metavar="EXT")
parser.add_option("--byjob",dest="byjob",metavar="byjob")
parser.add_option("--docopy",dest="docopy",metavar="docopy")
parser.add_option("--clean",dest="clean",metavar="clean")
parser.add_option("--jobs",dest="jobs",metavar="jobs")
parser.add_option("--nruns",dest="nruns",metavar="nruns")

(options, args) = parser.parse_args()

if options.switch != None:
    switch = int(options.switch)
else:
    switch = 0

MCname = copy.copy(options.ext)
sufMC = 'new'

if options.docopy != None:
    docopy = int(options.docopy)
else:
    docopy = 0
if options.clean != None:
    clean = int(options.clean)
else:
    clean = 0
if options.byjob != None:
    byjob = int(options.byjob)
else:
    byjob = 0

#switch = 1
#swpar = 1
#swsim = 1
#nbest = 30

# scratch directory
scrdir = '/scratch/users/s08sk8'

# -- Output directory
outdir = 'Outputs'

# -- Date set-up
tstart = datetime(2013,2,21) # start of 'useful' sim
leff = 1265 # useful sim length
spinup = 0#1095 # spinup length
lsim = leff + spinup # total sim length
simt = [tstart-timedelta(days=spinup)+timedelta(days=x) for x in range(lsim)]

# --- Definition ---

# Path to data
obsdir = '/users/s08sk8/Data/BB'
# Simulation subdirectories
#MCname = 'MC12'
# Observations
Vars = {}

Vars['Streamflow'] = {'sim':['Streamflow_all.tab'],
                      'obs':'/Q/BB_discharge_daily_01062011-01022017.csv','obscol':[1],
                      'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}

Vars['SWC_Peat'] = {'sim':['SWC_Peat.L1_all.tab'],
                   'obs':'/SWC/VSM_Peat_daily.csv','obscol':[1],'obsdph':[0.1],'loc':'P',
                   'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['SWC_Gley'] = {'sim':['SWC_Gley.L1_all.tab','SWC_Gley.L2_all.tab'],
                   'obs':'/SWC/VSM_Gley_daily.csv','obscol':[1,2,3],'obsdph':[0.1,0.2,0.4],'loc':'PG',
                   'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['SWC_Podzol'] = {'sim':['SWC_Podzol.L1_all.tab','SWC_Podzol.L2_all.tab'],
                     'obs':'/SWC/VSM_Podzol_daily.csv','obscol':[1,2,3],'obsdph':[0.1,0.2,0.4],'loc':'PP',
                     'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
#Vars['SWC_ForestA'] = {'sim':['SWC_ForestA.L1_all.tab','SWC_ForestA.L2_all.tab'],
#                      'obs':'/SWC/VSM_ForestA_daily.csv','obscol':[1,2,3],'obsdph':[0.1,0.2,0.4],'loc':'PP',
#                      'beg':datetime(2013,2,21),'end':datetime(2015,8,8)}
Vars['SWC_ForestB'] = {'sim':['SWC_ForestB.L1_all.tab','SWC_ForestB.L2_all.tab'],
                      'obs':'/SWC/VSM_ForestB_daily.csv','obscol':[1,2,3],'obsdph':[0.1,0.2,0.4],'loc':'PP',
                      'beg':datetime(2015,2,25),'end':datetime(2016,8,8)}

Vars['SWC_HeatherA'] = {'sim':['SWC_HeatherA.L1_all.tab','SWC_HeatherA.L2_all.tab'],
                      'obs':'/SWC/VSM_HeatherA_daily.csv','obscol':[1,2,3],'obsdph':[0.1,0.2,0.4],'loc':'PP',
                      'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}

Vars['T_ForestA'] = {'sim':['T_ForestA_all.tab'],
                      'obs':'/ET/T_ForestA_daily.csv','obscol':[1],
                      'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['T_ForestB'] = {'sim':['T_ForestB_all.tab'],
                      'obs':'/ET/T_ForestB_daily.csv','obscol':[1],
                      'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}

Vars['NetRad_WS1'] = {'sim':['NetRad_WS1_all.tab'],
                      'obs':'/NetRad/station1_chickencage_daily_17072014-08082016.csv','obscol':[1],
                      'beg':datetime(2014,10,9),'end':datetime(2016,8,8)}
Vars['NetRad_WS2'] = {'sim':['NetRad_WS2_all.tab'],
                      'obs':'/NetRad/station2_bog_daily_17072014-03082016.csv','obscol':[1],
                      'beg':datetime(2014,10,9),'end':datetime(2016,8,3)}
Vars['NetRad_WS3'] = {'sim':['NetRad_WS3_all.tab'],
                      'obs':'/NetRad/station3_hilltop_daily_17042015-12072016.csv','obscol':[1],
                      'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}


Vars['d2H_Stream'] = {'sim':['d2H_Stream_all.tab'],
                      'obs':'/Isotopes/d2H_Stream.csv','obscol':[1],
                      'beg':datetime(2013,2,21),'end':datetime(2015,2,20)}

Vars['d2H_DW1'] = {'sim':['d2H_DW1_all.tab'],
                   'obs':'/Isotopes/d2H_deeperwells.csv','obscol':[1],
                   'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['d2H_DW2'] = {'sim':['d2H_DW1_all.tab'], # uses same output as DW1 (same pixel at 100m)
                   'obs':'/Isotopes/d2H_deeperwells.csv','obscol':[2],
                   'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['d2H_DW3'] = {'sim':['d2H_DW3_all.tab'],
                   'obs':'/Isotopes/d2H_deeperwells.csv','obscol':[3],
                   'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['d2H_DW4'] = {'sim':['d2H_DW4_all.tab'],
                   'obs':'/Isotopes/d2H_deeperwells.csv','obscol':[4],
                   'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}



Vars['d18O_Stream'] = {'sim':['d18O_Stream_all.tab'],
                       'obs':'/Isotopes/d18O_Stream.csv','obscol':[1],
                       'beg':datetime(2013,2,21),'end':datetime(2015,2,20)}

Vars['d18O_DW1'] = {'sim':['d18O_DW1_all.tab'],
                    'obs':'/Isotopes/d18O_deeperwells.csv','obscol':[1],
                    'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['d18O_DW2'] = {'sim':['d18O_DW1_all.tab'], # uses same output as DW1 (same pixel at 100m)
                    'obs':'/Isotopes/d18O_deeperwells.csv','obscol':[2],
                    'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['d18O_DW3'] = {'sim':['d18O_DW3_all.tab'],
                    'obs':'/Isotopes/d18O_deeperwells.csv','obscol':[3],
                    'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}
Vars['d18O_DW4'] = {'sim':['d18O_DW4_all.tab'],
                    'obs':'/Isotopes/d18O_deeperwells.csv','obscol':[4],
                    'beg':datetime(2013,2,21),'end':datetime(2016,8,8)}


obsnames = sorted(Vars.keys(),key=str.lower)
#obsnames = ['Streamflow',
#            'SWC_Peat']#,'SWC_Gley','SWC_Podzol','SWC_ForestB',
#            'T_ForestA','T_ForestB',
#            'NetRad_WS1','NetRad_WS2','NetRad_WS3']


#print obsnames 
#obsnames = obsnames.sort()
nobs = len(obsnames)
#print obsnames

# Corresponding simulations outputs
simdir = os.getcwd()+'/'+MCname+'_sampling.'
pardir = os.getcwd()+'/Parameters_samples/'#+MCname+'_sampling_parameters.1.txt'

# Unit conversion sim --> obs
# simfct = [1 for i in len(obsnames)]

# Number of samples
# nit = 1500
nit = int(options.nruns)

# Number of parallel runs
#print options.jobs
if byjob == 1 and options.jobs != None:
    jobs = [int(x) for x in options.jobs.split(',')]
else:
    jobs = [13,14,15]
print(jobs)

print '========================================'
print '------ Calculating model-data fits -----'
print
print " Namely: KGE, MAE, RMSE, and Pearson's r"
print 
print 'Obsversations used here:'
for i in range(nobs):
    oname = obsnames[i]
    print oname+': from',Vars[oname]['beg'],'to',Vars[oname]['end']
print
print 'MC jobs included here:'
print ' '.join([str(i) for i in jobs]) 
print

#######################################################################
# Simulations outputs : saving the metrics
# if switch == 0:

#print
#print 'Summarizing the metrics....'
print

Obs={}
obst={}
nok = 0

if switch == 0:

    print 'Store measured datasets...'
    print

    for oname in obsnames:
        
        Obs[oname] = {}

        # -- Get the obs
        f_obs = obsdir+'/'+Vars[oname]['obs']
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
        tmpt = np.array([datetime.strptime(a, '%d/%m/%Y') for a in tmp])
        lobs = len(tmpt)
        # tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True)[Vars[oname]['obscol']]
        # tmp = pd.read_csv(f_obs).iloc[:,Vars[oname]['obscol']]
        parser = lambda date: pd.datetime.strptime(date,'%d/%m/%Y')
        tmp = pd.read_csv(f_obs, parse_dates=[0],
                          date_parser=parser).iloc[:,[0]+Vars[oname]['obscol']]
        # Just to be sure the first column is named DATE
        tmp.columns.values[0]='DATE'
        #tmp = tmp.set_index('DATE')
        #print oname

        # Crop obs between desired time frame
        fitbeg = Vars[oname]['beg']
        fitend = Vars[oname]['end']
        tmp = tmp.loc[(tmp.DATE>=fitbeg) & (tmp.DATE<=fitend)]
        #Obs[oname] = tmp.loc[(tmp.DATE>=fitbeg) & (tmp.DATE<=fitend)]
        Obs[oname] = tmp.dropna(how='all')#inplace=True,how='all')
        #Obs[oname]['v'] = [tmp.iloc[:,idx] for idx in range(lobs) if tmpt[idx]>=fitbeg and tmpt[idx]<=fitend]
        #Obs[oname]['t'] = np.array([tmpt[idx] for idx in range(lobs) if tmpt[idx]>=fitbeg and tmpt[idx]<=fitend])
        #print Obs[oname].DATE
    
    # -- Copy simulations outputs to scratch
    if docopy == 1:
        print 'Temporarily copy parameters to scratch...(outputs will done on-the-fly for each job)'
    
        tmpdir = os.path.abspath(os.path.join(scrdir,MCname))
        if len(glob.glob(tmpdir))==0:
            os.system('mkdir '+tmpdir)
        tmppar = tmpdir+'/'+MCname+'_parameters'
        if len(glob.glob(tmppar))==0:
            os.system('mkdir '+tmppar)
        os.system('cp -p '+pardir+'/'+MCname+'_sampling_parameters.*.txt '+tmppar+'/')


    print
    print 'Preparing summary files...'

    # -- All variables: for soil samples, take into account the one/two layers
    Out = {}
    Indf = {}
    for oname in obsnames:
        if len(Vars[oname]['sim']) > 1:
            for j in range(len(Vars[oname]['sim'])):
                Out[oname+'.'+Vars[oname]['sim'][j].split('_all.tab')[0].split('.')[1]] = oname
                Indf[oname+'.'+Vars[oname]['sim'][j].split('_all.tab')[0].split('.')[1]] = j
        else:
            Out[oname] = oname
            Indf[oname] = 0
    outnames = sorted(Out.keys(),key=str.lower)
    #print outnames

    # Parameters
    f_in = tmppar+'/'+MCname+'_sampling_parameters.char.txt'
    pnames2 = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20',unpack=True,skip_header=1)[0])
    #print tmp
    npar = len(pnames2)
    print
    print '(length of parameter vector: '+str(npar)+')'
    print
    pnames = ','.join([str(pnames2[idx]) for idx in range(npar)])
    # if byjob == 0:
    #     with open(tmpdir+'/'+MCname+'_parameters.txt','w') as f_out:
    #         f_out.write('Job,Iteration,Sample,'+pnames+'\n')
  
    itot = 1
    itot2 = 1

    print 'Calculate sample misfits...'

    for i in jobs:
        
        print 'job',i,'...'

        f_par = tmppar+'/'+MCname+'_sampling_parameters.'+str(i)+'.txt'
        #tmp_par = np.genfromtxt(f_par,delimiter=',')[::,1::]
        df_par = pd.read_csv(f_par,header=None).set_index(0)

        tmpdir2 = tmpdir+'/'+MCname+'_sampling.'+str(i)

        # Copy outputs from /users ?
        if docopy == 1:
            if len(glob.glob(tmpdir2))==0:
                os.system('mkdir '+tmpdir2)
            os.system('cp -p '+simdir+str(i)+'/*.tab '+tmpdir2+'/')

        # Initialize metrics
        KGE = pd.DataFrame(range(1,nit+1), columns=['Iteration'], index=range(1,nit+1))
        MAE = pd.DataFrame(range(1,nit+1), columns=['Iteration'], index=range(1,nit+1))
        RMSE = pd.DataFrame(range(1,nit+1), columns=['Iteration'], index=range(1,nit+1))
        corr = pd.DataFrame(range(1,nit+1), columns=['Iteration'], index=range(1,nit+1))

        # Loop over the datastream examined
        for oname in outnames:

            # corresponding "variable" name
            oname2 = Out[oname]

            f_sim = tmpdir2+'/'+Vars[oname2]['sim'][Indf[oname]]

            # ID the runs that worked
            js1 = pd.read_csv(f_sim)['Sample']
            nj1 = len(js1)
            #print nj1
            js1 = [np.int(js1[idx]) for idx in range(nj1) if js1[idx]<=nit]
            #print js1
            tmp = [js1[idx+1]-js1[idx] for idx in range(nj1-1)]
            # Index of non-doublons in js (if doublons-triplets-etc., take the last one)
            jid = [idx for idx in range(nj1-1) if tmp[idx]>0] + [nj1-1]
            #print jid
            # Corrected list of runs numbers
            js2 = [js1[ix] for ix in jid]
            nj = len(js2)
            #print nj
            #print js2
            print oname, ', runs OK =',nj,'(failed :',nit-nj,')'
            # Read all simulation for this obs / job
            df_sim = pd.read_csv(f_sim,header=None,skiprows=1).transpose().loc[1:lsim,]
            # Subset where there is an obs (even if NA)
            df_sim = df_sim.set_index(pd.date_range(tstart,periods=lsim))
            # Add runs actual number (with doublons)
            df_sim.columns = js1
            # Remove doublons and keep only observation days
            df_sim = df_sim.ix[Obs[oname2].DATE].iloc[:,jid]
            
            # If no transformation, direct fit quantification
            if len(Vars[oname2]['obscol'])==1:
                # KGE
                def func1(x):
                    return metrics.kling_gupta(x,Obs[oname2].iloc[:,1],method='2012')
                KGE[oname] = np.nan
                KGE.loc[js2,oname] = copy.copy(df_sim.apply(func1,axis=0))
                # MAE
                def func2(x):
                    return metrics.meanabs(np.asarray(x),Obs[oname2].iloc[:,1])
                MAE[oname] = np.nan
                MAE.loc[js2,oname] = df_sim.apply(func2,axis=0)
                # RMSE
                def func3(x):
                    return metrics.rmse(np.asarray(x),Obs[oname2].iloc[:,1])
                RMSE[oname] = np.nan
                RMSE.loc[js2,oname] = df_sim.apply(func3,axis=0)
                # Correlation
                def func4(x):
                    return metrics.corr(np.asarray(x),Obs[oname2].iloc[:,1])
                corr[oname] = np.nan
                corr.loc[js2,oname] = copy.copy(df_sim.apply(func4,axis=0))
                
            # Several observation: reconstruct soil profiles for the simulated layer
            else:
                # Simulation-wise for now
                # Convoluted df.apply is maybe possible, but complex
                KGE[oname] = np.nan
                RMSE[oname] = np.nan
                MAE[oname] = np.nan
                corr[oname] = np.nan

                for j in range(nj):

                    # -- Subset to runs that worked
                    sim = df_sim.iloc[:,j]#*simfct[iobs]

                    # Topsoil?
                    if oname.split('.')[1]=='L1':
                        d1 = 0
                        d2 = df_par.loc['HLayer1_'+Vars[oname2]['loc'],js2[j]]
                        #for k in range(npar) if pnames2[k]==]['loc']][0]
                    if oname.split('.')[1]=='L2':
                        d1 = df_par.loc['HLayer1_'+Vars[oname2]['loc'],js2[j]]
                        d2 = d1 + df_par.loc['HLayer2_'+Vars[oname2]['loc'],js2[j]]
                        #d1 = tmp_par[k,js2[j]-1] for k in range(npar) if pnames2[k]=='HLayer1_'+Vars[oname2]['loc']][0]
                        #d2 = d1 + [tmp_par[k,js2[j]-1] for k in range(npar) if pnames2[k]=='HLayer2_'+Vars[oname2]['loc']][0]

                    # Spline interpolation to then get layer-integrated value
                    d = np.asanyarray(Vars[oname2]['obsdph'])
                    obs = np.zeros(df_sim.shape[0])
                    obs[:] = np.nan

                    #for it in range(df_sim.shape[0]):
                    #    func = sp.interpolate.UnivariateSpline(d,Obs[oname2].iloc[it,1:],s=0)
                    #    obs[it] = sp.integrate.quad(func,d1,d2)[0]/(d2-d1)

                    def func(y):
                    #    return sp.integrate.quad(interp1d(d,y,kind='cubic'),d1,d2)/(d2-d1)
                    # return UnivariateSpline(d,y,s=0).integral(d1,d2)/(d2-d1)
                        return splint(d1,d2,splrep(d,y,k=min(len(d)-1,3),s=0))/(d2-d1)
                    # return sp.integrate.quad(interp1d.splrep(d,y,k=len(d)-1,s=0),d1,d2)[0]/(d2-d1)
              
                    obs = Obs[oname2].iloc[:,1:].apply(func,axis=1)
                    #print len(obs)

                    KGE.loc[js2[j],oname] = metrics.kling_gupta(sim,obs,method='2012')
                    MAE.loc[js2[j],oname] = metrics.meanabs(sim,obs)
                    RMSE.loc[js2[j],oname] = metrics.rmse(sim,obs)
                    corr.loc[js2[j],oname] = metrics.corr(sim,obs)       

            if oname == outnames[0]:
                itot+=nj

        # Clean the metrics dataframe to include only the successful runs common to all obs
        # Use MAE or RMSE, because KGE and corr can have NaN only for 'flat' succesful runs
        MAE.dropna(inplace=True)
        RMSE.dropna(inplace=True)

        js3 = MAE.index

        KGE = KGE.ix[js3]
        corr = corr.ix[js3]
        df_par = df_par.loc[:,js3]

        # ====================================================
        # -- Outputs
        #print j, js[j-1], idx, tmp_par.shape
        if byjob == 0:

            # Parameters
            if(i==jobs[0]):
                df_par.to_csv(tmpdir+'/'+MCname+sufMC+'_parameters.txt')
            else:
                df_par.to_csv(tmpdir+'/'+MCname+sufMC+'_parameters.txt',
                              mode='a',header=False)

            # Metrics
            if(i==jobs[0]):
                KGE.to_csv(tmpdir+'/'+MCname+sufMC+'_KGE2012.txt',index=False)
                MAE.to_csv(tmpdir+'/'+MCname+sufMC+'_MAE.txt',index=False)
                RMSE.to_csv(tmpdir+'/'+MCname+sufMC+'_RMSE.txt',index=False)
                corr.to_csv(tmpdir+'/'+MCname+sufMC+'_corr.txt',index=False)
            else:
                KGE.to_csv(tmpdir+'/'+MCname+sufMC+'_KGE2012.txt',index=False,
                           mode='a',header=False)
                MAE.to_csv(tmpdir+'/'+MCname+sufMC+'_MAE.txt',index=False,
                           mode='a',header=False)
                RMSE.to_csv(tmpdir+'/'+MCname+sufMC+'_RMSE.txt',index=False,
                            mode='a',header=False)
                corr.to_csv(tmpdir+'/'+MCname+sufMC+'_corr.txt',index=False,
                            mode='a',header=False)


        
        elif byjob==1:
            
            # Parameters
            df_par.to_csv(tmpdir+'/'+MCname+sufMC+'_job'+str(i)+'_parameters.txt')
            # with open(tmpdir+'/'+MCname+'_job'+str(i)+'_parameters.txt','a') as f_out:
            # #with open(tmpdir+'/'+MCname+'job_'+str(i)+'_parameters.txt','a') as f_out:
            #     for j in range(nj):
            #         f_out.write(str(js[j])+','+','.join([str(x) for x in tmp_par[::,js[j]-1]])+'\n')

            
            #print KGE

            # Metrics
            KGE.to_csv(tmpdir+'/'+MCname+sufMC+'_job'+str(i)+'_KGE2012.txt',index=False)
            MAE.to_csv(tmpdir+'/'+MCname+sufMC+'_job'+str(i)+'_MAE.txt',index=False)
            RMSE.to_csv(tmpdir+'/'+MCname+sufMC+'_job'+str(i)+'_RMSE.txt',index=False)
            corr.to_csv(tmpdir+'/'+MCname+sufMC+'_job'+str(i)+'_corr.txt',index=False)

            # -- Copy to /users
            print 'Copy back to home base...'
            print
            os.system('cp -p '+tmpdir+'/'+MCname+sufMC+'_job'+str(i)+'_*.txt '+os.getcwd()+'/'+outdir+'/')
            #os.system('cp -p '+tmpdir+'/'+MCname+'job_'+str(i)+'_*.txt '+os.getcwd()+'/'+outdir+'/')

            
            if clean == 1:
                print 'Cleaning up...'
                print
                os.system('rm -fr '+tmpdir2)

    if byjob == 0:
        print 'Number of complete samples:',itot, itot2
        print

        print 'Copy back to home base...'
        print
        os.system('cp -p '+tmpdir+'/'+MCname+sufMC+'_*.txt '+os.getcwd()+'/'+outdir+'/')

        if clean == 1:
            print 'Cleaning up...'
            print
            for i in range(1,njob+1):
                tmpdir2 = tmpdir+'/'+MCname+'_sampling.'+str(i)
                os.system('rm -fr '+tmpdir2)

    print 'Done !'

