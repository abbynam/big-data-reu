# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 09:49:39 2018

@author: paul.maggi
"""

import sys, os
import numpy as np
from matplotlib.pyplot import *
import matplotlib
from itertools import compress
import math
import timeit
import pandas as pdsingle
from openpyxl import load_workbook
from itertools import compress
import itertools
import bisect
from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse as sprs
from detSimFunctRepoErgInit import *
import numpy.random as random
import multiprocessing as mp
from pathlib import Path

if __name__ == '__main__':
    tStart = time.time()

    # python running dir: scripts/
    base_dir = Path(__file__).resolve().parent
  
#  #%%   
########################################
########### INPUT VARIABLES ############
########################################
    
    beamEne = 199 #MeV
    simProtons = 3.6E9 #total number of simulated protons
    DR = 100 # 0.001E3 #MU/min  20 kMU/min, 60 kMU/min, 100 kMU/min, 140 kMU/min, 180 kMU/min
    MUdelivered = 50E3
    useMultiprocessing = False
    makeEachCrystalAMod = False #Do not use with multiprocessing! unless you have some sort of 64 core computational beast
    simNum = 5
    simType = 'pjin_100DR_spot_3e9'  
    saveData = True
    saveLocation = 'saveMCDE'
    resTime = 1.5E-6 #abosulute value of active/coinicdence window 
    deadScl = 1.0 #this is a multiplicitive scaling factor for dead time. A value of 1 = 100%. 0.5 = 50%
    rReduction = 1 #this is a dividing scaling factor for electron cloud radius. 1 means use full electron cloud. If you want to 
                   #effectively disable electron cloud-multipix scaling, divide by 1000 or some other large number. 
########################################    

########################################
############## LOAD DATA ###############
########################################

    csv_dir = (base_dir/"../csv_in")

    if True:
        print("\nLoading 511")
        s511_path = '/umbc/rs/cybertrn/reu2026/team2/research/2-MCDE/csv_in/1-item1_PG_cat.csv'
        singlesANon511 = getSinglesErgInit(s511_path)
        print(f"\t{s511_path}")
        d511_path = '/umbc/rs/cybertrn/reu2026/team2/research/2-MCDE/csv_in/1-item1_PG_cat.csv'
        doublesANon511 = getDoublesErgInit(d511_path)
        print(f"\t{d511_path}")
        t511_path = '/umbc/rs/cybertrn/reu2026/team2/research/2-MCDE/csv_in/1-item1_PG_cat.csv'        
        triplesANon511 = getTriplesErgInit(t511_path)
        print(f"\t{t511_path}")
        print("\nLoading PG")
        sPG_path  = '/umbc/rs/cybertrn/reu2026/team2/research/2-MCDE/csv_in/1-item1_511.csv'
        singlesA511 = getSinglesErgInit(sPG_path)
        print(f"\t{sPG_path}")
        dPG_path  = '/umbc/rs/cybertrn/reu2026/team2/research/2-MCDE/csv_in/1-item1_511.csv'
        doublesA511 = getDoublesErgInit(dPG_path)
        print(f"\t{dPG_path}")
        tPG_path  = '/umbc/rs/cybertrn/reu2026/team2/research/2-MCDE/csv_in/1-item1_511.csv'
        triplesA511 = getTriplesErgInit(tPG_path)
        print(f"\t{tPG_path}")

    if False:
        print("\nLoading 511")
        s511_path = csv_dir/'compiled_PG.csv'
        singlesANon511 = getSinglesErgInit(s511_path)
        print(f"\t{s511_path}")
        d511_path = csv_dir/'compiled_PG.csv'
        doublesANon511 = getDoublesErgInit(d511_path)
        print(f"\t{d511_path}")
        t511_path = csv_dir/'compiled_PG.csv'
        triplesANon511 = getTriplesErgInit(t511_path)
        print(f"\t{t511_path}")
        print("\nLoading PG")
        sPG_path = csv_dir/'compiled_511.csv'
        singlesA511 = getSinglesErgInit(sPG_path)
        print(f"\t{sPG_path}")
        dPG_path = csv_dir/'compiled_511.csv'
        doublesA511 = getDoublesErgInit(dPG_path)
        print(f"\t{dPG_path}")
        tPG_path = csv_dir/'compiled_511.csv'
        triplesA511 = getTriplesErgInit(tPG_path)
        print(f"\t{tPG_path}")

    random.shuffle(singlesA511)
    random.shuffle(doublesA511)
    random.shuffle(triplesA511)
    
    print('\nSim Number: %i' % simNum)    
    print('Beam energy (MeV): %g' % beamEne)
    print('number of protons: %.2E' % simProtons)
    print('Dose Rate: (Mu/min): %.2E' % DR)
    print('MU Delivered: %.2E' % MUdelivered)


########################################
#### That's it! You're good to go! #####
########################################

    
    np.random.seed()
    mcThresh = 0.01
    maxEThresh = 5.0 

    print('Min proton Energy (MeV): %g' % mcThresh)
    print('Max proton Energy (MeV): %g' % maxEThresh)   


#    #%%
    
    singlesANon511 = singlesANon511[singlesANon511[:,0]>mcThresh,:].copy()
    doublesANon511 = doublesANon511[np.logical_and(doublesANon511[:,0]>mcThresh,doublesANon511[:,4]>mcThresh),:].copy()
    triplesANon511 = triplesANon511[np.logical_and(triplesANon511[:,0]>mcThresh,triplesANon511[:,4]>mcThresh),:].copy()
    triplesANon511 = triplesANon511[triplesANon511[:,8]>mcThresh,:].copy()
    
    singlesANon511 = singlesANon511[singlesANon511[:,0]<maxEThresh].copy()
    doublesANon511 = doublesANon511[np.logical_and(doublesANon511[:,0]<maxEThresh,doublesANon511[:,4]<maxEThresh),:].copy()
    triplesANon511 = triplesANon511[np.logical_and(triplesANon511[:,0]<maxEThresh,triplesANon511[:,4]<maxEThresh),:].copy()
    triplesANon511 = triplesANon511[triplesANon511[:,8]<maxEThresh,:].copy()

    nnS = singlesANon511.shape[0]
    nnD = doublesANon511.shape[0]
    nnT = triplesANon511.shape[0]
    numNon511 = nnS+nnD+nnT

    # calculate the eSum 13th column if needed
    # this makes MCDE 12 or 13 column agnostic
    
#    if singlesANon511.shape[1] == 12:
#        print("Adding 13th column to PG...")
#        eSum = singlesANon511[:,0]
#        singlesANon511 = np.hstack((singlesANon511, eSum[:,None]))
#    if doublesANon511.shape[1] == 12:
#        eSum = doublesANon511[:,0] + doublesANon511[:,4]
#        doublesANon511 = np.hstack((doublesANon511, eSum[:,None]))
#    if triplesANon511.shape[1] == 12:
#        eSum = triplesANon511[:,0] + triplesANon511[:,4] + triplesANon511[:,8]
#        triplesANon511 = np.hstack((triplesANon511, eSum[:,None]))


    singlesA511 = singlesA511[singlesA511[:,0]>mcThresh,:].copy()
    doublesA511 = doublesA511[np.logical_and(doublesA511[:,0]>mcThresh,doublesA511[:,4]>mcThresh),:].copy()
    triplesA511 = triplesA511[np.logical_and(triplesA511[:,0]>mcThresh,triplesA511[:,4]>mcThresh),:].copy()
    triplesA511 = triplesA511[triplesA511[:,8]>mcThresh,:].copy()
    
    singlesA511 = singlesA511[singlesA511[:,0]<maxEThresh].copy()
    doublesA511 = doublesA511[np.logical_and(doublesA511[:,0]<maxEThresh,doublesA511[:,4]<maxEThresh),:].copy()
    triplesA511 = triplesA511[np.logical_and(triplesA511[:,0]<maxEThresh,triplesA511[:,4]<maxEThresh),:].copy()
    triplesA511 = triplesA511[triplesA511[:,8]<maxEThresh,:].copy()
    
    nS = singlesA511.shape[0]
    nD = doublesA511.shape[0]
    nT = triplesA511.shape[0]
    tot = nS+nD+nT

    # calculate the eSum 13th column if needed
    # this makes MCDE 12 or 13 column agnostic
    
#    if singlesA511.shape[1] == 12:
#        print("Adding 13th column to 511...")
#        eSum = singlesA511[:,0]
#        singlesA511 = np.hstack((singlesA511, eSum[:,None]))
#    if doublesA511.shape[1] == 12:
#        eSum = doublesA511[:,0] + doublesA511[:,4]
#        doublesA511 = np.hstack((doublesA511, eSum[:,None]))
#    if triplesA511.shape[1] == 12:
#        eSum = triplesA511[:,0] + triplesA511[:,4] + triplesA511[:,8]
#        triplesA511 = np.hstack((triplesA511, eSum[:,None]))

    
    pPerS = DRtoPperS(beamEne,DR)
    irradTime = 60*MUdelivered/DR #in s
 
    print('\nprotons per second: %g' % pPerS)
    print('irradiation time: %g\n' % irradTime)
    
    tot511 = final511Count(beamEne,pPerS,irradTime,simProtons,num511=tot) #for 20k
    cRateMod = 0.4*(tot511 + numNon511) * pPerS / simProtons
    

    num511S = np.round(tot511*nS/tot).astype(int)
    num511D = np.round(tot511*nD/tot).astype(int)
    num511T = np.round(tot511*nT/tot).astype(int)
    
    print('\nTotal 511 included: %g' % tot511)
    print('single 511 included: %g' % num511S)
    print('double 511 included: %g' % num511D)
    print('triple 511 included: %g\n' % num511T)    
    
    singlesA = np.concatenate((singlesANon511,singlesA511[:num511S,:]))
    random.shuffle(singlesA)
    doublesA = np.concatenate((doublesANon511,doublesA511[:num511D,:]))
    random.shuffle(doublesA)
    triplesA = np.concatenate((triplesANon511,triplesA511[:num511T,:]))
    random.shuffle(triplesA)
    
    numSing = singlesA.shape[0]
    numDoub = doublesA.shape[0]
    numTrip = triplesA.shape[0]
    totEvents = numSing + numDoub + numTrip
    #totEvents = numDoub
    singL = np.zeros((numSing,))
    doubL = np.zeros((numDoub,))+1
    tripL = np.zeros((numTrip,))+2
    #nList = np.concatenate((singL,doubL))
    nList = np.concatenate((doubL,singL))
    nList = np.concatenate((nList,tripL))

    random.shuffle(nList)
    




 
#        #%%    
    
    print('Starting pixel indexing')
    print(f"\nsinglesA shape: {singlesA.shape}")
    print(f"doublesA shape: {doublesA.shape}")
    print(f"triplesA shape: {triplesA.shape}\n")

    outs = autoFullPixelateCommon(singlesA,doublesA,triplesA)
    if makeEachCrystalAMod:
        outs = makeIndivModReadout(outs)
    xRef = outs[9]
    zRef = outs[10]
    tSamp = -np.log(1-random.rand(totEvents))/cRateMod
    if DR == 1E-6:
        tSamp += 5 - tSamp
    print('Finished indexing')

#        #%%    
    
    if useMultiprocessing:
        if makeEachCrystalAMod:
            multiplicity = 64
        else:
            multiplicity = 16
    
        pr = []
        q = []
        z = []
        
        newSings = []
        newDoubs = []
        newTrips = []
        coincCount = 0
        deadCount = 0
        countOverflow = np.zeros((14,))
        
        for iii in range(multiplicity):
            q.append(mp.SimpleQueue())
            newSinglesA,newDoublesA,newTriplesA,newOuts,newtSamp,newNList = makeSingleMod(singlesA,doublesA,triplesA,outs,tSamp,nList,iii)
            pr.append(mp.Process(target=moduleSimulation,args=(newSinglesA,newDoublesA,newTriplesA,newtSamp,newNList,newOuts,resTime,deadScl,rReduction,q[iii])))
                #
            pr[iii].start()
            print(iii)
            
        for iii in range(multiplicity):
            z.append(q[iii].get())
            pr[iii].join()
    
        for iii in range(multiplicity):
            newSings.append(z[iii][0])
            newDoubs.append(z[iii][1])
            newTrips.append(z[iii][2])
            coincCount += z[iii][3]
            deadCount += z[iii][4]
            countOverflow += z[iii][5]
            
        newSings = list(itertools.chain(*newSings))
        newDoubs = list(itertools.chain(*newDoubs))
        newTrips = list(itertools.chain(*newTrips))
    else:
        newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow = moduleSimulation(singlesA,doublesA,triplesA,tSamp,nList,outs,resTime,deadScl,rReduction)

    tStop = time.time()
    print('It took %g seconds' % (tStop-tStart))

#    #%%    



    
#    

    newSings = pd.DataFrame(data=newSings)
    newSings.sort_values(5,inplace=True)      #ordering sngles according to timestamp
    newSings = np.array(newSings)
    tS = newSings[:,-2].copy()
    newSingsErgInit = newSings[:,-3].copy()
    peakTrackS = newSings[:,-1].copy()
    newSings = np.array(newSings[:,:4])
##    

    newDoubs = pd.DataFrame(data=newDoubs)
    newDoubs.sort_values(9,inplace=True)      #ordering doubles according to timestamp
    newDoubs = np.array(newDoubs)
    tD = newDoubs[:,-4].copy()
    newDoubsErgInit = newDoubs[:,-5].copy()    
    dType = newDoubs[:,-3].copy() #0 is true coincidence, 1 is false coincidence
    peakTrackD = newDoubs[:,-2:].copy()
    newDoubs = np.array(newDoubs[:,:8])
    eSum2 = newDoubs[:,0]+newDoubs[:,4]
#    
    #nts = newTrips.copy()
    newTrips = pd.DataFrame(data=newTrips)
    newTrips.sort_values(13,inplace=True)         #ordering doubles according to timestamp
    newTrips = np.array(newTrips)
    tT = newTrips[:,-5].copy()
    newTripsErgInit = newTrips[:,-6].copy()  
    tType = newTrips[:,-4].copy() #tType = 0: true coincidence. 1: false coinc; 2: doubToTrip, doub is first 2. 3: doubTotrip, single first
    peakTrackT = newTrips[:,-3:].copy()
    newTrips = np.array(newTrips[:,:12])
    eSum3 = newTrips[:,0]+newTrips[:,4]+newTrips[:,8]
    
    enes = np.linspace(0,5,501)
    
#    #%%
    print('Final pixel coincidence count was: %g' % coincCount)
    print('Final missed counts was: %g' % deadCount)
    print('Final count overflow was: %g' % countOverflow[0])
    print('Final false double was: %g' % countOverflow[1])
    print('Final false triple was: %g' % countOverflow[2])
    print('Final 2->3 was: %g' % countOverflow[3])
    print('Total Singles: %g' % newSings.shape[0])
    sInRange = np.logical_and(newSings[:,0]>0.05,newSings[:,0]<2.7).sum()
    print('Total Singles in range: %g' % sInRange)
    print('Total Doubles: %g' % newDoubs.shape[0])
    dInRange = np.logical_and(newDoubs[:,0]>0.05,newDoubs[:,0]<2.7)
    dInRange = np.logical_and(dInRange,newDoubs[:,4]<2.7)
    dInRange = np.logical_and(newDoubs[:,4]>0.05,dInRange).sum()
    print('Total Doubles in range: %g' % dInRange)
    print('Total Triples: %g' % newTrips.shape[0])
    tInRange = np.logical_and(newTrips[:,0]>0.05,newTrips[:,0]<2.7)
    tnRange = np.logical_and(tInRange,newTrips[:,4]<2.7)
    tInRange = np.logical_and(newTrips[:,4]>0.05,tInRange)
    tnRange = np.logical_and(tInRange,newTrips[:,8]<2.7)
    tInRange = np.logical_and(newTrips[:,8]>0.05,tInRange).sum()
    print('Total Triples in range: %g' % tInRange)
    print('Total Events: %g' % (newSings.shape[0]+newDoubs.shape[0]+newTrips.shape[0]))
    print('S:D Ratio: %g' % (newSings.shape[0]/newDoubs.shape[0]))
    print('S:T Ratio: %g' % (newSings.shape[0]/newTrips.shape[0]))
    print('Percent pix > 3: %g' % (countOverflow[0]/(countOverflow[0]+newSings.shape[0]+newDoubs.shape[0]+newTrips.shape[0])))

    print('Percent inter-module doubles detected: %g' % (countOverflow[6]/countOverflow[10]))
    print('Number inter-module doubles detected: %g' % countOverflow[6])
    print('Percent inter-module triples detected: %g' % (countOverflow[7]/countOverflow[11]))
    print('Number inter-module triples detected: %g' % countOverflow[7])
    
    if makeEachCrystalAMod:
        print('Dead time percentage: %g' % (countOverflow[4]/(tS.max()*64)))
    else:
        print('Dead time percentage: %g' % (countOverflow[4]/(tS.max()*16)))
        
#       #%%
    if saveData:
        os.makedirs(saveLocation, exist_ok=True)
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_singles.txt' % (beamEne,DR/1E3,simType),newSings,delimiter=',')
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_doubles.txt' % (beamEne,DR/1E3,simType),newDoubs,delimiter=',')
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_triples.txt' % (beamEne,DR/1E3,simType),newTrips,delimiter=',')
        
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_singlesInitErg.txt' %(beamEne,DR/1E3,simType),newSingsErgInit)        
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_doublesInitErg.txt' % (beamEne,DR/1E3,simType),newDoubsErgInit)  
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_triplesInitErg.txt' % (beamEne,DR/1E3,simType),newTripsErgInit)  

        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_dType.txt' % (beamEne,DR/1E3,simType,),dType)
        np.savetxt(saveLocation+'/%iMeV_%ikMUmin_%s_tType.txt' % (beamEne,DR/1E3,simType,),tType)
        #np.savetxt(saveLocation+'/%i_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tT.txt' % (simNum,beamEne,DR/1E3,simType,resTime,deadScl),tT)
        #np.savetxt(saveLocation+'/%i_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tD.txt' % (simNum,beamEne,DR/1E3,simType,resTime,deadScl),tD)
        #np.savetxt(saveLocation+'/%i_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tS.txt' % (simNum,beamEne,DR/1E3,simType,resTime,deadScl),tS)
#    


   #%%

  






