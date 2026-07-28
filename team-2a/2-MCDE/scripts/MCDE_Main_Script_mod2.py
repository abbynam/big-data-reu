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
import pandas as pd
from openpyxl import load_workbook
from itertools import compress
import itertools
import bisect
from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse as sprs
from detSimFunctRepo_mod import *
import numpy.random as random
import multiprocessing as mp
from pathlib import Path

def main(csv_511: str, csv_PG: str):
    tStart = time.time()
    
########################################
########### INPUT VARIABLES ############
########################################
    
    beamEne = 150 #MeV
    simProtons = 1E9 #total number of simulated protons
    DR = 180E3 #MU/min
    MUdelivered = 50E3
    useMultiprocessing = False
    makeEachCrystalAMod = False #Do not use with multiprocessing! unless you have some sort of 64 core computational beast
    simType = 'ref'  
    saveData = True
    resTime = 1.5E-6 #abosulute value of active/coinicdence window 
    deadScl = 1 #this is a multiplicitive scaling factor for dead time. A value of 1 = 100%. 0.5 = 50%
    rReduction = 1 #this is a dividing scaling factor for electron cloud radius. 1 means use full electron cloud. If you want to 
                   #effectively disable electron cloud-multipix scaling, divide by 1000 or some other large number. 
########################################    

########################################
############## LOAD DATA ###############
########################################
    
    # Script directory
    base_dir = Path(__file__).resolve().parent
    csv_dir = (base_dir/"../csv_in")

    yes511file = csv_dir/csv_511
    non511file = csv_dir/csv_PG
    
    # non511file = '21-item21_511(in).csv'
    # yes511file = '21-item21_PG_cat(in).csv'
    
    singlesANon511 = getSingles(non511file)
    doublesANon511 = getDoubles(non511file)
    triplesANon511 = getTriples(non511file)

    singlesA511 = getSingles(yes511file)
    random.shuffle(singlesA511)

    doublesA511 = getDoubles(yes511file)
    random.shuffle(doublesA511)

    triplesA511 = getTriples(yes511file)
    random.shuffle(triplesA511)
        
########################################
#### That's it! You're good to go! #####
########################################
    
    #%%
    
    np.random.seed()
    mcThresh = 0.005
    maxEThresh = 5    
    
    singlesANon511 = singlesANon511[singlesANon511[:,0]>mcThresh,:].copy()
    doublesANon511 = doublesANon511[np.logical_and(doublesANon511[:,0]>mcThresh,doublesANon511[:,4]>mcThresh),:].copy()
    triplesANon511 = triplesANon511[np.logical_and(triplesANon511[:,0]>mcThresh,triplesANon511[:,4]>mcThresh),:].copy()
    triplesANon511 = triplesANon511[triplesANon511[:,8]>mcThresh,:].copy()
    
    singlesANon511 = singlesANon511[singlesANon511[:,0]<maxEThresh].copy()
    doublesANon511 = doublesANon511[np.logical_and(doublesANon511[:,0]<maxEThresh,doublesANon511[:,4]<maxEThresh),:].copy()
    triplesANon511 = triplesANon511[np.logical_and(triplesANon511[:,0]<maxEThresh,triplesANon511[:,4]<maxEThresh),:].copy()
    triplesANon511 = triplesANon511[triplesANon511[:,8]<maxEThresh,:].copy()
    
    numNon511 = singlesANon511.shape[0]+doublesANon511.shape[0]+triplesANon511.shape[0]

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
    
    pPerS = DRtoPperS(beamEne,DR)
    irradTime = 60*MUdelivered/DR #in s
    
    tot511 = final511Count(beamEne,pPerS,irradTime,simProtons,num511=tot) #for 20k
    cRateMod = 0.4*(tot511 + numNon511) * pPerS / simProtons
    

    num511S = np.round(tot511*nS/tot).astype(int)
    num511D = np.round(tot511*nD/tot).astype(int)
    num511T = np.round(tot511*nT/tot).astype(int)
    
    singlesA = np.concatenate((singlesANon511,singlesA511[:num511S,:]))
    random.shuffle(singlesA)
    doublesA = np.concatenate((doublesANon511,doublesA511[:num511D,:]))
    random.shuffle(doublesA)
    triplesA = np.concatenate((triplesANon511,triplesA511[:num511T,:]))
    random.shuffle(triplesA)
    
    
    # FIX (Increased separation distance between front and back crystals)
    singlesA[:,2][singlesA[:,2]<-230] = singlesA[:,2][singlesA[:,2]<-230]-100
    doublesA[:,2][doublesA[:,2]<-230] = doublesA[:,2][doublesA[:,2]<-230]-100
    doublesA[:,6][doublesA[:,6]<-230] = doublesA[:,6][doublesA[:,6]<-230]-100
    triplesA[:,2][triplesA[:,2]<-230] = triplesA[:,2][triplesA[:,2]<-230]-100
    triplesA[:,6][triplesA[:,6]<-230] = triplesA[:,6][triplesA[:,6]<-230]-100
    triplesA[:,10][triplesA[:,10]<-230] = triplesA[:,10][triplesA[:,10]<-230]-100

    
    numSing = singlesA.shape[0]
    numDoub = doublesA.shape[0]
    numTrip = triplesA.shape[0]
    totEvents = numSing + numDoub + numTrip
    
    singL = np.zeros((numSing,))
    doubL = np.zeros((numDoub,))+1
    tripL = np.zeros((numTrip,))+2
    nList = np.concatenate((singL,doubL))
    nList = np.concatenate((nList,tripL))
    random.shuffle(nList)
    

    print('Starting pixel indexing')
    outs = autoFullPixelateCommon(singlesA,doublesA,triplesA)
    if makeEachCrystalAMod:
        outs = makeIndivModReadout(outs)
    xRef = outs[9]
    zRef = outs[10]
    tSamp = -np.log(1-random.rand(totEvents))/cRateMod
    if DR == 1E-6:
        tSamp += 5 - tSamp
    print('Finished indexing')
    
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
    
#    
    newSings = pd.DataFrame(data=newSings)
    newSings.sort_values(4,inplace=True)
    newSings = np.array(newSings)
    tS = newSings[:,-2].copy()
    peakTrackS = newSings[:,-1].copy()
    newSings = np.array(newSings[:,:4])
##    
    newDoubs = pd.DataFrame(data=newDoubs)
    newDoubs.sort_values(8,inplace=True)
    newDoubs = np.array(newDoubs)
    tD = newDoubs[:,-4].copy()
    dType = newDoubs[:,-3].copy() #0 is true coincidence, 1 is false coincidence
    peakTrackD = newDoubs[:,-2:].copy()
    newDoubs = np.array(newDoubs[:,:8])
    eSum2 = newDoubs[:,0]+newDoubs[:,4]
#    
    #nts = newTrips.copy()
    newTrips = pd.DataFrame(data=newTrips)
    newTrips.sort_values(12,inplace=True)
    newTrips = np.array(newTrips)
    tT = newTrips[:,-5].copy()
    tType = newTrips[:,-4].copy() #tType = 0: true coincidence. 1: false coinc; 2: doubToTrip, doub is first 2. 3: doubTotrip, single first
    peakTrackT = newTrips[:,-3:].copy()
    newTrips = np.array(newTrips[:,:12])
    eSum3 = newTrips[:,0]+newTrips[:,4]+newTrips[:,8]
    
    enes = np.linspace(0,5,501)
    
    #%%
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
        
        #%%
    if saveData:
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_singles.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),newSings)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_doubles.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),newDoubs)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_triples.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),newTrips)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_dType.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),dType)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tType.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),tType)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tT.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),tT)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tD.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),tD)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tS.txt' % (simProtons/1E9,beamEne,DR/1E3,simType,resTime,deadScl),tS)
#    

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} file_511.csv file_PG.csv")
        sys.exit(1)

    csv_511 = sys.argv[1]
    csv_PG = sys.argv[2]

    main(csv_511, csv_PG)
    
