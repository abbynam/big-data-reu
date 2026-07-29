# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 09:49:39 2018

@author: paul.maggi
"""

import sys, os
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
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
from detSimFunctRepo_UCT_3 import *
import numpy.random as random
import multiprocessing as mp
import matplotlib
import scipy.stats as stats

if __name__ == '__main__':
    tStart = time.time()
    
########################################
########### INPUT VARIABLES ############
########################################

    beamEne = 150 #MeV
    simProtons = 1e9 #total number of simulated protons
    DR = 20E3 #MU/min 
    MUdelivered = 50E3 # [FS] Still not sure what to set this to ...
    useMultiprocessing = False
    makeEachCrystalAMod = False #Do not use with multiprocessing! unless you have some sort of 64 core computational beast
    simType = 'ref'
    saveData = True
    resTime = 1.5E-6 #absolute value of active/coinicdence window 
    deadScl = 1 #this is a multiplicative scaling factor for dead time. A value of 1 = 100%. 0.5 = 50%
    rReduction = 1 #this is a dividing scaling factor for electron cloud radius. 1 means use full electron cloud. If you want to 
                   #effectively disable electron cloud-multipix scaling, divide by 1000 or some other large number. 

########################################
####### FRANK'S INPUT VARIABLES ########
########################################

    materialType = 'hdpe' # [FS] specifies material used for the source. Currently supports 'hdpe' and 'graphite'
    setup = 0 # [FS] 0: ortho setup 1: back-to-back with modules on XZ planes 2: alternate back-to-back (YZ plane orientation)
    dropTwoCrystals = True # [FS] sets number of crystals per module - True: 2 crystals arranged along z axis (two top crystals). False: 4 crystals in 2X2 square grid
    detectorWidth = 60 # [FS] Rough upper bound for the width of the sensor for each detector. Used for translations/rotations

########################################
############## LOAD DATA ###############
########################################

    # [FS] streamlined naming convention for annihilation/non-annihilation gammas
    filename = "Data/proton_HDPE_ortho_70_1e9"
    filename_annihil = filename + "_annihil.csv"
    filename_notAnnihil = filename + "_notAnnihil.csv"

    singlesANon511 = getSingles(filename_notAnnihil)
    doublesANon511 = getDoubles(filename_notAnnihil)
    triplesANon511 = getTriples(filename_notAnnihil)
    
    singlesA511 = getSingles(filename_annihil)
    random.shuffle(singlesA511)
    doublesA511 = getDoubles(filename_annihil)
    random.shuffle(doublesA511)
    triplesA511 = getTriples(filename_annihil)
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
    
    # [FS] counting the total number of intermod doubles and triples in the simulated dataset as a whole
    # -------------------------------------------------------------------------
    numSimIntermodDoubles = 0
    filt1 = doublesA511[:,2] < doublesA511[:,1]
    filt2 = doublesA511[:,6] >= doublesA511[:,5]
    intermodDoubles = doublesA511[filt1 & filt2]
    numSimIntermodDoubles += intermodDoubles.size
    
    filt1 = doublesA511[:,2] >= doublesA511[:,1]
    filt2 = doublesA511[:,6] < doublesA511[:,5]
    intermodDoubles = doublesA511[filt1 & filt2]
    numSimIntermodDoubles += intermodDoubles.size
    
    numSimIntermodTriples = 0
    filt1 = triplesA511[:,2] < triplesA511[:,1]
    filt2 = triplesA511[:,6] >= triplesA511[:,5]
    filt3 = triplesA511[:,10] >= triplesA511[:,9]
    intermodTriples = triplesA511[filt1 & (filt2 | filt3)]
    numSimIntermodTriples += intermodTriples.size
    
    filt1 = triplesA511[:,2] >= triplesA511[:,1]
    filt2 = triplesA511[:,6] < triplesA511[:,5]
    filt3 = triplesA511[:,10] < triplesA511[:,9]
    intermodTriples = triplesA511[filt1 & (filt2 | filt3)]
    numSimIntermodTriples += intermodTriples.size
    # -------------------------------------------------------------------------
    
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
    
    # [FS] counting the total number of intermod doubles and triples in the reduced dataset
    # -------------------------------------------------------------------------
    numRedIntermodDoubles = 0
    filt1 = doublesA511[:num511D,2] < doublesA511[:num511D,1]
    filt2 = doublesA511[:num511D,6] >= doublesA511[:num511D,5]
    intermodDoubles = doublesA511[:num511D,:][filt1 & filt2]
    numRedIntermodDoubles += intermodDoubles.size
    
    filt1 = doublesA511[:num511D,2] >= doublesA511[:num511D,1]
    filt2 = doublesA511[:num511D,6] < doublesA511[:num511D,5]
    intermodDoubles = doublesA511[:num511D,:][filt1 & filt2]
    numRedIntermodDoubles += intermodDoubles.size
    
    numRedIntermodTriples = 0
    filt1 = triplesA511[:num511T,2] < triplesA511[:num511T,1]
    filt2 = triplesA511[:num511T,6] >= triplesA511[:num511T,5]
    filt3 = triplesA511[:num511T,10] >= triplesA511[:num511T,9]
    intermodTriples = triplesA511[:num511T,:][filt1 & (filt2 | filt3)]
    numRedIntermodTriples += intermodTriples.size
    
    filt1 = triplesA511[:num511T,2] >= triplesA511[:num511T,1]
    filt2 = triplesA511[:num511T,6] < triplesA511[:num511T,5]
    filt3 = triplesA511[:num511T,10] < triplesA511[:num511T,9]
    intermodTriples = triplesA511[:num511T,:][filt1 & (filt2 | filt3)]
    numRedIntermodTriples += intermodTriples.size
    # -------------------------------------------------------------------------  
    # [FS] Here I am creating index lists of the same size as the full datasets
    # (singlesA, doublesA, triplesA), where a '1' at an index in the index list
    # corresponds to a 511 at the same index in the dataset, and '0' otherwise.
    # Notice that index lists are shuffled using the same random seed as their
    # corresponding dataset.
    # -------------------------------------------------------------------------
    singlesA = np.concatenate((singlesANon511,singlesA511[:num511S,:]))
    singlesA511Index = np.concatenate((np.zeros(singlesANon511.size),np.ones(num511S)))
    shuffle_in_unison(singlesA, singlesA511Index)
    
    doublesA = np.concatenate((doublesANon511,doublesA511[:num511D,:]))
    doublesA511Index = np.concatenate((np.zeros(doublesANon511.size),np.ones(num511D)))
    shuffle_in_unison(doublesA, doublesA511Index)
    
    triplesA = np.concatenate((triplesANon511,triplesA511[:num511T,:]))
    triplesA511Index = np.concatenate((np.zeros(triplesANon511.size),np.ones(num511T)))
    shuffle_in_unison(triplesA, triplesA511Index)
    # -------------------------------------------------------------------------    
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
    
    ########################################
    ######## SETUP FOR UCT DETECTOR ########
    ########################################
    
    
    # [FS] Translates and rotates the intersections with the second module to
    # lie on the same plane (parallel to XZ) as the first module, to the right
    # of the first module.
    # ------------------------------------------------------------------------- 
    if setup == 0:
        translate(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12],yT = -detectorWidth, filt_ = filter1)
        rotate(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12], angle = np.pi/2, filt_ = filter1)
    elif setup == 1:
        translate(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12],xT = detectorWidth, filt_ = filter1)
        rotate(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12], angle = np.pi, filt_ = filter1)
    elif setup == 2:
        translate(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12],yT = -detectorWidth)
        rotate(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12], angle = np.pi/2)
        rotateOther(singlesA, doublesA[:, :4], doublesA[:, 4:8], triplesA[:, :4], triplesA[:, 4:8], triplesA[:, 8:12], angle = -np.pi/2)
    # ------------------------------------------------------------------------- 
    # [FS] autoFullPixelateCommon(...) was simplified to handle the case of a
    # single module with two or four crystals.
    print('Starting pixel indexing')
    outs = autoFullPixelateCommon(dropTwoCrystals, singlesA,doublesA,triplesA)
    
    if makeEachCrystalAMod:
        outs = makeIndivModReadout(outs)
    xRef = outs[9]
    zRef = outs[10]
    tSamp = -np.log(1-random.rand(totEvents))/cRateMod # [FS] modifying to treat each module separately
    if DR == 1E-6:
        tSamp += 5 - tSamp
    print('Finished indexing')
    
    # [FS] I have ignored multiprocessing for now.    
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
        # [FS] Here I have added a lot of new parameters and outputs
        newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow,A511Counts,A511Filter, mod0count,mod1count =  moduleSimulation(singlesA,doublesA,triplesA,singlesA511Index,doublesA511Index,triplesA511Index,tSamp,nList,outs,resTime,deadScl,rReduction,dropTwoCrystals)
    
    tStop = time.time()
    print()
    print('It took %g seconds' % (tStop-tStart))
    
    # [FS] Here I am formatting the filters for the annihilation gamma stats
    # ------------------------------------------------------------------------- 
    A511singC = np.array(A511Filter[0])
    A511doubC = np.array(A511Filter[1])
    A511tripC = np.array(A511Filter[2])
    
    A511singTR = np.array(A511Filter[3])
    A511doubTR = np.array(A511Filter[4])
    A511tripTR = np.array(A511Filter[5])
    
    A511singD = np.array(A511Filter[6])
    A511doubD = np.array(A511Filter[7])
    A511tripD = np.array(A511Filter[8])

    A511singTD = np.array(A511Filter[9])
    A511doubTD = np.array(A511Filter[10])
    A511tripTD = np.array(A511Filter[11])
    # ------------------------------------------------------------------------- 
    
    # [FS] I removed the data sorting so the annihilation gamma labels would
    # correspond to the right scatter events. However, didn't seem to affect
    # the output anyway ...
    
    newSings = pd.DataFrame(data=newSings)
#    newSings.sort_values(4,inplace=True) 
    newSings = np.array(newSings)
    tS = newSings[:,-2].copy()
    peakTrackS = newSings[:,-1].copy()
    newSings = np.array(newSings[:,:4])
##    
    newDoubs = pd.DataFrame(data=newDoubs)
#    newDoubs.sort_values(8,inplace=True)
    newDoubs = np.array(newDoubs)
    tD = newDoubs[:,-4].copy()
    dType = newDoubs[:,-3].copy() #0 is true coincidence, 1 is false coincidence
    peakTrackD = newDoubs[:,-2:].copy()
    newDoubs = np.array(newDoubs[:,:8])
    eSum2 = newDoubs[:,0]+newDoubs[:,4]
#    
    #nts = newTrips.copy()
    newTrips = pd.DataFrame(data=newTrips)
#    newTrips.sort_values(12,inplace=True)
    newTrips = np.array(newTrips)
    tT = newTrips[:,-5].copy()
    tType = newTrips[:,-4].copy() #tType = 0: true coincidence. 1: false coinc; 2: doubToTrip, doub is first 2. 3: doubTotrip, single first
    peakTrackT = newTrips[:,-3:].copy()
    newTrips = np.array(newTrips[:,:12])
    eSum3 = newTrips[:,0]+newTrips[:,4]+newTrips[:,8]
    
    enes = np.linspace(0,5,501)
    
    #[FS] Translates and rotates the intersections with the second module to lie on the same plane (parallel to XZ) as the first module, to the right of the first module.
    # ------------------------------------------------------------------------- 
    if setup == 0:
        rotate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12], angle = -np.pi/2, filt_ = filter3)
        translate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12],yT = +detectorWidth, filt_ = filter3)
    elif setup == 1:
        rotate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12], angle = -np.pi, filt_ = filter3)
        translate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12],xT = -detectorWidth, filt_ = filter3)
    elif setup == 2:
        rotate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12], angle = -np.pi/2, filt_ = filter3)
        translate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12],yT = detectorWidth, filt_ = filter3)
        rotate(newSings, newDoubs[:, :4], newDoubs[:, 4:8], newTrips[:, :4], newTrips[:, 4:8], newTrips[:, 8:12], angle = -np.pi/2, filt_ = filter3)
    # ------------------------------------------------------------------------- 
    # [FS] a little scatter plot I used to check transformations
    #    plt.scatter(newTrips[:,1],newTrips[:,2])
    #%%
    print()
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

    # [FS] My output for the count statistics - see report for examples of the tables made from these
    # -------------------------------------------------------------------------
    print("-------------------------------------------------------------------")
    print(filename)
    print('Dose rate: %g kMU' % (DR/1000.0))
    print()
    print("TABLE 1")
    print('Singles per 10e7 p+ = %g' % (newSings.shape[0]/5))
    print('Doubles per 10e7 p+ = %g' % (newDoubs.shape[0]/5))
    print('Triples per 10e7 p+ = %g' % (newTrips.shape[0]/5))
    print('S:D Ratio: %g' % (newSings.shape[0]/newDoubs.shape[0]))
    print('S:T Ratio: %g' % (newSings.shape[0]/newTrips.shape[0]))
    print('Deadtime: %g %%' % (100*countOverflow[4]/(tS.max()*2))) # [FS] SOMETHING IS SUPER FISHY ABOUT THIS
    print()
    print("TABLE 2")
    print("False Doubles: %g %%" % (100*countOverflow[1]/newDoubs.shape[0]))
    print("False Triples: %g %%" % (100*countOverflow[2]/newTrips.shape[0]))
    print('D to T: %g %%' % (100*countOverflow[3]/newTrips.shape[0]))
    
    if countOverflow[6] != 0:
        print('Percent inter-module doubles detected: %g' % (countOverflow[6]/countOverflow[10]))
        print('Number inter-module doubles detected: %g' % countOverflow[6])
    else:
        print('No intermodule doubles detected.')
    if countOverflow[7] != 0:
        print('Percent inter-module triples detected: %g' % (countOverflow[7]/countOverflow[11]))
        print('Number inter-module triples detected: %g' % countOverflow[7])
    else:
        print('No intermodule triples detected.')
    if countOverflow[10] == 0:
        print('No ideal intermodule doubles.')
    if countOverflow[11] == 0:
        print('No ideal intermodule triples.')
    
    print("-------------------------------------------------------------------")
    print("511 INFORMATION")
    print()
    print("Total simulated 511s: %g" % (tot))
    print("Total simulated 511 singles: %g" % (nS) )
    print("Total simulated 511 doubles: %g" % (nD) )
    print("Total simulated 511 triples: %g" % nT )
    print("Total simulated intermod 511 doubles: %g" % numSimIntermodDoubles)
    print("Total simulated intermod 511 triples: %g" % numSimIntermodTriples)
    print()
    print("Total in-experiment 511s: %g" % np.round(tot511) )
    print("Total in-experiment 511 singles: %g" % num511S )
    print("Total in-experiment 511 doubles: %g" % num511D )
    print("Total in-experiment 511 triples: %g" % num511T )
    print("Total in-experiment intermod 511 doubles: %g" % numSimIntermodDoubles)
    print("Total in-experiment intermod 511 triples: %g" % numSimIntermodTriples)    
    print()
    print("Total contributing 511s: %g" % (A511Counts[0] + A511Counts[1] + A511Counts[2]) )
    print("Total contributing 511 singles: %g" % A511Counts[0] )
    print("Total contributing 511 doubles: %g" % A511Counts[1] )
    print("Total contributing 511 triples: %g" % A511Counts[2] )
    print()
    print("Total ideal 511 intermod doubles: %g" % A511Counts[3])
    print("Total ideal 511 intermod triples: %g" % A511Counts[4])
    print("Total detected 511 intermod doubles: %g" % A511Counts[5])
    print("Total detected 511 intermod triples: %g" % A511Counts[6])
    print()
    print("Total true recorded 511 singles %g" % A511Counts[7]) # "pure" 511 interactions recorded by the detector
    print("Total true recorded 511 doubles %g" % A511Counts[8])
    print("Total true recorded 511 triples %g" % A511Counts[9])
    print()
    print("Total detected 511 singles %g" % A511Counts[10])
    print("Total detected 511 doubles %g" % A511Counts[11])
    print("Total detected 511 triples %g" % A511Counts[12])
    print()
    print("Percentage true detected 511 singles %g" % (A511Counts[13]/A511Counts[10]))
    print("Percentage true detected 511 doubles %g" % (A511Counts[14]/A511Counts[11]))
    print("Percentage true detected 511 triples %g" % (A511Counts[15]/A511Counts[12]))
    print()
#    print("Total false 511 singles %g" % A511Counts[16])
#    print("Total false 511 doubles %g" % A511Counts[17])
#    print("Total false 511 triples %g" % A511Counts[18])
    # -------------------------------------------------------------------------
    # [FS] Lots of plots
    # -------------------------------------------------------------------------
    plt.rcParams["font.family"] = "serif"
    matplotlib.rcParams.update({"font.size": 14})
    
    # [FS] just adding a graph to visualise the energy spectrum of recorded measurements
    Ea = np.zeros(singlesA[:,0].size + doublesA[:,0].size + triplesA[:,0].size)
    Ea[:singlesA[:,0].size] = singlesA[:,0]
    Ea[singlesA[:,0].size:singlesA[:,0].size + doublesA[:,0].size] = doublesA[:,0] + doublesA[:,4]
    Ea[singlesA[:,0].size + doublesA[:,0].size:singlesA[:,0].size + doublesA[:,0].size + triplesA[:,0].size] = triplesA[:,0] + triplesA[:,4] + triplesA[:,8]   
    #hist(Ea, bins = 200, color = 'b', label = "Simulated (True) energy spectrum")
    Es = np.zeros(newSings[:,0].size + newDoubs[:,0].size + newTrips[:,0].size)
    EsSings = newSings[:,0]
    EsDoubs = newDoubs[:,0] + newDoubs[:,4]
    EsTrips = newTrips[:,0] + newTrips[:,4] + newTrips[:,8]   

    doubsfilter = dType.astype(int)
    tripsfilter = tType.astype(int)
    EsDoubsTrue = EsDoubs[doubsfilter==0]
    EsDoubsFalse = EsDoubs[doubsfilter==1]
    
    EsTripsTrue = EsTrips[tripsfilter==0]
    EsTripsFalse = EsTrips[tripsfilter==1]
    EsTripsDToT = EsTrips[tripsfilter>1]
    
    hist(EsDoubsTrue, bins = 200, color = 'orange', histtype = "step", label = "True")
    hist(EsDoubsFalse, bins = 200, color = 'blue', histtype = "step", label = "False")
    yscale("log")
    #density = stats.gaussian_kde(Es)
    #plot(np.sort(Es), density(np.sort(Es)), color = "orange", label = True)
    legend()
    xlim(0,5)
    xlabel("Energy (MeV)")
    ylabel("Relative intensity (counts)")
    show()
    
    hist(EsTripsTrue, bins = 200, color = 'green', histtype = "step", label = "True")
    hist(EsTripsFalse, bins = 200, color = 'blue', histtype = "step", label = "False")
    hist(EsTripsDToT, bins = 200, color = 'orange', histtype = "step", label = "D to T")
    yscale("log")
    #density = stats.gaussian_kde(Es)
    #plot(np.sort(Es), density(np.sort(Es)), color = "orange", label = True)
    legend()
    xlim(0,5)
    xlabel("Energy (MeV)")
    ylabel("Relative intensity (counts)")
    show()

    # [FS] Producing an event separation graph for doubles, as in Maggi's paper
    newDoubsTrue = newDoubs[dType.round() == 0]
    newDoubsFalse = newDoubs[dType.round() == 1]
    trueEventDist = np.sqrt((newDoubsTrue[:,1] - newDoubsTrue[:,5])**2 + (newDoubsTrue[:,2] - newDoubsTrue[:,6])**2 + (newDoubsTrue[:,3] - newDoubsTrue[:,7])**2)
    falseEventDist = np.sqrt((newDoubsFalse[:,1] - newDoubsFalse[:,5])**2 + (newDoubsFalse[:,2] - newDoubsFalse[:,6])**2 + (newDoubsFalse[:,3] - newDoubsFalse[:,7])**2)
    hist(trueEventDist, bins = 40, range= [0,40], color = "orange", alpha = 0.5, label = "True")
    xlim(0,40)
    hist(falseEventDist, bins = 40,range = [0,40], color = "blue", alpha = 0.5, label = "False")
    xlabel("Event Separation Distance [mm]")
    ylabel("Relative intensity (counts/mm)")
    legend()
    show()
    # [FS] Producing an event separation graph for triples, as in Maggi's paper
    newTripsTrue = newTrips[tType.round() == 0]
    newTripsFalse = newTrips[tType.round() == 1]
    newTripsDtoT = np.concatenate((newTrips[tType.round() == 2], newTrips[tType.round() == 3]), axis = 0) 
    trueEventDist = np.sqrt((newTripsTrue[:,1] - newTripsTrue[:,5])**2 + (newTripsTrue[:,2] - newTripsTrue[:,6])**2 + (newTripsTrue[:,3] - newTripsTrue[:,7])**2)
    falseEventDist = np.sqrt((newTripsFalse[:,1] - newTripsFalse[:,5])**2 + (newTripsFalse[:,2] - newTripsFalse[:,6])**2 + (newTripsFalse[:,3] - newTripsFalse[:,7])**2)
    dToTEventDist = np.sqrt((newTripsDtoT[:,1] - newTripsDtoT[:,5])**2 + (newTripsDtoT[:,2] - newTripsDtoT[:,6])**2 + (newTripsDtoT[:,3] - newTripsDtoT[:,7])**2)
    hist(trueEventDist,  bins = 40, range= [0,40],color = "green", alpha = 0.5, label = "True")
    hist(falseEventDist,  bins = 40, range= [0,40], color = "blue", alpha = 0.5, label = "False")
    hist(dToTEventDist,  bins = 40, range= [0,40], color = "orange", alpha = 0.5, label = "D to T")
    xlim(0,40)
    xlabel("Event Separation Distance [mm]")
    ylabel("Relative intensity (counts/mm)")
    legend()
    show()    
    # [FS] Producing a time difference graph for all data, as in Maggi's paper
    t1 = np.concatenate((tS[filter1(newSings)], tD[filter1(newDoubs[:,0:4])], tT[filter1(newTrips[:,0:4])]))
    t2 = np.concatenate((tS[filter1(newSings)], tD[filter1(newDoubs[:,0:4])], tT[filter1(newTrips[:,0:4])]))
    t1 = np.sort(t1)
    t2 = np.sort(t2)
    tDiff1 = np.diff(t1)
    tDiff2 = np.diff(t2)
    tDiff = np.concatenate((tDiff1, tDiff2))
    #t = np.concatenate((tS,tD,tT))
    #t = np.sort(t)
    #tDiff = np.diff(t)
    xscale("log")
    
    hist(tDiff, bins=numpy.logspace(start=np.log10(1e-5), stop=np.log10(numpy.max(tDiff)), num=40), color = "red", alpha = 0.5, label = "Simulated")
    plt.xscale("log")
    xlabel("Time Difference ($s$)")
    ylabel("Relative Intensity ($s^{-1}$)")
    legend()
    show()
    
    if makeEachCrystalAMod:
        print('Dead time percentage: %g' % (countOverflow[4]/(tS.max()*64)))
    else:
        print()
        #print('Dead time percentage: %g' % (countOverflow[4]/(tS.max()*16)))
        
        #%%
    if saveData:
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_singles.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),newSings)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_doubles.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),newDoubs)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_triples.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),newTrips)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_dType.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),dType)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tType.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),tType)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tT.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),tT)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),tD)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_tS.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),tS)
        # [FS] Exporting my 511 filters for the data ...
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511singC.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511singC)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511doubC.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511doubC)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511tripC.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511tripC)
        
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511singTR.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511singTR)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511doubTR.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511doubTR)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511tripTR.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511tripTR)
        
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511singD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511singD)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511doubD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511doubD)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511tripD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511tripD)
        
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511singTD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511singTD)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511doubTD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511doubTD)
        np.savetxt('Output/%iGEvent_%iMeV_%ikMUmin_%s_res%s_deadScl%s_511tripTD.txt' % (simProtons/1.0E9,beamEne,DR/1E3,simType,resTime,deadScl),A511tripTD)

#    