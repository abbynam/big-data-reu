import sys, os, subprocess, string
import numpy, scipy, scipy.interpolate
import numpy as np
import matplotlib
import pandas as pd
from openpyxl import load_workbook
from itertools import compress
import numpy.random as random
import copy
from scipy import fftpack, signal, ndimage, interpolate
import time
import scipy.sparse as sprs
import bisect
from pathlib import Path

cols = range(8)

class pixel:
    '''
    The pixel is the basic unit of an interaction. It stores all relevant information about the position, energy, and temporal evolution of an event.
    '''
    def __init__(self,module=0,crystal=0,xInd=0,zInd=0,edep=0,rawX=0,rawY=0,rawZ=0,ergInit=0,time=0,peakTrack=0,pulseShape=np.ones((1201)),eScl=1,eOut=0,pulseLength=5E-6,tShapeThresh=0.001,trigThresh=0.05,resTime = 1.5E-6):
        self.module = int(module)
        self.crystal = crystal
        self.xInd = xInd
        self.zInd = zInd
        self.edep = edep
        self.rawX = rawX
        self.rawY = rawY
        self.rawZ = rawZ
        self.ergInit = ergInit
        self.time = time
        self.peakTrack = 0
        self.eScl=eScl
        eOut = edep*eScl
        self.eOut=eOut
        self.pulseShape = pulseShape*edep*eScl
        tPulse = np.linspace(time,time+pulseLength,len(pulseShape))
        self.tPulse = tPulse
        self.tMax = tPulse[pulseShape==pulseShape.max()][0]
        if self.pulseShape.max() >= trigThresh:
            self.tTrig = tPulse[np.logical_and(tPulse<=self.tMax,pulseShape>=trigThresh)].min()
        else:
            self.tTrig = time
        if self.pulseShape.max()>= tShapeThresh and self.pulseShape[-1]< tShapeThresh:
            self.tStop = self.tTrig + resTime
        else:
            self.tStop = self.tMax

        
        
    def __len__(self):
        return 1
       
class module:
    '''
    The module class contains module status and timing information.
    '''
    def __init__(self,status=0,lastTime=0,nextTime=float('inf'),trigThresh=0.05):
        self.status = status #statuses are 0 = ready, 1 = triggered, 2 = dead
        self.lastTime = lastTime
        self.nextTime = nextTime
        self.trigThresh = trigThresh



def getDataPix( pixList,dataToGet):
    '''
    outData = getDataPix(pixList,dataToGet)
    This function take in a list of pixel objects and returns the requested atribute.
    Input:
        pixList : a list of pixel objects
        dataToGet : either a full string or corresponding letter of the data to get
    Output:
        outData : an array contatining the requested data for each item in the pixel list
    '''
    if (dataToGet == 'module') or (dataToGet == 'm'):
        outData = [temp.module for temp in pixList]
        
    elif (dataToGet == 'crystal') or (dataToGet == 'c'):
        outData = [temp.crystal for temp in pixList]
        
    elif (dataToGet == 'xInd') or (dataToGet == 'x'):
        outData = [temp.xInd for temp in pixList]
        
    elif (dataToGet == 'zInd') or (dataToGet == 'z'):
        outData = [temp.zInd for temp in pixList]
        
    elif (dataToGet == 'edep') or (dataToGet == 'e'):
        outData = [temp.edep for temp in pixList]
        
    elif (dataToGet == 'time') or (dataToGet == 't'):
        outData = [temp.time for temp in pixList]
        
    elif (dataToGet == 'ergInit') or (dataToGet == 'eo'):
        outData = [temp.ergInit for temp in pixList]
    
    elif (dataToGet == 'peakTrack') or (dataToGet == 'p'):
        outData = [temp.peakTrack for temp in pixList]
        
    elif (dataToGet == 'eOut') or ('o'):
        outData = [temp.eOut for temp in pixList]
        
    elif (dataToGet == 'eScl') or ('s'):
        outData = [temp.eScl for temp in pixList]
        
        
    elif (dataToGet == 'tPulse') or (dataToGet == 'u'):
        outData = [temp.tPulse.mean() for temp in pixList]
        
    else:
        outData = []

    return np.array(outData)

def getDataMod(modList,dataToGet):
    '''
    outData = getDataMod(modList,dataToGet)
    This function take in a list of module objects and returns the requested atribute.
    Input:
        modList : a list of module objects
        dataToGet : either a full string or corresponding letter of the data to get
    Output:
        outData : an array contatining the requested data for each item in the module list
    '''
    if (dataToGet == 'status') or (dataToGet == 's'):
        outData = [temp.status for temp in modList]
        
    elif (dataToGet == 'processed') or (dataToGet == 'p'):
        outData = [temp.processed for temp in modList]
        
    elif (dataToGet == 'lastTime') or (dataToGet == 'l'):
        outData = [temp.lastTime for temp in modList]
        
    elif (dataToGet == 'nextTime') or (dataToGet == 'n'):
        outData = [temp.nextTime for temp in modList]
        
    else:
        outData = []
        
    return np.array(outData)
    
def DRtoPperS(beamEne,DR):
    '''
    ionsPerS = DRtoPperS(beamEne,DR)
    This function calculates the number of protons output for a single pencil beam per second based on beam energy and dose rate.
    Input:
        beamEne : nominal beam energy in MeV
        DR : dose rate in MU/min
    Output:
        number of protons out of the nozzle per second
    '''
    #using limited log file data from 70 to 150 MeV
    #MUtoNA = (3.023210E-02)*np.power(beamEne,2) - (3.730140)*beamEne + 1.568968E+02
    #average trans factor good to <2%, 15% max at high ene, 2.2% stdev
    transFactor = np.power(10,(0.000000475)*(beamEne**3)-0.000213*(beamEne**2)+0.0422*beamEne-5.672)
    
    MUtoNA = 0.1302558*np.exp(-0.02447910*beamEne)
    estCurrent = DR*MUtoNA
    ionsPerS = transFactor*estCurrent/(1.602E-10)
    
    #using a different approach
    # ionsPerMU = -27.467*(beamEne**2)+25935*(beamEne)+484599 #this is a fit from raystation commisioning data
    # ionsPerS = ionsPerMU*DR/60
    
    # #beamRange = 0.0022*beamEne**1.7707#0.0004*(beamEne**2)+0.0627*beamEne-3.012
    # #proToPho = (-1.074102E-09)*(np.power(beamEne,3)) + (5.253571E-07)*(np.power(beamEne,2)) - (4.303668E-05)*beamEne + 1.645424E-03
    # #proToPho = (-0.01394*(beamEne**3) + 7.27633*(beamEne**2) - 429.36744*beamEne + 15131.75103)/20000000
    # estCurrent = ionsPerS*(1.602E-10)/transFactor
   
    if estCurrent > 600:
        decr = 600/estCurrent
        ionsPerS = decr*ionsPerS

    return ionsPerS#*0.425
    
def getSingles(csvName):
    '''
    c = getSingles(csvName)
    
    This function extacts the 'singles' data from a mixed-data CSV file. It assumes
    normal reconData tree format - that is a constant 12 columns for all singles, doubles, and triples.
    Singles only have the first 4 columns populated with edep, x, y, z, and the other 8 zeroed out.
    It does not use any special delimiter; just the python default.
    
    Input:
        csvName : file name of the data file
    Output:
        list of singles (edep,x,y,z)
    '''
    detM = pd.read_csv(csvName,names=range(12))
    detM = np.nan_to_num(detM)
    #detM = b
    sInd = detM[:,4] == 0
    newList = list(compress(range(len(sInd)),sInd))
    c = np.zeros((len(newList),4))
    c[range(len(newList)),:4]=detM[newList,:4]
    return c


def getSinglesErgInit(csvName):
    '''

    c = getSinglesErgInit(csvName)
    
    This function extacts the 'singles' data from a mixed-data CSV file. It assumes
    normal reconData tree format - that is a constant 12 columns for all singles, doubles, and triples.

    Singles only have the first 4 columns populated with edep, x, y, z, and the other 8 zeroed out.
    It does not use any special delimiter; just the python default.
    

    Input:
        csvName : file name of the data file
    Output:
        list of singles (edep,x,y,z)

    '''
    detM = pd.read_csv(csvName,names=range(13))
    detM = np.nan_to_num(detM)
    print(f"\n[getSinglesErgInit] Loaded singles with shape: {detM.shape}")
    print(f"One row: {detM[0]}")
    sInd = detM[:,4] == 0
    newList = list(compress(range(len(sInd)),sInd))
    c = np.zeros((len(newList),5))
    c[range(len(newList)),:4]=detM[newList,:4]
    c[range(len(newList)),4]=detM[newList,12]
    print(f"[getSinglesErgInit] Filtered singles to shape: {c.shape}, {100*c.shape[0]/detM.shape[0]} %")
    print(f"One row: {c[0]}")
    return c

    
def getDoubles(csvName):
    '''
    c = getDoubles(csvName)
    
    This function extacts the 'doubles' data from a mixed-data CSV file. It assumes
    normal reconData tree format - that is a constant 12 columns for all singles, doubles, and triples.
    Doubles have the first 8 columns populated with edep, x, y, z, and the other 4 zeroed out.
    It does not use any special delimiter; just the python default.
    
    Input:
        csvName : file name of the data file
    Output:
        list of doubles (edep1,x1,y1,z1,edep2,x2,y2,z2)
    '''
    detM = pd.read_csv(csvName,names=range(12))
    detM = np.nan_to_num(detM)
    #detM = b
    sInd = detM[:,4] == 0
    sorDind = detM[:,8] == 0
    dInd = np.logical_and(np.invert(sInd),sorDind)
    newList = list(compress(range(len(dInd)),dInd))
    c = np.zeros((len(newList),8))
    c[range(len(newList)),:8]=detM[newList,:8]

    return c

def getDoublesErgInit(csvName):
    '''
    c = getDoublesErgInit(csvName)

    
    This function extacts the 'doubles' data from a mixed-data CSV file. It assumes
    normal reconData tree format - that is a constant 12 columns for all singles, doubles, and triples.
    Doubles have the first 8 columns populated with edep, x, y, z, and the other 4 zeroed out.

    It does not use any special delimiter; just the python default.
    
    Input:
        csvName : file name of the data file

    Output:
        list of doubles (edep1,x1,y1,z1,edep2,x2,y2,z2)
    '''
    detM = pd.read_csv(csvName,names=range(13))
    detM = np.nan_to_num(detM)
    print(f"\n[getDoublesErgInit] Loaded doubles with shape: {detM.shape}")
    print(f"One row: {detM[0]}")
    sInd = detM[:,4] == 0
    sorDind = detM[:,9] == 0
    dInd = np.logical_and(np.invert(sInd),sorDind)
    newList = list(compress(range(len(dInd)),dInd))
    c = np.zeros((len(newList),9))
    c[range(len(newList)),:8]=detM[newList,:8]
    c[range(len(newList)),8]=detM[newList,12]
    print(f"[getDoublesErgInit] Filtered doubles to shape: {c.shape}, {100*c.shape[0]/detM.shape[0]} %")
    print(f"One row: {c[0]}")
    return c

    
def getTriples(csvName):
    '''
    c = getTriples(csvName)
    
    This function extacts the 'triples' data from a mixed-data CSV file. It assumes
    normal reconData tree format - that is a constant 12 columns for all singles, doubles, and triples.
    Triples only have the columns populated with edep, x, y, z.
    It does not use any special delimiter; just the python default.
    
    Input:
        csvName : file name of the data file
    Output:
        list of triples (edep1,x1,y1,z1,edep2,x2,y2,z2,edep3,x3,y3,z3)
    '''
    detM = pd.read_csv(csvName,names=range(12))
    detM = np.nan_to_num(detM)
    #detM = b
    #sInd = detM[:,4] == 0
    #sorDind = detM[:,9] == 0
    dInd = detM[:,9]!=0
    #dInd = np.logical_and(np.invert(sInd),sorDind)
    newList = list(compress(range(len(dInd)),dInd))
    c = np.zeros((len(newList),12))
    c[range(len(newList)),:12]=detM[newList,:12]
    return c

def getTriplesErgInit(csvName):
    '''
    c = getTriples(csvName)
    
    This function extacts the 'triples' data from a mixed-data CSV file. It assumes

    normal reconData tree format - that is a constant 12 columns for all singles, doubles, and triples.
    Triples only have the columns populated with edep, x, y, z.
    It does not use any special delimiter; just the python default.
    

    Input:
        csvName : file name of the data file
    Output:
        list of triples (edep1,x1,y1,z1,edep2,x2,y2,z2,edep3,x3,y3,z3)

    '''
    detM = pd.read_csv(csvName,names=range(13))
    detM = np.nan_to_num(detM)
    print(f"\n[getTriplesErgInit] Loaded triples with shape: {detM.shape}")
    print(f"One row: {detM[0]}")
    #sInd = detM[:,4] == 0
    #sorDind = detM[:,9] == 0
    dInd = detM[:,9]!=0
    #dInd = np.logical_and(np.invert(sInd),sorDind)
    newList = list(compress(range(len(dInd)),dInd))
    c = np.zeros((len(newList),13))
    c[range(len(newList)),:12]=detM[newList,:12]
    c[range(len(newList)),12]=detM[newList,12]
    print(f"[getTriplesErgInit] Filtered triples to shape: {c.shape}, {100*c.shape[0]/detM.shape[0]} %")
    print(f"One row: {c[0]}")
    return c

    
def protonStopPowerHDPE(beamEne):
    '''
    dEdX = protonStopPowerHDPE(beamEne)
    
    returns stopping power dE/dx [MeV/cm] in HDPE, assuming density = 0.97 g/cm^3
    
    Input:
        beamEne : proton beam energy, in [MeV]
    Output:
        dEdX : stopping power in [MeV/cm]
    '''
    
    if beamEne <= 0.1:
        msp = (8.2498E+08)*(beamEne**5) - (2.3969E+08)*(beamEne**4) + (2.6865E+07)*(beamEne**3) - (1.5449E+06)*(beamEne**2) + (5.1635E+04)*beamEne + 2.0330E+02
    if beamEne < 0.2:
        msp1 = (8.2498E+08)*(beamEne**5) - (2.3969E+08)*(beamEne**4) + (2.6865E+07)*(beamEne**3) - (1.5449E+06)*(beamEne**2) + (5.1635E+04)*beamEne + 2.0330E+02
        msp2 = (2.6826E+02)*np.power(beamEne,-0.75484)
        msp = (msp1+msp2)/2
    else:
        msp = (2.6826E+02)*np.power(beamEne,-0.75484)
    dEdX = msp*0.97
    return dEdX

def totalCrossSection(beamEne,mbList,eneList):
    '''
    interProb = totalCrossSection(beamEne,mbList,eneList)
    computes the total interaction probablity for a proton incident on HDPE based on the supplied
    cross sections and energies.
    input:
        beamEne : starting energy, should match a value in the eneList array
        mbList : list of differential (in energy) cross sections for a given reaction, in [mb]
        eneList : energy list in [MeV] corresponding to the data points in mbList
    output:
        interProb : total interaction probability, using a CSDA approach from beamEne to stop.
    '''
    densityHDPE = 0.97 #in units of g/cm^3
    molarMassHDPE = 28 #in units of g/mol
    nDense = 2*densityHDPE*(6.023E23)/molarMassHDPE #in units of atoms/cm^3
    interProb = 0
    mb = 1E-27 # 1 mb = 1E-27*cm^2
    dE = np.diff(eneList).mean()

    ### EDIT by Peter Jin 7/1/25 for DEBUGGING
    print(eneList)
    ###

    indStop = np.argmax(beamEne==eneList)
    for iii in range(indStop):
        dEdX  = protonStopPowerHDPE(eneList[iii]) #MeV/cm
        dX = dE/dEdX #cm
        interProb += mb*mbList[iii]*dX*nDense * np.exp(-dX * mb * mbList[iii])
    return interProb

def count511More(beamEne,pPerS,t,num511,numProtons):
    '''
    totalCounts, actB8, actC9, actC10, actC11, actN12 = count511More(beamEne,pPerS,t,effec)
    
    calculates how many 511 keV gammas should be detected for a given proton irradiation on HDPE, 
    as well as pathway ending activities. Looks at B8, C9, C10, C11, N12 pathways. Uses the full generation and decay,
    that is parent->daughter formalism and activation building approach.
    Note: it requires several cross section data files. Valid in the energy range of 10 to 250 MeV
    
    input:
        beamEne : energy of the proton beam [MeV]
        pPerS : number of protons per second incident on the phantom [p/s]
        t : total irradiation time [s]
        num511 : total number of 511 interactions recorded for a given MC setup. It is the total length of the data file.
            for example, for 700 singles, 200 doubles and 100 triples, num511 = 1000.
        numProtons : total number of protons simluated in monte carlo
    output:
        totalCounts : total number of 511 keV gammas that would be detected, assuming no detector timing effects
        COMMENTED OUT actXN : activities of the pathways of B8, C9, C10, C11, N12, N13
    '''

    eneList = np.linspace(10,250,241)
    base_dir = Path(__file__).resolve().parent
    txt_dir = (base_dir/"../txt_in")
    mb8 = np.loadtxt(txt_dir/'C12toB8crossSection.txt')
    mb9 = np.loadtxt(txt_dir/'C12toC9crossSection.txt')
    mb10 = np.loadtxt(txt_dir/'C12toC10crossSection.txt')
    mb11 = np.loadtxt(txt_dir/'C12toC11crossSection.txt')
    mb12 = np.loadtxt(txt_dir/'C12toN12crossSection.txt')
    mb13 = np.loadtxt(txt_dir/'C12toN13crossSection.txt')
    
    actLambda8 = np.log(2)/0.770 #770 ms is the approximate half life of B8
    actLambda9 = np.log(2)/0.1266 #126.6 ms is the approximate half life of C9
    actLambda10 = np.log(2)/19.3 #19.3 s is the approximate half life of C10
    actLambda11 = np.log(2)/1221.6 #1221.6 s is the approximate half life of C11
    actLambda12 = np.log(2)/0.011 #11 ms is the approximate half life of N12
    actLambda13 = np.log(2)/597.9 #597.9 s is the approximate half life of N13
    
    prodProb8 = totalCrossSection(beamEne,mb8,eneList)
    prodProb9 = totalCrossSection(beamEne,mb9,eneList)
    prodProb10 = totalCrossSection(beamEne,mb10,eneList)
    prodProb11 = totalCrossSection(beamEne,mb11,eneList)
    prodProb12 = totalCrossSection(beamEne,mb12,eneList)
    prodProb13 = totalCrossSection(beamEne,mb13,eneList)
    totalProdProb = (prodProb8 + prodProb9 + prodProb10 + prodProb11 + prodProb12 + prodProb13)
    
    ###EDIT FOR DEBUGGING by Peter Jin 7/1/25
    print(numProtons)
    print(totalProdProb)
    ###
    efficency = (num511 / numProtons) / totalProdProb
    
    prodRate8 = prodProb8 * pPerS
    prodRate9 = prodProb9 * pPerS
    prodRate10 = prodProb10 * pPerS
    prodRate11 = prodProb11 * pPerS
    prodRate12 = prodProb12 * pPerS
    prodRate13 = prodProb13 * pPerS
    
    act8 = prodRate8 * (1-np.exp(-actLambda8*t))
    act9 = prodRate9 * (1-np.exp(-actLambda9*t))
    act10 = prodRate10 * (1-np.exp(-actLambda10*t))
    act11 = prodRate11 * (1-np.exp(-actLambda11*t))
    act12 = prodRate12 * (1-np.exp(-actLambda12*t))
    act13 = prodRate13 * (1-np.exp(-actLambda13*t))
    
    totCounts8 = prodRate8 * ((np.exp(-actLambda8*t)-1)/actLambda8 + t)
    totCounts9 = prodRate9 * ((np.exp(-actLambda9*t)-1)/actLambda9 + t)
    totCounts10 = prodRate10 * ((np.exp(-actLambda10*t)-1)/actLambda10 + t)
    totCounts11 = prodRate11 * ((np.exp(-actLambda11*t)-1)/actLambda11 + t)
    totCounts12 = prodRate12 * ((np.exp(-actLambda12*t)-1)/actLambda12 + t)
    totCounts13 = prodRate13 * ((np.exp(-actLambda13*t)-1)/actLambda13 + t)
    
    totalCounts = (totCounts8 + totCounts9 + totCounts10 + totCounts11 + totCounts12 + totCounts13) * efficency

    return totalCounts#,act8,act9,act10,act11,act12,act13
 
def detPix(csvData,resXZ=2*0.86,resY=0.1):
    '''
    pixelized = detPix(csvData,resXZ,resY)
    
    a crude pixelation just based on data flooring. Accepts singles, doubles or triples in normal 4, 8, or 12 column format
    
    input:
        csvData : interaction data in normal (edep, x, y, z) format
        resXZ : pixelation resolution in X and Z; assumd symetric
        resY : y data 'resolution'. In reality this is incorrect; use this method to just approximate uncertainy
    output:
        pixelized : pixelated data.
    '''
    pixelized = csvData.copy()
    temp = np.floor(csvData[:,1]/resXZ)
    pixelized[:,1] = (temp*resXZ)


    temp = np.floor((csvData[:,2]/resY))
    pixelized[:,2] = (temp*resY)
    temp = np.floor((csvData[:,3]/resXZ))
    pixelized[:,3] = (temp*resXZ)
    if csvData.shape[1]>4:
        temp = np.floor((csvData[:,5]/resXZ))
        pixelized[:,5] = (temp*resXZ)
        temp = np.floor((csvData[:,6]/resY))
        pixelized[:,6] = (temp*resY)
        temp = np.floor((csvData[:,7]/resXZ))
        pixelized[:,7] = (temp*resXZ)
    if csvData.shape[1]>8:
        temp = np.floor((csvData[:,9]/resXZ))
        pixelized[:,9] = (temp*resXZ)
        temp = np.floor((csvData[:,10]/resY))
        pixelized[:,10] = (temp*resY)
        temp = np.floor((csvData[:,11]/resXZ))
        pixelized[:,11] = (temp*resXZ)

    return pixelized
    
def getModSings(singlesRaw):
    '''
    mdls = getModSings(singlesRaw)
    
    this function takes in the data from a single set of interactions (edep,x,y,z) and returns the module number of each interaction.
    Assumes 16 modules. It can be passed any interaction quartet, e.g. (edep2,x2,y2,z2) from a triple.

    input: 
        singlesRaw : list of (edep,x,y,z) values
    output: 
        mdls : module numbers [0,15]
    '''
    singles = detPix(singlesRaw,0.1,0.1)
    yMean = np.unique(singles[:,2]).mean() #two options
    xMean = np.unique(singles[:,1]).mean() #two options over 4 crystal locs
    zMean = np.unique(singles[:,3]).mean() #four options over 8 crystal locs
    zLeft = (singles[singles[:,3]<zMean,3].min()+singles[singles[:,3]<zMean,3].max())/2
    zRight = (singles[singles[:,3]>zMean,3].max()+singles[singles[:,3]>zMean,3].min())/2
    locs = np.zeros((singles.shape[0],3))
    locs[singles[:,2]<yMean,1] = -1
    locs[singles[:,2]>yMean,1] = 1
    locs[singles[:,1]<xMean,0] = -1
    locs[singles[:,1]>xMean,0] = 1
    locs[singles[:,3]<zLeft,2] = -2
    locs[np.logical_and(singles[:,3]>zLeft,singles[:,3]<zMean),2] = -1
    locs[np.logical_and(singles[:,3]>zMean,singles[:,3]<zRight),2] = 1
    locs[singles[:,3]>zRight,2] = 2
    
    mdls = modListGen(locs)
    
    return np.array((mdls)).astype(int)
       
def final511Count(beamEne,pPerS,t,simProtons,num511):
    '''
    final511 = final511Count(beamEne,pPerS,t,simProtons,num511)
    returns the expected number of 511s that should be detected for a given irradation of HDPE. There are several 'gotchas' to this:
        1) this is only accurate in HDPE
        2) it only generates an expected ratiometric number and cannont tell you about the temporal distribution of the 511s.
        3) it looks at the 6 most prominent positron generators in HDPE
        4) Only valid in the energy range of 10 to 250 MeV
    input:
        beamEne : energy of the proton beam [MeV]
        pPerS : number of protons per second incident on the phantom [p/s]
        t : total irradiation time [s]
        simProtons : total number of protons simluated in monte carlo
        num511 : total number of 511 interactions recorded for a given MC setup. It is the total length of the data file.
            for example, for 700 singles, 200 doubles and 100 triples, num511 = 1000.
    output:
        final511 : total number of 511 keV gammas that would be detected, assuming no detector timing effects
    '''
    
    totalCounts = count511More(beamEne,pPerS,t,num511,simProtons)
    non511Rat = pPerS*t/simProtons
    final511 = totalCounts/non511Rat
    
    return final511
    
def returnPixelLocs(events,numPixX=11,numPixZ=11):
    '''
    xS,zS = returnPixelLocs(events,numPixX,numPixZ)
    
    This function returns the center locations of where pixels should be based on the supplied data. 
    It assumes a 16 module system, with 4 crystals, and each crystal has the specified number of pixels in X and Z.
    
    input:
        events : set of singles-like (edep, x, y, z) data for all 16 modules
        numPixX, numPixZ : number of pixels per crystal in X, Z respectively
    output:
        xS,zS : location of each pixel in each module and crystal.
            size of xS and zS is (16,4,numPix), meaning indexing is (module number, crystal number, pixel number)
    
    '''
    wash = 0.0000001
    sideBuff = 0
    diffThresh = 0.2
    pixOut = events.copy()
    test = np.zeros(events.shape).copy()
    xLocsStore = np.zeros((16,4,numPixX))
    zLocsStore = np.zeros((16,4,numPixZ))
    modCrysList = np.zeros((events.shape[0],2))
    pixLocs = np.zeros((events.shape[0],2))
    modCrysList[:,0] = getModSings(events[:,:4])
    #singles
    pixInd = 0

    indMasterList = np.arange(modCrysList.shape[0])
    for iii in range(16):
        print(f"module {iii}")
        q = events[modCrysList[:,0]==iii,:4].copy()
        print(f"q shape: {q.shape}")
        print(f"q {q}")
        listPart = indMasterList[modCrysList[:,0]==iii]
        c = getCrystal(q)
        print(f"c {c}")
        for jjj in range(4):
            print(f"crystal {jjj}")
            indList = listPart[c==jjj]
            modCrysList[indList,1] = jjj
            temp = q[c==jjj,:].copy()
            print(f"temp shape: {temp.shape}")
            print(f"temp[:,1] = {temp[:,1] if temp.shape[1] > 1 else 'N/A'}")
            print(f"temp[:,3] = {temp[:,3] if temp.shape[1] > 3 else 'N/A'}")

            if temp.size == 0 or temp.shape[0] == 0 or temp.shape[1] < 4:
           # if temp.shape[0] == 0:
                print(f"Skipping module {iii}, crystal {jjj} due to empty data.")
                xLocsStore[iii,jjj,:] = np.linspace(0,0,11)
                zLocsStore[iii,jjj,:] = np.linspace(0,0,11)
                continue  # Skip this crystal
            xMin = temp[:,1].min()
            xMax = temp[:,1].max()
            xDiff = np.abs(xMin-xMax)
            zMin = temp[:,3].min()
            zMax = temp[:,3].max()
            zDiff = np.abs(zMin-zMax)
            xLocsStore[iii,jjj,:] = np.linspace(xMin+sideBuff,xMax-sideBuff,11)
            zLocsStore[iii,jjj,:] = np.linspace(zMin+sideBuff,zMax-sideBuff,11)
            
    xd = np.unique(xLocsStore)
    xx = np.diff(xd)
    xLen = len(xx[xx>diffThresh])+1
    xx = np.append(xx,xx.mean())
    xS = np.zeros((xLen,))
    jjj = 0
    for iii in range(xLen):
        temp = []
        temp.append(xd[jjj])
        
        while (xx[jjj]<diffThresh) and (jjj<(len(xd)-1)) and (jjj<(len(xx)-1)):
            if not(xd[jjj] in temp):
                temp.append(xd[jjj])
            jjj += 1
        xS[iii] = np.mean(np.array(temp))
    #    if (iii == (xLen-1)) and np.isnan(xS[-1]):
    #        xS[-1] = xd[-1]
        jjj += 1
        
        
    zd = np.unique(zLocsStore)
    zz = np.diff(zd)
    zLen = len(zz[zz>diffThresh])+1
    zz = np.append(zz,zz.mean())
    zS = np.zeros((zLen,))
    jjj = 0

        
    for iii in range(zLen):
        temp = []
        temp.append(zd[jjj])
        
        while (zz[jjj]<diffThresh) and (jjj<(len(zd)-1)) and (jjj<(len(zz)-1)):
            if not(zd[jjj] in temp):
                temp.append(zd[jjj])
            jjj += 1
        zS[iii] = np.mean(np.array(temp))
    #    if (iii == (xLen-1)) and np.isnan(xS[-1]):
    #        xS[-1] = xd[-1]
        jjj += 1
    return xS,zS

def getNewTime(modList,pixList,tCurrent,tSamp,tShape=0.1E-6):
    '''
    
    tNew,tSamp = getNewTime(modList,pixList,tCurrent,tSamp,tShape)
    
    this function takes in the modList, current (well, last) clock time, and the sampled deltaT
    and returns the new current time and the new sampled deltaT.
    
    inputs:
        modList, pixList : list containing the module and pixel objects, respectively
        tCurrent : current time
        tSamp : deltaT difference, essentially time till next photon event
        tShape : depricated, to be removed later
    outputs :
        tNew : new current time
        tSamp : new deltaT
    '''
    #tShape = 0.1E-6
    tList = getDataMod(modList,'n')
    #tList2 = getDataPix(pixList,'t')
    #tList2 = tList2[tList2>tCurrent] #this is handeling the shaping time for 
    #tList = np.concatenate((tList,tList2))
    if tCurrent == 0:
        tNew = tSamp.copy()
        tSamp = 0
    else:
        if tList.min()-tCurrent < tSamp: #means there is some processing that needs to happen before next event is popped off the stack
            tNew = tList.min()
            tSamp = tSamp - (tNew-tCurrent)
        else:
            tNew = tCurrent+tSamp.copy()
            tSamp = 0
    return tNew,tSamp
    
def processModStatusPulse(pixelList,modList,tNew,newSings,newDoubs,newTrips,countOverflow,trigThresh=0.05,resTime=1.5E-6,deadScl=1):

    '''
    pixelList,modList,newSings,newDoubs,newTrips,countOverflow = processModStatusPulse(pixelList,modList,tNew,newSings,newDoubs,newTrips,countOverflow,trigThresh=0.05,resTime=1.5E-6)
    

    The big boy! This is one of the the primary control functions for the detector simulation. It handles all module status-update calls, data cleanup, and output updating.
    Longer description TBD if I want to go into detail about the three module statuses, and what is done for each
    
    Pure input:
        tNew : current time that is checked against module status change times to determine which modules are changing
        trigThresh : trigger threshold for a pixel, in [MeV]
        resTime : active time after a trigger, [s]
    Input/output:
        pixelList, modList : list containing the module and pixel objects, respectively
        newSings, newDoubs, newTrips : list of output singles, doubles, and triples that have been processed, respectively. now includes gamma initial energy value
        countOverflow : an array that tracks special interaction types
            countOverflow[0] : number of read out events that had more than 3 interactions (quads). They are not output, just counted
            countOverflow[1] : number of false doubles
            countOverflow[2] : number of false triples
            countOverflow[3] : number of D->Ts
            countOverflow[4] : total time spent dead [s]
            countOverflow[5] : currently not used
            
    '''
    #trigThresh = 0.05
    #b = np.exp(152.39*trigThresh)-1
    #beta = 152.39
    #trigThresh = np.log(b*random.rand()+1)/beta
    tShape = 0.1E-6
    #resTime = 1.5E-6
    modInds = np.arange(len(modList))[getDataMod(modList,'n')==tNew]
    #print('len modList = %i' % len(modList))
    #print('len modInds = %i' % len(modInds))
    for mI in modInds:
        modStatus = modList[mI].status
        if modStatus == 0: #means was in ready, check if new event triggers
            pI = []
            #print('mod status = %i' % modStatus)
            
            for iii in np.arange(len(pixelList)):
                p = pixelList[iii]
                #if (p.module) == mI and (p.time + tShape == tNew):
                if (p.module) == mI and (p.tStop == tNew):
                    if p.eOut >= trigThresh:
                        #print('mod status = %i' % modStatus)
                        modList[mI].status = 1
                        modList[mI].lastTime = p.tTrig + random.rand()*1E-10
                        #modList[mI].lastTime = p.time
                        #modList[mI].trigThresh = trigThresh
                        #print('*** set mod newTime ***')
                        #modList[mI].nextTime = tNew.copy() + resTime
                        modList[mI].nextTime = p.tTrig + resTime
                    else:
                        pI.append(iii)
                if modList[mI].status == 0:
                    modList[mI].nextTime = float('inf')
            for iii in pI[::-1]:
                del(pixelList[iii])
                        
        elif modStatus == 1: #means was in triggered,
            #print('mod status = %i.  triggered mod = %i' % (modStatus, mI))
            #store data with edep >= trigThresh and lastTime<=tpix<=nextTime
            #calc deadtime add
            #clear data, set to dead
            pI = []
            pix = []
            for iii in np.arange(len(pixelList)):
                p = pixelList[iii]
                if (p.module == mI):
                    #trigThresh = modList[mI].trigThresh
                    pI.append(iii)
                    #print('p eOut is %0.5f' % p.eOut)
                    #print('p time is %0.9f' % p.time)
                    #print('p tMax is %0.9f' % p.tMax)
                    #print('mod last time is %0.9f' % modList[mI].lastTime)
                    #print('mod next time is %0.9f' % (modList[mI].nextTime))
                    #if (p.edep >= trigThresh) and (p.time >= modList[mI].lastTime) and (p.time <= modList[mI].nextTime - tShape):
                    if (p.eOut >= trigThresh) and (p.tMax >= modList[mI].lastTime) and (p.tMax <= modList[mI].nextTime):
                        #print('Found good event: len pixLIst = %i' % len(pixelList))
                        pix.append(p)
            numCrys = len(np.unique(getDataPix(pix,'c')))
            
            if len(pix) == 1 : #single
                sing = np.array((pix[0].edep,pix[0].rawX,pix[0].rawY,pix[0].rawZ,pix[0].ergInit,modList[mI].lastTime,pix[0].peakTrack))
                newSings.append(sing)
                #print('recorded single in mod %i' % mI)
            elif len(pix) == 2: #double
                doub = []
                pt = []
                #
                for p in pix:
                    doub.append(p.edep)
                    doub.append(p.rawX)
                    doub.append(p.rawY)
                    doub.append(p.rawZ)
                    pt.append(p.peakTrack)
                #doub = np.array((pix[0].edep,pix[0].rawX,pix[0].rawY,pix[0].rawZ,pix[0],pix[0].edep,pix[1].rawX,pix[1].rawY,pix[1].rawZ,pix[1].ergInit)
                doub.append(p.ergInit)
                doub.append(modList[mI].lastTime)
                if pix[0].time != pix[1].time: #false double, incr countOverflow[1]
                    countOverflow[1] += 1 
                    doub.append(1)
                else:
                    doub.append(0)

                doub.append(pt[0])
                doub.append(pt[1])
                    
                newDoubs.append(doub)
                #print('recorded Double in mod %i' % mI)
                #print(' new DOub number = %i ' % len(newDoubs))
            elif len(pix) == 3:  #triple
                trip = []
                pt = []
                for p in pix:
                    trip.append(p.edep)
                    trip.append(p.rawX)
                    trip.append(p.rawY)
                    trip.append(p.rawZ)
                    pt.append(p.peakTrack)
                    #print('initial energy = %g' % p.ergInit) 
                #trip.append(p.ergInit)                    
                #trip.append(modList[mI].lastTime)
                
                if (pix[0].time == pix[1].time) and (pix[0].time != pix[2].time):
                    countOverflow[3] += 1
                    trip.append(pix[0].ergInit)                    
                    trip.append(modList[mI].lastTime)                                 
                    trip.append(2)
                    #print('D to T -- double first') 
                elif (pix[0].time != pix[1].time) and (pix[1].time == pix[2].time):
                    countOverflow[3] += 1
                    trip.append(pix[1].ergInit)                    
                    trip.append(modList[mI].lastTime)                                     
                    trip.append(3)
                    #print('D to T -- single first') 
                elif (pix[0].time != pix[1].time) and (pix[1].time != pix[2].time) and (pix[0].time != pix[2].time): 
                    countOverflow[2] += 1
                    trip.append(pix[2].ergInit)                    
                    trip.append(modList[mI].lastTime)                                     
                    trip.append(1)
                    #print('false triple') 

                else:
                    trip.append(pix[0].ergInit)                    
                    trip.append(modList[mI].lastTime)    
                    trip.append(0)
                    #print('true triple')
                trip.append(pt[0])    
                trip.append(pt[1])
                trip.append(pt[2])
                newTrips.append(trip)  
                #print('new triple initial energy: %g' % trip[12])
                #print('recorded Triple in mod %i' % mI)
            elif len(pix)!=0:
                #true overflow, incr countOverflow[0]
                #print('Number of pix is %i' % len(pix))
                countOverflow[0] += 1
            deadAdd = deadScl * (np.sum(np.random.uniform(121*(40E-9),121*(50E-9),numCrys)) + (19E-6)*len(pix) - 9E-6 + np.random.uniform(0,5E-6))
            countOverflow[4] += deadAdd
            modList[mI].status = 2
            modList[mI].lastTime = tNew
            modList[mI].nextTime = tNew + deadAdd

            for pp in pI[::-1]:
                del(pixelList[pp])
                
        elif modStatus == 2: #means was in dead
            #print('mod status = %i ' % (modStatus))
            #print('Dead State mod = %i' % (mI))
            #change to ready
            modList[mI].status = 0
            modList[mI].lastTime = tNew
            modList[mI].nextTime = float('inf')
            
    return pixelList,modList,newSings,newDoubs,newTrips,countOverflow
    
def updateNextTimePulse(pixList,modList):
    '''
    modList = updateNextTimePulse(pixList,modList)
    
    This function determines which the next time a module's status changes.
    
    input/output:
        pixelList, modList : list containing the module and pixel objects, respectively
    '''
    #tShape = 0.1E-6
    for p in pixList:
        mI = int(p.module)
        if modList[mI].status == 0:
            #pt = p.time
            nT = np.min((p.tStop,modList[mI].nextTime))
            #nT = np.min((p.time + p.tShape,modList[mI].nextTime))
            modList[mI].nextTime = nT
    return modList

def eCloudRad(ene):
    '''
    r = eCloudRad(ene)
    
    This function determines the initial radius of an electron cloud in CZT. 
    
    input:
        ene : energy deposited, [MeV]
    output: 
        r : electron cloud radius, [mm?]
    '''
    
    #generated from fit of center points of figure 3.3 in feng's diss
    r = (0.00879*np.power(ene,2) + 0.02714*ene)/2
    return r
       
def getXZFromMandC(xL,zL,mNum,cNum):
    '''
    xOut,zOut = getXZFromMandC(xL,zL,mNum,cNum)
    returns the list of x's and z's from a given reference list based on module number and crystal number

    input:
        xL,zL : reference xList and zList for pixelated event locations; assumed unique list
        mNum,cNum : module number (0-16) and crystal number (0-4) 

    output:
        xOut,zOut : unique list of pixel locations in x and z, not grouped
    '''
    #xOut = np.zeros((len(mNum),))
    #zOut = np.zeros((len(mNum),))
    xCenter = xL.mean()
    zCenter = zL.mean()
    xP = xL[xL>xCenter].copy()
    xP03 = xP[xP>xP.mean()]
    xP12 = xP[xP<xP.mean()]
    xN = xL[xL<xCenter].copy()
    xN03 = xN[xN>xN.mean()]
    xN12 = xN[xN<xN.mean()]
    
    zP = zL[zL>zCenter].copy()
    zPP = zP[zP>zP.mean()]
    zPN = zP[zP<zP.mean()]
    zN = zL[zL<zCenter].copy()
    zNP = zN[zN>zN.mean()]
    zNN = zN[zN<zN.mean()]
    
    indXP = np.logical_or(mNum<4,np.logical_and(mNum>=8,mNum<=11)) #x is positive
    if indXP:
        if np.logical_or(cNum==0,cNum==3):
            xOut = xP03
        else:
            xOut = xP12
    else:
        if np.logical_or(cNum==0,cNum==3):
            xOut = xN03
        else:
            xOut = xN12
            
    
    indZPP = mNum in (0,7,11,12)
    indZPN = mNum in (1,6,10,13)
    indZNP = mNum in (2,5,9,14)
    indZNN = mNum in (3,4,8,15)
    
    if indZPP:
        if cNum<=1:
            zOut = zPP[zPP>zPP.mean()]
        else:
            zOut = zPP[zPP<zPP.mean()]
            
    if indZPN:
        if cNum<=1:
            zOut = zPN[zPN>zPN.mean()]
        else:
            zOut = zPN[zPN<zPN.mean()]
            
    if indZNP:
        if cNum<=1:
            zOut = zNP[zNP>zNP.mean()]
        else:
            zOut = zNP[zNP<zNP.mean()]
            
    if indZNN:
        if cNum<=1:
            zOut = zNN[zNN>zNN.mean()]
        else:
            zOut = zNN[zNN<zNN.mean()]

            
    return xOut,zOut
        
def findActivePixels(xTest,zTest,yTest,eDep,xRef,zRef,mInd,cInd,plane1Anode=295,plane2Anode=327.5,V=2000,crysThick1=10,crysThick2=15,rReduction=1):
    '''
    eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(xTest,zTest,yTest,eDep,xRef,zRef,mInd,cInd,plane1Anode,plane2Anode,V)
    
    This function determines which of the four pixels nearest to an interaction receive some collected charge. It then reports the proper amount of "energy depositied"
    based on the volume of the charge received. This does include electron diffusion towards the anode, but does not include charge trapping.
    
    input:
        xTest,zTest,yTest : x, y, and z position of the interaction, respectively
        eDep : energy deposited [MeV]
        xRef,zRef : reference pixel center locations in x and z, respectively
        mInd,cInd : module and crystal index of the interaction, respectively
        plane1Anode,plane2Anode : first plane and second plane anode locations, respectively
        V : applied detector voltage [V]
    output: 
        eDepPix : energy deposited in each pixel
        pixXInd,pixZInd : pixel indices for pixels that receive charge
        xC,zC : approximately locations of charge clouds for each pixel
    '''
    if mInd >= 16:
        cInd = np.mod(mInd-16,4)
        mInd = np.floor((mInd-16)/4)
    xO,zO = getXZFromMandC(xRef,zRef,mInd,cInd)
    xLims = np.zeros((len(xO)+1,))
    zLims = np.zeros((len(zO)+1,))
    xPix = np.mean(np.diff(xO))
    zPix = np.mean(np.diff(zO))
    xLims[:-1] = xO-xPix/2
    xLims[-1] = xO[-1]+xPix/2
    zLims[:-1] = zO-zPix/2
    zLims[-1] = zO[-1]+zPix/2

    #xTest = random.uniform(low=xO.min(),high=xO.max())
    #zTest = random.uniform(low=zO.min(),high=zO.max())
    if mInd <=7:
        D = crysThick1
        d = yTest - plane1Anode
    else:
        D = crysThick2
        d = plane2Anode - yTest

    r = eCloudRad(eDep)+1.175*eCloudSigma(d,D,V)
    r = r/rReduction

    xP = xTest + r
    xM = xTest - r
    zP = zTest + r
    zM = zTest - r

    # sanitize for <class 'numpy.ndarray'> (1,) when expecting scalar
    xP = xP.item() if isinstance(xP, np.ndarray) else xP 
    xM = xM.item() if isinstance(xM, np.ndarray) else xM
    zP = zP.item() if isinstance(zP, np.ndarray) else zP
    zM = zM.item() if isinstance(zM, np.ndarray) else zM

    xMI = bisect.bisect_right(xLims,xM)-1
    xPI = bisect.bisect_right(xLims,xP)

    zMI = bisect.bisect_right(zLims,zM)-1
    zPI = bisect.bisect_right(zLims,zP)

    pixXInd = np.arange(len(xO))[np.logical_and(xO<=xLims[xPI],xO>=xLims[xMI])]
    pixZInd = np.arange(len(zO))[np.logical_and(zO<=zLims[zPI],zO>=zLims[zMI])]
   
    #eDepPix is listed in increasing order
    if (len(pixXInd) == 1) and (len(pixZInd) == 1): #this means only a single pixel was activated
        eDepPix = eDep
        xC = xTest
        zC = zTest
        
    elif (len(pixXInd) == 1) and (len(pixZInd) != 1): #activated 2 in z direction
        eDepPix = np.zeros((1,len(pixZInd)))
        sVol = (4/3)*np.pi*np.power(r,3)
        zInd = zPI-1
        zDiv = zLims[zInd]
        xC = np.array((xTest,xTest))
        if zTest >= zDiv: #means subtract from below, interaction above dividing line
            vCap = np.pi*(np.power((zDiv-zM),2)/3)*(3*r-(zDiv-zM))
            eDepPix[0,0] = eDep*vCap/sVol
            eDepPix[0,1] = eDep - eDepPix[0,0]
            zC = np.array(((zM*eDepPix[0,0]/eDep)+(zDiv*eDepPix[0,1]/eDep),(zDiv*eDepPix[0,0]/eDep)+(zP*eDepPix[0,1]/eDep)))
        else:
            vCap = np.pi*(np.power((zP-zDiv),2)/3)*(3*r-(zP-zDiv))
            eDepPix[0,1] = eDep*vCap/sVol
            eDepPix[0,0] = eDep - eDepPix[0,1]
            zC = np.array(((zM*eDepPix[0,0]/eDep)+(zDiv*eDepPix[0,1]/eDep),(zDiv*eDepPix[0,0]/eDep)+(zP*eDepPix[0,1]/eDep)))
        #pixXInd = (np.zeros((pixZInd.shape))+pixXInd[0]).astype(int)
        #zC = np.array(((zM*eDepPix[0,0]/eDep)+(zDiv*eDepPix[0,1]/eDep),(zDiv*eDepPix[0,0]/eDep)+(zP*eDepPix[0,1]/eDep)))
    elif (len(pixXInd) != 1) and (len(pixZInd) == 1): #activated 2 in x direction
        eDepPix = np.zeros((len(pixXInd),1))
        sVol = (4/3)*np.pi*np.power(r,3)
        xInd = xPI-1
        xDiv = xLims[xInd]
        zC = np.array((zTest,zTest))
        if xTest >= xDiv: #means subtract from below, interaction above dividing line
            vCap = np.pi*(np.power((xDiv-xM),2)/3)*(3*r-(xDiv-xM))
            eDepPix[0,0] = eDep*vCap/sVol
            eDepPix[1,0] = eDep - eDepPix[0,0]
            print("xTest:", type(xTest), np.shape(xTest))
            print("xDiv:", type(xDiv), np.shape(xDiv))
            print("xP:", type(xP), np.shape(xP))
            print("eDep:", type(eDep), np.shape(eDep))
            print("eDepPix:", type(eDepPix), np.shape(eDepPix))
            xC = np.array(((xM*eDepPix[0,0]/eDep)+(xDiv*eDepPix[1,0]/eDep),(xDiv*eDepPix[0,0]/eDep)+(xTest*eDepPix[1,0]/eDep)))
        else:
            vCap = np.pi*(np.power((xP-xDiv),2)/3)*(3*r-(xP-xDiv))
            eDepPix[1,0] = eDep*vCap/sVol
            eDepPix[0,0] = eDep - eDepPix[1,0]

            xC = np.array(((xTest*eDepPix[0,0]/eDep)+(xDiv*eDepPix[1,0]/eDep),(xDiv*eDepPix[0,0]/eDep)+(xP*eDepPix[1,0]/eDep)))
        #pixZInd = (np.zeros((pixXInd.shape))+pixZInd[0]).astype(int)
        #xC = np.array(((xM*eDepPix[0,0]/eDep)+(xDiv*eDepPix[1,0]/eDep),(xDiv*eDepPix[0,0]/eDep)+(xP*eDepPix[1,0]/eDep)))
    else: #activated 4 pixels
        xInd = xPI-1
        zInd = zPI-1
        eDepPix = np.zeros((len(pixXInd),len(pixZInd)))
        xDiv = xLims[xInd]
        zDiv = zLims[zInd]
        sVol = (4/3)*np.pi*np.power(r,3)
        xC = xTest
        zC = zTest
        if np.sqrt(np.power(zLims[zInd]-zTest,2)+np.power(xLims[xInd]-xTest,2)) >= r: #means intersection outside of circle, only 3 pixels activated
            if (xTest >= xDiv) and (zTest >= zDiv): # centered xGreater,zGreater, eDepPix[0,0] = 0
                xC = np.zeros((2,2))
                zC = np.zeros((2,2))
                eDepPix[0,0] = 0 #lower left
                vCapX = np.pi*(np.power((xDiv-xM),2)/3)*(3*r-(xDiv-xM))
                vCapZ = np.pi*(np.power((zDiv-zM),2)/3)*(3*r-(zDiv-zM))
                
                eDepPix[0,1] = eDep*vCapX/sVol # upper left
                xC[0,1] = (xM*eDepPix[0,1]/eDep)+(xDiv*(1-eDepPix[0,1]/eDep))
                zC[0,1] = zTest
                
                eDepPix[1,0] = eDep*vCapZ/sVol # lower right
                xC[1,0] = xTest
                zC[1,0] = (zM*eDepPix[1,0]/eDep)+(zDiv*(1-eDepPix[1,0]/eDep))
                
                eDepPix[1,1] = eDep - np.sum(eDepPix)
                xC[1,1] = (xTest*eDepPix[1,1]/eDep)+(xDiv*(1-eDepPix[1,1]/eDep))
                zC[1,1] = (zTest*eDepPix[1,1]/eDep)+(zDiv*(1-eDepPix[1,1]/eDep))
                
            elif (xTest < xDiv) and (zTest < zDiv): # centered xLess,zLess, eDepPix[1,1] = 0
                xC = np.zeros((2,2))
                zC = np.zeros((2,2))
                eDepPix[1,1] = 0 #upper right
                vCapX = np.pi*(np.power((xP-xDiv),2)/3)*(3*r-(xP-xDiv))
                vCapZ = np.pi*(np.power((zP-zDiv),2)/3)*(3*r-(zP-zDiv))
                
                eDepPix[0,1] = eDep*vCapZ/sVol # upper left
                xC[0,1] = xTest
                zC[0,1] = (zP*eDepPix[0,1]/eDep)+(zDiv*(1-eDepPix[0,1]/eDep))
                
                eDepPix[1,0] = eDep*vCapX/sVol # lower right
                xC[1,0] = (xP*eDepPix[1,0]/eDep)+(xDiv*(1-eDepPix[1,0]/eDep))
                zC[1,0] = zTest
                
                eDepPix[0,0] = eDep - np.sum(eDepPix) # lower left
                xC[0,0] = (xTest*eDepPix[0,0]/eDep)+(xDiv*(1-eDepPix[0,0]/eDep))
                zC[0,0] = (zTest*eDepPix[0,0]/eDep)+(zDiv*(1-eDepPix[0,0]/eDep))
                
                
                
            elif (xTest >= xDiv) and (zTest < zDiv): # centered xGreater,zLess, eDepPix[0,1] = 0
                xC = np.zeros((2,2))
                zC = np.zeros((2,2))
                eDepPix[0,1] = 0 #upper left
                vCapX = np.pi*(np.power((xDiv-xM),2)/3)*(3*r-(xDiv-xM))
                vCapZ = np.pi*(np.power((zP-zDiv),2)/3)*(3*r-(zP-zDiv))
                
                eDepPix[0,0] = eDep*vCapX/sVol # lower left
                xC[0,0] = (xM*eDepPix[0,0]/eDep)+(xDiv*(1-eDepPix[0,0]/eDep))
                zC[0,0] = zTest
                
                eDepPix[1,1] = eDep*vCapZ/sVol # upper right
                xC[1,1] = xTest
                zC[1,1] = (zP*eDepPix[1,1]/eDep)+(zDiv*(1-eDepPix[1,1]/eDep))
                
                eDepPix[1,0] = eDep - np.sum(eDepPix) # lower right
                xC[1,0] = (xTest*eDepPix[1,0]/eDep)+(xDiv*(1-eDepPix[1,0]/eDep))
                zC[1,0] = (zTest*eDepPix[1,0]/eDep)+(zDiv*(1-eDepPix[1,0]/eDep))
                
            elif (xTest < xDiv) and (zTest >= zDiv): # centered xLess, zGreater, eDepPix[1,0] = 0
                xC = np.zeros((2,2))
                zC = np.zeros((2,2))
                eDepPix[1,0] = 0 #lower right
                vCapX =  np.pi*(np.power((xP-xDiv),2)/3)*(3*r-(xP-xDiv))
                vCapZ =  np.pi*(np.power((zDiv-zM),2)/3)*(3*r-(zDiv-zM))
                
                eDepPix[1,1] = eDep*vCapX/sVol # upper left
                xC[1,1] = (xP*eDepPix[1,1]/eDep)+(xDiv*(1-eDepPix[1,0]/eDep))
                zC[1,1] = zTest
                
                eDepPix[0,0] = eDep*vCapZ/sVol # lower right
                xC[0,0] = xTest
                zC[0,0] = (zM*eDepPix[0,0]/eDep)+(zDiv*(1-eDepPix[0,0]/eDep))
                
                eDepPix[0,1] = eDep - np.sum(eDepPix)
                xC[0,1] = (xTest*eDepPix[0,1]/eDep)+(xDiv*(1-eDepPix[0,1]/eDep))
                zC[0,1] = (zTest*eDepPix[0,1]/eDep)+(zDiv*(1-eDepPix[0,1]/eDep))
                
        else: #means all 4 pixels activated
            if (xTest >= xDiv) and (zTest >= zDiv): # centered xGreater,zGreater
                vCapX = np.pi*(np.power((xDiv-xM),2)/3)*(3*r-(xDiv-xM))
                vCapZ = np.pi*(np.power((zDiv-zM),2)/3)*(3*r-(zDiv-zM))
                
                eDepPix[0,0] = eDep*vCapX*vCapZ/np.power(sVol,2)
                eDepPix[1,0] = eDep*(1-vCapX/sVol)*vCapZ/sVol
                eDepPix[0,1] = eDep*(1-vCapZ/sVol)*vCapX/sVol
                eDepPix[1,1] = eDep - np.sum(eDepPix)

                
            elif (xTest < xDiv) and (zTest < zDiv): # centered xLess,zLess
                vCapX = np.pi*(np.power((xP-xDiv),2)/3)*(3*r-(xP-xDiv))
                vCapZ = np.pi*(np.power((zP-zDiv),2)/3)*(3*r-(zP-zDiv))
                
                eDepPix[1,1] = eDep*vCapX*vCapZ/np.power(sVol,2)
                eDepPix[1,0] = eDep*(1-vCapZ/sVol)*vCapX/sVol
                eDepPix[0,1] = eDep*(1-vCapX/sVol)*vCapZ/sVol
                eDepPix[0,0] = eDep - np.sum(eDepPix)

            elif (xTest >= xDiv) and (zTest < zDiv): # centered xGreater,zLess
                vCapX = np.pi*(np.power((xDiv-xM),2)/3)*(3*r-(xDiv-xM))
                vCapZ = np.pi*(np.power((zP-zDiv),2)/3)*(3*r-(zP-zDiv))
                
                eDepPix[0,1] = eDep*vCapX*vCapZ/np.power(sVol,2)
                eDepPix[1,1] = eDep*(1-vCapX/sVol)*vCapZ/sVol
                eDepPix[0,0] = eDep*(1-vCapZ/sVol)*vCapX/sVol
                eDepPix[1,0] = eDep - np.sum(eDepPix)

            elif (xTest < xDiv) and (zTest >= zDiv): # centered xLess, zGreater
                vCapX =  np.pi*(np.power((xP-xDiv),2)/3)*(3*r-(xP-xDiv))
                vCapZ =  np.pi*(np.power((zDiv-zM),2)/3)*(3*r-(zDiv-zM))
                eDepPix[1,0] = eDep*vCapX*vCapZ/np.power(sVol,2)
                eDepPix[0,0] = eDep*(1-vCapX/sVol)*vCapZ/sVol
                eDepPix[0,1] = eDep*(1-vCapZ/sVol)*vCapX/sVol
                eDepPix[1,1] = eDep - np.sum(eDepPix)
            
            xC = np.zeros((2,2))
            zC = np.zeros((2,2))

            xC[0,0] = (xM*eDepPix[0,0]/eDep)+(xDiv*(1-eDepPix[0,0]/eDep))
            zC[0,0] = (zM*eDepPix[0,0]/eDep)+(zDiv*(1-eDepPix[0,0]/eDep))

            xC[1,0] = (xP*eDepPix[1,0]/eDep)+(xDiv*(1-eDepPix[1,0]/eDep))
            zC[1,0] = (zM*eDepPix[1,0]/eDep)+(zDiv*(1-eDepPix[1,0]/eDep))

            xC[0,1] = (xM*eDepPix[0,1]/eDep)+(xDiv*(1-eDepPix[0,1]/eDep))
            zC[0,1] = (zP*eDepPix[0,1]/eDep)+(zDiv*(1-eDepPix[0,1]/eDep))

            xC[1,1] = (xP*eDepPix[1,1]/eDep)+(xDiv*(1-eDepPix[1,1]/eDep))
            zC[1,1] = (zP*eDepPix[1,1]/eDep)+(zDiv*(1-eDepPix[1,1]/eDep))
                
    return eDepPix,pixXInd,pixZInd,xC,zC
       
def calcWeightingPotential(y,p1A,p2A,crysThick1=10,crysThick2=15):
    '''
    edepScale = calcWeightingPotential(y,p1A,p2A)
    returns the energy depositied scaling factor based on the weighting potential method for signal induction.
    Weighting potential calculation is done using a fit calculated from the following reference:
    https://cztlab.engin.umich.edu/wp-content/uploads/sites/187/2015/03/zhangf.pdf
    
    The calculation assumes two opposed CZT planes (one 10 mm, one 15 mm) with anodes on the exterior surfaces.
    
    plane1 (10 mm)            plane 2 (15 mm)
    |                    +y2 <------ y2=p2A |
    |                                       |
    |                                       |
    |y1=p1A ------> +y1                     |
    
    input:
        y : raw interaction depth
        p1A, p2A : y location of the two anodes

    output:
        edepScale : multiplicitive scaling factor for deposited energy
    '''
    edepScale = 1
    if y <= p1A+crysThick1+0.1:
        edepScale = 1 - (0.9808 * np.exp(-23.11*(y-p1A)/crysThick1))
    else:
        edepScale = 1 - (0.9808 * np.exp(-23.11*(p2A-y)/crysThick1))
    return edepScale
    
def getPulseShape(t,d,eDep,enes,depth,muList,cList,aList,bList,dList):
    '''
    l = getPulseShape(t,d,eDep,enes,depth,muList,cList,aList,bList,dList)
    
    This function returns an expected temporal pulse shape based on a some precomupted fitting parameters. It approximately simulates the 
    detector preamp output.
    
    input:
        t : time list to compute pulse shape over
        d : interactoin depth 
        eDep : energy deposited [MeV]
        enes : energy values of the fitting parameters
        depth : depth values of the fitting paramters
        muList,cList,aList,bList,dList : precomputed fitting parameters
    output: 
        l : pulse shape value (normlaized to max of 1) at each input from t
    '''
    t1 = t.copy()
    t1 = t1-t1.min()
    #enes = np.linspace(0.001,5,501)
    #depth = np.linspace(0.01,1,501)

    eInd = bisect.bisect_left(enes,eDep)
    dInd = bisect.bisect_left(depth,d)
    
    if eInd != 0 and eInd != 501 and dInd != 0 and dInd != 501: 
        mu = (muList[dInd,eInd]+muList[dInd-1,eInd-1])/2
        c = (cList[dInd,eInd]+cList[dInd-1,eInd-1])/2
        aa = (aList[dInd,eInd]+aList[dInd-1,eInd-1])/2
        b = (bList[dInd,eInd]+bList[dInd-1,eInd-1])/2
        dd = (dList[dInd,eInd]+dList[dInd-1,eInd-1])/2
        l = landau(t1,mu,c,aa,b,dd)
        l = l/l.max()
    else:
        l = np.zeros((t1.shape))
    #l = l/l.max()
        
    return l
        
def makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,yTest,eInit,mInd,cInd,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,p1Anode=295,p2Anode=327.5,crysThick1=10,crysThick2=15):
    '''
    pList = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,yTest,eInit,mInd,cInd,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,p1Anode=295,p2Anode=327.5)
    
    This function creates a list of pixel objects based on the output of findActivePixels. 
    
    input:
        eDepPix : list of energy deposited in a pixel. 
        pixXInd,pixZInd : list of x, z pixel indicies for each interaction
        xC,zC : list of x, z values for each interaction
        xRef,zRef : reference pixel location lists
        yTest : y value of interaction
        eInit: gamma initial energy
        mInd,cInd : module and crystal indicies for the interaction
        t : time array to generate pulse shape over
        tCurrent : current time 
        gEnes,depths,muList,cList,aList,bList,dList : precomputed values for fitting temporal pulse shape
        p1Anode,p2Anode : anode locations
    output: 
        pList : list of pixel objects containing the input data. 
    '''
    
    pList = []
    mm = mInd
    cc = cInd
    rY = yTest
    #ergInit = eInit
    escl = calcWeightingPotential(rY,p1Anode,p2Anode,crysThick1,crysThick2)
    
    if mm>=16:
        if mm <= 47:
            dth = (rY - p1Anode)/crysThick1
        else:
            dth = (p2Anode - rY)/crysThick2
    else:
        
        if mm <= 7:
            dth = (rY - p1Anode)/crysThick1
        else:
            dth = (p2Anode - rY)/crysThick2
    if len(pixXInd) == 1 and len(pixZInd) == 1: #single pixel
        xx = pixXInd
        zz = pixZInd
        ed = eDepPix
        rX = xC
        rZ = zC
        
        l = getPulseShape(t,dth,ed,gEnes,depths,muList,cList,aList,bList,dList)
        pixNew = pixel(module=mm,crystal=cc,xInd=xx,zInd=zz,edep=ed,rawX=rX,rawY=rY,rawZ=rZ,ergInit=eInit,time=tCurrent,eScl=escl,pulseShape=l)
        pList.append(pixNew)
    elif ((len(pixXInd) == 1) and (len(pixZInd) != 1)) or ((len(pixXInd) != 1) and (len(pixZInd) == 1)): #means 2 pixel
        eDepPix = np.ravel(eDepPix)
        
        for iii in range(2):
            if len(pixXInd) == 1:
                xx = pixXInd
            else:
                xx = pixXInd[iii]
            if len(pixZInd) == 1:
                zz = pixZInd
            else:
                zz = pixZInd[iii]
            ed = eDepPix[iii]
            rX = xC[iii]
            rZ = zC[iii]

            l = getPulseShape(t,dth,ed,gEnes,depths,muList,cList,aList,bList,dList)
            pixNew = pixel(module=mm,crystal=cc,xInd=xx,zInd=zz,edep=ed,rawX=rX,rawY=rY,rawZ=rZ,ergInit=eInit,time=tCurrent,eScl=escl,pulseShape=l)
            pList.append(pixNew)
    else: #means 3 or 4
        for iii in range(2):
            for jjj in range(2):
                ed = eDepPix[iii,jjj]
                if ed > 0:
                    xx = pixXInd[iii]
                    zz = pixZInd[jjj]
                    rX = xC[iii,jjj]
                    rZ = zC[iii,jjj]

                    l = getPulseShape(t,dth,ed,gEnes,depths,muList,cList,aList,bList,dList)
                    pixNew = pixel(module=mm,crystal=cc,xInd=xx,zInd=zz,edep=ed,rawX=rX,rawY=rY,rawZ=rZ,ergInit=eInit,time=tCurrent,eScl=escl,pulseShape=l)
                    pList.append(pixNew)        
    return pList

def combinePix(pix1,pix2,tStep=10E-9):
    '''
    pixComb = combinePix(pix1,pix2,tStep)
    
    This function combines two pixel objects with the same pixel reference. It primarily changes the trigger flag of a pixel, and correctly adds the temporal pulse shapes
    
    input:
        pix1, pix2 : pixel objects to be added
        tStep : time stamp incriment
    output: 
        pixComb : combined pixel object
    '''
    
#    pix1 = pixNew
#    pix2 = pixNew1
    tMin = np.min((pix1.tPulse.min(),pix2.tPulse.min()))
    tMax = np.max((pix1.tPulse.max(),pix2.tPulse.max()))
    if pix1.tPulse.min() == pix2.tPulse.min():
        tNew = pix1.tPulse
    else:
        tNew = np.linspace(tMin,tMax,int((tMax-tMin)/tStep + 1))
    newPulse = np.zeros((len(tNew),))
    newPulse[:len(pix1.pulseShape)] += pix1.pulseShape
    newPulse[-len(pix2.pulseShape):] += pix2.pulseShape
    #newPulse -= 1
    eMax = newPulse.max()#-1
    #newPulse = newPulse/newPulse.max()
    #newPulse = newPulse/newPulse.max()
    
    
    xNew = (pix1.rawX + pix2.rawX)/2
    yNew = (pix1.rawY + pix2.rawY)/2
    zNew = (pix1.rawZ + pix2.rawZ)/2
    if np.sum(newPulse)>0:
        pixComb = pixel(module=pix1.module,crystal=pix1.crystal,xInd=pix1.xInd,zInd=pix1.zInd,\
                    edep=eMax,rawX=xNew,rawY=yNew,rawZ=zNew,time=tMin,eScl=1,\
                    pulseShape=newPulse/newPulse.max(),pulseLength=(tMax-tMin))
    else:
        pixComb = pixel(module=pix1.module,crystal=pix1.crystal,xInd=pix1.xInd,zInd=pix1.zInd,\
            edep=eMax,rawX=xNew,rawY=yNew,rawZ=zNew,time=tMin,eScl=1,\
            pulseShape=np.zeros((len(tNew))),pulseLength=(tMax-tMin))
    #pixComb.pulseShape -= 1
    return pixComb
        
def checkNewEventPulse(pixelList,modList,newEvent,coincCount,deadCount):
    '''  
    pixelList,coincCount,deadCount = checkNewEventPulse(pixelList,modList,newEvent,coincCount,deadCount)
    
    this function checks to see if the new event is acceptable to be put on the pixel stack
    
    input: 
        modList : list of module objects
        newEvent : pixel object of new event
    input/output:
        pixelList : list of pixel objects
        coincCount : count of how many interactions occur in the same pixel during a single readout cycle
        deadCount : number of interactions missed due to a module in the dead state

    '''
    modStatus = modList[newEvent.module].status
    if modStatus <= 1: #use shaping time because I haven't triggered yet
        #print('Mod status = 0')
        diffCount = 0   
        sameCount = 0
        for iii in range(len(pixelList)):
            #if (pixelList[iii].module == newEvent.module) and (pixelList[iii].crystal == newEvent.crystal) and (abs(pixelList[iii].xInd - newEvent.xInd) <= 1) and (abs(pixelList[iii].zInd - newEvent.zInd) <= 1) and (sameCount == 0):
            #if (pixelList[iii].module == newEvent.module) and (pixelList[iii].crystal == newEvent.crystal) and (((abs(pixelList[iii].xInd - newEvent.xInd) <= 1) and (pixelList[iii].zInd == newEvent.zInd)) or ((abs(pixelList[iii].zInd - newEvent.zInd) <= 1) and (pixelList[iii].xInd == newEvent.xInd))) and (sameCount == 0):
            if (pixelList[iii].module == newEvent.module) and (pixelList[iii].crystal == newEvent.crystal) and (pixelList[iii].xInd == newEvent.xInd) and (pixelList[iii].zInd == newEvent.zInd):# and (sameCount == 0):
                if abs(pixelList[iii].time - newEvent.time)>0:
                    coincCount += 1
                pComb = combinePix(pixelList[iii],newEvent)
                
                pixelList[iii] = pComb
                sameCount = 1
            else:
                diffCount += 1
        if diffCount == len(pixelList):
            pixelList.append(newEvent)
    else:
        deadCount += 1
    return pixelList,coincCount,deadCount#, modList

def eneUncCZT(ene):
    '''
    eneUnc = eneUncCZT(ene)
    returns the energy uncertainy for an edep event in CZT based on experimental relationship

    input:
        ene in MeV

    output:
        sigMeV in MeV
    '''
    eneKeV = ene*1000
    sigE = 4 + 0.129*np.sqrt(eneKeV)
    sigEMeV = sigE/1000
    return sigEMeV
    
def modListGen(loc):
    '''
    modList = modListGen(loc):
    
    looks at the 3 valued array from the input and translates it to a module number

    input:
        loc : a three valued array specificying x, y and z locations of a module
    outout:
        modList : module number from 0 to 15
    '''
    modList = np.zeros((loc.shape[0],))
    modList[np.sum((loc[:] == ((1,-1,2))),axis=1)==3] = 0
    modList[np.sum((loc[:] == ((1,-1,1))),axis=1)==3] = 1
    modList[np.sum((loc[:] == ((1,-1,-1))),axis=1)==3] = 2
    modList[np.sum((loc[:] == ((1,-1,-2))),axis=1)==3] = 3
    modList[np.sum((loc[:] == ((-1,-1,-2))),axis=1)==3] = 4
    modList[np.sum((loc[:] == ((-1,-1,-1))),axis=1)==3] = 5
    modList[np.sum((loc[:] == ((-1,-1,1))),axis=1)==3] = 6
    modList[np.sum((loc[:] == ((-1,-1,2))),axis=1)==3] = 7
    modList[np.sum((loc[:] == ((1,1,-2))),axis=1)==3] = 8
    modList[np.sum((loc[:] == ((1,1,-1))),axis=1)==3] = 9
    modList[np.sum((loc[:] == ((1,1,1))),axis=1)==3] = 10
    modList[np.sum((loc[:] == ((1,1,2))),axis=1)==3] = 11
    modList[np.sum((loc[:] == ((-1,1,2))),axis=1)==3] = 12
    modList[np.sum((loc[:] == ((-1,1,1))),axis=1)==3] = 13
    modList[np.sum((loc[:] == ((-1,1,-1))),axis=1)==3] = 14
    modList[np.sum((loc[:] == ((-1,1,-2))),axis=1)==3] = 15

    return modList
    
def getCrystal(eventsRaw):
    '''
    quad1[,quad2, quad3] = getCrystal(eventsRaw)
    
    returns which crystals an interaction took place in. Accepts singles, doubles, or triples. 
    NOTE: only input event data from a single module.
    
    input:
        eventsRaw : data to be analyzed in normal (edep,x,y,z) format. can be S, D or T. 
    outout:
        quadN : crystal of each interaction. Returns 1, 2, or 3 values.
    '''
    
    events = detPix(eventsRaw,0.1,0.1)
    if events.shape[1]==4:
        xMean = np.unique(events[:,1]).mean()
        zMean = np.unique(events[:,3]).mean()
        
        L1 = events[:,1]<xMean
        B1 = events[:,3]<zMean
        
        crysTL1 = np.logical_and(np.invert(B1),L1)
        crysBL1 = np.logical_and(B1,L1)
        crysTR1 = np.logical_and(np.invert(B1),np.invert(L1))
        crysBR1 = np.logical_and(B1,np.invert(L1))

        quad1 = np.zeros((len(L1),))
        
        quad1[crysTR1] = 0
        quad1[crysTL1] = 1
        quad1[crysBL1] = 2
        quad1[crysBR1] = 3
        
        return quad1.astype(int)
        
    elif events.shape[1]==8:
        xMean = np.unique(events[:,(1,5)]).mean()
        zMean = np.unique(events[:,(3,7)]).mean()
        
        L1 = events[:,1]<xMean
        L2 = events[:,5]<xMean
        B1 = events[:,3]<zMean
        B2 = events[:,7]<zMean
        
        crysTL1 = np.logical_and(np.invert(B1),L1)
        crysBL1 = np.logical_and(B1,L1)
        crysTR1 = np.logical_and(np.invert(B1),np.invert(L1))
        crysBR1 = np.logical_and(B1,np.invert(L1))
        
        crysTL2 = np.logical_and(np.invert(B2),L2)
        crysBL2 = np.logical_and(B2,L2)
        crysTR2 = np.logical_and(np.invert(B2),np.invert(L2))
        crysBR2 = np.logical_and(B2,np.invert(L2))    
        
        quad1 = np.zeros((len(L1),))
        quad2 = np.zeros((len(L1),))
        
        quad1[crysTR1] = 0
        quad1[crysTL1] = 1
        quad1[crysBL1] = 2
        quad1[crysBR1] = 3
        
        quad2[crysTR2] = 0
        quad2[crysTL2] = 1
        quad2[crysBL2] = 2
        quad2[crysBR2] = 3
        
        return quad1.astype(int), quad2.astype(int)

    elif events.shape[1]==12:
        xMean = np.unique(events[:,(1,5,9)]).mean()
        zMean = np.unique(events[:,(3,7,11)]).mean()

        L1 = events[:,1]<xMean
        L2 = events[:,5]<xMean
        L3 = events[:,9]<xMean
        B1 = events[:,3]<zMean
        B2 = events[:,7]<zMean
        B3 = events[:,11]<zMean

   
        crysTL1 = np.logical_and(np.invert(B1),L1)
        crysBL1 = np.logical_and(B1,L1)
        crysTR1 = np.logical_and(np.invert(B1),np.invert(L1))
        crysBR1 = np.logical_and(B1,np.invert(L1))
        
        crysTL2 = np.logical_and(np.invert(B2),L2)
        crysBL2 = np.logical_and(B2,L2)
        crysTR2 = np.logical_and(np.invert(B2),np.invert(L2))
        crysBR2 = np.logical_and(B2,np.invert(L2))    
        
        crysTL3 = np.logical_and(np.invert(B3),L3)
        crysBL3 = np.logical_and(B3,L3)
        crysTR3 = np.logical_and(np.invert(B3),np.invert(L3))
        crysBR3 = np.logical_and(B3,np.invert(L3))    
        
        quad1 = np.zeros((len(L1),))
        quad2 = np.zeros((len(L1),))
        quad3 = np.zeros((len(L1),))
        
        quad1[crysTR1] = 0
        quad1[crysTL1] = 1
        quad1[crysBL1] = 2
        quad1[crysBR1] = 3
        
        quad2[crysTR2] = 0
        quad2[crysTL2] = 1
        quad2[crysBL2] = 2
        quad2[crysBR2] = 3

        quad3[crysTR3] = 0
        quad3[crysTL3] = 1
        quad3[crysBL3] = 2            
        quad3[crysBR3] = 3
    
        return quad1.astype(int),quad2.astype(int), quad3.astype(int)
    else:
        return 0
        
def pixelateCommon(events,xRef,zRef,wash):
    '''
    pixOut,modCrysList,pixLocs = pixelateCommon(events,xRef,zRef,wash)
    
    pixelates the supplied event data to the reference data. 
    
    input:
        events : event data to be pixelated
        xRef, zRef : reference pixel centers
        wash : value needed for bit-rounding errors
    outout:
        pixOut : pixelated data
        modCrysList : module numbers of interactions
        pixLocs : pixel index (x,z) of the interaction
    '''
    print(f"\nLength of xRef: {len(xRef)}")
    print(f"Length of zRef: {len(zRef)}")

    numPixX = int(len(xRef)/4)
    numPixZ = int(len(zRef)/8)
    print(f"numPixX (pixels per crystal in X): {numPixX}")
    print(f"numPixZ (pixels per crystal in Z): {numPixZ}\n")
    wash += 0.000001
    sideBuff = 0
    pixOut = events.copy()
    test = np.zeros(events.shape).copy()
    #xL,zL = xL,zL = returnPixelLocs(events[:,:4])
    
    xLocsStore = np.zeros((16,4,numPixX))
    zLocsStore = np.zeros((16,4,numPixZ))
    if events.shape[1]<8:
        modCrysList = np.zeros((events.shape[0],2))
        pixLocs = np.zeros((events.shape[0],2))
        modCrysList[:,0] = getModSings(events[:,:4])
    elif events.shape[1]<12:
        modCrysList = np.zeros((events.shape[0],4))
        pixLocs = np.zeros((events.shape[0],4))
        modCrysList[:,0] = getModSings(events[:,:4])
        modCrysList[:,2] = getModSings(events[:,4:8])
    else:

        modCrysList = np.zeros((events.shape[0],6))
        pixLocs = np.zeros((events.shape[0],6))
        modCrysList[:,0] = getModSings(events[:,:4])
        modCrysList[:,2] = getModSings(events[:,4:8])
        modCrysList[:,4] = getModSings(events[:,8:12])
    #singles
    pixInd = 0

    indMasterList = np.arange(modCrysList.shape[0])
    for iii in range(16):
        for jjj in range(4):
            #xLocsStore[iii,jjj,:],zLocsStore[iii,jjj,:] = getXZFromMandC(xRef,zRef,iii,jjj)
            x,z = getXZFromMandC(xRef,zRef,iii,jjj)
            xLocsStore[iii,jjj,:] = x .copy()
            zLocsStore[iii,jjj,:] = z.copy()
            
    for iii in range(16):
        q = events[modCrysList[:,0]==iii,:4].copy()
        listPart = indMasterList[modCrysList[:,0]==iii]
        c = getCrystal(q)
        for jjj in range(4):
            indList = listPart[c==jjj]
            modCrysList[indList,1] = jjj
            x,z = getXZFromMandC(xRef,zRef,iii,jjj)
            xLocsStore[iii,jjj,:] = x .copy()
            zLocsStore[iii,jjj,:] = z.copy()
            
    for iii in range(16):
        for jjj in range(4):
            indList = indMasterList[np.logical_and(modCrysList[:,0]==iii,modCrysList[:,1]==jjj)]
            xLocs = xLocsStore[iii,jjj,:].copy()
            xRes = np.mean(np.diff(xLocs))/2
            zLocs = zLocsStore[iii,jjj,:].copy()
            zRes = np.mean(np.diff(zLocs))/2

            #x
            pixOut[indList[events[indList,1]<=xLocs[0]],1] = xLocs[0]
            pixOut[indList[events[indList,1]>=xLocs[-1]],1] = xLocs[-1]
            pixLocs[indList[events[indList,0]<=xLocs[0]],0] = 0
            pixLocs[indList[events[indList,0]>=xLocs[-1]],0] = numPixX-1
            
            indSmall = indList[np.logical_and(events[indList,1]>xLocs[0],events[indList,1]<xLocs[-1])]
            for kkk in indSmall:
                pixOut[kkk,1] = xLocs[np.abs(events[kkk,1]-xLocs)<=xRes+wash][0]
                pixLocs[kkk,0] = (np.arange(numPixX)[np.abs(events[kkk,1]-xLocs)<=xRes+wash])[0]

            #z
            pixOut[indList[events[indList,3]<=zLocs[0]],3] = zLocs[0]
            pixOut[indList[events[indList,3]>=zLocs[-1]],3] = zLocs[-1]
            pixLocs[indList[events[indList,3]<=zLocs[0]],1] = 0
            pixLocs[indList[events[indList,3]>=zLocs[-1]],1] = numPixZ-1
            indSmall = indList[np.logical_and(events[indList,3]>zLocs[0],events[indList,3]<zLocs[-1])]
            for kkk in indSmall:
#                if len(zLocs[np.abs(events[kkk,3]-zLocs)<=zRes+wash]) > 1:
#                    print('Huh')
                pixOut[kkk,3] = zLocs[np.abs(events[kkk,3]-zLocs)<=zRes+wash][0]
                pixLocs[kkk,1] = (np.arange(numPixZ)[np.abs(events[kkk,3]-zLocs)<=zRes+wash])[0]

    #doubles
    if np.logical_or(events.shape[1]==9,events.shape[1]==13):
        for iii in range(16):
            q = events[modCrysList[:,2]==iii,4:8].copy()
            listPart = indMasterList[modCrysList[:,2]==iii]
            c = getCrystal(q)
            for jjj in range(4):
                indList = listPart[c==jjj]
                modCrysList[indList,3] = jjj
                
        for iii in range(16):
            for jjj in range(4):
                indList = indMasterList[np.logical_and(modCrysList[:,2]==iii,modCrysList[:,3]==jjj)]
                xLocs = xLocsStore[iii,jjj,:].copy()
                xRes = np.mean(np.diff(xLocs))/2
                zLocs = zLocsStore[iii,jjj,:].copy()
                zRes = np.mean(np.diff(zLocs))/2

                #x
                pixOut[indList[events[indList,5]<=xLocs[0]],5] = xLocs[0]
                pixOut[indList[events[indList,5]>=xLocs[-1]],5] = xLocs[-1]
                pixLocs[indList[events[indList,5]<=xLocs[0]],2] = 0
                pixLocs[indList[events[indList,5]>=xLocs[-1]],2] = numPixX-1
                indSmall = indList[np.logical_and(events[indList,5]>xLocs[0],events[indList,5]<xLocs[-1])]
                for kkk in indSmall:
                    pixOut[kkk,5] = xLocs[np.abs(events[kkk,5]-xLocs)<=xRes+wash][0]
                    pixLocs[kkk,2] = (np.arange(numPixX)[np.abs(events[kkk,5]-xLocs)<=xRes+wash])[0]

                #z
                pixOut[indList[events[indList,7]<=zLocs[0]],7] = zLocs[0]
                pixOut[indList[events[indList,7]>=zLocs[-1]],7] = zLocs[-1]
                pixLocs[indList[events[indList,7]<=zLocs[0]],3] = 0
                pixLocs[indList[events[indList,7]>=zLocs[-1]],3] = numPixZ-1
                indSmall = indList[np.logical_and(events[indList,7]>zLocs[0],events[indList,7]<zLocs[-1])]
                for kkk in indSmall:
                    pixOut[kkk,7] = zLocs[np.abs(events[kkk,7]-zLocs)<=zRes+wash][0]
                    pixLocs[kkk,3] = (np.arange(numPixZ)[np.abs(events[kkk,7]-zLocs)<=zRes+wash])[0]

    if events.shape[1]==13:

        for iii in range(16):
            q = events[modCrysList[:,4]==iii,8:12].copy()
            listPart = indMasterList[modCrysList[:,4]==iii]
            c = getCrystal(q)
            for jjj in range(4):
                indList = listPart[c==jjj]
                modCrysList[indList,5] = jjj
                
        for iii in range(16):
            for jjj in range(4):
                indList = indMasterList[np.logical_and(modCrysList[:,4]==iii,modCrysList[:,5]==jjj)]
                xLocs = xLocsStore[iii,jjj,:].copy()
                xRes = np.mean(np.diff(xLocs))/2
                zLocs = zLocsStore[iii,jjj,:].copy()
                zRes = np.mean(np.diff(zLocs))/2

                #x
                pixOut[indList[events[indList,9]<=xLocs[0]],9] = xLocs[0]
                pixOut[indList[events[indList,9]>=xLocs[-1]],9] = xLocs[-1]
                pixLocs[indList[events[indList,9]<=xLocs[0]],4] = 0
                pixLocs[indList[events[indList,9]>=xLocs[-1]],4] = numPixX-1
                indSmall = indList[np.logical_and(events[indList,9]>xLocs[0],events[indList,9]<xLocs[-1])]
                for kkk in indSmall:
                    pixOut[kkk,9] = xLocs[np.abs(events[kkk,9]-xLocs)<=xRes+wash][0]
                    pixLocs[kkk,4] = (np.arange(numPixX)[np.abs(events[kkk,9]-xLocs)<=xRes+wash])[0]

                #z
                pixOut[indList[events[indList,11]<=zLocs[0]],11] = zLocs[0]
                pixOut[indList[events[indList,11]>=zLocs[-1]],11] = zLocs[-1]
                pixLocs[indList[events[indList,11]<=zLocs[0]],5] = 0
                pixLocs[indList[events[indList,11]>=zLocs[-1]],5] = numPixZ-1
                indSmall = indList[np.logical_and(events[indList,11]>zLocs[0],events[indList,11]<zLocs[-1])]
                for kkk in indSmall:
                    pixOut[kkk,11] = zLocs[np.abs(events[kkk,11]-zLocs)<=zRes+wash][0]
                    pixLocs[kkk,5] = (np.arange(numPixZ)[np.abs(events[kkk,11]-zLocs)<=zRes+wash])[0]
    #random.shuffle(pixOut)
    return pixOut,modCrysList,pixLocs
    
def eCloudSigma(d,D=10,V=2000):
    '''
    sigma = eCloudSigma(d,E)
    calculates the sigma (assuming Guassian) of the electron cloud due to drifting
    from Feng's disseration pg 33. Assumes 293K
    input:
        E = electric field strength [V/m]
        d : drift distance [mm]
        D = detector thickness [mm]
    output:
        sigma : sigma, [mm]
    '''
    kTe = 0.0253 #[V]
    sigma = np.sqrt(2*kTe*d*D/V)
    return sigma
    
def landau(x,mu,c,a,b,d):
    '''
    output = landau(x,mu,c,a,b,d)
    
    returns a landau function based on the inputs
    
    '''
    xp = (x-mu)/c
    return np.exp(-(xp*a + np.exp(-xp*b))/d)
    
def addEneUn(doubles):
    '''
    doubEU = addEneUn(doubles)
    
    Adds energy uncertainy based on expected values for CZT
    
    input:
        doubles : input doubles without uncertainy
    output:
        doubEU : doubles with energy uncertainy
    '''
    doubEU = doubles.copy()
    if doubles.ndim > 1:
        s = doubles.shape
        sigma1 = eneUncCZT(doubles[:,0])/(doubles[:,0])
        randList1 = abs(random.normal(1,sigma1/2.35,doubles.shape[0]))
        doubEU[:,0] = doubEU[:,0]*randList1
        if s[1]>4:
            sigma2 = eneUncCZT(doubles[:,4])/(doubles[:,4])
            randList2 = abs(random.normal(1,sigma2/2.35,doubles.shape[0]))
            doubEU[:,4] = doubEU[:,4]*randList2
        if doubles.shape[1] >8:
            sigma3 = eneUnc(doubles[:,8])/(doubles[:,8])
            randList3 = abs(random.normal(1,sigma3/2.35,doubles.shape[0]))
            doubEU[:,8] = doubEU[:,8]*randList3
    else:
        sigma1 = eneUncCZT(doubles[0])/(doubles[0])
        sigma2 = eneUncCZT(doubles[4])/(doubles[4])
        randList1 = abs(random.normal(1,sigma1/2.35))
        randList2 = abs(random.normal(1,sigma2/2.35))
        doubEU[0] = doubEU[0]*randList1
        doubEU[4] = doubEU[4]*randList2
        if len(doubles)>8:
            sigma3 = eneUnc(doubles[8])/(doubles[8])
            randList3 = abs(random.normal(1,sigma3/2.35))
            doubEU[8] = doubEU[8]*randList3
    return doubEU
    
def eneUnc(ene):
    '''
    hold over from unupdated code. See eneUncCZT
    '''
    return eneUncCZT(ene)
    
def autoFullPixelateCommon(*events):
    '''
    outs = autoFullPixelateCommon(events)
    
    This function takes in either an all module set of singles, doubles, or triples and determines what the best
    reference pixels should be. 
    
    input:
        events : singles, doubles, or triples for all 16 modules
    output :
        outs : set of reference data
            sP = outs[0] : pixelated singles
            mcS = outs[1] : module and crystal indices of singles
            pxS = outs[2] : pixel indicies of singles
            dP = outs[3] : pixelated doubles
            mcD = outs[4] : module and crystal indices of doubles
            pxD = outs[5] : pixel indicies of doubles
            tP = outs[6] : pixelated triples
            mcT = outs[7] : module and crystal indices of triples
            pxT = outs[8] : pixel indicies of triples
            xRef = outs[9] : reference x values for common pixelation
            zRef = outs[10] : reference z vales for common pixelation
    '''
    
    haveSings = -1
    haveDoubs = -1
    haveTrips = -1
    
    for evList in events:
        if evList.shape[1] == 5 and evList.shape[0] > 0:
            haveSings = 1
            singles = evList
        if evList.shape[1] == 9 and evList.shape[0] > 0:
            haveDoubs = 1
            doubles = evList
        if evList.shape[1] == 13 and evList.shape[0] > 0:
            haveTrips = 1
            triples = evList
    xRefList = []
    zRefList = []
    
    if haveSings != -1:

        
        xRefS,zRefS = returnPixelLocs(singles)
        xRefS = xRefS[np.isfinite(xRefS)]
        zRefS = zRefS[np.isfinite(zRefS)]
        xRefList.append(xRefS)
        zRefList.append(zRefS)

    if haveDoubs != -1:
        
        xRefD1,zRefD1 = returnPixelLocs(doubles[:,:4])
        xRefD1 = xRefD1[np.isfinite(xRefD1)]
        zRefD1 = zRefD1[np.isfinite(zRefD1)]
        xRefList.append(xRefD1)
        zRefList.append(zRefD1)

        xRefD2,zRefD2 = returnPixelLocs(doubles[:,4:])
        xRefD2 = xRefD2[np.isfinite(xRefD2)]
        zRefD2 = zRefD2[np.isfinite(zRefD2)]
        xRefList.append(xRefD2)
        zRefList.append(zRefD2)

    if haveTrips != -1:
    
        xRefT1,zRefT1 = returnPixelLocs(triples[:,:4])
        xRefT1 = xRefT1[np.isfinite(xRefT1)]
        zRefT1 = zRefT1[np.isfinite(zRefT1)]
        xRefList.append(xRefT1)
        zRefList.append(zRefT1)

        xRefT2,zRefT2 = returnPixelLocs(triples[:,4:8])
        xRefT2 = xRefT2[np.isfinite(xRefT2)]
        zRefT2 = zRefT2[np.isfinite(zRefT2)]
        xRefList.append(xRefT2)
        zRefList.append(zRefT2)

        xRefT3,zRefT3 = returnPixelLocs(triples[:,:8:])
        xRefT3 = xRefT3[np.isfinite(xRefT3)]
        zRefT3 = zRefT3[np.isfinite(zRefT3)]
        xRefList.append(xRefT3)
        zRefList.append(zRefT3)
    
    print("\nDEBUG xRefList")
    
    print("len =", len(xRefList))
    
    for i, x in enumerate(xRefList):
        try:
            print(i, type(x), np.shape(x))
        except:
            print(i, type(x), "no shape")
    
    xRefList = np.array(xRefList)
   
    zRefList = np.array(zRefList)
    #print(xRefList)  
    xRef = np.mean(xRefList,axis=0)
    #print(xRef)
    xRefMin = np.abs(np.min(xRefList,axis=0)-xRef).max()
    xRefMax = np.abs(np.max(xRefList,axis=0)-xRef).max()
    
    zRef = np.mean(zRefList,axis=0)
    #print(zRef)
    zRefMin = np.abs(np.min(zRefList,axis=0)-zRef).max()
    zRefMax = np.abs(np.max(zRefList,axis=0)-zRef).max()
    
    wash = np.max([xRefMin,xRefMax,zRefMin,zRefMax,1E-3])
    #print(wash)
    
    outList = []
    if haveSings != -1:
        sP,mcS,pxS = pixelateCommon(singles,xRef,zRef,wash)
        outList.append(sP)
        outList.append(mcS)
        outList.append(pxS)
    if haveDoubs != -1:
        dP,mcD,pxD = pixelateCommon(doubles,xRef,zRef,wash)
        outList.append(dP)
        outList.append(mcD)
        outList.append(pxD)
    if haveTrips != -1:
        tP,mcT,pxT = pixelateCommon(triples,xRef,zRef,wash)
        outList.append(tP)
        outList.append(mcT)
        outList.append(pxT)
    outList.append(xRef)
    outList.append(zRef)
    return outList
    
def makeSingleMod(singlesA,doublesA,triplesA,outs,tSamp,nList,modNum):

    '''
    singOut,doubOut,tripOut,newOut,finalTSamp,newNList = makeSingleMod(singlesA,doublesA,triplesA,outs,tSamp,nList,modNum)
    
    this function takes in data from all 16 modules, and outputs the data only for the specified module. Genreally used when multiprocessing.
    
    input:
        singlesA, doublesA, triplesA : S, D, or T for all 16 modules
        outs : output of autoFullPixelateCommon for all 16 modules
        tSamp : deltaT list for all 16 modules
        nList : interaction ordering for all 16 modules
        modNum : module number of data to be extracted
    output:
        same as input, only for single specified module
    '''
    sP = outs[0]
    mcS = outs[1]
    pxS = outs[2]
    dP = outs[3]
    mcD = outs[4]
    pxD = outs[5]
    tP = outs[6]
    mcT = outs[7]
    pxT = outs[8]
    xRef = outs[9]
    zRef = outs[10]
    
    indS = mcS[:,0]==modNum
    
    indD = mcD[:,0]==modNum
    indD = np.logical_and(indD,mcD[:,2]==modNum)
    
    indDtoS1 = np.logical_and(mcD[:,0]==modNum,mcD[:,2]!=modNum)
    indDtoS2 = np.logical_and(mcD[:,0]!=modNum,mcD[:,2]==modNum)
    
    indT = mcT[:,0]==modNum
    indT = np.logical_and(indT,mcT[:,2]==modNum)
    indT = np.logical_and(indT,mcT[:,4]==modNum)
    
    indTtoS1 = np.logical_and(np.logical_and(mcT[:,0]==modNum,mcT[:,2]!=modNum),mcT[:,4]!=modNum)
    indTtoS2 = np.logical_and(np.logical_and(mcT[:,0]!=modNum,mcT[:,2]==modNum),mcT[:,4]!=modNum)
    indTtoS3 = np.logical_and(np.logical_and(mcT[:,0]!=modNum,mcT[:,2]!=modNum),mcT[:,4]==modNum)
    
    indTtoD12 = np.logical_and(np.logical_and(mcT[:,0]==modNum,mcT[:,2]==modNum),mcT[:,4]!=modNum)
    indTtoD23 = np.logical_and(np.logical_and(mcT[:,0]!=modNum,mcT[:,2]==modNum),mcT[:,4]==modNum)
    indTtoD13 = np.logical_and(np.logical_and(mcT[:,0]!=modNum,mcT[:,2]==modNum),mcT[:,4]!=modNum)
    
    tempTInds = np.arange(len(nList))
    finalTInds = tempTInds[nList==0][indS]
    finalTInds = np.append(finalTInds,tempTInds[nList==1][indD])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indT])
    finalTInds = np.append(finalTInds,tempTInds[nList==1][indDtoS1])
    finalTInds = np.append(finalTInds,tempTInds[nList==1][indDtoS2])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indTtoS1])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indTtoS2])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indTtoS3])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indTtoD12])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indTtoD23])
    finalTInds = np.append(finalTInds,tempTInds[nList==2][indTtoD13])
    finalTInds = np.unique(finalTInds)
    
    finalTSamp = np.append(tSamp[finalTInds][0],np.diff(np.cumsum(tSamp)[finalTInds]))
    
    sPNew = sP[indS,:]
    mcSNew = mcS[indS,:]
    pxSNew = pxS[indS,:]
    
    dPNew = dP[indD,:]
    mcDNew = mcD[indD,:]
    pxDNew = pxD[indD,:]
    
    tPNew = tP[indT,:]
    mcTNew = mcT[indT,:]
    pxTNew = pxT[indT,:]

    sPNewDtoS1 = dP[indDtoS1,:4]
    singDtoS1 = doublesA[indDtoS1,:4]
    mcSNewDtoS1 = mcD[indDtoS1,:2]
    pxSNewDtoS1 = pxD[indDtoS1,:2]
    
    sPNewDtoS2 = dP[indDtoS2,4:]
    singDtoS2 = doublesA[indDtoS2,4:]
    mcSNewDtoS2 = mcD[indDtoS2,2:]
    pxSNewDtoS2 = pxD[indDtoS2,2:]
    
    sPNewTtoS1 = tP[indTtoS1,:4]
    singTtoS1 = triplesA[indTtoS1,:4]
    mcSNewTtoS1 = mcT[indTtoS1,:2]
    pxSNewTtoS1 = pxT[indTtoS1,:2]

    sPNewTtoS2 = tP[indTtoS2,4:8]
    singTtoS2 = triplesA[indTtoS2,4:8]
    mcSNewTtoS2 = mcT[indTtoS2,2:4]
    pxSNewTtoS2 = pxT[indTtoS2,2:4]
    
    sPNewTtoS3 = tP[indTtoS3,8:]
    singTtoS3 = triplesA[indTtoS3,8:]
    mcSNewTtoS3 = mcT[indTtoS3,4:]
    pxSNewTtoS3 = pxT[indTtoS3,4:]
    
    dPNewTtoD12 = tP[indTtoD12,:8]
    doubTtoD12 = triplesA[indTtoD12,:8]
    mcDNewTtoD12 = mcT[indTtoD12,:4]
    pxDNewTtoD12 = pxT[indTtoD12,:4]

    dPNewTtoD23 = tP[indTtoD23,4:]
    doubTtoD23 = triplesA[indTtoD23,4:]
    mcDNewTtoD23 = mcT[indTtoD23,2:]
    pxDNewTtoD23 = pxT[indTtoD23,2:]

    dPNewTtoD13 = np.concatenate((tP[indTtoD13,:4],tP[indTtoD13,8:]),axis=1)
    doubTtoD13 = np.concatenate((triplesA[indTtoD13,:4],triplesA[indTtoD13,8:]),axis=1)
    mcDNewTtoD13 = np.concatenate((pxT[indTtoD13,:2],pxT[indTtoD13,4:]),axis=1)
    
    sCount = 0
    dCount = 0
    tCount = 0
    dtos1Count = 0
    dtos2Count = 0
    ttos1Count = 0
    ttos2Count = 0
    ttos3Count = 0
    ttod12Count = 0
    ttod23Count = 0
    ttod13Count = 0
    
    sPOut = []
    mcSOut = []
    pxSOut = []
    dPOut = []
    mcDOut = []
    pxDOut = []
    tPOut = []
    mcTOut = []
    pxTOut = []
    singOut = []
    doubOut = []
    tripOut = []
    newNList = []
    
    for iii in range(len(nList)):
        if nList[iii] == 0: #single
            if indS[sCount]:
                sPOut.append(sP[sCount,:])
                mcSOut.append(mcS[sCount,:])
                pxSOut.append(pxS[sCount,:])
                singOut.append(singlesA[sCount,:])
                newNList.append(0)
            sCount += 1
            
        elif nList[iii] == 1: #double
            if indD[dCount]: #true double
                dPOut.append(dP[dCount,:])
                mcDOut.append(mcD[dCount,:])
                pxDOut.append(pxD[dCount,:])
                doubOut.append(doublesA[dCount,:])
                newNList.append(1)
            elif indDtoS1[dCount]: #d1->s
                sPOut.append(sPNewDtoS1[dtos1Count])
                mcSOut.append(mcSNewDtoS1[dtos1Count])
                pxSOut.append(pxSNewDtoS1[dtos1Count])
                singOut.append(singDtoS1[dtos1Count])
                dtos1Count += 1
                newNList.append(0)
            elif indDtoS2[dCount]: #d1->s
                sPOut.append(sPNewDtoS2[dtos2Count])
                mcSOut.append(mcSNewDtoS2[dtos2Count])
                pxSOut.append(pxSNewDtoS2[dtos2Count])
                singOut.append(singDtoS2[dtos2Count])
                dtos2Count += 1
                newNList.append(0)

            dCount += 1
            
        elif nList[iii] == 2: #triple
            if indT[tCount]: #true triple
                tPOut.append(tP[tCount,:])
                mcTOut.append(mcT[tCount,:])
                pxTOut.append(pxT[tCount,:])
                tripOut.append(triplesA[tCount,:])
                newNList.append(2)
            elif indTtoS1[tCount]: #t1->s
                sPOut.append(sPNewTtoS1[ttos1Count])
                mcSOut.append(mcSNewTtoS1[ttos1Count])
                pxSOut.append(pxSNewTtoS1[ttos1Count])
                singOut.append(singTtoS1[ttos1Count])
                ttos1Count += 1
                newNList.append(0)
            elif indTtoS2[tCount]: #t2->s
                sPOut.append(sPNewTtoS2[ttos2Count])
                mcSOut.append(mcSNewTtoS2[ttos2Count])
                pxSOut.append(pxSNewTtoS2[ttos2Count])
                singOut.append(singTtoS2[ttos2Count])
                ttos2Count += 1
                newNList.append(0)
            elif indTtoS3[tCount]: #t3->s
                sPOut.append(sPNewTtoS3[ttos3Count])
                mcSOut.append(mcSNewTtoS3[ttos3Count])
                pxSOut.append(pxSNewTtoS3[ttos3Count])
                singOut.append(singTtoS3[ttos3Count])
                ttos3Count += 1
                newNList.append(0)
                
            elif indTtoD12[tCount]: #t12->d
                dPOut.append(dPNewTtoD12[ttod12Count,:])
                mcDOut.append(mcDNewTtoD12[ttod12Count,:])
                pxDOut.append(pxDNewTtoD12[ttod12Count,:])
                doubOut.append(doubTtoD12[ttod12Count,:])
                ttod12Count += 1
                newNList.append(1)
            elif indTtoD23[tCount]: #t23->d
                dPOut.append(dPNewTtoD23[ttod23Count,:])
                mcDOut.append(mcDNewTtoD23[ttod23Count,:])
                pxDOut.append(pxDNewTtoD23[ttod23Count,:])
                doubOut.append(doubTtoD23[ttod23Count,:])
                ttod23Count += 1
                newNList.append(1)
            elif indTtoD13[tCount]: #t13->d
                dPOut.append(dPNewTtoD13[ttod13Count,:])
                mcDOut.append(mcDNewTtoD13[ttod13Count,:])
                pxDOut.append(pxDNewTtoD13[ttod13Count,:])
                doubOut.append(doubTtoD13[ttod13Count,:])
                ttod13Count += 1
                newNList.append(1)
                
            tCount += 1
            

    sPOut = np.array(sPOut)
    mcSOut = np.array(mcSOut)
    pxSOut = np.array(pxSOut)
    dPOut = np.array(dPOut)
    mcDOut = np.array(mcDOut)
    pxDOut = np.array(pxDOut)
    tPOut = np.array(tPOut)
    mcTOut = np.array(mcTOut)
    pxTOut = np.array(pxTOut)
    singOut = np.array(singOut)
    doubOut = np.array(doubOut)
    tripOut = np.array(tripOut)
    newNList = np.array(newNList)
    
    newOut = [sPOut,mcSOut,pxSOut,dPOut,mcDOut,pxDOut,tPOut,mcTOut,pxTOut,xRef,zRef]
    
    return singOut,doubOut,tripOut,newOut,finalTSamp,newNList
    
def moduleSimulation(singles,doubles,triples,tSamp,nList,outs,resTime,deadScl,rReduction,q=0):
    
    '''
    newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow =  moduleSimulation(singles,doubles,triples,tSamp,nList,outs,q=0)
    
    The control function to rule them all! This is the main function that does the module simulation. There are a few specif
    
    '''

    #singles = singlesA
    #doubles = doublesA
    #triples = triplesA
    
    #q = 0

    #resTime = 1.5E-6
    trigThresh = 0.05
    tShape = 0.40E-6
    base_dir = Path(__file__).resolve().parent
    txt_dir = (base_dir/"../txt_in")
    muList = np.loadtxt(txt_dir/'muList_501x1001.txt')
    cList = np.loadtxt(txt_dir/'cList_501x1001.txt')    
    aList = np.loadtxt(txt_dir/'aList_501x1001.txt')
    bList = np.loadtxt(txt_dir/'bList_501x1001.txt')
    dList = np.loadtxt(txt_dir/'dList_501x1001.txt')
    gEnes = np.linspace(0.001,10,1001) #for 1001 enes
    
    depths = np.linspace(0.01,1,501)
    t = np.linspace(0,5E-6,5001)
    
    print('Determining outs')
   
    sP = outs[0]
    mcS = outs[1]
    pxS = outs[2]
    dP = outs[3]
    mcD = outs[4]
    pxD = outs[5]
    tP = outs[6]
    mcT = outs[7]
    pxT = outs[8]
    xRef = outs[9]
    zRef = outs[10]
 
    print('Finished dealing with outs')
    newDoubs = []
    newSings = []
    newTrips = []



    sCount = 0
    dCount = 0
    tCount = 0
    
    numSing = singles.shape[0]
    numDoub = doubles.shape[0]
    numTrip = triples.shape[0]
    totEvents = numSing + numDoub + numTrip
    #totEvents = 100000
    modNums = range(int(1+mcS[:,0].max()))
    modList = [module() for i in modNums]
    
    pixList = []
    iii = 0
    
    tDiff = 0
    coincCount = 0
    deadCount = 0
    tCurrent = tSamp[0]
    countOverflow = np.zeros((14,)) #[0] is true overflows from any source, [1] is fake doubles, [2] is fake triples, [3] is double to triple [4] is total dead time
    #[6] is detected interMod doubles, [7] is detected intermod triples, [8] is detected intercrys doubles, [9] is detected intercrys triples
    #[10] is ideal interMod doubles, [11] is ideal intermod triples, [12] is ideal intercrys doubles, [13] is ideal intercrys triples
    pixScl = []
    
    t = np.linspace(0,5E-6,501)
    tStep = np.mean(np.diff(t))
    tau = 200E-9
    num = [1/tau,0]
    den = [1,2/tau,1/(tau**2)]
    #a = signal.TransferFunction(num,den)
    
    p1Anode = singles[:,2].min()
    p2Anode = singles[:,2].max()
    
    if q!=0: #means multiprocessing
        crysThick1 = np.abs(p2Anode-p1Anode)
        crysThick2 = np.abs(p2Anode-p1Anode)
    else:
        unqSings = np.unique(singles[:,2])
        diffUS = np.diff(unqSings)
        thInd = diffUS>diffUS.mean()
    
        crysThick1 = unqSings[:-1][thInd]-p1Anode
        crysThick2 = p2Anode-unqSings[unqSings>unqSings[:-1][thInd]].min()
    
    tStart = time.time()
    plane1Anode = p1Anode
    plane2Anode = p2Anode
    V = 2000
    eLim = 10
    eChange = 10
    
  
    #totEvents = 10
    print('Total events = %i ' % totEvents)
    print(' ')
    
    while iii<totEvents:
        #print('*** ')
        #print('********** event number %i' % iii)
        #print(' ')        
        while tDiff > 0:
            tCurrent,tDiff = getNewTime(modList,pixList,tCurrent,tDiff,tShape)
            #print('*** event %i checkModStatus' % iii)
            pixList,modList,newSings,newDoubs,newTrips,countOverflow = processModStatusPulse(pixList,modList,tCurrent,newSings,newDoubs,newTrips,countOverflow,trigThresh,resTime,deadScl)
            modList = updateNextTimePulse(pixList,modList)
              
        tDiff = tSamp[iii]    
        if nList[iii] == 0: #single

            #print('*** event %i is a single' % iii)
            eInit = singles[sCount,4]
            mm = mcS[sCount,0]
            interMod = False
            if modList[mm.astype(int)].status <= 1:
                cc = mcS[sCount,1]
                xx = pxS[sCount,0]
                zz = pxS[sCount,1]
                ed = singles[sCount,0]
                rX = singles[sCount,1]
                rY = singles[sCount,2]
                rZ = singles[sCount,3]

                if ed > eLim:
                    ed = eChange
                
                if mm<=7:
                    V = 2300
                else:
                    V = 2600
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,eInit,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
                
                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
            else:
                deadCount += 1
                interMod = False
                interCrys = False
    
            sCount += 1
        elif nList[iii] == 1: #double


            interMod = mcD[dCount,0]!=mcD[dCount,2]
            interCrys = not(interMod) and (mcD[dCount,1]!=mcD[dCount,3])
            
            countOverflow[10] += int(interMod)
            countOverflow[12] += int(interCrys)
            
            eInit = doubles[dCount,8]
            mm = mcD[dCount,0]
            if modList[mm.astype(int)].status <= 1:
                cc = mcD[dCount,1]
                xx = pxD[dCount,0]
                zz = pxD[dCount,1]
                ed = doubles[dCount,0]
                rX = doubles[dCount,1]
                rY = doubles[dCount,2]
                rZ = doubles[dCount,3]
        
                if ed > eLim:
                    ed = eChange
                    
                if mm<=7:
                    V = 2300
                else:
                    V = 2600
                    
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,eInit,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
                #print('interaction 1 = %g' % eInit)                

                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
                if len(pOut) > 1 or (getDataPix(pOut,'eOut').min()<trigThresh):
                    interMod = False
            else:
                deadCount += 1
                interMod = False
                interCrys = False
            
            mm = mcD[dCount,2]
            if modList[mm.astype(int)].status <= 1:
                cc = mcD[dCount,3]
                xx = pxD[dCount,2]
                zz = pxD[dCount,3]
                ed = doubles[dCount,4]
                rX = doubles[dCount,5]
                rY = doubles[dCount,6]
                rZ = doubles[dCount,7]
                if ed > eLim:
                    ed = eChange   
                    
                if mm<=7:
                    V = 2300
                else:
                    V = 2600
    
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,eInit,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
                #print('interaction 2 = %g' % eInit)                 

                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
                if len(pOut) > 1 or (getDataPix(pOut,'eOut').min()<trigThresh):
                    interMod = False
            else:
                deadCount += 1
                interMod = False
                interCrys = False
                
            #print('len(pOut) %i' % len(pOut))
            dCount += 1
            countOverflow[6] += int(interMod)
            countOverflow[8] += int(interCrys)
        
            
        elif nList[iii] == 2: #triple

            #print('*** event %i is a triple' % iii)
            interMod = (mcT[tCount,0]!=mcT[tCount,2]) or (mcT[tCount,0]!=mcT[tCount,4]) or (mcT[tCount,2]!=mcT[tCount,4])
            ic12 = (mcT[tCount,1]!=mcT[tCount,3]) and (mcT[tCount,0]==mcT[tCount,2])
            ic13 = (mcT[tCount,1]!=mcT[tCount,5]) and (mcT[tCount,0]==mcT[tCount,4])
            ic23 = (mcT[tCount,3]!=mcT[tCount,5]) and (mcT[tCount,2]==mcT[tCount,4])
            interCrys = ic12 or ic13 or ic23
   
            countOverflow[11] += int(interMod)
            countOverflow[13] += int(interCrys)
            
            eInit = triples[tCount,12]
            mm = mcT[tCount,0]
            if modList[mm.astype(int)].status <= 1:
                cc = mcT[tCount,1]
                xx = pxT[tCount,0]
                zz = pxT[tCount,1]
                ed = triples[tCount,0]
                rX = triples[tCount,1]
                rY = triples[tCount,2]
                rZ = triples[tCount,3]
                if ed > eLim:
                    ed = eChange
                    
                if mm<=7:
                    V = 2300
                else:
                    V = 2600
                
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,eInit,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
                
                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
                if len(pOut) > 1 or (getDataPix(pOut,'eOut').min()<trigThresh):
                    interMod = False
            else:
                deadCount += 1
                interMod = False
                interCrys = False
                
            
            mm = mcT[tCount,2]
            if modList[mm.astype(int)].status <= 1:
                cc = mcT[tCount,3]
                xx = pxT[tCount,2]
                zz = pxT[tCount,3]
                ed = triples[tCount,4]
                rX = triples[tCount,5]
                rY = triples[tCount,6]
                rZ = triples[tCount,7]
                if ed > eLim:
                    ed = eChange
                    
                if mm<=7:
                    V = 2300
                else:
                    V = 2600
    
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,eInit,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
                
                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
                if len(pOut) > 1 or (getDataPix(pOut,'eOut').min()<trigThresh):
                    interMod = False
            else:
                deadCount += 1
                interMod = False
                interCrys = False
            
            mm = mcT[tCount,4]
            if modList[mm.astype(int)].status <= 1:
                cc = mcT[tCount,5]
                xx = pxT[tCount,4]
                zz = pxT[tCount,5]
                ed = triples[tCount,8]
                rX = triples[tCount,9]
                rY = triples[tCount,10]
                rZ = triples[tCount,11]
                if ed > eLim:
                    ed = eChange
                    
                if mm<=7:
                    V = 2300
                else:
                    V = 2600
                
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,eInit,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
                
                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
                if len(pOut) > 1 or (getDataPix(pOut,'eOut').min()<trigThresh):
                    interMod = False
            else:
                deadCount += 1
                interMod = False
                interCrys = False
                
            tCount += 1

            countOverflow[7] += int(interMod)
            countOverflow[9] += int(interCrys)
            
            
        modList = updateNextTimePulse(pixList,modList)
        
        iii += 1
        if iii % 1000000 == 0:
            print('starting event %i ' % iii)
        
        
    print('Putting variables')
    if q!= 0:   
        asdf = [newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow]
        q.put(asdf)
    print('Done with put')
    return newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow
    
def makeIndivModReadout(outs):
    '''
    
    '''

    mcS = outs[1]
    mcD = outs[4]
    mcT = outs[7]

    
    mcS[:,0] = 4*mcS[:,0]+mcS[:,1]+16
    mcS[:,1] = 0
    
    mcD[:,0] = 4*mcD[:,0]+mcD[:,1]+16
    mcD[:,2] = 4*mcD[:,2]+mcD[:,3]+16
    mcD[:,1] = 0
    mcD[:,3] = 0
    
    mcT[:,0] = 4*mcT[:,0]+mcT[:,1]+16
    mcT[:,2] = 4*mcT[:,2]+mcT[:,3]+16
    mcT[:,4] = 4*mcT[:,4]+mcT[:,5]+16
    mcT[:,1] = 0
    mcT[:,3] = 0
    mcT[:,5] = 0
    
    outs[1] = mcS
    outs[4] = mcD
    outs[7] = mcT
    
    return outs

# def findTripFromSandD(sings,tS,doubs,tD):
    # newT = []
    # n = 0
    # dInd = 0
    # sInd = 0
    # if tS[0]<tD[0]:
        # while tS[sInd]<tD[0]-1:
            # sInd += 1
        # sInd -= 1
    # else:
        # while tD[dInd]<tS[0]-1:
            # dInd += 1
# dInd -= 1
    # a,b = np.unique(tS,return_counts=True)
    # while (dInd < len(tD)-2) and (sInd < len(tS)-2):
        # if tS[sInd] <= tD[dInd]:
            # if tD[dInd] - tS[sInd] < 1: 
                # if b[tS[sInd]==a]==1:
                    # n += 1
                    # print('Found the trip %i at double %i, single %i' % (n,dInd,sInd))
                    # newT.append(np.concatenate((sings[sInd,:],doubs[dInd,:])))
            # sInd += 1
        # elif tS[sInd] > tD[dInd]:
            # if tS[sInd] - tD[dInd] < 1:
                # if b[tS[sInd]==a]==1:
                    # n += 1
                    # print('Found the trip %i at double %i, single %i' % (n,dInd,sInd))
                    # newT.append(np.concatenate((doubs[dInd,:],sings[sInd,:])))
            # dInd +=1
        # if tS[sInd] == tD[dInd]:
            # sInd += 1

    # return np.array(newT)
    
# def findDoubFromS(sings,tS,tThresh = 1):
    # newD = []
    # n = 0
    # tDiff = np.diff(tS)
    # sInd = 0
    # goodCount = 0
    # while sInd < len(tS)-1:
        # if tDiff[sInd] > tThresh:
            # if goodCount == 1:
                # newD.append(np.concatenate((sings[sInd,:],sings[sInd+1,:])))
                # n += 1
                # print('Found the doub %i at single %i' % (n,sInd))
            # goodCount = 0

        # else:
            # if tDiff[sInd] <tThresh:
                # goodCount += 1
        # sInd += 1
    # return np.array(newD)
    
    
# def findTripFromS(sings,tS,tThresh = 1):
    # newT = []
    # n = 0
    # tdiff = np.diff(tS)
    # for iii in range(1,len(tS)-1):
        # if np.logical_and(np.abs(tS[iii]-tS[iii-1])<tThresh,np.abs(tS[iii]-tS[iii+1])<tThresh):
            # newT.append(np.concatenate((sings[iii-1,:],sings[iii,:],sings[iii+1,:])))
            # n += 1
            # print('Found the trip %i at single %i' % (n,iii))
    # return np.array(newT)
    
def getDistance(x1,y1,z1,x2,y2,z2):
    return np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
    
def findIntermodDoubles(timeStamps,ene,xList,yList,zList,cutOffClk,peakMin,peakMax):
  ''' timeDiff,doubles = findIntermodDoubles(timeStamps,ene,xList,yList,zList,cutOffClk,peakMin,peakMax)
  
  This function takes in data from singles and tries to find a combination of singles that 
  occured within the cutoffclk window with an energy sum between peakMin and peakMax.
  It then outputs the time difference and the resultant double.
  
  input:
      timestamps in [s]
      ene in [MeV]
      xList,yList,zList self explainatory
      cutOffClk in [s]
      peakMin, peakMax in [MeV]
  output:
      timeDiff : time between singles that were used for inter-module event
      doubles : inter-mod doubles
  '''
  
  maxPts = np.floor(len(timeStamps)/2).astype(int)
  deltaT = np.zeros((maxPts,))
  yList1 = np.zeros((maxPts,))
  yList2 = np.zeros((maxPts,))
  ed1 = np.zeros((maxPts,))
  ed2 = np.zeros((maxPts,))
  xList1 = np.zeros((maxPts,))
  xList2 = np.zeros((maxPts,))
  zList1 = np.zeros((maxPts,))
  zList2 = np.zeros((maxPts,))
  jjj = 0
  timeDiff = timeStamps[1:]-timeStamps[0:-1]
  distThresh = 4
  foundLast = False
  for iii in range(len(timeDiff)-1):
    checkNum = 1
    buffTime = timeDiff[iii]
    if int(iii/100000)==float(iii/100000):
        print(iii)
    if (ene[iii] == 0) or (foundLast):
        foundLast = False
        continue
    while (buffTime <= cutOffClk) & ((checkNum+iii)<len(timeStamps)-3):
      #print(timeDiff[iii])
      eneSum = ene[iii]+ene[checkNum+iii]
    #  print('increased time')
    #  if sum(loc[iii+checkNum,:] == loc[iii,:])==3:
    #    coincNum = coincNum + 1
      if (eneSum <= peakMax) & (eneSum >= peakMin) & (getDistance(xList[iii],yList[iii],zList[iii],xList[iii+checkNum],yList[iii+checkNum],zList[iii+checkNum])>distThresh):
          deltaT[jjj] = buffTime
          ed1[jjj] = ene[iii]
          ed2[jjj] = ene[iii+checkNum]
          yList1[jjj] = yList[iii]
          yList2[jjj] = yList[checkNum+iii]
          xList1[jjj]=xList[iii]
          xList2[jjj]=xList[checkNum+iii]
          zList1[jjj]=zList[iii]
          zList2[jjj]=zList[checkNum+iii]
          jjj += 1
          foundLast = True
          
      checkNum = checkNum + 1
      buffTime = buffTime + timeDiff[iii+checkNum]
      
    #print(timeStamps[iii+checkNum] -timeStamps[iii])
    if jjj > maxPts-2:
      break
  wasdf = np.array((np.trim_zeros(ed1),np.trim_zeros(xList1),np.trim_zeros(yList1),np.trim_zeros(zList1),np.trim_zeros(ed2),np.trim_zeros(xList2),np.trim_zeros(yList2),np.trim_zeros(zList2))).T
  return np.trim_zeros(deltaT),wasdf
  
def findTripFromSandD(sings,tS,doubs,tD,peakMin,peakMax,tThresh=1E-6):
    newT = np.zeros((doubs.shape[0],12))
    #tInd = 0
    n = 0
    dInd = 0
    sInd = 0
    eSum2 = doubs[:,0]+doubs[:,4]
    eSum = 0
    if tS[0]<tD[0]:
        while tS[sInd]<tD[0]-tThresh:
            sInd += 1
        sInd -= 1
    else:
        while tD[dInd]<tS[0]-tThresh:
            dInd += 1
#dInd -= 1
    a,b = np.unique(tS,return_counts=True)
    while (dInd < len(tD)-2) and (sInd < len(tS)-2):
        if tS[sInd] <= tD[dInd]:
            if tD[dInd] - tS[sInd] < tThresh: 
                if b[tS[sInd]==a]==1:
                    eSum = eSum2[dInd] + sings[sInd,0]
                    if np.logical_and(eSum>=peakMin,eSum<=peakMax):
                        newT[n,:4] = sings[sInd,:].copy()
                        newT[n,4:] = doubs[dInd,:].copy()
                        n += 1
                        print('Found the trip %i at double %i, single %i' % (n,dInd,sInd))
                        #newT.append(np.concatenate((sings[sInd,:],doubs[dInd,:])))
            sInd += 1
        elif tS[sInd] > tD[dInd]:
            if tS[sInd] - tD[dInd] < tThresh:
                if b[tS[sInd]==a]==1:
                    eSum = eSum2[dInd] + sings[sInd,0]
                    if np.logical_and(eSum>=peakMin,eSum<=peakMax):
                        newT[n,:4] = sings[sInd,:].copy()
                        newT[n,4:] = doubs[dInd,:].copy()
                        n += 1
                        print('Found the trip %i at double %i, single %i' % (n,dInd,sInd))
                        #newT.append(np.concatenate((sings[sInd,:],doubs[dInd,:])))
            dInd +=1
        # if tS[sInd] == tD[dInd]:
            # sInd += 1

    return newT[~np.all(newT==0,axis=1)]
