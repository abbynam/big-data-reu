#@author: paul.maggi
#@edited: FS
#@adapted for chip: hlewis3

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
import scipy as sp
import time
import scipy.sparse as sprs
import bisect
from pathlib import Path

#np.random.seed(0) # [FS] Fixing randomization

cols = range(8)
## [FS] TODO: MAKE PIXEL CLASS TAKE 511 INDEX AS INPUT (AND NONE IF IT'S NOT A 511) (INPUT AS LIST)
# [FS] I modified the pixel class to take the 511 index as input (and NONE if
# not a 511). This can be a list if a event is composed of multiple 511 interactions.
class pixel: 
    '''
    The pixel is the basic unit of an interaction. It stores all relevant information about the position, energy, and temporal evolution of an event.
    '''
    def __init__(self,module=0,crystal=0,xInd=0,zInd=0,edep=0,rawX=0,rawY=0,rawZ=0,time=0,peakTrack=0,pulseShape=np.ones((501)),eScl=1,eOut=0,pulseLength=5E-6,tShapeThresh=0.001,trigThresh=0.05,resTime = 1.5E-6,A511Index=[]):
        self.module = int(module)
        self.crystal = crystal
        self.xInd = xInd
        self.zInd = zInd
        self.edep = edep
        self.rawX = rawX
        self.rawY = rawY
        self.rawZ = rawZ
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
        # [FS] 511 index
        if 'nay' in A511Index: # [FS] Keeping track of 511s
            self.A511Index = []
        else:
            self.A511Index = A511Index
        
        
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

# [FS] using my own bad file importer because I'm getting annoyed with numpy
def file_import(fName, skipToLine, delimiter = None):
    '''
    Imports files written in a tabular format (.csv, .tsv etc.). Can specify the delimiter for columns, and number of header lines to skip.
    
    Input:
        fName: self-explanatory
        skipToLine: number of lines to skip for file header.
        delimiter: specify character or short string separating individual values in a line of the file.
    Output:
        A 2D numpy array storing all columns of data in parallel 1D arrays.
    '''
    f = open(fName,"r" )
    i = 0
    j = 0
    arrays = []
    for line in f : # loopoverlines
        if j < skipToLine-1:
            j+=1
            continue
        if j == skipToLine-1:
            lineCheck = line.strip()
            columnCheck = line.split(delimiter)
            #print(len(columnCheck))
            arrayNum = len(columnCheck)
            for k in range(arrayNum):
                arrays.append([])
            j+=1
        line = line.strip()
        columns = line.split(delimiter)
        for k in range(len(columns)):
            arrays[k].append(float(columns[k]))
        i +=1
    #print([len(x) for x in arrays])
    #print([x for x in arrays])
    return np.asarray(arrays)

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
    sorDind = detM[:,9] == 0
    dInd = np.logical_and(np.invert(sInd),sorDind)
    newList = list(compress(range(len(dInd)),dInd))
    c = np.zeros((len(newList),8))
    c[range(len(newList)),:8]=detM[newList,:8]
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

# [FS] modified to use spline interpolant instead of fit function, so can take broader range of carbon-based materials
def protonStopPower(beamEne, stoppingPowerFunction,density):
    '''
    dEdX = protonStopPowerHDPE(beamEne)
    
    returns stopping power dE/dx [MeV/cm] in HDPE, assuming density = 0.97 g/cm^3
    
    Input:
        beamEne : proton beam energy, in [MeV]
        stoppingPowerFunction: stopping power function of target material (here expected to be a cubic spline interpolant to raw data)
        density: material density of target
    Output:
        dEdX : stopping power in [MeV/cm]
    ''' 
    msp = stoppingPowerFunction(beamEne)
    dEdX = msp*density
    return dEdX

# [FS] passing the stopping power function through, along with density and molar mass
def totalCrossSection(beamEne,mbList,eneList, stoppingPowerFunction, density, molarMass):
    '''
    interProb = totalCrossSection(beamEne,mbList,eneList)
    computes the total interaction probablity for a proton incident on HDPE based on the supplied
    cross sections and energies.
    input:
        beamEne : starting energy, should match a value in the eneList array
        mbList : list of differential (in energy) cross sections for a given reaction, in [mb]
        eneList : energy list in [MeV] corresponding to the data points in mbList
        stoppingPowerFunction: stopping power function of target material (here expected to be a cubic spline interpolant to raw data)
        density: material density of target
        molarMass: molar mass of target
    output:
        interProb : total interaction probability, using a CSDA approach from beamEne to stop.
    '''
    nDense = 2*density*(6.023E23)/molarMass #in units of atoms/cm^3
    interProb = 0
    mb = 1E-27 # 1 mb = 1E-27*cm^2
    dE = np.diff(eneList).mean()
    indStop = np.argmax(beamEne==eneList)
    for iii in range(indStop):
        dEdX  = protonStopPower(eneList[iii], stoppingPowerFunction,density) #MeV/cm
        dX = dE/dEdX #cm
        interProb += mb*mbList[iii]*dX*nDense * np.exp(-dX * mb * mbList[iii])
    return interProb

# [FS] I am modifying this method to calculate stopping power for graphite as well!
def count511More(beamEne,pPerS,t,num511,numProtons, materialType):
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
        materialType: string argument for the type of target material used. Currently, 'HDPE' and "graphite' are options.
    output:
        totalCounts : total number of 511 keV gammas that would be detected, assuming no detector timing effects
        COMMENTED OUT actXN : activities of the pathways of B8, C9, C10, C11, N12, N13
    '''
    
    # Importing stopping power data, as given by NIST. Also density (g/cm3) and molar masses (g/mol) of materials
    #-------------------------------------------------------------------------
    base_dir = Path(__file__).resolve().parent
    txt_dir = (base_dir/"../txt_in")
    if materialType == 'hdpe':
        fName = txt_dir/'stopping_power_hdpe.txt'
        density = 0.97
        molarMass = 28
    elif materialType == 'graphite':
        fName = txt_dir/'stopping_power_graphite.txt'
        density = 1.7
        molarMass = 12.011
    energy_StoppingPower = file_import(fName,skipToLine = 9) # imports stopping power as function of energy in MeV
    stoppingPowerFunction = sp.interpolate.interp1d(energy_StoppingPower[0], energy_StoppingPower[1], kind = 'cubic')
     #-------------------------------------------------------------------------
   
    eneList = np.linspace(10,250,241)
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
    
    prodProb8 = totalCrossSection(beamEne,mb8,eneList, stoppingPowerFunction, density, molarMass)
    prodProb9 = totalCrossSection(beamEne,mb9,eneList, stoppingPowerFunction, density, molarMass)
    prodProb10 = totalCrossSection(beamEne,mb10,eneList, stoppingPowerFunction, density, molarMass)
    prodProb11 = totalCrossSection(beamEne,mb11,eneList, stoppingPowerFunction, density, molarMass)
    prodProb12 = totalCrossSection(beamEne,mb12,eneList, stoppingPowerFunction, density, molarMass)
    prodProb13 = totalCrossSection(beamEne,mb13,eneList, stoppingPowerFunction, density, molarMass)
    totalProdProb = (prodProb8 + prodProb9 + prodProb10 + prodProb11 + prodProb12 + prodProb13)
    
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

# [FS] Modified this for two modules in standard configuration
def getModSings(singlesRaw): 
    '''
    mdls = getModSings(singlesRaw)
    
    this function takes in the data from a single set of interactions (edep,x,y,z) and returns the module number of each interaction.
    Assumes 2 modules. It can be passed any interaction quartet, e.g. (edep2,x2,y2,z2) from a triple.

    input: 
        singlesRaw : list of (edep,x,y,z) values
    output: 
        mdls : module numbers [0,15]
    '''
    singles = detPix(singlesRaw,0.1,0.1)
    xMean = np.unique(singles[:,1]).mean() #two options over 4 crystal locs
    zwischenzug = np.copy(singles[:,1])
    mods = np.where(zwischenzug < xMean, 0, 1)
    return mods
       
def final511Count(beamEne,pPerS,t,simProtons,num511, materialType = 'hdpe'):
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
        materialType: string argument for the type of target material used. Currently, 'HDPE' and "graphite' are options.
    output:
        final511 : total number of 511 keV gammas that would be detected, assuming no detector timing effects
    '''
    
    totalCounts = count511More(beamEne,pPerS,t,num511,simProtons, materialType)
    non511Rat = pPerS*t/simProtons
    final511 = totalCounts/non511Rat
    
    return final511
    
# [FS] returnPixelLocs() was modified to handle a single module with 2 or 4 crystals. All outer for loops - used to parse through the modules - were removed.
def returnPixelLocs(dropTwoCrystals, events, numPixX=11,numPixZ=11):
    
    '''
    xS,zS = returnPixelLocs(events,numPixX,numPixZ)
    
    This function returns the center locations of where pixels should be based on the supplied data. 
    It assumes a 2 module system, with 2 or 4 crystals, and each crystal has the specified number of pixels in X and Z.
    
    input:
        dropTwoCrystals: If False, assumes 2x2 grid of crystals. Else, assumes only top two crystals are being used.
        events : set of singles-like (edep, x, y, z) data for all 16 modules
        numPixX, numPixZ : number of pixels per crystal in X, Z respectively
    output:
        xS,zS : location of each pixel in each module and crystal.
            size of xS and zS is (16,4,numPix), meaning indexing is (module number, crystal number, pixel number)
    
    '''
    
    if dropTwoCrystals:
        crystalCount = 2
    else:
        crystalCount = 4
    
    wash = 0.0000001
    sideBuff = 0
    diffThresh = 0.2
    pixOut = events.copy()
    test = np.zeros(events.shape).copy()
    xLocsStore = np.zeros((2,crystalCount,numPixX)) # [FS] Changed to two modules /w 2 or 4 crystals
    zLocsStore = np.zeros((2,crystalCount,numPixZ)) # [FS] Changed to two modules /w 2 or 4 crystals
    modCrysList = np.zeros((events.shape[0],2))
    pixLocs = np.zeros((events.shape[0],2))
    modCrysList[:,0] = getModSings(events[:,:4]) # [FS] getModSings was modified for my standard two-module configuration
    #singles
    pixInd = 0

    indMasterList = np.arange(modCrysList.shape[0])
    for iii in range(2):    
        q = events[modCrysList[:,0]==iii,:4].copy()
        listPart = indMasterList[modCrysList[:,0]==iii]
        c = getCrystal(q, dropTwoCrystals)
        for jjj in range(crystalCount):
            indList = listPart[c==jjj]
            modCrysList[indList,1] = jjj
            temp = q[c==jjj,:].copy()

            #DEBUG
            print("i,j,shape,Xmin,Xmax,Zmin,Zmax")
            print(
                iii,
                jjj,
                temp.shape[0],
                temp[:,1].min(),
                temp[:,1].max(),
                temp[:,3].min(),
                temp[:,3].max()
            )

            xMin = temp[:,1].min()
            xMax = temp[:,1].max()
            xDiff = np.abs(xMin-xMax)
            zMin = temp[:,3].min()
            zMax = temp[:,3].max()
            zDiff = np.abs(zMin-zMax)
            xLocsStore[iii,jjj,:] = np.linspace(xMin+sideBuff,xMax-sideBuff,numPixX) # [FS] I replaced 11 with numPixX since I assume this is what they were going for
            zLocsStore[iii,jjj,:] = np.linspace(zMin+sideBuff,zMax-sideBuff,numPixZ) # [FS] I replaced 11 with numPixZ.
    
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

# [FS] processModStatusPulse was mostly modified to record annihilation gamma statistics, and record counts of interections in each module (including count overflows)
def processModStatusPulse(pixelList,modList,tNew,newSings,newDoubs,newTrips,A511Counts,countOverflow,A511Filter, mod0count,mod1count, trigThresh=0.05,resTime=1.5E-6,deadScl=1, eWindow = 0.015):
    '''
    pixelList,modList,newSings,newDoubs,newTrips,countOverflow = processModStatusPulse(pixelList,modList,tNew,newSings,newDoubs,newTrips,countOverflow,trigThresh=0.05,resTime=1.5E-6)
    
    The big boy! This is one of the the primary control functions for the detector simulation. It handles all module status-update calls, data cleanup, and output updating.
    Longer description TBD if I want to go into detail about the three module statuses, and what is done for each
    
    Pure input:
        tNew : current time that is checked against module status change times to determine which modules are changing
        trigThresh : trigger threshold for a pixel, in [MeV]
        resTime : active time after a trigger, [s]
        eWindow: simulated energy window width for identifying annihilation gammas (needs an better implementation)
    Input/output:
        pixelList, modList : list containing the module and pixel objects, respectively
        newSings, newDoubs, newTrips : list of output singles, doubles, and triples that have been processed, respectively
        A511Counts: an array that tracks annihilation gamma statistics:
            [0] is true overflows from any source, [1] is fake doubles, [2] is fake triples, [3] is double to triple [4] is total dead time
            [6] is detected interMod doubles, [7] is detected intermod triples, [8] is detected intercrys doubles, [9] is detected intercrys triples
            [10] is ideal interMod doubles, [11] is ideal intermod triples, [12] is ideal intercrys doubles, [13] is ideal intercrys triples
        A511Filter: boolean arrays parallel to newSings, newDoubs and newTrips that will specify which members of newSings, newDoubs and newTrips have the following associations with 511s:
            [0], [1], [2] give recorded events where 511s contributed (sing, doub, trip).
            [3], [4], [5] give true recorded 511 events (sing, doub, trip)
            [6], [7], [8] give detected 511 events (sing, doub, trip)
            [9], [10], [11] give true detected 511 events (sing, doub, trip)
        mod0count, mod1count: counters for the total number of interactions registered by the CC in each module.
        countOverflow : an array that tracks special interaction types
            countOverflow[0] : number of read out events that had more than 3 interactions (quads). They are not output, just counted
            countOverflow[1] : number of false doubles
            countOverflow[2] : number of false triples
            countOverflow[3] : number of D->Ts
            countOverflow[4] : total time spent dead [s]
            countOverflow[5] : currently not used
    '''
    
    modInds = np.arange(len(modList))[getDataMod(modList,'n')==tNew]
    for mI in modInds:
        modStatus = modList[mI].status
        if modStatus == 0: #means was in ready, check if new event triggers
            pI = []
            
            for iii in np.arange(len(pixelList)):
                p = pixelList[iii]
                if (p.module) == mI and (p.tStop == tNew):
                    if p.eOut >= trigThresh:
                        modList[mI].status = 1
                        #np.random.seed(0)
                        modList[mI].lastTime = p.tTrig + random.rand()*1E-10
                        modList[mI].nextTime = p.tTrig + resTime
                    else:
                        pI.append(iii)
                if modList[mI].status == 0:
                    modList[mI].nextTime = float('inf')
            for iii in pI[::-1]:
                del(pixelList[iii])
                        
        elif modStatus == 1: #means was in triggered,
            #store data with edep >= trigThresh and lastTime<=tpix<=nextTime
            #calc deadtime add
            #clear data, set to dead
            pI = []
            pix = []
            
            # [FS] Here, determining whether an interaction was a pure 511
            true511 = True
            if len(pixelList) == 1 and len(pixelList[0].A511Index) == 1:
                A511Counts[7] += 1
            elif len(pixelList) == 2 and len(pixelList[0].A511Index) == 1 and len(pixelList[1].A511Index) == 1 and pixelList[0].A511Index[0] == pixelList[1].A511Index[0]:
                A511Counts[8] += 1
            elif len(pixelList) == 3 and len(pixelList[0].A511Index) == 1 and len(pixelList[1].A511Index) == 1 and len(pixelList[2].A511Index) == 1 and pixelList[0].A511Index[0] == pixelList[1].A511Index[0] and pixelList[0].A511Index[0] == pixelList[2].A511Index[0]:
                A511Counts[9] += 1
            else:
                true511 = False
            
            for iii in np.arange(len(pixelList)):
                p = pixelList[iii]
                if (p.module == mI):
                    pI.append(iii)
                    if (p.eOut >= trigThresh) and (p.tMax >= modList[mI].lastTime) and (p.tMax <= modList[mI].nextTime):
                        pix.append(p)
                        
                        # [FS] Here I am counting the single, double and triple 511s that contribute to the final readings
                        
                        for item in p.A511Index:
                            if item[0]=='1':
                                A511Counts[0]+=1
                            elif item[0]=='2':
                                A511Counts[1]+=1
                            elif item[0]=='3':
                                A511Counts[2]+=1
                            for pixo in pixelList:
                                if item in pixo.A511Index:
                                    pixo.A511Index.remove(item)
                        
                        
            numCrys = len(np.unique(getDataPix(pix,'c')))
            if len(pix) > 0:
                if pix[0].module == 0:
                    mod0count += len(pix)
                elif pix[0].module == 1:
                    mod1count += len(pix)
            
            if len(pix) == 1 : #single                
                sing = np.array((pix[0].edep,pix[0].rawX,pix[0].rawY,pix[0].rawZ,modList[mI].lastTime,pix[0].peakTrack))
                newSings.append(sing)
                e = pix[0].edep
                # [FS] Checks whether energy adds up to 511keV (with specified energy window)
                if 0.511 - eWindow <= e and e <= 0.511  + eWindow:
                    A511Counts[10] += 1
                    A511Filter[6].append(1)
                    if true511:
                        A511Counts[13] += 1
                        A511Filter[9].append(1)
                    else:
                        A511Counts[16] += 1
                        A511Filter[9].append(0)
                else:
                    A511Filter[6].append(0)
                    A511Filter[9].append(0)
                # [FS] Adds to contributing 511 index
                if len(p.A511Index) != 0:
                    A511Filter[0].append(1)
                else:
                    A511Filter[0].append(0)
                if true511:
                    A511Filter[3].append(1)
                else:
                    A511Filter[3].append(0)
                
            elif len(pix) == 2: #double
                doub = []
                pt = []
                e = 0 # [FS] Finding total energy
                for p in pix:
                    doub.append(p.edep)
                    doub.append(p.rawX)
                    doub.append(p.rawY)
                    doub.append(p.rawZ)
                    pt.append(p.peakTrack)
                    e+=p.edep
                #[FS] Checks whether energy adds up to 511keV (with specified energy window)
                if 0.511 - eWindow <= e and e <= 0.511 + eWindow:
                    A511Counts[11] += 1
                    A511Filter[7].append(1)
                    if true511:
                        A511Counts[14] += 1
                        A511Filter[10].append(1)
                    else:
                        A511Counts[17] += 1
                        A511Filter[10].append(0)
                else:
                    A511Filter[7].append(0)
                    A511Filter[10].append(0)
                
                if len(p.A511Index) != 0:
                    A511Filter[1].append(1)
                else:
                    A511Filter[1].append(0)
                if true511:
                    A511Filter[4].append(1)
                else:
                    A511Filter[4].append(0)                
                doub.append(modList[mI].lastTime)
                if pix[0].time != pix[1].time: #false double, incr countOverflow[1]
                    countOverflow[1] += 1 
                    doub.append(1)
                else:
                    doub.append(0)

                doub.append(pt[0])
                doub.append(pt[1])
                    
                newDoubs.append(doub)
                
            elif len(pix) == 3:
                trip = []
                pt = []
                e = 0 # [FS] Finding total energy
                for p in pix:
                    trip.append(p.edep)
                    trip.append(p.rawX)
                    trip.append(p.rawY)
                    trip.append(p.rawZ)
                    pt.append(p.peakTrack)
                    e+=p.edep

                #[FS] Checks whether energy adds up to 511keV (with specified energy window)
                if 0.511 - eWindow <= e and e <= 0.511 + eWindow:
                    A511Counts[12] += 1
                    A511Filter[8].append(1)
                    if true511:
                        A511Counts[15] += 1
                        A511Filter[11].append(1)
                    else:
                        A511Counts[18] += 1
                        A511Filter[11].append(0)
                else:
                    A511Filter[8].append(0)
                    A511Filter[11].append(0)

                if len(p.A511Index) != 0:
                    A511Filter[2].append(1)
                else:
                    A511Filter[2].append(0)
                if true511:
                    A511Filter[5].append(1)
                else:
                    A511Filter[5].append(0)

                trip.append(modList[mI].lastTime)
                
                if (pix[0].time == pix[1].time) and (pix[0].time != pix[2].time):
                    countOverflow[3] += 1
                    trip.append(2)
                elif (pix[0].time != pix[1].time) and (pix[1].time == pix[2].time):
                    countOverflow[3] += 1
                    trip.append(3)
                elif (pix[0].time != pix[1].time) and (pix[1].time != pix[2].time) and (pix[0].time != pix[2].time): 
                    countOverflow[2] += 1
                    trip.append(1)
                else:
                    trip.append(0)
                trip.append(pt[0])    
                trip.append(pt[1])
                trip.append(pt[2])
                newTrips.append(trip)    
            elif len(pix)!=0:
                #true overflow, incr countOverflow[0]
                #print('Number of pix is %i' % len(pix))
                countOverflow[0] += 1
            #np.random.seed(0)
#            deadAdd = deadScl*(np.sum(121*(45E-9) + (19E-6)*len(pix) - 9E-6 + 2.5E6)) # [FS] just setting deadAdd to the expectation value, in the hopes of minimizing variance ...
            deadAdd = (np.sum(np.random.uniform(121*(40E-9),121*(50E-9),numCrys)) + (19E-6)*len(pix) + deadScl*1E-6 + np.random.uniform(0,5E-6)) # [FS] just doing with the tuning factor instead of the coefficient ...
            countOverflow[4] += deadAdd
            modList[mI].status = 2
            modList[mI].lastTime = tNew
            modList[mI].nextTime = tNew + deadAdd

            for pp in pI[::-1]:
                del(pixelList[pp])
                
        elif modStatus == 2: #means was in dead
            #change to ready
            modList[mI].status = 0
            modList[mI].lastTime = tNew
            modList[mI].nextTime = float('inf')
            
    return pixelList,modList,newSings,newDoubs,newTrips,countOverflow,A511Counts, A511Filter, mod0count, mod1count
    
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
       
# [FS] Modified for a 2-module system in standard configuration.
def getXZFromMandC(xL,zL,mNum,cNum, dropTwoCrystals):
    '''
    xOut,zOut = getXZFromMandC(xL,zL,mNum,cNum)
    returns the list of x's and z's from a given reference list based on module number and crystal number

    input:
        xL,zL : reference xList and zList for pixelated event locations; assumed unique list
        mNum,cNum : module number (0-2) and crystal number (0-2 or 0-4) 
        dropTwoCrystals: If False, assumes 2x2 grid of crystals. Else, assumes only top two crystals are being used.
    output:
        xOut,zOut : unique list of pixel locations in x and z, not grouped
    '''

    # [FS] I am replacing this with my own simplified version. Returns xOut and zOut corresponding to a crystal number of a module. Assumes 2 or 4 crystals.
    ##########################################################################
    xCenter = xL.mean()
    zCenter = zL.mean()
    
    xL1 = xL[xL < xCenter]
    xL2 = xL[xL >= xCenter]
    xCenter1 = xL1.mean()
    xCenter2 = xL2.mean()
    
    
    if dropTwoCrystals:
        zOut = zL
        if mNum == 0:
            if cNum == 0: # [FS] top 
                xOut = xL1[xL1 >= xCenter1]
            elif cNum == 1: # [FS] bottom
                xOut = xL1[xL1 < xCenter1]
        elif mNum == 1:
            if cNum == 0: # [FS] top
                xOut = xL2[xL2 >= xCenter2]
            elif cNum == 1: # [FS] bottom
                xOut = xL2[xL2 < xCenter2]
    else:
        if mNum == 0:
            if cNum == 0: #top right
                xOut = xL1[xL1 >= xCenter1]
                zOut = zL[zL >= zCenter]
            elif cNum == 1: #top left
                xOut = xL1[xL1 < xCenter1]
                zOut = zL[zL >= zCenter]
            elif cNum == 2: # bottom left
                xOut = xL1[xL1 < xCenter1]
                zOut = zL[zL < zCenter]
            elif cNum == 3: # bottom right
                xOut = xL1[xL1 >= xCenter1]
                zOut = zL[zL < zCenter]
        elif mNum == 1:
            if cNum == 0: #top right
                xOut = xL2[xL2 >= xCenter2]
                zOut = zL[zL >= zCenter]
            elif cNum == 1: #top left
                xOut = xL2[xL2 < xCenter2]
                zOut = zL[zL >= zCenter]
            elif cNum == 2: # bottom left
                xOut = xL2[xL2 < xCenter2]
                zOut = zL[zL < zCenter]
            elif cNum == 3: # bottom right
                xOut = xL2[xL2 >= xCenter2]
                zOut = zL[zL < zCenter]
    
    return xOut,zOut
    ##########################################################################
  
# [FS] Modified to only take one anode      
def findActivePixels(xTest,zTest,yTest,eDep,xRef,zRef,mInd,cInd,plane1Anode=60,V=2000,crysThick1=10,rReduction=1, dropTwoCrystals=4):
    '''
    eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(xTest,zTest,yTest,eDep,xRef,zRef,mInd,cInd,plane1Anode,plane2Anode,V)
    
    This function determines which of the four pixels nearest to an interaction receive some collected charge. It then reports the proper amount of "energy deposited"
    based on the volume of the charge received. This does include electron diffusion towards the anode, but does not include charge trapping.
    
    input:
        xTest,zTest,yTest : x, y, and z position of the interaction, respectively
        eDep : energy deposited [MeV]
        xRef,zRef : reference pixel center locations in x and z, respectively
        mInd,cInd : module and crystal index of the interaction, respectively
        plane1Anode: plane anode location
        V : applied detector voltage [V]
    output: 
        eDepPix : energy deposited in each pixel
        pixXInd,pixZInd : pixel indices for pixels that receive charge
        xC,zC : approximately locations of charge clouds for each pixel
    '''
    
    #if mInd >= 16:
    #    cInd = np.mod(mInd-16,4)
    #    mInd = np.floor((mInd-16)/4)
    xO,zO = getXZFromMandC(xRef,zRef,mInd,cInd, dropTwoCrystals)

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

    D = abs(crysThick1) # [FS] Modified
    d = plane1Anode-yTest # [FS] Modified  
    
    r = eCloudRad(eDep)+1.175*eCloudSigma(d,D,V)
    r = r/rReduction

    xP = xTest + r
    xM = xTest - r
    zP = zTest + r
    zM = zTest - r   

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

# [FS] Changed this to work for a SINGLE anode
def calcWeightingPotential(y,pA,crysThick=10):
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
    edepScale = 1 - (0.9808 * np.exp(-23.11*(pA-y)/crysThick))
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
        
# [FS] Modified to handle the case of a single anode plane, as well as annihilation gamma index information for the Pixel instance made in the method.
def makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,yTest,mInd,cInd,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList, I511, p1Anode=295,crysThick1=10): 
    '''
    pList = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,yTest,mInd,cInd,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,p1Anode=295,p2Anode=327.5)
    
    This function creates a list of pixel objects based on the output of findActivePixels. 
    
    input:
        eDepPix : list of energy deposited in a pixel. 
        pixXInd,pixZInd : list of x, z pixel indicies for each interaction
        xC,zC : list of x, z values for each interaction
        xRef,zRef : reference pixel location lists
        yTest : y value of interaction
        mInd,cInd : module and crystal indicies for the interaction
        t : time array to generate pulse shape over
        tCurrent : current time 
        gEnes,depths,muList,cList,aList,bList,dList : precomputed values for fitting temporal pulse shape
        p1Anode: anode location
        I511: list containing 511 indexing info
        crysThick1: crystal plane thickness
    output: 
        pList : list of pixel objects containing the input data. 
    '''
    
    pList = []
    mm = mInd
    cc = cInd
    rY = yTest
    # [FS] Changed weighting potential for case of one anode only
    escl = calcWeightingPotential(rY,p1Anode,crysThick1)
    
    # [FS] Modified for case of two modules, and one anode plane (actually not sure what this should be)
    dth = (p1Anode-rY)/crysThick1
    
    if len(pixXInd) == 1 and len(pixZInd) == 1: #single pixel
        xx = pixXInd
        zz = pixZInd
        ed = eDepPix
        rX = xC
        rZ = zC
        
        l = getPulseShape(t,dth,ed,gEnes,depths,muList,cList,aList,bList,dList)
        pixNew = pixel(module=mm,crystal=cc,xInd=xx,zInd=zz,edep=ed,rawX=rX,rawY=rY,rawZ=rZ,time=tCurrent,eScl=escl,pulseShape=l, A511Index=I511)
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
            pixNew = pixel(module=mm,crystal=cc,xInd=xx,zInd=zz,edep=ed,rawX=rX,rawY=rY,rawZ=rZ,time=tCurrent,eScl=escl,pulseShape=l, A511Index=I511)
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
                    pixNew = pixel(module=mm,crystal=cc,xInd=xx,zInd=zz,edep=ed,rawX=rX,rawY=rY,rawZ=rZ,time=tCurrent,eScl=escl,pulseShape=l,A511Index=I511)
                    pList.append(pixNew)        
    return pList

# [FS] modified to combine lists of 511 index lists when Pixel objects are combined.
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
    
    I511 = pix1.A511Index + pix2.A511Index #[FS] Here I am combining 511 IDs from both pixels, so we can keep track of all of them.
    
    if np.sum(newPulse)>0:
        pixComb = pixel(module=pix1.module,crystal=pix1.crystal,xInd=pix1.xInd,zInd=pix1.zInd,\
                    edep=eMax,rawX=xNew,rawY=yNew,rawZ=zNew,time=tMin,eScl=1,\
                    pulseShape=newPulse/newPulse.max(),pulseLength=(tMax-tMin),A511Index=I511)
    else:
        pixComb = pixel(module=pix1.module,crystal=pix1.crystal,xInd=pix1.xInd,zInd=pix1.zInd,\
            edep=eMax,rawX=xNew,rawY=yNew,rawZ=zNew,time=tMin,eScl=1,\
            pulseShape=np.zeros((len(tNew))),pulseLength=(tMax-tMin),A511Index=I511)
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
    
# [FS] this method is not needed in my implementation
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

# [FS] Adding a new section in this code to deal with the case of 2 crystals arranged horisontally.    
def getCrystal(eventsRaw, dropTwoCrystals):
    '''
    quad1[,quad2, quad3] = getCrystal(eventsRaw)
    
    returns which crystals an interaction took place in. Accepts singles, doubles, or triples. 
    NOTE: only input event data from a single module.
    
    input:
        eventsRaw : data to be analyzed in normal (edep,x,y,z) format. can be S, D or T.
        dropTwoCrystals: If False, assumes 2x2 grid of crystals. Else, assumes only top two crystals are being used.
    output:
        quadN : crystal of each interaction. Returns 1, 2, or 3 values.
    '''
    
    events = detPix(eventsRaw,0.1,0.1)
    ###########################################################################
    if dropTwoCrystals:
        if events.shape[1]==4:
            xMean = np.unique(events[:,1]).mean()
            
            R1 = events[:,1]>=xMean # [FS] Top crystal
            L1 = events[:,1]<xMean # [FS] Bottom crystal

            dual1 = np.zeros((len(R1),))
            
            dual1[R1] = 0
            dual1[L1] = 1
            
            return dual1.astype(int)
        
        elif events.shape[1]==8:
            xMean = np.unique(events[:,(1,5)]).mean()

            R1 = events[:,1]>=xMean # [FS] Top crystal
            L1 = events[:,1]<xMean # [FS] Bottom crystal
            R2 = events[:,5]>=xMean
            L2 = events[:,5]<xMean
            
            dual1 = np.zeros((len(R1),))
            dual2 = np.zeros((len(R1),))
            
            dual1[R1] = 0
            dual1[L1] = 1
            
            dual2[R2] = 0
            dual2[L2] = 1
            
            return dual2.astype(int)
 
        elif events.shape[1]==12:
            xMean = np.unique(events[:,(1,5,9)]).mean()
            
            R1 = events[:,1]>=xMean # [FS] Top crystal
            L1 = events[:,1]<xMean # [FS] Bottom crystal
            R2 = events[:,5]>=xMean
            L2 = events[:,5]<xMean
            R3 = events[:,9]>=xMean
            L3 = events[:,9]<xMean
            
            dual1 = np.zeros((len(R1),))
            dual2 = np.zeros((len(R1),))
            dual3 = np.zeros((len(R1),))
            
            dual1[R1] = 0
            dual1[L1] = 1
            
            dual2[R2] = 0
            dual2[L2] = 1
            
            dual3[R3] = 0
            dual3[L3] = 1
            
            return dual3.astype(int)
    ###########################################################################
    else:
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

# [FS] pixelateCommon was modified identically to returnPixelLocs.
def pixelateCommon(dropTwoCrystals, events,xRef,zRef,wash):
    '''
    pixOut,modCrysList,pixLocs = pixelateCommon(events,xRef,zRef,wash)
    
    pixelates the supplied event data to the reference data. 
    
    input:
        dropTwoCrystals: If False, assumes 2x2 grid of crystals. Else, assumes only top two crystals are being used.
        events : event data to be pixelated
        xRef, zRef : reference pixel centers
        wash : value needed for bit-rounding errors
    output:
        pixOut : pixelated data
        modCrysList : module numbers of interactions
        pixLocs : pixel index (x,z) of the interaction
    '''
    if dropTwoCrystals: #[FS] edited numPix for the cases of a single module with four or two crystals
        crystalCount = 2
        numPixX = int(len(xRef)/4) 
        numPixZ = int(len(zRef)/1)
    else:
        crystalCount = 4
        numPixX = int(len(xRef)/4) 
        numPixZ = int(len(zRef)/2)
    
    wash += 0.000001
    sideBuff = 0
    pixOut = events.copy()
    test = np.zeros(events.shape).copy()
    #xL,zL = xL,zL = returnPixelLocs(events[:,:4])
    
    xLocsStore = np.zeros((2,crystalCount,numPixX)) # [FS] reduced from 16 module down to 1 module.
    zLocsStore = np.zeros((2,crystalCount,numPixZ))
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
        modCrysList[:,4] = getModSings(events[:,8:])
    #singles
    pixInd = 0

    indMasterList = np.arange(modCrysList.shape[0])
    for iii in range(2):
        for jjj in range(crystalCount):
            
            x,z = getXZFromMandC(xRef,zRef,iii,jjj, dropTwoCrystals)
            xLocsStore[iii,jjj,:] = x.copy()
            zLocsStore[iii,jjj,:] = z.copy()
    
    for iii in range(2):
        q = events[modCrysList[:,0]==iii,:4].copy()
        listPart = indMasterList[modCrysList[:,0]==iii]
        c = getCrystal(q,dropTwoCrystals)
        for jjj in range(crystalCount):
            indList = listPart[c==jjj]
            modCrysList[indList,1] = jjj
            x,z = getXZFromMandC(xRef,zRef,iii,jjj, dropTwoCrystals)
            xLocsStore[iii,jjj,:] = x.copy()
            zLocsStore[iii,jjj,:] = z.copy()
    
    for iii in range(2):
        for jjj in range(crystalCount):
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
    if np.logical_or(events.shape[1]==8,events.shape[1]==12):
        for iii in range(2):
            q = events[modCrysList[:,2]==iii,4:8].copy()
            listPart = indMasterList[modCrysList[:,2]==iii]
            c = getCrystal(q, dropTwoCrystals)
            for jjj in range(crystalCount):
                indList = listPart[c==jjj]
                modCrysList[indList,3] = jjj
        
        for iii in range(2):
            for jjj in range(crystalCount):
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

    if events.shape[1]==12:
        
        for iii in range(2):
            q = events[modCrysList[:,4]==iii,8:].copy()
            listPart = indMasterList[modCrysList[:,4]==iii]
            c = getCrystal(q, dropTwoCrystals)
            for jjj in range(crystalCount):
                indList = listPart[c==jjj]
                modCrysList[indList,5] = jjj
        
        for iii in range(2):
            for jjj in range(crystalCount):
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
    calculates the sigma (assuming Gaussian) of the electron cloud due to drifting
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
        #np.random.seed(0)
        randList1 = abs(random.normal(1,sigma1/2.35,doubles.shape[0]))
        doubEU[:,0] = doubEU[:,0]*randList1
        if s[1]>4:
            sigma2 = eneUncCZT(doubles[:,4])/(doubles[:,4])
            #np.random.seed(0)
            randList2 = abs(random.normal(1,sigma2/2.35,doubles.shape[0]))
            doubEU[:,4] = doubEU[:,4]*randList2
        if doubles.shape[1] >8:
            sigma3 = eneUnc(doubles[:,8])/(doubles[:,8])
            #np.random.seed(0)
            randList3 = abs(random.normal(1,sigma3/2.35,doubles.shape[0]))
            doubEU[:,8] = doubEU[:,8]*randList3
    else:
        sigma1 = eneUncCZT(doubles[0])/(doubles[0])
        sigma2 = eneUncCZT(doubles[4])/(doubles[4])
        #np.random.seed(0)
        randList1 = abs(random.normal(1,sigma1/2.35))
        #np.random.seed(0)
        randList2 = abs(random.normal(1,sigma2/2.35))
        doubEU[0] = doubEU[0]*randList1
        doubEU[4] = doubEU[4]*randList2
        if len(doubles)>8:
            sigma3 = eneUnc(doubles[8])/(doubles[8])
            #np.random.seed(0)
            randList3 = abs(random.normal(1,sigma3/2.35))
            doubEU[8] = doubEU[8]*randList3
    return doubEU
    
def eneUnc(ene):
    '''
    hold over from unupdated code. See eneUncCZT
    '''
    return eneUncCZT(ene)

# [FS] Modified to handle the 2-module standard configuration    
def autoFullPixelateCommon(dropTwoCrystals, *events):
    '''
    outs = autoFullPixelateCommon(events)
    
    This function takes in either an all module set of singles, doubles, or triples and determines what the best
    reference pixels should be. 
    
    input:
        dropTwoCrystals: If False, assumes 2x2 grid of crystals. Else, assumes only top two crystals are being used.
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
        if evList.shape[1] == 4:
            haveSings = 1
            singles = evList
        if evList.shape[1] == 8:
            haveDoubs = 1
            doubles = evList
        if evList.shape[1] == 12:
            haveTrips = 1
            triples = evList
    xRefList = []
    zRefList = []
    
    if haveSings != -1:

        
        xRefS,zRefS = returnPixelLocs(dropTwoCrystals, singles) # [FS] These are edited to handle the standard two-module configuration
        print(f" xRefS len: {len(xRefS)}  zRefS len:{len(zRefS)}\n")
        xRefS = xRefS[np.isfinite(xRefS)]
        zRefS = zRefS[np.isfinite(zRefS)]
        xRefList.append(xRefS)
        zRefList.append(zRefS)

    if haveDoubs != -1:
        
        xRefD1,zRefD1 = returnPixelLocs(dropTwoCrystals, doubles[:,:4]) # [FS]
        print(f"xRefD1 len: {len(xRefD1)} zRefD1 len:{len(zRefD1)}\n")
        xRefD1 = xRefD1[np.isfinite(xRefD1)]
        zRefD1 = zRefD1[np.isfinite(zRefD1)]
        xRefList.append(xRefD1)
        zRefList.append(zRefD1)

        xRefD2,zRefD2 = returnPixelLocs(dropTwoCrystals, doubles[:,4:]) # [FS]
        print(f"xRefD2 len: {len(xRefD2)} zRefD2 len:{len(zRefD2)}\n")
        xRefD2 = xRefD2[np.isfinite(xRefD2)]
        zRefD2 = zRefD2[np.isfinite(zRefD2)]
        xRefList.append(xRefD2)
        zRefList.append(zRefD2)

    if haveTrips != -1:
        
        xRefT1,zRefT1 = returnPixelLocs(dropTwoCrystals, triples[:,:4]) # [FS]
        print(f"xRefT1 len: {len(xRefT1)} zRefT1 len:{len(zRefT1)}\n")
        xRefT1 = xRefT1[np.isfinite(xRefT1)]
        zRefT1 = zRefT1[np.isfinite(zRefT1)]
        xRefList.append(xRefT1)
        zRefList.append(zRefT1)

        xRefT2,zRefT2 = returnPixelLocs(dropTwoCrystals, triples[:,4:8]) # [FS]
        print(f"xRefT2 len: {len(xRefT2)} zRefT2 len:{len(zRefT2)}\n")
        xRefT2 = xRefT2[np.isfinite(xRefT2)]
        zRefT2 = zRefT2[np.isfinite(zRefT2)]
        xRefList.append(xRefT2)
        zRefList.append(zRefT2)

        xRefT3,zRefT3 = returnPixelLocs(dropTwoCrystals, triples[:,:8:]) # [FS]
        print(f"xRefT3 len: {len(xRefT3)} zRefT3 len:{len(zRefT3)}\n")
        xRefT3 = xRefT3[np.isfinite(xRefT3)]
        zRefT3 = zRefT3[np.isfinite(zRefT3)]
        xRefList.append(xRefT3)
        zRefList.append(zRefT3)
    
    print("\nDEBUG xRefList")

    print(f"\nxref list len: {len(xRefList)}\n")

    for i, x in enumerate(xRefList):
        try:
            print(i, type(x), np.shape(x))
        except:
            print(i, type(x), "no shape")

    print("\n", xRefList)

    xRefList = np.array(xRefList)
    
    print("\nDEBUG zRefList")

    print(f"\nzref list len: {len(zRefList)}\n")

    for j, z in enumerate(zRefList):
        try:
            print(j, type(z), np.shape(z))
        except:
            print(j, type(z), "no shape",)

    print("\n", zRefList)

    zRefList = np.array(zRefList)
    
    print(xRefList)
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
        sP,mcS,pxS = pixelateCommon(dropTwoCrystals, singles,xRef,zRef,wash) # As with returnPixelLocs(), the function was simplified, and module number was added as a parameter
        outList.append(sP)
        outList.append(mcS)
        outList.append(pxS)
    if haveDoubs != -1:
        dP,mcD,pxD = pixelateCommon(dropTwoCrystals, doubles,xRef,zRef,wash) # [FS]
        outList.append(dP)
        outList.append(mcD)
        outList.append(pxD)
    if haveTrips != -1:
        tP,mcT,pxT = pixelateCommon(dropTwoCrystals, triples,xRef,zRef,wash) # [FS]
        outList.append(tP)
        outList.append(mcT)
        outList.append(pxT)
    outList.append(xRef)
    outList.append(zRef)
    return outList
    
# [FS] I haven't touched this.
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

# Calls the following functions (note to self):
# getNewTime(modList,pixList,tCurrent,tDiff,tShape)
# processModStatusPulse(pixList,modList,tCurrent,newSings,newDoubs,newTrips,countOverflow,trigThresh,resTime,deadScl)
# updateNextTimePulse(pixList,modList)
# findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,plane2Anode,V,crysThick1,crysThick2,rReduction)
# makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,plane1Anode,plane2Anode,crysThick1,crysThick2)    
# checkNewEventPulse(pixList,modList,p,coincCount,deadCount)

# [FS] Modified through the processModStatusPulse() method, and some code to record intermodule events for annihilation gammas
def moduleSimulation(singles,doubles,triples,singlesA511Index,doublesA511Index,triplesA511Index,tSamp,nList,outs,resTime,deadScl,rReduction, dropTwoCrystals, q=0, printOut = True, sizeReduction = 1): # [FS] TODO: MAKE 511 INDICES, include input for makePixFromCloud
    
    '''
    newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow =  moduleSimulation(singles,doubles,triples,tSamp,nList,outs,q=0)
    
    The control function to rule them all! This is the main function that does the module simulation. There are a few specif
    
    '''
    #print("Running ModuleSim")
    #np.random.seed(0) # [FS] Fixing randomization
    #resTime = 1.5E-6
    trigThresh = 0.05
    tShape = 0.40E-6
    
    muList = np.loadtxt('muList_501x1001.txt')
    cList = np.loadtxt('cList_501x1001.txt')
    aList = np.loadtxt('aList_501x1001.txt')
    bList = np.loadtxt('bList_501x1001.txt')
    dList = np.loadtxt('dList_501x1001.txt')
    gEnes = np.linspace(0.001,10,1001) #for 1001 enes
    
    depths = np.linspace(0.01,1,501)
    t = np.linspace(0,5E-6,5001)
    
    sP = outs[0][:int(round(singles.shape[0]))]
    mcS = outs[1][:int(round(singles.shape[0]))]
    pxS = outs[2][:int(round(singles.shape[0]))]
    dP = outs[3][:int(round(doubles.shape[0]))]
    mcD = outs[4][:int(round(doubles.shape[0]))]
    pxD = outs[5][:int(round(doubles.shape[0]))]
    tP = outs[6][:int(round(triples.shape[0]))]
    mcT = outs[7][:int(round(triples.shape[0]))]
    pxT = outs[8][:int(round(triples.shape[0]))]
    xRef = outs[9]
    zRef = outs[10]

    newDoubs = []
    newSings = []
    newTrips = []
    
    sCount = 0
    dCount = 0
    tCount = 0
    
    s511Count = 0
    d511Count = 0
    t511Count = 0
    
    numSing = singles.shape[0]
    numDoub = doubles.shape[0]
    numTrip = triples.shape[0]
    totEvents = numSing + numDoub + numTrip
    percPoints = numpy.round(np.linspace(0, totEvents, 100))
    
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
    
    # [FS] 
    A511Counts = np.zeros((19,)) # [FS] means of keeping track of 511s. [0] is total contributing singles, [1] is total contributing doubles, [2] is total contributing triples, [3] is ideal intermod doubles, [4] is ideal intermod triples, [5] is detected intermod doubles, [6] is detected intermod triples,
                                # [7], [8] and [9] total true recorded 511 singles, doubles and triples, [10], [11], [12] detected 511 singles, doubles and triples, [13], [14], [15] true detected 511 singles, doubles and triples, [16], [17], [18] false detected singles, doubles and triples.
    A511Filters = [[],[],[],[],[],[],[],[],[],[],[],[]] # [FS] Here I initialise boolean arrays parallel to newSings, newDoubs and newTrips that will specify which members of newSings, newDoubs and newTrips have the following associations with 511s:
                # [0], [1], [2] give recorded events where 511s contributed (sing, doub, trip).
                # [3], [4], [5] give true recorded 511 events (sing, doub, trip)
                # [6], [7], [8] give detected 511 events (sing, doub, trip)
                # [9], [10], [11] give true detected 511 events (sing, doub, trip)
    mod0count = 0   
    mod1count = 0 # [FS] I'm getting a bit lazy, using this to indicate raw singles counts associated with each module
    
    t = np.linspace(0,5E-6,501)
    tStep = np.mean(np.diff(t))
    tau = 200E-9
    num = [1/tau,0]
    den = [1,2/tau,1/(tau**2)]
    #a = signal.TransferFunction(num,den)
    
    # [FS] Here, changed the code to only use p1Anode and crysThick1 since we only have one plane
    
    p1Anode = outs[0][:,2].max()
    depthref = outs[0][:,2].min()
    
    if q!=0: #means multiprocessing
        crysThick1 = np.abs(depthref-p1Anode)
    else:
        unqSings = np.unique(singles[:,2])
        diffUS = np.diff(unqSings)
        thInd = diffUS>diffUS.mean()
    
        crysThick1 = abs(depthref-p1Anode)
    
    tStart = time.time()
    plane1Anode = p1Anode
    
    V = 2000
    eLim = 10
    eChange = 10
    
    perc = -1 # [FS] keep track of percentage points
    
    if printOut: print("Module simulation progress:")
    
    while iii<totEvents:
        if iii % round(totEvents/100) == 0 and printOut:
            perc += 1
            sys.stdout.write("\r%i%%" % perc)
            sys.stdout.flush()
        
        while tDiff > 0:
            tCurrent,tDiff = getNewTime(modList,pixList,tCurrent,tDiff,tShape)
            pixList,modList,newSings,newDoubs,newTrips,countOverflow,A511Counts, A511Filter, mod0count,mod1count = processModStatusPulse(pixList,modList,tCurrent,newSings,newDoubs,newTrips,A511Counts,countOverflow,A511Filters,mod0count,mod1count,trigThresh,resTime,deadScl)
            modList = updateNextTimePulse(pixList,modList)
        tDiff = tSamp[iii]
        if nList[iii] == 0: #single
            #[FS] I check whether the interaction is a 511 using my index list, and then give the 511 a unique id
            #where the first digit gives the type (1 for single, 2 for double or 3 for triple) and the current count
            #is concatenated after. If it is not a 511, I give it a 'nay' label used when generating pixel objects 
            if singlesA511Index[sCount] == 1:
                A511Index = [('1' + str(sCount))]
            else:
                A511Index = ['nay']
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
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,V,crysThick1,rReduction, dropTwoCrystals) # [FS] Modified to handle the case of one plane, one anode only
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,A511Index,plane1Anode,crysThick1) # [FS] Modified to handle the case of one plane, one anode only. Also creates pixel objects that are labeled with the ID of the 511 they were triggered by.

                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount) # When pixel objects are combined in checkNewEventPulse, their 511 labels are also combined
            else:
                deadCount += 1
                interMod = False
                interCrys = False
    
            sCount += 1
        elif nList[iii] == 1: #double
            interMod = mcD[dCount,0]!=mcD[dCount,2]
            interCrys = not(interMod) and (mcD[dCount,1]!=mcD[dCount,3])

            if doublesA511Index[dCount] == 1: #[FS] Naming the 511
                A511Index = [("2" + str(dCount))]
            else:
                A511Index = ['nay']
            
            countOverflow[10] += int(interMod)
            if A511Index != 'nay': #[FS] counting intermods that are 511s
                A511Counts[3] += int(interMod)
            
            countOverflow[12] += int(interCrys)

                
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
                    
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,V,crysThick1,rReduction, dropTwoCrystals) # [FS] Modified to handle the case of one plane, one anode only
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,A511Index,plane1Anode,crysThick1) # [FS] Modified to handle the case of one plane, one anode only   
                
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
    
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,V,crysThick1,rReduction, dropTwoCrystals) # [FS] Modified to handle the case of one plane, one anode only
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,A511Index,plane1Anode,crysThick1) # [FS] Modified to handle the case of one plane, one anode only    
                
                for p in pOut:
                    pixList,coincCount,deadCount = checkNewEventPulse(pixList,modList,p,coincCount,deadCount)
                if len(pOut) > 1 or (getDataPix(pOut,'eOut').min()<trigThresh):
                    interMod = False
            else:
                deadCount += 1
                interMod = False
                interCrys = False
                
            dCount += 1
            countOverflow[6] += int(interMod)
            if A511Index != 'nay': #[FS] counting detected intermod doubles that are 511s
                A511Counts[5] += int(interMod)
            countOverflow[8] += int(interCrys)
            
        elif nList[iii] == 2: #triple
            interMod = (mcT[tCount,0]!=mcT[tCount,2]) or (mcT[tCount,0]!=mcT[tCount,4]) or (mcT[tCount,2]!=mcT[tCount,4])
            ic12 = (mcT[tCount,1]!=mcT[tCount,3]) and (mcT[tCount,0]==mcT[tCount,2])
            ic13 = (mcT[tCount,1]!=mcT[tCount,5]) and (mcT[tCount,0]==mcT[tCount,4])
            ic23 = (mcT[tCount,3]!=mcT[tCount,5]) and (mcT[tCount,2]==mcT[tCount,4])
            interCrys = ic12 or ic13 or ic23

            if triplesA511Index[tCount] == 1: #[FS] Naming the 511
                A511Index = [("3" + str(tCount))]
            else:
                A511Index = ['nay']    

            countOverflow[11] += int(interMod)
            if A511Index != 'nay': #[FS] counting ideal intermod triples that are 511s
                A511Counts[4] += int(interMod)
            countOverflow[13] += int(interCrys)
                
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
                
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,V,crysThick1,rReduction, dropTwoCrystals) # [FS] Modified to handle the case of one plane, one anode only
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,A511Index,plane1Anode,crysThick1) # [FS] Modified to handle the case of one plane, one anode only    
                
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
    
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,V,crysThick1,rReduction, dropTwoCrystals) # [FS] Modified to handle the case of one plane, one anode only
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,A511Index,plane1Anode,crysThick1) # [FS] Modified to handle the case of one plane, one anode only   
                
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
                
                eDepPix,pixXInd,pixZInd,xC,zC = findActivePixels(rX,rZ,rY,ed,xRef,zRef,mm,cc,plane1Anode,V,crysThick1,rReduction, dropTwoCrystals) # [FS] Modified to handle the case of one plane, one anode only
                pOut = makePixFromCloud(eDepPix,pixXInd,pixZInd,xC,zC,xRef,zRef,rY,mm,cc,t,tCurrent,gEnes,depths,muList,cList,aList,bList,dList,A511Index,plane1Anode,crysThick1) # [FS] Modified to handle the case of one plane, one anode only    
                
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
            if A511Index != 'nay': #[FS] counting ideal intermod triples that are 511s
                A511Counts[6] += int(interMod)
            countOverflow[9] += int(interCrys)
            
        modList = updateNextTimePulse(pixList,modList)
        
        iii += 1
    if q!= 0:   
        asdf = [newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow]
        q.put(asdf)
    return newSings,newDoubs,newTrips,coincCount,deadCount,countOverflow,A511Counts,A511Filter, mod0count, mod1count
    
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

# [FS] this filter separates data using the line y = x
def filter1(arg):
    return arg[:, 2] < 1.5*arg[:, 1]
# [FS] this filter separates data using the line y = x, but other way round
def filter2(arg):
    return arg[:, 2] >= 1.5*arg[:, 1]
# filter for rotating back
def filter3(arg, detectorwidth = 25):
    return arg[:,1] > detectorwidth

def filter4(arg, detectorwidth = 25):
    return np.invert(filter3(arg))

# [fS] Procedure for translating interaction data by a given amount in the x direction (can generalize later) (Takes singles-type data)
def translate(*args, xT = 0, yT = 0, zT = 0, filt_ = filter1):
    '''
    [fS] Procedure for translating interaction data by a given amount in the x direction (can generalize later) (Takes singles-type data)

    Parameters
    ----------
    *args : nparray
        takes an arbitrary number of singles-type data arrays
    xT : float, optional
        x translation. The default is 0.
    yT : TYPE, optional
        y translation. The default is 0.
    zT : TYPE, optional
        z translation. The default is 0.
    filt_ : function, optional
        data filter. The default is filter1.

    Returns
    -------
    None.

    '''
    for arg in args:
        filt = filt_(arg)
        arg[:, 1][filt] += xT
        arg[:, 2][filt] += yT
        arg[:, 3][filt] += zT

# [FS] Obsolete method, should be able to use translate to do the same thing
def translateOther(*args, xT = 0, yT = 0, zT = 0):
    for arg in args:
        filt = arg[:, 2] >= arg[:, 1]
        arg[:, 1][filt] += xT
        arg[:, 2][filt] += yT
        arg[:, 3][filt] += zT

# [FS] New rotation procedure, to rotate stuff in-place. Uses y=x line separator as condition (Takes singles-type data)
def rotate(*args, angle = 0, filt_ = filter1):
    '''
    [FS] Procedure to rotate interaction data in-place. Uses y=x line separator as condition (Takes singles-type data)

    Parameters
    ----------
    *args : nparray
        takes an arbitrary number of singles-type data arrays
    filt_ : function, optional
        data filter. The default is filter1.
    
    Returns
    -------
    None.        
    '''
    
    cosA = np.cos(angle)
    sinA = np.sin(angle)
    
    for arg in args:
        filt = filt_(arg)
        argX = np.copy(arg[:, 1][filt])
        argY = np.copy(arg[:, 2][filt])
        
        arg[:, 1][filt] = cosA*argX - sinA*argY
        arg[:, 2][filt] = sinA*argX + cosA*argY

# [FS] Obsolete method, should be able to use rotate to do the same thing
def rotateOther(*args, angle = 0):
    cosA = np.cos(angle)
    sinA = np.sin(angle)
    
    for arg in args:
        filt = arg[:, 2] >= 1.2*arg[:, 1]
        argX = np.copy(arg[:, 1][filt])
        argY = np.copy(arg[:, 2][filt])
        
        arg[:, 1][filt] = cosA*argX - sinA*argY
        arg[:, 2][filt] = sinA*argX + cosA*argY

def shuffle_in_unison(*arrs):
    '''
    [FS] Procedure for shuffling arbitrary number of arrays in the same order (https://stackoverflow.com/questions/4601373/better-way-to-shuffle-two-numpy-arrays-in-unison)
    
    Parameters
    ----------
    *arrs : nparray
        takes an arbitrary number of singles-type data arrays

    Returns
    -------
    None.

    '''
    #np.random.seed(0)
    rng_state = numpy.random.get_state()
    for arr in arrs:
        #np.random.seed(0)
        np.random.set_state(rng_state)
        #np.random.seed(0)
        np.random.shuffle(arr)
