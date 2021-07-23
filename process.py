#!/usr/bin/env python3
import pandas as pd
import numpy as np
from rotate import rotate_to_sector_0
import matplotlib.pyplot as plt
import ROOT
import time
import itertools
import random
import sys
from root_numpy import hist2array
import ctypes
import re

np.set_printoptions(threshold=sys.maxsize)
pd.set_option('display.max_rows', None)

infiles = []

def loadDataFile(MappingFile):

    column_names=['layer', 'u', 'v', 'density', 'shape', 'nDAQ', 'nTPG','DAQId1','nDAQeLinks1','DAQId2','nDAQeLinks2','TPGId1','nTPGeLinks1','TPGId2','nTPGeLinks2']

    #Read in data file
    data = pd.read_csv(MappingFile,delim_whitespace=True,names=column_names)

    #For a given module you need to know the total number of e-links
    data['TPGeLinkSum'] = data[['nTPGeLinks1','nTPGeLinks2']].sum(axis=1).where( data['nTPG'] == 2 , data['nTPGeLinks1'])

    #And find the fraction taken by the first or second lpgbt
    data['TPGeLinkFrac1'] = np.where( (data['nTPG'] > 0), data['nTPGeLinks1']/data['TPGeLinkSum'], 0  )
    data['TPGeLinkFrac2'] = np.where( (data['nTPG'] > 1), data['nTPGeLinks2']/data['TPGeLinkSum'], 0  )
    return data

def getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator):

    #Set number of TCs by hand for now (In reality this number is taken per event from CMSSW)
    column_names=[ 'u', 'v', 'layer', 'nTCs', 'nWords' ]
    data = pd.read_csv(CMSSW_Silicon,names=column_names)
    data_scin = pd.read_csv(CMSSW_Scintillator,names=column_names)

    return data,data_scin

def loadModuleTowerMappingFile(MappingFile):
    #File detailing which towers (in eta and phi coordinates)
    #are overlapped by which modules
    module_towermap = {}
    
    with open(MappingFile) as towerfile:

        for m,module in enumerate(towerfile):
            module_towers = []
            if (m==0): continue #header line
            modulesplit = module.split()
            moduleType = int(modulesplit[0])
            nTowers = modulesplit[4]
            towermaps = modulesplit[5:]

            countType = 0 #CE-E towers are separate from CE-H/Scintillator towers
            if moduleType > 0:
                countType = 1
            for tower in range (int(nTowers)):
                module_towers.append([countType, int(towermaps[3*tower]), int(towermaps[3*tower+1])])

            if moduleType == 2:
                module_towermap[1, int(modulesplit[1]), int(modulesplit[2]) , int(modulesplit[3])] = module_towers
            else:
                module_towermap[0, int(modulesplit[1]), int(modulesplit[2]) , int(modulesplit[3])] = module_towers

    return module_towermap

#Legacy function to read ROverZHistograms file with 1D histograms
def getModuleHists1D(HistFile):

    module_hists = []
    inclusive_hists = []
    
    infiles.append(ROOT.TFile.Open(HistFile,"READ"))

    inclusive = {}
    phi60 = {}

    list_of_hists = [hist.GetName() for hist in infiles[-1].GetListOfKeys()]

    for hist in list_of_hists:
        if "Inclusive" in hist or "ROverZ" not in hist:
            continue
        
        ijk = hist.split("_")
        
        if "silicon" in hist:
            if "Phi60" in hist:
                phi60[0,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = infiles[-1].Get(hist)
            else:
                inclusive[0,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = infiles[-1].Get(hist)
        elif "scintillator" in hist:
            if "Phi60" in hist:
                phi60[1,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = infiles[-1].Get(hist)
            else:
                inclusive[1,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = infiles[-1].Get(hist)

    inclusive_hists.append(infiles[-1].Get("ROverZ_Inclusive" ))
    inclusive_hists.append(infiles[-1].Get("ROverZ_Inclusive_Phi60" ))
                
    module_hists.append(inclusive)
    module_hists.append(phi60)
            
    return inclusive_hists,module_hists

#Function to read TC histograms file with 1D histograms
def getModuleTCHists(HistFile):

    module_hists = {}
    
    infiles.append(ROOT.TFile.Open(HistFile,"READ"))
    list_of_hists = [hist.GetName() for hist in infiles[-1].GetListOfKeys()]
    
    for hist in list_of_hists:
        if "nTCs" not in hist:
            continue

        ijk = hist.split("_")

        if "silicon" in hist:
            module_hists[0,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = infiles[-1].Get(hist)
        elif "scintillator" in hist:
            module_hists[1,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = infiles[-1].Get(hist)

    return module_hists


def getPhiSplitIndices( PhiVsROverZ, split = "fixed", fixvalue = 55):
    #If split is a fixed value in each R/Z bin, the fixvalue can be set, otherwise it is ignored
    #Fix value should be a multiple of 5/3 degrees, but if not is anyway forced to the closest multiple

    PhiVsROverZ_profile = PhiVsROverZ.ProfileX()
    
    if ( split == "per_roverz_bin" ):
        phi_divisions = np.empty(0)
        for b in range(1, PhiVsROverZ_profile.GetNbinsX()+1):#begin from first bin
            phi_divisions = np.append(phi_divisions,PhiVsROverZ_profile.GetBinContent(b))
    elif ( split == "fixed" ):
        phi_divisions = np.full(PhiVsROverZ_profile.GetNbinsX(), np.radians(fixvalue))

    #The Zeroth element is the bin 1 low edge
    edges = np.zeros(PhiVsROverZ.GetNbinsY()+1,dtype='double')
    PhiVsROverZ.GetYaxis().GetLowEdge(edges) 

    split_indices = []

    for value in phi_divisions:
        split_indices.append( (np.abs(edges - value)).argmin() + 1 )# +1 so that an element is converted into a ROOT bin number

    return split_indices
    
def getModuleHists(HistFile, split = "fixed", phidivisionX_fixvalue_min = 55, phidivisionY_fixvalue_max = None):
    #phidivisionX and phidivisionY refer to the traditional regions of phi>60 degrees and phi<60 degrees respectively.
    #Renamed such that their definitions can be more flexible

    #Split in phi is either 'fixed' (set at fixvalue in degrees)
    #or 'per_roverz_bin' in which case the split will be taken
    #at the mean point in phi in each R/Z bin (inclusive over all events)

    #If split is 'fixed' phidivisionX is lower-bounded by 'phidivisionX_fixvalue_min', with no upper bound
    #phidivisionY is upper-bounded by 'phidivisionY_fixvalue_max'
    #(or 'phidivisionX_fixvalue_min' if 'phidivisionY_fixvalue_max' is not set), with no lower bound
    if ( phidivisionY_fixvalue_max == None ):
        phidivisionY_fixvalue_max = phidivisionX_fixvalue_min
    
    module_hists = []
    inclusive_hists = []
    
    infiles.append(ROOT.TFile.Open(HistFile,"READ"))

    if not infiles[-1]:
        raise EnvironmentError
    
    PhiVsROverZ = infiles[-1].Get("ROverZ_Inclusive" )

    nBinsPhi = PhiVsROverZ.GetNbinsY()    

    #Get the phi indices where the split between "high" and "low" phi should be
    #Could either be constant with R/Z, or different for each R/Z bin 

    #Gives the bin number, whose low edge is the splitting point
    split_indices_DivisionX = getPhiSplitIndices( PhiVsROverZ, split = split, fixvalue = phidivisionX_fixvalue_min)
    split_indices_DivisionY = getPhiSplitIndices( PhiVsROverZ, split = split, fixvalue = phidivisionY_fixvalue_max)

    projectionX_PhiDivisionX = PhiVsROverZ.ProjectionX( "ROverZ_PhiDivisionX" )
    projectionX_PhiDivisionY = PhiVsROverZ.ProjectionX( "ROverZ_PhiDivisionY" )
    projectionX_PhiDivisionX.Reset()
    projectionX_PhiDivisionY.Reset()

    #Get an independent projection for each R/Z bin
    for x in range(1,PhiVsROverZ.GetNbinsX()+1):
        error = ctypes.c_double(-1)
        projectionX_PhiDivisionX.SetBinContent(x,PhiVsROverZ.IntegralAndError(x,x,int(split_indices_DivisionX[x-1]),int(nBinsPhi),error))
        projectionX_PhiDivisionX.SetBinError(x,error.value)
        projectionX_PhiDivisionY.SetBinContent(x,PhiVsROverZ.IntegralAndError(x,x,1,int(split_indices_DivisionY[x-1]-1),error))
        projectionX_PhiDivisionY.SetBinError(x,error.value)

    inclusive_hists.append(projectionX_PhiDivisionX)
    inclusive_hists.append(projectionX_PhiDivisionY)
    
    phiDivisionX = {}
    phiDivisionY = {}

    list_of_hists = [hist.GetName() for hist in infiles[-1].GetListOfKeys()]

    for hist in list_of_hists:
        if "Inclusive" in hist or "ROverZ" not in hist:
            continue

        ijk = hist.split("_") #To get u (ieta), v (iphi), and layer
        PhiVsROverZ = infiles[-1].Get(hist)
        nBinsPhi = PhiVsROverZ.GetNbinsY()        

        if "silicon" in hist:                            
            projectionX_PhiDivisionX = PhiVsROverZ.ProjectionX( "ROverZ_silicon_"+ijk[-3]+"_"+ijk[-2]+"_"+ijk[-1]+"_PhiDivisionX" )
            projectionX_PhiDivisionY = PhiVsROverZ.ProjectionX( "ROverZ_silicon_"+ijk[-3]+"_"+ijk[-2]+"_"+ijk[-1]+"_PhiDivisionY" )
        elif "scintillator" in hist:
            projectionX_PhiDivisionX = PhiVsROverZ.ProjectionX( "ROverZ_scintillator_"+ijk[-3]+"_"+ijk[-2]+"_"+ijk[-1]+"_PhiDivisionX" )
            projectionX_PhiDivisionY = PhiVsROverZ.ProjectionX( "ROverZ_scintillator_"+ijk[-3]+"_"+ijk[-2]+"_"+ijk[-1]+"_PhiDivisionY" )

        projectionX_PhiDivisionX.Reset()
        projectionX_PhiDivisionY.Reset()

        #Get an independent projection for each R/Z bin
        for x in range(1,PhiVsROverZ.GetNbinsX()+1):
            error = ctypes.c_double(-1)
            projectionX_PhiDivisionX.SetBinContent(x,PhiVsROverZ.IntegralAndError(x,x,int(split_indices_DivisionX[x-1]),int(nBinsPhi),error))
            projectionX_PhiDivisionX.SetBinError(x,error.value)
            projectionX_PhiDivisionY.SetBinContent(x,PhiVsROverZ.IntegralAndError(x,x,1,int(split_indices_DivisionY[x-1]-1),error))
            projectionX_PhiDivisionY.SetBinError(x,error.value)
        
        if "silicon" in hist:
            phiDivisionX[0,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = projectionX_PhiDivisionX
            phiDivisionY[0,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = projectionX_PhiDivisionY
        elif "scintillator" in hist:
            phiDivisionX[1,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = projectionX_PhiDivisionX
            phiDivisionY[1,int(ijk[-3]),int(ijk[-2]),int(ijk[-1])] = projectionX_PhiDivisionY

    module_hists.append(phiDivisionX)
    module_hists.append(phiDivisionY)
            
    return inclusive_hists,module_hists

def getHistsPerLayer(module_hists):

    f = ROOT.TFile.Open("HistsPerLayer.root","RECREATE")
    histsPerLayer = []
    
    #Get binning from module hists
    example_hist = next(iter(module_hists[0].values()))
    nROverZBins = example_hist.GetNbinsX()
    rOverZMin = example_hist.GetXaxis().GetBinLowEdge(1)
    rOverZMax = example_hist.GetXaxis().GetBinLowEdge(nROverZBins + 1)

    for layer in range(1,53):
        layer_hist = ROOT.TH1D( ("layer_ROverZ" + "_" + str(layer)),"",nROverZBins,rOverZMin,rOverZMax)
        histsPerLayer.append(layer_hist)

    print (len(histsPerLayer))

    #phi > 60
    for key,hist in module_hists[0].items():
        #key[3] is the layer number (which goes from 1-52)
        histsPerLayer[key[3]-1].Add( hist )

    #phi < 60
    for key,hist in module_hists[1].items():
        #key[3] is the layer number (which goes from 1-52)
        histsPerLayer[key[3]-1].Add( hist )

    for hist in histsPerLayer:
        hist.Write()

    f.Close()
        
    return histsPerLayer

def getlpGBTLoadInfo(data,data_tcs_passing,data_tcs_passing_scin):
    #Loop over all lpgbts
        
    lpgbt_loads_tcs=[]
    lpgbt_loads_words=[]
    layers=[]

    maxLpGBTNumber = data[['TPGId1','TPGId2']].max().max()
    
    for lpgbt in range(0,maxLpGBTNumber+1):

        tc_load = 0.
        word_load = 0.
        lpgbt_layer = -1.
        lpgbt_found = False
        
        for tpg_index in ['TPGId1','TPGId2']:#lpgbt may be in the first or second position in the file

            for index, row in (data[data[tpg_index]==lpgbt]).iterrows():  
                passing = data_tcs_passing
                if (row['density']==2):#Scintillator
                    passing = data_tcs_passing_scin

                lpgbt_found = True
                lpgbt_layer = row['layer']
                #Get the number of TCs passing threshold
                line = passing[(passing['layer']==row['layer'])&(passing['u']==row['u'])&(passing['v']==row['v'])]
                nwords = -1
                ntcs = -1
                if not line.empty: 
                    nwords = line['nWords'].values[0]
                    ntcs = line['nTCs'].values[0]

                    linkfrac = 'TPGeLinkFrac1'
                    if ( tpg_index == 'TPGId2' ):
                        linkfrac = 'TPGeLinkFrac2'
                    tc_load += row[linkfrac] * ntcs #Add the number of trigger cells from a given module to the lpgbt
                    word_load += row[linkfrac] * nwords #Add the number of trigger cells from a given module to the lpgbt
                    
        if (lpgbt_found):
            lpgbt_loads_tcs.append(tc_load)

            lpgbt_loads_words.append(word_load)
            layers.append( lpgbt_layer )

    return lpgbt_loads_tcs,lpgbt_loads_words,layers

def getHexModuleLoadInfo(data,data_tcs_passing,data_tcs_passing_scin,print_modules_no_tcs=False):

    module_loads_words=[]
    layers = []
    u_list = []
    v_list = []

    for index, row in data.iterrows():
        if ( row['TPGeLinkSum'] < 0 ):
            continue
        passing = data_tcs_passing
        if (row['density']==2):
            passing = data_tcs_passing_scin
        
        module_layer = row['layer']
        line = passing[(passing['layer']==row['layer'])&(passing['u']==row['u'])&(passing['v']==row['v'])]

        nwords = -1
        ntcs = -1
        #if not line.empty: 
        nwords = line['nWords'].values[0]
        ntcs = line['nTCs'].values[0]
        if (print_modules_no_tcs):
            if ( ntcs == 0 ):
                if (row['density']<2):
                    print ("sili",row['layer'],row['u'],row['v'])
                if (row['density']==2):
                    print ("scin",row['layer'],row['u'],row['v'])
        
        word_load = nwords / (2 * row['TPGeLinkSum'] )

        module_loads_words.append(word_load)
        layers.append(module_layer)
        u_list.append(row['u'])
        v_list.append(row['v'])

    return module_loads_words,layers,u_list,v_list

#Minimal custom implementation of root_numpy hist2array returning errors
#See https://github.com/scikit-hep/root_numpy/blob/master/root_numpy/_hist.py
# for original
def TH1D2array(hist, include_overflow=False, return_error_squares=False):

    # Determine dimensionality and shape
    shape = (hist.GetNbinsX() + 2)
    dtype = np.dtype('f8')
    array = np.ndarray(shape=shape, dtype=dtype, buffer=hist.GetArray())
    errors = np.ndarray(shape=shape, dtype=dtype)

    for e in range (hist.GetNbinsX() + 1):
        errors[e] = np.power(hist.GetBinError(e),2)
    
    if not include_overflow:
        # Remove overflow and underflow bins
        array = array[slice(1, -1)]
        errors = errors[slice(1, -1)]

    array = np.copy(array)
    errors = np.copy(errors)

    if return_error_squares:
        array = np.column_stack((array,errors))

    return array
    

def getMiniGroupHists(lpgbt_hists, minigroups_swap, root=False, return_error_squares=False):
    
    minigroup_hists = []
    minigroup_hists_errors = []

    minigroup_hists_phiGreater60 = {}
    minigroup_hists_phiLess60 = {}

    #Get binning from lpgbt hists
    example_hist = lpgbt_hists[0][0]
    nROverZBins = example_hist.GetNbinsX()
    rOverZMin = example_hist.GetXaxis().GetBinLowEdge(1)
    rOverZMax = example_hist.GetXaxis().GetBinLowEdge(nROverZBins + 1)

    for minigroup, lpgbts in minigroups_swap.items():
        
        phiGreater60 = ROOT.TH1D( "minigroup_ROverZ_silicon_" + str(minigroup) + "_0","",nROverZBins,rOverZMin,rOverZMax)
        phiLess60    = ROOT.TH1D( "minigroup_ROverZ_silicon_" + str(minigroup) + "_1","",nROverZBins,rOverZMin,rOverZMax)

        for lpgbt in lpgbts:

            phiGreater60.Add( lpgbt_hists[0][lpgbt] )
            phiLess60.Add( lpgbt_hists[1][lpgbt] )
            
        phiGreater60_array = TH1D2array(phiGreater60,return_error_squares=return_error_squares)
        phiLess60_array = TH1D2array(phiLess60,return_error_squares=return_error_squares) 

        if ( root ):
            minigroup_hists_phiGreater60[minigroup] = phiGreater60
            minigroup_hists_phiLess60[minigroup] = phiLess60
        else:
            minigroup_hists_phiGreater60[minigroup] = phiGreater60_array
            minigroup_hists_phiLess60[minigroup] = phiLess60_array
            
    minigroup_hists.append(minigroup_hists_phiGreater60)
    minigroup_hists.append(minigroup_hists_phiLess60)

    return minigroup_hists
    
def getlpGBTHists(data, module_hists):

    lpgbt_hists = []

    for p,phiselection in enumerate(module_hists):#phi > 60 and phi < 60

        temp = {}

        #Get binning from module hists
        example_hist = next(iter(phiselection.values()))
        nROverZBins = example_hist.GetNbinsX()
        rOverZMin = example_hist.GetXaxis().GetBinLowEdge(1)
        rOverZMax = example_hist.GetXaxis().GetBinLowEdge(nROverZBins + 1)

        maxLpGBTNumber = data[['TPGId1','TPGId2']].max().max()
        
        for lpgbt in range(0,maxLpGBTNumber+1):
            lpgbt_found = False

            lpgbt_hist = ROOT.TH1D( ("lpgbt_ROverZ_silicon_" + str(lpgbt) + "_" + str(p)),"",nROverZBins,rOverZMin,rOverZMax)
            
            for tpg_index in ['TPGId1','TPGId2']:#lpgbt may be in the first or second position in the file

                for index, row in (data[data[tpg_index]==lpgbt]).iterrows():  
                    lpgbt_found = True
                    
                    if (row['density']==2):#Scintillator

                        hist = phiselection[ 1, row['u'], row['v'], row['layer'] ] # get module hist
                    else:
                        hist = phiselection[ 0, row['u'], row['v'], row['layer'] ] # get module hist        

                    linkfrac = 'TPGeLinkFrac1'
                    if ( tpg_index == 'TPGId2' ):
                        linkfrac = 'TPGeLinkFrac2'

                    lpgbt_hist.Add( hist, row[linkfrac] ) # add module hist with the correct e-link weight

            if lpgbt_found:
                temp[lpgbt] = lpgbt_hist

        lpgbt_hists.append(temp)
    
    return lpgbt_hists

def getMiniModuleGroups(data,minigroups_swap):

    minigroups_modules = {}
    
    for minigroup_id,lpgbts in minigroups_swap.items():

        module_list = []

        for lpgbt in lpgbts:

            data_list = data[ ((data['TPGId1']==lpgbt) | (data['TPGId2']==lpgbt)) ]

            for index, row in data_list.iterrows():

                if ( row['density']==2 ):
                    mod = [1, row['u'], row['v'], row['layer']]
                else:
                    mod = [0, row['u'], row['v'], row['layer']]

                if mod not in module_list:
                    module_list.append(mod)

        minigroups_modules[minigroup_id] = module_list
        
    return minigroups_modules
    
def getMiniTowerGroups(data, minigroups_modules):

    minigroups_towers = {}
    
    for mg,module_list in minigroups_modules.items():

        towerlist = []
        for module in module_list:
            if ( (module[0],module[3],module[1],module[2]) in data ):
                #(silicon CE-E/CE-H/scintillator), layer, u, v
                towerlist.append(data[module[0],module[3],module[1],module[2]]);
            else:
                print ( "Module missing from tower file (layer, u, v): ",module[3],module[1],module[2])

        #Remove duplicates
        towerlist = [item for sublist in towerlist for item in sublist]
        if len(towerlist) > 0:
            towerlist = np.unique(towerlist,axis=0).tolist()
        minigroups_towers[mg] = towerlist

    return minigroups_towers

def getTowerBundles(minigroups_towers, bundles, phisplit=None):

    #For each bundle get a list of the unique tower bins touched by the constituent modules
    #phisplit is a list indicating the bin numbers in phi where a split must be made
    all_bundles_towers = []
    
    for bundle in bundles:
        bundle_towers = []
        for minigroup in bundle:
            bundle_towers.append(minigroups_towers[minigroup])
        #Remove duplicates
        bundle_towers = [item for sublist in bundle_towers for item in sublist]
        bundle_towers = np.unique(bundle_towers,axis=0)

        bundle_towers_phi_split = []
        
        if phisplit != None and len(phisplit)>0:
            bundle_towers_phi_split.append(bundle_towers[(bundle_towers[:,2]<phisplit[0])].tolist())
            for i in range( (len(phisplit)-1) ):
                bundle_towers_phi_split.append(bundle_towers[(bundle_towers[:,2]>=phisplit[i])&(bundle_towers[:,2]<phisplit[i+1])].tolist())
            bundle_towers_phi_split.append(bundle_towers[(bundle_towers[:,2]>=phisplit[-1])].tolist())
        else:
            bundle_towers_phi_split.append(bundle_towers)
        
        all_bundles_towers.append(bundle_towers_phi_split)
        
    return all_bundles_towers

def getMinilpGBTGroups(data, minigroup_type="minimal"):

    minigroups = {}
    counter = 0
    minigroups_swap = {}
    n_lpgbt = 1 + np.max([np.max(data["TPGId1"]), np.max(data["TPGId2"])])
    
    if(minigroup_type=='minimal'):
        for index, row in data.iterrows():
            if (row['nTPG']==0): continue
            elif (row['nTPG']==1): 
                #If there is only one lpgbt attached to the module
                #Check if it is attached to another module
                #If it is not create a new minigroup
                #which is labelled with 'counter'
                if not (row['TPGId1'] in minigroups):

                    minigroups[row['TPGId1']] = counter
                    counter+=1

                
            elif (row['nTPG']==2): 
                
                #If there are two lpgbts attached to the module
                #The second one is "new"

                if (row['TPGId2'] in minigroups):
                    print ("should not be the case?")

                if (row['TPGId1'] in minigroups):
                    minigroups[row['TPGId2']] = minigroups[row['TPGId1']]
                else:
                    minigroups[row['TPGId1']] = counter
                    minigroups[row['TPGId2']] = counter

                    counter+=1

    elif (minigroup_type=='bylayer_silicon_seprated'):
        layer = 1 
        for lpgbt in range(n_lpgbt):
            layers_connected_to_lpgbt = list(data["layer"][(data["TPGId1"]==lpgbt) | (data["TPGId2"]==lpgbt)])
            if(len(dict.fromkeys(layers_connected_to_lpgbt))!=1):
                print("LpGBT connected to more than one layer. Change algorithm!!")
            
            if (layer !=layers_connected_to_lpgbt[0]):
                counter+=1
            minigroups[lpgbt] = counter
            layer = layers_connected_to_lpgbt[0]

    elif (minigroup_type=='bylayer'):
        layer = 1
        for lpgbt in range(n_lpgbt):
            layers_connected_to_lpgbt = list(data["layer"][(data["TPGId1"]==lpgbt) | (data["TPGId2"]==lpgbt)])
            if(len(dict.fromkeys(layers_connected_to_lpgbt))!=1):
                print("LpGBT connected to more than one layer. Change algorithm!!")

            if (layer !=layers_connected_to_lpgbt[0]):
                counter+=1
            if(layer ==50 and layers_connected_to_lpgbt[0]==37):
                counter=22
            minigroups[lpgbt] = counter
            layer = layers_connected_to_lpgbt[0]    

    for lpgbt, minigroup in minigroups.items(): 
        if minigroup in minigroups_swap: 
            minigroups_swap[minigroup].append(lpgbt) 
        else: 
            minigroups_swap[minigroup]=[lpgbt] 

    return minigroups,minigroups_swap

def find_nearest(array, values):
    indices = []
    for value in values:
        array_subtract = array - value
        index = (np.abs(array_subtract)).argmin()
        if array_subtract[index] > 0: index+=1
        indices.append(index)
        
    return indices

def getBundles(minigroups_swap,combination,nBundles=24,maxInputs=72):
    #Reimplemented in externals/mlrose_mod/opt_probs.py
    
    #Need to divide the minigroups into nBundles groups (24 by default) taking into account their different size

    #The weights are the numbers of lpgbts in each mini-groups
    try:
        weights = np.array([ len(minigroups_swap[x]) for x in combination ])
    except KeyError:
        print ("Requested minigroup does not exist in input mapping file")
        exit()
    cumulative_arr = weights.cumsum() / weights.sum()
    #Calculate the indices where to perform the split

    #Method 1
    #idx = np.searchsorted(cumulative_arr, np.linspace(0, 1, nBundles, endpoint=False)[1:])
    #Method 2 (improved)
    idx = find_nearest(cumulative_arr, np.linspace(0, 1, nBundles, endpoint=False)[1:])

    bundles = np.array_split(combination,idx)

    for bundle in bundles:
        weight_bundles = np.array([ len(minigroups_swap[x])  for x in bundle ])
        if (weight_bundles.sum() > maxInputs ):
            print ( "Error: more than " + str(maxInputs) + " lpgbts in bundle")
            
    return bundles
            
def getBundledlpgbtHistsRoot(minigroup_hists,bundles):

    bundled_lpgbthists = []
    bundled_lpgbthists_list = []

    #Get binning from minigroup hists
    example_hist = minigroup_hists[0][0]
    nROverZBins = example_hist.GetNbinsX()
    rOverZMin = example_hist.GetXaxis().GetBinLowEdge(1)
    rOverZMax = example_hist.GetXaxis().GetBinLowEdge(nROverZBins + 1)

    for p,phiselection in enumerate(minigroup_hists):

        temp = {}

        for i in range(len(bundles)):#loop over bundles

            #Create one lpgbt histogram per bundle
            lpgbt_hist = ROOT.TH1D( ("lpgbt_ROverZ_bundled_" + str(i) + "_" + str(p)),"",nROverZBins,rOverZMin,rOverZMax)
            
            for minigroup in bundles[i]:#loop over each lpgbt in the bundle
                lpgbt_hist.Add( phiselection[minigroup] )

            temp[i] = lpgbt_hist

        bundled_lpgbthists_list.append(temp)

    return bundled_lpgbthists_list


def getBundledlpgbtHists(minigroup_hists,bundles):

    bundled_lpgbthists_list = []

    for p,phiselection in enumerate(minigroup_hists):

        temp_list = {}

        for i in range(len(bundles)):#loop over bundles

            #Create one lpgbt histogram per bundle
            #Use minigroup_hists[0][0] shape as  the minigroup_hists
            #might have been produced with additional squared error information    
            lpgbt_hist_list = np.zeros(minigroup_hists[0][0].shape)

            for minigroup in bundles[i]:#loop over each lpgbt in the bundle
                lpgbt_hist_list+= phiselection[minigroup] 

            temp_list[i] = lpgbt_hist_list

        bundled_lpgbthists_list.append(temp_list)

    return bundled_lpgbthists_list

def getNumberOfModulesInEachBundle(minigroups_modules,bundles):

    data = []
    for bundle in bundles:
        size_of_bundle = 0
        for minigroup in bundle:
            size_of_bundle += len(minigroups_modules[minigroup])
        data.append(size_of_bundle)

    return data

def getMaximumNumberOfModulesInABundle(minigroups_modules,bundles):

    maximum = max(getNumberOfModulesInEachBundle( minigroups_modules,bundles ))
    return maximum

def calculateChiSquared(inclusive,grouped,nBundles=24,max_modules=None,weight_max_modules=1000,max_towers=None,weight_max_towers=[1000,1,180],weight_proportionally=True):

    use_error_squares = False
    #Check if the minigroup_hists were produced
    #with additional squared error information    

    if ( grouped[0][0].ndim == 2 ):
        use_error_squares = True

    #If optimisation of the number of modules in a bundle is performed
    #Aim for the maximum to be as low as possible -
    #i.e. for the number of modules in each bundle to be similar
    use_max_modules = False
    if ( max_modules != None):
        use_max_modules = True

    #If optimisation of the number of towers touched in a bundle is performed
    #Either aim for the maximum to be as low as possible (option 1) or
    #use a function that penalises very high values, but allows lower values with no penalty (option 2)
    use_max_towers = False
    if ( max_towers != None):
        use_max_towers = True

    chi2_total = 0
    
    for i in range(len(inclusive)):

        for key,hist in grouped[i].items():

            for b in range(len(hist)):

                if not use_error_squares:
                    squared_diff = np.power(hist[b]-inclusive[i].GetBinContent(b+1)/nBundles, 2 )
                    chi2_total+=squared_diff
                else:
                    weight_p = 1.0
                    if weight_proportionally:
                        if hist[b][0] > 0:
                            weight_p = hist[b][0]/500.

                    squared_diff = np.power(hist[b][0]-inclusive[i].GetBinContent(b+1)/nBundles, 2 ) / weight_p

                    squared_error = hist[b][1]

                    if ( squared_error == 0 ):
                        squared_error = inclusive[i].GetBinError(b+1)/nBundles
                    
                    if ( squared_error != 0 ):
                        chi2_total+=(squared_diff/squared_error)


    if use_max_modules:
        chi2_total += weight_max_modules * max_modules

    if use_max_towers:
        if weight_max_towers[1] == 1:
            chi2_total += weight_max_towers[0] * max_towers
        elif weight_max_towers[1] == 2:
            if max_towers > weight_max_towers[2]:
                chi2_total += weight_max_towers[0] * pow((max_towers-weight_max_towers[2]),2)
            
    return chi2_total
