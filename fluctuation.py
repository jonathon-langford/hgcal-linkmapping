#!/usr/bin/env python3
import random
random.seed(202008)
import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import math
import pickle
from scipy import optimize
from process import loadDataFile,loadConfiguration
from process import getPhiSplitIndices
from process import getMinilpGBTGroups,getBundles,getBundledlpgbtHists,getMiniModuleGroups
from rotate import rotate_to_sector_0
from geometryCorrections import applyGeometryCorrectionsNumpy,loadSiliconNTCCorrectionFile,applyGeometryCorrectionsTCPtRawData
from fluctuation_postprocess import plotMeanMax, plotTruncation, studyTruncationOptions, plot_Truncation_tc_Pt
import time
import yaml
import sys, os

def getMiniGroupHistsNumpy(module_hists, minigroups_modules):
    
    minigroup_hists = []

    minigroup_hists_phidivisionX = {}
    minigroup_hists_phidivisionY = {}

    #Get binning from module hists
    example_hist = next(iter(module_hists[0].values()))
    nROverZBins = len(example_hist)

    for minigroup, modules in minigroups_modules.items():
        
        phidivisionX = np.zeros(nROverZBins)
        phidivisionY = np.zeros(nROverZBins)
        
        for module in modules:
            phidivisionX = phidivisionX + module_hists[0][module[0],module[1],module[2],module[3]]
            phidivisionY = phidivisionY + module_hists[1][module[0],module[1],module[2],module[3]]
        
        minigroup_hists_phidivisionX[minigroup] = phidivisionX.copy()
        minigroup_hists_phidivisionY[minigroup] = phidivisionY.copy()

        
    minigroup_hists.append(minigroup_hists_phidivisionX)
    minigroup_hists.append(minigroup_hists_phidivisionY)

    return minigroup_hists


def getMiniGroupTCPtRawData(module_rawdata, minigroups_modules):
    
    minigroup_rawdata = []
    
    minigroup_rawdata_phidivisionX = {}
    minigroup_rawdata_phidivisionY = {}

    for minigroup, modules in minigroups_modules.items():

        phidivisionX = []
        phidivisionY = []

        for module in modules:
        
            phidivisionX += module_rawdata[0][module[0],module[1],module[2],module[3]] 
            phidivisionY += module_rawdata[1][module[0],module[1],module[2],module[3]]

        minigroup_rawdata_phidivisionX[minigroup] = phidivisionX.copy()
        minigroup_rawdata_phidivisionY[minigroup] = phidivisionY.copy()

    minigroup_rawdata.append(minigroup_rawdata_phidivisionX)
    minigroup_rawdata.append(minigroup_rawdata_phidivisionY)
    
    return minigroup_rawdata


def getBundledTCPtRawData(minigroup_rawdata,bundles):

    bundled_rawdata = []

    for phiselection in minigroup_rawdata:

        phi_region_rawdata = {}

        for i in range(len(bundles)):#loop over bundles
            
            one_bundle_rawdata = []

            for minigroup in bundles[i]:#loop over each minigroup in the bundle
                one_bundle_rawdata += phiselection[minigroup] 
            
            phi_region_rawdata[i] = one_bundle_rawdata.copy()

        bundled_rawdata.append(phi_region_rawdata)

    return bundled_rawdata


def getROverZPhi(x, y, z, sector = 0):

    if (z > 0):
        x = x*-1
    
    r = math.sqrt( x*x + y*y  )
    phi = np.arctan2(y,x)
    
    if (sector == 1):
        if ( phi < np.pi and phi > 0):
            phi = phi-(2*np.pi/3)
        else:
            phi = phi+(4*np.pi/3)
    elif (sector == 2):
        phi = phi+(2*np.pi/3)

    roverz_phi = [r/z,phi]
    return roverz_phi

def etaphiMapping(layer, etaphi, mappingFile):

    if (etaphi[1] > 24 and etaphi[1] <= 72):
        sector = 0
    elif (etaphi[1] > 72 and etaphi[1] <= 120):
        sector = 2
    else:
        sector = 1
        
    if (sector==0):
        pp=etaphi[1]-24
    elif (sector==2):
        pp=etaphi[1]-72
    elif (sector==1):
        if (etaphi[1]<=24):
            etaphi[1] = etaphi[1]+144
        pp = etaphi[1]-120
  
    pp = (pp-1)//4# //Phi index 1-12

    if "FeMappingV7" in mappingFile:
        if ( etaphi[0] <= 3 ):
            ep = 0
        elif ( etaphi[0] <= 9 ):
            ep = 1
        elif ( etaphi[0] <= 13 ):
            ep = 2
        elif ( etaphi[0] <= 17 ):
            ep = 3
        else:
            ep = 4
    elif "FeMappingTpgV7" in mappingFile:
        split = 12
        if layer > 40:
            split = 8
        if ( etaphi[0] <= split ):
            ep = 0
        else:
            ep = 1
    else:
        print( "Expected config file version to be either V7 or TpgV7" )
        
    return [ep,pp],sector

def applyTruncationAndGetPtSums(bundled_tc_Pt_rawdata, truncation_options, roverzBinning):

    #truncation_options is a list containing the truncation_options to study.
    #an element is inserted such that the sum without truncation is also available

    TCratio = []
    predetermined_values = []
    regionADefinitions = []
    regionBDefinitions = []
    definitionOptions = ['X', 'Y', 'X+Y']
    
    for option in truncation_options:
        if 'regionADefinition' in option.keys():
            if option['regionADefinition'] in definitionOptions:
                regionADefinitions.append(option['regionADefinition'])
            else:
                print ("Not a valid option for regionADefinition, assuming regionADefinition==X")
                regionADefinitions.append('X')
        else:
            print ("regionADefinition not given, assuming regionADefinition==X")
            regionADefinitions.append('X')

        if 'regionBDefinition' in option.keys():
            if option['regionBDefinition'] in definitionOptions:
                regionBDefinitions.append(option['regionBDefinition'])
            else:
                print ("Not a valid option for regionBDefinition, assuming regionBDefinition==Y")
                regionBDefinitions.append('Y')
        else:
            print ("regionBDefinition not given, assuming regionBDefinition==Y")
            regionBDefinitions.append('Y')
            
        TCratio.append(option['maxTCsA']/option['maxTCsB'])
        predetermined_values.append(option['predetermined_values'])
        
    truncation_max = np.full(len(predetermined_values[0]),1000)

    truncation_values_loop = predetermined_values.copy()
    TCratio_loop = TCratio.copy()
    regionADefinitions_loop = regionADefinitions.copy()
    regionBDefinitions_loop = regionBDefinitions.copy()
    truncation_values_loop.insert( 0, truncation_max )
    TCratio_loop.insert( 0, 1. )
    #The two phi divisions are saved independently, so both options can be recreated later
    regionADefinitions_loop.insert( 0, "X" )
    regionBDefinitions_loop.insert( 0, "Y" )
    
    alldata = [] #Output list

    for a,(truncation,ratio,adef,bdef) in enumerate(zip(truncation_values_loop,TCratio_loop,regionADefinitions_loop,regionBDefinitions_loop)):
        
        bundled_pt_hists_truncated = []

        #Need to consider phi regions at the same time, as the relevant "A" or "B" region might be inclusive in phi
        regionA_truncated_summed = np.zeros(len(roverzBinning)-1)
        regionB_truncated_summed = np.zeros(len(roverzBinning)-1)
        
        #Loop over each bundle
        for b in range(len(bundled_tc_Pt_rawdata[0])):

            #Get lists of (r/z, pt) pairs
            phidivisionX = np.asarray(bundled_tc_Pt_rawdata[0][b])
            phidivisionY = np.asarray(bundled_tc_Pt_rawdata[1][b])
            inclusive = np.asarray(bundled_tc_Pt_rawdata[0][b] + bundled_tc_Pt_rawdata[1][b])
                       
            #Find out how many TCs should be truncated
            #Bin the raw pT data
            if adef == "X":
                regionA = phidivisionX[:]
            elif adef == "Y":
                regionA = phidivisionY[:]
            elif  adef == "X+Y":
                regionA = inclusive[:]
            #Should not be any further possibilities
            
            if bdef == "X":
                regionB = phidivisionX[:]
            elif bdef == "Y":
                regionB = phidivisionY[:]
            elif  bdef == "X+Y":
                regionB = inclusive[:]
            #Should not be any further possibilities

            digitised_regionA_rawdata = np.digitize(regionA[:,0],roverzBinning)
            digitised_regionB_rawdata = np.digitize(regionB[:,0],roverzBinning)                
            
            sumPt_truncated_regionA = np.zeros(len(roverzBinning)-1)
            sumPt_truncated_regionB = np.zeros(len(roverzBinning)-1)

            #Loop over R/Z bins

            for roverz in range(len(roverzBinning)-1):
                #Get the pT values for the relevant R/Z bin
                pt_values_regionA = regionA[digitised_regionA_rawdata==roverz+1][:,1] #roverz+1 to convert from index to digitised bin number
                pt_values_regionB = regionB[digitised_regionB_rawdata==roverz+1][:,1] 
                #Get the number to be truncated in each region in an R/Z bin
                number_truncated_regionA = int(max(0,len(pt_values_regionA)-truncation[roverz]))
                if ( ratio.is_integer() ):
                    number_truncated_regionB = int(max(0,len(pt_values_regionB)-(np.ceil(truncation[roverz]/ratio))))#ceil rather than round in the cases ratio=1 or ratio=2 to make sure 0.5 goes to 1.0 (not default in python). 
                else:
                    number_truncated_regionB = int(max(0,len(pt_values_regionB)-(np.round(truncation[roverz]/ratio))))
                #Find the lowest 'n' values (number_truncated[roverz]), in order to truncate these

                sum_truncated_regionA = 0
                sum_truncated_regionB = 0
                if number_truncated_regionA > 0:
                    #Sort the pt_values array such that lower pT values are below number_truncated_regionA-1 and higher pT values are above
                    partition_regionA = pt_values_regionA[np.argpartition(pt_values_regionA, number_truncated_regionA-1)]#-1 to convert to index number 
                    #Then sum the lower values
                    sum_truncated_regionA = np.cumsum(partition_regionA)[number_truncated_regionA-1]
                if number_truncated_regionB > 0:
                    partition_regionB = pt_values_regionB[np.argpartition(pt_values_regionB, number_truncated_regionB-1)]#-1 to convert to index number 
                    sum_truncated_regionB = np.cumsum(partition_regionB)[number_truncated_regionB-1]

                #Save the sum pT after truncation
                total_sum_regionA = np.sum(pt_values_regionA)
                total_sum_regionB = np.sum(pt_values_regionB)

                sumPt_truncated_regionA[roverz] = total_sum_regionA - sum_truncated_regionA
                sumPt_truncated_regionB[roverz] = total_sum_regionB - sum_truncated_regionB

            #Sum all bundles together at this stage
            regionA_truncated_summed+=sumPt_truncated_regionA
            regionB_truncated_summed+=sumPt_truncated_regionB

        bundled_pt_hists_truncated.append( regionA_truncated_summed.copy() )
        bundled_pt_hists_truncated.append( regionB_truncated_summed.copy() )

        alldata.append(bundled_pt_hists_truncated.copy())
        
    return alldata


def checkFluctuations(initial_state, cmsswNtuple, outputName="alldata", tcPtConfig=None, truncationConfig = None, binningConfig = None, save_ntc_hists=False, beginEvent = -1, endEvent = -1):

    if ( binningConfig != None ):        
        nROverZBins = binningConfig["nROverZBins"]
        rOverZMin = binningConfig["rOverZMin"]
        rOverZMax = binningConfig["rOverZMax"]
    else:
        #Set defaults
        nROverZBins = 42
        rOverZMin = 0.07587128
        rOverZMax = 0.55563514

    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = nROverZBins, range = (rOverZMin,rOverZMax) )
    roverzBinning = inclusive_hists[1]

    #Load allocation information
    info = loadConfiguration(initial_state)
    data = info['data']
    init_state = info['mapping']
    nBundles = info['nBundles']
    maxInputs = info['maxInputs']
    mappingFile = info['configuration']['MappingFile']
    phisplitConfig = info['phisplitConfig']
    correctionConfig = info['correctionConfig']
    CMSSW_ModuleHists = info['CMSSW_ModuleHists']
    
    #Load the truncation options, if need to truncate based on E_T when running over ntuple (save_sum_tcPt == True)
    truncation_options = []
    save_sum_tcPt = False
    if ( tcPtConfig != None ):
        save_sum_tcPt = tcPtConfig['save_sum_tcPt']
        if save_sum_tcPt:
            options_to_study = tcPtConfig['options_to_study']
            if ( truncationConfig != None ):
                for option in options_to_study:
                    truncation_options.append(truncationConfig['option'+str(option)])

    #Load the CMSSW ntuple to get per event and per trigger cell information
    rootfile = ROOT.TFile.Open( cmsswNtuple , "READ" )
    tree = rootfile.Get("HGCalTriggerNtuple")

    max_ieta = 2 #Definition of a scintillator module is different between V7 and TpgV7 mapping files 
    if "FeMappingV7" in mappingFile:
        max_ieta = 5
    
    #Load geometry corrections
    if correctionConfig['nTCCorrectionFile'] != None:
        modulesToCorrect = loadSiliconNTCCorrectionFile( correctionConfig['nTCCorrectionFile'] )
    else:
        modulesToCorrect = pd.DataFrame()
        
    #Get list of which lpgbts are in each minigroup
    minigroups,minigroups_swap = getMinilpGBTGroups(data)

    #Get list of which modules are in each minigroup
    minigroups_modules = getMiniModuleGroups(data,minigroups_swap)
    bundles = getBundles(minigroups_swap,init_state,nBundles,maxInputs)

    bundled_lpgbthists_allevents = []
    bundled_pt_hists_allevents = []
    
    ROverZ_per_module_phidivisionX = {} #traditionally phi > 60 degrees
    ROverZ_per_module_phidivisionY = {} #traditionally phi < 60 degrees
    ROverZ_per_module_phidivisionX_tcPt = {} 
    ROverZ_per_module_phidivisionY_tcPt = {} 

    nTCs_per_module = {}

    #Value of split in phi (traditionally 60 degrees)
    if phisplitConfig == None:
        phi_split_phidivisionX = np.full( nROverZBins, np.pi/3 )
        phi_split_phidivisionY = np.full( nROverZBins, np.pi/3 )
    else:
        if phisplitConfig['type'] == "fixed":
            phi_split_phidivisionX = np.full( nROverZBins, np.radians(phisplitConfig['phidivisionX_fixvalue_min']) )
            phi_split_phidivisionY = np.full( nROverZBins, np.radians(phisplitConfig['phidivisionY_fixvalue_max']) )
        else:
            file_roverz_inclusive = ROOT.TFile(CMSSW_ModuleHists,"READ")
            PhiVsROverZ_Total = file_roverz_inclusive.Get("ROverZ_Inclusive")
            split_indices_phidivisionX = getPhiSplitIndices( PhiVsROverZ_Total, split = "per_roverz_bin")
            split_indices_phidivisionY = getPhiSplitIndices( PhiVsROverZ_Total, split = "per_roverz_bin")
            phi_split_phidivisionX = np.zeros( nROverZBins )
            phi_split_phidivisionY = np.zeros( nROverZBins )
            for i,(idxX,idxY) in enumerate(zip(split_indices_phidivisionX, split_indices_phidivisionY)):
                phi_split_phidivisionX[i] = PhiVsROverZ_Total.GetYaxis().GetBinLowEdge(int(idxY))
                phi_split_phidivisionY[i] = PhiVsROverZ_Total.GetYaxis().GetBinLowEdge(int(idxY))

    if save_ntc_hists:
        for i in range (15):
            for j in range (15):
                for k in range (1,53):
                    if  k < 28 and k%2 == 0:
                        continue
                    nTCs_per_module[0,i,j,k] = ROOT.TH1D( "nTCs_silicon_" + str(i) + "_" + str(j) + "_" + str(k), "", 49, -0.5, 48.5 )

        for i in range (max_ieta):
            for j in range (12):
                for k in range (37,53):
                    nTCs_per_module[1,i,j,k] = ROOT.TH1D( "nTCs_scintillator_" + str(i) + "_" + str(j) + "_" + str(k), "", 49, -0.5, 48.5 )


    for z in (-1,1):
        for sector in (0,1,2):
            key1 = (z,sector)
            ROverZ_per_module_phidivisionX[key1] = {}
            ROverZ_per_module_phidivisionY[key1] = {}
            if save_sum_tcPt:
                ROverZ_per_module_phidivisionX_tcPt[key1] = {}
                ROverZ_per_module_phidivisionY_tcPt[key1] = {}

            for i in range (15):
                for j in range (15):
                    for k in range (1,53):
                        if  k < 28 and k%2 == 0:
                            continue
                        ROverZ_per_module_phidivisionX[key1][0,i,j,k] = np.empty(0)
                        ROverZ_per_module_phidivisionY[key1][0,i,j,k] = np.empty(0)
                        if save_sum_tcPt:
                            ROverZ_per_module_phidivisionX_tcPt[key1][0,i,j,k] = [] #np.empty(0)
                            ROverZ_per_module_phidivisionY_tcPt[key1][0,i,j,k] = [] #np.empty(0)

            for i in range (max_ieta):
                for j in range (12):
                    for k in range (37,53):
                        ROverZ_per_module_phidivisionX[key1][1,i,j,k] = np.empty(0)
                        ROverZ_per_module_phidivisionY[key1][1,i,j,k] = np.empty(0)
                        if save_sum_tcPt:
                            ROverZ_per_module_phidivisionX_tcPt[key1][0,i,j,k] = [] #np.empty(0)
                            ROverZ_per_module_phidivisionY_tcPt[key1][0,i,j,k] = [] #np.empty(0)
    
    try:
        for entry,event in enumerate(tree):

            if ( beginEvent != -1 and entry < beginEvent ):
                if ( entry == 0 ):
                    print ("Event number less than " + str(beginEvent) + ", continue" )
                continue;

            if ( endEvent != -1 and entry > endEvent ):
                print ("Event number greater than " + str(endEvent) + ", break" )
                break;

            # if entry > 10:
            #     break
            print ("Event number " + str(entry))

            for key1 in ROverZ_per_module_phidivisionX.keys():
                for key2 in ROverZ_per_module_phidivisionX[key1].keys():
                    ROverZ_per_module_phidivisionX[key1][key2] = np.empty(0)
                    ROverZ_per_module_phidivisionY[key1][key2] = np.empty(0)
                    if save_sum_tcPt:
                        ROverZ_per_module_phidivisionX_tcPt[key1][key2] = [] #np.empty(0)
                        ROverZ_per_module_phidivisionY_tcPt[key1][key2] = [] #np.empty(0)
                        
            #Loop over list of trigger cells in a particular
            #event and fill R/Z histograms for each module
            #for phidivisionX and phidivisionY (traditionally phi > 60 degrees and phi < 60 degrees respectively)

            #Check if tc_pt exists (needed if we want to save the sum of (truncated) TC's pT)
            eventzip = zip(event.tc_waferu,event.tc_waferv,event.tc_layer,event.tc_x,event.tc_y,event.tc_z,event.tc_cellu,event.tc_cellv)
            if ( save_sum_tcPt ):
                if hasattr(event, 'tc_pt'):
                    eventzip = zip(event.tc_waferu,event.tc_waferv,event.tc_layer,event.tc_x,event.tc_y,event.tc_z,event.tc_cellu,event.tc_cellv,event.tc_pt)
                else:
                    print ('tc_pt not found in TTree - switching to non-save_sum_pt mode')
                    save_sum_tcPt = False

            for variables in eventzip:
                u,v,layer,x,y,z,cellu,cellv = variables[:8]
                if save_sum_tcPt: pt = variables[8]
                
                if ( u > -990 ): #Silicon
                    uv,sector = rotate_to_sector_0(u,v,layer)
                    roverz_phi = getROverZPhi(x,y,z,sector)
                    roverz_bin = np.argmax( roverzBinning > abs(roverz_phi[0]) )

                    if (roverz_phi[1] >= phi_split_phidivisionX[roverz_bin-1]):
                        #There should be no r/z values lower than rOverZMin (around 0.076)
                        ROverZ_per_module_phidivisionX[np.sign(z),sector][0,uv[0],uv[1],layer] = np.append(ROverZ_per_module_phidivisionX[np.sign(z),sector][0,uv[0],uv[1],layer],abs(roverz_phi[0]))
                        if save_sum_tcPt:
                            ROverZ_per_module_phidivisionX_tcPt[np.sign(z),sector][0,uv[0],uv[1],layer].append( [abs(roverz_phi[0]),pt] )
                    if (roverz_phi[1] < phi_split_phidivisionY[roverz_bin-1]):
                        ROverZ_per_module_phidivisionY[np.sign(z),sector][0,uv[0],uv[1],layer] = np.append(ROverZ_per_module_phidivisionY[np.sign(z),sector][0,uv[0],uv[1],layer],abs(roverz_phi[0]))
                        if save_sum_tcPt:
                            ROverZ_per_module_phidivisionY_tcPt[np.sign(z),sector][0,uv[0],uv[1],layer].append( [abs(roverz_phi[0]),pt] )
                        
                else: #Scintillator  
                    eta = cellu
                    phi = cellv
                    etaphi,sector = etaphiMapping(layer,[eta,phi],mappingFile)
                    roverz_phi = getROverZPhi(x,y,z,sector)
                    roverz_bin = np.argmax( roverzBinning > abs(roverz_phi[0]) )
                    
                    if (roverz_phi[1] >= phi_split_phidivisionX[roverz_bin-1]):
                        ROverZ_per_module_phidivisionX[np.sign(z),sector][1,etaphi[0],etaphi[1],layer] = np.append(ROverZ_per_module_phidivisionX[np.sign(z),sector][1,etaphi[0],etaphi[1],layer],abs(roverz_phi[0]))
                        if save_sum_tcPt:
                            ROverZ_per_module_phidivisionX_tcPt[np.sign(z),sector][1,etaphi[0],etaphi[1],layer].append( [abs(roverz_phi[0]),pt] )
                    if (roverz_phi[1] < phi_split_phidivisionY[roverz_bin-1]):
                        ROverZ_per_module_phidivisionY[np.sign(z),sector][1,etaphi[0],etaphi[1],layer] = np.append(ROverZ_per_module_phidivisionY[np.sign(z),sector][1,etaphi[0],etaphi[1],layer],abs(roverz_phi[0]))
                        if save_sum_tcPt:
                            ROverZ_per_module_phidivisionY_tcPt[np.sign(z),sector][1,etaphi[0],etaphi[1],layer].append( [abs(roverz_phi[0]),pt] )
            #Bin the TC module data
            module_hists_phidivisionX = {}
            module_hists_phidivisionY = {}

            for key1,value1 in ROverZ_per_module_phidivisionX.items():
                module_hists_phidivisionX[key1] = {}
                for key2,value2 in value1.items():
                    module_hists_phidivisionX[key1][key2] = np.histogram( value2, bins = nROverZBins, range = (rOverZMin,rOverZMax) )[0]

            for key1,value1 in ROverZ_per_module_phidivisionY.items():
                module_hists_phidivisionY[key1] = {}
                for key2,value2 in value1.items():
                    module_hists_phidivisionY[key1][key2] = np.histogram( value2, bins = nROverZBins, range = (rOverZMin,rOverZMax) )[0]

            for z in (-1,1):
                for sector in (0,1,2):
                        
                    #the module hists are a numpy array of size nROverZBins (42 by default)
                    module_hists = [module_hists_phidivisionX[z,sector],module_hists_phidivisionY[z,sector]]
                    
                    #Apply geometry corrections
                    applyGeometryCorrectionsNumpy( module_hists, modulesToCorrect )

                    #Save the integral of module_hists, per event 
                    if save_ntc_hists:
                        for module,hist in nTCs_per_module.items():
                            hist.Fill( np.round(np.sum(module_hists[0][module]) + np.sum(module_hists[1][module])) )
                    
                    #Sum the individual module histograms to get the minigroup histograms
                    minigroup_hists = getMiniGroupHistsNumpy(module_hists,minigroups_modules)

                    #Sum the minigroup histograms to get the bundle histograms
                    bundled_lpgbthists = getBundledlpgbtHists(minigroup_hists,bundles)

                    bundled_lpgbthists_allevents.append(bundled_lpgbthists)
                    
                    #Collect the individual TC Pt values for a given minigroup, with the view to truncate and sum
                    if ( save_sum_tcPt ):
                        tc_Pt_rawdata = [ROverZ_per_module_phidivisionX_tcPt[z,sector],ROverZ_per_module_phidivisionY_tcPt[z,sector]]

                        #Apply geometry corrections
                        applyGeometryCorrectionsTCPtRawData( tc_Pt_rawdata, modulesToCorrect )

                        #Get lists of (r/z, pt) pairs, first for minigroups and then for bundles
                        minigroup_tc_Pt_rawdata = getMiniGroupTCPtRawData(tc_Pt_rawdata,minigroups_modules)
                        bundled_tc_Pt_rawdata = getBundledTCPtRawData(minigroup_tc_Pt_rawdata,bundles)

                        #Get histograms of (truncated) sum pT per r/z bin
                        bundled_pt_hists = applyTruncationAndGetPtSums(bundled_tc_Pt_rawdata, truncation_options, roverzBinning)

                        bundled_pt_hists_allevents.append(bundled_pt_hists)

    except KeyboardInterrupt:
        print("interrupt received, stopping and saving")

    finally:

        #Write all data to file for later analysis (Pickling)
        if ( beginEvent != -1 ):
            outputName = outputName + "_from" + str(beginEvent)
        if ( endEvent != -1 ):
            outputName = outputName + "_to" + str(endEvent)
        
        with open( outputName + ".txt", "wb") as filep:
            pickle.dump(bundled_lpgbthists_allevents, filep)

        if save_sum_tcPt:
            with open( outputName + "_sumpt.txt", "wb") as filep:
                pickle.dump(bundled_pt_hists_allevents, filep)

        if save_ntc_hists:
            outfile = ROOT.TFile(outputName + "_nTCs.root","RECREATE")
            for hist in nTCs_per_module.values():
                hist.Write()

def main():

    try:
        config_file = sys.argv[1]
    except IndexError:
         print ("Please give valid config file")
         exit()
    try:
        with open(config_file,'r') as file:
            config = yaml.load(file,Loader=yaml.FullLoader)
    except EnvironmentError:
        print ("Please give valid config file")
        exit()
        
    #Code to process the input root file,
    #and to get the bundle R/Z histograms per event
    truncationConfig = None
    if 'truncationConfig' in config.keys():
        truncationConfig = config['truncationConfig']

    binningConfig = None
    if 'binningConfig' in config.keys():
        binningConfig = config['binningConfig']

    if (config['function']['checkFluctuations']):
        subconfig = config['checkFluctuations']
        
        tcPtConfig = None
        if 'tcPtConfig' in subconfig.keys():
            tcPtConfig = subconfig['tcPtConfig']

        checkFluctuations(initial_state=subconfig['initial_state'], cmsswNtuple=subconfig['cmsswNtuple'], outputName=subconfig['outputName'], tcPtConfig = tcPtConfig, truncationConfig = truncationConfig, binningConfig = binningConfig, save_ntc_hists=subconfig['save_ntc_hists'], beginEvent = subconfig['beginEvent'], endEvent = subconfig['endEvent'])

    #Plotting functions
    
    if (config['function']['plot_MeanMax']):
        subconfig = config['plot_MeanMax']
        plotMeanMax(eventData = config['eventData'], outdir = config['output_dir'], xyTreatment = subconfig['xyTreatment'], binningConfig = binningConfig, plotIndividualEvents = subconfig['plotIndividualEvents'])

    if (config['function']['plot_Truncation']):
        subconfig = config['plot_Truncation']
        plotTruncation(eventData = config['eventData'],outdir = config['output_dir'], useMaximumXY = subconfig['useMaximumXY'], binningConfig = binningConfig )
        
    if (config['function']['studyTruncationOptions']):
        subconfig = config['studyTruncationOptions']
        studyTruncationOptions(eventData = config['eventData'], options_to_study = subconfig['options_to_study'], truncation_values_method = subconfig['truncation_values_method'], truncationConfig = config['truncationConfig'], binningConfig = binningConfig, outdir = config['output_dir'] )
        
    if (config['function']['plot_Truncation_tc_Pt']):
        subconfig = config['plot_Truncation_tc_Pt']
        plot_Truncation_tc_Pt(eventData = subconfig['eventData_Pt'], options_to_study = subconfig['options_to_study'], truncationConfig = config['truncationConfig'], binningConfig = binningConfig, outdir = config['output_dir'] )

    
main()
