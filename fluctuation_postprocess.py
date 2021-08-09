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
import time
import yaml
import sys, os

def plotMeanMax(eventData, outdir = ".", xyTreatment = "maximum", binningConfig = None, plotIndividualEvents = False):

    #Load pickled per-event bundle histograms
    phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents = loadFluctuationData(eventData)

    #Plotting Max, mean and standard deviation per bundle:

    #For each r/z bin take either the sum of the phidivisionX and phidivisionY arrays (inclusive) or the maximum of the two
    if ( xyTreatment == 'maximum' ):
        bundled_lpgbthists_allevents = np.maximum(phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents)
        plotMeanMax_MakePlots(bundled_lpgbthists_allevents, binningConfig = binningConfig,
                              savename = "_xy_maximum", outdir = outdir, plotIndividualEvents = plotIndividualEvents)
    elif ( xyTreatment == 'inclusive' ):
        bundled_lpgbthists_allevents = phidivisionX_bundled_lpgbthists_allevents + phidivisionY_bundled_lpgbthists_allevents
        plotMeanMax_MakePlots(bundled_lpgbthists_allevents, binningConfig = binningConfig,
                              savename = "_xy_inclusive", outdir = outdir, plotIndividualEvents = plotIndividualEvents)
    elif ( xyTreatment == 'separate' ):
        plotMeanMax_MakePlots(phidivisionX_bundled_lpgbthists_allevents, binningConfig = binningConfig,
                              savename = "_phiDivisionX", outdir = outdir, plotIndividualEvents = plotIndividualEvents)
        plotMeanMax_MakePlots(phidivisionY_bundled_lpgbthists_allevents, binningConfig = binningConfig,
                              savename = "_phiDivisionY", outdir = outdir, plotIndividualEvents = plotIndividualEvents)

def plotMeanMax_MakePlots(bundled_lpgbthists_allevents, binningConfig = None, outdir = ".", savename = "", plotIndividualEvents = False):

    if ( binningConfig != None ):        
        rOverZMin = binningConfig["rOverZMin"]
        rOverZMax = binningConfig["rOverZMax"]
    else:
        #Set defaults
        rOverZMin = 0.07587128
        rOverZMax = 0.55563514

    nROverZBins = len(bundled_lpgbthists_allevents[0][0]) #42 by default
    nBundles = len(bundled_lpgbthists_allevents[0]) #24 by default
    
    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = nROverZBins, range = (rOverZMin,rOverZMax) )

    os.system("mkdir -p " + outdir)
    
    hists_max = [] 
    
    for bundle in range(nBundles):
        bundle_list = bundled_lpgbthists_allevents[:,bundle,:]

        hist_max = np.amax(bundle_list,axis=0)
        hist_mean = np.mean(bundle_list, axis=0)
        hist_std = np.std(bundle_list, axis=0)

        for s,std in enumerate(hist_std):
            hist_std[s] = std + hist_mean[s]

        pl.bar((inclusive_hists[1])[:-1], hist_max, width=0.012,align='edge')
        pl.bar((inclusive_hists[1])[:-1], hist_std, width=0.012,align='edge')
        pl.bar((inclusive_hists[1])[:-1], hist_mean, width=0.012,align='edge')

        #Plot all events for a given bundle on the same plot
        #(just for one bundle by default)
        if plotIndividualEvents and bundle == 0:
            for e,event in enumerate(bundle_list):
                pl.bar((inclusive_hists[1])[:-1], event, width=0.012,fill=False)
                
        pl.ylim((0,71))
        pl.savefig( outdir + "/bundle_" + str(bundle) + "max" + savename + ".png" )        
        pl.clf()

        hists_max.append(hist_max)
        
    #Plot maxima for all bundles on the same plot
    for hist in hists_max:
        pl.bar((inclusive_hists[1])[:-1], hist, width=0.012,align='edge')
    pl.ylim((0,71))
    pl.xlabel('r/z')
    pl.ylabel('Maximum number of TCs per bin')
    pl.savefig( outdir + "/maxima" + savename + ".png" )
    pl.clf()

def sumMaximumOverAllEventsAndBundles(truncation,data):
    #Solve for the truncation factor, given two datasets, A and B (data[0] and data[1])
    #And the maximum number of trigger cells allowed in each dataset (data[2] and data[3] respectively)
    #The truncation parameter must be less than or equal to 1
    maximum_A = np.amax(data[0],axis=(1,0))
    maximum_B = np.amax(data[1],axis=(1,0))

    Bscaling_factor = data[2]/data[3]
    
    maxAB = np.maximum(maximum_A,maximum_B*Bscaling_factor)

    #Use ceiling rather than round to get worst case
    overallsum_A = np.sum(np.amax(np.where(data[0]<truncation*maxAB,data[0],np.where(truncation*maxAB<maximum_A,truncation*maxAB,maximum_A)),axis=(1,0)))
    overallsum_B = np.sum(np.amax(np.where(data[1]<truncation*(maxAB/Bscaling_factor),data[1],np.where(truncation*(maxAB/Bscaling_factor)<maximum_B,truncation*maxAB/Bscaling_factor,maximum_B)),axis=(1,0)))

    valA = data[2] - overallsum_A
    valB = data[3] - overallsum_B

    #Give a preference that the sum is less than the maximum allowed
    if ( valA < 0 ):
         valA = valA * -1.5
    if ( valB < 0 ):
         valB = valB * -1.5

        
    optval = valA + valB

    return optval

def plot_NTCs_Vs_ROverZ(inputdata,axis,savename,truncation_curves=None,scaling=None):

    #Fill a 2D histogram per bunch-crossing with N_TCs (maximum over bundles) 

    #Each row represents the r/z bins in a bundle, there are n_bundles*n_events rows
    data = inputdata.reshape(-1,inputdata.shape[-1])

    #Swap axes, such that each row represents an r/z bin, there are n_roverz_bins rows (later flattened)
    data_swap = np.swapaxes(data,0,1)

    #Get the r/z bin axis indices, n_bundles*n_events*[0]+n_bundles*n_events*[1]+...n_bundles*n_events*[n_roverz_bins]
    axis_indices = np.where(data_swap==data_swap)[0]
    #Then get the roverz bin values corresponding to the indices
    roverz = np.take(axis,axis_indices)

    #Plot the 2D histogram
    nTCbins = int(50*(42/(len(axis)-1)))
    pl.clf()
    pl.hist2d( roverz , data_swap.flatten() , bins = (len(axis)-1,nTCbins),range=[[axis[0],axis[-1]], [0, nTCbins]],norm=LogNorm())
    pl.colorbar().set_label("Number of Events")
    #Plot the various 1D truncation curves
    colours = ['red','orange','cyan','green','teal','darkviolet']

    if ( truncation_curves is not None ):
        for t,truncation_option in enumerate(truncation_curves):
            scale = 1.
            if (scaling is not None):
                scale=scaling[t]
            plotted_truncation_curve = np.append(truncation_option,truncation_option[-1])/scale
            pl.step(axis,plotted_truncation_curve+1, where = 'post' , color=colours[t],linewidth='3')
            #plotted_truncation_curve+1 so that bin 'n' is visually included if the truncation value is 'n'
            #Note because of the geometric corrections the number of trigger cells might be fractional,
            #in which case the +1 is not correct (but only applies to the bins at low and high r/z)
            
    pl.xlabel('r/z')
    pl.ylabel('Number of TCs')
    pl.savefig( savename + ".png" )
    pl.clf()

def plot_frac_Vs_ROverZ( dataA, dataB, truncation_curve, TCratio, axis, title, savename ):

    #Sum over all events and bundles of TCs (in each R/Z bin) 
    totalsumA = np.sum( dataA , axis=(0,1) )
    totalsumB = np.sum( dataB , axis=(0,1) )

    #Sum over all events and bundles of truncated TCs (in each R/Z bin) 
    truncatedsum_A = np.sum(np.where(dataA<truncation_curve, dataA, truncation_curve),axis=(0,1))

    if ( TCratio.is_integer() ):
        truncation_curveB = np.ceil(truncation_curve/TCratio) #ceil rather than round in the cases ratio=1 or ratio=2 to make sure 0.5 goes to 1.0 (not default in python). 
    else:
        truncation_curveB = np.round(truncation_curve/TCratio)
    truncatedsum_B = np.sum(np.where(dataB<truncation_curveB, dataB, truncation_curveB),axis=(0,1))

    #Divide to get the fraction, taking into account division by zero
    ratioA = np.divide(   truncatedsum_A, totalsumA , out=np.ones_like(truncatedsum_A), where=totalsumA!=0 )
    ratioB = np.divide(   truncatedsum_B, totalsumB , out=np.ones_like(truncatedsum_B), where=totalsumB!=0 )

    pl.clf()
    pl.step(axis,np.append( ratioA , ratioA[-1] ),color='red',linewidth='1', where = 'post', label='data A')
    pl.step(axis,np.append( ratioB , ratioB[-1] ),color='orange',linewidth='1', where = 'post', label='data B')
    pl.xlabel('r/z')
    pl.ylabel('Sum TCs with truncation / Sum all TCs')
    pl.title(title)
    pl.ylim((0.6,1.05))
    pl.legend(loc='lower left')
    pl.savefig( savename + ".png" )
    pl.clf()    

def getTruncationValuesForGivenScalar(scalar,data):

    #For the get reverse truncation values method
    #For a given "Scalar" (i.e. efficiency or sum of truncated TCs / sum of all TCs)

    #Calculate the B scaling factor such that it has the same weight in the fit as A
    Bscaling_factor = data[2]/data[3]

    #Sum over all events and bundles of TCs (in each R/Z bin) 
    totalsumA = np.sum( data[0] , axis=1 )
    totalsumB = np.sum( data[1] , axis=1 )
    totalsum = (totalsumA + Bscaling_factor * totalsumB) / 2
        
    #Largest possible truncation value (end of loop)
    largestTruncationMax = 100

    #Loop over bins
    truncation_values = []

    for b in range(len(data[0])):
        scalarset = False
        for x in range (0,largestTruncationMax,1):
            
            truncation_sumA = np.sum(np.where(data[0][b]<x, data[0][b], x))
            truncation_sumB = np.sum(np.where(data[1][b]*Bscaling_factor<x, data[1][b]*Bscaling_factor, x))
            truncation_sum = (truncation_sumA + truncation_sumB) / 2
            
            if ( scalar*totalsum[b] == 0):
                truncation_values.append(0)
                scalarset = True
                break
            
            if truncation_sum > scalar*totalsum[b]:#break when the truncated sum is below the desired value
                
                #There is the option to get the truncation sum for x-1 and check whether this or the sum for x is closer to giving the desired efficiency.
                #Or by default we assume we always want greater than the desired efficiency
                use_higher_efficiency=True

                if use_higher_efficiency:
                    truncation_values.append(x)
                else:
                    truncation_sum_x_minus_1A = np.sum(np.where(data[0][b]<x-1, data[0][b], x-1))
                    truncation_sum_x_minus_1B = np.sum(np.where(data[1][b]*Bscaling_factor<x-1, data[1][b]*Bscaling_factor, x-1))
                    truncation_sum_x_minus_1 = (truncation_sum_x_minus_1A + truncation_sum_x_minus_1B) / 2
                    if ( abs( scalar*totalsum[b] - truncation_sum ) < abs( scalar*totalsum[b] - truncation_sum_x_minus_1 )):
                        truncation_values.append(x)
                    else:
                        truncation_values.append(x-1)
                        
                scalarset = True
                break

        if not scalarset:
            print ("Scalar not set for bin", b)
            print ("This should not generally happen")
            truncation_values.append(1)
    
    return truncation_values
    
def getOptimalScalingValue(scalar,data):

    #Calculate the B scaling factor such that it has the same weight in the fit as A
    Bscaling_factor = data[2]/data[3]

    truncation_values = getTruncationValuesForGivenScalar(scalar,data)
    difference = abs(data[2]-np.sum(truncation_values)) + Bscaling_factor*abs(data[3]-np.sum(truncation_values))
        
    return difference
    
def getReverseTruncationValues( dataA, dataB, maxTCsA, maxTCsB ):
    #Obtain the truncation values by demanding that the efficiency be constant as a function of r/z
    
    nbins = len(dataA[0][0]) #42 by default

    #Option to read out means
    print_means = True
    mean_A = np.mean(np.sum(dataA, axis=(2)))
    mean_B = np.mean(np.sum(dataB, axis=(2)))
    mean_A2 = np.mean(dataA, axis=(0,1))
    mean_B2 = np.mean(dataB, axis=(0,1))
    if print_means:
        print ("mean region A = "  + str(mean_A))
        print ("mean region B = "  + str(mean_B))
    
    #reorganise data to get bin b values across all events and bundles
    reorganisedDataA = []
    reorganisedDataB = []
    for b in range(nbins):
        reorganisedDataA.append(np.array([j[b] for i in dataA for j in i]))
        reorganisedDataB.append(np.array([j[b] for i in dataB for j in i]))
    reorganisedDataA = np.asarray(reorganisedDataA)
    reorganisedDataB = np.asarray(reorganisedDataB)
            
    #Get optimal scalar (efficiency), i.e. the highest possible, whilst also keeping constant as a function of r/z
    result = optimize.minimize_scalar(getOptimalScalingValue,args=[reorganisedDataA,reorganisedDataB,maxTCsA,maxTCsB],bounds=(0.8,1.0),method='bounded')

    #Get first iteration of truncation values
    truncation_values = getTruncationValuesForGivenScalar(result.x,[reorganisedDataA,reorganisedDataB,maxTCsA,maxTCsB])

    #Give 'spare' TCs to low bins (or take away from higher bins)
    truncation_values = np.array(truncation_values)

    print ("Sum of truncation values before redistribution = " + str(np.sum(truncation_values)) + "; max = " + str(maxTCsA))
    print (truncation_values)

    n_spare = maxTCsA - np.sum(truncation_values)
    if ( n_spare != 0):
        give_to_zero = False #Redistribute to bins where the truncation value is 0 (normally the last two in r/z)
        if give_to_zero:
            argordered = np.argsort(truncation_values)
        else:                   
            non_zero_indices = np.nonzero(truncation_values)[0]
            arg_ordered = non_zero_indices[np.argsort(truncation_values[non_zero_indices])]

        arg = arg_ordered[::np.sign(n_spare)][:int(abs(n_spare))]
        truncation_values[arg]+=np.sign(n_spare)
    print ("Sum of truncation values after redistribution = " + str(np.sum(truncation_values)) + "; max = " + str(maxTCsA))
    
    return truncation_values
    
def getTruncationValuesRoverZ(data_A, data_B, maxtcs_A, maxtcs_B):
    #Get an array of size nROverZbins, which indicates the maximum number of TCs allowed in each RoverZ bin
    
    #'scalar' is the value by which the maximum (of data_A or data_B x TCratio) is multiplied to get the truncation values
    # result = optimize.minimize_scalar(sumMaximumOverAllEventsAndBundles,args=[data_A, data_B, maxtcs_A, maxtcs_B],bounds=(-1,1.0),method='bounded')
    # scalar2 = result.x

    maximum_A = np.amax(data_A,axis=(1,0))
    maximum_B = np.amax(data_B,axis=(1,0))

    #Option to read out means
    print_means = False
    mean_A = np.mean(np.sum(data_A, axis=(2)))
    mean_B = np.mean(np.sum(data_B, axis=(2)))
    mean_A2 = np.mean(data_A, axis=(0,1))
    mean_B2 = np.mean(data_B, axis=(0,1))
    if print_means:
        print ("mean region A = "  + str(mean_A))
        print ("mean region B = "  + str(mean_B))

    pc_A = np.percentile(data_A, 99, axis=(0,1))
    pc_B = np.percentile(data_B, 99,axis=(0,1))

    Bscaling_factor = maxtcs_A/maxtcs_B
    #maxAB = np.maximum(maximum_A,maximum_B*Bscaling_factor)
    maxABpc = np.maximum( pc_A, pc_B * Bscaling_factor )#percentile

    #Find the floating-point truncation values and round down to make these integers.
    #This will be less than the allowed total due to rounding down.

    scalar = min( maxtcs_A/np.sum(maxABpc), maxtcs_B*Bscaling_factor/np.sum(maxABpc ))
    truncation_float = maxABpc * scalar
    truncation_floor = np.floor(truncation_float)

    #The integer difference from the allowed total gives the number of bins that should have their limit incremented by 1.
    integer_difference = np.floor(np.sum(truncation_float)-np.sum(truncation_floor))
    
    #Find the N bins which have the smallest floor values
    #and add 1 to these. This gives limits for A, and for B (divided by TC ratio)

    arg = np.argsort(truncation_floor)[:int(integer_difference)]
    truncation_floor[arg]+=1

    #Reassign from highest to lowest
    nTimesToReassign = 1
    nReassign = 10
    for n in range (nTimesToReassign):
        argReduce = np.argsort(truncation_floor)[-int(nReassign):]
        argIncrease = np.argsort(truncation_floor)[:int(nReassign)]
        truncation_floor[argIncrease]+=1
        truncation_floor[argReduce]-=1

    #Reduce the maximum bins of truncation_floor if sum of truncation_floor is larger than that allowed by maxtcs_A and maxtcs_B
    #Done consecutively in A and B so as not to overcorrect
    diffA = np.sum(truncation_floor) - maxtcs_A
    if ( diffA < 0 ): diffA = 0
    arg = np.argsort(truncation_floor)[:int(diffA)]
    #arg is a list of the indices of the highest bins
    truncation_floor[arg]-=1

    diffB = np.sum(truncation_floor)*(maxtcs_B/maxtcs_A) - maxtcs_B
    if ( diffB < 0 ): diffB = 0
    arg = np.argsort(truncation_floor)[:int(diffB)]
    truncation_floor[arg]-=1
    
    return truncation_floor

def loadFluctuationData(eventData):
    #Load the per-event flucation data produced using 'checkFluctuations'
    #Return two arrays (for phi divisions X and Y) containing for each event and
    #bundle, the number of TCs in each R/Z bin]
    
    with open(eventData, "rb") as filep:   
        bundled_lpgbthists_allevents = pickle.load(filep)
    
    #Names for phi > split_value and phi < split_value indices
    dataX = 0
    dataY = 1

    nBundles = len(bundled_lpgbthists_allevents[0][0]) #24 by default
    nbins = len(bundled_lpgbthists_allevents[0][0][0]) #42 by default
    
    dataX_bundled_lpgbthists_allevents = np.empty((len(bundled_lpgbthists_allevents),nBundles,nbins))
    dataY_bundled_lpgbthists_allevents = np.empty((len(bundled_lpgbthists_allevents),nBundles,nbins))

    for e,event in enumerate(bundled_lpgbthists_allevents):        
        dataX_bundled_lpgbthists_allevents[e] = np.array(list(event[dataX].values()))
        dataY_bundled_lpgbthists_allevents[e] = np.array(list(event[dataY].values()))

    return dataX_bundled_lpgbthists_allevents, dataY_bundled_lpgbthists_allevents
    
def studyTruncationOptions(eventData, options_to_study, truncation_values_method, truncationConfig, binningConfig = None, outdir = "."):

    if ( binningConfig != None ):        
        rOverZMin = binningConfig["rOverZMin"]
        rOverZMax = binningConfig["rOverZMax"]
    else:
        #Set defaults
        rOverZMin = 0.07587128
        rOverZMax = 0.55563514

    #Load pickled per-event bundle histograms
    phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents = loadFluctuationData(eventData)

    os.system("mkdir -p " + outdir)

    nROverZBins = len(phidivisionX_bundled_lpgbthists_allevents[0][0]) #42 by default
    
    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = nROverZBins, range = (rOverZMin, rOverZMax) )

    inclusive_bundled_lpgbthists_allevents = phidivisionX_bundled_lpgbthists_allevents + phidivisionY_bundled_lpgbthists_allevents

    #RegionA and Region B are defined in terms of phidivisionX and phidivisionY in the configuration (either X, Y or the sum of the two)
    
    truncation_values = []
    truncation_options = []
    regionA_bundled_lpgbthists_allevents = []
    regionB_bundled_lpgbthists_allevents = []

    for option in options_to_study:
        print ("Get truncation value for option " + str(option))

        truncation_options.append(truncationConfig['option'+str(option)])

        if 'regionADefinition' in truncation_options[-1].keys():
            if truncation_options[-1]['regionADefinition'] == "X":
                regionA_bundled_lpgbthists_allevents.append(phidivisionX_bundled_lpgbthists_allevents)
            elif truncation_options[-1]['regionADefinition'] == "Y":
                regionA_bundled_lpgbthists_allevents.append(phidivisionY_bundled_lpgbthists_allevents)
            elif truncation_options[-1]['regionADefinition'] == "X+Y":
                regionA_bundled_lpgbthists_allevents.append(inclusive_bundled_lpgbthists_allevents)
            else:
                print ("Not a valid option for regionADefinition, assuming regionADefinition==X")
                regionA_bundled_lpgbthists_allevents.append(phidivisionX_bundled_lpgbthists_allevents)
        else:
            print ("regionADefinition not given, assuming regionADefinition==X")
            regionA_bundled_lpgbthists_allevents.append(phidivisionX_bundled_lpgbthists_allevents)

                
        if 'regionBDefinition' in truncation_options[-1].keys():
            if truncation_options[-1]['regionBDefinition'] == "X":
                regionB_bundled_lpgbthists_allevents.append(phidivisionX_bundled_lpgbthists_allevents)
            elif truncation_options[-1]['regionBDefinition'] == "Y":
                regionB_bundled_lpgbthists_allevents.append(phidivisionY_bundled_lpgbthists_allevents)
            elif truncation_options[-1]['regionBDefinition'] == "X+Y":
                regionB_bundled_lpgbthists_allevents.append(inclusive_bundled_lpgbthists_allevents)
            else:
                print ("Not a valid option for regionBDefinition, assuming regionBDefinition==Y")
                regionB_bundled_lpgbthists_allevents.append(phidivisionY_bundled_lpgbthists_allevents)
        else:
            print ("regionBDefinition not given, assuming regionBDefinition==Y")
            regionB_bundled_lpgbthists_allevents.append(phidivisionY_bundled_lpgbthists_allevents)
            
        if truncation_values_method == "original":
            truncation_values.append( getTruncationValuesRoverZ(regionA_bundled_lpgbthists_allevents[-1],regionB_bundled_lpgbthists_allevents[-1],truncation_options[-1]['maxTCsA'],truncation_options[-1]['maxTCsB']) )

        if truncation_values_method == "reverse":
            truncation_values.append( getReverseTruncationValues( regionA_bundled_lpgbthists_allevents[-1], regionB_bundled_lpgbthists_allevents[-1], truncation_options[-1]['maxTCsA'], truncation_options[-1]['maxTCsB']) )


    print ("Using truncation values method " + truncation_values_method)
    for option,truncation in zip(options_to_study,truncation_values):
        print ("Truncation Option " + str(option) + " = ")
        print ( repr(truncation) )
    
    #Once we have the truncation values, need to find how many TCs are lost
    print ("Plotting histograms")
    #Fill a 2D histogram per bunch-crossing with N_TCs (maximum over bundles) 
    TCratios = []
    for option in truncation_options:
        TCratios.append(option['maxTCsA']/option['maxTCsB'])

    if ( len(truncation_options) > 0 ):                
        plot_NTCs_Vs_ROverZ(regionA_bundled_lpgbthists_allevents[-1],inclusive_hists[1],outdir + "/NTCs_Vs_ROverZ_A",truncation_values)
        plot_NTCs_Vs_ROverZ(regionB_bundled_lpgbthists_allevents[-1],inclusive_hists[1],outdir + "/NTCs_Vs_ROverZ_B",truncation_values,TCratios)

    #Plot sum of truncated TCs over the sum of all TCs
    for (study_num,option,values,regionA,regionB) in zip(options_to_study,truncation_options,truncation_values,regionA_bundled_lpgbthists_allevents,regionB_bundled_lpgbthists_allevents):
        plot_frac_Vs_ROverZ( regionA, regionB, values, option['maxTCsA']/option['maxTCsB'], inclusive_hists[1], "Sum No. TCs Option " + str(study_num), outdir + "/frac_option_"+str(study_num))


def plotTruncation(eventData, outdir = ".", useMaximumXY = True, binningConfig = None):
    os.system("mkdir -p " + outdir)

    if ( binningConfig != None ):        
        rOverZMin = binningConfig["rOverZMin"]
        rOverZMax = binningConfig["rOverZMax"]
    else:
        #Set defaults
        rOverZMin = 0.07587128
        rOverZMax = 0.55563514
    
    #Load pickled per-event bundle histograms
    phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents = loadFluctuationData(eventData)

    nBundles = len(phidivisionX_bundled_lpgbthists_allevents[0]) #24 by default
    nROverZBins = len(phidivisionX_bundled_lpgbthists_allevents[0][0]) #42 by default
    
    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = nROverZBins, range = (rOverZMin,rOverZMax) )
    roverzBinning = inclusive_hists[1]
    
    #For each r/z bin take either the sum of the phidivisionX and phidivisionY arrays (inclusive) or the maximum of the two
    if ( useMaximumXY ):
        bundled_lpgbthists_allevents = np.maximum(phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents)
    else:
        bundled_lpgbthists_allevents = phidivisionX_bundled_lpgbthists_allevents + phidivisionY_bundled_lpgbthists_allevents

    #Per event take the maximum in each r/z bin across the nBundles bundles 
    hists_max = np.amax(bundled_lpgbthists_allevents,axis=1)
            
    #Find the maximum per bin over all events,
    #Then find the 99th percentile for a 1% truncation

    overall_max = np.amax(hists_max, axis=0)    
    
    overall_max99p = np.round(np.percentile(hists_max,99,axis=0))
    overall_max95p = np.round(np.percentile(hists_max,95,axis=0))
    overall_max90p = np.round(np.percentile(hists_max,90,axis=0))

    #Loop back over events, counting the maximum wait time
    #for each bin, with and without truncation
    total_per_event = []
    total_per_event99 = []
    total_per_event95 = []
    total_per_event90 = []

    max_per_event_perbin = []
    max_per_event_perbin99 = []
    max_per_event_perbin90 = []
    max_per_event_perbin95 = []
    
    for bundle_hists_phidivisionX,bundle_hists_phidivisionY in zip(phidivisionX_bundled_lpgbthists_allevents, phidivisionY_bundled_lpgbthists_allevents):

        bundle_hists_inclusive = bundle_hists_phidivisionX + bundle_hists_phidivisionY
        bundle_hists_maximum = np.maximum(bundle_hists_inclusive,bundle_hists_phidivisionY)
        #nBundles arrays (24 by default), with length of nROverZBins (42 by default)

        sum99 = []
        sum95 = []
        sum90 = []

        if ( useMaximumXY ):
            bundle_hists = bundle_hists_maximum
        else:
            bundle_hists = bundle_hists_inclusive

        for bundle in bundle_hists:
            
            #If a given r/z bin is greater than the maximum allowed by truncation then set to the truncated value
            sum99.append ( np.where( np.less( overall_max99p, bundle ), overall_max99p, bundle )  )
            sum95.append ( np.where( np.less( overall_max95p, bundle ), overall_max95p, bundle )  )
            sum90.append ( np.where( np.less( overall_max90p, bundle ), overall_max90p, bundle )  )
        

        total_per_event.append( np.sum(bundle_hists, axis=1 )) #array with length of nBundles (24 by default), (sum over the bins - 42 by default)
        total_per_event99.append( np.sum(sum99, axis=1 ))
        total_per_event95.append( np.sum(sum95, axis=1 ))
        total_per_event90.append( np.sum(sum90, axis=1 ))

        max_per_event_perbin.append( np.amax(bundle_hists, axis=0 ) )
        max_per_event_perbin99.append( np.amax(sum99, axis=0 ) )
        max_per_event_perbin95.append( np.amax(sum95, axis=0 ) )
        max_per_event_perbin90.append( np.amax(sum90, axis=0 ) )

    #Calculating the best possible given the per-event fluctuations

    #For a given r/z bin calculate the mean over all events
    #and add 2.5x the average of the nBundles (24 by default) bundles' RMS
    best_likely = np.mean(bundled_lpgbthists_allevents,axis=(0,1))+2.5*np.mean(np.std(bundled_lpgbthists_allevents,axis=(0)),axis=0)
    ratio_to_best = np.divide(overall_max99p,best_likely,out=np.zeros_like(overall_max99p),where=best_likely!=0)
    
    print ("Maximum TC in any bundle in any event (per bin) = ", np.round(np.amax(max_per_event_perbin,axis=0)))
    print ("Sum of per-bin maximum TC (over bundles and events) = ",  np.round(np.sum(np.amax(max_per_event_perbin,axis=0))))
    print ("Sum of per-bin maximum TC (over bundles and events) with 1% truncation =", np.round(np.sum(np.amax(max_per_event_perbin99,axis=0))))
    print ("Sum of per-bin maximum TC (over bundles and events) with 5% truncation = ", np.round(np.sum(np.amax(max_per_event_perbin95,axis=0))))
    print ("Sum of per-bin maximum TC (over bundles and events) with 10% truncation = ", np.round(np.sum(np.amax(max_per_event_perbin90,axis=0))))

    pl.hist(np.sum(np.array(total_per_event)-np.array(total_per_event99),axis=1)/(nBundles),50,(0,5),histtype='step',log=True,label='1% truncation')
    pl.hist(np.sum(np.array(total_per_event)-np.array(total_per_event95),axis=1)/(nBundles),50,(0,5),histtype='step',log=True,label='5% truncation')
    pl.hist(np.sum(np.array(total_per_event)-np.array(total_per_event90),axis=1)/(nBundles),50,(0,5),histtype='step',log=True,label='10% truncation')    
    pl.xlabel('Number of TCs truncated on average per bundle')
    
    pl.ylabel('Number of Events')
    pl.legend()
    pl.savefig( outdir + "/truncation.png" )

    pl.clf()
    pl.step(roverzBinning, np.append(ratio_to_best,ratio_to_best[-1]), where='post')
    pl.axhline(y=1, color='r', linestyle='--')
    pl.xlabel('r/z')
    pl.ylabel('Ratio of 1% truncation to likely best')
    pl.ylim((0,10))
    pl.savefig( outdir + "/ratio_to_best.png" )

    #As a cross check plot the bundle R/Z histograms integrated over all events.
    #These should be the same as those produced by plotbundles.py
    pl.clf()
    for bundle in np.sum(phidivisionX_bundled_lpgbthists_allevents,axis=0):
        pl.step(roverzBinning, np.append(bundle,bundle[-1]), where='post')
    pl.ylim((0,1100000))
    pl.savefig( outdir + "/phidivisionXIntegrated.png" )
    pl.clf()
    for bundle in np.sum(phidivisionY_bundled_lpgbthists_allevents,axis=0):
        pl.step(roverzBinning, np.append(bundle,bundle[-1]), where='post')
    pl.ylim((0,1100000))
    pl.savefig( outdir + "/phidivisionYIntegrated.png" )

def plot_Truncation_tc_Pt(eventData, options_to_study, truncationConfig = None, binningConfig = None, outdir = "."):
    os.system("mkdir -p " + outdir)
    
    if ( binningConfig != None ):        
        rOverZMin = binningConfig["rOverZMin"]
        rOverZMax = binningConfig["rOverZMax"]
    else:
        #Set defaults
        rOverZMin = 0.07587128
        rOverZMax = 0.55563514
    
    #Load the per-event flucation data produced using 'checkFluctuations'
    with open(eventData, "rb") as filep:   
        data = pickle.load(filep)

    #Identifiers for regionA and regionB
    dataA = 0
    dataB = 1

    nROverZBins = len(data[0][0][0]) #42 by default
    inclusive_hists =  np.histogram( np.empty(0), bins = nROverZBins, range = (rOverZMin,rOverZMax) )[1]

    truncation_options_regionA = []
    truncation_options_regionB = []

    #Get Np arrays for regions A and B and for each truncation option
    #loop over number of truncation options

    for t in range(len(data[0])):
        
        dataA_allevents = np.empty((len(data),nROverZBins))
        dataB_allevents = np.empty((len(data),nROverZBins)) 

        for e,event in enumerate(data):        
            dataA_allevents[e] = np.asarray( event[t][dataA] )
            dataB_allevents[e] = np.asarray( event[t][dataB] )
        
        truncation_options_regionA.append(dataA_allevents)
        truncation_options_regionB.append(dataB_allevents)

    #Sum over all events and bundles of TCs (in each R/Z bin) 
    totalsumA = np.sum( truncation_options_regionA[0] , axis=0 )
    totalsumB = np.sum( truncation_options_regionB[0] , axis=0 )

    #Loop over truncation options
    for t,(truncationA,truncationB) in enumerate(zip(truncation_options_regionA,truncation_options_regionB)):
        if t == 0: continue #total sum
        if (t > len(options_to_study) ):
            break
        
        #Sum over all events of truncated TCs (in each R/Z bin) 
        truncatedsum_A = np.sum( truncationA, axis=0 )
        truncatedsum_B = np.sum( truncationB, axis=0 )

        #Divide to get the fraction, taking into account division by zero
        #We assume that the order of options in the input-data is the same as
        #the order of options provided in the config
        ratioA = np.divide(   truncatedsum_A, totalsumA , out=np.ones_like(truncatedsum_A), where=totalsumA!=0 )
        ratioB = np.divide(   truncatedsum_B, totalsumB , out=np.ones_like(truncatedsum_B), where=totalsumB!=0 )
        
        pl.clf()
        pl.step(inclusive_hists,np.append( ratioA , ratioA[-1] ),color='red',linewidth='1', where = 'post', label='data A')
        pl.step(inclusive_hists,np.append( ratioB , ratioB[-1] ),color='orange',linewidth='1', where = 'post', label='data B')

        pl.xlabel('r/z')
        pl.ylabel('pT sum TCs with truncation / pT sum all TCs')
        pl.title("Sum pT TCs Option " + str(options_to_study[t-1]) )
        pl.ylim((0.6,1.05))
        pl.legend()
        pl.savefig( outdir + "/pt_truncation_option_" + str(options_to_study[t-1]) + ".png" )
        pl.clf()
    
