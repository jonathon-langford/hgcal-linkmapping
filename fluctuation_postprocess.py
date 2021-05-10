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

def plotMeanMax(eventData, outdir = ".", includePhi60 = True):
    #Load pickled per-event bundle histograms
    with open(eventData, "rb") as filep:   
        bundled_lpgbthists_allevents = pickle.load(filep)
    os.system("mkdir -p " + outdir)
    
    nbins = 42
    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = nbins, range = (0.076,0.58) )

    #Names for inclusive and phi < 60 indices
    inclusive = 0
    phi60 = 1

    hists_max = [] 
    
    #Plotting Max, mean and standard deviation per bundle:

    for bundle in range(24):

        list_over_events_inclusive = np.empty(((len(bundled_lpgbthists_allevents)),nbins))
        list_over_events_phi60 = np.empty(((len(bundled_lpgbthists_allevents)),nbins))
        
        for e,event in enumerate(bundled_lpgbthists_allevents):
            list_over_events_inclusive[e] = np.array(event[inclusive][bundle])/6
            list_over_events_phi60[e] = np.array(event[phi60][bundle])/6

        list_over_events_maximum = np.maximum(list_over_events_inclusive, list_over_events_phi60*2 )

        if ( includePhi60 ):
            list_over_events = list_over_events_maximum
        else:
            list_over_events = list_over_events_inclusive

        hist_max = np.amax(list_over_events,axis=0)
        hist_mean = np.mean(list_over_events, axis=0)
        hist_std = np.std(list_over_events, axis=0)

        for s,std in enumerate(hist_std):
            hist_std[s] = std + hist_mean[s]

        pl.bar((inclusive_hists[1])[:-1], hist_max, width=0.012,align='edge')
        pl.bar((inclusive_hists[1])[:-1], hist_std, width=0.012,align='edge')
        pl.bar((inclusive_hists[1])[:-1], hist_mean, width=0.012,align='edge')

        #Plot all events for a given bundle on the same plot
        # for e,event in enumerate(list_over_events):
        #     pl.bar((inclusive_hists[1])[:-1], event, width=0.012,fill=False)
        #     #if (e>200): break

        pl.ylim((0,31))
        pl.savefig( outdir + "/bundle_" + str(bundle) + "max.png" )        
        pl.clf()

        hists_max.append(hist_max)
        
    #Plot maxima for all bundles on the same plot
    for hist in hists_max:
        pl.bar((inclusive_hists[1])[:-1], hist, width=0.012,align='edge')
    pl.ylim((0,31))
    pl.xlabel('r/z')
    pl.ylabel('Maximum number of TCs per bin')
    pl.savefig( outdir + "/maxima.png" )
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
    pl.clf()
    pl.hist2d( roverz , data_swap.flatten() , bins = (len(axis)-1,50),range=[[0.076,0.58], [0, 50]],norm=LogNorm())
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
    pl.ylabel('Sum truncated TCs / Sum all TCs')
    pl.title(title)
    pl.ylim((0.6,1.05))
    pl.legend()
    pl.savefig( savename + ".png" )
    pl.clf()    

def getTruncationValuesForGivenScalar(scalar,data):

    #For the get reverse truncation values method
    #For a given "Scalar" (i.e. efficiency or sum of truncated TCs / sum of all TCs)
    
    #Sum over all events and bundles of TCs (in each R/Z bin) 
    totalsumA = np.sum( data[0] , axis=1 )
    totalsumB = np.sum( data[1] , axis=1 )

    #Largest possible truncation value (end of loop)
    largestTruncationMax = 100

    #Loop over bins
    truncation_values = []

    for b in range(len(data[0])):
        scalarset = False
        for x in range (0,largestTruncationMax,1):
            truncation_sum = np.sum(np.where(data[0][b]<x, data[0][b], x))

            if ( scalar*totalsumA[b] == 0):
                truncation_values.append(0)
                scalarset = True
                break

            if truncation_sum > scalar*totalsumA[b]:#break when the truncated sum is below the desired value

                #Get the truncation sum for x-1 and check whether this or the sum for x is closer to giving the desired efficiency
                truncation_sum_x_minus_1 = np.sum(np.where(data[0][b]<x-1, data[0][b], x-1))
                if ( abs( scalar*totalsumA[b] - truncation_sum ) < abs( scalar*totalsumA[b] - truncation_sum_x_minus_1 )):
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

    truncation_values = getTruncationValuesForGivenScalar(scalar,data)
    difference = abs(data[2]-np.sum(truncation_values))
        
    return difference
    
def getReverseTruncationValues( dataA, dataB, maxTCsA, maxTCsB ):
    #Obtain the truncation values by demanding that the efficiency be constant as a function of r/z
    
    nbins = len(dataA[0][0]) #42

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

    #Give 'spare' TCs to low bins
    truncation_values = np.array(truncation_values)
    n_spare = maxTCsA - np.sum(truncation_values)
    if ( n_spare != 0):
        arg = np.argsort(truncation_values)[::np.sign(n_spare)][:int(abs(n_spare))]
        truncation_values[arg]+=np.sign(n_spare)
        
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
    #For 3-link options RegionA = phidivisionX+phidivisionY, for the 4-link
    #options RegionA = phidivisionX. In both cases RegionB = phidivisionY
    
    with open(eventData, "rb") as filep:   
        bundled_lpgbthists_allevents = pickle.load(filep)
    
    #Names for phi > split_value and phi < split_value indices
    dataX = 0
    dataY = 1

    nbundles = len(bundled_lpgbthists_allevents[0][0]) #24
    nbins = len(bundled_lpgbthists_allevents[0][0][0]) #42
    
    dataX_bundled_lpgbthists_allevents = np.empty((len(bundled_lpgbthists_allevents),nbundles,nbins))
    dataY_bundled_lpgbthists_allevents = np.empty((len(bundled_lpgbthists_allevents),nbundles,nbins))

    for e,event in enumerate(bundled_lpgbthists_allevents):        
        dataX_bundled_lpgbthists_allevents[e] = np.array(list(event[dataX].values()))
        dataY_bundled_lpgbthists_allevents[e] = np.array(list(event[dataY].values()))

    return dataX_bundled_lpgbthists_allevents, dataY_bundled_lpgbthists_allevents
    
def studyTruncationOptions(eventData, options_to_study, truncation_values_method, truncationConfig, outdir = "."):
    
    #Load pickled per-event bundle histograms
    phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents = loadFluctuationData(eventData)

    os.system("mkdir -p " + outdir)

    nbinsROverZ = len(phidivisionX_bundled_lpgbthists_allevents[0][0]) #42
    
    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = nbinsROverZ, range = (0.076,0.58) )

    inclusive_bundled_lpgbthists_allevents = phidivisionX_bundled_lpgbthists_allevents + phidivisionY_bundled_lpgbthists_allevents

    #RegionA is either phidivisionX (in the case of options 4 and 5) or phidivisionX+phidivisionY (in the cases of 1, 2 and 3)
    
    truncation_values = []
    truncation_options = []
    regionA_bundled_lpgbthists_allevents = []
    regionB_bundled_lpgbthists_allevents = []

    for option in options_to_study:
        print ("Get truncation value for option " + str(option))

        truncation_options.append(truncationConfig['option'+str(option)])
    
        if truncation_options[-1]['nLinks'] == 3:            
            regionA_bundled_lpgbthists_allevents.append(inclusive_bundled_lpgbthists_allevents)
        elif truncation_options[-1]['nLinks'] == 4:
            regionA_bundled_lpgbthists_allevents.append(phidivisionX_bundled_lpgbthists_allevents)
            
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

    #If options 1,2 or 3 (3-links) or 4,5 (4-links) are included in options_to_study
    options_3links = []
    options_3links_TCratio = []
    options_4links = []
    options_4links_TCratio = []
    for option,truncation in zip(truncation_options,truncation_values):
        if option['nLinks'] == 3:
            options_3links.append(truncation)
            options_3links_TCratio.append(option['maxTCsA']/option['maxTCsB'])
        elif option['nLinks'] == 4:
            options_4links.append(truncation)
            options_4links_TCratio.append(option['maxTCsA']/option['maxTCsB'])

    if ( len(options_3links) > 0 ):                
        plot_NTCs_Vs_ROverZ(inclusive_bundled_lpgbthists_allevents,inclusive_hists[1],outdir + "/NTCs_Vs_ROverZ_A_3links",options_3links)
        plot_NTCs_Vs_ROverZ(phidivisionY_bundled_lpgbthists_allevents,inclusive_hists[1],outdir + "/NTCs_Vs_ROverZ_B_3links",options_3links,options_3links_TCratio)

    if ( len(options_4links) > 0 ):                
        #Don't want inclusive here to be region A, rather phidivisionX
        plot_NTCs_Vs_ROverZ(phidivisionX_bundled_lpgbthists_allevents,inclusive_hists[1],outdir + "/NTCs_Vs_ROverZ_A_4links",options_4links)
        plot_NTCs_Vs_ROverZ(phidivisionY_bundled_lpgbthists_allevents,inclusive_hists[1],outdir + "/NTCs_Vs_ROverZ_B_4links",options_4links,options_4links_TCratio)

    #Plot sum of truncated TCs over the sum of all TCs
    for (study_num,option,values,regionA,regionB) in zip(options_to_study,truncation_options,truncation_values,regionA_bundled_lpgbthists_allevents,regionB_bundled_lpgbthists_allevents):
        plot_frac_Vs_ROverZ( regionA, regionB, values, option['maxTCsA']/option['maxTCsB'], inclusive_hists[1], "Sum No. TCs Option " + str(study_num), outdir + "/frac_option_"+str(study_num))


def plotTruncation(eventData, outdir = ".", includePhi60 = True):
    os.system("mkdir -p " + outdir)
    
    #Load pickled per-event bundle histograms
    phidivisionX_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents = loadFluctuationData(eventData)

    #To get binning for r/z histograms
    inclusive_hists = np.histogram( np.empty(0), bins = 42, range = (0.076,0.58) )
    roverzBinning = inclusive_hists[1]
    
    #Form the intersection of the inclusive and phi60 arrays,
    #taking for each bin the maximum of the inclusive and phidivisionY x 2
    inclusive_bundled_lpgbthists_allevents = phidivisionX_bundled_lpgbthists_allevents + phidivisionY_bundled_lpgbthists_allevents
    maximum_bundled_lpgbthists_allevents = np.maximum(inclusive_bundled_lpgbthists_allevents,phidivisionY_bundled_lpgbthists_allevents*2)
    
    if ( includePhi60 ):
        hists_max = np.amax(maximum_bundled_lpgbthists_allevents,axis=1)
    else:
        hists_max = np.amax(inclusive_bundled_lpgbthists_allevents,axis=1)
            
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
        bundle_hists_maximum = np.maximum(bundle_hists_inclusive,bundle_hists_phidivisionY*2)
        #24 arrays, with length of 42        

        sum99 = []
        sum95 = []
        sum90 = []

        if ( includePhi60 ):
            bundle_hists = bundle_hists_maximum
        else:
            bundle_hists = bundle_hists_inclusive

        for bundle in bundle_hists:
            
            #If a given r/z bin is greater than the maximum allowed by truncation then set to the truncated value
            sum99.append ( np.where( np.less( overall_max99p, bundle ), overall_max99p, bundle )  )
            sum95.append ( np.where( np.less( overall_max95p, bundle ), overall_max95p, bundle )  )
            sum90.append ( np.where( np.less( overall_max90p, bundle ), overall_max90p, bundle )  )
        

        total_per_event.append( np.sum(bundle_hists, axis=1 ))        #array with length of 24 (sum over the 42 bins)
        total_per_event99.append( np.sum(sum99, axis=1 ))
        total_per_event95.append( np.sum(sum95, axis=1 ))
        total_per_event90.append( np.sum(sum90, axis=1 ))

        max_per_event_perbin.append( np.amax(bundle_hists, axis=0 ) )
        max_per_event_perbin99.append( np.amax(sum99, axis=0 ) )
        max_per_event_perbin95.append( np.amax(sum95, axis=0 ) )
        max_per_event_perbin90.append( np.amax(sum90, axis=0 ) )

    #Calculating the best possible given the per-event fluctuations

    #For a given r/z bin calculate the mean over all events
    #and add 2.5x the average of the 24 bundles' RMS
    best_likely = np.mean(inclusive_bundled_lpgbthists_allevents,axis=(0,1))+2.5*np.mean(np.std(inclusive_bundled_lpgbthists_allevents,axis=(0)),axis=0)
    ratio_to_best = np.divide(overall_max99p,best_likely,out=np.zeros_like(overall_max99p),where=best_likely!=0)
    
    print ("Maximum TC in any bundle in any event (per bin) = ", np.round(np.amax(max_per_event_perbin,axis=0)))
    print ("Sum of per-bin maximum TC (over bundles and events) = ",  np.round(np.sum(np.amax(max_per_event_perbin,axis=0))))
    print ("Sum of per-bin maximum TC (over bundles and events) with 1% truncation =", np.round(np.sum(np.amax(max_per_event_perbin99,axis=0))))
    print ("Sum of per-bin maximum TC (over bundles and events) with 5% truncation = ", np.round(np.sum(np.amax(max_per_event_perbin95,axis=0))))
    print ("Sum of per-bin maximum TC (over bundles and events) with 10% truncation = ", np.round(np.sum(np.amax(max_per_event_perbin90,axis=0))))

    pl.hist(np.sum(np.array(total_per_event)-np.array(total_per_event99),axis=1)/(24),50,(0,5),histtype='step',log=True,label='1% truncation')
    pl.hist(np.sum(np.array(total_per_event)-np.array(total_per_event95),axis=1)/(24),50,(0,5),histtype='step',log=True,label='5% truncation')
    pl.hist(np.sum(np.array(total_per_event)-np.array(total_per_event90),axis=1)/(24),50,(0,5),histtype='step',log=True,label='10% truncation')    
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

def plot_Truncation_tc_Pt(eventData, options_to_study, truncationConfig = None, outdir = ".",  ):
    os.system("mkdir -p " + outdir)
    
    #Load the per-event flucation data produced using 'checkFluctuations'
    with open(eventData, "rb") as filep:   
        data = pickle.load(filep)

    #Load the number of links for the truncation options
    nLinks = []
    if ( truncationConfig != None ):
        for option in options_to_study:
            nLinks.append(truncationConfig['option'+str(option)]['nLinks'])

    nbinsROverZ = len(data[0][0][0]) #42
    axis =  np.histogram( np.empty(0), bins = nbinsROverZ, range = (0.076,0.58) )[1]

    truncation_options_regionA = []
    truncation_options_regionB = []

    #Get Np arrays for regions A and B and for each truncation option
    #loop over number of truncation options

    for t in range(len(data[0])):
        
        dataA_allevents = np.empty((len(data),nbinsROverZ))
        dataB_allevents = np.empty((len(data),nbinsROverZ)) 

        for e,event in enumerate(data):        
            #if e>4: continue
            dataA_allevents[e] = np.asarray( event[t][0] )
            dataB_allevents[e] = np.asarray( event[t][1] )
        
        truncation_options_regionA.append(dataA_allevents)
        truncation_options_regionB.append(dataB_allevents)


    #Sum over all events and bundles of TCs (in each R/Z bin) 
    totalsumA = np.sum( truncation_options_regionA[0] , axis=0 )
    totalsumB = np.sum( truncation_options_regionB[0] , axis=0 )
    totalsumInclusive = totalsumA + totalsumB

    #Loop over truncation options
    for t,(truncationA,truncationB) in enumerate(zip(truncation_options_regionA,truncation_options_regionB)):
        if t == 0: continue #total sum
        
        #Sum over all events of truncated TCs (in each R/Z bin) 
        truncatedsum_A = np.sum( truncationA, axis=0 )
        truncatedsum_B = np.sum( truncationB, axis=0 )

        #Divide to get the fraction, taking into account division by zero
        #We assume that the order of options in the input-data is the same as
        #the order of options provided in the config
        if (nLinks[t-1] == 3 ):
            ratioA = np.divide(   truncatedsum_A, totalsumInclusive , out=np.ones_like(truncatedsum_A), where=totalsumInclusive!=0 )
        elif (nLinks[t-1] == 4 ):
            ratioA = np.divide(   truncatedsum_A, totalsumA , out=np.ones_like(truncatedsum_A), where=totalsumA!=0 )
        ratioB = np.divide(   truncatedsum_B, totalsumB , out=np.ones_like(truncatedsum_B), where=totalsumB!=0 )
        
        pl.clf()
        pl.step(axis,np.append( ratioA , ratioA[-1] ),color='red',linewidth='1', where = 'post', label='data A')
        pl.step(axis,np.append( ratioB , ratioB[-1] ),color='orange',linewidth='1', where = 'post', label='data B')

        pl.xlabel('r/z')
        pl.ylabel('pT sum truncated TCs / pT sum all TCs')
        pl.title("Sum pT TCs Option " + str(options_to_study[t-1]) )
        pl.ylim((0.6,1.05))
        pl.legend()
        pl.savefig( outdir + "/pt_truncation_option_" + str(options_to_study[t-1]) + ".png" )
        pl.clf()
    
