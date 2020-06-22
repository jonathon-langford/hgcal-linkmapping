#!/usr/bin/env python3
import sys
sys.path.insert(1, './externals')
import ROOT
import numpy as np
import mlrose_mod as mlrose
import time
import yaml
import signal

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from sklearn.metrics import accuracy_score

from process import getModuleHists,getlpGBTHists,getMiniGroupHists,getMinilpGBTGroups,getBundles, getBundledlpgbtHists,getBundledlpgbtHistsRoot,calculateChiSquared
from process import loadDataFile,getTCsPassing,getlpGBTLoadInfo,getHexModuleLoadInfo
from process2D import getModuleHists2D, getlpGBTHists2D, getMiniGroupHists2D, getMiniGroupHists2D, getBundledlpgbtHists2D, getBundledlpgbtHistsRoot2D, calculateChiSquared_modif
from plotting import plot, plot2D
from bestchi2 import bestsofar

chi2_min = 50000000000000000000000
combbest = []

def handler(signum, frame):
    raise Exception()    

def plot_lpGBTLoads(MappingFile,CMSSW_Silicon,CMSSW_Scintillator):

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW
    
    lpgbt_loads_tcs,lpgbt_loads_words,lpgbt_layers = getlpGBTLoadInfo(data,data_tcs_passing,data_tcs_passing_scin)

    plot(lpgbt_loads_tcs,"loads_tcs.png",binwidth=0.1,xtitle='Number of TCs on a single lpGBT')
    plot(lpgbt_loads_words,"loads_words.png",binwidth=0.1,xtitle='Number of words on a single lpGBT')
    plot2D(lpgbt_loads_tcs,lpgbt_layers,"tcs_vs_layer.png",xtitle='Number of TCs on a single lpGBT')
    plot2D(lpgbt_loads_words,lpgbt_layers,"words_vs_layer.png",xtitle='Number of words on a single lpGBT')

def plot_ModuleLoads(MappingFile,CMSSW_Silicon,CMSSW_Scintillator):

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW
    
    lpgbt_loads_tcs,lpgbt_loads_words,lpgbt_layers = getHexModuleLoadInfo(data,data_tcs_passing,data_tcs_passing_scin)

    plot(lpgbt_loads_tcs,"loads_tcs.png",binwidth=0.1,xtitle='Number of TCs on a single lpGBT')
    plot(lpgbt_loads_words,"loads_words.png",binwidth=0.1,xtitle='Number of words on a single lpGBT')
    plot2D(lpgbt_loads_tcs,lpgbt_layers,"tcs_vs_layer.png",xtitle='Number of TCs on a single lpGBT')
    plot2D(lpgbt_loads_words,lpgbt_layers,"words_vs_layer.png",xtitle='Number of words on a single lpGBT')

def check_for_missing_modules_inMappingFile(MappingFile,CMSSW_Silicon,CMSSW_Scintillator):

    #Check for modules missing in the mapping file
    
    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW
    
    mappingfile_sil = data[data['density']<2][['layer', 'u', 'v']]
    mappingfile_scin = data[data['density']==2][['layer', 'u', 'v']]

    cmssw_sil = data_tcs_passing[['u','v','layer','nTCs']]
    cmssw_scin = data_tcs_passing_scin[['u','v','layer','nTCs']]

    #onlymapping_sil = mappingfile.merge(cmssw.drop_duplicates(), on=['u','v','layer'],how='left', indicator=True)
    onlycmssw_sil = cmssw_sil.merge(mappingfile_sil.drop_duplicates(), on=['u','v','layer'],how='left', indicator=True)
    onlycmssw_scin = cmssw_scin.merge(mappingfile_scin.drop_duplicates(), on=['u','v','layer'],how='left', indicator=True)

    onlycmssw_sil = onlycmssw_sil[onlycmssw_sil['_merge'] == 'left_only']
    onlycmssw_scin = onlycmssw_scin[onlycmssw_scin['_merge'] == 'left_only']

    print ("Silicon")
    print (onlycmssw_sil[onlycmssw_sil['nTCs']>0][['layer','u','v']].to_string(index=False))
    print ("Scintillator")
    print (onlycmssw_scin[onlycmssw_scin['nTCs']>0][['layer','u','v']].to_string(index=False))

def check_for_missing_modules_inCMSSW(MappingFile,CMSSW_Silicon,CMSSW_Scintillator):

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW
    getHexModuleLoadInfo(data,data_tcs_passing,data_tcs_passing_scin,True)
    
    
def study_mapping(MappingFile,CMSSW_ModuleHists,algorithm="random_hill_climb",initial_state="best_so_far",random_seed=1,max_iterations=100000,output_dir=".",print_level=0, minigroup_type="minimal",chi2_coef=0.1):

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
   # inclusive_hists,module_hists = getModuleHists(CMSSW_ModuleHists)
   # 
   # #Form hists corresponding to each lpGBT from module hists
   # lpgbt_hists = getlpGBTHists(data, module_hists)
   # minigroups, minigroups_swap = getMinilpGBTGroups(data, minigroup_type)
   # minigroup_hists = getMiniGroupHists(lpgbt_hists,minigroups_swap)
   # minigroup_hists_root = getMiniGroupHists(lpgbt_hists,minigroups_swap,root=True)

    inclusive_hists, module_hists2D = getModuleHists2D(CMSSW_ModuleHists)
    lpgbt_hists2D = getlpGBTHists2D(data, module_hists2D)
    minigroups, minigroups_swap = getMinilpGBTGroups(data, minigroup_type)
    minigroup_hists2D = getMiniGroupHists2D(lpgbt_hists2D, minigroups_swap)
    minigroup_hists_root2D = getMiniGroupHists2D(lpgbt_hists2D,minigroups_swap,root=True)


    def mapping_max(state):
        global chi2_min
        global combbest

        chi2 = 0
    
        bundles = getBundles(minigroups_swap,state)
        bundled_lpgbthists = getBundledlpgbtHists2D(minigroup_hists2D,bundles)

        chi2 = calculateChiSquared_modif(inclusive_hists,bundled_lpgbthists,chi2_coef,root=False)

        typicalchi2 = 600000000000
        if (chi2<chi2_min):
            chi2_min = chi2
            combbest = np.copy(state)
            if ( print_level > 0 ):
                print (algorithm," ", chi2_min, " ", chi2_min/typicalchi2)
            if ( print_level > 1 ):
                print (repr(combbest))

        return chi2

    
    init_state = []
    if (initial_state == "best_so_far"):
        init_state = bestsofar
    if (initial_state[-4:] == ".npy"):
        print (initial_state)
        init_state = np.hstack(np.load(initial_state,allow_pickle=True))
    elif (initial_state == "random"):
       init_state = np.arange(len(minigroups_swap))
       np.random.shuffle(init_state)

       
    fitness_cust = mlrose.CustomFitness(mapping_max)
    # Define optimization problem object
    problem_cust = mlrose.DiscreteOpt(length = len(init_state), fitness_fn = fitness_cust, maximize = False, max_val = len(minigroups_swap), minigroups = minigroups_swap)

    # Define decay schedule
    schedule = mlrose.ExpDecay()
    #schedule = mlrose.ArithDecay()

    filename = "bundles_job_"
    filenumber = ""
    if ( len(sys.argv) > 2 ):
        filenumber = str(sys.argv[2])
    else:
        filenumber = "default"
    filename+=filenumber
        
    if ( algorithm == "save_root" ):
        #Save best combination so far into a root file
        bundles = getBundles(minigroups_swap,init_state)

        bundled_hists = getBundledlpgbtHistsRoot2D(minigroup_hists_root2D,bundles)
        chi2 = calculateChiSquared_modif(inclusive_hists,bundled_hists,chi2_coef,root=True)
        newfile = ROOT.TFile("lpgbt_10.root","RECREATE")
        np.save(output_dir + "/" + filename + ".npy",bundles)
        for sector in bundled_hists:
            for key, value in sector.items():
                value.Write()
        for sector in inclusive_hists:
            sector.Scale(1./24.)
            sector.Write()
        newfile.Close()
        print ("Chi2:",chi2)
        print ("List of Bundles:")
        for b,bundle in enumerate(bundles):
            print ("" )
            print ("bundle" + str(b) )
            for minigroup in bundle:
                #print (minigroup)
                lpgbts = minigroups_swap[minigroup]
                for lpgbt in lpgbts:
                    print (str(lpgbt) + ", "  , end = '')

        
    elif (algorithm == "random_hill_climb"):
        try:

            best_state, best_fitness = mlrose.random_hill_climb(problem_cust, max_attempts=10000, max_iters=max_iterations, restarts=0, init_state=init_state, random_state=random_seed)
            print (repr(best_state))

        except:
            print("interrupt received, stopping and saving")
            

        finally:
            bundles = getBundles(minigroups_swap,combbest)
            np.save(output_dir + "/" + filename + ".npy",bundles)
            file1 = open(output_dir + "/chi2.txt","a")
            file1.write( "bundles[" + filenumber + "] = " + str(chi2_min) + "\n" )
            file1.close( )

            
    elif (algorithm == "genetic_alg"):
        best_state, best_fitness = mlrose.genetic_alg(problem_cust, pop_size=200, mutation_prob=0.1, max_attempts=1000, max_iters=10000000, curve=False, random_state=random_seed)
    elif (algorithm == "mimic"):
        best_state, best_fitness = mlrose.mimic(problem_cust, pop_size=200,  keep_pct=0.2, max_attempts=10, max_iters=np.inf, curve=False, random_state=random_seed)
    elif (algorithm == "simulated_annealing"):
        best_state, best_fitness = mlrose.simulated_annealing(problem_cust, schedule = schedule, max_attempts = 100000, max_iters = 10000000, init_state = init_state, random_state = 1)
    else:
        print("Algorithm "+ algorithm + " not known" )


    
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

    #Catch possible exceptions from batch system
    signal.signal(signal.SIGINT,handler)
    signal.signal(signal.SIGUSR1,handler)
    signal.signal(signal.SIGXCPU,handler)

    if ( config['function']['study_mapping'] ):
        subconfig = config['study_mapping']
        study_mapping(subconfig['MappingFile'],subconfig['CMSSW_ModuleHists'],algorithm=subconfig['algorithm'],initial_state=subconfig['initial_state'],random_seed=subconfig['random_seed'],max_iterations=subconfig['max_iterations'],output_dir=config['output_dir'],print_level=config['print_level'],minigroup_type=subconfig['minigroup_type'],chi2_coef=subconfig['chi2_coef'])

    if ( config['function']['check_for_missing_modules'] ):
        subconfig = config['check_for_missing_modules']
        if ( subconfig['inMappingFile'] ):
            print("Missing modules in mapping file: "+ subconfig['MappingFile'] + "\n")
            check_for_missing_modules_inMappingFile(subconfig['MappingFile'],subconfig['CMSSW_Silicon'],subconfig['CMSSW_Scintillator'])
        if ( subconfig['inCMSSW'] ):
            print("\nMissing modules in CMSSW\n")
            check_for_missing_modules_inCMSSW(subconfig['MappingFile'],subconfig['CMSSW_Silicon'],subconfig['CMSSW_Scintillator'])

    if ( config['function']['plot_lpGBTLoads'] ):
        subconfig = config['plot_lpGBTLoads']
        plot_lpGBTLoads(subconfig['MappingFile'],subconfig['CMSSW_Silicon'],subconfig['CMSSW_Scintillator'])

    if ( config['function']['plot_ModuleLoads'] ):
        subconfig = config['plot_ModuleLoads']
        plot_ModuleLoads(subconfig['MappingFile'],subconfig['CMSSW_Silicon'],subconfig['CMSSW_Scintillator'])

    
main()
