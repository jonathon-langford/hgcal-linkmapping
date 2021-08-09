#!/usr/bin/env python3
import sys
sys.path.insert(1, './externals')
import ROOT
import numpy as np
import mlrose_mod as mlrose # Author: Genevieve Hayes https://github.com/gkhayes/mlrose/tree/master/mlrose
import yaml
import signal
import pickle
import json
import re
import subprocess
import random

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from sklearn.metrics import accuracy_score
from _ctypes import PyObj_FromPtr

from process import getModuleHists, getCMSSWNtupleName, getlpGBTHists, getMiniGroupHists, getMinilpGBTGroups, getMiniModuleGroups, getBundles, getBundledlpgbtHists, getBundledlpgbtHistsRoot, calculateChiSquared, getMaximumNumberOfModulesInABundle
from process import loadDataFile, loadModuleTowerMappingFile, loadConfiguration, getTCsPassing, getlpGBTLoadInfo, getHexModuleLoadInfo, getModuleTCHists, getMiniTowerGroups, getMaxTowersList
from plotting import plot, plot2D

from geometryCorrections import applyGeometryCorrections

chi2_min = 50000000000000000000000
combbest = []
nCallsToMappingMax = 0

class exitProgramSignal(LookupError):
    pass
 
def handler(signum, frame):
    raise exitProgramSignal()    

def dummy_handler(signum, frame):
    pass

def plot_lpGBTLoads(Configuration):

    subconfig = Configuration['plot_lpGBTLoads']
    output_dir = Configuration['output_dir']
    MappingFile = subconfig['MappingFile']
    CMSSW_Silicon = subconfig['CMSSW_Silicon']
    CMSSW_Scintillator = subconfig['CMSSW_Scintillator']

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW
    
    lpgbt_loads_tcs,lpgbt_loads_words,lpgbt_layers = getlpGBTLoadInfo(data,data_tcs_passing,data_tcs_passing_scin)

    plot(lpgbt_loads_tcs,"loads_tcs.png",binwidth=0.1,xtitle='Number of TCs on a single lpGBT',outdir=output_dir)
    plot(lpgbt_loads_words,"loads_words.png",binwidth=0.1,xtitle='Number of words on a single lpGBT',outdir=output_dir)
    plot2D(lpgbt_loads_tcs,lpgbt_layers,"tcs_vs_layer.png",xtitle='Number of TCs on a single lpGBT',outdir=output_dir)
    plot2D(lpgbt_loads_words,lpgbt_layers,"words_vs_layer.png",xtitle='Number of words on a single lpGBT',outdir=output_dir)

def plot_ModuleLoads(Configuration):

    subconfig = Configuration['plot_ModuleLoads']
    output_dir = Configuration['output_dir']
    MappingFile = subconfig['MappingFile']
    CMSSW_Silicon = subconfig['CMSSW_Silicon']
    CMSSW_Scintillator = subconfig['CMSSW_Scintillator']

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW    
    module_loads_words,module_layers,u,v = getHexModuleLoadInfo(data,data_tcs_passing,data_tcs_passing_scin)

    plot(module_loads_words,"module_loads_words.png",binwidth=0.01,xtitle=r'Average number of words on a single module / $2 \times N_{e-links}$',outdir=output_dir)
    plot2D(module_loads_words,module_layers,"module_words_vs_layer.png",binwidthx=0.05,binwidthy=1,xtitle=r'Average number of words on a single module / $2 \times N_{e-links}$',outdir=output_dir)
    
def produce_AllocationFile(Configuration):

    subconfig = Configuration['produce_AllocationFile']
    allocation = subconfig['allocation']
    file_name = subconfig['file_name']
    
    #Load allocation information
    info = loadConfiguration(allocation)
    data = info['data']
    mapping = info['mapping']
    minigroup_type = info['minigroup_type']
    nBundles = info['nBundles']
    maxInputs = info['maxInputs']

    #Get minigroups
    minigroups,minigroups_swap = getMinilpGBTGroups(data, minigroup_type)

    #Bundle together minigroup configuration
    bundles = getBundles(minigroups_swap,mapping,nBundles,maxInputs)

    #Open output file
    fileout = open(file_name, 'w')
    fileout.write( '(lpGBT_number) (number_modules) (sil=0scin=1) (layer) (u/eta) (v/phi) (number_elinks)\n' )
    for b,bundle in enumerate(bundles):
        fileout.write(str(b) + "\n")
        for minigroup in bundle:

            #list lpgbts in minigroup:
            for lpgbt in minigroups_swap[minigroup]:
                fileout.write(str(lpgbt) + " ")
                
                #Get modules associated to each lpgbt:
                data_list = data[ ((data['TPGId1']==lpgbt) | (data['TPGId2']==lpgbt)) ]
                fileout.write(str(len(data_list)) + " ")
                for index, row in data_list.iterrows():
                    if ( row['density']==2 ):
                        fileout.write("1 " + str(row['layer']) + " " + str(row['u']) + " " + str(row['v']) + " " + str(row['TPGeLinkSum']) + " " )
                    else:
                        fileout.write("0 " + str(row['layer']) + " " + str(row['u']) + " " + str(row['v']) + " " + str(row['TPGeLinkSum']) + " " )
                fileout.write("\n")
                
    fileout.close()

#Code necessary for indentation:
#from https://stackoverflow.com/questions/13249415/how-to-implement-custom-indentation-when-pretty-printing-with-the-json-module
class NoIndent(object):
    """ Value wrapper. """
    def __init__(self, value):
        self.value = value

class MyEncoder(json.JSONEncoder):
    FORMAT_SPEC = '@@{}@@'
    regex = re.compile(FORMAT_SPEC.format(r'(\d+)'))

    def __init__(self, **kwargs):
        # Save copy of any keyword argument values needed for use here.
        self.__sort_keys = kwargs.get('sort_keys', None)
        super(MyEncoder, self).__init__(**kwargs)

    def default(self, obj):
        return (self.FORMAT_SPEC.format(id(obj)) if isinstance(obj, NoIndent)
                else super(MyEncoder, self).default(obj))

    def encode(self, obj):
        format_spec = self.FORMAT_SPEC  # Local var to expedite access.
        json_repr = super(MyEncoder, self).encode(obj)  # Default JSON.

        # Replace any marked-up object ids in the JSON repr with the
        # value returned from the json.dumps() of the corresponding
        # wrapped Python object.
        for match in self.regex.finditer(json_repr):
            # see https://stackoverflow.com/a/15012814/355230
            id = int(match.group(1))
            no_indent = PyObj_FromPtr(id)
            json_obj_repr = json.dumps(no_indent.value, sort_keys=self.__sort_keys)

            # Replace the matched id string with json formatted representation
            # of the corresponding Python object.
            json_repr = json_repr.replace(
                            '"{}"'.format(format_spec.format(id)), json_obj_repr)

        return json_repr
    
def produce_JsonMappingFile(Configuration):

    subconfig = Configuration['produce_JsonMappingFile']
    allocation = subconfig['allocation']
    output_name = subconfig['output_name']
    disconnected_modules = None
    if 'disconnected_modules' in subconfig.keys():
        disconnected_modules = subconfig['disconnected_modules']
    
    #Load allocation information
    info = loadConfiguration(allocation)
    data = info['data']
    mapping = info['mapping']
    minigroup_type = info['minigroup_type']
    nBundles = info['nBundles']
    maxInputs = info['maxInputs']
    configuration = info['configuration']
    random_seed = info['random_seed']
    nCallsToMappingMax = info['nCallsToMappingMax']
    max_modules = info['max_modules']
    cmsswNtuple = info['cmsswNtuple']
    max_towers_phi_region = max(info['max_towers'])
    git = info['git']

    #Get minigroups
    minigroups,minigroups_swap = getMinilpGBTGroups(data, minigroup_type)
    
    #Bundle together minigroup configuration
    bundles = getBundles(minigroups_swap,mapping,nBundles,maxInputs)

    #Open output file
    json_main = {}

    stage2list = []
    stage1linkslist = []
    stage1list = []
    #intialise empty list with number of minigroups
    lpgbtlist = [None]*len(minigroups)
    modulelist = []
    config_output = {}
    
    #1a) Stage 1 links to Stage 2 mapping (still preliminary)
    #Assume for now that the Stage 2 FPGA is attached to
    #two links from the current sector and one from the next sector

    nStage2Boards = 1
    nStage1Boards = len(bundles)

    for two in range(nStage2Boards):
        stage2dict = {}
        stage2_stage1links_list = []

        for one in range(nStage1Boards):
            link_dict = {}
            link_dict['SameSector'] = True
            stage2_stage1links_list.append(NoIndent(link_dict))
            link_dict = {}
            link_dict['SameSector'] = True
            stage2_stage1links_list.append(NoIndent(link_dict))
            link_dict = {}
            link_dict['SameSector'] = False
            stage2_stage1links_list.append(NoIndent(link_dict))
            
        stage2dict['Stage1Links'] = stage2_stage1links_list

        stage2list.append(stage2dict)

    #1b) Stage 1 FPGAs to Stage 1 links (still preliminary)
    #Assume for now each that each Stage 1 FPGA is connected to two links 
    #for the current sector and one link, which will go to the previous sector
    
    for stage1 in range(nStage1Boards):

        for i in range(2):
            stage1linkdict = {}
            stage1linkdict['Stage1'] = stage1
            stage1linkdict['Stage2SameSector'] = True
            stage1linkslist.append(NoIndent(stage1linkdict))
        stage1linkdict = {}
        stage1linkdict['Stage1'] = stage1
        stage1linkdict['Stage2SameSector'] = False
        stage1linkslist.append(NoIndent(stage1linkdict))
    
    #2) LpGBT mapping to Stage 1 and modules
    for b,bundle in enumerate(bundles):
       
        stage1dict = {}
        stage1dict["Stage1Links"] = {}
        stage1dict["lpgbts"] = []

        stage2_stage1links_list = [b*3,b*3+1,b*3+2]
        stage1dict["Stage1Links"] = stage2_stage1links_list
        
        for minigroup in bundle:

            #list lpgbts in minigroup:
            for lpgbt in minigroups_swap[minigroup]:
                stage1dict["lpgbts"].append(lpgbt)

                lpgbtdict = {}
                lpgbtdict['Stage1'] = b
                lpgbtdict['Modules'] = []
                
                #Get modules associated to each lpgbt:
                data_list = data[ ((data['TPGId1']==lpgbt) | (data['TPGId2']==lpgbt)) ]

                for index, row in data_list.iterrows():
                    lpgbt_moddict = {}
                    if ( row['density']==2 ):
                        lpgbt_moddict['isSilicon'] = False
                    else:
                        lpgbt_moddict['isSilicon'] = True
                    lpgbt_moddict['u'] = row['u']
                    lpgbt_moddict['v'] = row['v']
                    lpgbt_moddict['layer'] = row['layer']

                    lpgbtdict['Modules'].append(NoIndent(lpgbt_moddict))

                lpgbtlist[lpgbt] = lpgbtdict
                
        stage1list.append(NoIndent(stage1dict))

    #3) Get the module mapping information directly from the input mapping file
    for index, row in data.iterrows():

        module_lpgbtlist = []
        moduledict = {}

        if ( row['density']==2 ):
            moduledict['isSilicon'] = False
        else:
            moduledict['isSilicon'] = True
        moduledict['u'] = row['u']
        moduledict['v'] = row['v']
        moduledict['layer'] = row['layer']

        lpgbt_list = []
        for lpgbt in range(row['nTPG']):
            lpgbt_dict = {}
            if lpgbt == 0:
                lpgbt_dict['id'] = row['TPGId1']
                lpgbt_dict['nElinks'] = row['nTPGeLinks1']
            elif lpgbt == 1:
                lpgbt_dict['id'] = row['TPGId2']
                lpgbt_dict['nElinks'] = row['nTPGeLinks2']
            else:
                print ("Number of lpGBTs is limited to two per module")
            lpgbt_list.append(lpgbt_dict)
            
        moduledict['lpgbts'] = lpgbt_list
        modulelist.append(NoIndent(moduledict))

    #4) Optionally add "disconnected" modules, i.e. modules that may exist in the latest geometry, but not in the input mapping file used to produce the json output
    #Assumes as input a ROOT tree, produced from the CMSSW geometry tester (HGCalTriggerGeomTesterV9Imp3)
    heOffset = 28
    if disconnected_modules != None:
        disconnected_file = ROOT.TFile.Open(disconnected_modules,"READ")
        disconnected_tree = disconnected_file.Get("hgcaltriggergeomtester/TreeModuleErrors")

        for entry,event in enumerate(disconnected_tree):
            moduledict = {}
            if ( event.subdet==1 or event.subdet==2 ):
                moduledict['isSilicon'] = True
                if event.subdet==1:
                    moduledict['layer'] = event.layer
                else:
                    moduledict['layer'] = event.layer + heOffset
            else:
                moduledict['isSilicon'] = False
                moduledict['layer'] = event.layer + heOffset
                
            moduledict['u'] = event.waferu
            moduledict['v'] = event.waferv
            moduledict['lpgbts'] = [] #empty, i.e. not connected
            modulelist.append(NoIndent(moduledict))


    #Save input configuration
    config_output["hgcal-linkmapping-configuration"] = configuration
    config_output["hgcal-linkmapping-randomseed"] = random_seed
    config_output["hgcal-linkmapping-ncalls"] = nCallsToMappingMax
    config_output["hgcal-linkmapping-maxmodules"] = max_modules
    config_output["hgcal-linkmapping-maxtowersphiregion"] = max_towers
    config_output["hgcal-linkmapping-cmsswntuple"] = cmsswNtuple
    config_output["hgcal-linkmapping-gitrevision"] = git
    
    json_main['Config'] = config_output            
    json_main['Stage2'] = stage2list
    json_main['Stage1Links'] = stage1linkslist
    json_main['Stage1'] = stage1list
    json_main['lpgbt'] = lpgbtlist
    json_main['Module'] = modulelist
    
    #Write to file
    with open(output_name, 'w') as fp:
        #data = json.dumps(json_main, indent=2, ensure_ascii=False)
        data = json.dumps(json_main, ensure_ascii=False, cls=MyEncoder, indent=4)
        fp.write(data)
    
def produce_nTCsPerModuleHists(Configuration):

    subconfig = Configuration['produce_nTCsPerModuleHists']
    allocation = subconfig['allocation']
    CMSSW_ModuleHists = subconfig['CMSSW_ModuleHists']
    
    #Load allocation information
    info = loadConfiguration(allocation)
    data = info['data']
    mapping = info['mapping']
    minigroup_type = info['minigroup_type']
    nBundles = info['nBundles']
    maxInputs = info['maxInputs']

    #Get minigroups
    minigroups,minigroups_swap = getMinilpGBTGroups(data, minigroup_type)

    #Get list of which modules are in each minigroup
    minigroups_modules = getMiniModuleGroups(data,minigroups_swap)
    
    #Bundle together minigroup configuration
    bundles = getBundles(minigroups_swap,mapping,nBundles,maxInputs)

    #Get nTC hists per module
    module_hists = getModuleTCHists(CMSSW_ModuleHists)

    #Open output file
    outfile = ROOT.TFile.Open("hists_per_bundle.root","RECREATE")
    for b,bundle in enumerate(bundles):
        outfile.mkdir("bundle_" + str(b))
        outfile.cd("bundle_" + str(b)) 
        for minigroup in bundle:

            for module in minigroups_modules[minigroup]:

                module_hists[tuple(module)].Write()

        outfile.cd()

def check_for_missing_modules(Configuration):

    subconfig = Configuration['check_for_missing_modules']
    
    if ( subconfig['inMappingFile'] ):
        print("Missing modules in mapping file: "+ subconfig['MappingFile'] + "\n")
        check_for_missing_modules_inMappingFile(subconfig)
    if ( subconfig['inCMSSW'] ):
        print("\nMissing modules in CMSSW\n")
        check_for_missing_modules_inCMSSW(subconfig)
            
def check_for_missing_modules_inMappingFile(subconfig):

    #Check for modules missing in the mapping file
    MappingFile = subconfig['MappingFile']
    CMSSW_Silicon = subconfig['CMSSW_Silicon']
    CMSSW_Scintillator = subconfig['CMSSW_Scintillator']
    
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

def check_for_missing_modules_inCMSSW(subconfig):

    MappingFile = subconfig['MappingFile']
    CMSSW_Silicon = subconfig['CMSSW_Silicon']
    CMSSW_Scintillator = subconfig['CMSSW_Scintillator']

    #Load external data
    data = loadDataFile(MappingFile) #dataframe    
    data_tcs_passing,data_tcs_passing_scin = getTCsPassing(CMSSW_Silicon,CMSSW_Scintillator) #from CMSSW
    getHexModuleLoadInfo(data,data_tcs_passing,data_tcs_passing_scin,True)
    
def study_mapping(Configuration):

    #Get Git revision
    try:
        git = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode("utf-8")
    except:
        git = "no git repository detected"
        
    #Global settings
    output_dir = Configuration['output_dir']
    print_level = Configuration['print_level']
    
    #Function specific settings
    subconfig = Configuration['study_mapping']
    MappingFile = subconfig['MappingFile']
    CMSSW_ModuleHists = subconfig['CMSSW_ModuleHists']
    algorithm = subconfig['algorithm']
    initial_state = subconfig['initial_state']
    max_iterations = subconfig['max_iterations']
    max_calls = subconfig['max_calls']
    minigroup_type = subconfig['minigroup_type']
    TowerMappingFile = subconfig['TowerMappingFile']
    TowerPhiSplit = subconfig['TowerPhiSplit']

    random_seed = subconfig['random_seed']
    if random_seed == None:
        random_seed = random.randrange(2**32-1)
    
    fpgaConfig = None
    correctionConfig = None
    chi2Config = None
    phisplitConfig = None
    cmsswNtuple = ""
    
    if 'fpgas' in subconfig.keys():
        fpgaConfig = subconfig['fpgas']
    if 'corrections' in subconfig.keys():
        correctionConfig = subconfig['corrections']
    if 'chi2' in subconfig.keys():
        chi2Config = subconfig['chi2']
    if 'phisplit' in subconfig.keys():
        phisplitConfig = subconfig['phisplit']            

    #Load external data
    data = loadDataFile(MappingFile) #dataframe

    #Load FPGA Information
    if ( fpgaConfig != None ):
        nBundles = fpgaConfig["nBundles"]
        maxInputs = fpgaConfig["maxInputs"]
    else:
        #Set defaults
        nBundles = 14
        maxInputs = 120
    
    try:

        #Configuration for how to divide TCs into phidivisionX and phidivisionY (traditionally phi > 60 and phi < 60)
        split = "per_roverz_bin"
        phidivisionX_fixvalue_min = 55
        phidivisionY_fixvalue_max = None
        
        if phisplitConfig != None:
            split = phisplitConfig['type']
            if 'phidivisionX_fixvalue_min' in phisplitConfig.keys():
                phidivisionX_fixvalue_min = phisplitConfig['phidivisionX_fixvalue_min']
            if 'phidivisionY_fixvalue_max' in phisplitConfig.keys():
                phidivisionY_fixvalue_max = phisplitConfig['phidivisionY_fixvalue_max']

        inclusive_hists,module_hists = getModuleHists(CMSSW_ModuleHists, split = split, phidivisionX_fixvalue_min = phidivisionX_fixvalue_min, phidivisionY_fixvalue_max = phidivisionY_fixvalue_max)
        cmsswNtuple = getCMSSWNtupleName(CMSSW_ModuleHists)

    except EnvironmentError:
        print ( "File " + CMSSW_ModuleHists + " does not exist" )
        exit()
    # Apply various corrections to r/z distributions from CMSSW

    if correctionConfig != None:
        print ( "Applying geometry corrections" )
        applyGeometryCorrections( inclusive_hists, module_hists, correctionConfig )

    include_errors_in_chi2 = False
    include_max_modules_in_chi2 = False
    include_max_towers_in_chi2 = False
    max_modules_weighting_factor = 1000
    max_towers_weighting_factor = 30000
    max_towers_weighting_option = 2
    max_towers_step_point = 180
    weight_bins_proportionally = True
    if chi2Config != None:
        if 'include_errors_in_chi2' in chi2Config.keys():
            include_errors_in_chi2 = chi2Config['include_errors_in_chi2']
        if 'include_max_modules_in_chi2' in chi2Config.keys():
            include_max_modules_in_chi2 = chi2Config['include_max_modules_in_chi2']
        if 'max_modules_weighting_factor' in chi2Config.keys():
            max_modules_weighting_factor = chi2Config['max_modules_weighting_factor']
        if 'include_max_towers_in_chi2' in chi2Config.keys():
            include_max_towers_in_chi2 = chi2Config['include_max_towers_in_chi2']
        if 'max_modules_weighting_factor' in chi2Config.keys():
            max_towers_weighting_factor = chi2Config['max_towers_weighting_factor']
        if 'max_towers_weighting_option' in chi2Config.keys():
            max_towers_weighting_option = chi2Config['max_towers_weighting_option']
        if 'max_towers_step_point' in chi2Config.keys():
            max_towers_step_point = chi2Config['max_towers_step_point']
        if 'weight_bins_proportionally' in chi2Config.keys():
            weight_bins_proportionally = chi2Config['weight_bins_proportionally']
            
    #Load tower data if required
    if include_max_towers_in_chi2:
        try:
            towerdata = loadModuleTowerMappingFile(TowerMappingFile)
        except EnvironmentError:
            print ( "File " + TowerMappingFile + " does not exist" )
            exit()
            
    #Form hists corresponding to each lpGBT from module hists
    lpgbt_hists = getlpGBTHists(data, module_hists)

    minigroups,minigroups_swap = getMinilpGBTGroups(data, minigroup_type)
    minigroup_hists = getMiniGroupHists(lpgbt_hists,minigroups_swap,return_error_squares=include_errors_in_chi2)
    minigroup_hists_root = getMiniGroupHists(lpgbt_hists,minigroups_swap,root=True)
    #Get list of which modules are in each minigroup
    minigroups_modules = getMiniModuleGroups(data,minigroups_swap)

    #Get list of which towers are in each minigroup
    if include_max_towers_in_chi2:
        minigroups_towers = getMiniTowerGroups(towerdata, minigroups_modules)

    
    def mapping_max(state):
        global chi2_min
        global combbest
        global nCallsToMappingMax
    
        max_modules = None
        max_towers = None
        chi2 = 0
    
        bundles = getBundles(minigroups_swap,state,nBundles,maxInputs)
        bundled_lpgbthists = getBundledlpgbtHists(minigroup_hists,bundles)

        if include_max_modules_in_chi2:
            max_modules = getMaximumNumberOfModulesInABundle(minigroups_modules,bundles)
        if include_max_towers_in_chi2:
            max_towers_list = getMaxTowersList(minigroups_towers, bundles, TowerPhiSplit)
            max_towers = max(max_towers_list)

        chi2 = calculateChiSquared(inclusive_hists,bundled_lpgbthists,nBundles,max_modules,max_modules_weighting_factor,max_towers,[max_towers_weighting_factor,max_towers_weighting_option,max_towers_step_point], weight_bins_proportionally)

        typicalchi2 = 600000000000
        if include_errors_in_chi2:
            typicalchi2 = 10000000
        if (chi2<chi2_min):
            chi2_min = chi2
            combbest = np.copy(state)
            if ( print_level > 0 ):
                print (algorithm," ", chi2_min, " ", chi2_min/typicalchi2)
                if include_max_towers_in_chi2:
                    print ("max_towers = ", max_towers)
                if include_max_modules_in_chi2:
                    print ("max_modules = ", max_modules)
            if ( print_level > 1 ):
                print (repr(combbest))

        nCallsToMappingMax += 1
        if max_calls != None:
            if max_calls == nCallsToMappingMax:
                raise exitProgramSignal
        
        return chi2

    init_state = []
    if (initial_state[-4:] == ".npy"):
        print (initial_state)
        previousConfig = loadConfiguration(initial_state)
        init_state = previousConfig['mapping']
        #Get the previous random seed, and number of calls
        #to restart where previous run finished
        #This does assume the rest of the config is the same
        previousnCalls = previousConfig['nCallsToMappingMax']
        previousRandomSeed = previousConfig['random_seed']    
        global nCallsToMappingMax
        nCallsToMappingMax += previousnCalls
        np.random.seed(previousRandomSeed)
        
        if ( len(init_state) != len(minigroups_swap) ):
            print ( "Initial state should be the same length as the number of mini groups")
            exit()
    elif (initial_state == "random"):
        np.random.seed(random_seed)
        init_state = np.arange(len(minigroups_swap))
        np.random.shuffle(init_state)

    
    fitness_cust = mlrose.CustomFitness(mapping_max)
    # Define optimization problem object
    problem_cust = mlrose.DiscreteOpt(length = len(init_state), fitness_fn = fitness_cust, maximize = False, max_val = len(minigroups_swap), minigroups = minigroups_swap, nBundles = nBundles)

    # Define decay schedule
    decay_schedule = "ExponentialDecay"
    if decay_schedule == "ExponentialDecay":
        schedule = mlrose.ExpDecay()
    elif decay_schedule == "ArithmeticDecay":
        schedule = mlrose.ArithDecay()
    else:
        print ("Unknown decay schedule")
        exit()
        
    filename = "bundles_job_"
    filenumber = ""
    if ( len(sys.argv) > 2 ):
        filenumber = str(sys.argv[2])
    else:
        filenumber = "default"
    filename+=filenumber
    
    if ( algorithm == "save_root" ):
        #Save best combination so far into a root file
        bundles = getBundles(minigroups_swap,init_state,nBundles,maxInputs)

        bundled_hists_root = getBundledlpgbtHistsRoot(minigroup_hists_root,bundles)
        bundled_hists = getBundledlpgbtHists(minigroup_hists,bundles)

        chi2 = calculateChiSquared(inclusive_hists,bundled_hists,nBundles,max_modules,max_modules_weighting_factor,max_towers,[max_towers_weighting_factor,max_towers_weighting_option,max_towers_step_point], weight_bins_proportionally)
        newfile = ROOT.TFile("bundles_roverz.root","RECREATE")
        for sector in bundled_hists_root:
            for key, value in sector.items():
                value.Write()
        for sector in inclusive_hists:
            sector.Scale(1./float(nBundles))
            sector.Write()
        newfile.Close()
        print ("Chi2:",chi2)
        print ("List of Bundles:")
        for b,bundle in enumerate(bundles):
            print ("" )
            print ("bundle" + str(b) )
            for minigroup in bundle:
                lpgbts = minigroups_swap[minigroup]
                for lpgbt in lpgbts:
                    print (str(lpgbt) + ", "  , end = '')

    elif algorithm == "random_hill_climb" or algorithm == "simulated_annealing":

        try:
            if (algorithm == "random_hill_climb"):
                best_state, best_fitness = mlrose.random_hill_climb(problem_cust, max_attempts=10000, max_iters=max_iterations, restarts=0, init_state=init_state, random_state=random_seed)
            elif (algorithm == "simulated_annealing"):
                best_state, best_fitness = mlrose.simulated_annealing(problem_cust, schedule = schedule, max_attempts = 100000, max_iters = 10000000, init_state = init_state, random_state=random_seed)
                

        except exitProgramSignal:
            print("interrupt received, stopping and saving")

        finally:
            signal.signal(signal.SIGUSR1,dummy_handler) # avoid any interrupt when finalising
            bundles = getBundles(minigroups_swap,combbest,nBundles,maxInputs)
            if include_max_modules_in_chi2:
                max_modules = getMaximumNumberOfModulesInABundle(minigroups_modules,bundles)
            else:
                max_modules = 'Not used in chi2'
            if include_max_towers_in_chi2:
                max_towers_list = getMaxTowersList(minigroups_towers, bundles, TowerPhiSplit)
            else:
                max_towers_list = 'Not used in chi2'
            with open( output_dir + "/" + filename + ".npy", "wb") as filep:
                result_config_git = [bundles,subconfig,random_seed,nCallsToMappingMax,max_modules,max_towers_list,cmsswNtuple,git]
                pickle.dump(result_config_git, filep)
            file1 = open(output_dir + "/chi2_"+filenumber+".txt","a")
            file1.write( "bundles[" + filenumber + "] = " + str(chi2_min) + "\n" )
            file1.close( )

    else:
        print("Algorithm "+ algorithm + " currently not implemented" )

    
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

    ROOT.TH1.SetDefaultSumw2()

    function_choice = config['function']
    
    if ( function_choice['study_mapping'] ):                    
        study_mapping(config)

    if ( function_choice['check_for_missing_modules'] ):
        check_for_missing_modules(config)

    if ( function_choice['plot_lpGBTLoads'] ):
        plot_lpGBTLoads(config)

    if ( function_choice['plot_ModuleLoads'] ):
        plot_ModuleLoads(config)

    if ( function_choice['produce_AllocationFile'] ):
        produce_AllocationFile(config)

    if ( function_choice['produce_nTCsPerModuleHists'] ):
        produce_nTCsPerModuleHists(config)

    if ( function_choice['produce_JsonMappingFile'] ):
        produce_JsonMappingFile(config)
    
main()
