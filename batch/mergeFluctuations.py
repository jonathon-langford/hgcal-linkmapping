#!/usr/bin/python3
import pickle
import numpy as np

start = np.arange(0,9000,50)
stop = np.arange(49,9000,50)

bundled_lpgbthists_allevents = []

for sta,sto in zip(start,stop):
    
    #filename = "outputFluctuationsNG-phi3060/alldata_from"+str(sta)+"_to"+str(sto)+".txt"
    #filename = "outputFluctuationsTT-phi3060/alldata_from"+str(sta)+"_to"+str(sto)+".txt"
    #filename = "outputFluctuationsTT-210505/alldata_from"+str(sta)+"_to"+str(sto)+"_sumpt.txt"
    #filename = "outputFluctuationsTT-210511-14bins/alldata_from"+str(sta)+"_to"+str(sto)+"_sumpt.txt"
    filename = "outputFluctuationsTT-210617/alldata_from"+str(sta)+"_to"+str(sto)+".txt"
    
    with open(filename, "rb") as filep:   
        bundled_lpgbthists_partial = pickle.load(filep)
        bundled_lpgbthists_allevents.extend(bundled_lpgbthists_partial)
        

#with open( "alldata_NG-phi3060_merged.txt", "wb") as filep:
#with open( "alldata_TT.txt", "wb") as filep:
#with open( "alldata_TT-phi3060_merged.txt", "wb") as filep:
#with open( "alldata_ttbar_200911_sumpt.txt", "wb") as filep:
#with open( "alldata_ttbar_200911_phi3090_sumpt.txt", "wb") as filep:
#with open( "alldata_ttbar_210505_pt.txt", "wb") as filep:
#with open( "alldata_ttbar_210511_14bins_pt.txt", "wb") as filep:
with open( "alldata_ttbar_210617.txt", "wb") as filep:
    pickle.dump(bundled_lpgbthists_allevents, filep)


