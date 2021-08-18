#!/usr/bin/env python
import os, re, subprocess
import math, time
import sys
import numpy as np

#OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputFluctuationsNG-phi3060/'
#OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputFluctuationsTT-phi3060/'
#OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputFluctuationsTT-1812-7545/'
#OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputFluctuationsTT-210511-14bins/'
#OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputFluctuationsTT-210602/'
OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputFluctuationsTT-210617/'

start = np.arange(0,9000,50)
stop = np.arange(49,9000,50)

#NumberOfJobs = 200

def main():
   print ()
   print ('START')
   print ()

   os.system("mkdir -p " + OutputDir + "/tmp")
   os.chdir(OutputDir + "/tmp/")

   ##### loop for creating and sending jobs #####
   #for x in range(NumberOfJobs):
   for sta,sto in zip(start,stop):
   ##### creates directory and file list for job #######

      ##### create config file #######
      configfile = open('/vols/cms/snwebb/hgcal/analysis/fluctuations_config_template.yaml', 'r')
      string_list = configfile.readlines()
      for l,line in enumerate(string_list):
         if line.find("START")!=-1:
            string_list[l] = line.replace("START",str(sta))
         if line.find("FINISH")!=-1:
            string_list[l] = line.replace("FINISH",str(sto))
         if line.find("OUTNAME")!=-1:
            string_list[l] = line.replace("OUTNAME",OutputDir+"alldata")

      configfile = open('config_'+str(sta)+"_"+str(sto)+'.yaml', 'w')
      configfile.writelines(string_list)

      ##### creates jobs #######
      with open('job_'+str(sta)+"_"+str(sto)+'.sh', 'w') as fout:                  

         fout.write("#!/bin/sh\n")
         fout.write("echo\n")
         fout.write("echo\n")
         fout.write("ulimit -c 0\n")
         fout.write("echo 'START---------------'\n")
         fout.write("echo 'WORKDIR ' ${PWD}\n")

         fout.write("trap \"echo SIGINT seen\"  SIGINT\n")
         fout.write("trap \"echo SIGUSR1 seen\" SIGUSR1\n")
         fout.write("trap \"echo SIGUSR2 seen\" SIGUSR2\n")
         fout.write("trap \"echo SIGTERM seen\" SIGTERM\n")
         fout.write("trap \"echo SIGXCPU seen\" SIGXCPU\n")

         fout.write("cd /vols/cms/snwebb/hgcal/analysis/\n")
         fout.write("source start_mapping_env.sh \n")
         fout.write("cd hgcal-linkmapping\n")
         fout.write("./fluctuation.py " + OutputDir + "/tmp/config_"+str(sta)+"_"+str(sto)+".yaml \n")

         fout.write("echo 'STOP---------------'\n")
         fout.write("echo\n")
         fout.write("echo\n")

      os.system("chmod 755 job_"+str(sta)+"_"+str(sto)+".sh")

      ###### sends jobs ######
      os.system("qsub -cwd -q hep.q -l h_vmem=6G -l s_vmem=5.8G -l h_rt=1:00:0 -l s_rt=0:57:0 job_"+str(sta)+"_"+str(sto)+".sh")


      print ("job nr " + str(sta)+"_"+str(sto) + " submitted")

   print ()
   print ("your jobs:")
   os.system("qstat")
   print ()
   print ('END')
   print ()
 

if __name__ == '__main__':
  main()
