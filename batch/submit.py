#!/usr/bin/env python
import os, re, subprocess
import math, time
import sys
import itertools
from datetime import datetime

today = datetime.today().strftime('%Y%m%d')
OutputDir = '/vols/cms/snwebb/hgcal/analysis/outputs/' + today + "-1"
NumberOfJobs = 100
UniqueSeed = 0

def getUniqueSeed():
   global UniqueSeed
   UniqueSeed = UniqueSeed+1
   return UniqueSeed

def createAndSubmitJob(jobnumber,modules,towers):

   randomseed = getUniqueSeed()
   #jobnumber+1000

   ##### create config file #######
   configfile = open('/vols/cms/snwebb/hgcal/analysis/config_template.yaml', 'r')
   string_list = configfile.readlines()
   for l,line in enumerate(string_list):
      if line.find("RANDOMSEED_TEMPLATE")!=-1:
         string_list[l] = line.replace("RANDOMSEED_TEMPLATE",str(randomseed))
      if line.find("MAXMODWEIGHT")!=-1:
         string_list[l] = line.replace("MAXMODWEIGHT",str(modules))
      if line.find("MAXTOWERWEIGHT")!=-1:
         string_list[l] = line.replace("MAXTOWERWEIGHT",str(towers))

   #configfile = open('config_'+str(randomseed)+'.yaml', 'w')
   jobname = 'module'+str(modules)+ '_tower' + str(towers) + '_' + str(jobnumber) 
   configfile = open('config_' + jobname + '.yaml', 'w')
   configfile.writelines(string_list)

   ##### creates jobs #######
   with open('job_'+ str(jobname)+'.sh', 'w') as fout:                  

      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo\n")
      fout.write("ulimit -c 0\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'WORKDIR ' ${PWD}\n")
      fout.write("echo 'TEMPDIR ' ${TMPDIR}\n")

      fout.write("trap \"echo SIGINT seen\"  SIGINT\n")
      fout.write("trap \"echo SIGUSR1 seen\" SIGUSR1\n")
      fout.write("trap \"echo SIGUSR2 seen\" SIGUSR2\n")
      fout.write("trap \"echo SIGTERM seen\" SIGTERM\n")
      fout.write("trap \"echo SIGXCPU seen\" SIGXCPU\n")

      fout.write("cd /vols/cms/snwebb/hgcal/analysis/\n")
      fout.write("source start_mapping_env.sh \n")

      fout.write("cd $TMPDIR \n")

      #Copy relevant files
      #fout.write("mkdir data \n")
      #fout.write("cp /vols/cms/snwebb/hgcal/analysis/hgcal-linkmapping/*py . \n")
      fout.write("cp -r /vols/cms/snwebb/hgcal/analysis/hgcal-linkmapping . \n")
      fout.write("cd hgcal-linkmapping \n")
      # fout.write("cp /vols/cms/snwebb/hgcal/analysis/hgcal-linkmapping/data/ROverZ*root data/. \n")
      # fout.write("cp -r /vols/cms/snwebb/hgcal/analysis/hgcal-linkmapping/externals . \n")
      # fout.write("cp /vols/cms/snwebb/hgcal/analysis/hgcal-linkmapping/data/*txt data/. \n")
      fout.write("cp " + OutputDir + "/tmp/config_" + jobname + ".yaml . \n")

      #fout.write("cd hgcal-linkmapping\n")
      #fout.write("./main.py " + OutputDir + "/tmp/config_" + jobname + ".yaml " + str(randomseed)+" \n")
      fout.write("./main.py config_" + jobname + ".yaml " + str(randomseed)+" \n")

      fout.write("echo 'copying output:'\n")
      fout.write("cp *npy " + OutputDir + "/. \n")
      fout.write("cp *txt " + OutputDir + "/. \n")
 
      fout.write("echo 'WAIT 120 seconds'\n")
      fout.write("sleep 120\n")
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")

   os.system("chmod 755 job_"+str(jobname)+".sh")

   ###### sends jobs ######
   os.system("qsub -cwd -q hep.q -l h_vmem=12G -l s_vmem=11.8G -l h_rt=3:5:0 -l s_rt=2:55:0 job_"+str(jobname)+".sh")

   print ("job nr " + str(randomseed) + " submitted")


def main():
   print ()
   print ('START')
   print ()

   os.system("mkdir -p " + OutputDir + "/tmp")
   os.chdir(OutputDir + "/tmp/")

   ##### loop for creating and sending jobs #####
   
#   for modules,towers in itertools.product([0,50000,100000,250000,1000000],repeat=2):
   for modules,towers in itertools.product([30000],repeat=2):
      
      for x in range(NumberOfJobs):
         ##### creates directory and file list for job #######
         createAndSubmitJob(x,modules,towers)


   #TESTING
   # for modules,towers in itertools.product([0,50000,100000,250000,1000000],repeat=2):
   #    for x in range(NumberOfJobs):
   #       print (modules,towers,x,getUniqueSeed())

   print ()
   print ("your jobs:")
   os.system("qstat")
   print ()
   print ('END')
   print ()
 

if __name__ == '__main__':
  main()
