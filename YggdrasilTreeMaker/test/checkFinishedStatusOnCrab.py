# /usr/bin/python

#Author: Ben Tannenwald
#Date: March 13, 2019
#Purpose: Script to check % finished status for CRAB3 jobs

import os,sys,subprocess

datasetTag = 'crab_preApproval2017_v2_'
crabSamplesToCheck = d= [y for y in [x[1] for x in os.walk('./')][0] if datasetTag in y] # auto-generate samples

# *** 0. Preliminaries

# ** A. Check if grid proxy exists and exit if not
gridProxyExists = True if os.path.isfile('/tmp/x509up_u40094') else False
if gridProxyExists is False:
    print('No grid proxy detected. Setup voms first plz. EXITING\n')
    quit()

# ** B. Get date+time and make output file name
raw_date = (subprocess.check_output("date", shell=True)).rstrip('\n').split(' ')
outName = 'crabFinishedStatus_' + raw_date[1] + '-' + raw_date[2] + '-' + raw_date[5] + '_' + raw_date[3]+ '.txt'
outFile = open( outName, 'w')


# *** 1. Get output fom crab status -d
#for sample in range(0
iDataset = 0
for iDataset in crabSamplesToCheck:
    raw_statusOutputLines = (subprocess.check_output("crab status -d {0}".format( iDataset ), shell=True)).split('\n')
    passesQuickCheck = ['finished' in x for x in raw_statusOutputLines].count(True)

    if passesQuickCheck:
        finishedLine = [ x for x in [line.replace('\t','').split(' ') for line in raw_statusOutputLines if 'finished' in line][0] if x !='' ]
        if len(finishedLine) == 5: # 100% finished!!
            finishedStatus = iDataset.lstrip(datasetTag) + '\t\t'+finishedLine[2]+'\t'+finishedLine[3]+'\t'+finishedLine[4]
        elif len(finishedLine) == 3: # not 100% finished
            finishedStatus = iDataset.lstrip(datasetTag)+'\t\t'+finishedLine[0]+'\t'+finishedLine[1]+'\t'+finishedLine[2]
        elif len(finishedLine) == 0:
            finishedStatus = 'CONFUSED:\t\t'+iDataset.lstrip(datasetTag)+'\t\t'+finishedLine
    else:
        finishedStatus = iDataset.lstrip(datasetTag)+'\t\t\t\t RUNNING'

    print ( finishedStatus )
    outFile.write( finishedStatus+'\n' )
