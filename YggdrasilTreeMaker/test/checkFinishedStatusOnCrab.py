# /usr/bin/python

#Author: Ben Tannenwald
#Date: March 13, 2019
#Purpose: Script to check % finished status for CRAB3 jobs

import os,sys,subprocess


def getOutputLine( _statusWord, _inputLines, _dataset ):

    _passesQuickCheck    = [_statusWord in x for x in _inputLines].count(True)
    _runningWord = 'running'
    _stillHasRunningJobs = [_runningWord in x for x in _inputLines].count(True)
    _stillHasTransferringJobs = ['transferring' in x for x in _inputLines].count(True)
    _stillHasIdleJobs = ['idle' in x for x in _inputLines].count(True)

    _statusOutput = ''


    if _passesQuickCheck:
        _statusLine = [ x for x in [line.replace('\t','').split(' ') for line in _inputLines if _statusWord in line][0] if x !='' ]
        
        if len(_statusLine) == 5: # 100% finished!!
            _statusOutput = _dataset.lstrip(datasetTag) + '\t\t'+_statusLine[2]+'\t'+_statusLine[3]+'\t'+_statusLine[4]
        if len(_statusLine) == 6: # only just started to complete jobs
            if _stillHasRunningJobs or _stillHasTransferringJobs or _stillHasIdleJobs:
                _statusOutput = _dataset.lstrip(datasetTag) + '\t\t'+_runningWord+'\t'+_statusLine[3]+'\t('+_statusLine[5]
            else:
                _statusOutput = _dataset.lstrip(datasetTag) + '\t\t'+_statusWord+'\t'+_statusLine[3]+'\t('+_statusLine[5]

        elif len(_statusLine) == 3: # not 100% finished
            if _stillHasRunningJobs or _stillHasTransferringJobs or _stillHasIdleJobs:
                _statusOutput = _dataset.lstrip(datasetTag)+'\t\t'+_runningWord+'\t'+_statusLine[1]+'\t'+_statusLine[2]
            else:
                _statusOutput = _dataset.lstrip(datasetTag)+'\t\t'+_statusWord+'\t'+_statusLine[1]+'\t'+_statusLine[2]

        elif len(_statusLine) == 0:
            _statusOutput = 'CONFUSED:\t\t'+_dataset.lstrip(datasetTag)+'\t\t'+_statusLine
    else:
        _statusOutput = _dataset.lstrip(datasetTag)+'\t\t\t\t NO FINISHED INFO'

    return _statusOutput

######################################################################################################################
######################################################################################################################


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
iDataset = 0
for iDataset in crabSamplesToCheck:
    raw_statusOutputLines = (subprocess.check_output("crab status -d {0}".format( iDataset ), shell=True)).split('\n')

    finishedStatus = getOutputLine( 'finished', raw_statusOutputLines, iDataset )
    print ( finishedStatus )
    outFile.write( finishedStatus+'\n' )

    #runningStatus  = getOutputLine( 'running', raw_statusOutputLines, iDataset )
    #print ( runningStatus )
    #outFile.write( runningStatus+'\n' )
       
    #failedStatus   = getOutputLine( 'failed', raw_statusOutputLines, iDataset )
    #print ( failedStatus )
    #outFile.write( failedStatus+'\n' )
