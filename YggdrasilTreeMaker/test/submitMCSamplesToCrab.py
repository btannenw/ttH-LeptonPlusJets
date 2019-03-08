# /usr/bin/python

#Author: Ben Tannenwald
#Date: March 8, 2019
#Purpose: Script to submit CRAB3 jobs for Yggdrasil production over subsets of necessary MC samples

import os,sys, argparse

# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--sampleSet", help="subset of MC samples user specifies to submit")
args = parser.parse_args()

if (len(vars(args)) != 1): # 1 --> three for options
    os.system('python submitMCSamplesToCrab.py -h')
    quit()

# ** A. Test sampleSet value and exit if bad
if args.sampleSet is None:
    print ("#### Need to specify sample subset file using --sampleSet <1/2/3/4> ####\nEXITING")
    quit()
else:
    if (args.sampleSet).isdigit():
        sampleSet = int(args.sampleSet)
        if sampleSet > 0 and sampleSet <=4:
            print ("#### Setting sampleSet = {0} ####\n".format(sampleSet))
        else:
            print ("#### Passed sampleSet = {0} is not in proper range <1/2/3/4>. Try again. EXITING ####\n".format(sampleSet))
            quit()
    else:
        print ("#### Passed sampleSet={0} is not a number. Need to specify sample subset file using --sampleSet <1/2/3/4> ####\nEXITING".format(args.sampleSet))
        quit()

# ** B. Exit if no grid proxy
#if ( not os.path.exists(os.path.expandvars("$X509_USER_PROXY")) ):
#    print "#### No GRID PROXY detected. Please do voms-proxy-init -voms cms before submitting CRAB jobs ####.\nEXITING"
#    quit()


# *** 1. Do some stuff
