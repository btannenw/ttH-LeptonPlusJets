# Instructions for running

## for csh shell
source /cvmfs/cms.cern.ch/cmsset_default.csh

setenv SCRAM_ARCH slc6_amd64_gcc630

## for sh shell
export SCRAM_ARCH="slc6_amd64_gcc630"

## common
cmsrel CMSSW_9_4_9

cd CMSSW_9_4_9/src

cmsenv

git cms-init

git clone https://github.com/btannenw/MiniAOD.git -b CMSSW_940 

git clone https://github.com/hsatoshi/GenParticleTopOriginChargedleptonFilter.git ttHAnalysisSubprogram/GenParticleTopOriginChargedleptonFilter -b  master 

git clone https://github.com/hsatoshi/PuppiLeptonIsolationhelper.git -b  master  

git clone https://github.com/btannenw/ttH-LeptonPlusJets.git  -b 94x

git cms-merge-topic cms-egamma:EgammaPostRecoTools_940

git cms-merge-topic guitargeek:EgammaID_9_4_X 

git cms-merge-topic cms-met:METFixEE2017_949_v2

git cms-merge-topic yrath:deterministicSeeds 

scram b -j 8

## Correcting github origins
cd ttH-LeptonPlusJets
git remote set-url origin git@github.com:btannenw/ttH-LeptonPlusJets.git
cd ../MiniAOD
git remote set-url origin git@github.com:btannenw/MiniAOD.git

# Adding BDT/DNN/MEM code

git clone https://gitlab.cern.ch/ttH/CommonClassifier.git TTH/CommonClassifier

source TTH/CommonClassifier/setup/install_mem.sh

source TTH/CommonClassifier/setup/install_recoLikelihood.sh

scram b -j 8

## Crab submission
Adding TTH/CommonClassifier results in a crab tarball that's too large for standard submission. Is there a smart way around this? Probably. In the meantime, you can do the following to clean out some of the older BDT/MEM weight collections that should be unnecessary

> cd $CMSSW_BASE

> source src/ttH-LeptonPlusJets/cleanup.sh 

