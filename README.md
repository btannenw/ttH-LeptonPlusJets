ttH-LeptonPlusJets
==================

### 74x
    source /cvmfs/cms.cern.ch/cmsset_default.csh
    setenv SCRAM_ARCH slc6_amd64_gcc491
    cmsrel CMSSW_7_4_4_patch4
    cd CMSSW_7_4_4_patch4/src/
    cmsenv

    git cms-addpkg PhysicsTools/JetMCAlgos
    cp /afs/cern.ch/user/l/lwming/public/74X/GenHFHadronMatcher.cc PhysicsTools/JetMCAlgos/plugins/GenHFHadronMatcher.cc
    cp /afs/cern.ch/user/l/lwming/public/74X/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/GenTtbarCategorizer.cc
    cp /afs/cern.ch/user/l/lwming/public/74X/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/GenTtbarCategorizer_cfi.py

    git cms-merge-topic gkasieczka:htt-v2-74X
    git clone https://github.com/cms-ttH/BoostedTTH.git
    ln -s $CMSSW_RELEASE_BASE/src/RecoJets/JetProducers/plugins/VirtualJetProducer.h BoostedTTH/BoostedProducer/plugins/VirtualJetProducer.h
    ln -s $CMSSW_RELEASE_BASE/src/RecoJets/JetProducers/plugins/VirtualJetProducer.cc BoostedTTH/BoostedProducer/plugins/VirtualJetProducer.cc

    git clone https://github.com/cms-ttH/MiniAOD.git

\#   git clone https://github.com/cms-ttH/ttH-LeptonPlusJets.git
    git clone https://github.com/hsatoshi/ttH-LeptonPlusJets.git
    cd ttH-LeptonPlusJets
    git fetch 
    git checkout -b satoshi origin/satoshi
    cd ../

    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git checkout 74x-root6
    cd -
    git clone https://github.com/cms-ttH/ttH-Limits.git ttH/Limits

    scram b -j 32


### 72X
    source /cvmfs/cms.cern.ch/cmsset_default.csh
    setenv SCRAM_ARCH slc6_amd64_gcc481
    cmsrel CMSSW_7_2_3
    cd CMSSW_7_2_3/src/
    cmsenv

    git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720

    git cms-addpkg PhysicsTools/JetMCAlgos
    cp /afs/cern.ch/user/l/lwming/public/74X/GenHFHadronMatcher.cc PhysicsTools/JetMCAlgos/plugins/GenHFHadronMatcher.cc
    cp /afs/cern.ch/user/l/lwming/public/74X/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/GenTtbarCategorizer.cc
    cp /afs/cern.ch/user/l/lwming/public/74X/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/GenTtbarCategorizer_cfi.py

    git cms-merge-topic gkasieczka:htt-v2-74X
    git clone https://github.com/cms-ttH/BoostedTTH.git
    ln -s $CMSSW_RELEASE_BASE/src/RecoJets/JetProducers/plugins/VirtualJetProducer.h BoostedTTH/BoostedProducer/plugins/VirtualJetProducer.h
    ln -s $CMSSW_RELEASE_BASE/src/RecoJets/JetProducers/plugins/VirtualJetProducer.cc BoostedTTH/BoostedProducer/plugins/VirtualJetProducer.cc

    git clone https://github.com/cms-ttH/MiniAOD.git
    git clone https://github.com/cms-ttH/ttH-LeptonPlusJets.git

    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    git clone https://github.com/cms-ttH/ttH-Limits.git ttH/Limits

    scram b -j 32


#############
 
