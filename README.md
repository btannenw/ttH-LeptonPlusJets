ttH-LeptonPlusJets
==================

### 76x

    cmsrel CMSSW_7_6_3
    cd CMSSW_7_6_3/src/
    cmsenv
    
    
    git clone https://github.com/cms-ttH/MiniAOD.git 
    cd MiniAOD/
    git checkout ccf0ec8
    cd ..
    
    wget https://raw.githubusercontent.com/hsatoshi/MiniAOD/09a6e2a2d20314afee38c6ccbcae5279217c7a90/BoostedObjects/src/BoostedUtils.cpp
    mv BoostedUtils.cpp MiniAOD/BoostedObjects/src
    
    wget https://raw.githubusercontent.com/hsatoshi/MiniAOD/09a6e2a2d20314afee38c6ccbcae5279217c7a90/BoostedObjects/interface/BoostedUtils.hpp
    mv BoostedUtils.hpp MiniAOD/BoostedObjects/interface
    
    
    # For the moment, do below instead of git@github.com:cms-ttH/ttH-LeptonPlusJets.git -b 76x  
    git clone git@github.com:hsatoshi/ttH-LeptonPlusJets.git -b 76x
    
    scram b 
