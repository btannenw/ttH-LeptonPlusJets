ttH-LeptonPlusJets
==================

+++++++++++++++++
+++  Install  +++
+++++++++++++++++

# csh 
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc530
# bash 
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530

cmsrel CMSSW_8_0_24_patch1
cd CMSSW_8_0_24_patch1/src/
cmsenv

export CMSSW_SRC="$( pwd )"

# updated MET tools
git cms-init
git cms-merge-topic cms-met:METRecipe_8020
git cms-merge-topic cms-met:fromCMSSW_8_0_20_postICHEPfilter
git cms-merge-topic ahinzmann:METRecipe_8020_Moriond17


scram b -j 10 
# use "scam b" at least once to create directory "../external" which we use below.

# update electron MVA tools and data
git cms-merge-topic ikrav:egm_id_80X_v2
cd ../external/$SCRAM_ARCH
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout egm_id_80X_v1
cd $CMSSW_SRC


scram b -j 10 ;
# Let's check if setup CMSSW environment can be compiled before setting up our Yggdrasil code.
# You may need to repeat scram several times (not only twice)...


# - - - Our yggdrasil code - - - 

git clone https://github.com/hsatoshi/MiniAOD.git -b satoshi_8020
git clone https://github.com/hsatoshi/GenParticleTopOriginChargedleptonFilter.git ttHAnalysisSubprogram/GenParticleTopOriginChargedleptonFilter
git clone https://github.com/hsatoshi/PuppiLeptonIsolationhelper.git

scram b -j 10 ;
# You may need to repeat scram several times (not only twice)...

git clone https://github.com/cms-ttH/ttH-LeptonPlusJets.git
# or, use your branch, such as 
#   git clone git@github.com:hsatoshi/ttH-LeptonPlusJets.git

cd ttH-LeptonPlusJets/YggdrasilTreeMaker
scram b -j 10 ;
# You may need to repeat scram several times (not only twice)...

# If you submit crab job,
# 
cd $CMSSW_SRC
cp -rv ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1 RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1




+++++++++++++++++
+++ Selection +++
+++++++++++++++++


Selection of electron
 - slimmedElectrons from MiniAOD
 - pt > 20 GeV
 - |eta| < 2.4

Selection of muon
 - slimmedMuons from MiniAOD
 - POG loose ID
 - pt > 20 GeV
 - |eta| < 2.4

Selection of jets
 - slimmedJets from MiniAOD 
 - POG loose Jet ID (jetID::jetLoose) with MiniAODHelper
 - Un-correct with MiniAODHelper
 - Correct    with MiniAODHelper
 - PT > 20 GeV, |eta|<5.
 - Stored in PT-order.


Selection of gen jets
 - this is independent from reco jet.
 - Input can be set with option "genjet", but usually 
 - PT cut 8 GeV.
 - no Eta cut


+++++++++++++++++
+++ Variables +++
+++++++++++++++++

Event

 run_  : run number
 lumi_ : lumi block
 evt_  : event number

 numTruePV_ : getTrueNumInteractions()
 numGenPV_  : getPU_NumInteractions()

 GoodFirstPV_ : if the LV is good (it checks !isFake, ndof, z and rho.)
 numPVs_ : number of good vertex.

 numSys_ : number of systematics in the tree.

 additionalJetEventId_;

 higgsDecayType_;
 ttbarDecayType_;


Trigger 

passHLT_XXXX : trigger for the path. 1=pass



Lepton.

lepton_pt_      : pt
lepton_eta_     : eta
lepton_phi_     : phi
lepton_e_       : energy
lepton_isMuon_  : If muon, 1. If electron, 0. 
lepton_relIso_  : Muon, delta_beta corrected relative isolation. cone 04. Calculated with MiniAODHelper.
                : Ele , EffectiveArea-corrected isolation, cone 0.3 (Calculated with MiniAODHelper with effAreaType::spring15)
lepton_puppirelIso_  : Muon , PUPPI isolation
                     : Ele  , (same as lepton_relIso_ for the moment.)
lepton_isTight_ : POG Tight ID. 1=pass. [*1] 
lepton_isLoose_ : POG Loose ID. 1=pass. [*1]
lepton_scEta_   : Electron Super Cluster. For muon, -99 is filled.

[*1] POG lepton ID : 
 Muon     POG Tight ID : = Tight ID.
 Muon     POG loose ID is not used at the moment. Always 1.
 Electron POG Tight = GeneralPurposeMVA2016 WP80
 Electron POG Loose = GeneralPurposeMVA2016 WP90

 lepnums_  : number of tight/loose leptons [*2]
[2*] not used in YggdrasilTupleMaker(and not filled) but needed by "AnalysisCode/macros/Yggdrasil_Slim.C"

Reconstructed ak4 jet.
 
jet_pt_  : pt
jet_phi_ : phi
jet_eta_ : eta
jet_m_   : mass
jet_combinedMVABJetTags_ : "pfCombinedMVAV2BJetTags"
jet_combinedInclusiveSecondaryVertexV2BJetTags_ : "pfCombinedInclusiveSecondaryVertexV2BJetTags"

jet_combinedMVABJetTags_HIP_ : "pfCombinedMVAV2BJetTags" HIP mitigation
jet_combinedInclusiveSecondaryVertexV2BJetTags_HIP_ : "pfCombinedInclusiveSecondaryVertexV2BJetTags" HIP mitigation

jet_partonflavour_ : result of jet->partonFlavour()
jet_flavour_       : result of jet->hadronFlavour()

jet_AssociatedGenJet_pt_ : genJet associated with the reco jet
jet_AssociatedGenJet_eta_: genJet associated with the reco jet
jet_AssociatedGenJet_phi_: genJet associated with the reco jet
jet_AssociatedGenJet_m_  : genJet associated with the reco jet

jet_genId_
jet_genParentId_
jet_genGrandParentId_

jet_vtxMass_
jet_vtxNtracks_
jet_vtx3DVal_
jet_vtx3DSig_


GenJet ak4 (independent from reco jet):

genjet_pt_ : pt
genjet_eta_: eta
genjet_phi_: phi
genjet_m_  : mass
genjet_BhadronMatch_ : =1 if this gen jet is regarded as a B-hadron with by GenHFHadronMatcher --[*2]

[*2] GenHFHadronMatcher in PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff
