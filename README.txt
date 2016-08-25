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

cmsrel CMSSW_8_0_12
cd CMSSW_8_0_12/src/
cmsenv

#[skip this for the moment]  # Apply bug fix of met correction. 
#[skip this for the moment]  #  - See https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription?rev=48#Instructions_for_8_0_X_X_10
#[skip this for the moment]  #  - Combination with CMSSW_8_0_13 did not work (June 30).
#[skip this for the moment]  git cms-init
#[skip this for the moment]  echo /PhysicsTools/PatUtils/ >> .git/info/sparse-checkout
#[skip this for the moment]  git cms-merge-topic cms-met:metTool80X   # ---[*1]

git clone https://github.com/cms-ttH/ttH-LeptonPlusJets.git

git clone https://github.com/hsatoshi/MiniAOD.git -b satoshi__CMSSW_8_0_8
git clone https://github.com/hsatoshi/GenParticleTopOriginChargedleptonFilter.git ttHAnalysisSubprogram/GenParticleTopOriginChargedleptonFilter
git clone https://github.com/hsatoshi/PuppiLeptonIsolationhelper.git

scram b -j 10 ; scram b -j 10
# Yes,.. compile twice for the moment. First compile claims errors while the second works. Need to be solved.




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
 Muon POG loose ID is not used at the moment. Always 1.
 Electron POG tight=MVA80. Loose=MVA90.

lepnums_  : number of tight/loose leptons [*2]

[2*] not used in YggdrasilTupleMaker(and not filled) but needed by "AnalysisCode/macros/Yggdrasil_Slim.C"

Reconstructed ak4 jet.
 
jet_pt_  : pt
jet_phi_ : phi
jet_eta_ : eta
jet_m_   : mass
jet_combinedMVABJetTags_ : "pfCombinedMVAV2BJetTags"
jet_combinedInclusiveSecondaryVertexV2BJetTags_ : "pfCombinedInclusiveSecondaryVertexV2BJetTags"

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
