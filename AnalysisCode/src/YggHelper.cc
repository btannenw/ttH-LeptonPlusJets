#include "../interface/YggHelper.h"

using namespace std;

YggHelper::YggHelper(){

	h_lep_wgt_ele_trig = NULL;
	h_lep_wgt_ele_iso  = NULL;
	h_lep_wgt_ele_id   = NULL;
	h_lep_wgt_mu_trig1 = NULL;
	h_lep_wgt_mu_iso   = NULL;
	h_lep_wgt_mu_id    = NULL;
	h_lep_wgt_mu_trig2 = NULL;
	
	csvwp = 0.8;


}

YggHelper::~YggHelper(){

}

void YggHelper::fillCSVhistos(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;
}



double YggHelper::get_csv_wgt( std::vector<TLorentzVector> jets, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

  int iSysHF = 0;
  switch(iSys){
  case 7: iSysHF=1; break;
  case 8: iSysHF=2; break;
  case 9: iSysHF=3; break;
  case 10: iSysHF=4; break;
  case 13: iSysHF=5; break;
  case 14: iSysHF=6; break;
  case 15: iSysHF=7; break;
  case 16: iSysHF=8; break;
  default : iSysHF = 0; break;
  }

  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }

  int iSysLF = 0;
  switch(iSys){
  case 7: iSysLF=1; break;
  case 8: iSysLF=2; break;
  case 11: iSysLF=3; break;
  case 12: iSysLF=4; break;
  case 17: iSysLF=5; break;
  case 18: iSysLF=6; break;
  case 19: iSysLF=7; break;
  case 20: iSysLF=8; break;
  default : iSysLF = 0; break;
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;


  for( int iJet=0; iJet<int(jets.size()); iJet++ ){
    TLorentzVector myJet;
    myJet = jets[iJet];
  
	  
    double csv = jetCSV[iJet];
    double jetPt = myJet.Pt();
    double jetAbsEta = fabs( myJet.Eta() );
    int flavor = jetFlavor[iJet];

    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100)iPt=4;
   // else if (jetPt >=100 && jetPt<160) iPt = 4;
    //else if (jetPt >=160 && jetPt<10000) iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //      flavor, jetPt, iJet->eta, csv, iCSVWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;

      // if( iSysLF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtLF );
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;
}

void YggHelper::FillLeptonWgtHistos(TFile *fileEleTrig, TFile *fileEleID, TFile *fileEleIso, TFile *fileMuonTrig, TFile *fileMuonID, TFile *fileMuonIso){



    h_lep_wgt_ele_trig   = (TH2D*)fileEleTrig->Get( Form("h_eleTrig_SF") );
    h_lep_wgt_ele_id     = (TH2D*)fileEleID->Get( Form("EGamma_SF2D"));
    h_lep_wgt_ele_iso    = (TH2D*)fileEleIso->Get( Form("IsolationSF"));
    h_lep_wgt_mu_trig1   = (TH2D*)fileMuonTrig->Get( Form("runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio"));
    h_lep_wgt_mu_id      = (TH2D*)fileMuonID->Get(Form("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio"));
    h_lep_wgt_mu_iso     = (TH2D*)fileMuonIso->Get(Form("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio"));
    h_lep_wgt_mu_trig2   = (TH2D*)fileMuonTrig->Get(Form(" runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio"));

   
  return;
}

double YggHelper::getLeptonWgt(vint lepnums, std::vector<TLorentzVector> vlepton){

  double LeptonWgtTotal =1;
  int leptonType = -1;
  
  if(lepnums[1]==1)leptonType=1; //Lepton is Electron
  if(lepnums[3]==1)leptonType=2; //Lepton is Muon

  
  

  if(leptonType==1){
  double elePt = vlepton[0].Pt();
  double eleEta = vlepton[0].Eta();

	
	//ElectronSF ID
	
	int iptEID=-1;
	int ietaEID=-1;
	double eleIDWgt = 1;

	if(elePt>=14.99 && elePt<25)iptEID=1;
	if(elePt>=25 && elePt<35)iptEID=2;
	if(elePt>=35 && elePt<45)iptEID=3;
	if(elePt>=35 && elePt<55)iptEID=4;
	if(elePt>=55)iptEID=5;

	if(eleEta>=-2.5 && eleEta<-1.57)ietaEID=1;
	if(eleEta>=-1.57 && eleEta<-1.44)ietaEID=2;
	if(eleEta>=-1.44 && eleEta<-0.8)ietaEID=3;
	if(eleEta>=-0.8 && eleEta<0)ietaEID=4;
	if(eleEta>=0 && eleEta<0.8)ietaEID=5;
	if(eleEta>=0.8 && eleEta<1.44)ietaEID=6;
	if(eleEta>=1.44 && eleEta<1.57)ietaEID=7;
	if(eleEta>=1.57 && eleEta<=2.5)ietaEID=8;

	if(iptEID>0 && ietaEID>0)eleIDWgt = h_lep_wgt_ele_id->GetBinContent(ietaEID,iptEID);
	else cout<<"Could not find eleID bin for electron with pt, eta = "<<elePt<<" "<<eleEta<<endl;
	
	//
	
	//ElectronSF Trigger
	
	int iptETrig=-1;
	int ietaETrig=-1;
	double eleTrigWgt = 1;

	//if(elePt>=14.99 && elePt<20)iptETrig=1;
	//if(elePt>=20 && elePt<25)iptETrig=2;
	if(elePt>=25 && elePt<30)iptETrig=1;
	if(elePt>=30 && elePt<35)iptETrig=2;
	if(elePt>=35 && elePt<40)iptETrig=3;
	if(elePt>=40 && elePt<50)iptETrig=4;
	if(elePt>=50 && elePt<60)iptETrig=5;
	if(elePt>=60 && elePt<100)iptETrig=6;
	if(elePt>=100)iptETrig=7;

	if(eleEta>=-2.5 && eleEta<-2.1)ietaETrig=1;
	if(eleEta>=-2.1 && eleEta<-1.57)ietaETrig=2;
	if(eleEta>=-1.57 && eleEta<-1.44)ietaETrig=3;
	if(eleEta>=-1.44 && eleEta<-0.8)ietaETrig=4;
	if(eleEta>=-0.8 && eleEta<0)ietaETrig=5;
	if(eleEta>=0 && eleEta<0.8)ietaETrig=6;
	if(eleEta>=0.8 && eleEta<1.44)ietaETrig=7;
	if(eleEta>=1.44 && eleEta<1.57)ietaETrig=8;
	if(eleEta>=1.57 && eleEta<2.1)ietaETrig=9;
	if(eleEta>=2.1 && eleEta<=2.5)ietaETrig=10;
	
	//cout<<" $ "<<iptETrig<<" "<<ietaETrig<<" $ ";
	
	//if(ietaETrig ==5 || ietaETrig ==6)eleTrigWgt=1;

	if(iptETrig>0 && ietaETrig>0)eleTrigWgt = h_lep_wgt_ele_trig->GetBinContent(iptETrig,ietaETrig);
	else cout<<"Could not find eleTrig bin for electron with pt, eta = "<<elePt<<" "<<eleEta<<endl;


  
	
	//ElectronSF Iso
	
	int iptEIso=-1;
	int ietaEIso=-1;
	double eleIsoWgt = 1;

	if(elePt>=14.99 && elePt<20)iptEIso=1;
	if(elePt>=20 && elePt<30)iptEIso=2;
	if(elePt>=30 && elePt<40)iptEIso=3;
	if(elePt>=40 && elePt<50)iptEIso=4;
	if(elePt>=50)iptEIso=5;
	
	
	if(eleEta>=-2.5 && eleEta<-1.57)ietaEIso=1;
	if(eleEta>=-1.57 && eleEta<-1.44)ietaEIso=2;
	if(eleEta>=-1.44 && eleEta<-0.8)ietaEIso=3;
	if(eleEta>=-0.8 && eleEta<0)ietaEIso=4;
	if(eleEta>=0 && eleEta<0.8)ietaEIso=5;
	if(eleEta>=0.8 && eleEta<1.44)ietaEIso=6;
	if(eleEta>=1.44 && eleEta<1.57)ietaEIso=7;
	if(eleEta>=1.57 && eleEta<=2.5)ietaEIso=8;
	
	if(iptEIso>0 && ietaEIso>0)eleIsoWgt = h_lep_wgt_ele_iso->GetBinContent(ietaEIso,ietaEIso);
	else cout<<"Could not find eleIso bin for electron with pt, eta = "<<elePt<<" "<<eleEta<<endl;

	
	double eleWgt = eleTrigWgt * eleIDWgt * eleIsoWgt;
	LeptonWgtTotal = eleWgt;
  }

  if(leptonType==2){
  double muPt = vlepton[0].Pt();
  double muEta = fabs(vlepton[0].Eta());
  
  	int iptM = -1;
	int ietaM =-1;
  	
	
	double muTrigWgt  = 1;
	//double muTrigWgt1 = 1;
	//double muTrigWgt2 = 1;
	double muIDWgt    = 1;
	double muIsoWgt   = 1;
	
	
	if(muPt>=24.99 && muPt<25)iptM=1;
	if(muPt>=25 && muPt<30)iptM=2;
	if(muPt>=30 && muPt<40)iptM=3;
	if(muPt>=40 && muPt<50)iptM=4;
	if(muPt>=50 && muPt<60)iptM=5;
	if(muPt>=60)iptM=6;
	
	if(muEta>=0 && muEta<0.9)ietaM=1;
	if(muEta>=0.9 && muEta<1.2)ietaM=2;
	if(muEta>=1.2 && muEta<2.1)ietaM=3;
	if(muEta>=2.1)ietaM=4;
	
	if(iptM>0 && ietaM>0){
		//muTrigWgt1 = h_lep_wgt_mu_trig1->GetBinContent(ietaM,iptM);
		//muTrigWgt2 = h_lep_wgt_mu_trig2->GetBinContent(ietaM,iptM);
		//muIsoWgt   = h_lep_wgt_mu_iso->GetBinContent(ietaM,iptM);
		//muIDWgt = h_lep_wgt_mu_id->GetBinContent(ietaM,iptM);
		//muTrigWgt = (muTrigWgt1*754.394 + muTrigWgt2*1908.10)/(754.394+1908.10);
	}
	else cout<<"Could not find muSF bin for muon with pt, eta = "<<muPt<<" "<<muEta<<endl;
	
	double muWgt = muTrigWgt * muIDWgt * muIsoWgt;
	LeptonWgtTotal = muWgt;
  }
  
  return LeptonWgtTotal;
}


//Feed in Loose Jets and get a new Jet vecTLV. For PT/ETA optimization studies.

void YggHelper::remakeJetV(std::vector<TLorentzVector> JetV, vdouble csvV, double JetPTcut, double JetEtacut, std::vector<TLorentzVector> &NewJetV, vdouble &NewJetcsvV, int &NewnumJets, int &NewtagJets ){
	
	NewnumJets = 0;
	NewtagJets = 0;
	for( int iJet=0; iJet<int(JetV.size()); iJet++ ){
	  TLorentzVector myJet;
	  myJet = JetV[iJet];
	  double myJetPT = myJet.Pt();
	  double myCSV = csvV[iJet];
	  double myJetEta = fabs(myJet.Eta());
	  if(myJetPT<JetPTcut)continue;
	  if(myJetEta>JetEtacut)continue;
	  NewJetV.push_back(myJet);
	  NewnumJets++;
	  NewJetcsvV.push_back(myCSV);	
	  if(myCSV>csvwp)NewtagJets++;	  
	}
		
	return;	
		
}

//Feed in Loose leptons and get a new lepton vecTLV. For PT/ETA optimization studies.

void YggHelper::remakeLeptonV(std::vector<TLorentzVector> LeptonV, double eleCut, double muonCut, std::vector<TLorentzVector> &NewLeptonV){
	
	for(int iLep=0; iLep<int(LeptonV.size()); iLep++){
		TLorentzVector myLep;
		myLep = LeptonV[iLep];
		double myLepPT = myLep.Pt();
		double isMu = LeptonV[iLep][4];
		if(isMu && myLepPT<muonCut)continue;
		if(!isMu && myLepPT<eleCut)continue;
  	    NewLeptonV.push_back(myLep);	
	}
	
	return;
}





int YggHelper::passCuts(bool goodfirstPV, int ElectronTrigger, int MuonTrigger, vint lepnums, int numJet, int numTag){
	
	int cutsPassed=1;
	
	if(goodfirstPV){;
	cutsPassed++;
	
	
	int numLooseElectrons = lepnums[0];
	int numTightElectrons = lepnums[1];
	int numLooseMuons = lepnums[2];
	int numTightMuons = lepnums[3];
	
	if(ElectronTrigger==1 && MuonTrigger==0){ //singleElectronTrigger has Passed
		cutsPassed++;
		if(numLooseElectrons==1){
			cutsPassed++;
			if(numTightElectrons==1){
				cutsPassed++;
				if(numLooseMuons==0){
					cutsPassed++;
					if(numJet>=4 && numTag>=2){
						cutsPassed++;
					}
				}
			}
		}
	}
	if(ElectronTrigger==0 && MuonTrigger==1){ //singleMuonTrigger has Passed
		cutsPassed++;
		if(numLooseMuons==1){
			cutsPassed++;
			if(numTightMuons==1){
				cutsPassed++;
				if(numLooseElectrons==0){
					cutsPassed++;
					if(numJet>=4 && numTag>=2){
						cutsPassed++;
					}
				}
			}
		}
	}
	if(ElectronTrigger==1 && MuonTrigger==1){ //singleElectronTrigger and singleMuonTrigger have Passed
		cutsPassed++;
		if(numLooseElectrons==1 || numLooseMuons==1){
			cutsPassed++;
			if((numLooseElectrons==1&&numTightElectrons==1)||(numLooseMuons==1&&numTightMuons==1)){
				cutsPassed++;
				if((numLooseElectrons==1&&numTightElectrons==1&&numLooseMuons==0)||(numLooseMuons==1&&numTightMuons==1&&numLooseElectrons==0)){
					cutsPassed++;
					if(numJet>=4 && numTag>=2){
						cutsPassed++;
					}
				}
			}
		}
	}
	}
	
	return cutsPassed;
	
	
}



void YggHelper::makeVTLVjets(std::vector<TLorentzVector> &NewJetV, JetCollection inputJets, int &NewNumJets, int &NewNumTags, ThresholdList cuts){
	
	NewNumJets = 0;
	NewNumTags = 0;
	for( int iJet=0; iJet<int(inputJets.jet_pt.size()); iJet++ ){
	  TLorentzVector vjets_new;
	  vjets_new.SetPtEtaPhiM(inputJets.jet_pt.at(iJet),inputJets.jet_eta.at(iJet),inputJets.jet_phi.at(iJet),inputJets.jet_m.at(iJet));
	  NewJetV.push_back(vjets_new);
	  NewNumJets++;
	  if(inputJets.jet_bDiscriminant.at(iJet)>cuts.Thre_Jet_Btag)NewNumTags++;	  
	}
		
	return;	
		
}



void YggHelper::makeVTLVleptons(std::vector<TLorentzVector> &NewLeptonV, vint &lepnums, LeptonCollection inputTightLeptons, LeptonCollection inputLooseLeptons){

	int numLooseElectrons=0;
	int numTightElectrons=0;
	int numLooseMuons=0;
	int numTightMuons=0;
	
	
	for(unsigned int iLep=0; iLep<inputTightLeptons.lep_pt.size(); iLep++){
		int isMu =inputTightLeptons.lep_isMuon.at(iLep);
		if(isMu)numTightMuons++;
		if(!isMu)numTightElectrons++;
		TLorentzVector vlep_new;
		vlep_new.SetPtEtaPhiE(inputTightLeptons.lep_pt.at(iLep),inputTightLeptons.lep_eta.at(iLep),inputTightLeptons.lep_phi.at(iLep),inputTightLeptons.lep_e.at(iLep));
		NewLeptonV.push_back(vlep_new);
	}
	
	for(unsigned int iLep=0; iLep<inputLooseLeptons.lep_pt.size(); iLep++){
		int isMu = inputLooseLeptons.lep_isMuon.at(iLep);
		if(isMu)numLooseMuons++;
		if(!isMu)numLooseElectrons++;
	}
		
	
	
	lepnums.push_back(numLooseElectrons);
	lepnums.push_back(numTightElectrons);
	lepnums.push_back(numLooseMuons);
	lepnums.push_back(numTightMuons);
	
	return;
}



void YggHelper::SetCuts(ThresholdList &cuts, double tightmupt, double loosemupt, double tightmueta, double loosemueta, double tightmuiso, double loosemuiso, double tightelpt, double looseelpt, double tighteleta, double looseeleta, double tighteliso, double looseeliso, double jetpt, double jeteta, double jetbtag){

	cuts.Thre_TightMu_PT = tightmupt;
	cuts.Thre_LooseMu_PT = loosemupt;
	cuts.Thre_TightMu_Eta = tightmueta;
	cuts.Thre_LooseMu_Eta = loosemueta;
	cuts.Thre_TightMu_Iso = tightmuiso;
	cuts.Thre_LooseMu_Iso = loosemuiso;
	cuts.Thre_TightEl_PT = tightelpt;
	cuts.Thre_LooseEl_PT = looseelpt;
	cuts.Thre_TightEl_Eta = tighteleta;
	cuts.Thre_LooseEl_Eta = looseeleta;
	cuts.Thre_TightEl_Iso = tighteliso;
	cuts.Thre_LooseEl_Iso = looseeliso;
	cuts.Thre_Jet_PT = jetpt;
	cuts.Thre_Jet_Eta = jeteta;
	cuts.Thre_Jet_Btag = jetbtag;


	return;
}



void YggHelper::SelectLeptons(LeptonCollection inputLeptons, LeptonCollection &outputTightLeptons, LeptonCollection &outputLooseLeptons, ThresholdList cuts){


	ElectronSelection(inputLeptons, outputTightLeptons, outputLooseLeptons, cuts);
	MuonSelection(inputLeptons, outputTightLeptons, outputLooseLeptons, cuts);
	outputTightLeptons.SortLeptons();
	outputLooseLeptons.SortLeptons();

	return;
}


void YggHelper::ElectronSelection(LeptonCollection inputLeptons, LeptonCollection &outputTightLeptons, LeptonCollection &outputLooseLeptons, ThresholdList cuts){

	for( unsigned int idx = 0 ; idx < inputLeptons.lep_pt.size() ; idx++){
		
		if(  inputLeptons.lep_isMuon.at(idx) !=0)continue;
		
		if(  inputLeptons.lep_POGTight.at(idx) !=1)continue;
		if(  inputLeptons.lep_pt.at(idx)        < cuts.Thre_LooseEl_PT)  continue;
		if( fabs (inputLeptons.lep_eta.at(idx)) > cuts.Thre_LooseEl_Eta)  continue;
		if( inputLeptons.lep_relIso.at(idx)     > cuts.Thre_LooseEl_Iso) continue;
		
		outputLooseLeptons.AddLep(inputLeptons,idx);
		
		if( inputLeptons.lep_POGTight.at(idx) !=1)continue;
		if( inputLeptons.lep_pt.at(idx)           < cuts.Thre_TightEl_PT ) continue;
		if( fabs(inputLeptons.lep_eta.at(idx))    > cuts.Thre_TightEl_Eta) continue;
		if( fabs(inputLeptons.lep_relIso.at(idx)) > cuts.Thre_TightEl_Iso) continue;
		
		outputTightLeptons.AddLep(inputLeptons,idx);
		
		
	}


	return;
}


void YggHelper::MuonSelection(LeptonCollection inputLeptons, LeptonCollection &outputTightLeptons, LeptonCollection &outputLooseLeptons, ThresholdList cuts){


	for( unsigned int idx = 0; idx < inputLeptons.lep_pt.size() ; idx++ ){
	
		if( inputLeptons.lep_isMuon.at(idx) !=1 ) continue;
		
		if(  inputLeptons.lep_POGTight.at(idx) !=1) continue;
		if(  inputLeptons.lep_pt.at(idx)        < cuts.Thre_LooseMu_PT ) continue;
		if(  fabs(inputLeptons.lep_eta.at(idx)) > cuts.Thre_LooseMu_Eta) continue;
		if(  inputLeptons.lep_relIso.at(idx)   > cuts.Thre_LooseMu_Iso ) continue;
		
		outputLooseLeptons.AddLep(inputLeptons,idx);
		
		if(  inputLeptons.lep_POGTight.at(idx) !=1) continue;
		if(  inputLeptons.lep_pt.at(idx)           < cuts.Thre_TightMu_PT ) continue;
		if(  fabs(inputLeptons.lep_eta.at(idx))    > cuts.Thre_TightMu_Eta) continue;
		if(  fabs(inputLeptons.lep_relIso.at(idx)) > cuts.Thre_TightMu_Iso) continue;
		
		outputTightLeptons.AddLep(inputLeptons,idx);
	}



	return;
}



void YggHelper::SelectJets(JetCollection inputJets, JetCollection &outputJets, LeptonCollection inputLooseLeptons, ThresholdList cuts){

	for(unsigned int idx = 0; idx < inputJets.jet_pt.size() ; idx++){
	
	
		if( inputJets.jet_pt.at(idx) < cuts.Thre_Jet_PT) continue;
		if( fabs(inputJets.jet_eta.at(idx)) > cuts.Thre_Jet_Eta) continue;
		
		if(OverlapWithLooseLeptons(inputJets.jet_eta.at(idx),inputJets.jet_phi.at(idx),inputLooseLeptons)) continue;
		
		outputJets.AddJet(inputJets,idx);
		
	}


	return;
}



bool YggHelper::OverlapWithLooseLeptons(double jeta, double jphi, LeptonCollection inputLoose){
		bool overlap = false;
		for(unsigned int iLep=0; iLep<inputLoose.lep_pt.size();iLep++){
			double d_eta = jeta-inputLoose.lep_eta.at(iLep);
			double d_phi = fabs(jphi - inputLoose.lep_phi.at(iLep));
			d_phi = (d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi;
			double DR2 = d_eta*d_eta + d_phi*d_phi;
			if(DR2 < 0.4*0.4)overlap = true;
		}
	return overlap;
}

double YggHelper::reweightPU( int nPU, int iSys ){

  double PUscale[50];

  if( nPU>49 ) nPU = 49;

  if( iSys==0 ){
  PUscale[0] = 0.544141;
  PUscale[1] = 0.683152;
  PUscale[2] = 1.0588;
  PUscale[3] = 1.36648;
  PUscale[4] = 1.62214;
  PUscale[5] = 1.93965;
  PUscale[6] = 1.45916;
  PUscale[7] = 1.29301;
  PUscale[8] = 1.39344;
  PUscale[9] = 1.37393;
  PUscale[10] = 1.26638;
  PUscale[11] = 1.1637;
  PUscale[12] = 1.05777;
  PUscale[13] = 0.901478;
  PUscale[14] = 0.695434;
  PUscale[15] = 0.48665;
  PUscale[16] = 0.323436;
  PUscale[17] = 0.233179;
  PUscale[18] = 0.204042;
  PUscale[19] = 0.179695;
  PUscale[20] = 0.117906;
  PUscale[21] = 0.0535705;
  PUscale[22] = 0.0179761;
  PUscale[23] = 0.00486845;
  PUscale[24] = 0.00120473;
  PUscale[25] = 0.000306356;
  PUscale[26] = 8.44097e-05;
  PUscale[27] = 2.35712e-05;
  PUscale[28] = 6.33309e-06;
  PUscale[29] = 1.62779e-06;
  PUscale[30] = 3.92723e-07;
  PUscale[31] = 8.94453e-08;
  PUscale[32] = 1.90178e-08;
  PUscale[33] = 3.84203e-09;
  PUscale[34] = 7.38767e-10;
  PUscale[35] = 1.28456e-10;
  PUscale[36] = 2.08504e-11;
  PUscale[37] = 3.31182e-12;
  PUscale[38] = 4.63316e-13;
  PUscale[39] = 4.72854e-14;
  PUscale[40] = 0;
  PUscale[41] = 0;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;
  }
  else if( iSys==1 ){
  PUscale[0] = 0.458261;
  PUscale[1] = 0.626711;
  PUscale[2] = 0.971785;
  PUscale[3] = 1.21776;
  PUscale[4] = 1.34764;
  PUscale[5] = 1.4248;
  PUscale[6] = 0.922526;
  PUscale[7] = 0.744979;
  PUscale[8] = 0.913723;
  PUscale[9] = 1.06132;
  PUscale[10] = 1.11114;
  PUscale[11] = 1.12292;
  PUscale[12] = 1.13425;
  PUscale[13] = 1.10829;
  PUscale[14] = 1.00558;
  PUscale[15] = 0.83577;
  PUscale[16] = 0.639871;
  PUscale[17] = 0.469319;
  PUscale[18] = 0.361301;
  PUscale[19] = 0.306711;
  PUscale[20] = 0.237734;
  PUscale[21] = 0.142901;
  PUscale[22] = 0.0656738;
  PUscale[23] = 0.0238299;
  PUscale[24] = 0.00722283;
  PUscale[25] = 0.00197294;
  PUscale[26] = 0.000550982;
  PUscale[27] = 0.000165108;
  PUscale[28] = 5.14307e-05;
  PUscale[29] = 1.60005e-05;
  PUscale[30] = 4.76288e-06;
  PUscale[31] = 1.35321e-06;
  PUscale[32] = 3.62213e-07;
  PUscale[33] = 9.29401e-08;
  PUscale[34] = 2.28987e-08;
  PUscale[35] = 5.14682e-09;
  PUscale[36] = 1.08953e-09;
  PUscale[37] = 2.27672e-10;
  PUscale[38] = 4.24244e-11;
  PUscale[39] = 6.6564e-12;
  PUscale[40] = 1.21726e-12;
  PUscale[41] = 1.25024e-13;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;
  }
  else if( iSys==-1 ){
  PUscale[0] = 0.637352;
  PUscale[1] = 0.745524;
  PUscale[2] = 1.16977;
  PUscale[3] = 1.55999;
  PUscale[4] = 2.01896;
  PUscale[5] = 2.79883;
  PUscale[6] = 2.50005;
  PUscale[7] = 2.22714;
  PUscale[8] = 2.00474;
  PUscale[9] = 1.67907;
  PUscale[10] = 1.37446;
  PUscale[11] = 1.13578;
  PUscale[12] = 0.899819;
  PUscale[13] = 0.646461;
  PUscale[14] = 0.412783;
  PUscale[15] = 0.243661;
  PUscale[16] = 0.153684;
  PUscale[17] = 0.125858;
  PUscale[18] = 0.116543;
  PUscale[19] = 0.0860169;
  PUscale[20] = 0.0414331;
  PUscale[21] = 0.0133115;
  PUscale[22] = 0.00326222;
  PUscale[23] = 0.000723524;
  PUscale[24] = 0.000169442;
  PUscale[25] = 4.20509e-05;
  PUscale[26] = 1.04015e-05;
  PUscale[27] = 2.41401e-06;
  PUscale[28] = 5.19674e-07;
  PUscale[29] = 1.0525e-07;
  PUscale[30] = 1.97832e-08;
  PUscale[31] = 3.47419e-09;
  PUscale[32] = 5.63765e-10;
  PUscale[33] = 8.60403e-11;
  PUscale[34] = 1.23701e-11;
  PUscale[35] = 1.59276e-12;
  PUscale[36] = 1.89108e-13;
  PUscale[37] = 1.98468e-14;
  PUscale[38] = 0;
  PUscale[39] = 0;
  PUscale[40] = 0;
  PUscale[41] = 0;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;
  }

  return PUscale[nPU];
}


void JetCollection::SetJetCollection(vdouble j_pt, vdouble j_eta, vdouble j_phi, vdouble j_m, vdouble j_csv){
	jet_pt=j_pt;
	jet_eta=j_eta;
	jet_phi=j_phi;
	jet_m = j_m;
	jet_bDiscriminant = j_csv;
}

void JetCollection::AddJet(JetCollection JetCollectionB, int iJet){
	jet_pt.push_back(JetCollectionB.jet_pt.at(iJet));
	jet_eta.push_back(JetCollectionB.jet_eta.at(iJet));
	jet_phi.push_back(JetCollectionB.jet_phi.at(iJet));
	jet_m.push_back(JetCollectionB.jet_m.at(iJet));
	jet_bDiscriminant.push_back(JetCollectionB.jet_bDiscriminant.at(iJet));
}

void LeptonCollection::SetLeptonCollection(vdouble l_pt, vdouble l_eta, vdouble l_phi, vdouble l_e, vint l_charge, vint l_isMuon, vdouble l_relIso, vint l_isTight, vint l_isLoose){
	lep_pt = l_pt;
	lep_eta = l_eta;
	lep_phi = l_phi;
	lep_e = l_e;
	lep_charge = l_charge;
	lep_isMuon = l_isMuon;
	lep_relIso = l_relIso;
	lep_POGLoose = l_isLoose;
	lep_POGTight = l_isTight;
}

void LeptonCollection::AddLep(LeptonCollection LeptonCollectionB,int iLep){
	lep_pt.push_back(LeptonCollectionB.lep_pt.at(iLep));
	lep_eta.push_back(LeptonCollectionB.lep_eta.at(iLep));
	lep_phi.push_back(LeptonCollectionB.lep_phi.at(iLep));
	lep_e.push_back(LeptonCollectionB.lep_e.at(iLep));
	lep_charge.push_back(LeptonCollectionB.lep_charge.at(iLep));
	lep_isMuon.push_back(LeptonCollectionB.lep_isMuon.at(iLep));
	lep_relIso.push_back(LeptonCollectionB.lep_relIso.at(iLep));
	lep_POGLoose.push_back(LeptonCollectionB.lep_POGLoose.at(iLep));
	lep_POGTight.push_back(LeptonCollectionB.lep_POGTight.at(iLep));
}

void LeptonCollection::ClearLep(){
	lep_pt.clear();
	lep_eta.clear();
	lep_phi.clear();
	lep_e.clear();
	lep_charge.clear();
	lep_isMuon.clear();
	lep_relIso.clear();
	lep_POGLoose.clear();
	lep_POGTight.clear();
}

void LeptonCollection::SortLeptons(){

	//LeptonCollection dummy;
	
	for(unsigned int i=0; i< lep_pt.size(); i++){
	for(unsigned int j=0; j< lep_pt.size(); j++){
		if( lep_pt.at(i) > lep_pt.at(j) )continue;
		
		double dpt = lep_pt.at(i);
		double deta = lep_eta.at(i);
		double dphi = lep_phi.at(i);
		double de = lep_e.at(i);
		int dcharge = lep_charge.at(i);
		int dmuon = lep_isMuon.at(i);
		double diso = lep_relIso.at(i);
		int dloose = lep_POGLoose.at(i);
		int dtight = lep_POGTight.at(i);
		
		lep_pt.at(i) = lep_pt.at(j);
		lep_eta.at(i) = lep_eta.at(j);
		lep_phi.at(i) = lep_phi.at(j);
		lep_e.at(i) = lep_e.at(j);
		lep_charge.at(i) = lep_charge.at(j);
		lep_isMuon.at(i) = lep_isMuon.at(j);
		lep_relIso.at(i) = lep_relIso.at(j);
		lep_POGLoose.at(i) = lep_POGLoose.at(j);
		lep_POGTight.at(i) = lep_POGTight.at(j);
		
		lep_pt.at(j) = dpt;
		lep_eta.at(j) = deta;
		lep_phi.at(j) = dphi;
		lep_e.at(j) = de;
		lep_charge.at(j) = dcharge;
		lep_isMuon.at(j) = dmuon;
		lep_relIso.at(j) = diso;
		lep_POGLoose.at(j) = dloose;
		lep_POGTight.at(j) = dtight;
	}
	}	
}

