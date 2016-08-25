{

  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy1 = makeshared.ReplaceAll("-W ", "");
  gSystem->SetMakeSharedLib(makeshared);
  //TString dummy = makeshared.ReplaceAll("-Wshadow ", "");
  TString dummy = makeshared.ReplaceAll("-Wshadow ", "-std=c++0x ");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libMiniAODMiniAODHelper.so");
  gSystem->Load("libttH-LeptonPlusJetsAnalysisCode.so");

  TString path = gSystem->GetIncludePath();
  path.Append(" -I$ROOTSYS/include/root -I./include  ");
  //path.Append(" -I$CMSSW_BASE/src/BEAN/BEANmaker/interface/ ");
  gSystem->SetIncludePath(path);


  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  gStyle->SetPalette(1);

}
