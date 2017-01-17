


echo ""
echo ""
echo "# - - - - - - - "
echo remove large files from directory to submit jobs.
rm -rf ../../AnalysisCode/data/*
echo " Removing large file is done..."
echo "# - - - - - - "

cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=False|isMC=True|g" > __yggdrasil_treeMaker_MC_cfg.py 
diff yggdrasil_treeMaker_cfg.py __yggdrasil_treeMaker_MC_cfg.py

cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=False|isMC=True|g" | sed "s|isTTBARMC=False|isTTBARMC=True|g" > __yggdrasil_treeMaker_MCTTBAR_cfg.py 
diff yggdrasil_treeMaker_cfg.py __yggdrasil_treeMaker_MCTTBAR_cfg.py

cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=True|isMC=False|g" > __yggdrasil_treeMaker_DATA_cfg.py 
diff yggdrasil_treeMaker_cfg.py __yggdrasil_treeMaker_DATA_cfg.py

nickname="Satoshi_Moriond17_Yggdra_20170117_submit2"


JobIndexList=""

ds[1]=/SingleElectron/Run2016B-23Sep2016-v3/MINIAOD
ds[2]=/SingleElectron/Run2016C-23Sep2016-v1/MINIAOD
ds[3]=/SingleElectron/Run2016D-23Sep2016-v1/MINIAOD
ds[4]=/SingleElectron/Run2016E-23Sep2016-v1/MINIAOD
ds[5]=/SingleElectron/Run2016F-23Sep2016-v1/MINIAOD
ds[6]=/SingleElectron/Run2016G-23Sep2016-v1/MINIAOD
ds[7]=/SingleElectron/Run2016H-PromptReco-v2/MINIAOD
ds[8]=/SingleElectron/Run2016H-PromptReco-v3/MINIAOD


name[1]=DataElB
name[2]=DataElC
name[3]=DataElD
name[4]=DataElE
name[5]=DataElF
name[6]=DataElG
name[7]=DataElHv2
name[8]=DataElHv3
ismc[1]=DATA
ismc[2]=DATA
ismc[3]=DATA
ismc[4]=DATA
ismc[5]=DATA
ismc[6]=DATA
ismc[7]=DATA
ismc[8]=DATA

JobIndexList=${JobIndexList}" 1 2 3 4 5 6 7 8 "

ds[11]=/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD
ds[12]=/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD
ds[13]=/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD
ds[14]=/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD
ds[15]=/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD
ds[16]=/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD
ds[17]=/SingleMuon/Run2016H-PromptReco-v2/MINIAOD
ds[18]=/SingleMuon/Run2016H-PromptReco-v3/MINIAOD

name[11]=DataMuB
name[12]=DataMuC
name[13]=DataMuD
name[14]=DataMuE
name[15]=DataMuF
name[16]=DataMuG
name[17]=DataMuHv2
name[18]=DataMuHv3
ismc[11]=DATA
ismc[12]=DATA
ismc[13]=DATA
ismc[14]=DATA
ismc[15]=DATA
ismc[16]=DATA
ismc[17]=DATA
ismc[18]=DATA

JobIndexList=${JobIndexList}" 11 12 13 14 15 16 17 18 "



ds[101]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT10-28028af67189b3de7224b79195bd0e1d/USER
ds[102]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT11-28028af67189b3de7224b79195bd0e1d/USER
ds[103]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT12-28028af67189b3de7224b79195bd0e1d/USER
ds[104]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT13-28028af67189b3de7224b79195bd0e1d/USER
ds[105]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT14-28028af67189b3de7224b79195bd0e1d/USER
ds[106]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT15-28028af67189b3de7224b79195bd0e1d/USER
ds[107]=/privateMCProductionAODSIMMiniAOD/friese-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT9-28028af67189b3de7224b79195bd0e1d/USER
ds[108]=/privateMCProductionAODSIMMiniAOD/kelmorab-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT1-28028af67189b3de7224b79195bd0e1d/USER
ds[109]=/privateMCProductionAODSIMMiniAOD/kelmorab-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT2-28028af67189b3de7224b79195bd0e1d/USER
ds[110]=/privateMCProductionAODSIMMiniAOD/kelmorab-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT3-28028af67189b3de7224b79195bd0e1d/USER
ds[111]=/privateMCProductionAODSIMMiniAOD/kelmorab-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT4-28028af67189b3de7224b79195bd0e1d/USER
ds[112]=/privateMCProductionAODSIMMiniAOD/kelmorab-eventAODSIMMiniAOD-TTToSemiLepton_hvq_ttHtranche3-KIT5-28028af67189b3de7224b79195bd0e1d/USER

name[101]=ttbar-KIT10
name[102]=ttbar-KIT11
name[103]=ttbar-KIT12
name[104]=ttbar-KIT13
name[105]=ttbar-KIT14
name[106]=ttbar-KIT15
name[107]=ttbar-KIT9

name[108]=ttbar-KIT1
name[109]=ttbar-KIT2
name[110]=ttbar-KIT3
name[111]=ttbar-KIT4
name[112]=ttbar-KIT5


ismc[101]=MC
ismc[102]=MC
ismc[103]=MC
ismc[104]=MC
ismc[105]=MC
ismc[106]=MC
ismc[107]=MC
ismc[108]=MC
ismc[109]=MC
ismc[110]=MC
ismc[111]=MC
ismc[112]=MC

# JobIndexList=${JobIndexList}" 101  102 103 104 105 106 107 108 109 110 111 112 "


ds[30]=/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
name[30]=ttbar
ismc[30]=MC

JobIndexList=${JobIndexList}" 30 "

# ds[21]=/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[21]=ttbar_dilep
# ismc[21]=MC


# ds[21]=/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
# name[21]=tthbb
# ismc[21]=MC
# 
# 
# 
# ds[22]=/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
# name[22]=tthnonbb
# ismc[22]=MC
# 

# ds[3]=/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# name[3]=DYlow
# ismc[3]=MC
# 

#ds[24]=/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
#name[24]=DYhigh
#ismc[24]=MC
#
#        
#
#ds[25]=/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/MINIAODSIM
#name[25]=WjetIncl
#ismc[25]=MC
# 
#ds[26]=/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/MINIAODSIM
#name[26]=ttbar
#ismc[26]=MC
#   # HLT = HLT, not HLT2
#   ds[5]=/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
#   name[5]=ttbar
#   ismc[5]=MC
#   
#   # HLT = HLT, not HLT2
#   ds[5]=/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
#   name[5]=ttbar
#   ismc[5]=MC
#   
#   
#   # HLT = HLT, not HLT2
#   ds[5]=/TTToSemiLeptonic_TuneCUETP8M1T4_alphaS01108_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
#   name[5]=ttbar
#   ismc[5]=MC
#   
#   # HLT = HLT, not HLT2
#   ds[5]=/TTTo2L2Nu_TuneCUETP8M1T4_alphaS01108_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
#   name[5]=ttbar
#   ismc[5]=MC
#   
#   # HLT = HLT, not HLT2
#   ds[5]=/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM
#   name[5]=ttbar
#   ismc[5]=MC
#   
#   # HLT = HLT, not HLT2
#   ds[5]=/TTToSemiLeptonic_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v2/MINIAODSIM
#   name[5]=ttbar
#   ismc[5]=MC


# 
# 
# ds[8]=/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# name[8]=ttbarmadgraphpythia8
# ismc[8]=MC
# 
# ds[9]=/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# name[9]=ttbaramcatnlopythia8
# ismc[9]=MC
# 
# 
# name[10]=SingleTop
# ds[10]=/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ismc[10]=MC
# 
# 
# name[11]=SingleAntiTop
# ds[11]=/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ismc[11]=MC
# 
# name[12]=ww
# ds[12]=/WW_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ismc[12]=MC
# 
# 
# name[13]=wz
# ds[13]=/WZ_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ismc[13]=MC
# 
# name[14]=zz
# ds[14]=/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ismc[14]=MC
# 
# 
# name[15]=powhegherwig
# ds[15]=/TT_TuneEE5C_13TeV-powheg-herwigpp/RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ismc[15]=MC
# 
# 
# 
# 
# 
# 
# ds[18]=/TTToSemiLeptonic_13TeV-powheg/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/MINIAODSIM
# name[18]=ttbarsemilep
# ismc[18]=MC
# 
# 
# 
# ds[19]=/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ds[20]=/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ds[21]=/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ds[22]=/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# name[19]=QCDHT100to200
# name[20]=QCDHT200to300
# name[21]=QCDHT300to500
# name[22]=QCDHT500to700
# ismc[19]=MC
# ismc[20]=MC
# ismc[21]=MC
# ismc[22]=MC
# 
# 
# ds[23]=/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ds[24]=/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ds[25]=/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# ds[26]=/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
# name[23]=QCDHT700to1000
# name[24]=QCDHT1000to1500
# name[25]=QCDHT1500to2000
# name[26]=QCDHT2000toInf
# ismc[23]=MC
# ismc[24]=MC
# ismc[25]=MC
# ismc[26]=MC



for idx in ` echo $JobIndexList `
do

echo ${name[${idx}]}


if [ ${ismc[$idx]} == "MC" ]
then

if [ ` echo ${name[$idx]} | grep ttbar ` ] 
then
cat crabconfig_template.py | sed "s|XXXXX|${name[${idx}]}|g" | sed "s|YYYYY|${nickname}|g" | sed "s|ZZZZZ|${ds[$idx]}|g" | \
                             sed "s|QQQQQ|MCTTBAR|g"> __CRABCONFIG__${name[${idx}]}.py
else
cat crabconfig_template.py | sed "s|XXXXX|${name[${idx}]}|g" | sed "s|YYYYY|${nickname}|g" | sed "s|ZZZZZ|${ds[$idx]}|g" | \
                             sed "s|QQQQQ|${ismc[$idx]}|g"> __CRABCONFIG__${name[${idx}]}.py

fi

else
cat crabconfig_template.py | sed "s|XXXXX|${name[${idx}]}|g" | sed "s|YYYYY|${nickname}|g" | sed "s|ZZZZZ|${ds[$idx]}|g" | \
                             sed "s|QQQQQ|${ismc[$idx]}|g" | \
   sed "s|#PPPPP#|config.Data.lumiMask=\"../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt\"|g" \
    > __CRABCONFIG__${name[${idx}]}.py
fi



done