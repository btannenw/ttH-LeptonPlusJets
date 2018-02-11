


echo ""
echo ""
echo "# - - - - - - - "
echo remove large files from directory to submit jobs.
rm -rf ../../AnalysisCode/data/*
echo " Removing large file is done..."
echo " Removing files ${CMSSW_BASE}/lib/slc6_amd64_gcc530/proclib/*openloop*so because they are too large."
rm -rf ${CMSSW_BASE}/lib/slc6_amd64_gcc530/proclib/*openloop*so
echo "# - - - - - - "

cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=False|isMC=True|g" > __yggdrasil_treeMaker_MC_cfg.py 
diff yggdrasil_treeMaker_cfg.py __yggdrasil_treeMaker_MC_cfg.py

cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=False|isMC=True|g" | sed "s|isTTBARMC=False|isTTBARMC=True|g" > __yggdrasil_treeMaker_MCTTBAR_cfg.py 
diff yggdrasil_treeMaker_cfg.py __yggdrasil_treeMaker_MCTTBAR_cfg.py


for P in B C D E F
do
cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=True|isMC=False|g" | sed "s|XXXPERIODXXX|2016${P}|g"> __yggdrasil_treeMaker_DATA_${P}_cfg.py 
done


nickname="Satoshi_jobsubmit_mkdir_2018_01_22__test3"

JobIndexList=""

ds[1]=/SingleElectron/Run2017B-17Nov2017-v1/MINIAOD
ds[2]=/SingleElectron/Run2017C-17Nov2017-v1/MINIAOD
ds[3]=/SingleElectron/Run2017D-17Nov2017-v1/MINIAOD
ds[4]=/SingleElectron/Run2017E-17Nov2017-v1/MINIAOD
ds[5]=/SingleElectron/Run2017F-17Nov2017-v1/MINIAOD



name[1]=DataElB
name[2]=DataElC
name[3]=DataElD
name[4]=DataElE
name[5]=DataElF

ismc[1]=DATA_B
ismc[2]=DATA_C
ismc[3]=DATA_D
ismc[4]=DATA_E
ismc[5]=DATA_F

JobIndexList=${JobIndexList}" 1 2 3  4   5  "



ds[11]=/SingleMuon/Run2017B-17Nov2017-v1/MINIAOD
ds[12]=/SingleMuon/Run2017C-17Nov2017-v1/MINIAOD
ds[13]=/SingleMuon/Run2017D-17Nov2017-v1/MINIAOD
ds[14]=/SingleMuon/Run2017E-17Nov2017-v1/MINIAOD
ds[15]=/SingleMuon/Run2017F-17Nov2017-v1/MINIAOD


name[11]=DataMuB
name[12]=DataMuC
name[13]=DataMuD
name[14]=DataMuE
name[15]=DataMuF

ismc[11]=DATA_B
ismc[12]=DATA_C
ismc[13]=DATA_D
ismc[14]=DATA_E
ismc[15]=DATA_F

JobIndexList=${JobIndexList}" 11 12 13 14 15 "



ds[30]=/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#      /TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
name[30]=ttbar
ismc[30]=MC

# 
# ds[31]=/TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[31]=ttbarISRUp
# ismc[31]=MC
# 
# ds[32]=/TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[32]=ttbarISRDown
# ismc[32]=MC
# 
# ds[33]=/TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[33]=ttbarFSRUp
# ismc[33]=MC
# 
# ds[34]=/TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[34]=ttbarFSRDown
# ismc[34]=MC
# 
# 
# 
# ds[35]=/TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# 
# name[35]=ttbarUEup
# ismc[35]=MC
# 
# ds[36]=/TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[36]=ttbarUEdown
# ismc[36]=MC
# 

JobIndexList=${JobIndexList}" 30 "

ds[40]=/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
name[40]=WjetIncl
ismc[40]=MC


JobIndexList=${JobIndexList}" 40 "


ds[50]=/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
name[50]=ZjetIncl
ismc[50]=MC

# ismc[51]=MC
# ismc[52]=MC
# ismc[53]=MC
# ismc[54]=MC
# ismc[55]=MC
# ismc[56]=MC
# ismc[57]=MC
# ismc[58]=MC
# 
# name[51]=DYM50HT70to100
# name[52]=DYM50HT100to200
# name[53]=DYM50HT200to400
# name[54]=DYM50HT400to600
# name[55]=DYM50HT600to800
# name[56]=DYM50HT800to1200
# name[57]=DYM50HT1200to2500
# name[58]=DYM50HT2500toInf
# 
# ds[51]=/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[52]=/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[53]=/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[54]=/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[55]=/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM
# ds[56]=/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[57]=/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[58]=/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# 


ds[59]=/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM
name[59]=ZjetLowMass
ismc[59]=MC


#JobIndexList=${JobIndexList}" 50 51 52 53 54 55 56 57 58 "
JobIndexList=${JobIndexList}" 50 59 "




name[60]=ww
ds[60]=/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[60]=MC

name[61]=wz
ds[61]=/WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[61]=MC

name[62]=zz
ds[62]=/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[62]=MC

JobIndexList=${JobIndexList}" 60 61 62 "



name[70]=tchan_top
ds[70]=/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[70]=MC

name[71]=tchan_tbar
ds[71]=/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[71]=MC

name[72]=tW
ds[72]=/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[72]=MC

name[73]=tbarW
ds[73]=/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[73]=MC

name[74]=schan_both
ds[74]=/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[74]=MC

JobIndexList=${JobIndexList}" 70 71 72 73 74  "




#   ds[80]=/QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#   ds[81]=/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
#   ds[82]=/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#   ds[83]=/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#   ds[84]=/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM
#   ds[85]=/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#   ds[86]=/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#   ds[87]=/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#   ds[88]=/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
#   
#   name[80]=QCDHT50to100
#   name[81]=QCDHT100to200
#   name[82]=QCDHT200to300
#   name[83]=QCDHT300to500
#   name[84]=QCDHT500to700
#   name[85]=QCDHT700to1000
#   name[86]=QCDHT1000to1500
#   name[87]=QCDHT1500to2000
#   name[88]=QCDHT2000toInf
#   
#   ismc[80]=MC
#   ismc[81]=MC
#   ismc[82]=MC
#   ismc[83]=MC
#   ismc[84]=MC
#   ismc[85]=MC
#   ismc[86]=MC
#   ismc[87]=MC
#   ismc[88]=MC
#   
#   JobIndexList=${JobIndexList}" 80 81 82 83 84 85 86 87 88 " 
#   




ds[100]=/ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
name[100]=tthbb
ismc[100]=MC

ds[101]=/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
name[101]=tthnonbb
ismc[101]=MC

JobIndexList=${JobIndexList}" 100 101" 

ds[109]=/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
name[109]=tttofullhad
ismc[109]=MC


ds[110]=/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM
name[110]=ttto2l2nu
ismc[110]=MC

ds[111]=/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
name[111]=tttosemilep
ismc[111]=MC


JobIndexList=${JobIndexList}" 109 110 111 "


# ismc[112]=MC
# ismc[113]=MC
# ismc[114]=MC
# ismc[115]=MC
# ismc[116]=MC
# ismc[117]=MC
# ismc[118]=MC
# ismc[119]=MC
# ismc[120]=MC
# ismc[121]=MC
# ismc[122]=MC
# ismc[123]=MC
# 
# 
# ds[112]=/TT_TuneCUETP8M2T4_mtop1735_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# ds[113]=/TT_TuneCUETP8M2T4_mtop1715_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[112]=ttbarmass1735
# name[113]=ttbarmass1715
# 
# ds[114]=/TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# ds[115]=/TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[114]=ttbartuneup
# name[115]=ttbartunedown
# 
# ds[116]=/TT_hdampUP_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[116]=ttbarhdampup
# 
# ds[117]=/TT_hdampDOWN_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[117]=ttbarhdampdown
# 
# ds[118]=/TT_widthx0p2_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[118]=ttbarwidth0p2
# 
# ds[119]=/TT_widthx0p5_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[119]=ttbarwidth0p5
# 
# ds[120]=/TT_widthx4_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[120]=ttbarwithd4
# 
# ds[121]=/TT_widthx8_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[121]=ttbarwidth8
# 
# ds[122]=/TT_TuneCUETP8M2T4_13TeV-powheg-colourFlip-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
# name[122]=ttbarcolorflip
# 
# ds[123]=/TT_TuneCUETP8M2T4_erdON_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM
# name[123]=ttbarerdon
# 
# 
# JobIndexList=${JobIndexList}" 112 113 114 115 116 117 118 119 120 121 122 123  "



ds[150]=/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ds[152]=/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ds[153]=/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ds[156]=/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM

name[150]=TTWJetsToLNu_1
name[151]=TTWJetsToLNu_2
name[152]=TTWJetsToQQ
name[153]=TTZToLLNuNu_1
name[154]=TTZToLLNuNu_2
name[155]=TTZToLLNuNu_3
name[156]=TTZToQQ

ismc[150]=MC
ismc[151]=MC
ismc[152]=MC
ismc[153]=MC
ismc[154]=MC
ismc[155]=MC
ismc[156]=MC

JobIndexList=${JobIndexList}" 150 152 153 156 " 




for idx in ` echo $JobIndexList `
do

echo ${name[${idx}]}


if [ ${ismc[$idx]} == "MC" ]
then

if [ ` echo ${name[$idx]} | grep 'ttbar\|ttto' ` ] 
then
cat crabconfig_template.py | sed "s|XXXXX|${name[${idx}]}|g" | sed "s|YYYYY|${nickname}|g" | sed "s|ZZZZZ|${ds[$idx]}|g" | \
                             sed "s|QQQQQ|MCTTBAR|g"> __CRABCONFIG__${name[${idx}]}.py
else
cat crabconfig_template.py | sed "s|XXXXX|${name[${idx}]}|g" | sed "s|YYYYY|${nickname}|g" | sed "s|ZZZZZ|${ds[$idx]}|g" | \
                             sed "s|QQQQQ|${ismc[$idx]}|g"> __CRABCONFIG__${name[${idx}]}.py

fi

else

# Data : 

JsonFile=Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt




cat crabconfig_template.py | sed "s|XXXXX|${name[${idx}]}|g" | sed "s|YYYYY|${nickname}|g" | sed "s|ZZZZZ|${ds[$idx]}|g" | \
                             sed "s|QQQQQ|${ismc[$idx]}|g" | \
   sed "s|#PPPPP#|config.Data.lumiMask=\"../data/${JsonFile}\"|g" \
    > __CRABCONFIG__${name[${idx}]}.py
fi



done
