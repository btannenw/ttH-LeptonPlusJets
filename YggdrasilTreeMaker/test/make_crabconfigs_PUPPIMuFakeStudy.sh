


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
cat yggdrasil_treeMaker_cfg.py | sed "s|isMC=True|isMC=False|g" | sed "s|XXXPERIODXXX|2017${P}|g"> __yggdrasil_treeMaker_DATA_${P}_cfg.py 
done


nickname="Yggdrasil_PUPPIMuonIsolationFakeStudy_2018_04_30"

JobIndexList=""

ds[1]=/DoubleMuon/Run2017B-17Nov2017-v1/MINIAOD
ds[2]=/DoubleMuon/Run2017C-17Nov2017-v1/MINIAOD
ds[3]=/DoubleMuon/Run2017D-17Nov2017-v1/MINIAOD
ds[4]=/DoubleMuon/Run2017E-17Nov2017-v1/MINIAOD
ds[5]=/DoubleMuon/Run2017F-17Nov2017-v1/MINIAOD



name[1]=DataMuB
name[2]=DataMuC
name[3]=DataMuD
name[4]=DataMuE
name[5]=DataMuF

ismc[1]=DATA_B
ismc[2]=DATA_C
ismc[3]=DATA_D
ismc[4]=DATA_E
ismc[5]=DATA_F

JobIndexList=${JobIndexList}" 1 2 3  4   5  "


ds[50]=/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
name[50]=ZjetIncl
ismc[50]=MC

JobIndexList=${JobIndexList}" 50 "

name[61]=wz
ds[61]=/WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[61]=MC

name[62]=zz
ds[62]=/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
ismc[62]=MC

JobIndexList=${JobIndexList}" 61 62 "



ds[110]=/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM
name[110]=ttto2l2nu
ismc[110]=MC


JobIndexList=${JobIndexList}" 110 "



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
