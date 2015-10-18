#!/usr/local/bin/perl

####make DIRs
@dir = ("Condor", "Script", "Log", "Error", "Output");   

for $dir (@dir){
    if(! -d $dir ){
	mkdir $dir
    }
}

####################
#### Modify $script and $workingDir if necessary
###################
$script = "yggdrasil_treeReader_13TeV";       # Name of job

#$workingDir = "/uscms_data/d2/lwming/CMSSW_7_2_3/src/ttH-LeptonPlusJets/AnalysisCode"; #### change this
#$workingDir = "/uscms_data/d3/flowers1/May2015/CMSSW_7_2_3/src/ttH-LeptonPlusJets/AnalysisCode"; #### change this
$workingDir = "/uscms_data/d1/satoshi/ohio/CMSSW_7_4_4_patch4/src/ttH-LeptonPlusJets/AnalysisCode"; #### change this

$SampleType = 2;
$SplitType   = -1;
$Njobs  = 10;

$num = @ARGV;

if( $num >= 1 ){
    $SampleType = $ARGV[0];

    if( $num>=2 ){
	$SplitType = $ARGV[1];

	if( $num>=3 ){
	    $Njobs = $ARGV[2];
	}
    }
}



open SHFILE, "> Script/condor\_$script.sh";
print SHFILE "#!/bin/sh\n";
print SHFILE "\n";
print SHFILE "echo $1\n";
print SHFILE "\n";
print SHFILE "echo \"\"\n";
print SHFILE "echo \"Using ROOT on Condor\"\n";
print SHFILE "echo \"\"\n";
print SHFILE "cd \${_CONDOR_SCRATCH_DIR}\n";
print SHFILE "\n";
print SHFILE "SampleType=\$1\n";
print SHFILE "SplitType=\$2\n";
print SHFILE "NumEvents=\$3\n";
print SHFILE "NumJobs=\$4\n";
print SHFILE "jobN=\$5+1\n";
print SHFILE "\n";
print SHFILE "root -b -q $workingDir/macros/head13TeV.C '$workingDir/macros/$script.C+('\$SampleType','\$SplitType','\$NumEvents','\$NumJobs','\$jobN')'\n";
print SHFILE "\n";
close SHFILE;


#%samples = ("ttH" => 1, "ttJets" => 2, "ttW" => 3, "ttZ" => 4);
#%splits = ("ttb" => 1, "ttbb" => 2, "tt2b" => 3, "ttcc" => 4, "ttlf" => 5);

for $iSample (1, 2, 3, 4){
    $Njobs = 1;
    $SampleType = $iSample; 
    if ($iSample == 2){
	$Njobs = 10;
    }
    for $iSplit( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11){
	$SplitType = $iSplit;
	if ($iSample >2  && $iSplit !=1){
	    next;
	    }
	    if($iSample >1 && $iSplit>5){
	    next
	}

open CONDORFILE, "> Condor/condor\_$script.jdl";
print CONDORFILE "# A Condor submission file\n";
print CONDORFILE "Executable              = Script/condor\_$script.sh\n";
print CONDORFILE "Universe                = vanilla\n";
print CONDORFILE "Getenv                  = true\n";
print CONDORFILE "\n";
print CONDORFILE "Arguments               = $SampleType $SplitType -1 $Njobs \$(Process)\n";
print CONDORFILE "Output                  = Output/condor\_$SampleType\_$SplitType\_\$(Process).out\n";
print CONDORFILE "Error                   = Error/condor\_$SampleType\_$SplitType\_\$(Process).err\n";
print CONDORFILE "Log                     = Log/condor\_$SampleType\_$SplitType\_\$(Process).log\n";
print CONDORFILE "\n";
print CONDORFILE "use_x509userproxy = true\n";
print CONDORFILE "Should_Transfer_Files   = YES\n";
print CONDORFILE "When_To_Transfer_Output = ON_EXIT\n";
#print CONDORFILE "Transfer_Input_Files = filelist_DYJetsToLL_M_50_13TeV_madgraph_pythia8_PU20bx25_POSTLS170_V5_v1.txt, filelist_TTJets_MSDecaysCKM_central_Tune4C_13TeV_madgraph_PU20bx25_POSTLS170_V5_v1.txt\n";
print CONDORFILE "\n";
print CONDORFILE "#+IsLocalJob             = true\n";
print CONDORFILE "#Rank                    = TARGET.IsLocalSlot\n";
print CONDORFILE "\n";
print CONDORFILE "Queue $Njobs\n";
print CONDORFILE "\n";
close CONDORFILE;

system("chmod a+x Script/condor\_$script.sh");
print "submitting: condor_submit Condor/condor\_$script.jdl\n";
system("condor_submit Condor/condor\_$script.jdl");

    }
}
